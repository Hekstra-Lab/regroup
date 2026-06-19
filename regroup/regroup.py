#!/usr/bin/env python3
"""Minimal regroup CLI for the PDZ3 workflow.

Determine the angles between crystal facets and an input field
vector.

Also print the field vector in the crystal frame as well as 
the field-symmetry alignment (FSA) for each broken symop, 
where the FSA is defined as the dot of the unit field vector
and the symop-transformed unit field vector. 

All 1D arrays are treated as 1D row vectors. 
"""

import argparse
import itertools
import pandas as pd
import gemmi
import numpy as np
from regroup import FrameGeometry
from regroup.framegeometry import get_orthogonalization_matrix
from regroup import ExptList
from regroup.geom_utils import *
from cctbx import sgtbx
from cctbx.sgtbx import subgroups
from regroup.low_sym import cctbx_cb_op_to_rs_op
from regroup.write_subgroup_json import (
    extract_basis_change_op,
    subgroup_record,
    write_subgroup_json,
)

def print_fsa_table(parent_sg, vec=None, O=None, file=None, den=None, opnums=None):
    """
    Print FSA table using Gemmi operation ordering and 
    crystal Euclidean space metric tensor.
    """
    sg = gemmi.SpaceGroup(str(parent_sg))
    ops = sg.operations().sym_ops

    field = np.asarray(vec, dtype=float)

    G = O.T @ O

    field_display = np.array(field, dtype=float)/np.linalg.norm(field)
    
    if opnums is None:
        opnums = list(range(len(ops)))
    elif isinstance(opnums, str):
        opnums = [int(x) for x in opnums.replace(",", " ").split()]
    else:
        opnums = [int(x) for item in opnums for x in str(item).replace(",", " ").split()]

    print(file=file)
    print(f"Space group: {parent_sg}", file=file)
    print(f"Field vector, crystal frame: {np.round(field_display, 6)}", file=file)
    print(file=file)

    header = (
        f"{'fsa':>10s}  "
        f"{'angle_deg':>10s}  "
        f"{'broken?':>8s}  "
        f"{'opnum':>6s}  "
        f"{'symop':25s}"
    )
    print(header, file=file)
    print("-" * len(header), file=file)

    for opnum in opnums:
        if opnum < 0 or opnum >= len(ops):
            print(
                f"Skipping invalid opnum {opnum}; valid range is 0 to {len(ops) - 1}.",
                file=file,
            )
            continue

        op = ops[opnum]
        R = gemmi_to_rot(op, den=den)

        rotated = R @ field

        # Metric-aware version of fsa = -v dot Rv
        fsa = -cosine_metric(field, rotated, G)
        angle_deg = angle_metric(field, rotated, G)

        # not broken means the rotated fractional vector represents the same
        # physical direction as the original. Use metric norm for tolerance.
        broken = not metric_close(rotated, field, G, atol=1e-5)

        print(
            f"{fsa:10.5f}  "
            f"{angle_deg:10.3f}  "
            f"{str(broken):>8s}  "
            f"{opnum:6d}  "
            f"{op.triplet():25s}",
            file=file,
        )

def _op_rot(op):
    return np.array(op.rot, dtype=float) / float(getattr(op, "DEN", 24))


def get_spacegroup(facet, parent_sg, O, subgroup_index=None, return_candidates=False):
    """
    Get reduced symmetry spacegroup for crystal aligned on `facet`.

    If subgroup_index is supplied, choose that exact cctbx subgroup index, but
    only if it preserves the field/facet vector. The returned subgroup operations
    are in the parent setting, which is the right thing to serialize for cosets.
    """

    guess_e_field_unit_vector = facet_normal_to_crystal_frame(facet, O)

    parent = sgtbx.space_group_info(parent_sg)
    subgrs = subgroups.subgroups(parent).groups_parent_setting()

    possible = []
    candidate_records = []
    for idx, subgroup in enumerate(subgrs):
        subgroup_info = sgtbx.space_group_info(group=subgroup)
        valid = True

        for op in subgroup.smx():
            rot_mat = np.array(op.r().as_double()).reshape((3, 3))
            valid &= np.allclose(rot_mat @ guess_e_field_unit_vector, guess_e_field_unit_vector)

        if valid:
            sg_symbol = subgroup_info.symbol_and_number()
            cb_op = extract_basis_change_op(sg_symbol)
            record = subgroup_record(idx, subgroup, subgroup_info=subgroup_info, facet=facet)
            possible.append([subgroup.n_smx(), str(sg_symbol), subgroup, cb_op, idx, record])
            candidate_records.append(record)

    if not possible:
        raise ValueError(f"No valid subgroups found for facet {facet}")

    if subgroup_index is not None:
        matches = [item for item in possible if item[4] == subgroup_index]
        if not matches:
            valid_indices = ", ".join(str(item[4]) for item in possible)
            raise ValueError(
                f"Requested subgroup index {subgroup_index} is not valid for facet {facet}. "
                f"Valid indices preserving the vector are: {valid_indices}"
            )
        chosen = matches[0]
    else:
        possible = sorted(possible, key=lambda item: (item[0], item[4]))
        chosen = possible[-1]

    result = [chosen[1], chosen[0], chosen[2], chosen[3], chosen[4], chosen[5]]
    if return_candidates:
        result.append(candidate_records)
    return result

def mean_vec(values, ndigits=4):
    """
    Average a pandas Series of tuple/list/array vectors.
    """
    arr = np.array([np.array(v, dtype=float) for v in values])
    v = arr.mean(axis=0)
    return fmt_vec(v, ndigits=ndigits)

def run_regroup(inp, spacegroup, hmax=1, efvector=(0, -1, 0), filename=None, fsa=False, opnums=None, subgroup_index=None, export_subgroup_json=None, export_candidates_json=None):
    """
    Computes A matrix and angle between vector and facet normals.
    We deal with four coordinate frames: 
        - the lab Cartesian frame. 
        - the crystal Cartesian frame.
        - the crystal fractional coordinate frame, 
          which we utilize for symmetry breaking and visualization. 
        - the reciprocal lattice. 
    Given these, we rely on the following formal statements: 
        - Astar transforms facet normals of reciprocal lattice hkls 
          into the lab Cartesian frame.
            - This allows us to calculate angles between the E field
              and facet normals in the lab frame. 
        - The metric tensor O^TO transforms reciprocal lattice points, i.e.
          Miller plane normals, into fractional coordinates.  
            - This allows us to rotate/translate facet normals and 
              determine preserved/broken symmetries. 
              
    """
    facets = list(itertools.product(np.arange(-hmax, hmax + 1), repeat=3))
    facets.remove((0, 0, 0))

    l_facets = []
    l_images = []
    l_angles = []
    l_ef_cryst = []

    if spacegroup is None:
        raise ValueError("Please provide parent spacegroup with -sg / --spacegroup")

    precog = False
    dials = False
    Astars = []
    images = []

    if inp[0][-4:] == ".inp":
        precog = True
        for _inp in inp:
            geometry = FrameGeometry(_inp)
            Astars.append(geometry.get_reciprocal_Amatrix())
            images.append(_inp)

    dials_expts = None
    if inp[0][-5:] == ".expt":
        dials = True
        dials_expts = ExptList(inp[0])
        Astars.extend(dials_expts.get_reciprocal_Amatrices())
        images.extend(dials_expts.get_image_filenames())

    if not dials and not precog:
        raise ValueError("File extension unrecognized. Please enter .inp or .expt files.")

    efvector = np.array(efvector, dtype=float)
    
    if precog:
        geom = FrameGeometry(inp[0])
        O = geom.get_orthogonalization_matrix().T
    elif dials:
        O = dials_expts.get_orthogonalization_matrix().T
        
    for facet in facets:
        hkl = np.array(facet)

        if np.gcd.reduce(hkl) > 1:
            continue

        for i, Astar in enumerate(Astars):
            normal = get_normal_vector(hkl, Astar)
            theta = np.rad2deg(angle(normal, efvector))
            n_frac = facet_normal_to_crystal_frame(facet, O)

            #check that Euclidean coordinate angle computation 
            #matches fractional coordinate space.
            ef_cryst = lab_vec_to_crystal(efvector, Astar)
            theta1 = angle_metric(ef_cryst, n_frac, O.T @ O)
            assert np.allclose(theta, theta1, rtol=1e-2)

            l_facets.append(facet)
            l_images.append(images[i])
            l_angles.append(theta)
            l_ef_cryst.append(fmt_vec(ef_cryst))

    df = pd.DataFrame(
        {
            "Facet": l_facets,
            "Image": l_images,
            "Angle": l_angles,
            "ef_crystal": l_ef_cryst,
        }
    )

    results = df.groupby("Facet").agg(
        {
            "Angle": ["mean", "std", "count"],
            "ef_crystal": mean_vec,
        }
    )

    results.sort_values(("Angle", "mean"), inplace=True)
    results.reset_index(inplace=True)
    mv = results.loc[0, "ef_crystal"][0]
    print("Fractional coordinates of field vector:", mv)

    results["facet_normal_crystal"] = results.Facet.apply(
        lambda hkl: fmt_vec(facet_normal_to_crystal_frame(hkl, O))
    )

    sg_results = results.Facet.apply(get_spacegroup, parent_sg=spacegroup, O=O, subgroup_index=subgroup_index).tolist()
    results["spacegroup"] = [item[0] for item in sg_results]
    results["n_symops"] = [item[1] for item in sg_results]
    results["basis_change_op"] = [item[3] for item in sg_results]
    results["cctbx_subgroup_index"] = [item[4] for item in sg_results]

    if export_subgroup_json:
        best_facet = results.loc[0, "Facet"][0]
        best_result = get_spacegroup(
            best_facet,
            parent_sg=spacegroup,
            O=O,
            subgroup_index=subgroup_index,
        )
        write_subgroup_json(
            export_subgroup_json,
            parent_sg=spacegroup,
            selected_record=best_result[5],
        )
        print(f"Wrote selected subgroup JSON: {export_subgroup_json}")
    
    fsa_vec = np.array(mv)

    #compute and display change-of-basis op and transformed hkl.
    best_cb_op = sg_results[0][3]
    cb = cctbx_cb_op_to_rs_op(best_cb_op)
    fv = results.loc[0, "Facet"][0]
    new_fv_hkl = fv @ _op_rot(cb)

    results = results.drop(columns=["ef_crystal", "basis_change_op"], level=0)

    print("Best-match basis-change op:", best_cb_op)
    print("Best-match facet:", fv)
    print("Transformed best-match facet:", new_fv_hkl)

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more 
        print(results)
        if fsa:
            print_fsa_table(spacegroup, fsa_vec, O=O, opnums=opnums)

        if filename:
            with open(filename, 'w') as fname:
                print("Best-match basis-change op:", best_cb_op)
                print("Best-match facet:", fv)
                print("Transformed best-match facet:", new_fv_hkl)
                print(results, file=fname)
                if fsa:
                    print_fsa_table(spacegroup, fsa_vec, O=O, opnums=opnums, file=fname)

    return results

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "inp",
        nargs="+",
        help="Precognition geometry file, .mccd.inp, DIALS experiment file, .expt, or facet, if known",
    )
    parser.add_argument("-sg", "--spacegroup", type=int, help="Parent spacegroup")
    parser.add_argument(
        "--hmax",
        default=1,
        help="Maximal index in candidate Miller planes",
        type=int,
    )
    parser.add_argument(
        "-ef",
        "--efvector",
        nargs=3,
        type=float,
        default=(0, 1, 0),
        metavar=("efx", "efy", "efz"),
        help="field vector in lab frame",
    )
    parser.add_argument(
        "--fsa",
        action="store_true",
        help="Print field-symmetry alignment table for broken symops in the top-rated low-symmetry spacegroup",
    )
    parser.add_argument(
        "--filename",
        default=None,
        help="Filename for saved output",
        type=str,
    )
    parser.add_argument(
        "--export-subgroup-json",
        default=None,
        help=(
            "Write the selected cctbx subgroup as JSON, including exact parent-setting symops "
            "for basis_change_pdb.py --child-subgroup-json."
        ),
    )

    args = parser.parse_args()

    run_regroup(
        inp=args.inp,
        spacegroup=args.spacegroup,
        hmax=args.hmax,
        efvector=args.efvector,
        filename=args.filename,
        fsa=args.fsa,
        export_subgroup_json=args.export_subgroup_json,
    )

if __name__ == "__main__":
    main()