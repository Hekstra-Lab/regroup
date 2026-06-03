#!/usr/bin/env python
"""
Determine the angles between the crystal facets and an input field 
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
from regroup import ExptList
from regroup.geom_utils import *
from cctbx import sgtbx
from cctbx.sgtbx import subgroups
from regroup.low_sym import _has_any, cctbx_cb_op_to_rs_op

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

def _extract_basis_change_op(symbol_and_number):
    """
    Extract the regroup/cctbx basis-change operator embedded in cctbx's
    symbol_and_number() display. If no operator is present, use identity.
    """
    s = str(symbol_and_number)
    if " (No." in s:
        s = s.split(" (No.", 1)[0]
    if "(" in s and ")" in s:
        return s.rsplit("(", 1)[1].split(")", 1)[0].strip()
    return "x,y,z"

def _op_rot(op):
    return np.array(op.rot, dtype=float) / float(getattr(op, "DEN", 24))


def get_spacegroup(facet, parent_sg, O):
    """
    Get reduced symmetry spacegroup for crystal aligned on `facet`.
    """

    #generate a fractional E field guess, on which to apply rotation matrices
    #for checking subgroup validity. 
    guess_e_field_unit_vector = facet_normal_to_crystal_frame(facet, O)

    parent = sgtbx.space_group_info(parent_sg)
    subgrs = subgroups.subgroups(parent).groups_parent_setting()

    possible = []
    for subgroup in subgrs:
        subgroup_info = sgtbx.space_group_info(group=subgroup)
        valid = True

        #check that all the subgroup operations preserve the vector orientation. 
        for op in subgroup.smx():
            rot_mat = np.array(op.r().as_double()).reshape((3, 3))
            valid &= np.allclose(rot_mat @ guess_e_field_unit_vector, guess_e_field_unit_vector)
            # is the operation above correct?
        if valid:
            sg_symbol = subgroup_info.symbol_and_number()
            cb_op = _extract_basis_change_op(sg_symbol)
            possible.append([subgroup.n_smx(), sg_symbol, subgroup, cb_op])

    possible = sorted(possible)
    return possible[-1][1], possible[-1][0], possible[-1][2], possible[-1][3]

def mean_vec(values, ndigits=4):
    """
    Average a pandas Series of tuple/list/array vectors.
    """
    arr = np.array([np.array(v, dtype=float) for v in values])
    v = arr.mean(axis=0)
    return fmt_vec(v, ndigits=ndigits)

def run_regroup(inp, spacegroup, hmax=1, efvector=(0, -1, 0), filename=None, fsa=False, opnums=None):
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

    sg_results = results.Facet.apply(get_spacegroup, parent_sg=spacegroup, O=O).tolist()
    results["spacegroup"] = [item[0] for item in sg_results]
    results["n_symops"] = [item[1] for item in sg_results]
    results["basis_change_op"] = [item[3] for item in sg_results]
    
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
        help="Precognition geometry file, .mccd.inp, or DIALS experiment file, .expt",
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
        "--opnums",
        nargs="+",
        default=None,
        help="FSA parent symop indices. Accepts space- or comma-separated values.",
    )
    parser.add_argument(
        "--filename",
        default=None,
        help="Filename for saved output",
        type=str,
    )

    args = parser.parse_args()

    run_regroup(
        inp=args.inp,
        spacegroup=args.spacegroup,
        hmax=args.hmax,
        efvector=args.efvector,
        filename=args.filename,
        fsa=args.fsa,
        opnums=args.opnums,
    )

if __name__ == "__main__":
    main()