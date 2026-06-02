#!/usr/bin/env python
"""
Determine the angles between the crystal facets and the electric field 
vector.

Also print the electric-field vector in the crystal frame.
"""

import argparse
import itertools
import pandas as pd
import numpy as np
from regroup import FrameGeometry
from regroup import ExptList
from cctbx import sgtbx
from cctbx.sgtbx import subgroups


def angle(v1, v2):
    """Compute angle between two vectors"""
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))


def get_normal_vector(hkl, Astar):
    """
    Get normal vector to real-space Miller plane, hkl.

    Note
    ----
        The normal vector to a real-space Miller plane is collinear with
        the reciprocal dHKL vector. As such, we can use this simpler 
        formula in the reciprocal lattice basis to get the correct 
        orientation in real-space. For a graphical explanation of this,
        please look at Rupp, p238.
    """
    return hkl @ Astar.T


def lab_vec_to_crystal(v_lab, Astar):
    """
    Convert a lab-frame vector into direct crystal fractional coordinates.

    Astar maps reciprocal fractional hkl -> lab Cartesian reciprocal vector.
    Therefore A = inv(Astar).T maps direct fractional -> lab Cartesian,
    and v_frac = inv(A) @ v_lab = Astar.T @ v_lab.
    """
    v_lab = np.array(v_lab, dtype=float)
    v_cryst = Astar.T @ v_lab
    return v_cryst / np.linalg.norm(v_cryst)


def facet_normal_to_crystal_frame(hkl, O):
    """
    Convert a reciprocal-lattice facet normal hkl into a direct crystal-frame
    vector parallel to the real-space plane normal.

    This is the same transform already used in get_spacegroup().
    """
    hkl = np.array(hkl, dtype=float)
    v = np.linalg.inv(O) @ (np.linalg.inv(O).T @ hkl)
    return v / np.linalg.norm(v)


def fmt_vec(v, ndigits=4):
    """Pretty-print a vector as a rounded tuple."""
    v = np.array(v, dtype=float)
    return tuple(np.round(v, ndigits))


def get_spacegroup(facet, parent_sg, O):
    """
    Get reduced symmetry spacegroup for crystal aligned on `facet`.
    """
    e_field_unit_vector = facet_normal_to_crystal_frame(facet, O)

    parent = sgtbx.space_group_info(parent_sg)
    subgrs = subgroups.subgroups(parent).groups_parent_setting()

    possible = []
    for subgroup in subgrs:
        subgroup_info = sgtbx.space_group_info(group=subgroup)
        valid = True
        for op in subgroup.smx():
            rot_mat = np.array(op.r().as_double()).reshape((3, 3))
            valid &= np.allclose(rot_mat @ e_field_unit_vector, e_field_unit_vector)

        if valid:
            possible.append([subgroup.n_smx(), subgroup_info.symbol_and_number(), subgroup])

    possible = sorted(possible)
    return possible[-1][1], possible[-1][0]

def mean_vec(values, ndigits=4):
    """
    Average a pandas Series of tuple/list/array vectors.
    """
    arr = np.array([np.array(v, dtype=float) for v in values])
    v = arr.mean(axis=0)
    return fmt_vec(v, ndigits=ndigits)

def run_regroup(inp, spacegroup, hmax=1, efvector=(0, -1, 0), filename=None):
    """
    Core regroup logic, usable from Python without argparse.
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

    for facet in facets:
        hkl = np.array(facet)

        if np.gcd.reduce(hkl) > 1:
            continue

        for i, Astar in enumerate(Astars):
            normal = get_normal_vector(hkl, Astar)
            theta = np.rad2deg(angle(normal, efvector))

            ef_cryst = lab_vec_to_crystal(efvector, Astar)

            l_facets.append(facet)
            l_images.append(images[i])
            l_angles.append(theta)
            l_ef_cryst.append(fmt_vec(ef_cryst))

    if precog:
        geom = FrameGeometry(inp[0])
        O = geom.get_orthogonalization_matrix().T
    elif dials:
        O = dials_expts.get_orthogonalization_matrix().T

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

    results["facet_normal_crystal"] = results.Facet.apply(
        lambda hkl: fmt_vec(facet_normal_to_crystal_frame(hkl, O))
    )

    results["spacegroup"], results["n_symops"] = zip(
        *results.Facet.apply(get_spacegroup, parent_sg=spacegroup, O=O)
    )

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more 
        print(results)
        if filename:
            with open(filename, 'w') as fname:
                print(results, file=fname)

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
        help="Maximal number to include in a Miller plane",
        type=int,
    )
    parser.add_argument(
        "-ef",
        "--efvector",
        nargs=3,
        type=float,
        default=(0, -1, 0),
        metavar=("efx", "efy", "efz"),
        help="EF vector in lab frame",
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
    )

if __name__ == "__main__":
    main()