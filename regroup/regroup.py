#!/usr/bin/env python
"""
Determine the angles between the crystal facets and the electric field 
vector.
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
    return hkl@Astar.T

def get_spacegroup(facet, parent_sg, O):
    """
    Get reduced symmetry spacegroup for crystal aligned on `facet`.
    """
    e_field_direction = np.linalg.inv(O)@(np.linalg.inv(O).T@facet)
    e_field_unit_vector = np.array(e_field_direction)/np.linalg.norm(e_field_direction)
    
    parent = sgtbx.space_group_info(parent_sg)
    subgrs = subgroups.subgroups(parent).groups_parent_setting()

    possible = []
    for subgroup in subgrs:
        subgroup_info = sgtbx.space_group_info(group=subgroup)
        valid = True
        for op in subgroup.smx():
            rot_mat = np.array(op.r().as_double()).reshape((3, 3))
            valid &= np.allclose(rot_mat@e_field_unit_vector, e_field_unit_vector)

        if valid:
            possible.append([subgroup.n_smx(), subgroup_info.symbol_and_number(), subgroup])

    possible = sorted(possible)
    return possible[-1][1], possible[-1][0]

def main():

    # CLI
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=__doc__)
    parser.add_argument("inp", nargs="+",
                        help="Precognition geometry file (suffixed with .mccd.inp) or DIALS experiment file (suffixed with .expt)")
    parser.add_argument("-sg", "--spacegroup", type=int, help="Parent spacegroup")
    parser.add_argument("--hmax", default=1, help="Maximal number to include in a Miller plane", type=int)
    parser.add_argument("-ef", "--efvector", nargs=3, type=int, default=(0, -1, 0), metavar=("efx", "efy", "efz"),
                        help="EF vector")
    parser.add_argument("--filename", default=None, help="Filename for saved output", type=str)
    
    args = parser.parse_args()

    # Relevant Miller planes
    facets = list(itertools.product(np.arange(-args.hmax, args.hmax+1), repeat=3))
    facets.remove((0, 0, 0))

    # Initialize results lists
    l_facets = []
    l_images = []
    l_angles = []

    # Load geometry
    if args.spacegroup:
        spacegroup = args.spacegroup

    # Initialize A* matrix list and filetypes
    precog = False
    dials = False
    Astars = []
    images = []

    # Process Precognition files
    if args.inp[0][-4:] == '.inp':
        precog = True
        for inp in args.inp:
            geometry = FrameGeometry(inp)
            Astars.append(geometry.get_reciprocal_Amatrix())
            images.append(inp)

    # Process DIALS files
    dials_expts = None
    if args.inp[0][-5:] == '.expt':
        dials = True
        dials_expts = ExptList(args.inp[0])
        Astars.extend(dials_expts.get_reciprocal_Amatrices())
        images.extend(dials_expts.get_image_filenames())

    if not dials and not precog:
        print('ERROR: File extension unrecognized. Please enter .inp or .expt files.')
        exit()

    for facet in facets:
        hkl = np.array(facet)

        # GH#1: Remove extra parallel Miller planes from output
        if np.gcd.reduce(hkl) > 1:
            continue
        
        # Loop over orientation matrices
        for i, Astar in enumerate(Astars):
            normal = get_normal_vector(hkl, Astar)
            theta = np.rad2deg(angle(normal, np.array(args.efvector)))
            l_facets.append(facet)
            l_images.append(images[i])
            l_angles.append(theta)

    # Use first frame for orthogonalization matrix
    O = None # Instantiate variable
    if precog:
        geom = FrameGeometry(args.inp[0])
        O = geom.get_orthogonalization_matrix().T
    elif dials:
        O = dials_expts.get_orthogonalization_matrix().T
    
    # Format output
    df = pd.DataFrame({"Facet": l_facets, "Image": l_images, "Angle": l_angles})
    results = df.groupby("Facet").agg({"Angle":["mean", "std", "count"]})
    results.sort_values(("Angle", "mean"), inplace=True)
    results.reset_index(inplace=True)
    results["spacegroup"], results["n_symops"] = zip(*results.Facet.apply(get_spacegroup, parent_sg=spacegroup, O=O))
    print(results)
    if args.filename:
        pd.options.display.max_rows = len(results)
        with open(args.filename, 'w') as fname:
            fname.write(results.__repr__())
    
if __name__ == "__main__":
    main()
