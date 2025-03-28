# regroup

Determine new space groups for analyzing pump-probe crystallography experiments

## Installation  

This program relies on `sgtbx` for the hierarchical grouping of different crystallographic space groups. Currently, it does not seem that `cctbx` can be easily installed with `pip`, so this dependency must be installed separately. Similarly, for processing DIALS files, `dxtbx` is required, which must also be installed. The following snippet should install the `regroup` command-line program to your current environment:

```shell
conda install -c conda-forge cctbx
conda install -c conda-forge dxtbx
pip install git+https://github.com/Hekstra-Lab/regroup.git
```

## Features  

`regroup` requires knowledge of the experimental geometry of the crystal in the lab reference frame in order to determine the new space group based on the orientation of the crystal relative to the "pump" perturbation. Since much of our group's work is conducted at the BioCARS Laue beamline (APS 14-ID-B), this program currently supports Precognition geometry files (`.inp` format) or DIALS experiment files (`.expt` format). Only DIALS stills can be processed -- scans are not handled currently.

If anyone is interested in support for additional file formats, please reach out by filing an [issue](https://github.com/Hekstra-Lab/regroup/issues).

For a full list of options and parameters, type `regroup --help` into your terminal.

## Conventions
Since `regroup` can process both DIALS experiment files and Precognition input files, it is important to make use of the correct lab frame convention. The program will automatically infer the lab frame convention given the type of file input, but the user must be sure to provide the electric field direction (via the `-ef` flag) in the correct convention. The DIALS and Precognition convention documentation can be found in the following:

[DIALS Convention ("Laboratory Frame" section)](https://dials.github.io/documentation/conventions.html)

[Precognition Convention (Section 3.2 - Goniometer Setting)](https://biocars.uchicago.edu/wp-content/uploads/2019/06/PrecognitionUserGuide_5-0.pdf)

Note that for data collected with a "vertical" electric field, the correct usage is likely `-ef 0 -1 0`; accordingly, this is the default value for `-ef` and can be omitted. If the `-ef` convention differs, it is essential that it be specified, e.g. as `-ef 1 0 0` with "horizontal" electric field.

## What does it do?
The user supplies a series of `.inp` files or an `.expt` file, which (among other things) describe the crystal orientation with an A matrix. The user also supplies the orientation of the electric field in the lab frame (via the `-ef` flag).  Given the crystal orientation and the EF direction, `regroup` walks through the possible crystal facets to see which will have a facet-normal parallel to the electric field. Then, `regroup` walks through the possible subgroups of the spacegroup to find those that will preserve the direction of the facet-normal (as an idealized electric field vector) for all of their subgroup symops.

## What do I do with the outputs?
`regroup` will print out something like the following:
```
           Facet       Angle                                    spacegroup n_symops
                        mean       std count                                       
0      (1, 1, 1)    2.231812  0.147201    46      P 1 (a+b,a-b,-c) (No. 1)        1
1      (1, 1, 0)   31.680535  0.145917    46      P 1 (a+b,a-b,-c) (No. 1)        1
2      (1, 0, 1)   34.806997  0.098994    46      P 1 (a+b,a-b,-c) (No. 1)        1
3      (0, 1, 1)   39.268555  0.117977    46      P 1 (a+b,a-b,-c) (No. 1)        1
4      (1, 0, 0)   50.767813  0.114569    46               C 2 1 1 (No. 5)        2
```
Each row corresponds to a facet of the crystal, ranked by how close the normal vector of that facet is to the electric field vector. In the above example, you can see that the electric field is very nearly normal to the (1, 1, 1) facet of the crystal. Across the 46 `.inp` files included, this facet is just ~2.2 degrees off from the EF vector with standard deviation <0.15 degrees.  
  
For each facet, `regroup` also tells you the new spacegroup caused by symmetry breaking along the EF vector. In the case of row 4 above, the new spacegroup would be C211 with no change in unit cell or indexing. However, rows 0-3 above require the change-of-basis `(a+b,a-b,-c)`. This can be accomplished via the following code (or similar):
```python
import reciprocalspaceship as rs
import gemmi

def change_spacegroup(filename, opstring, sg='P1'):
    """
    Function for reading in an unmerged, high-symmetry MTZ, changing basis and reindexing, and writing out the new version
    
    Parameters
    ----------
    filename : str
        mtz to be read in
    opstring : str
        String to feed to gemmi.Op constructor. Used as change of basis, and the inverse of which is used as the reindexing operation
    sg : str
        Optionally, provide a non-P1 reduced-symmetry spacegroup for the output.
    """
    
    mtz = rs.read_mtz(filename)
    mtz.remove_absences(inplace=True)
    op = gemmi.Op(opstring).inverse() # can confirm visually that this gives the intended rhombus

    # note that this works equivalently to gemmi's "unit cell reduction" algorithm
    mtz.cell = mtz.cell.changed_basis_forward(op, False) 
    mtz.spacegroup = sg
    mtz.apply_symop(op.inverse(), inplace=True) # reindexing via inverse operation

    mtz.write_mtz(filename.removesuffix('.mtz') + '_new_spacegroup.mtz')
    
    return

```
