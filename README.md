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

`regroup.low_sym` converts high-symmetry MTZ files to a lower-symmetry spacegroup using a regroup change-of-basis operator, if any. This command saves the original high-symmetry HKLs as Hh, Kh, Lh, changes the unit cell and space group, and if there is a change of basis, applies the basis-change operation and checks for correct basis change.


## Conventions
Since `regroup` can process both DIALS experiment files and Precognition input files, it is important to make use of the correct lab frame convention. The program will automatically infer the lab frame convention given the type of file input, but the user must be sure to provide the electric field direction (via the `-ef` flag) in the correct convention. The DIALS and Precognition convention documentation can be found in the following:

[DIALS Convention ("Laboratory Frame" section)](https://dials.github.io/documentation/conventions.html)

[Precognition Convention (Section 3.2 - Goniometer Setting)](https://biocars.uchicago.edu/wp-content/uploads/2019/06/PrecognitionUserGuide_5-0.pdf)

Note that for data collected with a positive electric field pointing to the lab floor, the correct usage is likely `-ef 0 1 0`; accordingly, this is the default value for `-ef` and can be omitted. If the `-ef` convention differs, it is essential that it be specified, e.g. as `-ef 1 0 0` with "horizontal" electric field, pointed to the center of the storage ring, in the case of APS. Regroup does not yet support cases where the electric field moves with the crystal.

## What does it do?
The user supplies a series of `.inp` files or an `.expt` file, which (among other things) describe the crystal orientation with an A matrix. The user also supplies the orientation of the electric field in the lab frame (via the `-ef` flag).  Given the crystal orientation and the EF direction, `regroup` walks through the possible crystal facets to see which will have a facet-normal parallel to the electric field. Then, `regroup` walks through the possible subgroups of the spacegroup to find those that will preserve the direction of the facet-normal (as an idealized electric field vector) for all of their subgroup symops. The following discussion will assume familiarity with the contents of the accompanying manuscript. 

Here is an example `regroup` function call. 

`regroup ../INDEXING/precog_index/e047d_off_*.mccd.inp --hmax 1 -sg 213 -ef 0 1 0 --filename figures/regroup.log --fsa --export-subgroup-json ./scripts/subgroup_PDZ3.json --ideal_vector`

The first argument is a set of Precognition `.inp` files. The high-symmetry space group is specified with the `-sg` flag. More information on flags can be found with `regroup -h` or the below. 

## What should I do with the outputs?

```
Best-match basis-change op: y-z,-x+y+z,x+z+1/2
Best-match facet: (1, -1, -1)
Transformed best-match facet: [ 0.  0. -3.]
```

Regroup first provides a best-match change-of-basis operation for use later (`y-z,-x+y+z,x+z+1/2`), as well as a best match facet normal (here, `(1,-1-1)`) . The facet normal is also provided as transformed into the low-symmetry space group (here, `(0, 0,-3)`). In other words, in the low-symmetry space group, the best facet normal is `(0,0,-3)`. This can be used later as well. 

```
           Facet       Angle                                             spacegroup n_symops
                        mean       std count                                                
0    (1, -1, -1)   15.166536  0.234400    46  R 3 :H (y-z,-x+y+z,x+z+1/2) (No. 146)        3
1    (0, -1, -1)   27.256358  0.209043    46   C 1 2 1 (z-3/8,x+y-1/4,-x+y) (No. 5)        2
2     (1, -1, 0)   33.393002  0.256032    46    C 1 2 1 (x-y,x+y-1/4,z-1/8) (No. 5)        2
3     (0, -1, 0)   39.984649  0.275708    46          P 41 (b+1/4,c+1/2,a) (No. 76)        4

```
Each row corresponds to a facet of the crystal, ranked by how close the normal vector of that facet is to the electric field vector. The first row is considered the "best". In the above example, you can see that the electric field is close to the crystal's (1, -1, -1) facet normal. Across the 46 `.inp` files included, this facet is ~15 degrees off from the EF vector with standard deviation ~0.23 degrees. `n_symops` indicates that 3 symops (including the identity) are preserved.

Regroup also tells you the field-symmetry alignment (FSA) of each symmetry operation. This is how much each symmetry operation "reorients" the electric field. If the FSA is 1, then the symmetry operation does not affect the electric field. If the FSA is -1, then the symmetry operation reverses the electric field. This can be used for identifying broken and preserved symmetry operations. More information about the FSA can be found in the accompanying manuscript. 

```
       fsa   angle_deg   broken?   opnum  symop                    
-------------------------------------------------------------------
   1.00000       0.000     False       0  x,y,z                    
   0.33333      70.529      True       1  -y+1/4,x+3/4,z+1/4       
  -0.33333     109.471      True       2  -x+1/2,-y,z+1/2          
   0.33333      70.529      True       3  y+1/4,-x+1/4,z+3/4       
  -0.33333     109.471      True       4  x+1/2,-y+1/2,-z          
   0.33333      70.529      True       5  -y+3/4,-x+3/4,-z+3/4     
  -0.33333     109.471      True       6  -x,y+1/2,-z+1/2          
  -1.00000     180.000      True       7  y+3/4,x+1/4,-z+1/4       
  -0.33333     109.471      True       8  z,x,y                    
  -1.00000     180.000      True       9  z+1/4,-y+1/4,x+3/4       
  -0.33333     109.471      True      10  z+1/2,-x+1/2,-y          
   0.33333      70.529      True      11  z+3/4,y+1/4,-x+1/4       
  -0.33333     109.471      True      12  -z,x+1/2,-y+1/2          
   0.33333      70.529      True      13  -z+3/4,-y+3/4,-x+3/4     
   1.00000       0.000     False      14  -z+1/2,-x,y+1/2          
   0.33333      70.529      True      15  -z+1/4,y+3/4,x+1/4       
   0.33333      70.529      True      16  -x+1/4,z+3/4,y+1/4       
  -0.33333     109.471      True      17  y,z,x                    
   0.33333      70.529      True      18  x+3/4,z+1/4,-y+1/4       
   1.00000       0.000     False      19  -y,z+1/2,-x+1/2          
  -1.00000     180.000      True      20  -x+3/4,-z+3/4,-y+3/4     
  -0.33333     109.471      True      21  y+1/2,-z+1/2,-x          
   0.33333      70.529      True      22  x+1/4,-z+1/4,y+3/4       
  -0.33333     109.471      True      23  -y+1/2,-z,x+1/2
```

Operations 1, 14 and 19 are not broken when using the best `(1,-1,-1)` facet normal. We can turn on this option with `--ideal_vector`. Else, operations 14 and 19 are broken, but their FSA are close to 1: 

```
       fsa   angle_deg   broken?   opnum  symop                    
-------------------------------------------------------------------
   1.00000       0.000     False       0  x,y,z                    
...    
   0.89736      26.187      True      14  -z+1/2,-x,y+1/2          
...
   0.89736      26.187      True      19  -y,z+1/2,-x+1/2          
...
```

We can use this later to pick the ideal CCsym.


For each facet, `regroup` also tells you the new spacegroup caused by symmetry breaking along the EF vector. In the case of row 0 above, the new spacegroup would be R3:H with change in unit cell and indexing using the "y-z,-x+y+z,x+z+1/2" operation. 
For example, we can change the basis of any unmerged MTZ files with observed Miller indices. We can do this via the following function call: 

```
regroup.low_sym ./mtzs/cdef_hs_e047_{off,200ns}.mtz
regroup.low_sym ./mtzs/cdef_hs_e047_{off,200ns}.mtz --op="y-z,-x+y+z,x+z+1/2" --ls_sg="R3:H"
```

The first call of `regroup.low_sym` is for adding `Hh, Kh, Lh` flags to `./mtzs/cdef_hs_e047_{off,200ns}.mtz` without changing space group, outputting the following MTZs: 
```
  ./mtzs/cdef_hs_e047_off_sg213.mtz
  ./mtzs/cdef_hs_e047_200ns_sg213.mtz
```

The second call of `regroup.low_sym` changes the `HKLobs` of `./mtzs/cdef_hs_e047_{off,200ns}.mtz`, outputs 
```

  ./mtzs/cdef_hs_e047_off_sg146.mtz 
  ./mtzs/cdef_hs_e047_200ns_sg146.mtz 
```
and the following logging:

```
old cell: <gemmi.UnitCell(90.5, 90.5, 90.5, 90, 90, 90)>
new cell: <gemmi.UnitCell(127.986, 127.986, 156.751, 90, 90, 120)>
The new cell has a volume 0.33 times smaller than the old cell.
op for careless:  x/3-y/3+2/3*z-1/3,2/3*x+y/3+z/3-1/6,-x/3+y/3+z/3-1/6
```

This is what we expect for the primitive-to-centered unit cell change of basis. We also get an op for careless, so that is what we will use. 

To run careless, we use, for example: 

```
expt=e047
in=./mtzs
out=./careless_runs/careless-
careless poly \
  <other keys>
  --double-wilson-parents=None,0,0 \
  --double-wilson-r=0,0.998,0.998 \
   --double-wilson-reindexing-ops="x,y,z;x/3-y/3+2/3*z-1/3,2/3*x+y/3+z/3-1/6,-x/3+y/3+z/3-1/6;x/3-y/3+2/3*z-1/3,2/3*x+y/3+z/3-1/6,-x/3+y/3+z/3-1/6" \
  "Hh,Kh,Lh,X,Y,file_id,Wavelength,dHKL,BATCH" \
  $in/cdef_hs_${expt}_off_sg213.mtz \
  $in/cdef_hs_${expt}_off_sg146.mtz \
  $in/cdef_hs_${expt}_200ns_sg146.mtz \
  $out/merged_${expt}
```

Full careless scripts can be found in the accompanying Zenodo deposition. We can then compute CCsyms from the sg146 Careless outputs, as in the accompanying Zenodo deposition. 

To change the basis of a low-symmetry model:

```
regroup.basis_change_pdb ./refinement/rbr/pdz3_cript_ds4-2_refine_85_chain-merge.pdb ./refinement/rbr/ls_cb.pdb --child-subgroup-json=./scripts/subgroup_PDZ3.json
```

where we have obtained the `./scripts/subgroup_PDZ3.json` from using the `--export-subgroup-json ./scripts/subgroup_PDZ3.json` in the original regroup function call. `./refinement/rbr/pdz3_cript_ds4-2_refine_85_chain-merge.pdb` is the input pdb file with the high-symmetry asymmetric unit, and `./refinement/rbr/ls_cb.pdb` is the desired output pdb destination. We recommend, however, first building your own pdb file using molecular replacement to get a sense of the symmetry. 

Once careless finishes runnning, you are ready to refine your `./refinement/rbr/ls_cb.pdb` against the Careless output (with or without extrapolation). 


