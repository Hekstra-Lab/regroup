# regroup

Determine new space groups for analyzing pump-probe crystallography experiments

## Installation  

This program relies on `sgtbx` for the hierarchical grouping of different crystallographic space groups. Currently, it does not seem that `cctbx` can be easily installed with `pip`, so this dependency must be installed separately. The following snippet should install the `regroup` command-line program to your current environment:

```shell
git clone https://github.com/Hekstra-Lab/regroup.git
cd regroup
conda install -c conda-forge cctbx
python setup.py install
```

## Features  

`regroup` requires knowledge of the experimental geometry of the crystal in the lab reference frame in order to determine the new space group based on the orientation of the crystal relative to the "pump" perturbation. Since much of our group's work is conducted at the BioCARS Laue beamline (APS 14-ID-B), this program currently only supports Precognition geometry files (`.inp` format).

If anyone is interested in support for additional file formats, please reach out by filing an [issue](https://github.com/Hekstra-Lab/regroup/issues).
