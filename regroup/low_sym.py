"""
Convert high-symmetry MTZ files to a lower-symmetry spacegroup using a regroup
change-of-basis operator, if any.

This command saves the original high-symmetry HKLs as Hh, Kh, Lh, changes the
unit cell and space group, and if there is a change of basis, applies the 
basis-change operation and checks for correct basis change.

Example
-------
regroup.low_sym e020_{off,200ns}.mtz --op 'a-b,b-c,a+b+c' --ls_sg P1
"""



import argparse
from pathlib import Path

import reciprocalspaceship as rs
import gemmi
import numpy as np


def _transpose(op):
    op1 = gemmi.Op(op.triplet()) #  clone the rotation part of the operation
    op1.rot = np.array(op.rot).T
    return op1
    
def _new_op(rot, tran):
    op = gemmi.Op()
    op.rot = [list(map(int, row)) for row in rot]
    op.tran = list(map(int, tran))
    return op

def _has_any(s, letters):
    letters = letters + letters.upper()
    return any(ch in letters for ch in s)


def cctbx_cb_op_to_rs_op(s):
    """
    gemmi.Op into rs.apply_symop treats hkl, xyz, and abc equivalently, 
    but cctbx does not treat these equivalently.
    this function converts a cctbx change_of_basis_op string 
    to a gemmi.Op suitable for reciprocalspaceship.DataSet.apply_symop(). 
    if the cctbx change_of_basis_op is in xyz notation, then rs can use it directly.
    if the cctbx change_of_basis_op is in abc or hkl notation, then we transpose.

    Parameters
    ----------
    s : str
        cctbx cb_op string, e.g.
        "a-b,b-c,a+b+c"
        "a+b,a-b,-c"
        "-x+y+z,x+z+1/2,y-z"

    Returns
    -------
    gemmi.Op
        Op such that reciprocalspaceship apply_to_hkl performs
        H_new = H_old @ op.rot correctly.
        Translation is not set with the phase-shift sign expected by
        reciprocalspaceship, but this can be implemented as needed. 
    """
    # cctbx xyz notation is consistent with gemmi/rs. 
    if _has_any(s, "xyz"):
        return gemmi.Op(s)
    if _has_any(s, "hkl"):
        return gemmi.Op(s)

    # cctbx abc/hkl notation stores c_inv.T, while gemmi/rs store c.
    if _has_any(s, "abc"):
        c = _transpose(gemmi.Op(s)).inverse()
        return c 
    raise ValueError(f"Could not determine cctbx notation in {s!r}")


def _valid_cell_volume(mtz, op, _ratio, verbose=True):
    #cell1 = mtz.cell.changed_basis_forward(op,True) # with the right transformation, 
    cell1 = mtz.cell.changed_basis_forward(op.inverse(), True)
    old_new_smaller = mtz.cell.volume / cell1.volume
    if verbose:
        print("\n")
        print("old cell:", mtz.cell)
        print("new cell:" , cell1)
        print(f"The new cell has a volume {mtz.cell.volume/cell1.volume:0.2f} times smaller than the old cell.")
        
    return np.isclose(old_new_smaller, _ratio)

def _add_Hhs(mtz_in):
    mtz = mtz_in.copy()
    mtz["Hh"], mtz["Kh"], mtz["Lh"] = mtz.get_hkls().T
    mtz.Hh = mtz.Hh.astype("MTZInt")
    mtz.Kh = mtz.Kh.astype("MTZInt")
    mtz.Lh = mtz.Lh.astype("MTZInt")
    return mtz

def mtz_regroup_basis_change(mtz_path, op_from_regroup, lowsym, verbose=True):

    # we save copies of the high-symmetry HKLs. 
    mtz = rs.read_mtz(mtz_path)
    mtz = _add_Hhs(mtz)
    if lowsym is None: 
        return mtz
    op1 = cctbx_cb_op_to_rs_op(op_from_regroup)
    
    #we check that unit cell scaling is correct. 
    if isinstance(lowsym, str):
        lowsym = gemmi.SpaceGroup(lowsym)

    nl_cenops = len(lowsym.operations().cen_ops)

    nh_cenops = len(mtz.spacegroup.operations().cen_ops)
    _ratio = nh_cenops / nl_cenops

    if ~_valid_cell_volume(mtz, op1, _ratio):
        raise RuntimeError(
            f"Unit cell does not scale by the expected amount {_ratio}, (see documentation). please check your operation."
        )

    mtz = mtz.compute_dHKL()
    mtz["dHKL_old"] = mtz["dHKL"]

    #change spacegroup and unit cell 
    mtz.cell = mtz.cell.changed_basis_forward(op1.inverse(), True)
    mtz.spacegroup = lowsym
        
    #apply symop 
    if verbose: 
        print(f"op for careless: ", op1.inverse().triplet())
        #print(f"final operation: {op1}")
    try:
        mtz = mtz.apply_symop(op1)
    except:
        raise RuntimeError(
            f"Fractional Miller indices when reindexing HKLs, suggesting incorrect change-of-basis operation."
        )
    mtz = mtz.compute_dHKL()

    if not np.allclose(mtz["dHKL"], mtz["dHKL_old"], atol=1e-3):
        raise RuntimeError(
            f"old unit cell and new unit cell do not provide consistent dHKLs."
        )
    return mtz.drop(columns=["dHKL", "dHKL_old"])


def main():

    # CLI
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "hs_mtz",
        nargs="+",
        help="Input high-symmetry MTZ file",
    )
    parser.add_argument(
        "--op",
        default="x,y,z",
        help="Change-of-basis operator from regroup, e.g. 'a-b,b-c,a+b+c'",
    )
    parser.add_argument(
        "--ls_sg",
        default=None,
        help="Low-symmetry space group, e.g. 'P1'",
    )
    parser.add_argument(
        "-o", "--out",
        default=None,
        help="Output MTZ filename. Default: <input>_sg<ls_sg>.mtz",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress verbose cell-volume output",
    )

    args = parser.parse_args()

    if args.out is not None and len(args.hs_mtz) > 1:
        raise ValueError("--out can only be used with a single input MTZ.")

    for hs_mtz in args.hs_mtz:
        out = args.out
        if args.ls_sg is None:
            sg_num = rs.read_mtz(hs_mtz).spacegroup.number
        else:
            sg_num = gemmi.SpaceGroup(args.ls_sg).number
        if out is None:
            p = Path(hs_mtz)
            out = str(p.with_name(f"{p.stem}_sg{sg_num}{p.suffix}"))

        mtz = mtz_regroup_basis_change(
            hs_mtz,
            args.op,
            args.ls_sg,
            verbose=not args.quiet,
        )
        mtz.write_mtz(out)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()