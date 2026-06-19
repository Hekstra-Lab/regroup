#!/usr/bin/env python3
"""Minimal subgroup JSON export helpers for regroup.py.

The JSON written here keeps entries needed for basis_change_pdb.py
to perform compact expansion of the high-symmetry ASU, placed
into the right origin.

"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import gemmi
from cctbx import sgtbx


KEEP_FIELDS = (
    "label",
    "index",
    "n_ops",
    "n_ops_parent_setting_compact",
    "symbol_and_number",
    "basis_change_op",
    "parent_sg",
    "parent_spacegroup",
    "ops_parent_setting",
    "ops_parent_setting_compact",
    "parent_ops_parent_setting",
    "parent_ops_parent_setting_compact",
    "n_parent_ops",
    "n_parent_ops_compact",
    "centering_ops_parent_setting",
    "output_spacegroup_symbol",
    "ops_output_setting_compact",
    "n_ops_output_setting_compact",
    "centering_ops_output_setting",
)


def extract_basis_change_op(symbol_and_number: Any) -> str:
    """Extract the cctbx basis-change op from symbol_and_number()."""
    text = str(symbol_and_number)
    if " (No." in text:
        text = text.split(" (No.", 1)[0]
    if "(" in text and ")" in text:
        return text.rsplit("(", 1)[1].split(")", 1)[0].strip()
    return "x,y,z"


def rt_mx_to_xyz(op: Any) -> str:
    """output a compact xyz string for a cctbx sgtbx.rt_mx."""
    for name in ("as_xyz", "as_xyz_string"):
        fn = getattr(op, name, None)
        if callable(fn):
            try:
                return str(fn()).replace(" ", "")
            except TypeError:
                pass
    return str(op).replace(" ", "")


def extract_output_spacegroup_symbol(symbol_and_number: Any) -> str:
    """Return the conventional output SG symbol from cctbx symbol_and_number()."""
    text = str(symbol_and_number).strip()
    if " (No." in text:
        text = text.split(" (No.", 1)[0].strip()
    if "(" in text and ")" in text:
        text = text.rsplit("(", 1)[0].strip()
    return " ".join(text.split())


def gemmi_op(op: Any) -> gemmi.Op:
    """Return a Gemmi operator from either an xyz string or op-like object."""
    if isinstance(op, gemmi.Op):
        return op
    if isinstance(op, str):
        return gemmi.Op(op.replace(" ", ""))
    return gemmi.Op(rt_mx_to_xyz(op))


def op_key(op: gemmi.Op) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]:
    """canonical key for a Gemmi Op, with translations modulo one unit cell."""
    wrapped = op.wrap()
    rot = tuple(tuple(int(v) for v in row) for row in wrapped.rot)
    tran = tuple(int(t) % wrapped.DEN for t in wrapped.tran)
    return rot, tran


def gemmi_op_strings(ops: list[Any]) -> list[str]:
    """Return unique canonical Gemmi op triplets, preserving order."""
    out: list[str] = []
    seen: set[tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]] = set()
    for op in ops:
        gop = gemmi_op(op).wrap()
        key = op_key(gop)
        if key not in seen:
            seen.add(key)
            out.append(gop.triplet().replace(" ", ""))
    return out


def tran_part(t: int) -> str:
    """Format a Gemmi translation integer as an xyz-triplet fraction suffix."""
    from fractions import Fraction

    t = int(t) % gemmi.Op.DEN
    if t == 0:
        return ""
    frac = Fraction(t, gemmi.Op.DEN)
    return f"+{frac.numerator}/{frac.denominator}"


def clean_spacegroup_symbol(spacegroup: Any) -> str:
    """Return a Gemmi-readable Hermann-Mauguin symbol from a cctbx or  SG object."""
    try:
        symbol = str(sgtbx.space_group_info(spacegroup).symbol_and_number())
    except Exception:
        symbol = str(spacegroup)
    return extract_output_spacegroup_symbol(symbol)


def centering_ops_for_spacegroup(spacegroup: Any) -> list[gemmi.Op]:
    """Return centering translation operators for a space group setting."""
    symbol = clean_spacegroup_symbol(spacegroup)
    if not symbol:
        return [gemmi.Op("x,y,z")]
    ops = gemmi.SpaceGroup(symbol).operations()
    centering = []
    for tran in ops.cen_ops:
        triplet = ",".join(axis + tran_part(t) for axis, t in zip("xyz", tran))
        centering.append(gemmi.Op(triplet).wrap())
    return [gemmi.Op(op) for op in gemmi_op_strings(centering)]


def expand_ops_by_centering(ops: list[Any], centering_ops: list[gemmi.Op]) -> list[str]:
    """Expand compact representatives by centering translations.

    cctbx's smx() list omits centering translations.  For coset arithmetic in
    the parent setting we need the full crystallographic operation group, so
    each compact representative R*x+v is written once for every parent-lattice
    centering translation c, i.e. R*x+v+c.
    """
    expanded: list[gemmi.Op] = []
    for op in ops:
        base = gemmi_op(op).wrap()
        for cen in centering_ops:
            expanded.append(cen.combine(base).wrap())
    return gemmi_op_strings(expanded)


def output_setting_compact_fields(symbol_and_number: Any) -> dict[str, Any]:
    """Return output-setting symops used only for compact placement."""
    output_sg = extract_output_spacegroup_symbol(symbol_and_number)
    fields: dict[str, Any] = {"output_spacegroup_symbol": output_sg}
    if not output_sg:
        return fields

    ops = gemmi.SpaceGroup(output_sg).operations()
    fields["ops_output_setting_compact"] = gemmi_op_strings(list(ops))
    fields["n_ops_output_setting_compact"] = len(fields["ops_output_setting_compact"])

    centering = []
    for tran in ops.cen_ops:
        triplet = ",".join(axis + tran_part(t) for axis, t in zip("xyz", tran))
        centering.append(gemmi.Op(triplet).triplet().replace(" ", ""))
    fields["centering_ops_output_setting"] = gemmi_op_strings(centering)
    return fields


def subgroup_record(index: int, subgroup: Any, subgroup_info: Any | None = None, facet: Any = None) -> dict[str, Any]:
    """Serialize one cctbx subgroup in the parent setting."""
    if subgroup_info is None:
        subgroup_info = sgtbx.space_group_info(group=subgroup)
    symbol = str(subgroup_info.symbol_and_number())
    compact_ops = [rt_mx_to_xyz(op) for op in subgroup.smx()]
    record: dict[str, Any] = {
        "label": f"cctbx_subgroup_{index}",
        "index": int(index),
        # This is updated to the full count in write_subgroup_json(), where the
        # parent centering translations are known.
        "n_ops": int(subgroup.n_smx()),
        "n_ops_parent_setting_compact": int(subgroup.n_smx()),
        "symbol_and_number": symbol,
        "basis_change_op": extract_basis_change_op(symbol),
        "ops_parent_setting": compact_ops,
        "ops_parent_setting_compact": compact_ops,
    }
    if facet is not None:
        record["facet"] = list(map(int, facet))
    record.update(output_setting_compact_fields(symbol))
    return record


def write_subgroup_json(filename: str | Path, parent_sg: Any, selected_record: dict[str, Any]) -> None:
    """Write only the subgroup fields consumed by basis_change_pdb.py."""
    selected = dict(selected_record)
    selected.setdefault("parent_sg", str(parent_sg))
    selected.setdefault("parent_spacegroup", str(parent_sg))

    parent_group = sgtbx.space_group_info(parent_sg).group()
    parent_compact_ops = [rt_mx_to_xyz(op) for op in parent_group.smx()]
    child_compact_ops = list(
        selected.get("ops_parent_setting_compact", selected.get("ops_parent_setting", []))
    )

    parent_centering_ops = centering_ops_for_spacegroup(parent_sg)
    selected["centering_ops_parent_setting"] = gemmi_op_strings(parent_centering_ops)

    # basis_change_pdb.py can form H\G cosets directly without knowing about
    # lattice-centering conventions.
    selected["ops_parent_setting_compact"] = gemmi_op_strings(child_compact_ops)
    selected["n_ops_parent_setting_compact"] = len(selected["ops_parent_setting_compact"])
    selected["ops_parent_setting"] = expand_ops_by_centering(child_compact_ops, parent_centering_ops)
    selected["n_ops"] = len(selected["ops_parent_setting"])

    selected["parent_ops_parent_setting_compact"] = gemmi_op_strings(parent_compact_ops)
    selected["n_parent_ops_compact"] = len(selected["parent_ops_parent_setting_compact"])
    selected["parent_ops_parent_setting"] = expand_ops_by_centering(parent_compact_ops, parent_centering_ops)
    selected["n_parent_ops"] = len(selected["parent_ops_parent_setting"])

    selected = {key: selected[key] for key in KEEP_FIELDS if key in selected}
    payload = {
        "parent_sg": str(parent_sg),
        "selected_subgroup": selected,
    }
    with Path(filename).open("w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
