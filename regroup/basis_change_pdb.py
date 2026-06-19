#!/usr/bin/env python3
"""Basis-change and symmetry-expand a PDB using a regroup subgroup JSON."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np
from regroup.low_sym import cctbx_cb_op_to_rs_op as cb_op_to_rs_op


CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
SOLVENT_RESNAMES = {"HOH", "WAT", "DOD", "SOL", "H2O"}
WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O"}
ION_RESNAMES = {
    "LI", "NA", "K", "RB", "CS", "BE", "MG", "CA", "SR", "BA",
    "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG", "AL", "GA",
    "IN", "TL", "PB", "AG", "AU", "PT", "PD", "F", "CL", "BR", "I", "IOD",
}
WATER_CHAIN = "W"
WATER_DUPLICATE_CUTOFF = 0.4


@dataclass
class AtomRecord:
    atom: gemmi.Atom
    frac0: gemmi.Fractional
    chain_id: str
    resname: str
    seqid: gemmi.SeqId
    het_flag: str
    segment: str
    source_residue_id: tuple[int, str, str, str, str]
    is_atom: bool

    @property
    def atom_name(self) -> str:
        return self.atom.name.strip()


@dataclass
class CopyBlock:
    rep_index: int
    rep_op: gemmi.Transform
    fracs: list[gemmi.Fractional]
    centroid: np.ndarray


@dataclass
class CopyCandidate:
    block: CopyBlock
    choice_index: int
    shift: np.ndarray
    point: np.ndarray


@dataclass
class ImportedSubgroup:
    label: str
    n_ops: int
    symbol_and_number: str
    basis_change_op: str
    ops_parent_setting: list[str]
    parent_ops_parent_setting: list[str]
    output_spacegroup_symbol: str
    ops_output_setting_compact: list[str]


OutputEntry = tuple[AtomRecord, gemmi.Position, str, int, gemmi.Transform]


def frac_arr(f: gemmi.Fractional | gemmi.Position | gemmi.Vec3) -> np.ndarray:
    return np.array([f.x, f.y, f.z], dtype=float)


def to_frac(a: np.ndarray) -> gemmi.Fractional:
    return gemmi.Fractional(float(a[0]), float(a[1]), float(a[2]))


def cell_matrix(cell: gemmi.UnitCell) -> np.ndarray:
    return np.array([frac_arr(cell.orthogonalize(to_frac(v))) for v in np.eye(3)], dtype=float).T


def op_to_transform(op: gemmi.Op) -> gemmi.Transform:
    s = op.float_seitz()
    return gemmi.Transform(gemmi.Mat33([row[:3] for row in s[:3]]), gemmi.Vec3(s[0][3], s[1][3], s[2][3]))


def translation_transform(shift: np.ndarray | tuple[float, float, float]) -> gemmi.Transform:
    return gemmi.Transform(gemmi.Mat33([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), gemmi.Vec3(*map(float, shift)))


def apply_transform(t: gemmi.Transform, f: gemmi.Fractional) -> gemmi.Fractional:
    v = t.apply(gemmi.Vec3(f.x, f.y, f.z))
    return gemmi.Fractional(v.x, v.y, v.z)


def transform_key(t: gemmi.Transform, ndigits: int = 7) -> tuple:
    m = np.round(np.array(t.mat.tolist(), dtype=float), ndigits)
    v = np.round(np.mod([t.vec.x, t.vec.y, t.vec.z], 1.0), ndigits)
    v[np.isclose(v, 1.0)] = 0.0
    return tuple(map(tuple, m)), tuple(v)


def unique_transforms(ops: list[gemmi.Transform]) -> list[gemmi.Transform]:
    out, seen = [], set()
    for op in ops:
        key = transform_key(op)
        if key not in seen:
            seen.add(key)
            out.append(op)
    return out


def transform_to_text(t: gemmi.Transform) -> str:
    m = np.array(t.mat.tolist(), dtype=float)
    v = np.array([t.vec.x, t.vec.y, t.vec.z], dtype=float)
    return f"mat={np.array2string(m, precision=5, suppress_small=True)}, vec={np.array2string(v, precision=5, suppress_small=True)}"


def transform_aniso_cartesian(aniso: gemmi.SMat33f, from_cell: gemmi.UnitCell, to_cell: gemmi.UnitCell, op: gemmi.Transform) -> gemmi.SMat33f:
    if not aniso.nonzero():
        return gemmi.SMat33f(0, 0, 0, 0, 0, 0)
    b_from, b_to = cell_matrix(from_cell), cell_matrix(to_cell)
    a = b_to @ np.array(op.mat.tolist(), dtype=float) @ np.linalg.inv(b_from)
    u = np.array([[aniso.u11, aniso.u12, aniso.u13], [aniso.u12, aniso.u22, aniso.u23], [aniso.u13, aniso.u23, aniso.u33]], dtype=float)
    u = 0.5 * (a @ u @ a.T + (a @ u @ a.T).T)
    return gemmi.SMat33f(float(u[0, 0]), float(u[1, 1]), float(u[2, 2]), float(u[0, 1]), float(u[0, 2]), float(u[1, 2]))


def compactness_key(points: np.ndarray, cell: gemmi.UnitCell) -> tuple[float, float, float, float]:
    xyz = points @ cell_matrix(cell).T
    center = np.mean(xyz, axis=0)
    rg = float(np.sum((xyz - center) ** 2))
    if len(xyz) <= 1:
        return rg, 0.0, 0.0, 0.0
    deltas = xyz[:, None, :] - xyz[None, :, :]
    extent = np.max(xyz, axis=0) - np.min(xyz, axis=0)
    return rg, float(np.sqrt(np.max(np.sum(deltas * deltas, axis=2)))), float(np.sum(extent)), float(np.prod(extent))


def rg_from_sums(sum_xyz: np.ndarray, sum_sq: float, n: int) -> float:
    return float(sum_sq - float(np.dot(sum_xyz, sum_xyz)) / float(n))


def operator_for_model_and_cell(op_string: str) -> tuple[gemmi.Op, gemmi.Op]:
    hkl_op = cb_op_to_rs_op(op_string)
    return hkl_op.inverse(), hkl_op


def op_key(op: gemmi.Op) -> tuple:
    op = op.wrap()
    return tuple(tuple(int(v) for v in row) for row in op.rot), tuple(int(t) % op.DEN for t in op.tran)


def unique_ops(ops: list[gemmi.Op]) -> list[gemmi.Op]:
    out, seen = [], set()
    for op in ops:
        key = op_key(op)
        if key not in seen:
            seen.add(key)
            out.append(op.wrap())
    return out


def load_imported_subgroup(path: str | Path) -> ImportedSubgroup:
    sel = json.loads(Path(path).read_text())["selected_subgroup"]
    return ImportedSubgroup(
        label=str(sel.get("label", Path(path).name)),
        n_ops=int(sel["n_ops"]),
        symbol_and_number=str(sel["symbol_and_number"]),
        basis_change_op=str(sel["basis_change_op"]),
        ops_parent_setting=[str(op).replace(" ", "") for op in sel["ops_parent_setting"]],
        parent_ops_parent_setting=[str(op).replace(" ", "") for op in sel["parent_ops_parent_setting"]],
        output_spacegroup_symbol=str(sel["output_spacegroup_symbol"]),
        ops_output_setting_compact=[str(op).replace(" ", "") for op in sel["ops_output_setting_compact"]],
    )


def conjugate_ops_to_output_basis(ops: list[gemmi.Op], basis_t: gemmi.Transform) -> list[gemmi.Transform]:
    basis_inv = basis_t.inverse()
    return unique_transforms([basis_t.combine(op_to_transform(op)).combine(basis_inv) for op in ops])


def output_origin_shift_to_declared_setting(basis_t: gemmi.Transform, child_parent_ops: list[gemmi.Op], output_ops: list[gemmi.Op]) -> np.ndarray:
    raw_child = conjugate_ops_to_output_basis(child_parent_ops, basis_t)
    target = {transform_key(op_to_transform(op)) for op in output_ops}
    if {transform_key(op) for op in raw_child}.issubset(target):
        return np.zeros(3, dtype=float)

    for q in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 96):
        for ijk in itertools.product(range(q), repeat=3):
            shift = np.array(ijk, dtype=float) / q
            t, tinv = translation_transform(shift), translation_transform(-shift)
            if {transform_key(t.combine(op).combine(tinv)) for op in raw_child}.issubset(target):
                return shift

    print("No output-origin correction found on the float grid; continuing with the raw basis.", file=sys.stderr)
    return np.zeros(3, dtype=float)


def parent_setting_cosets(imported: ImportedSubgroup) -> tuple[list[list[gemmi.Op]], list[gemmi.Op], list[gemmi.Op]]:
    h_parent = unique_ops([gemmi.Op(op) for op in imported.ops_parent_setting])
    g_parent = unique_ops([gemmi.Op(op) for op in imported.parent_ops_parent_setting])
    print(
        f"Using imported subgroup {imported.label}: n_ops={imported.n_ops}, {imported.symbol_and_number}, cb_op={imported.basis_change_op}",
        file=sys.stderr,
    )
    print(f"Using {len(g_parent)} parent operations from subgroup JSON for cosets", file=sys.stderr)

    cosets, seen = [], set()
    for g in g_parent:
        if op_key(g) in seen:
            continue
        coset = unique_ops([h.combine(g).wrap() for h in h_parent])
        cosets.append(coset)
        seen.update(op_key(product) for product in coset)
    return cosets, g_parent, h_parent


def child_parent_coset_rep_choices(imported: ImportedSubgroup, basis_t: gemmi.Transform) -> tuple[list[list[gemmi.Transform]], list[gemmi.Transform], list[gemmi.Transform]]:
    cosets_parent, g_parent, h_parent = parent_setting_cosets(imported)
    return (
        [conjugate_ops_to_output_basis(coset, basis_t) for coset in cosets_parent],
        conjugate_ops_to_output_basis(g_parent, basis_t),
        conjugate_ops_to_output_basis(h_parent, basis_t),
    )


def parse_atom_records(structure: gemmi.Structure, basis_t: gemmi.Transform, new_cell: gemmi.UnitCell) -> list[AtomRecord]:
    records: list[AtomRecord] = []
    for model_i, model in enumerate(structure):
        for chain in model:
            chain_id = chain.name or " "
            for residue in chain:
                resname = residue.name.strip()
                rid = (model_i, chain_id, str(residue.seqid), resname, residue.segment or "")
                het = residue.het_flag or "A"
                for atom in residue:
                    atom0 = atom.clone()
                    atom0.aniso = transform_aniso_cartesian(atom.aniso, structure.cell, new_cell, basis_t)
                    records.append(AtomRecord(atom0, apply_transform(basis_t, structure.cell.fractionalize(atom.pos)), chain_id, resname, residue.seqid, het, residue.segment or "", rid, het == "A"))
    return records


def build_copy_blocks(atom_records: list[AtomRecord], reps: list[gemmi.Transform]) -> list[CopyBlock]:
    non_solvent = [i for i, rec in enumerate(atom_records) if rec.resname not in SOLVENT_RESNAMES]
    anchors = set([i for i, rec in enumerate(atom_records) if rec.is_atom and rec.atom_name == "CA"] or non_solvent or range(len(atom_records)))
    blocks: list[CopyBlock] = []
    for rep_i, rep in enumerate(reps):
        fracs, centroid_fracs = [], []
        for atom_i, rec in enumerate(atom_records):
            frac = apply_transform(rep, rec.frac0)
            fracs.append(frac)
            if atom_i in anchors:
                centroid_fracs.append(frac_arr(frac))
        blocks.append(CopyBlock(rep_i, rep, fracs, np.mean(centroid_fracs, axis=0)))
    return blocks


def expand_choice_blocks_with_output_child_symmetry(choice_blocks: list[list[CopyBlock]], output_ops: list[gemmi.Op]) -> list[list[CopyBlock]]:
    output_ts = unique_transforms([op_to_transform(op.wrap()) for op in output_ops])
    expanded: list[list[CopyBlock]] = []
    for blocks in choice_blocks:
        new_blocks, seen = [], set()
        for block in blocks:
            for op_t in output_ts:
                rep = op_t.combine(block.rep_op)
                key = transform_key(rep)
                if key in seen:
                    continue
                seen.add(key)
                fracs = [apply_transform(op_t, frac) for frac in block.fracs]
                new_blocks.append(CopyBlock(block.rep_index, rep, fracs, frac_arr(apply_transform(op_t, to_frac(block.centroid)))))
        expanded.append(new_blocks)
    return expanded


def greedy_compact_block_translations(centroids: np.ndarray, cell: gemmi.UnitCell, search_radius: int) -> np.ndarray:
    candidates = np.array(list(itertools.product(range(-search_radius, search_radius + 1), repeat=3)), dtype=float)
    shifts = np.zeros((len(centroids), 3), dtype=float)
    placed = [centroids[0].copy()]
    b = cell_matrix(cell)
    for i in range(1, len(centroids)):
        best_shift, best_score = candidates[0], float("inf")
        for shift in candidates:
            c = centroids[i] + shift
            score = sum(float(np.dot((c - p) @ b.T, (c - p) @ b.T)) for p in placed)
            if score < best_score:
                best_shift, best_score = shift, score
        shifts[i] = best_shift
        placed.append(centroids[i] + best_shift)
    return shifts


def refine_compact_shifts(centroids: np.ndarray, cell: gemmi.UnitCell, start_shifts: np.ndarray, search_radius: int, max_sweeps: int = 20) -> np.ndarray:
    shifts = np.array(start_shifts, dtype=float) - start_shifts[0]
    n = len(centroids)
    if n <= 1:
        return shifts

    frac_to_cart = cell_matrix(cell).T
    centroid_cart = centroids @ frac_to_cart
    points_cart = (centroids + shifts) @ frac_to_cart
    sum_xyz = np.sum(points_cart, axis=0)
    sum_sq = float(np.sum(points_cart * points_cart))
    best_score = rg_from_sums(sum_xyz, sum_sq, n)
    local_moves = np.array(list(itertools.product(range(-search_radius, search_radius + 1), repeat=3)), dtype=float)
    pair_moves = np.array(list(itertools.product((-1, 0, 1), repeat=3)), dtype=float)

    for _ in range(max_sweeps):
        changed = False
        for i in range(1, n):
            old_shift, old_point = shifts[i].copy(), points_cart[i].copy()
            local = (best_score, old_shift, old_point, sum_xyz, sum_sq)
            for move in local_moves:
                trial_shift = old_shift + move
                if np.array_equal(trial_shift, old_shift):
                    continue
                trial_point = centroid_cart[i] + trial_shift @ frac_to_cart
                trial_sum = sum_xyz - old_point + trial_point
                trial_sum_sq = sum_sq - float(np.dot(old_point, old_point)) + float(np.dot(trial_point, trial_point))
                score = rg_from_sums(trial_sum, trial_sum_sq, n)
                if score < local[0] - 1e-7:
                    local = (score, trial_shift.copy(), trial_point.copy(), trial_sum, trial_sum_sq)
            if local[0] < best_score - 1e-7:
                best_score, shifts[i], points_cart[i], sum_xyz, sum_sq = local
                changed = True

        pair_changed = True
        while pair_changed:
            pair_changed = False
            for i in range(1, n):
                for j in range(i + 1, n):
                    old_i_shift, old_j_shift = shifts[i].copy(), shifts[j].copy()
                    old_i_point, old_j_point = points_cart[i].copy(), points_cart[j].copy()
                    local = (best_score, old_i_shift, old_j_shift, old_i_point, old_j_point, sum_xyz, sum_sq)
                    for move_i in pair_moves:
                        trial_i_shift = old_i_shift + move_i
                        trial_i_point = centroid_cart[i] + trial_i_shift @ frac_to_cart
                        for move_j in pair_moves:
                            if np.all(move_i == 0) and np.all(move_j == 0):
                                continue
                            trial_j_shift = old_j_shift + move_j
                            trial_j_point = centroid_cart[j] + trial_j_shift @ frac_to_cart
                            trial_sum = sum_xyz - old_i_point - old_j_point + trial_i_point + trial_j_point
                            trial_sum_sq = (
                                sum_sq
                                - float(np.dot(old_i_point, old_i_point))
                                - float(np.dot(old_j_point, old_j_point))
                                + float(np.dot(trial_i_point, trial_i_point))
                                + float(np.dot(trial_j_point, trial_j_point))
                            )
                            score = rg_from_sums(trial_sum, trial_sum_sq, n)
                            if score < local[0] - 1e-7:
                                local = (score, trial_i_shift.copy(), trial_j_shift.copy(), trial_i_point.copy(), trial_j_point.copy(), trial_sum, trial_sum_sq)
                    if local[0] < best_score - 1e-7:
                        best_score, shifts[i], shifts[j], points_cart[i], points_cart[j], sum_xyz, sum_sq = local
                        changed = pair_changed = True
        if not changed:
            break
    return shifts - shifts[0]


def compact_block_translations(blocks: list[CopyBlock], cell: gemmi.UnitCell, search_radius: int = 2) -> list[np.ndarray]:
    if len(blocks) <= 1:
        return [np.zeros(3, dtype=float) for _ in blocks]
    centroids = np.array([block.centroid for block in blocks], dtype=float)
    starts = [np.zeros_like(centroids), greedy_compact_block_translations(centroids, cell, search_radius)]
    starts.extend((-np.round(centroids - centroids[i])).astype(float) for i in range(len(centroids)))

    best_shifts, best_key = None, None
    for start in starts:
        refined = refine_compact_shifts(centroids, cell, start, search_radius)
        key = compactness_key(centroids + refined, cell)
        if best_key is None or key < best_key:
            best_shifts, best_key = refined.copy(), key
    return [shift - np.floor(np.mean(centroids + best_shifts, axis=0)) for shift in best_shifts]


def make_copy_candidates(choice_blocks: list[list[CopyBlock]], search_radius: int) -> list[list[CopyCandidate]]:
    translations = [np.array(t, dtype=float) for t in itertools.product(range(-search_radius, search_radius + 1), repeat=3)]
    return [[CopyCandidate(block, choice_i, shift.copy(), block.centroid + shift) for choice_i, block in enumerate(blocks) for shift in translations] for blocks in choice_blocks]


def nearest_candidate_index(candidates: list[CopyCandidate], target: np.ndarray, cell: gemmi.UnitCell) -> int:
    b = cell_matrix(cell)
    distances = [float(np.dot((candidate.point - target) @ b.T, (candidate.point - target) @ b.T)) for candidate in candidates]
    return int(np.argmin(distances))


def refine_candidate_selection(candidates: list[list[CopyCandidate]], cell: gemmi.UnitCell, start: list[int], max_sweeps: int = 30) -> list[int]:
    selection = list(start)
    n = len(candidates)
    frac_to_cart = cell_matrix(cell).T
    candidate_points_frac = [np.array([c.point for c in row], dtype=float) for row in candidates]
    candidate_points_cart = [pts @ frac_to_cart for pts in candidate_points_frac]
    points_frac = np.array([candidate_points_frac[i][selection[i]] for i in range(n)], dtype=float)
    points_cart = np.array([candidate_points_cart[i][selection[i]] for i in range(n)], dtype=float)
    sum_xyz, sum_sq = np.sum(points_cart, axis=0), float(np.sum(points_cart * points_cart))
    best_primary, best_key = rg_from_sums(sum_xyz, sum_sq, n), compactness_key(points_frac, cell)

    for _ in range(max_sweeps):
        changed = False
        for i in range(n):
            old, old_frac, old_cart = selection[i], points_frac[i].copy(), points_cart[i].copy()
            local = (best_primary, best_key, old, old_frac, old_cart, sum_xyz, sum_sq)
            for j in range(len(candidates[i])):
                if j == old:
                    continue
                new_frac, new_cart = candidate_points_frac[i][j], candidate_points_cart[i][j]
                trial_sum = sum_xyz - old_cart + new_cart
                trial_sum_sq = sum_sq - float(np.dot(old_cart, old_cart)) + float(np.dot(new_cart, new_cart))
                primary = rg_from_sums(trial_sum, trial_sum_sq, n)
                key = None
                take = primary < local[0] - 1e-7
                if not take and abs(primary - local[0]) <= 1e-7:
                    trial_points = points_frac.copy()
                    trial_points[i] = new_frac
                    key = compactness_key(trial_points, cell)
                    take = key < local[1]
                if take:
                    local = (primary, key or (primary, 0.0, 0.0, 0.0), j, new_frac.copy(), new_cart.copy(), trial_sum, trial_sum_sq)
            if local[2] != old:
                best_primary, _key, selection[i], points_frac[i], points_cart[i], sum_xyz, sum_sq = local
                best_key = compactness_key(points_frac, cell)
                changed = True
        if not changed:
            break
    return selection


def refine_candidate_selection_pairs(candidates: list[list[CopyCandidate]], cell: gemmi.UnitCell, start: list[int], max_pair_candidates: int = 32, max_sweeps: int = 8) -> list[int]:
    selection = list(start)
    n = len(candidates)
    frac_to_cart = cell_matrix(cell).T
    candidate_points_frac = [np.array([c.point for c in row], dtype=float) for row in candidates]
    candidate_points_cart = [pts @ frac_to_cart for pts in candidate_points_frac]
    points_frac = np.array([candidate_points_frac[i][selection[i]] for i in range(n)], dtype=float)
    points_cart = np.array([candidate_points_cart[i][selection[i]] for i in range(n)], dtype=float)
    sum_xyz, sum_sq = np.sum(points_cart, axis=0), float(np.sum(points_cart * points_cart))
    best_primary, best_key = rg_from_sums(sum_xyz, sum_sq, n), compactness_key(points_frac, cell)

    def ranked_for(i: int) -> list[int]:
        if len(candidates[i]) <= max_pair_candidates:
            return list(range(len(candidates[i])))
        old_cart = points_cart[i].copy()
        scored = []
        for j, new_cart in enumerate(candidate_points_cart[i]):
            trial_sum = sum_xyz - old_cart + new_cart
            trial_sum_sq = sum_sq - float(np.dot(old_cart, old_cart)) + float(np.dot(new_cart, new_cart))
            scored.append((rg_from_sums(trial_sum, trial_sum_sq, n), j))
        keep = {selection[i], *(j for _score, j in sorted(scored)[: max_pair_candidates - 1])}
        return sorted(keep)

    for _ in range(max_sweeps):
        ranked = [ranked_for(i) for i in range(n)]
        best_move = None
        for i in range(n):
            for j in range(i + 1, n):
                old_i, old_j = selection[i], selection[j]
                old_i_cart, old_j_cart = points_cart[i], points_cart[j]
                for ci in ranked[i]:
                    new_i_cart = candidate_points_cart[i][ci]
                    for cj in ranked[j]:
                        if ci == old_i and cj == old_j:
                            continue
                        new_j_cart = candidate_points_cart[j][cj]
                        trial_sum = sum_xyz - old_i_cart - old_j_cart + new_i_cart + new_j_cart
                        trial_sum_sq = (
                            sum_sq
                            - float(np.dot(old_i_cart, old_i_cart))
                            - float(np.dot(old_j_cart, old_j_cart))
                            + float(np.dot(new_i_cart, new_i_cart))
                            + float(np.dot(new_j_cart, new_j_cart))
                        )
                        primary = rg_from_sums(trial_sum, trial_sum_sq, n)
                        key = None
                        take = primary < best_primary - 1e-7
                        if not take and abs(primary - best_primary) <= 1e-7:
                            trial_points = points_frac.copy()
                            trial_points[i] = candidate_points_frac[i][ci]
                            trial_points[j] = candidate_points_frac[j][cj]
                            key = compactness_key(trial_points, cell)
                            take = key < best_key
                        if take and (best_move is None or (primary, key or (primary, 0.0, 0.0, 0.0)) < (best_move[0], best_move[1])):
                            best_move = (primary, key or (primary, 0.0, 0.0, 0.0), i, j, ci, cj, trial_sum, trial_sum_sq)
        if best_move is None:
            break
        best_primary, _key, i, j, ci, cj, sum_xyz, sum_sq = best_move
        selection[i], selection[j] = ci, cj
        points_frac[i], points_frac[j] = candidate_points_frac[i][ci], candidate_points_frac[j][cj]
        points_cart[i], points_cart[j] = candidate_points_cart[i][ci], candidate_points_cart[j][cj]
        best_key = compactness_key(points_frac, cell)
    return selection


def select_compact_copy_choices(choice_blocks: list[list[CopyBlock]], cell: gemmi.UnitCell, search_radius: int = 2) -> tuple[list[CopyBlock], list[np.ndarray], list[int]]:
    if all(len(blocks) == 1 for blocks in choice_blocks):
        blocks = [blocks[0] for blocks in choice_blocks]
        return blocks, compact_block_translations(blocks, cell, search_radius), [0 for _ in blocks]

    candidates = make_copy_candidates(choice_blocks, min(search_radius, 1))
    starts = [[nearest_candidate_index(candidates[i], choice_blocks[i][0].centroid, cell) for i in range(len(candidates))]]
    for anchor_i, anchor_candidates in enumerate(candidates):
        best_anchor: dict[int, tuple[int, CopyCandidate]] = {}
        for idx, candidate in enumerate(anchor_candidates):
            old = best_anchor.get(candidate.choice_index)
            if old is None or tuple(np.abs(candidate.shift)) < tuple(np.abs(old[1].shift)):
                best_anchor[candidate.choice_index] = (idx, candidate)
        for anchor_idx, anchor in best_anchor.values():
            starts.append([anchor_idx if i == anchor_i else nearest_candidate_index(candidates[i], anchor.point, cell) for i in range(len(candidates))])

    best_selection, best_key = starts[0], None
    for start in starts:
        selection = refine_candidate_selection(candidates, cell, start)
        key = compactness_key(np.array([candidates[i][selection[i]].point for i in range(len(candidates))]), cell)
        if best_key is None or key < best_key:
            best_selection, best_key = selection, key
    best_selection = refine_candidate_selection_pairs(candidates, cell, best_selection)

    selected_blocks, shifts, choice_indices = [], [], []
    placed = np.array([candidates[i][best_selection[i]].point for i in range(len(candidates))])
    global_shift = -np.floor(np.mean(placed, axis=0))
    for i, sel in enumerate(best_selection):
        candidate = candidates[i][sel]
        selected_blocks.append(candidate.block)
        shifts.append(candidate.shift + global_shift)
        choice_indices.append(candidate.choice_index)

    translation_shifts = compact_block_translations(selected_blocks, cell, search_radius)
    if compactness_key(np.array([b.centroid + s for b, s in zip(selected_blocks, translation_shifts)]), cell) < compactness_key(np.array([b.centroid + s for b, s in zip(selected_blocks, shifts)]), cell):
        shifts = translation_shifts
    return selected_blocks, shifts, choice_indices


def assign_chain_ids(atom_records: list[AtomRecord], n_copies: int) -> dict[tuple[int, str], str]:
    old_chain_ids = []
    for rec in atom_records:
        if rec.resname not in WATER_RESNAMES and rec.chain_id not in old_chain_ids:
            old_chain_ids.append(rec.chain_id)
    available = (chain for chain in CHAIN_IDS if chain != WATER_CHAIN)
    return {(copy_i, old): next(available) for copy_i in range(n_copies) for old in old_chain_ids}


def build_output_entries(atom_records: list[AtomRecord], blocks: list[CopyBlock], shifts: list[np.ndarray], new_cell: gemmi.UnitCell) -> list[OutputEntry]:
    chain_map = assign_chain_ids(atom_records, len(blocks))
    entries: list[OutputEntry] = []
    for copy_i, block in enumerate(blocks):
        for rec, frac in zip(atom_records, block.fracs):
            chain_id = WATER_CHAIN if rec.resname in WATER_RESNAMES else chain_map[(copy_i, rec.chain_id)]
            pos = new_cell.orthogonalize(to_frac(frac_arr(frac) + shifts[copy_i]))
            entries.append((rec, pos, chain_id, copy_i, block.rep_op))
    return entries


def remove_close_waters_and_ions(entries: list[OutputEntry], cell: gemmi.UnitCell, output_ops: list[gemmi.Op]) -> tuple[list[OutputEntry], int]:
    cutoff2 = WATER_DUPLICATE_CUTOFF * WATER_DUPLICATE_CUTOFF
    ops = unique_ops([gemmi.Op("x,y,z"), *[op.wrap() for op in output_ops]])
    kept_sites: list[tuple[str, np.ndarray]] = []
    drop: set[tuple[int, tuple[int, str, str, str, str]]] = set()

    for rec, pos, _chain_id, copy_i, _op in entries:
        if rec.resname in WATER_RESNAMES:
            if rec.atom_name not in {"O", "OH2", "OW"}:
                continue
            site_key = "WATER"
        elif rec.resname.upper() in ION_RESNAMES:
            site_key = f"ION:{rec.resname.upper()}"
        else:
            continue

        frac_arr0 = np.mod(frac_arr(cell.fractionalize(pos)), 1.0)
        duplicate = False
        for kept_key, kept_frac in kept_sites:
            if kept_key != site_key:
                continue
            for op in ops:
                sym = np.array(op.apply_to_xyz([float(x) for x in frac_arr0]), dtype=float)
                delta = sym - kept_frac - np.round(sym - kept_frac)
                d = cell.orthogonalize(to_frac(delta))
                if d.x * d.x + d.y * d.y + d.z * d.z <= cutoff2:
                    duplicate = True
                    break
            if duplicate:
                break
        if duplicate:
            drop.add((copy_i, rec.source_residue_id))
        else:
            kept_sites.append((site_key, frac_arr0))

    return [e for e in entries if (e[3], e[0].source_residue_id) not in drop], len(drop)


def water_numbers(entries: list[OutputEntry]) -> dict[tuple[int, tuple[int, str, str, str, str]], int]:
    numbers: dict[tuple[int, tuple[int, str, str, str, str]], int] = {}
    for rec, _pos, chain_id, copy_i, _op in entries:
        if chain_id == WATER_CHAIN and rec.resname in WATER_RESNAMES:
            numbers.setdefault((copy_i, rec.source_residue_id), len(numbers) + 1)
    return numbers


def write_pdb_with_gemmi(entries: list[OutputEntry], outfile: str | Path, new_cell: gemmi.UnitCell, spacegroup: str) -> int:
    structure = gemmi.Structure()
    structure.cell = new_cell
    structure.spacegroup_hm = spacegroup
    model = gemmi.Model("1")
    chains: dict[str, gemmi.Chain] = {}
    residues: dict[tuple[str, object], gemmi.Residue] = {}
    residue_order: list[tuple[str, object]] = []
    nums = water_numbers(entries)

    for serial, (rec, pos, chain_id, copy_i, op) in enumerate(entries, start=1):
        chains.setdefault(chain_id, gemmi.Chain(chain_id))
        water_num = nums.get((copy_i, rec.source_residue_id))
        key = (chain_id, water_num if water_num is not None else (copy_i, rec.source_residue_id))
        if key not in residues:
            residue = gemmi.Residue()
            residue.name = rec.resname
            residue.seqid = gemmi.SeqId(water_num, " ") if water_num is not None else rec.seqid
            residue.het_flag = rec.het_flag
            residue.segment = rec.segment
            residues[key] = residue
            residue_order.append(key)
        atom = rec.atom.clone()
        atom.pos = pos
        atom.serial = serial
        atom.aniso = transform_aniso_cartesian(rec.atom.aniso, new_cell, new_cell, op)
        residues[key].add_atom(atom)

    for key in residue_order:
        chains[key[0]].add_residue(residues[key])
    for chain in chains.values():
        model.add_chain(chain)
    structure.add_model(model)
    if hasattr(structure, "write_minimal_pdb"):
        structure.write_minimal_pdb(str(outfile))
    else:
        structure.write_pdb(str(outfile))
    return len(nums)


def write_expanded_low_sym_asu(infile: str | Path, outfile: str | Path, child_subgroup_json: str | Path) -> tuple:
    input_structure = gemmi.read_structure(str(infile), merge_chain_parts=False)
    old_cell = input_structure.cell
    imported = load_imported_subgroup(child_subgroup_json)
    basis_op_raw, hkl_op = operator_for_model_and_cell(imported.basis_change_op)
    new_cell = old_cell.changed_basis_forward(basis_op_raw, True)
    output_child_ops = unique_ops([gemmi.Op(op).wrap() for op in imported.ops_output_setting_compact])

    basis_t = op_to_transform(basis_op_raw)
    origin_shift = output_origin_shift_to_declared_setting(basis_t, unique_ops([gemmi.Op(op) for op in imported.ops_parent_setting]), output_child_ops)
    if np.any(origin_shift):
        basis_t = translation_transform(origin_shift).combine(basis_t)
        print(f"Applied output-origin correction {tuple(origin_shift)}", file=sys.stderr)

    rep_choices, g_output, h_output = child_parent_coset_rep_choices(imported, basis_t)
    atom_records = parse_atom_records(input_structure, basis_t, new_cell)
    choice_blocks = expand_choice_blocks_with_output_child_symmetry([build_copy_blocks(atom_records, choices) for choices in rep_choices], output_child_ops)
    blocks, shifts, _choice_indices = select_compact_copy_choices(choice_blocks, new_cell, search_radius=2)
    output_entries, removed = remove_close_waters_and_ions(build_output_entries(atom_records, blocks, shifts, new_cell), new_cell, output_child_ops)
    renumbered = write_pdb_with_gemmi(output_entries, outfile, new_cell, imported.output_spacegroup_symbol)

    if removed:
        print(f"Removed {removed} symmetry-equivalent waters/ions within {WATER_DUPLICATE_CUTOFF} A", file=sys.stderr)
    if renumbered:
        print(f"Renumbered {renumbered} water residues on chain {WATER_CHAIN}", file=sys.stderr)
    return old_cell, new_cell, basis_op_raw, hkl_op, [block.rep_op for block in blocks], g_output, h_output, imported


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Change basis and expand the ASU from a regroup subgroup JSON.")
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("--child-subgroup-json", required=True)
    args = parser.parse_args(argv)

    old_cell, new_cell, basis_op, hkl_op, reps, g_ops, h_ops, imported = write_expanded_low_sym_asu(args.infile, args.outfile, args.child_subgroup_json)
    print(f"input operator: {imported.basis_change_op}")
    print(f"HKL/low_sym operator: {hkl_op.triplet()}")
    print(f"model/cell operator used for unit cell: {basis_op.triplet()}")
    if reps:
        print(f"first output coset rep: {transform_to_text(reps[0])}")
    print(f"old cell: {old_cell}")
    print(f"new cell: {new_cell}")
    print(f"output SG: {imported.output_spacegroup_symbol}")
    print(f"|G|={len(g_ops)}, |H|={len(h_ops)}, coset reps={len(reps)}")
    print(f"wrote: {args.outfile}")


if __name__ == "__main__":
    main()
