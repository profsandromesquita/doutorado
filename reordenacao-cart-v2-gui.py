# -*- coding: utf-8 -*-
"""
================================================================================
PIPELINE UNIFICADO DE REORDENAÇÃO DE ESTRUTURA CAR-T - VERSÃO 2.0 COM GUI
================================================================================

Versão com interface gráfica (tkinter) para processamento de trajetórias
de dinâmica molecular. Pode ser empacotado como executável standalone.

AUTOR: Pipeline consolidado automaticamente
DATA: 2025
================================================================================
"""

import math
import sys
import os
import threading
import queue
from typing import List, Dict, Set, Tuple, Optional, Any
from collections import defaultdict

# Interface gráfica
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext


# ==============================================================================
# SEÇÃO 1: FUNÇÕES GEOMÉTRICAS FUNDAMENTAIS
# ==============================================================================

def dist(a: Dict, b: Dict) -> float:
    return math.sqrt(
        (a['x'] - b['x'])**2 +
        (a['y'] - b['y'])**2 +
        (a['z'] - b['z'])**2
    )


def dist_tuple(p1: Tuple[float, float, float], p2: Tuple[float, float, float]) -> float:
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def angle(a: Dict, b: Dict, c: Dict) -> float:
    v1 = (a['x'] - b['x'], a['y'] - b['y'], a['z'] - b['z'])
    v2 = (c['x'] - b['x'], c['y'] - b['y'], c['z'] - b['z'])
    n1 = math.sqrt(sum(v*v for v in v1))
    n2 = math.sqrt(sum(v*v for v in v2))
    if n1 < 1e-6 or n2 < 1e-6:
        return 0.0
    dot = sum(v1[i]*v2[i] for i in range(3))
    cosang = max(-1.0, min(1.0, dot/(n1*n2)))
    return math.degrees(math.acos(cosang))


def angle_tuple(p1: Tuple[float, float, float], p2: Tuple[float, float, float],
                p3: Tuple[float, float, float]) -> float:
    v1 = (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])
    v2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
    n1 = math.sqrt(sum(v*v for v in v1))
    n2 = math.sqrt(sum(v*v for v in v2))
    if n1 < 1e-6 or n2 < 1e-6:
        return 0.0
    dot = sum(v1[i]*v2[i] for i in range(3))
    cosang = max(-1.0, min(1.0, dot/(n1*n2)))
    return math.degrees(math.acos(cosang))


def same_coords(a: Dict, b: Dict, tol: float = 1e-3) -> bool:
    return (
        abs(a['x'] - b['x']) < tol and
        abs(a['y'] - b['y']) < tol and
        abs(a['z'] - b['z']) < tol
    )


def swap_coords(a: Dict, b: Dict) -> None:
    for coord in ('x', 'y', 'z'):
        a[coord], b[coord] = b[coord], a[coord]


# ==============================================================================
# SEÇÃO 2: FUNÇÕES DE PARSING E ESCRITA DE PDB - MULTI-FRAME
# ==============================================================================

def split_trajectory_into_frames(content: str) -> List[str]:
    frames = []
    current_frame_lines = []
    lines = content.splitlines()

    for line in lines:
        current_frame_lines.append(line)
        if line.strip() == "ENDMDL":
            frames.append("\n".join(current_frame_lines))
            current_frame_lines = []

    if current_frame_lines:
        remaining = "\n".join(current_frame_lines).strip()
        if remaining:
            frames.append(remaining)

    return frames


def parse_frame_content(content: str) -> Tuple[List[Dict], List[str]]:
    lines = content.splitlines()
    atoms = []

    for i, line in enumerate(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                rec = {
                    'line_idx': i,
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'altloc': line[16] if len(line) > 16 else '',
                    'resname': line[17:20].strip(),
                    'chain': line[21].strip() if len(line) > 21 else '',
                    'resseq': int(line[22:26]),
                    'icode': line[26] if len(line) > 26 else '',
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'raw_line': line
                }
                atoms.append(rec)
            except (ValueError, IndexError):
                continue

    return atoms, lines


def update_pdb_lines(lines: List[str], atoms: List[Dict]) -> List[str]:
    for at in atoms:
        i = at['line_idx']
        line = lines[i]
        if len(line) < 54:
            line = line.ljust(54)
        new_line = (
            line[:30] +
            f"{at['x']:8.3f}{at['y']:8.3f}{at['z']:8.3f}" +
            line[54:]
        )
        lines[i] = new_line
    return lines


# ==============================================================================
# SEÇÃO 3: FUNÇÕES DE BUSCA E INDEXAÇÃO DE ÁTOMOS
# ==============================================================================

def find_residue_atoms(atoms: List[Dict], resname: str, resseq: int, chain: Optional[str] = None) -> List[int]:
    idxs = []
    for i, a in enumerate(atoms):
        if a['resname'] == resname and a['resseq'] == resseq:
            if chain is None or a['chain'] == chain:
                idxs.append(i)
    return idxs


def find_atom_by_name(atoms: List[Dict], resseq: int, atom_name: str, chain: Optional[str] = None) -> Optional[int]:
    for i, a in enumerate(atoms):
        if a['resseq'] == resseq and a['name'] == atom_name:
            if chain is None or a['chain'] == chain:
                return i
    return None


def build_name_index(atoms: List[Dict], residue_idxs: List[int]) -> Dict[str, int]:
    name2idx = {}
    for i in residue_idxs:
        nm = atoms[i]['name']
        if nm not in name2idx:
            name2idx[nm] = i
    return name2idx


def build_serial_index(atoms: List[Dict]) -> Dict[int, Dict]:
    return {a['serial']: a for a in atoms}


def get_scope_serials(atoms: List[Dict], res_start: int, res_end: int, chain: str) -> List[int]:
    return [
        a['serial'] for a in atoms
        if a['chain'] == chain and res_start <= a['resseq'] <= res_end
    ]


def get_scope_indices(atoms: List[Dict], res_start: int, res_end: int, chain: str) -> List[int]:
    return [
        i for i, a in enumerate(atoms)
        if a['chain'] == chain and res_start <= a['resseq'] <= res_end
    ]


# ==============================================================================
# SEÇÃO 4: FUNÇÕES DE IDENTIFICAÇÃO DE BACKBONE
# ==============================================================================

def identify_backbone_serial_range(atoms: List[Dict], start_serial: int, end_serial: int) -> List[int]:
    bb_names = {'N', 'CA', 'C'}
    return sorted([
        a['serial'] for a in atoms
        if start_serial <= a['serial'] <= end_serial and a['name'] in bb_names
    ])


def get_backbone_residues(backbone_ids: List[int], atom_by_serial: Dict[int, Dict]) -> List[Dict]:
    seen_res = set()
    residues = []
    for serial in backbone_ids:
        a = atom_by_serial[serial]
        key = (a['chain'], a['resseq'])
        if key not in seen_res:
            seen_res.add(key)
            residues.append({'chain': a['chain'], 'resseq': a['resseq'], 'resname': a['resname']})
    return residues


# ==============================================================================
# SEÇÃO 5: FUNÇÕES DE SELEÇÃO DINÂMICA DE VIZINHOS
# ==============================================================================

def escolher_vizinhos_dinamico(atoms: List[Dict], locked: Set[int], idx_ref: int,
                                target_count: int, dmin_init: float, dmax_init: float,
                                delta: float = 0.005, max_iter: int = 200,
                                label: str = "H?", scope_indices: List[int] = None) -> Tuple[List[int], float, float]:
    ref = atoms[idx_ref]
    dmin, dmax = dmin_init, dmax_init

    indices_to_search = scope_indices if scope_indices is not None else range(len(atoms))

    for iteration in range(max_iter):
        candidates = []
        for i in indices_to_search:
            if i in locked:
                continue
            d = dist(ref, atoms[i])
            if dmin <= d <= dmax:
                candidates.append((i, d))

        if len(candidates) >= target_count:
            candidates.sort(key=lambda x: x[1])
            chosen = [c[0] for c in candidates[:target_count]]
            return chosen, dmin, dmax

        dmin = max(0.0, dmin - delta)
        dmax = dmax + delta

    candidates = []
    for i in indices_to_search:
        if i in locked:
            continue
        d = dist(ref, atoms[i])
        if dmin <= d <= dmax:
            candidates.append((i, d))

    if candidates:
        candidates.sort(key=lambda x: x[1])
        chosen = [c[0] for c in candidates[:target_count]]
        return chosen, dmin, dmax

    raise RuntimeError(f"[{label}] Não foi possível encontrar {target_count} vizinho(s) para átomo ref={idx_ref}")


def atribuir_coord_alvos(atoms: List[Dict], locked: Set[int], idx_ref: int, alvo_idxs: List[int],
                          dmin_init: float, dmax_init: float, target_count: int, label: str,
                          scope_indices: List[int] = None) -> None:
    if not alvo_idxs:
        return
    try:
        vizinhos, dmin_final, dmax_final = escolher_vizinhos_dinamico(
            atoms, locked, idx_ref, target_count, dmin_init, dmax_init, label=label,
            scope_indices=scope_indices
        )
        for i, alvo_idx in enumerate(alvo_idxs):
            if i < len(vizinhos):
                viz_idx = vizinhos[i]
                if not same_coords(atoms[alvo_idx], atoms[viz_idx]):
                    swap_coords(atoms[alvo_idx], atoms[viz_idx])
                locked.add(alvo_idx)
    except RuntimeError:
        pass


def escolher_pesado_com_angulo(atoms: List[Dict], locked: Set[int], idx_ref: int, idx_ref2: int,
                                dmin: float, dmax: float, ang1_min: float, ang1_max: float,
                                ang2_min: float, ang2_max: float, alvo_angulo: float,
                                label: str = "", expand_distance: bool = False,
                                delta: float = 0.005, max_iter: int = 200,
                                idx_angle_a: int = None, idx_angle_b: int = None,
                                scope_indices: List[int] = None) -> Tuple[int, float, float]:
    ref = atoms[idx_ref]
    ref2 = atoms[idx_ref2]

    indices_to_search = scope_indices if scope_indices is not None else range(len(atoms))

    current_dmin, current_dmax = dmin, dmax

    for iteration in range(max_iter if expand_distance else 1):
        candidates = []
        for i in indices_to_search:
            if i in locked:
                continue
            a = atoms[i]
            if a['name'].startswith('H'):
                continue
            d = dist(ref, a)
            if current_dmin <= d <= current_dmax:
                ang = angle(ref2, ref, a)
                if ang1_min <= ang <= ang1_max:
                    if idx_angle_a is not None and idx_angle_b is not None:
                        ang2 = angle(atoms[idx_angle_a], atoms[idx_angle_b], a)
                        if ang2_min <= ang2 <= ang2_max:
                            candidates.append((i, d, ang))
                    else:
                        candidates.append((i, d, ang))

        if candidates:
            candidates.sort(key=lambda x: abs(x[2] - alvo_angulo))
            best = candidates[0]
            return best[0], best[1], best[2]

        if expand_distance:
            current_dmin = max(0.0, current_dmin - delta)
            current_dmax = current_dmax + delta

    raise RuntimeError(f"[{label}] Nenhum candidato pesado encontrado para ref={idx_ref}")


# ==============================================================================
# SEÇÃO 6: BUSCA DFS PARA CAMINHO DO BACKBONE
# ==============================================================================

def find_all_backbone_paths_dfs(atoms: List[Dict], backbone_ids: List[int],
                                 coords: Dict[int, Tuple[float, float, float]],
                                 verbose: bool = False) -> Tuple[List[int], List[Dict]]:
    if len(backbone_ids) < 2:
        return backbone_ids, []

    atom_by_serial = {a['serial']: a for a in atoms}

    start_atom = atom_by_serial[backbone_ids[0]]
    end_atom = atom_by_serial[backbone_ids[-1]]
    res_start = start_atom['resseq']
    res_end = end_atom['resseq']
    chain_alvo = start_atom['chain']
    all_serials = get_scope_serials(atoms, res_start, res_end, chain_alvo)

    def get_coord(serial):
        return coords.get(serial, (atom_by_serial[serial]['x'],
                                   atom_by_serial[serial]['y'],
                                   atom_by_serial[serial]['z']))

    def get_neighbors(serial, used, dmin, dmax):
        p1 = get_coord(serial)
        neighbors = []
        for s in all_serials:
            if s in used:
                continue
            p2 = get_coord(s)
            d = dist_tuple(p1, p2)
            if dmin <= d <= dmax:
                neighbors.append((s, d))
        neighbors.sort(key=lambda x: x[1])
        return [n[0] for n in neighbors]

    target_length = len(backbone_ids)
    start_serial = backbone_ids[0]

    DMIN_INITIAL = 1.20
    DMAX_INITIAL = 1.60
    DELTA = 0.02
    MIN_DMIN = 1.10
    MAX_EXPANSIONS = 15

    expansoes = []

    def dfs_search(dmin, dmax):
        stack = [(start_serial, [start_serial], {start_serial})]

        while stack:
            current, path, used = stack.pop()

            if len(path) == target_length:
                return path

            neighbors = get_neighbors(current, used, dmin, dmax)

            for neighbor in neighbors:
                new_used = used | {neighbor}
                new_path = path + [neighbor]
                stack.append((neighbor, new_path, new_used))

        return None

    dmin, dmax = DMIN_INITIAL, DMAX_INITIAL
    result = dfs_search(dmin, dmax)

    if result:
        return result, expansoes

    for expansion in range(MAX_EXPANSIONS):
        new_dmin = max(MIN_DMIN, dmin - DELTA)
        new_dmax = dmax + DELTA

        if new_dmin == dmin and new_dmax == dmax:
            break

        dmin, dmax = new_dmin, new_dmax

        exp_info = {
            'expansion': expansion + 1,
            'dmin': dmin,
            'dmax': dmax
        }
        expansoes.append(exp_info)

        result = dfs_search(dmin, dmax)
        if result:
            return result, expansoes

    return backbone_ids, expansoes


def formatar_relatorio_expansoes(expansoes: List[Dict], frame_num: int = None) -> str:
    if not expansoes:
        return ""

    lines = []
    header = f"Frame {frame_num}" if frame_num else "Backbone"
    lines.append(f"\n{'='*60}")
    lines.append(f"EXPANSÕES DE JANELA - {header}")
    lines.append(f"{'='*60}")

    for exp in expansoes:
        lines.append(f"  Expansão {exp['expansion']}: dmin={exp['dmin']:.3f}, dmax={exp['dmax']:.3f}")

    return "\n".join(lines)


# ==============================================================================
# SEÇÃO 7: FASE 1 - RECONSTRUÇÃO DO BACKBONE
# ==============================================================================

def fase1_backbone(atoms: List[Dict], coords: Dict[int, Tuple[float, float, float]],
                   start_serial: int, end_serial: int,
                   verbose: bool = False) -> Tuple[Dict[int, Tuple[float, float, float]], List[int], List[Dict]]:
    backbone_ids = identify_backbone_serial_range(atoms, start_serial, end_serial)

    if not backbone_ids:
        return coords, [], []

    ordered_path, expansoes = find_all_backbone_paths_dfs(atoms, backbone_ids, coords, verbose=verbose)

    new_coords = coords.copy()
    for i, target_serial in enumerate(backbone_ids):
        source_serial = ordered_path[i]
        if target_serial != source_serial:
            new_coords[target_serial], new_coords[source_serial] = \
                new_coords[source_serial], new_coords[target_serial]

    return new_coords, backbone_ids, expansoes


# ==============================================================================
# SEÇÃO 8: FASE 2 - ATRIBUIÇÃO DE HN, HA, O
# ==============================================================================

def fase2_hn_ha_o(atoms: List[Dict], coords: Dict[int, Tuple[float, float, float]],
                  backbone_ids: List[int], verbose: bool = False) -> Dict[int, Tuple[float, float, float]]:
    atom_by_serial = {a['serial']: a for a in atoms}
    residues = get_backbone_residues(backbone_ids, atom_by_serial)

    start_atom = atom_by_serial[backbone_ids[0]]
    end_atom = atom_by_serial[backbone_ids[-1]]
    res_start = start_atom['resseq']
    res_end = end_atom['resseq']
    chain_alvo = start_atom['chain']
    all_serials = get_scope_serials(atoms, res_start, res_end, chain_alvo)

    locked = set(backbone_ids)
    new_coords = coords.copy()

    def get_coord(serial):
        return new_coords.get(serial, (atom_by_serial[serial]['x'],
                                       atom_by_serial[serial]['y'],
                                       atom_by_serial[serial]['z']))

    for res in residues:
        chain = res['chain']
        resseq = res['resseq']

        N_serial = None
        CA_serial = None
        C_serial = None

        for serial in backbone_ids:
            a = atom_by_serial[serial]
            if a['chain'] == chain and a['resseq'] == resseq:
                if a['name'] == 'N':
                    N_serial = serial
                elif a['name'] == 'CA':
                    CA_serial = serial
                elif a['name'] == 'C':
                    C_serial = serial

        # HN
        HN_serial = None
        for a in atoms:
            if a['chain'] == chain and a['resseq'] == resseq and a['name'] == 'HN':
                HN_serial = a['serial']
                break

        if HN_serial and N_serial:
            N_coord = get_coord(N_serial)
            best_cand = None
            best_dist = float('inf')

            for s in all_serials:
                if s in locked:
                    continue
                p = get_coord(s)
                d = dist_tuple(N_coord, p)
                if 0.90 <= d <= 1.15 and d < best_dist:
                    best_dist = d
                    best_cand = s

            if best_cand and best_cand != HN_serial:
                new_coords[HN_serial], new_coords[best_cand] = new_coords[best_cand], new_coords[HN_serial]
            if HN_serial:
                locked.add(HN_serial)

        # HA
        HA_serial = None
        for a in atoms:
            if a['chain'] == chain and a['resseq'] == resseq and a['name'] == 'HA':
                HA_serial = a['serial']
                break

        if HA_serial and CA_serial:
            CA_coord = get_coord(CA_serial)
            best_cand = None
            best_dist = float('inf')

            for s in all_serials:
                if s in locked:
                    continue
                p = get_coord(s)
                d = dist_tuple(CA_coord, p)
                if 0.95 <= d <= 1.15 and d < best_dist:
                    best_dist = d
                    best_cand = s

            if best_cand and best_cand != HA_serial:
                new_coords[HA_serial], new_coords[best_cand] = new_coords[best_cand], new_coords[HA_serial]
            if HA_serial:
                locked.add(HA_serial)

        # O
        O_serial = None
        for a in atoms:
            if a['chain'] == chain and a['resseq'] == resseq and a['name'] == 'O':
                O_serial = a['serial']
                break

        if O_serial and C_serial:
            C_coord = get_coord(C_serial)
            best_cand = None
            best_dist = float('inf')

            for s in all_serials:
                if s in locked:
                    continue
                a = atom_by_serial.get(s)
                if a and a['name'].startswith('H'):
                    continue
                p = get_coord(s)
                d = dist_tuple(C_coord, p)
                if 1.15 <= d <= 1.35 and d < best_dist:
                    best_dist = d
                    best_cand = s

            if best_cand and best_cand != O_serial:
                new_coords[O_serial], new_coords[best_cand] = new_coords[best_cand], new_coords[O_serial]
            if O_serial:
                locked.add(O_serial)

    return new_coords


# ==============================================================================
# SEÇÃO 9: FASE 3 - ATRIBUIÇÃO DE CB
# ==============================================================================

def fase3_ca_cb(atoms: List[Dict], coords: Dict[int, Tuple[float, float, float]],
                backbone_ids: List[int], locked_serials: Set[int],
                verbose: bool = False) -> Tuple[Dict[int, Tuple[float, float, float]], Set[int]]:
    atom_by_serial = {a['serial']: a for a in atoms}
    residues = get_backbone_residues(backbone_ids, atom_by_serial)

    start_atom = atom_by_serial[backbone_ids[0]]
    end_atom = atom_by_serial[backbone_ids[-1]]
    res_start = start_atom['resseq']
    res_end = end_atom['resseq']
    chain_alvo = start_atom['chain']
    all_serials = get_scope_serials(atoms, res_start, res_end, chain_alvo)

    new_coords = coords.copy()
    locked = locked_serials.copy()

    def get_coord(serial):
        return new_coords.get(serial, (atom_by_serial[serial]['x'],
                                       atom_by_serial[serial]['y'],
                                       atom_by_serial[serial]['z']))

    for res in residues:
        chain = res['chain']
        resseq = res['resseq']
        resname = res['resname']

        if resname == 'GLY':
            continue

        CA_serial = None
        N_serial = None
        C_serial = None

        for serial in backbone_ids:
            a = atom_by_serial[serial]
            if a['chain'] == chain and a['resseq'] == resseq:
                if a['name'] == 'CA':
                    CA_serial = serial
                elif a['name'] == 'N':
                    N_serial = serial
                elif a['name'] == 'C':
                    C_serial = serial

        CB_serial = None
        for a in atoms:
            if a['chain'] == chain and a['resseq'] == resseq and a['name'] == 'CB':
                CB_serial = a['serial']
                break

        if not all([CA_serial, N_serial, C_serial, CB_serial]):
            continue

        CA_coord = get_coord(CA_serial)
        N_coord = get_coord(N_serial)
        C_coord = get_coord(C_serial)

        candidates = []
        for s in all_serials:
            if s in locked:
                continue
            a = atom_by_serial.get(s)
            if a and a['name'].startswith('H'):
                continue

            p = get_coord(s)
            d = dist_tuple(CA_coord, p)

            if 1.45 <= d <= 1.65:
                ang_N = angle_tuple(N_coord, CA_coord, p)
                ang_C = angle_tuple(C_coord, CA_coord, p)

                if 105.0 <= ang_N <= 115.0 and 105.0 <= ang_C <= 115.0:
                    candidates.append((s, d, ang_N, ang_C))

        chosen_serial = None
        if candidates:
            if len(candidates) > 1:
                working = [(s, d, tN, tC) for s, d, tN, tC in candidates]
                def score(item):
                    _, _, tN, tC = item
                    return abs(tN - 110.0) + abs(tC - 111.0)
                working.sort(key=score)
                chosen_serial = working[0][0]
            else:
                chosen_serial = candidates[0][0]

            if chosen_serial != CB_serial:
                tmp = new_coords[CB_serial]
                new_coords[CB_serial] = new_coords[chosen_serial]
                new_coords[chosen_serial] = tmp

            locked.add(CB_serial)

    return new_coords, locked


# ==============================================================================
# SEÇÃO 10: FUNÇÕES AUXILIARES PARA CADEIAS LATERAIS
# ==============================================================================

def aplicar_coords_para_atoms(atoms: List[Dict], coords: Dict[int, Tuple[float, float, float]]) -> None:
    serial_to_idx = {a['serial']: i for i, a in enumerate(atoms)}
    for serial, (x, y, z) in coords.items():
        if serial in serial_to_idx:
            idx = serial_to_idx[serial]
            atoms[idx]['x'] = x
            atoms[idx]['y'] = y
            atoms[idx]['z'] = z


def extrair_coords_de_atoms(atoms: List[Dict]) -> Dict[int, Tuple[float, float, float]]:
    return {a['serial']: (a['x'], a['y'], a['z']) for a in atoms}


def bloquear_residuo(atoms: List[Dict], locked: Set[int], resname: str, resseq: int, chain: str) -> None:
    for i, a in enumerate(atoms):
        if a['resname'] == resname and a['resseq'] == resseq and a['chain'] == chain:
            locked.add(i)


def inicializar_locked_base(atoms: List[Dict]) -> Set[int]:
    locked = set()
    base_names = {"N", "CA", "C", "HN", "HA", "O", "CB"}
    for i, a in enumerate(atoms):
        if a['name'] in base_names:
            locked.add(i)
    return locked


# ==============================================================================
# SEÇÃO 11: PROCESSAMENTO DE CADEIAS LATERAIS POR RESÍDUO
# ==============================================================================

def processar_ile(atoms: List[Dict], locked: Set[int], resseq: int = 576, chain: str = "B",
                  verbose: bool = True, scope_indices: List[int] = None) -> None:
    res_idxs = find_residue_atoms(atoms, "ILE", resseq, chain)
    if not res_idxs:
        return
    name2idx = build_name_index(atoms, res_idxs)
    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_CG1 = name2idx.get("CG1")
    idx_CG2 = name2idx.get("CG2")
    idx_CD = name2idx.get("CD")
    if idx_CD is None:
        idx_CD = name2idx.get("CD1")

    if idx_CG2 is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.40, dmax=1.65,
                ang1_min=101.0, ang1_max=119.0, ang2_min=104.0, ang2_max=116.0, alvo_angulo=110.0, label="CG2",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CG2], atoms[cand_idx]):
                swap_coords(atoms[idx_CG2], atoms[cand_idx])
            locked.add(idx_CG2)
        except RuntimeError:
            pass

    alvo_hg2 = [name2idx.get(n) for n in ("1HG2", "2HG2", "3HG2") if name2idx.get(n) is not None]
    if alvo_hg2 and idx_CG2 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CG2, alvo_hg2, dmin_init=0.990, dmax_init=1.150,
                            target_count=len(alvo_hg2), label="HG2", scope_indices=scope_indices)

    if idx_CG1 is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.40, dmax=1.65,
                ang1_min=101.0, ang1_max=119.0, ang2_min=104.0, ang2_max=116.0, alvo_angulo=110.0, label="CG1",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CG1], atoms[cand_idx]):
                swap_coords(atoms[idx_CG1], atoms[cand_idx])
            locked.add(idx_CG1)
        except RuntimeError:
            pass

    idx_1HG1 = name2idx.get("1HG1")
    idx_2HG1 = name2idx.get("2HG1")
    alvo_cg1 = [x for x in [idx_1HG1, idx_2HG1, idx_CD] if x is not None]
    if alvo_cg1 and idx_CG1 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CG1, alvo_cg1, dmin_init=0.990, dmax_init=1.150,
                            target_count=len(alvo_cg1), label="CG1_viz", scope_indices=scope_indices)

    alvo_hd = [name2idx.get(n) for n in ("HD1", "HD2", "HD3") if name2idx.get(n) is not None]
    if alvo_hd and idx_CD is not None:
        atribuir_coord_alvos(atoms, locked, idx_CD, alvo_hd, dmin_init=0.990, dmax_init=1.150,
                            target_count=len(alvo_hd), label="HD", scope_indices=scope_indices)


def processar_thr(atoms: List[Dict], locked: Set[int], resseq: int = 577, chain: str = "B",
                  verbose: bool = True, scope_indices: List[int] = None) -> None:
    res_idxs = find_residue_atoms(atoms, "THR", resseq, chain)
    if not res_idxs:
        return
    name2idx = build_name_index(atoms, res_idxs)
    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_OG1 = name2idx.get("OG1")
    idx_HG1 = name2idx.get("HG1")
    idx_HB = name2idx.get("HB")
    idx_CG2 = name2idx.get("CG2")

    if idx_OG1 is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.30, dmax=1.50,
                ang1_min=101.0, ang1_max=119.0, ang2_min=104.0, ang2_max=116.0, alvo_angulo=110.0, label="OG1",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_OG1], atoms[cand_idx]):
                swap_coords(atoms[idx_OG1], atoms[cand_idx])
            locked.add(idx_OG1)
        except RuntimeError:
            pass

    if idx_HG1 is not None and idx_OG1 is not None:
        atribuir_coord_alvos(atoms, locked, idx_OG1, [idx_HG1], dmin_init=0.90, dmax_init=1.10,
                            target_count=1, label="HG1", scope_indices=scope_indices)

    if idx_HB is not None and idx_CB is not None:
        atribuir_coord_alvos(atoms, locked, idx_CB, [idx_HB], dmin_init=1.00, dmax_init=1.20,
                            target_count=1, label="HB", scope_indices=scope_indices)

    if idx_CG2 is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0, ang2_min=104.0, ang2_max=116.0, alvo_angulo=110.0, label="CG2",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CG2], atoms[cand_idx]):
                swap_coords(atoms[idx_CG2], atoms[cand_idx])
            locked.add(idx_CG2)
        except RuntimeError:
            pass

    alvo_hg2 = [name2idx.get(n) for n in ("1HG2", "2HG2", "3HG2") if name2idx.get(n) is not None]
    if alvo_hg2 and idx_CG2 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CG2, alvo_hg2, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hg2), label="HG2", scope_indices=scope_indices)


def processar_leu(atoms: List[Dict], locked: Set[int], resseq: int = 578, chain: str = "B",
                  verbose: bool = True, scope_indices: List[int] = None) -> None:
    res_idxs = find_residue_atoms(atoms, "LEU", resseq, chain)
    if not res_idxs:
        return
    name2idx = build_name_index(atoms, res_idxs)
    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_CG = name2idx.get("CG")
    idx_HG = name2idx.get("HG")
    idx_CD1 = name2idx.get("CD1")
    idx_CD2 = name2idx.get("CD2")
    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")

    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(atoms, locked, idx_CB, alvo_hb, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hb), label="HB", scope_indices=scope_indices)

    if idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0, ang2_min=104.0, ang2_max=116.0, alvo_angulo=110.0, label="CG",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CG], atoms[cand_idx]):
                swap_coords(atoms[idx_CG], atoms[cand_idx])
            locked.add(idx_CG)
        except RuntimeError:
            pass

    if idx_HG is not None and idx_CG is not None:
        atribuir_coord_alvos(atoms, locked, idx_CG, [idx_HG], dmin_init=0.995, dmax_init=1.115,
                            target_count=1, label="HG", scope_indices=scope_indices)

    if idx_CD1 is not None and idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CG, idx_CB, dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0, ang2_min=104.0, ang2_max=116.0, alvo_angulo=110.0, label="CD1",
                idx_angle_a=idx_CA, idx_angle_b=idx_CB, scope_indices=scope_indices)
            if not same_coords(atoms[idx_CD1], atoms[cand_idx]):
                swap_coords(atoms[idx_CD1], atoms[cand_idx])
            locked.add(idx_CD1)
        except RuntimeError:
            pass

    if idx_CD2 is not None and idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CG, idx_CB, dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0, ang2_min=104.0, ang2_max=116.0, alvo_angulo=110.0, label="CD2",
                idx_angle_a=idx_CA, idx_angle_b=idx_CB, scope_indices=scope_indices)
            if not same_coords(atoms[idx_CD2], atoms[cand_idx]):
                swap_coords(atoms[idx_CD2], atoms[cand_idx])
            locked.add(idx_CD2)
        except RuntimeError:
            pass

    alvo_hd1 = [name2idx.get(n) for n in ("1HD1", "2HD1", "3HD1") if name2idx.get(n) is not None]
    if alvo_hd1 and idx_CD1 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CD1, alvo_hd1, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hd1), label="HD1", scope_indices=scope_indices)

    alvo_hd2 = [name2idx.get(n) for n in ("1HD2", "2HD2", "3HD2") if name2idx.get(n) is not None]
    if alvo_hd2 and idx_CD2 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CD2, alvo_hd2, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hd2), label="HD2", scope_indices=scope_indices)


def processar_tyr(atoms: List[Dict], locked: Set[int], resseq: int = 579, chain: str = "B",
                  verbose: bool = True, scope_indices: List[int] = None) -> None:
    res_idxs = find_residue_atoms(atoms, "TYR", resseq, chain)
    if not res_idxs:
        return
    name2idx = build_name_index(atoms, res_idxs)
    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_CG = name2idx.get("CG")
    idx_CD1 = name2idx.get("CD1")
    idx_CD2 = name2idx.get("CD2")
    idx_CE1 = name2idx.get("CE1")
    idx_CE2 = name2idx.get("CE2")
    idx_CZ = name2idx.get("CZ")
    idx_OH = name2idx.get("OH")
    idx_HH = name2idx.get("HH")
    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")

    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(atoms, locked, idx_CB, alvo_hb, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hb), label="HB", scope_indices=scope_indices)

    if idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0, ang2_min=104.0, ang2_max=116.0, alvo_angulo=110.0, label="CG",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CG], atoms[cand_idx]):
                swap_coords(atoms[idx_CG], atoms[cand_idx])
            locked.add(idx_CG)
        except RuntimeError:
            pass

    if idx_CD1 is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CG, idx_CB, dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0, ang2_min=117.0, ang2_max=123.0, alvo_angulo=120.0, label="CD1",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CD1], atoms[cand_idx]):
                swap_coords(atoms[idx_CD1], atoms[cand_idx])
            locked.add(idx_CD1)
        except RuntimeError:
            pass

    idx_HD1 = name2idx.get("HD1")
    if idx_HD1 is not None and idx_CD1 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CD1, [idx_HD1], dmin_init=0.995, dmax_init=1.115,
                            target_count=1, label="HD1", scope_indices=scope_indices)

    if idx_CD2 is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CG, idx_CB, dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0, ang2_min=117.0, ang2_max=123.0, alvo_angulo=120.0, label="CD2",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CD2], atoms[cand_idx]):
                swap_coords(atoms[idx_CD2], atoms[cand_idx])
            locked.add(idx_CD2)
        except RuntimeError:
            pass

    idx_HD2 = name2idx.get("HD2")
    if idx_HD2 is not None and idx_CD2 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CD2, [idx_HD2], dmin_init=0.995, dmax_init=1.115,
                            target_count=1, label="HD2", scope_indices=scope_indices)

    if idx_CE1 is not None and idx_CD1 is not None and idx_CG is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CD1, idx_CG, dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0, ang2_min=117.0, ang2_max=123.0, alvo_angulo=120.0, label="CE1",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CE1], atoms[cand_idx]):
                swap_coords(atoms[idx_CE1], atoms[cand_idx])
            locked.add(idx_CE1)
        except RuntimeError:
            pass

    idx_HE1 = name2idx.get("HE1")
    if idx_HE1 is not None and idx_CE1 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CE1, [idx_HE1], dmin_init=0.995, dmax_init=1.115,
                            target_count=1, label="HE1", scope_indices=scope_indices)

    if idx_CE2 is not None and idx_CD2 is not None and idx_CG is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CD2, idx_CG, dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0, ang2_min=117.0, ang2_max=123.0, alvo_angulo=120.0, label="CE2",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CE2], atoms[cand_idx]):
                swap_coords(atoms[idx_CE2], atoms[cand_idx])
            locked.add(idx_CE2)
        except RuntimeError:
            pass

    idx_HE2 = name2idx.get("HE2")
    if idx_HE2 is not None and idx_CE2 is not None:
        atribuir_coord_alvos(atoms, locked, idx_CE2, [idx_HE2], dmin_init=0.995, dmax_init=1.115,
                            target_count=1, label="HE2", scope_indices=scope_indices)

    if idx_CZ is not None and idx_CE1 is not None and idx_CD1 is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CE1, idx_CD1, dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0, ang2_min=117.0, ang2_max=123.0, alvo_angulo=120.0, label="CZ",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_CZ], atoms[cand_idx]):
                swap_coords(atoms[idx_CZ], atoms[cand_idx])
            locked.add(idx_CZ)
        except RuntimeError:
            pass

    if idx_OH is not None and idx_CZ is not None and idx_CE1 is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CZ, idx_CE1, dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0, ang2_min=117.0, ang2_max=123.0, alvo_angulo=120.0, label="OH",
                scope_indices=scope_indices)
            if not same_coords(atoms[idx_OH], atoms[cand_idx]):
                swap_coords(atoms[idx_OH], atoms[cand_idx])
            locked.add(idx_OH)
        except RuntimeError:
            pass

    if idx_HH is not None and idx_OH is not None:
        atribuir_coord_alvos(atoms, locked, idx_OH, [idx_HH], dmin_init=0.85, dmax_init=1.05,
                            target_count=1, label="HH", scope_indices=scope_indices)


def processar_cys(atoms: List[Dict], locked: Set[int], resseq: int = 580, chain: str = "B",
                  verbose: bool = True, scope_indices: List[int] = None) -> None:
    res_idxs = find_residue_atoms(atoms, "CYS", resseq, chain)
    if not res_idxs:
        return
    name2idx = build_name_index(atoms, res_idxs)
    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_SG = name2idx.get("SG")
    idx_HG1 = name2idx.get("HG1")
    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")

    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(atoms, locked, idx_CB, alvo_hb, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hb), label="HB", scope_indices=scope_indices)

    if idx_SG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.60, dmax=1.90,
                ang1_min=100.0, ang1_max=125.0, ang2_min=104.0, ang2_max=121.0, alvo_angulo=112.5, label="SG",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_SG], atoms[cand_idx]):
                swap_coords(atoms[idx_SG], atoms[cand_idx])
            locked.add(idx_SG)
        except RuntimeError:
            pass

    if idx_HG1 is not None and idx_SG is not None:
        atribuir_coord_alvos(atoms, locked, idx_SG, [idx_HG1], dmin_init=1.20, dmax_init=1.45,
                            target_count=1, label="HG1", scope_indices=scope_indices)


def processar_lys(atoms: List[Dict], locked: Set[int], resseq: int = 581, chain: str = "B",
                  verbose: bool = True, scope_indices: List[int] = None) -> None:
    res_idxs = find_residue_atoms(atoms, "LYS", resseq, chain)
    if not res_idxs:
        return
    name2idx = build_name_index(atoms, res_idxs)
    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_CG = name2idx.get("CG")
    idx_CD = name2idx.get("CD")
    idx_CE = name2idx.get("CE")
    idx_NZ = name2idx.get("NZ")

    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")
    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(atoms, locked, idx_CB, alvo_hb, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hb), label="HB", scope_indices=scope_indices)

    if idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0, ang2_min=103.0, ang2_max=124.0, alvo_angulo=113.6, label="CG",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_CG], atoms[cand_idx]):
                swap_coords(atoms[idx_CG], atoms[cand_idx])
            locked.add(idx_CG)
        except RuntimeError:
            pass

    idx_HG1 = name2idx.get("HG1")
    idx_HG2 = name2idx.get("HG2")
    alvo_hg = [x for x in [idx_HG1, idx_HG2] if x is not None]
    if alvo_hg and idx_CG is not None:
        atribuir_coord_alvos(atoms, locked, idx_CG, alvo_hg, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hg), label="HG", scope_indices=scope_indices)

    if idx_CD is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CG, idx_CB, dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0, ang2_min=103.0, ang2_max=124.0, alvo_angulo=113.6, label="CD",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_CD], atoms[cand_idx]):
                swap_coords(atoms[idx_CD], atoms[cand_idx])
            locked.add(idx_CD)
        except RuntimeError:
            pass

    idx_HD1 = name2idx.get("HD1")
    idx_HD2 = name2idx.get("HD2")
    alvo_hd = [x for x in [idx_HD1, idx_HD2] if x is not None]
    if alvo_hd and idx_CD is not None:
        atribuir_coord_alvos(atoms, locked, idx_CD, alvo_hd, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hd), label="HD", scope_indices=scope_indices)

    if idx_CE is not None and idx_CD is not None and idx_CG is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CD, idx_CG, dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0, ang2_min=103.0, ang2_max=124.0, alvo_angulo=113.6, label="CE",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_CE], atoms[cand_idx]):
                swap_coords(atoms[idx_CE], atoms[cand_idx])
            locked.add(idx_CE)
        except RuntimeError:
            pass

    idx_HE1 = name2idx.get("HE1")
    idx_HE2 = name2idx.get("HE2")
    alvo_he = [x for x in [idx_HE1, idx_HE2] if x is not None]
    if alvo_he and idx_CE is not None:
        atribuir_coord_alvos(atoms, locked, idx_CE, alvo_he, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_he), label="HE", scope_indices=scope_indices)

    if idx_NZ is not None and idx_CE is not None and idx_CD is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CE, idx_CD, dmin=1.38, dmax=1.58,
                ang1_min=95.0, ang1_max=125.0, ang2_min=100.0, ang2_max=120.0, alvo_angulo=110.0, label="NZ",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_NZ], atoms[cand_idx]):
                swap_coords(atoms[idx_NZ], atoms[cand_idx])
            locked.add(idx_NZ)
        except RuntimeError:
            pass

    idx_HZ1 = name2idx.get("HZ1")
    idx_HZ2 = name2idx.get("HZ2")
    idx_HZ3 = name2idx.get("HZ3")
    alvo_hz = [x for x in [idx_HZ1, idx_HZ2, idx_HZ3] if x is not None]
    if alvo_hz and idx_NZ is not None:
        atribuir_coord_alvos(atoms, locked, idx_NZ, alvo_hz, dmin_init=0.90, dmax_init=1.10,
                            target_count=len(alvo_hz), label="HZ", scope_indices=scope_indices)


def processar_arg(atoms: List[Dict], locked: Set[int], resseq: int = 582, chain: str = "B",
                  verbose: bool = True, scope_indices: List[int] = None) -> None:
    res_idxs = find_residue_atoms(atoms, "ARG", resseq, chain)
    if not res_idxs:
        return
    name2idx = build_name_index(atoms, res_idxs)
    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_CG = name2idx.get("CG")
    idx_CD = name2idx.get("CD")
    idx_NE = name2idx.get("NE")
    idx_CZ = name2idx.get("CZ")
    idx_NH1 = name2idx.get("NH1")
    idx_NH2 = name2idx.get("NH2")

    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")
    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(atoms, locked, idx_CB, alvo_hb, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hb), label="HB", scope_indices=scope_indices)

    if idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CB, idx_CA, dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0, ang2_min=103.0, ang2_max=124.0, alvo_angulo=113.6, label="CG",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_CG], atoms[cand_idx]):
                swap_coords(atoms[idx_CG], atoms[cand_idx])
            locked.add(idx_CG)
        except RuntimeError:
            pass

    idx_HG1 = name2idx.get("HG1")
    idx_HG2 = name2idx.get("HG2")
    alvo_hg = [x for x in [idx_HG1, idx_HG2] if x is not None]
    if alvo_hg and idx_CG is not None:
        atribuir_coord_alvos(atoms, locked, idx_CG, alvo_hg, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hg), label="HG", scope_indices=scope_indices)

    if idx_CD is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CG, idx_CB, dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0, ang2_min=103.0, ang2_max=124.0, alvo_angulo=113.6, label="CD",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_CD], atoms[cand_idx]):
                swap_coords(atoms[idx_CD], atoms[cand_idx])
            locked.add(idx_CD)
        except RuntimeError:
            pass

    idx_HD1 = name2idx.get("HD1")
    idx_HD2 = name2idx.get("HD2")
    alvo_hd = [x for x in [idx_HD1, idx_HD2] if x is not None]
    if alvo_hd and idx_CD is not None:
        atribuir_coord_alvos(atoms, locked, idx_CD, alvo_hd, dmin_init=0.995, dmax_init=1.115,
                            target_count=len(alvo_hd), label="HD", scope_indices=scope_indices)

    if idx_NE is not None and idx_CD is not None and idx_CG is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CD, idx_CG, dmin=1.40, dmax=1.60,
                ang1_min=103.0, ang1_max=133.0, ang2_min=108.0, ang2_max=128.0, alvo_angulo=118.0, label="NE",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_NE], atoms[cand_idx]):
                swap_coords(atoms[idx_NE], atoms[cand_idx])
            locked.add(idx_NE)
        except RuntimeError:
            pass

    idx_HE = name2idx.get("HE")
    if idx_HE is not None and idx_NE is not None:
        atribuir_coord_alvos(atoms, locked, idx_NE, [idx_HE], dmin_init=0.90, dmax_init=1.10,
                            target_count=1, label="HE", scope_indices=scope_indices)

    if idx_CZ is not None and idx_NE is not None and idx_CD is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_NE, idx_CD, dmin=1.25, dmax=1.45,
                ang1_min=100.0, ang1_max=130.0, ang2_min=114.0, ang2_max=126.0, alvo_angulo=120.0, label="CZ",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_CZ], atoms[cand_idx]):
                swap_coords(atoms[idx_CZ], atoms[cand_idx])
            locked.add(idx_CZ)
        except RuntimeError:
            pass

    if idx_NH1 is not None and idx_CZ is not None and idx_NE is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CZ, idx_NE, dmin=1.25, dmax=1.45,
                ang1_min=100.0, ang1_max=130.0, ang2_min=114.0, ang2_max=126.0, alvo_angulo=120.0, label="NH1",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_NH1], atoms[cand_idx]):
                swap_coords(atoms[idx_NH1], atoms[cand_idx])
            locked.add(idx_NH1)
        except RuntimeError:
            pass

    idx_1HH1 = name2idx.get("1HH1")
    idx_2HH1 = name2idx.get("2HH1")
    alvo_hh1 = [x for x in [idx_1HH1, idx_2HH1] if x is not None]
    if alvo_hh1 and idx_NH1 is not None:
        atribuir_coord_alvos(atoms, locked, idx_NH1, alvo_hh1, dmin_init=0.90, dmax_init=1.10,
                            target_count=len(alvo_hh1), label="HH1", scope_indices=scope_indices)

    if idx_NH2 is not None and idx_CZ is not None and idx_NE is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(atoms, locked, idx_CZ, idx_NE, dmin=1.25, dmax=1.45,
                ang1_min=100.0, ang1_max=130.0, ang2_min=114.0, ang2_max=126.0, alvo_angulo=120.0, label="NH2",
                expand_distance=True, delta=0.005, max_iter=200, scope_indices=scope_indices)
            if not same_coords(atoms[idx_NH2], atoms[cand_idx]):
                swap_coords(atoms[idx_NH2], atoms[cand_idx])
            locked.add(idx_NH2)
        except RuntimeError:
            pass

    idx_1HH2 = name2idx.get("1HH2")
    idx_2HH2 = name2idx.get("2HH2")
    alvo_hh2 = [x for x in [idx_1HH2, idx_2HH2] if x is not None]
    if alvo_hh2 and idx_NH2 is not None:
        atribuir_coord_alvos(atoms, locked, idx_NH2, alvo_hh2, dmin_init=0.90, dmax_init=1.10,
                            target_count=len(alvo_hh2), label="HH2", scope_indices=scope_indices)


# ==============================================================================
# SEÇÃO 12: PROCESSAMENTO DE UM ÚNICO FRAME
# ==============================================================================

def processar_frame(frame_content: str, frame_num: int, start_serial: int = 8641,
                    end_serial: int = 8776, chain: str = "B", verbose: bool = False) -> Tuple[str, bool, str, List[Dict]]:
    atoms, lines = parse_frame_content(frame_content)

    if not atoms:
        return frame_content, False, "Frame vazio ou sem átomos válidos", []

    coords = extrair_coords_de_atoms(atoms)
    erro_msg = ""
    expansoes_backbone = []

    try:
        # FASE 1: Backbone
        coords, backbone_ids, expansoes_backbone = fase1_backbone(atoms, coords, start_serial, end_serial, verbose=False)
        aplicar_coords_para_atoms(atoms, coords)

        # FASE 2: HN, HA, O
        coords = extrair_coords_de_atoms(atoms)
        coords = fase2_hn_ha_o(atoms, coords, backbone_ids, verbose=False)
        aplicar_coords_para_atoms(atoms, coords)

        # FASE 3: CA-CB
        coords = extrair_coords_de_atoms(atoms)
        locked_serials = set(backbone_ids)
        atom_by_serial = {a['serial']: a for a in atoms}
        residues = get_backbone_residues(backbone_ids, atom_by_serial)
        for res in residues:
            for aname in ('HN', 'HA', 'O'):
                for a in atoms:
                    if a['chain'] == res['chain'] and a['resseq'] == res['resseq'] and a['name'] == aname:
                        locked_serials.add(a['serial'])
        coords, locked_serials = fase3_ca_cb(atoms, coords, backbone_ids, locked_serials, verbose=False)
        aplicar_coords_para_atoms(atoms, coords)

        # FASES 4-10: Cadeias laterais
        locked = inicializar_locked_base(atoms)

        # Calcular scope_indices para restringir buscas aos resíduos 576-583
        scope_indices = get_scope_indices(atoms, 576, 583, chain)

        processar_ile(atoms, locked, resseq=576, chain=chain, verbose=False, scope_indices=scope_indices)
        bloquear_residuo(atoms, locked, "ILE", 576, chain)

        processar_thr(atoms, locked, resseq=577, chain=chain, verbose=False, scope_indices=scope_indices)
        bloquear_residuo(atoms, locked, "THR", 577, chain)

        processar_leu(atoms, locked, resseq=578, chain=chain, verbose=False, scope_indices=scope_indices)
        bloquear_residuo(atoms, locked, "LEU", 578, chain)

        processar_tyr(atoms, locked, resseq=579, chain=chain, verbose=False, scope_indices=scope_indices)
        bloquear_residuo(atoms, locked, "TYR", 579, chain)

        processar_cys(atoms, locked, resseq=580, chain=chain, verbose=False, scope_indices=scope_indices)
        bloquear_residuo(atoms, locked, "CYS", 580, chain)

        processar_lys(atoms, locked, resseq=581, chain=chain, verbose=False, scope_indices=scope_indices)
        bloquear_residuo(atoms, locked, "LYS", 581, chain)

        processar_arg(atoms, locked, resseq=582, chain=chain, verbose=False, scope_indices=scope_indices)
        bloquear_residuo(atoms, locked, "ARG", 582, chain)

        # Atualizar linhas
        lines = update_pdb_lines(lines, atoms)
        return "\n".join(lines), True, "", expansoes_backbone

    except Exception as e:
        erro_msg = str(e)
        return "\n".join(lines), False, erro_msg, expansoes_backbone


# ==============================================================================
# SEÇÃO 13: PIPELINE PRINCIPAL MULTI-FRAME (para GUI)
# ==============================================================================

def executar_pipeline_multiframe_gui(pdb_entrada: str, pdb_saida: str, start_serial: int = 8641,
                                      end_serial: int = 8776, chain: str = "B",
                                      callback_log=None, callback_progress=None) -> dict:
    """
    Executa o pipeline completo para todos os frames de uma trajetória.
    Retorna um dicionário com estatísticas do processamento.
    """
    resultado = {
        'sucesso': False,
        'total_frames': 0,
        'frames_sucesso': 0,
        'frames_falha': 0,
        'mensagem': '',
        'arquivo_saida': pdb_saida
    }

    def log(msg):
        if callback_log:
            callback_log(msg)

    def progress(valor, maximo):
        if callback_progress:
            callback_progress(valor, maximo)

    log("="*60)
    log("PIPELINE CAR-T MULTI-FRAME v2.0")
    log("="*60)
    log(f"Arquivo de entrada: {pdb_entrada}")
    log(f"Arquivo de saída: {pdb_saida}")
    log("="*60)

    try:
        # Ler arquivo completo
        with open(pdb_entrada, 'r') as f:
            content = f.read()

        # Dividir em frames
        frames = split_trajectory_into_frames(content)
        total_frames = len(frames)
        resultado['total_frames'] = total_frames
        log(f"\nTotal de frames detectados: {total_frames}")

        if total_frames == 0:
            resultado['mensagem'] = "Nenhum frame encontrado no arquivo!"
            log("ERRO: Nenhum frame encontrado no arquivo!")
            return resultado

        # Processar cada frame
        processed_frames = []
        frames_com_falha = []
        frames_sucesso = 0
        todas_expansoes = []

        for i, frame_content in enumerate(frames):
            frame_num = i + 1
            progress(frame_num, total_frames)
            log(f"Processando frame {frame_num}/{total_frames}...")

            processed_frame, sucesso, erro_msg, expansoes = processar_frame(
                frame_content, frame_num, start_serial, end_serial, chain, verbose=False
            )

            processed_frames.append(processed_frame)

            if expansoes:
                todas_expansoes.append((frame_num, expansoes))

            if sucesso:
                frames_sucesso += 1
            else:
                frames_com_falha.append((frame_num, erro_msg))

        log("\nProcessamento concluído!")

        # Escrever arquivo de saída
        with open(pdb_saida, 'w') as f:
            f.write("\n".join(processed_frames))

        # Relatório final
        log("\n" + "="*60)
        log("RELATÓRIO DE PROCESSAMENTO")
        log("="*60)
        log(f"Total de frames:      {total_frames}")
        log(f"Frames com sucesso:   {frames_sucesso}")
        log(f"Frames com falha:     {len(frames_com_falha)}")

        resultado['frames_sucesso'] = frames_sucesso
        resultado['frames_falha'] = len(frames_com_falha)

        if frames_com_falha:
            log("\n" + "-"*60)
            log("FRAMES COM FALHA:")
            log("-"*60)

            diag_file = pdb_saida.replace(".pdb", "_diagnostico.txt")
            with open(diag_file, 'w') as f_diag:
                f_diag.write("="*80 + "\n")
                f_diag.write("RELATÓRIO DETALHADO - FRAMES COM FALHA\n")
                f_diag.write("="*80 + "\n")
                f_diag.write(f"Arquivo de entrada: {pdb_entrada}\n")
                f_diag.write(f"Total de frames: {total_frames}\n")
                f_diag.write(f"Frames com falha: {len(frames_com_falha)}\n")
                f_diag.write("="*80 + "\n\n")

                for frame_num, erro in frames_com_falha:
                    primeira_linha = erro.split('\n')[0] if '\n' in erro else erro
                    log(f"  Frame {frame_num}: {primeira_linha[:60]}...")
                    f_diag.write(f"\n{'#'*80}\n")
                    f_diag.write(f"# FRAME {frame_num}\n")
                    f_diag.write(f"{'#'*80}\n\n")
                    f_diag.write(erro)
                    f_diag.write("\n\n")

            log(f"\nDiagnóstico salvo em: {diag_file}")

        # Relatório de expansões
        if todas_expansoes:
            expansao_file = pdb_saida.replace(".pdb", "_expansoes.txt")
            with open(expansao_file, 'w') as f_exp:
                f_exp.write("=" * 80 + "\n")
                f_exp.write("RELATÓRIO DE EXPANSÕES DE JANELA DO BACKBONE\n")
                f_exp.write("=" * 80 + "\n")
                f_exp.write(f"Frames com expansão: {len(todas_expansoes)}\n")
                f_exp.write("=" * 80 + "\n\n")

                for frame_num, expansoes in todas_expansoes:
                    relatorio = formatar_relatorio_expansoes(expansoes, frame_num)
                    f_exp.write(relatorio)
                    f_exp.write("\n\n")

            log(f"Expansões salvas em: {expansao_file}")

        log("\n" + "="*60)
        log(f"Arquivo de saída: {pdb_saida}")
        log("="*60)
        log("PIPELINE CONCLUÍDO COM SUCESSO!")
        log("="*60)

        resultado['sucesso'] = True
        resultado['mensagem'] = "Processamento concluído com sucesso!"
        return resultado

    except Exception as e:
        resultado['mensagem'] = f"Erro durante processamento: {str(e)}"
        log(f"\nERRO: {str(e)}")
        return resultado


# ==============================================================================
# SEÇÃO 14: INTERFACE GRÁFICA (TKINTER)
# ==============================================================================

class CARTReorderApp:
    """Aplicação GUI para reordenação de estruturas CAR-T."""

    def __init__(self, root):
        self.root = root
        self.root.title("CAR-T Reorder Pipeline v2.0")
        self.root.geometry("800x700")
        self.root.minsize(700, 600)

        # Variáveis
        self.arquivo_entrada = tk.StringVar()
        self.arquivo_saida = tk.StringVar()
        self.start_serial = tk.StringVar(value="8641")
        self.end_serial = tk.StringVar(value="8776")
        self.chain = tk.StringVar(value="B")

        # Queue para comunicação entre threads
        self.log_queue = queue.Queue()
        self.progress_queue = queue.Queue()

        # Flag de processamento
        self.processando = False

        self.criar_interface()
        self.iniciar_atualizacao_log()

    def criar_interface(self):
        """Cria todos os widgets da interface."""

        # Frame principal com padding
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # ===== TÍTULO =====
        titulo = ttk.Label(main_frame, text="Pipeline de Reordenação CAR-T",
                          font=('Helvetica', 16, 'bold'))
        titulo.pack(pady=(0, 10))

        subtitulo = ttk.Label(main_frame, text="Versão 2.0 - Interface Gráfica",
                             font=('Helvetica', 10))
        subtitulo.pack(pady=(0, 20))

        # ===== SELEÇÃO DE ARQUIVOS =====
        arquivos_frame = ttk.LabelFrame(main_frame, text="Arquivos", padding="10")
        arquivos_frame.pack(fill=tk.X, pady=(0, 10))

        # Arquivo de entrada
        ttk.Label(arquivos_frame, text="Arquivo PDB de Entrada:").grid(row=0, column=0, sticky=tk.W, pady=5)
        entrada_frame = ttk.Frame(arquivos_frame)
        entrada_frame.grid(row=0, column=1, sticky=tk.EW, pady=5)
        ttk.Entry(entrada_frame, textvariable=self.arquivo_entrada, width=50).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(entrada_frame, text="Procurar...", command=self.selecionar_entrada).pack(side=tk.LEFT, padx=(5, 0))

        # Arquivo de saída
        ttk.Label(arquivos_frame, text="Arquivo PDB de Saída:").grid(row=1, column=0, sticky=tk.W, pady=5)
        saida_frame = ttk.Frame(arquivos_frame)
        saida_frame.grid(row=1, column=1, sticky=tk.EW, pady=5)
        ttk.Entry(saida_frame, textvariable=self.arquivo_saida, width=50).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(saida_frame, text="Procurar...", command=self.selecionar_saida).pack(side=tk.LEFT, padx=(5, 0))

        arquivos_frame.columnconfigure(1, weight=1)

        # ===== PARÂMETROS =====
        params_frame = ttk.LabelFrame(main_frame, text="Parâmetros", padding="10")
        params_frame.pack(fill=tk.X, pady=(0, 10))

        # Start Serial
        ttk.Label(params_frame, text="Serial Inicial:").grid(row=0, column=0, sticky=tk.W, pady=5, padx=(0, 10))
        ttk.Entry(params_frame, textvariable=self.start_serial, width=15).grid(row=0, column=1, sticky=tk.W, pady=5)

        # End Serial
        ttk.Label(params_frame, text="Serial Final:").grid(row=0, column=2, sticky=tk.W, pady=5, padx=(20, 10))
        ttk.Entry(params_frame, textvariable=self.end_serial, width=15).grid(row=0, column=3, sticky=tk.W, pady=5)

        # Chain
        ttk.Label(params_frame, text="Chain:").grid(row=0, column=4, sticky=tk.W, pady=5, padx=(20, 10))
        ttk.Entry(params_frame, textvariable=self.chain, width=5).grid(row=0, column=5, sticky=tk.W, pady=5)

        # ===== BARRA DE PROGRESSO =====
        progress_frame = ttk.Frame(main_frame)
        progress_frame.pack(fill=tk.X, pady=(0, 10))

        self.progress_label = ttk.Label(progress_frame, text="Aguardando...")
        self.progress_label.pack(anchor=tk.W)

        self.progress_bar = ttk.Progressbar(progress_frame, mode='determinate', length=400)
        self.progress_bar.pack(fill=tk.X, pady=5)

        # ===== BOTÕES =====
        botoes_frame = ttk.Frame(main_frame)
        botoes_frame.pack(fill=tk.X, pady=(0, 10))

        self.btn_executar = ttk.Button(botoes_frame, text="Executar Pipeline",
                                        command=self.executar_pipeline, style='Accent.TButton')
        self.btn_executar.pack(side=tk.LEFT, padx=(0, 10))

        self.btn_limpar = ttk.Button(botoes_frame, text="Limpar Log", command=self.limpar_log)
        self.btn_limpar.pack(side=tk.LEFT)

        # ===== ÁREA DE LOG =====
        log_frame = ttk.LabelFrame(main_frame, text="Log de Execução", padding="5")
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.log_text = scrolledtext.ScrolledText(log_frame, height=15, font=('Consolas', 9))
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.config(state=tk.DISABLED)

        # ===== RODAPÉ =====
        rodape = ttk.Label(main_frame, text="Pipeline CAR-T v2.0 | 2025",
                          font=('Helvetica', 8))
        rodape.pack(pady=(10, 0))

    def selecionar_entrada(self):
        """Abre diálogo para selecionar arquivo de entrada."""
        filename = filedialog.askopenfilename(
            title="Selecionar arquivo PDB de entrada",
            filetypes=[("Arquivos PDB", "*.pdb"), ("Todos os arquivos", "*.*")]
        )
        if filename:
            self.arquivo_entrada.set(filename)
            # Sugerir nome de saída
            base = os.path.splitext(filename)[0]
            self.arquivo_saida.set(f"{base}_reordenado.pdb")

    def selecionar_saida(self):
        """Abre diálogo para selecionar arquivo de saída."""
        filename = filedialog.asksaveasfilename(
            title="Salvar arquivo PDB de saída",
            defaultextension=".pdb",
            filetypes=[("Arquivos PDB", "*.pdb"), ("Todos os arquivos", "*.*")]
        )
        if filename:
            self.arquivo_saida.set(filename)

    def adicionar_log(self, mensagem):
        """Adiciona mensagem à queue de log."""
        self.log_queue.put(mensagem)

    def atualizar_progresso(self, valor, maximo):
        """Atualiza a barra de progresso."""
        self.progress_queue.put((valor, maximo))

    def iniciar_atualizacao_log(self):
        """Inicia o loop de atualização do log."""
        self.processar_queues()
        self.root.after(100, self.iniciar_atualizacao_log)

    def processar_queues(self):
        """Processa as filas de log e progresso."""
        # Processar logs
        while not self.log_queue.empty():
            try:
                msg = self.log_queue.get_nowait()
                self.log_text.config(state=tk.NORMAL)
                self.log_text.insert(tk.END, msg + "\n")
                self.log_text.see(tk.END)
                self.log_text.config(state=tk.DISABLED)
            except queue.Empty:
                break

        # Processar progresso
        while not self.progress_queue.empty():
            try:
                valor, maximo = self.progress_queue.get_nowait()
                self.progress_bar['maximum'] = maximo
                self.progress_bar['value'] = valor
                self.progress_label.config(text=f"Frame {valor}/{maximo}")
            except queue.Empty:
                break

    def limpar_log(self):
        """Limpa a área de log."""
        self.log_text.config(state=tk.NORMAL)
        self.log_text.delete(1.0, tk.END)
        self.log_text.config(state=tk.DISABLED)
        self.progress_bar['value'] = 0
        self.progress_label.config(text="Aguardando...")

    def validar_parametros(self):
        """Valida os parâmetros de entrada."""
        if not self.arquivo_entrada.get():
            messagebox.showerror("Erro", "Selecione o arquivo PDB de entrada!")
            return False

        if not os.path.exists(self.arquivo_entrada.get()):
            messagebox.showerror("Erro", "Arquivo de entrada não encontrado!")
            return False

        if not self.arquivo_saida.get():
            messagebox.showerror("Erro", "Especifique o arquivo PDB de saída!")
            return False

        try:
            int(self.start_serial.get())
            int(self.end_serial.get())
        except ValueError:
            messagebox.showerror("Erro", "Os valores de Serial devem ser números inteiros!")
            return False

        if not self.chain.get():
            messagebox.showerror("Erro", "Especifique a Chain!")
            return False

        return True

    def executar_pipeline(self):
        """Inicia a execução do pipeline em uma thread separada."""
        if self.processando:
            messagebox.showwarning("Aviso", "Processamento já em andamento!")
            return

        if not self.validar_parametros():
            return

        self.processando = True
        self.btn_executar.config(state=tk.DISABLED)
        self.limpar_log()

        # Executar em thread separada
        thread = threading.Thread(target=self.executar_pipeline_thread)
        thread.daemon = True
        thread.start()

    def executar_pipeline_thread(self):
        """Executa o pipeline (chamado pela thread)."""
        try:
            resultado = executar_pipeline_multiframe_gui(
                pdb_entrada=self.arquivo_entrada.get(),
                pdb_saida=self.arquivo_saida.get(),
                start_serial=int(self.start_serial.get()),
                end_serial=int(self.end_serial.get()),
                chain=self.chain.get(),
                callback_log=self.adicionar_log,
                callback_progress=self.atualizar_progresso
            )

            # Mostrar resultado final
            self.root.after(0, lambda: self.finalizar_processamento(resultado))

        except Exception as e:
            self.adicionar_log(f"\nERRO CRÍTICO: {str(e)}")
            self.root.after(0, lambda: self.finalizar_processamento({'sucesso': False, 'mensagem': str(e)}))

    def finalizar_processamento(self, resultado):
        """Finaliza o processamento e reativa a interface."""
        self.processando = False
        self.btn_executar.config(state=tk.NORMAL)

        if resultado['sucesso']:
            messagebox.showinfo("Sucesso",
                f"Processamento concluído!\n\n"
                f"Frames processados: {resultado.get('frames_sucesso', 0)}/{resultado.get('total_frames', 0)}\n"
                f"Arquivo salvo em:\n{resultado.get('arquivo_saida', '')}")
        else:
            messagebox.showerror("Erro", resultado.get('mensagem', 'Erro desconhecido'))


# ==============================================================================
# SEÇÃO 15: PONTO DE ENTRADA
# ==============================================================================

def main():
    """Função principal - inicia a aplicação GUI."""
    root = tk.Tk()

    # Tentar aplicar tema moderno se disponível
    try:
        root.tk.call('source', 'azure.tcl')
        root.tk.call('set_theme', 'light')
    except:
        pass

    app = CARTReorderApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
