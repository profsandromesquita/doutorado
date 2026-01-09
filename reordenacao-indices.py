# -*- coding: utf-8 -*-
"""
================================================================================
PIPELINE UNIFICADO DE REORDENAÇÃO DE ESTRUTURA CAR-T - MULTI-FRAME (v2 CORRIGIDO)
================================================================================

Versão que processa múltiplos frames de uma trajetória de dinâmica molecular.
Cada frame é separado por TER/ENDMDL.

AUTOR: Pipeline consolidado automaticamente
DATA: 2025
VERSÃO: 2.0 - Correção da função fase3_ca_cb com expansão dinâmica de janela
================================================================================
"""

import math
import sys
from typing import List, Dict, Set, Tuple, Optional, Any
from collections import defaultdict


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
    """
    Calcula o ângulo em graus formado por três pontos (p1-p2-p3).
    O ângulo é medido no vértice p2.
    """
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
    """
    Divide o conteúdo de uma trajetória em frames individuais.
    Cada frame termina com TER seguido de ENDMDL.
    """
    frames = []
    current_frame_lines = []
    lines = content.splitlines()

    for line in lines:
        current_frame_lines.append(line)
        if line.strip() == "ENDMDL":
            frames.append("\n".join(current_frame_lines))
            current_frame_lines = []

    # Se sobrou conteúdo sem ENDMDL (último frame incompleto)
    if current_frame_lines:
        remaining = "\n".join(current_frame_lines).strip()
        if remaining:
            frames.append(remaining)

    return frames


def parse_frame_content(content: str) -> Tuple[List[Dict], List[str]]:
    """
    Faz o parsing do conteúdo de um único frame.
    """
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
    """
    Atualiza as coordenadas nas linhas do PDB.
    """
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
    """
    Retorna lista de seriais dos átomos que estão dentro do escopo definido.

    Args:
        atoms: Lista de todos os átomos
        res_start: Número do primeiro resíduo do escopo (ex: 576)
        res_end: Número do último resíduo do escopo (ex: 583)
        chain: Cadeia alvo (ex: "B")

    Returns:
        Lista de seriais dos átomos no escopo
    """
    return [
        a['serial'] for a in atoms
        if a['chain'] == chain and res_start <= a['resseq'] <= res_end
    ]


def get_scope_indices(atoms: List[Dict], res_start: int, res_end: int, chain: str) -> List[int]:
    """
    Retorna lista de índices dos átomos que estão dentro do escopo definido.

    Args:
        atoms: Lista de todos os átomos
        res_start: Número do primeiro resíduo do escopo (ex: 576)
        res_end: Número do último resíduo do escopo (ex: 583)
        chain: Cadeia alvo (ex: "B")

    Returns:
        Lista de índices dos átomos no escopo
    """
    return [
        i for i, a in enumerate(atoms)
        if a['chain'] == chain and res_start <= a['resseq'] <= res_end
    ]


# ==============================================================================
# SEÇÃO 4: CONSTRUÇÃO DO BACKBONE CANÔNICO
# ==============================================================================

def build_backbone_order(atoms: List[Dict], atom_by_serial: Dict[int, Dict],
                         start_serial: int, end_serial: int) -> List[int]:
    if start_serial not in atom_by_serial:
        raise ValueError(f"Start_serial {start_serial} não encontrado no PDB.")
    if end_serial not in atom_by_serial:
        raise ValueError(f"End_serial {end_serial} não encontrado no PDB.")

    start_atom = atom_by_serial[start_serial]
    end_atom = atom_by_serial[end_serial]
    chain = start_atom["chain"]
    res_start = start_atom["resseq"]
    res_end = end_atom["resseq"]

    if chain != end_atom["chain"]:
        raise ValueError("Start e end estão em cadeias diferentes.")

    backbone_ids = []
    for res in range(res_start, res_end + 1):
        for aname in ("N", "CA", "C"):
            candidates = [
                a for a in atoms
                if a["chain"] == chain and a["resseq"] == res and a["name"] == aname
            ]
            if not candidates:
                raise ValueError(f"Não encontrei átomo {aname} no resíduo {res} cadeia {chain}.")
            backbone_ids.append(candidates[0]["serial"])

    if backbone_ids[0] != start_serial:
        raise ValueError(f"Primeiro átomo do backbone ({backbone_ids[0]}) não é o start_serial ({start_serial}).")
    if backbone_ids[-1] != end_serial:
        raise ValueError(f"Último átomo do backbone ({backbone_ids[-1]}) não é o end_serial ({end_serial}).")

    return backbone_ids


def get_backbone_residues(backbone_ids: List[int], atom_by_serial: Dict[int, Dict]) -> List[Dict]:
    residues = []
    if len(backbone_ids) % 3 != 0:
        raise ValueError("backbone_ids não está múltiplo de 3 (N, CA, C).")

    for i in range(0, len(backbone_ids), 3):
        N_id, CA_id, C_id = backbone_ids[i:i+3]
        N_atom = atom_by_serial[N_id]
        CA_atom = atom_by_serial[CA_id]
        C_atom = atom_by_serial[C_id]

        if not (N_atom["resseq"] == CA_atom["resseq"] == C_atom["resseq"] and
                N_atom["chain"] == CA_atom["chain"] == C_atom["chain"]):
            raise ValueError("Inconsistência N/CA/C no mesmo resíduo.")

        residues.append({
            "resseq": CA_atom["resseq"],
            "resname": CA_atom["resname"],
            "chain": CA_atom["chain"],
            "N": N_id,
            "CA": CA_id,
            "C": C_id,
        })

    return residues


# ==============================================================================
# SEÇÃO 5: LIMITES DE LIGAÇÃO QUÍMICA
# ==============================================================================

def bond_limits_backbone(name1: str, name2: str) -> Optional[Tuple[float, float]]:
    """
    Retorna os limites de distância (min, max) para ligações do backbone.

    Limites ajustados:
    - N-CA / CA-N: 1.30 - 1.60 Å
    - CA-C / C-CA: 1.30 - 1.70 Å
    - C-N / N-C:   1.15 - 1.45 Å
    """
    pair = (name1, name2)
    if pair in (("N", "CA"), ("CA", "N")):
        return (1.30, 1.60)
    if pair in (("CA", "C"), ("C", "CA")):
        return (1.30, 1.70)
    if pair in (("C", "N"), ("N", "C")):
        return (1.15, 1.45)
    return None


# ==============================================================================
# SEÇÃO 6: ALGORITMOS DE SELEÇÃO DE ÁTOMOS
# ==============================================================================

def escolher_vizinhos_dinamico(atoms: List[Dict], locked: Set[int], idx_ref: int,
                                target_count: int, dmin_init: float, dmax_init: float,
                                delta: float = 0.005, max_iter: int = 200,
                                label: str = "H?", scope_indices: List[int] = None) -> Tuple[List[int], float, float]:
    """
    Busca vizinhos dinâmicos com restrição de escopo.

    Args:
        scope_indices: Se fornecido, restringe a busca apenas a estes índices.
                      Se None, busca em todos os átomos (comportamento legado).
    """
    ref = atoms[idx_ref]
    dmin = dmin_init
    dmax = dmax_init

    # CORREÇÃO: Se scope_indices fornecido, usar apenas esses índices
    indices_to_search = scope_indices if scope_indices is not None else range(len(atoms))

    for _ in range(max_iter):
        cand = []
        for i in indices_to_search:
            if i in locked or i == idx_ref:
                continue
            a = atoms[i]
            d = dist(ref, a)
            if dmin <= d <= dmax:
                cand.append((i, d))
        if len(cand) == target_count:
            return [i for (i, _) in cand], dmin, dmax
        elif len(cand) > target_count:
            dmin += delta
            dmax -= delta
            if dmin >= dmax:
                break
        else:
            dmin = max(0.0, dmin - delta)
            dmax += delta

    cand_all = []
    for i in indices_to_search:
        if i in locked or i == idx_ref:
            continue
        a = atoms[i]
        d = dist(ref, a)
        cand_all.append((i, d))
    cand_all.sort(key=lambda t: t[1])
    chosen = cand_all[:target_count]
    return [i for (i, _) in chosen], dmin, dmax


def atribuir_coord_alvos(atoms: List[Dict], locked: Set[int], idx_ref: int,
                         alvo_idxs: List[int], dmin_init: float, dmax_init: float,
                         target_count: int, delta: float = 0.005, max_iter: int = 200,
                         label: str = "H?", verbose: bool = False,
                         scope_indices: List[int] = None) -> Dict:
    """
    Atribui coordenadas aos átomos alvo com restrição de escopo.

    Args:
        scope_indices: Se fornecido, restringe a busca apenas a estes índices.
    """
    ref = atoms[idx_ref]
    cand_idxs, dmin_final, dmax_final = escolher_vizinhos_dinamico(
        atoms, locked, idx_ref, target_count, dmin_init, dmax_init, delta, max_iter, label,
        scope_indices=scope_indices
    )
    cand_coords = {i: (atoms[i]['x'], atoms[i]['y'], atoms[i]['z']) for i in cand_idxs}
    used_cands = set()
    detalhes = {'ref': ref, 'ref_idx': idx_ref, 'janela_final': (dmin_final, dmax_final), 'mapeamentos': []}

    for idx_alvo in alvo_idxs:
        alvo = atoms[idx_alvo]
        match = None
        for i_cand in cand_idxs:
            if i_cand in used_cands:
                continue
            cx, cy, cz = cand_coords[i_cand]
            if abs(alvo['x'] - cx) < 1e-3 and abs(alvo['y'] - cy) < 1e-3 and abs(alvo['z'] - cz) < 1e-3:
                match = i_cand
                break
        if match is not None:
            used_cands.add(match)
            locked.add(idx_alvo)
            d = dist(ref, alvo)
            detalhes['mapeamentos'].append({'alvo_idx': idx_alvo, 'alvo': alvo, 'cand_idx': match,
                                            'cand': atoms[match], 'dist_ref_alvo': d, 'swap_feito': False})

    for idx_alvo in alvo_idxs:
        if idx_alvo in locked:
            continue
        alvo = atoms[idx_alvo]
        cand_rest = [i for i in cand_idxs if i not in used_cands]
        if not cand_rest:
            locked.add(idx_alvo)
            detalhes['mapeamentos'].append({'alvo_idx': idx_alvo, 'alvo': alvo, 'cand_idx': None,
                                            'cand': None, 'dist_ref_alvo': dist(ref, alvo), 'swap_feito': False})
            continue
        cand_rest.sort(key=lambda i: dist(ref, atoms[i]))
        i_cand = cand_rest[0]
        used_cands.add(i_cand)
        cand_atom = atoms[i_cand]
        d_before = dist(ref, alvo)
        if not same_coords(alvo, cand_atom):
            swap_coords(alvo, cand_atom)
            swap_feito = True
        else:
            swap_feito = False
        locked.add(idx_alvo)
        d_after = dist(ref, alvo)
        detalhes['mapeamentos'].append({'alvo_idx': idx_alvo, 'alvo': alvo, 'cand_idx': i_cand,
                                        'cand': cand_atom, 'dist_ref_alvo_antes': d_before,
                                        'dist_ref_alvo_depois': d_after, 'swap_feito': swap_feito})
    return detalhes


def escolher_pesado_com_angulo(atoms: List[Dict], locked: Set[int], idx_center: int, idx_prev: int,
                                dmin: float, dmax: float, ang1_min: float = 101.0, ang1_max: float = 119.0,
                                ang2_min: float = 104.0, ang2_max: float = 116.0, alvo_angulo: float = 110.0,
                                label: str = "PESO", idx_angle_a: int = None, idx_angle_b: int = None,
                                expand_distance: bool = False, delta: float = 0.005,
                                max_iter: int = 200, scope_indices: List[int] = None) -> Tuple[int, float, float]:
    """
    Escolhe átomo pesado com critério de ângulo e restrição de escopo.

    Args:
        scope_indices: Se fornecido, restringe a busca apenas a estes índices.
    """
    ref_dist = atoms[idx_center]
    if idx_angle_a is None:
        idx_angle_a = idx_prev
    if idx_angle_b is None:
        idx_angle_b = idx_center
    ref_angle_a = atoms[idx_angle_a]
    ref_angle_b = atoms[idx_angle_b]

    dmin_cur = dmin
    dmax_cur = dmax
    candidatos = []

    # CORREÇÃO: Se scope_indices fornecido, usar apenas esses índices
    indices_to_search = scope_indices if scope_indices is not None else range(len(atoms))

    for iteration in range(max_iter if expand_distance else 1):
        candidatos = []
        for i in indices_to_search:
            if i in locked or i == idx_center:
                continue
            a = atoms[i]
            d = dist(ref_dist, a)
            if dmin_cur <= d <= dmax_cur:
                ang = angle(ref_angle_a, ref_angle_b, a)
                candidatos.append((i, d, ang))
        if candidatos:
            break
        if expand_distance:
            dmin_cur = max(0.0, dmin_cur - delta)
            dmax_cur = dmax_cur + delta

    if not candidatos:
        raise RuntimeError(f"Nenhum candidato encontrado para {label}.")

    cand1 = [(i, d, ang) for (i, d, ang) in candidatos if ang1_min <= ang <= ang1_max]
    if len(cand1) == 1:
        return cand1[0]
    if len(cand1) > 1:
        cand2 = [(i, d, ang) for (i, d, ang) in cand1 if ang2_min <= ang <= ang2_max]
        if len(cand2) == 1:
            return cand2[0]
        if len(cand2) > 1:
            cand2.sort(key=lambda t: abs(t[2] - alvo_angulo))
            return cand2[0]
        cand1.sort(key=lambda t: abs(t[2] - alvo_angulo))
        return cand1[0]
    candidatos.sort(key=lambda t: abs(t[2] - alvo_angulo))
    return candidatos[0]


# ==============================================================================
# SEÇÃO 7: FASE 1 - RECONSTRUÇÃO DO BACKBONE (DFS)
# ==============================================================================

def get_backbone_ideal_angle(prev_name: str, curr_name: str, next_name: str) -> float:
    """
    Retorna o ângulo ideal (em graus) para a sequência de átomos do backbone.

    Ângulos típicos do backbone:
    - N-CA-C:  ~111° (carbono alfa tetraédrico)
    - CA-C-N:  ~116° (ligação peptídica planar)
    - C-N-CA:  ~121° (ligação peptídica planar)
    """
    triplet = (prev_name, curr_name, next_name)

    # Mapeamento de ângulos ideais para diferentes sequências
    angles = {
        ("N", "CA", "C"): 111.0,
        ("CA", "C", "N"): 116.0,
        ("C", "N", "CA"): 121.0,
    }

    return angles.get(triplet, 110.0)  # Default para 110° se não encontrado


def find_candidates_dynamic_window(
    all_serials: List[int],
    locked: Set[int],
    coords: Dict[int, Tuple[float, float, float]],
    curr_coord: Tuple[float, float, float],
    prev_coord: Optional[Tuple[float, float, float]],
    d_min_base: float,
    d_max_base: float,
    ideal_angle: float,
    delta: float = 0.02
) -> Tuple[List[Tuple[int, float]], float, float, int]:
    """
    Busca candidatos com janela dinâmica de distância.

    Lógica de expansão:
    1. Primeiro tenta com limites originais
    2. Se não encontrou, expande o máximo em +delta até 15 tentativas
    3. Se ainda não encontrou, reduz o mínimo em -delta até 15 tentativas
    4. Se ainda não encontrou, continua expandindo ambos indefinidamente

    Args:
        all_serials: Lista de todos os seriais de átomos
        locked: Conjunto de seriais já utilizados
        coords: Dicionário de coordenadas atuais
        curr_coord: Coordenada do átomo atual
        prev_coord: Coordenada do átomo anterior (para cálculo de ângulo)
        d_min_base: Distância mínima inicial (pode ser reduzida)
        d_max_base: Distância máxima inicial (pode ser expandida)
        ideal_angle: Ângulo ideal para refinamento
        delta: Incremento de expansão (default: 0.02 Å)

    Returns:
        Tupla (candidatos, d_min_final, d_max_final, num_expansoes)
        candidatos: Lista de tuplas (serial, distância)
    """
    d_min = d_min_base
    d_max = d_max_base
    num_expansoes = 0
    max_expansoes_por_direcao = 15

    def buscar_na_janela(d_min_busca, d_max_busca):
        cands = []
        for cand_id in all_serials:
            if cand_id in locked:
                continue
            cand_coord = coords[cand_id]
            d = dist_tuple(curr_coord, cand_coord)
            if d_min_busca <= d <= d_max_busca:
                cands.append((cand_id, d))
        return cands

    # 1. Tentar com limites originais
    candidatos = buscar_na_janela(d_min, d_max)

    # 2. Se não encontrou, expandir MAX até max_expansoes_por_direcao tentativas
    if not candidatos:
        for _ in range(max_expansoes_por_direcao):
            d_max += delta
            num_expansoes += 1
            candidatos = buscar_na_janela(d_min, d_max)
            if candidatos:
                break

    # 3. Se ainda não encontrou, reduzir MIN até max_expansoes_por_direcao tentativas
    if not candidatos:
        for _ in range(max_expansoes_por_direcao):
            d_min = max(0.0, d_min - delta)
            num_expansoes += 1
            candidatos = buscar_na_janela(d_min, d_max)
            if candidatos:
                break

    # 4. Se AINDA não encontrou, continuar expandindo ambos indefinidamente
    while not candidatos:
        d_max += delta
        d_min = max(0.0, d_min - delta)
        num_expansoes += 1
        candidatos = buscar_na_janela(d_min, d_max)

    # Aplicar refinamento angular se houver múltiplos candidatos
    if len(candidatos) > 1 and prev_coord is not None:
        candidatos_com_angulo = []
        for cand_id, d in candidatos:
            cand_coord = coords[cand_id]
            ang = angle_tuple(prev_coord, curr_coord, cand_coord)
            diff_angulo = abs(ang - ideal_angle)
            candidatos_com_angulo.append((cand_id, d, ang, diff_angulo))

        # Ordenar por diferença de ângulo (menor primeiro)
        candidatos_com_angulo.sort(key=lambda x: x[3])
        candidatos = [(c[0], c[1]) for c in candidatos_com_angulo]

    return candidatos, d_min, d_max, num_expansoes


def find_all_backbone_paths_dfs(atoms: List[Dict], atom_by_serial: Dict[int, Dict],
                                 backbone_ids: List[int], coords: Dict[int, Tuple[float, float, float]],
                                 max_steps: int = 23, delta: float = 0.02) -> Tuple[List[Dict], Dict]:
    """
    Encontra caminho válido para o backbone usando DFS completo.

    Testa TODAS as possibilidades de candidatos em cada passo, usando backtracking
    para explorar diferentes combinações e encontrar a melhor solução global.

    Lógica de expansão de janela:
    1. Primeiro tenta com limites originais
    2. Se não encontrou, expande o máximo em +delta até 15 tentativas
    3. Se ainda não encontrou, reduz o mínimo em -delta até 15 tentativas

    Args:
        atoms: Lista de átomos
        atom_by_serial: Dicionário de átomos por serial
        backbone_ids: Sequência canônica de IDs do backbone
        coords: Coordenadas atuais
        max_steps: Número máximo de passos (default: 23)
        delta: Incremento de expansão (default: 0.02 Å)

    Returns:
        Tupla (solutions, stats)
    """
    # CORREÇÃO: Filtrar apenas átomos do escopo (576-583)
    # Determinar escopo a partir do backbone_ids
    start_atom = atom_by_serial[backbone_ids[0]]
    end_atom = atom_by_serial[backbone_ids[-1]]
    res_start = start_atom['resseq']
    res_end = end_atom['resseq']
    chain_alvo = start_atom['chain']

    all_serials = get_scope_serials(atoms, res_start, res_end, chain_alvo)

    L = len(backbone_ids)
    if L - 1 != max_steps:
        raise ValueError(f"Número de passos ({L-1}) não bate com max_steps={max_steps}.")

    start_id = backbone_ids[0]
    stats = {
        "total_caminhos_testados": 0,
        "caminhos_completos_encontrados": 0,
        "melhor_desvio_angular": float('inf'),
        "expansoes_janela": []
    }

    # Estado inicial
    initial_state = {
        "step": 0,
        "coords": coords.copy(),
        "locked": {start_id},
        "total_dist": 0.0,
        "total_angular_deviation": 0.0,
        "edges": [],
        "expansoes": []
    }

    # Lista para armazenar todas as soluções válidas
    all_solutions = []

    # Pilha para DFS (usando lista como pilha)
    stack = [initial_state]

    while stack:
        state = stack.pop()
        k = state["step"]
        stats["total_caminhos_testados"] += 1

        # Se completou todos os passos, é uma solução válida
        if k == max_steps:
            stats["caminhos_completos_encontrados"] += 1
            if state["total_angular_deviation"] < stats["melhor_desvio_angular"]:
                stats["melhor_desvio_angular"] = state["total_angular_deviation"]
            all_solutions.append(state)
            continue

        # Informações do passo atual
        curr_id = backbone_ids[k]
        next_id = backbone_ids[k + 1]
        curr_atom = atom_by_serial[curr_id]
        next_atom = atom_by_serial[next_id]

        # Obter limites de distância para este par
        limits = bond_limits_backbone(curr_atom['name'], next_atom['name'])
        if limits is None:
            continue  # Pular se par inválido
        d_min, d_max_base = limits

        # Coordenada atual
        curr_coord = state["coords"][curr_id]

        # Coordenada anterior para cálculo de ângulo
        prev_coord = None
        ideal_angle = 110.0
        if k > 0:
            prev_id = backbone_ids[k - 1]
            prev_coord = state["coords"][prev_id]
            prev_name = atom_by_serial[prev_id]['name']
            ideal_angle = get_backbone_ideal_angle(prev_name, curr_atom['name'], next_atom['name'])

        # Buscar TODOS os candidatos com expansão dinâmica
        d_min_base = d_min  # Guardar limite mínimo original
        d_max = d_max_base
        d_min_atual = d_min_base
        candidatos = []
        num_expansoes = 0
        max_expansoes_por_direcao = 15  # Limite de tentativas por direção

        # Função auxiliar para buscar candidatos na janela atual
        def buscar_candidatos_na_janela(d_min_busca, d_max_busca):
            cands = []
            for cand_id in all_serials:
                if cand_id in state["locked"]:
                    continue
                cand_coord = state["coords"][cand_id]
                d = dist_tuple(curr_coord, cand_coord)
                if d_min_busca <= d <= d_max_busca:
                    # Calcular desvio angular
                    if prev_coord is not None:
                        ang = angle_tuple(prev_coord, curr_coord, cand_coord)
                        angular_dev = abs(ang - ideal_angle)
                    else:
                        ang = ideal_angle
                        angular_dev = 0.0
                    cands.append((cand_id, d, ang, angular_dev))
            return cands

        # 1. Tentar com limites originais
        candidatos = buscar_candidatos_na_janela(d_min_atual, d_max)

        # 2. Se não encontrou, expandir MAX até max_expansoes_por_direcao tentativas
        if not candidatos:
            for _ in range(max_expansoes_por_direcao):
                d_max += delta
                num_expansoes += 1
                candidatos = buscar_candidatos_na_janela(d_min_atual, d_max)
                if candidatos:
                    break

        # 3. Se ainda não encontrou, reduzir MIN até max_expansoes_por_direcao tentativas
        if not candidatos:
            for _ in range(max_expansoes_por_direcao):
                d_min_atual = max(0.0, d_min_atual - delta)
                num_expansoes += 1
                candidatos = buscar_candidatos_na_janela(d_min_atual, d_max)
                if candidatos:
                    break

        # 4. Se AINDA não encontrou, continuar expandindo ambos indefinidamente
        while not candidatos:
            d_max += delta
            d_min_atual = max(0.0, d_min_atual - delta)
            num_expansoes += 1
            candidatos = buscar_candidatos_na_janela(d_min_atual, d_max)

        # Atualizar d_min para o relatório de expansões
        d_min = d_min_atual

        # Ordenar candidatos por desvio angular (melhor primeiro para eficiência)
        candidatos.sort(key=lambda x: x[3])

        # Adicionar TODOS os candidatos à pilha (DFS explora todos)
        for cand_id, cand_dist, cand_ang, angular_dev in reversed(candidatos):
            # Criar novo estado com este candidato
            new_coords = state["coords"].copy()

            # Fazer troca se necessário
            if cand_id != next_id:
                tmp = new_coords[next_id]
                new_coords[next_id] = new_coords[cand_id]
                new_coords[cand_id] = tmp

            new_locked = set(state["locked"])
            new_locked.add(next_id)

            new_edges = list(state["edges"])
            new_edges.append({
                "from": curr_id,
                "to": next_id,
                "donor": cand_id,
                "distance": cand_dist,
                "angulo": cand_ang,
                "desvio_angular": angular_dev,
                "janela_expandida": num_expansoes > 0
            })

            new_expansoes = list(state["expansoes"])
            if num_expansoes > 0:
                new_expansoes.append({
                    "passo": k + 1,
                    "passo_total": max_steps,
                    "atomo_atual": {
                        "serial": curr_id,
                        "nome": curr_atom['name'],
                        "residuo": curr_atom['resname'],
                        "res_seq": curr_atom['resseq']
                    },
                    "atomo_proximo": {
                        "serial": next_id,
                        "nome": next_atom['name'],
                        "residuo": next_atom['resname'],
                        "res_seq": next_atom['resseq']
                    },
                    "janela_original": {"min": d_min_base, "max": d_max_base},
                    "janela_expandida": {"min": d_min, "max": d_max},
                    "expansoes": num_expansoes,
                    "min_reduzido": d_min < d_min_base,
                    "max_expandido": d_max > d_max_base,
                    "candidato_encontrado": {
                        "serial": cand_id,
                        "nome": atom_by_serial.get(cand_id, {}).get('name', '?'),
                        "residuo": atom_by_serial.get(cand_id, {}).get('resname', '?'),
                        "res_seq": atom_by_serial.get(cand_id, {}).get('resseq', '?'),
                        "distancia": cand_dist
                    },
                    "houve_troca": cand_id != next_id
                })

            new_state = {
                "step": k + 1,
                "coords": new_coords,
                "locked": new_locked,
                "total_dist": state["total_dist"] + cand_dist,
                "total_angular_deviation": state["total_angular_deviation"] + angular_dev,
                "edges": new_edges,
                "expansoes": new_expansoes
            }
            stack.append(new_state)

    # Selecionar a melhor solução (menor desvio angular total)
    if all_solutions:
        best_solution = min(all_solutions, key=lambda s: s["total_angular_deviation"])
        stats["expansoes_janela"] = best_solution["expansoes"]
        return [best_solution], stats

    return [], stats


def find_all_backbone_paths(atoms: List[Dict], atom_by_serial: Dict[int, Dict],
                            backbone_ids: List[int], coords: Dict[int, Tuple[float, float, float]],
                            max_steps: int = 23, min_total: float = 30.0,
                            max_total: float = 40.0,
                            delta: float = 0.02) -> Tuple[List[Dict], Dict]:
    """
    Wrapper que chama a busca DFS completa para o backbone.

    Lógica de expansão de janela:
    1. Primeiro tenta com limites originais
    2. Se não encontrou, expande o máximo em +delta até 15 tentativas
    3. Se ainda não encontrou, reduz o mínimo em -delta até 15 tentativas
    """
    solutions, stats = find_all_backbone_paths_dfs(
        atoms, atom_by_serial, backbone_ids, coords, max_steps, delta
    )

    # Adicionar estatísticas de distância total
    if solutions:
        total_dist = solutions[0]["total_dist"]
        if total_dist < min_total:
            stats["dist_status"] = "menor_que_minimo"
        elif total_dist <= max_total:
            stats["dist_status"] = "dentro_faixa"
        else:
            stats["dist_status"] = "maior_que_maximo"

    return solutions, stats


def formatar_relatorio_expansoes(expansoes: List[Dict], frame_num: int = None) -> str:
    """
    Formata um relatório detalhado de todas as expansões de janela realizadas.

    Args:
        expansoes: Lista de dicionários com informações de cada expansão
        frame_num: Número do frame (opcional, para relatórios multi-frame)

    Returns:
        String formatada com o relatório completo
    """
    if not expansoes:
        return ""

    lines = []
    lines.append("=" * 90)
    if frame_num is not None:
        lines.append(f"RELATÓRIO DE EXPANSÕES DE JANELA - FRAME {frame_num}")
    else:
        lines.append("RELATÓRIO DE EXPANSÕES DE JANELA DO BACKBONE")
    lines.append("=" * 90)
    lines.append(f"\nTotal de átomos que necessitaram expansão: {len(expansoes)}")
    lines.append("")

    for i, exp in enumerate(expansoes, 1):
        lines.append("-" * 90)
        lines.append(f"EXPANSÃO #{i}")
        lines.append("-" * 90)

        # Informações do passo
        lines.append(f"Passo: {exp['passo']}/{exp['passo_total']}")

        # Átomo atual
        curr = exp['atomo_atual']
        lines.append(f"\nÁtomo atual:     {curr['nome']:4s} (serial {curr['serial']}) - {curr['residuo']} {curr['res_seq']}")

        # Átomo próximo esperado
        prox = exp['atomo_proximo']
        lines.append(f"Próximo esperado: {prox['nome']:4s} (serial {prox['serial']}) - {prox['residuo']} {prox['res_seq']}")

        # Janela de distância
        jan_orig = exp['janela_original']
        jan_exp = exp['janela_expandida']
        lines.append(f"\nJanela original:  {jan_orig['min']:.3f} - {jan_orig['max']:.3f} Å")
        lines.append(f"Janela expandida: {jan_exp['min']:.3f} - {jan_exp['max']:.3f} Å")
        lines.append(f"Expansões necessárias: {exp['expansoes']}")

        # Indicar tipo de expansão
        min_reduzido = exp.get('min_reduzido', jan_exp['min'] < jan_orig['min'])
        max_expandido = exp.get('max_expandido', jan_exp['max'] > jan_orig['max'])

        if min_reduzido and max_expandido:
            delta_min = jan_orig['min'] - jan_exp['min']
            delta_max = jan_exp['max'] - jan_orig['max']
            lines.append(f"  - Mínimo reduzido em: -{delta_min:.3f} Å")
            lines.append(f"  - Máximo expandido em: +{delta_max:.3f} Å")
        elif min_reduzido:
            delta_min = jan_orig['min'] - jan_exp['min']
            lines.append(f"  - Mínimo reduzido em: -{delta_min:.3f} Å")
        elif max_expandido:
            delta_max = jan_exp['max'] - jan_orig['max']
            lines.append(f"  - Máximo expandido em: +{delta_max:.3f} Å")

        # Candidato encontrado
        cand = exp['candidato_encontrado']
        lines.append(f"\nCandidato encontrado:")
        lines.append(f"  Serial: {cand['serial']}")
        lines.append(f"  Nome: {cand['nome']}")
        lines.append(f"  Resíduo: {cand['residuo']} {cand['res_seq']}")
        lines.append(f"  Distância: {cand['distancia']:.4f} Å")

        # Houve troca?
        if exp['houve_troca']:
            lines.append(f"\n  ⚠ HOUVE TROCA DE COORDENADAS: {cand['serial']} ↔ {prox['serial']}")
        else:
            lines.append(f"\n  ✓ Candidato já estava na posição correta")

        lines.append("")

    lines.append("=" * 90)
    return "\n".join(lines)


def fase1_backbone(atoms: List[Dict], coords: Dict[int, Tuple[float, float, float]],
                   start_serial: int, end_serial: int, verbose: bool = True) -> Tuple[Dict[int, Tuple[float, float, float]], List[int], List[Dict]]:
    """
    Executa a Fase 1: Reconstrução do backbone.

    Returns:
        Tupla (coords_atualizadas, backbone_ids, expansoes_janela)
    """
    atom_by_serial = {a['serial']: a for a in atoms}
    backbone_ids = build_backbone_order(atoms, atom_by_serial, start_serial, end_serial)
    max_steps = len(backbone_ids) - 1
    solutions, stats = find_all_backbone_paths(atoms, atom_by_serial, backbone_ids, coords, max_steps, 30.0, 40.0)

    # Com janela dinâmica, sempre encontra solução
    solution = solutions[0]
    expansoes = stats.get("expansoes_janela", [])

    # Mostrar estatísticas de expansão se verbose
    if verbose and expansoes:
        print(f"  Backbone: {len(expansoes)} passos precisaram de expansão de janela")

    return solution["coords"], backbone_ids, expansoes


# ==============================================================================
# SEÇÃO 8: FASE 2 - ATRIBUIÇÃO DE HN, HA, O
# ==============================================================================

def fase2_hn_ha_o(atoms: List[Dict], coords: Dict[int, Tuple[float, float, float]],
                  backbone_ids: List[int], verbose: bool = True) -> Dict[int, Tuple[float, float, float]]:
    atom_by_serial = {a['serial']: a for a in atoms}

    # CORREÇÃO: Filtrar apenas átomos do escopo (576-583)
    start_atom = atom_by_serial[backbone_ids[0]]
    end_atom = atom_by_serial[backbone_ids[-1]]
    res_start = start_atom['resseq']
    res_end = end_atom['resseq']
    chain_alvo = start_atom['chain']
    all_serials = get_scope_serials(atoms, res_start, res_end, chain_alvo)

    backbone_set = set(backbone_ids)
    locked_sidechain = set()

    for serial_ref in backbone_ids:
        ref_atom = atom_by_serial[serial_ref]
        ref_name = ref_atom['name']
        chain = ref_atom['chain']
        resseq = ref_atom['resseq']

        if ref_name == "N":
            target_type = "HN"
        elif ref_name == "CA":
            target_type = "HA"
        elif ref_name == "C":
            target_type = "O"
        else:
            continue

        ligado_atom = None
        for a in atoms:
            if a['chain'] == chain and a['resseq'] == resseq and a['name'] == target_type:
                ligado_atom = a
                break
        if ligado_atom is None:
            continue

        ligado_serial = ligado_atom['serial']
        ref_coord = coords[serial_ref]
        best_serial = None
        best_dist = None

        for s in all_serials:
            if s in backbone_set or s in locked_sidechain:
                continue
            cand_coord = coords[s]
            d = dist_tuple(ref_coord, cand_coord)
            if best_dist is None or d < best_dist:
                best_dist = d
                best_serial = s

        if best_serial is None:
            continue

        if best_serial != ligado_serial:
            tmp = coords[ligado_serial]
            coords[ligado_serial] = coords[best_serial]
            coords[best_serial] = tmp

        locked_sidechain.add(ligado_serial)

    return coords


# ==============================================================================
# SEÇÃO 9: FASE 3 - CONEXÃO CA-CB
# ==============================================================================

def fase3_ca_cb(atoms: List[Dict], coords: Dict[int, Tuple[float, float, float]],
                backbone_ids: List[int], locked: Set[int], verbose: bool = True) -> Tuple[Dict[int, Tuple[float, float, float]], Set[int]]:
    """
    VERSÃO CORRIGIDA - Conexão CA-CB com expansão dinâmica de janela.
    
    CORREÇÕES:
    1. Expansão de janela quando não encontra candidatos
    2. Fallback para átomo mais próximo em caso extremo
    3. Nunca ignora um CB silenciosamente
    """
    atom_by_serial = {a['serial']: a for a in atoms}

    # Obter escopo
    start_atom = atom_by_serial[backbone_ids[0]]
    end_atom = atom_by_serial[backbone_ids[-1]]
    res_start = start_atom['resseq']
    res_end = end_atom['resseq']
    chain_alvo = start_atom['chain']
    all_serials = get_scope_serials(atoms, res_start, res_end, chain_alvo)

    residues = get_backbone_residues(backbone_ids, atom_by_serial)

    # Parâmetros de expansão de janela
    d_min_base = 1.40
    d_max_base = 1.60
    delta = 0.02
    max_expansoes = 30

    for res in residues:
        resname = res['resname']
        resseq = res['resseq']
        chain = res['chain']

        if resname == "GLY":
            continue

        N_id = res['N']
        CA_id = res['CA']
        C_id = res['C']

        cb_atom = None
        for a in atoms:
            if a['chain'] == chain and a['resseq'] == resseq and a['name'] == 'CB':
                cb_atom = a
                break
        if cb_atom is None:
            continue

        CB_id = cb_atom['serial']
        CA_coord = coords[CA_id]
        N_coord = coords[N_id]
        C_coord = coords[C_id]

        # CORREÇÃO: Busca com expansão dinâmica de janela
        d_min = d_min_base
        d_max = d_max_base
        candidates = []

        for iteration in range(max_expansoes):
            candidates = []
            for s in all_serials:
                if s in locked:
                    continue
                cand_coord = coords[s]
                d = dist_tuple(CA_coord, cand_coord)
                if d_min <= d <= d_max:
                    ca_dict = {'x': CA_coord[0], 'y': CA_coord[1], 'z': CA_coord[2]}
                    n_dict = {'x': N_coord[0], 'y': N_coord[1], 'z': N_coord[2]}
                    c_dict = {'x': C_coord[0], 'y': C_coord[1], 'z': C_coord[2]}
                    cand_dict = {'x': cand_coord[0], 'y': cand_coord[1], 'z': cand_coord[2]}
                    theta_N = angle(n_dict, ca_dict, cand_dict)
                    theta_C = angle(c_dict, ca_dict, cand_dict)
                    candidates.append((s, d, theta_N, theta_C))

            if candidates:
                break

            # Expandir janela
            d_min = max(0.5, d_min - delta)
            d_max = min(3.0, d_max + delta)

        # CORREÇÃO: Fallback quando não encontra candidatos
        if not candidates:
            if verbose:
                print(f"  AVISO: Fallback para CB do resíduo {resseq} ({resname})")

            all_distances = []
            for s in all_serials:
                if s in locked:
                    continue
                cand_coord = coords[s]
                d = dist_tuple(CA_coord, cand_coord)
                all_distances.append((s, d))

            if all_distances:
                all_distances.sort(key=lambda x: x[1])
                chosen_serial = all_distances[0][0]

                if chosen_serial != CB_id:
                    tmp = coords[CB_id]
                    coords[CB_id] = coords[chosen_serial]
                    coords[chosen_serial] = tmp

                locked.add(CB_id)
            else:
                locked.add(CB_id)
            continue

        # Seleção do melhor candidato
        if len(candidates) > 1:
            enriched = [(s, d, tN, tC) for (s, d, tN, tC) in candidates if tN is not None and tC is not None]
            if not enriched:
                candidates.sort(key=lambda x: x[1])
                chosen_serial = candidates[0][0]
            else:
                filtered = [(s, d, tN, tC) for (s, d, tN, tC) in enriched if 104.0 <= tN <= 116.0 and 105.0 <= tC <= 118.0]
                if filtered:
                    working = filtered
                else:
                    working = enriched
                def score(item):
                    _, _, tN, tC = item
                    return abs(tN - 110.0) + abs(tC - 111.0)
                working.sort(key=score)
                chosen_serial = working[0][0]
        else:
            chosen_serial = candidates[0][0]

        if chosen_serial != CB_id:
            tmp = coords[CB_id]
            coords[CB_id] = coords[chosen_serial]
            coords[chosen_serial] = tmp

        locked.add(CB_id)

    return coords, locked




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
    """Processa cadeia lateral de ILE com restrição de escopo."""
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
    """Processa cadeia lateral de THR com restrição de escopo."""
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
    """Processa cadeia lateral de LEU com restrição de escopo."""
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
    """
    Processa um único frame da trajetória.

    Retorna:
        Tupla (conteúdo_processado, sucesso, mensagem_erro, expansoes_backbone)
        - conteúdo_processado: String com frame (processado ou original se falhou)
        - sucesso: True se processou sem erros, False caso contrário
        - mensagem_erro: Descrição do erro (vazia se sucesso)
        - expansoes_backbone: Lista com informações de expansões de janela
    """
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

        # CORREÇÃO: Calcular scope_indices para restringir buscas aos resíduos 576-583
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
        if verbose:
            print(f"  AVISO Frame {frame_num}: {erro_msg}")
        # Retorna frame original em caso de erro (com expansões parciais se existirem)
        return "\n".join(lines), False, erro_msg, expansoes_backbone


# ==============================================================================
# SEÇÃO 13: PIPELINE PRINCIPAL MULTI-FRAME
# ==============================================================================

def executar_pipeline_multiframe(pdb_entrada: str, pdb_saida: str, start_serial: int = 8641,
                                  end_serial: int = 8776, chain: str = "B", verbose: bool = True) -> None:
    """
    Executa o pipeline completo para todos os frames de uma trajetória.
    """
    print("="*80)
    print("PIPELINE UNIFICADO DE REORDENAÇÃO CAR-T - MULTI-FRAME")
    print("="*80)
    print(f"Arquivo de entrada: {pdb_entrada}")
    print(f"Arquivo de saída: {pdb_saida}")
    print("="*80)

    # Ler arquivo completo
    with open(pdb_entrada, 'r') as f:
        content = f.read()

    # Dividir em frames
    frames = split_trajectory_into_frames(content)
    total_frames = len(frames)
    print(f"\nTotal de frames detectados: {total_frames}")

    if total_frames == 0:
        print("ERRO: Nenhum frame encontrado no arquivo!")
        return

    # Processar cada frame
    processed_frames = []
    frames_com_falha = []  # Lista de tuplas (numero_frame, mensagem_erro)
    frames_sucesso = 0
    todas_expansoes = []  # Lista de tuplas (numero_frame, lista_expansoes)

    for i, frame_content in enumerate(frames):
        frame_num = i + 1
        if verbose:
            print(f"\rProcessando frame {frame_num}/{total_frames}...", end="", flush=True)

        processed_frame, sucesso, erro_msg, expansoes = processar_frame(
            frame_content, frame_num, start_serial, end_serial, chain, verbose=False
        )
        processed_frames.append(processed_frame)

        # Registrar expansões (mesmo para frames com sucesso)
        if expansoes:
            todas_expansoes.append((frame_num, expansoes))

        if sucesso:
            frames_sucesso += 1
        else:
            frames_com_falha.append((frame_num, erro_msg))

    print(f"\rProcessamento concluído!                                    ")

    # Escrever arquivo de saída
    with open(pdb_saida, 'w') as f:
        f.write("\n".join(processed_frames))

    # Relatório final
    print("\n" + "="*80)
    print("RELATÓRIO DE PROCESSAMENTO")
    print("="*80)
    print(f"Total de frames:      {total_frames}")
    print(f"Frames com sucesso:   {frames_sucesso}")
    print(f"Frames com falha:     {len(frames_com_falha)}")

    if frames_com_falha:
        print("\n" + "-"*80)
        print("FRAMES COM FALHA (mantidos com coordenadas originais):")
        print("-"*80)

        # Gerar arquivo de diagnóstico detalhado
        diag_file = pdb_saida.replace(".pdb", "_diagnostico.txt")
        with open(diag_file, 'w') as f_diag:
            f_diag.write("="*80 + "\n")
            f_diag.write("RELATÓRIO DETALHADO DE DIAGNÓSTICO - FRAMES COM FALHA\n")
            f_diag.write("="*80 + "\n")
            f_diag.write(f"Arquivo de entrada: {pdb_entrada}\n")
            f_diag.write(f"Total de frames: {total_frames}\n")
            f_diag.write(f"Frames com falha: {len(frames_com_falha)}\n")
            f_diag.write("="*80 + "\n\n")

            for frame_num, erro in frames_com_falha:
                # Extrair resumo curto para console
                primeira_linha = erro.split('\n')[0] if '\n' in erro else erro
                print(f"  Frame {frame_num}: {primeira_linha[:100]}...")

                # Escrever diagnóstico completo no arquivo
                f_diag.write(f"\n{'#'*80}\n")
                f_diag.write(f"# FRAME {frame_num}\n")
                f_diag.write(f"{'#'*80}\n\n")
                f_diag.write(erro)
                f_diag.write("\n\n")

        print("\n" + "-"*80)
        print(f"Diagnóstico detalhado salvo em: {diag_file}")

    # Gerar relatório de expansões de janela
    if todas_expansoes:
        total_atomos_com_expansao = sum(len(exp) for _, exp in todas_expansoes)
        frames_com_expansao = len(todas_expansoes)

        print("\n" + "-"*80)
        print("ESTATÍSTICAS DE EXPANSÃO DE JANELA DO BACKBONE:")
        print("-"*80)
        print(f"Frames que precisaram de expansão: {frames_com_expansao}/{total_frames}")
        print(f"Total de átomos com expansão:      {total_atomos_com_expansao}")

        # Gerar arquivo de relatório de expansões
        expansao_file = pdb_saida.replace(".pdb", "_expansoes_backbone.txt")
        with open(expansao_file, 'w') as f_exp:
            f_exp.write("=" * 90 + "\n")
            f_exp.write("RELATÓRIO COMPLETO DE EXPANSÕES DE JANELA DO BACKBONE\n")
            f_exp.write("=" * 90 + "\n")
            f_exp.write(f"Arquivo de entrada: {pdb_entrada}\n")
            f_exp.write(f"Total de frames: {total_frames}\n")
            f_exp.write(f"Frames com expansão: {frames_com_expansao}\n")
            f_exp.write(f"Total de átomos com expansão: {total_atomos_com_expansao}\n")
            f_exp.write("=" * 90 + "\n\n")

            for frame_num, expansoes in todas_expansoes:
                # Escrever relatório formatado para cada frame
                relatorio = formatar_relatorio_expansoes(expansoes, frame_num)
                f_exp.write(relatorio)
                f_exp.write("\n\n")

        print(f"Relatório de expansões salvo em:   {expansao_file}")

    print("\n" + "="*80)
    print(f"Arquivo de saída: {pdb_saida}")
    print("="*80)
    print("PIPELINE CONCLUÍDO!")
    print("="*80)


# ==============================================================================
# SEÇÃO 14: INTERFACE DE LINHA DE COMANDO
# ==============================================================================

def main():
    """Função principal para execução via linha de comando."""
    try:
        from google.colab import files
        print("="*60)
        print("PIPELINE CAR-T MULTI-FRAME")
        print("="*60)
        print("\nSelecione o arquivo PDB/trajetória de entrada:")
        uploaded = files.upload()
        if not uploaded:
            raise RuntimeError("Nenhum arquivo foi enviado.")
        pdb_entrada = list(uploaded.keys())[0]
        pdb_saida = "trajetoria_cart_reordenada.pdb"
        executar_pipeline_multiframe(pdb_entrada, pdb_saida, verbose=True)
        files.download(pdb_saida)
    except ImportError:
        if len(sys.argv) < 2:
            print("Uso: python reordenacao-cart-v1.py <arquivo_entrada.pdb> [arquivo_saida.pdb]")
            print("\nExemplo:")
            print("  python reordenacao-cart-v1.py trajetoria.pdb trajetoria_reordenada.pdb")
            sys.exit(1)
        pdb_entrada = sys.argv[1]
        pdb_saida = sys.argv[2] if len(sys.argv) > 2 else "trajetoria_cart_reordenada.pdb"
        executar_pipeline_multiframe(pdb_entrada, pdb_saida, verbose=True)


if __name__ == "__main__":
    main()