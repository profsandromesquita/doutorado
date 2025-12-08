# -*- coding: utf-8 -*-
"""
================================================================================
PIPELINE UNIFICADO DE REORDENAÇÃO DE ESTRUTURA CAR-T
================================================================================

Este módulo consolida todos os 10 scripts de reordenação da estrutura CAR-T
em um único pipeline automatizado que processa internamente sem necessidade
de arquivos intermediários.

FLUXO DE PROCESSAMENTO:
1. GERANDO UM PDB COM BACKBONE REESTRUTURADO
2. BUSCA E REESTRUTURAÇÃO DOS ÁTOMOS HN, HA E C
3. BUSCA E REORDENAÇÃO DO SIDE CHAIN - CONEXAO CA-CB
4. REESTRUTURA 576 (ILE)
5. REESTRUTURA 577 (THR)
6. REORDENAÇÃO 578 (LEU)
7. REESTRUTURA 579 (TYR)
8. REESTRUTURAR 580 (CYS)
9. REESTRUTURAR 581 (LYS)
10. REESTRUTURA 582 (ARG)

AUTOR: Pipeline consolidado automaticamente
DATA: 2025
================================================================================
"""

import math
from typing import List, Dict, Set, Tuple, Optional, Any
from collections import defaultdict


# ==============================================================================
# SEÇÃO 1: FUNÇÕES GEOMÉTRICAS FUNDAMENTAIS
# ==============================================================================

def dist(a: Dict, b: Dict) -> float:
    """
    Calcula a distância euclidiana 3D entre dois átomos.

    Args:
        a: Dicionário do primeiro átomo com chaves 'x', 'y', 'z'
        b: Dicionário do segundo átomo com chaves 'x', 'y', 'z'

    Returns:
        Distância em Angstroms
    """
    return math.sqrt(
        (a['x'] - b['x'])**2 +
        (a['y'] - b['y'])**2 +
        (a['z'] - b['z'])**2
    )


def dist_tuple(p1: Tuple[float, float, float], p2: Tuple[float, float, float]) -> float:
    """
    Calcula a distância euclidiana 3D entre duas coordenadas em tupla.

    Args:
        p1: Tupla (x, y, z) do primeiro ponto
        p2: Tupla (x, y, z) do segundo ponto

    Returns:
        Distância em Angstroms
    """
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def angle(a: Dict, b: Dict, c: Dict) -> float:
    """
    Calcula o ângulo A-B-C em graus.

    O ângulo é medido no vértice B, formado pelos vetores BA e BC.

    Args:
        a: Dicionário do átomo A com chaves 'x', 'y', 'z'
        b: Dicionário do átomo B (vértice) com chaves 'x', 'y', 'z'
        c: Dicionário do átomo C com chaves 'x', 'y', 'z'

    Returns:
        Ângulo em graus (0-180)
    """
    v1 = (a['x'] - b['x'], a['y'] - b['y'], a['z'] - b['z'])
    v2 = (c['x'] - b['x'], c['y'] - b['y'], c['z'] - b['z'])

    n1 = math.sqrt(sum(v*v for v in v1))
    n2 = math.sqrt(sum(v*v for v in v2))

    if n1 < 1e-6 or n2 < 1e-6:
        return 0.0

    dot = sum(v1[i]*v2[i] for i in range(3))
    cosang = max(-1.0, min(1.0, dot/(n1*n2)))
    return math.degrees(math.acos(cosang))


def same_coords(a: Dict, b: Dict, tol: float = 1e-3) -> bool:
    """
    Verifica se dois átomos possuem as mesmas coordenadas dentro de uma tolerância.

    Args:
        a: Dicionário do primeiro átomo
        b: Dicionário do segundo átomo
        tol: Tolerância em Angstroms

    Returns:
        True se as coordenadas são iguais dentro da tolerância
    """
    return (
        abs(a['x'] - b['x']) < tol and
        abs(a['y'] - b['y']) < tol and
        abs(a['z'] - b['z']) < tol
    )


def swap_coords(a: Dict, b: Dict) -> None:
    """
    Troca as coordenadas entre dois átomos (in-place).

    Args:
        a: Dicionário do primeiro átomo
        b: Dicionário do segundo átomo
    """
    for coord in ('x', 'y', 'z'):
        a[coord], b[coord] = b[coord], a[coord]


# ==============================================================================
# SEÇÃO 2: FUNÇÕES DE PARSING E ESCRITA DE PDB
# ==============================================================================

def parse_pdb_content(content: str) -> Tuple[List[Dict], List[str]]:
    """
    Faz o parsing do conteúdo de um arquivo PDB.

    Args:
        content: Conteúdo completo do arquivo PDB como string

    Returns:
        Tupla contendo:
        - Lista de dicionários de átomos
        - Lista de linhas originais do arquivo
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
                # Linha mal formatada, ignorar
                continue

    return atoms, lines


def parse_pdb_file(filename: str) -> Tuple[List[Dict], List[str]]:
    """
    Lê e faz o parsing de um arquivo PDB.

    Args:
        filename: Caminho do arquivo PDB

    Returns:
        Tupla contendo:
        - Lista de dicionários de átomos
        - Lista de linhas originais do arquivo
    """
    with open(filename, 'r') as f:
        content = f.read()
    return parse_pdb_content(content)


def update_pdb_lines(lines: List[str], atoms: List[Dict]) -> List[str]:
    """
    Atualiza as coordenadas nas linhas do PDB com base nos átomos modificados.

    Args:
        lines: Lista de linhas originais do PDB
        atoms: Lista de dicionários de átomos com coordenadas atualizadas

    Returns:
        Lista de linhas atualizadas
    """
    for at in atoms:
        i = at['line_idx']
        line = lines[i]

        # Garante que a linha tenha pelo menos 54 colunas
        if len(line) < 54:
            line = line.ljust(54)

        # Atualiza colunas 31-54 (coordenadas X, Y, Z)
        new_line = (
            line[:30] +
            f"{at['x']:8.3f}{at['y']:8.3f}{at['z']:8.3f}" +
            line[54:]
        )
        lines[i] = new_line

    return lines


def write_pdb_file(filename: str, lines: List[str]) -> None:
    """
    Escreve um arquivo PDB a partir de uma lista de linhas.

    Args:
        filename: Caminho do arquivo de saída
        lines: Lista de linhas do PDB
    """
    with open(filename, 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"PDB escrito em: {filename}")


# ==============================================================================
# SEÇÃO 3: FUNÇÕES DE BUSCA E INDEXAÇÃO DE ÁTOMOS
# ==============================================================================

def find_residue_atoms(
    atoms: List[Dict],
    resname: str,
    resseq: int,
    chain: Optional[str] = None
) -> List[int]:
    """
    Encontra todos os índices de átomos de um resíduo específico.

    Args:
        atoms: Lista de dicionários de átomos
        resname: Nome do resíduo (ex: "ILE", "THR")
        resseq: Número do resíduo
        chain: Identificador da cadeia (opcional)

    Returns:
        Lista de índices dos átomos encontrados
    """
    idxs = []
    for i, a in enumerate(atoms):
        if a['resname'] == resname and a['resseq'] == resseq:
            if chain is None or a['chain'] == chain:
                idxs.append(i)
    return idxs


def find_atom_by_name(
    atoms: List[Dict],
    resseq: int,
    atom_name: str,
    chain: Optional[str] = None
) -> Optional[int]:
    """
    Encontra o índice de um átomo específico por nome e resíduo.

    Args:
        atoms: Lista de dicionários de átomos
        resseq: Número do resíduo
        atom_name: Nome do átomo (ex: "CA", "CB")
        chain: Identificador da cadeia (opcional)

    Returns:
        Índice do átomo ou None se não encontrado
    """
    for i, a in enumerate(atoms):
        if a['resseq'] == resseq and a['name'] == atom_name:
            if chain is None or a['chain'] == chain:
                return i
    return None


def build_name_index(atoms: List[Dict], residue_idxs: List[int]) -> Dict[str, int]:
    """
    Constrói um mapa de nome de átomo para índice dentro de um resíduo.

    Args:
        atoms: Lista de dicionários de átomos
        residue_idxs: Lista de índices de átomos do resíduo

    Returns:
        Dicionário {nome_átomo: índice}
    """
    name2idx = {}
    for i in residue_idxs:
        nm = atoms[i]['name']
        if nm not in name2idx:
            name2idx[nm] = i
    return name2idx


def build_serial_index(atoms: List[Dict]) -> Dict[int, Dict]:
    """
    Constrói um mapa de número serial para dicionário do átomo.

    Args:
        atoms: Lista de dicionários de átomos

    Returns:
        Dicionário {serial: átomo}
    """
    return {a['serial']: a for a in atoms}


# ==============================================================================
# SEÇÃO 4: CONSTRUÇÃO DO BACKBONE CANÔNICO
# ==============================================================================

def build_backbone_order(
    atoms: List[Dict],
    atom_by_serial: Dict[int, Dict],
    start_serial: int = 8641,
    end_serial: int = 8776
) -> List[int]:
    """
    Constrói a ordem canônica do backbone (N, CA, C) para os resíduos especificados.

    Args:
        atoms: Lista de dicionários de átomos
        atom_by_serial: Mapa serial -> átomo
        start_serial: Serial do N inicial (default: 8641 = N ILE 576)
        end_serial: Serial do C final (default: 8776 = C GLY 583)

    Returns:
        Lista de seriais do backbone na ordem [N576, CA576, C576, N577, ...]
    """
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
                if a["chain"] == chain
                and a["resseq"] == res
                and a["name"] == aname
            ]
            if not candidates:
                raise ValueError(
                    f"Não encontrei átomo {aname} no resíduo {res} cadeia {chain}."
                )
            backbone_ids.append(candidates[0]["serial"])

    # Checagens de consistência
    if backbone_ids[0] != start_serial:
        raise ValueError(
            f"Primeiro átomo do backbone ({backbone_ids[0]}) "
            f"não é o start_serial ({start_serial})."
        )
    if backbone_ids[-1] != end_serial:
        raise ValueError(
            f"Último átomo do backbone ({backbone_ids[-1]}) "
            f"não é o end_serial ({end_serial})."
        )

    return backbone_ids


def get_backbone_residues(
    backbone_ids: List[int],
    atom_by_serial: Dict[int, Dict]
) -> List[Dict]:
    """
    Constrói uma lista de resíduos do backbone com seus átomos N, CA, C.

    Args:
        backbone_ids: Lista de seriais do backbone
        atom_by_serial: Mapa serial -> átomo

    Returns:
        Lista de dicionários com informações de cada resíduo
    """
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
    Retorna os limites de distância para ligações do backbone.

    Args:
        name1: Nome do primeiro átomo
        name2: Nome do segundo átomo

    Returns:
        Tupla (d_min, d_max) em Angstroms ou None se par inválido
    """
    pair = (name1, name2)

    if pair in (("N", "CA"), ("CA", "N")):
        return (1.40, 1.60)
    if pair in (("CA", "C"), ("C", "CA")):
        return (1.40, 1.70)
    if pair in (("C", "N"), ("N", "C")):
        return (1.25, 1.45)

    return None


# Constantes de distância para diferentes tipos de ligação
BOND_LIMITS = {
    # Backbone
    'N-CA': (1.40, 1.60),
    'CA-C': (1.40, 1.70),
    'C-N': (1.25, 1.45),

    # Hidrogênios
    'C-H': (0.995, 1.115),       # Carbono sp3 - Hidrogênio
    'N-H_backbone': (0.90, 1.10), # Nitrogênio backbone - HN
    'N-H_amino': (0.90, 1.10),   # Nitrogênio amino - HZ
    'O-H': (0.85, 1.05),         # Oxigênio - Hidrogênio (OH)
    'S-H': (1.20, 1.45),         # Enxofre - Hidrogênio

    # Carbonos alifáticos (sp3)
    'CA-CB': (1.40, 1.60),
    'CB-CG': (1.40, 1.60),
    'CG-CD': (1.40, 1.60),
    'CD-CE': (1.40, 1.60),

    # Carbonos aromáticos (sp2)
    'C-C_arom': (1.30, 1.50),

    # Heteroátomos
    'CB-OG': (1.30, 1.50),       # THR: CB-OG1
    'CB-SG': (1.60, 1.90),       # CYS: CB-SG
    'CE-NZ': (1.38, 1.58),       # LYS: CE-NZ
    'CD-NE': (1.40, 1.60),       # ARG: CD-NE
    'NE-CZ': (1.25, 1.45),       # ARG: NE-CZ
    'CZ-NH': (1.25, 1.45),       # ARG: CZ-NH1/NH2
    'CZ-OH': (1.30, 1.50),       # TYR: CZ-OH
}


# ==============================================================================
# SEÇÃO 6: ALGORITMOS DE SELEÇÃO DE ÁTOMOS
# ==============================================================================

def escolher_vizinhos_dinamico(
    atoms: List[Dict],
    locked: Set[int],
    idx_ref: int,
    target_count: int,
    dmin_init: float,
    dmax_init: float,
    delta: float = 0.005,
    max_iter: int = 200,
    label: str = "H?"
) -> Tuple[List[int], float, float]:
    """
    Encontra exatamente N vizinhos usando janela dinâmica de distância.

    Se há mais candidatos que o necessário, estreita a janela.
    Se há menos candidatos que o necessário, alarga a janela.

    Args:
        atoms: Lista de dicionários de átomos
        locked: Conjunto de índices de átomos bloqueados
        idx_ref: Índice do átomo de referência
        target_count: Número de vizinhos desejados
        dmin_init: Distância mínima inicial
        dmax_init: Distância máxima inicial
        delta: Incremento/decremento por iteração
        max_iter: Número máximo de iterações
        label: Rótulo para mensagens de debug

    Returns:
        Tupla (lista de índices, dmin_final, dmax_final)
    """
    ref = atoms[idx_ref]
    dmin = dmin_init
    dmax = dmax_init

    for _ in range(max_iter):
        cand = []
        for i, a in enumerate(atoms):
            if i in locked or i == idx_ref:
                continue
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
        else:  # len(cand) < target_count
            dmin = max(0.0, dmin - delta)
            dmax += delta

    # Fallback: pega os mais próximos
    cand_all = []
    for i, a in enumerate(atoms):
        if i in locked or i == idx_ref:
            continue
        d = dist(ref, a)
        cand_all.append((i, d))

    cand_all.sort(key=lambda t: t[1])
    chosen = cand_all[:target_count]
    return [i for (i, _) in chosen], dmin, dmax


def atribuir_coord_alvos(
    atoms: List[Dict],
    locked: Set[int],
    idx_ref: int,
    alvo_idxs: List[int],
    dmin_init: float,
    dmax_init: float,
    target_count: int,
    delta: float = 0.005,
    max_iter: int = 200,
    label: str = "H?",
    verbose: bool = False
) -> Dict:
    """
    Atribui coordenadas de candidatos para átomos alvo usando janela dinâmica.

    Primeiro verifica se algum alvo já coincide com um candidato.
    Depois faz swap para os alvos restantes.

    Args:
        atoms: Lista de dicionários de átomos
        locked: Conjunto de índices de átomos bloqueados (será modificado)
        idx_ref: Índice do átomo de referência
        alvo_idxs: Lista de índices dos átomos alvo
        dmin_init: Distância mínima inicial
        dmax_init: Distância máxima inicial
        target_count: Número de candidatos a buscar
        delta: Incremento/decremento por iteração
        max_iter: Número máximo de iterações
        label: Rótulo para mensagens
        verbose: Se True, imprime detalhes

    Returns:
        Dicionário com detalhes da atribuição
    """
    ref = atoms[idx_ref]
    cand_idxs, dmin_final, dmax_final = escolher_vizinhos_dinamico(
        atoms, locked, idx_ref, target_count,
        dmin_init, dmax_init, delta, max_iter, label
    )

    cand_coords = {i: (atoms[i]['x'], atoms[i]['y'], atoms[i]['z']) for i in cand_idxs}
    used_cands = set()

    detalhes = {
        'ref': ref,
        'ref_idx': idx_ref,
        'janela_final': (dmin_final, dmax_final),
        'mapeamentos': []
    }

    # 1) Alvos que já coincidem com algum candidato
    for idx_alvo in alvo_idxs:
        alvo = atoms[idx_alvo]
        match = None
        for i_cand in cand_idxs:
            if i_cand in used_cands:
                continue
            cx, cy, cz = cand_coords[i_cand]
            if (abs(alvo['x'] - cx) < 1e-3 and
                abs(alvo['y'] - cy) < 1e-3 and
                abs(alvo['z'] - cz) < 1e-3):
                match = i_cand
                break

        if match is not None:
            used_cands.add(match)
            locked.add(idx_alvo)
            d = dist(ref, alvo)
            detalhes['mapeamentos'].append({
                'alvo_idx': idx_alvo,
                'alvo': alvo,
                'cand_idx': match,
                'cand': atoms[match],
                'dist_ref_alvo': d,
                'swap_feito': False
            })

    # 2) Alvos restantes com swap
    for idx_alvo in alvo_idxs:
        if idx_alvo in locked:
            continue

        alvo = atoms[idx_alvo]
        cand_rest = [i for i in cand_idxs if i not in used_cands]

        if not cand_rest:
            locked.add(idx_alvo)
            detalhes['mapeamentos'].append({
                'alvo_idx': idx_alvo,
                'alvo': alvo,
                'cand_idx': None,
                'cand': None,
                'dist_ref_alvo': dist(ref, alvo),
                'swap_feito': False
            })
            continue

        # Escolhe candidato mais próximo do ref
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

        detalhes['mapeamentos'].append({
            'alvo_idx': idx_alvo,
            'alvo': alvo,
            'cand_idx': i_cand,
            'cand': cand_atom,
            'dist_ref_alvo_antes': d_before,
            'dist_ref_alvo_depois': d_after,
            'swap_feito': swap_feito
        })

    return detalhes


def escolher_pesado_com_angulo(
    atoms: List[Dict],
    locked: Set[int],
    idx_center: int,
    idx_prev: int,
    dmin: float,
    dmax: float,
    ang1_min: float = 101.0,
    ang1_max: float = 119.0,
    ang2_min: float = 104.0,
    ang2_max: float = 116.0,
    alvo_angulo: float = 110.0,
    delta: float = 0.005,
    max_iter: int = 200,
    label: str = "PESO"
) -> Tuple[int, float, float]:
    """
    Seleciona um átomo pesado (C, N, O, S) usando critérios de distância e ângulo.

    Usa refinamento em 3 etapas:
    - REF1: Filtra por ângulo amplo [ang1_min, ang1_max]
    - REF2: Filtra por ângulo estreito [ang2_min, ang2_max]
    - REF3: Seleciona mais próximo do ângulo alvo

    Se REF1 ou REF2 eliminam todos candidatos vindos de >=2, volta à lista anterior.

    Args:
        atoms: Lista de dicionários de átomos
        locked: Conjunto de índices de átomos bloqueados
        idx_center: Índice do átomo central (de onde mede distância)
        idx_prev: Índice do átomo anterior (para medir ângulo)
        dmin: Distância mínima
        dmax: Distância máxima
        ang1_min, ang1_max: Limites REF1
        ang2_min, ang2_max: Limites REF2
        alvo_angulo: Ângulo alvo para REF3
        delta: Incremento para expansão de janela
        max_iter: Iterações máximas para expansão
        label: Rótulo para mensagens

    Returns:
        Tupla (índice do candidato escolhido, distância, ângulo)
    """
    center = atoms[idx_center]
    prev = atoms[idx_prev]
    dmin_curr = dmin
    dmax_curr = dmax

    # Busca candidatos, expandindo janela se necessário
    for _ in range(max_iter):
        base = []
        for i, a in enumerate(atoms):
            if i in locked or i == idx_center:
                continue
            d = dist(center, a)
            if dmin_curr <= d <= dmax_curr:
                ang = angle(prev, center, a)
                base.append((i, d, ang))
        if base:
            break
        dmin_curr = max(0.0, dmin_curr - delta)
        dmax_curr += delta

    if not base:
        raise RuntimeError(f"Nenhum candidato encontrado para {label} em {dmin:.2f}-{dmax:.2f} Å.")

    if len(base) == 1:
        return base[0]

    # REF1
    prev_list = base
    prev_count = len(prev_list)
    ref1 = [(i, d, ang) for (i, d, ang) in prev_list if ang1_min <= ang <= ang1_max]

    if len(ref1) == 0 and prev_count >= 2:
        ref3 = sorted(prev_list, key=lambda t: abs(t[2] - alvo_angulo))
        return ref3[0]
    if len(ref1) == 1:
        return ref1[0]

    # REF2
    if ref1:
        prev_list = ref1
        prev_count = len(prev_list)
    ref2 = [(i, d, ang) for (i, d, ang) in prev_list if ang2_min <= ang <= ang2_max]

    if len(ref2) == 0 and prev_count >= 2:
        ref3 = sorted(prev_list, key=lambda t: abs(t[2] - alvo_angulo))
        return ref3[0]
    if len(ref2) == 1:
        return ref2[0]

    # REF3
    if ref2:
        ref3 = sorted(ref2, key=lambda t: abs(t[2] - alvo_angulo))
    else:
        ref3 = sorted(prev_list, key=lambda t: abs(t[2] - alvo_angulo))

    return ref3[0]


# ==============================================================================
# SEÇÃO 7: FASE 1 - RECONSTRUÇÃO DO BACKBONE
# ==============================================================================

def fase1_backbone(
    atoms: List[Dict],
    coords: Dict[int, Tuple[float, float, float]],
    start_serial: int = 8641,
    end_serial: int = 8776,
    verbose: bool = True
) -> Tuple[Dict[int, Tuple[float, float, float]], List[int]]:
    """
    FASE 1: Reconstrução do backbone N-CA-C.

    Explora caminhos possíveis e encontra uma solução válida para o backbone.

    Args:
        atoms: Lista de dicionários de átomos
        coords: Dicionário serial -> (x, y, z) de coordenadas
        start_serial: Serial do N inicial
        end_serial: Serial do C final
        verbose: Se True, imprime progresso

    Returns:
        Tupla (coordenadas atualizadas, lista de seriais do backbone)
    """
    if verbose:
        print("\n" + "="*80)
        print("FASE 1: RECONSTRUÇÃO DO BACKBONE")
        print("="*80)

    atom_by_serial = {a['serial']: a for a in atoms}
    backbone_ids = build_backbone_order(atoms, atom_by_serial, start_serial, end_serial)

    if verbose:
        print(f"Backbone canônico construído: {len(backbone_ids)} átomos")
        print(f"Resíduos: {backbone_ids[0]} (N) até {backbone_ids[-1]} (C)")

    # Algoritmo de busca de caminho válido (DFS simplificado - pega primeiro válido)
    all_serials = [a['serial'] for a in atoms]
    max_steps = len(backbone_ids) - 1

    # Estado inicial
    start_id = backbone_ids[0]
    locked = {start_id}

    # Processa cada passo do backbone
    for k in range(max_steps):
        curr_id = backbone_ids[k]
        next_id = backbone_ids[k + 1]

        curr_atom = atom_by_serial[curr_id]
        next_atom = atom_by_serial[next_id]

        limits = bond_limits_backbone(curr_atom['name'], next_atom['name'])
        if limits is None:
            raise ValueError(f"Par de backbone inesperado: {curr_atom['name']}-{next_atom['name']}")

        d_min, d_max = limits
        curr_coord = coords[curr_id]

        # Busca melhor candidato
        best_serial = None
        best_dist = None

        for cand_id in all_serials:
            if cand_id in locked:
                continue

            cand_coord = coords[cand_id]
            d = dist_tuple(curr_coord, cand_coord)

            if d_min <= d <= d_max:
                if best_dist is None or d < best_dist:
                    best_dist = d
                    best_serial = cand_id

        if best_serial is None:
            raise RuntimeError(f"Sem candidato para passo {k}: {curr_atom['name']} -> {next_atom['name']}")

        # Swap de coordenadas se necessário
        if best_serial != next_id:
            tmp = coords[next_id]
            coords[next_id] = coords[best_serial]
            coords[best_serial] = tmp

        locked.add(next_id)

        if verbose and (k + 1) % 5 == 0:
            print(f"  Passo {k+1}/{max_steps} concluído")

    if verbose:
        print(f"Backbone reconstruído com sucesso!")

    return coords, backbone_ids


# ==============================================================================
# SEÇÃO 8: FASE 2 - ATRIBUIÇÃO DE HN, HA, O
# ==============================================================================

def fase2_hn_ha_o(
    atoms: List[Dict],
    coords: Dict[int, Tuple[float, float, float]],
    backbone_ids: List[int],
    verbose: bool = True
) -> Dict[int, Tuple[float, float, float]]:
    """
    FASE 2: Atribuição dos átomos HN, HA e O do backbone.

    Para cada N do backbone, encontra o HN mais próximo.
    Para cada CA do backbone, encontra o HA mais próximo.
    Para cada C do backbone, encontra o O mais próximo.

    Args:
        atoms: Lista de dicionários de átomos
        coords: Dicionário serial -> (x, y, z) de coordenadas
        backbone_ids: Lista de seriais do backbone
        verbose: Se True, imprime progresso

    Returns:
        Coordenadas atualizadas
    """
    if verbose:
        print("\n" + "="*80)
        print("FASE 2: ATRIBUIÇÃO DE HN, HA E O")
        print("="*80)

    atom_by_serial = {a['serial']: a for a in atoms}
    all_serials = [a['serial'] for a in atoms]

    backbone_set = set(backbone_ids)
    locked_sidechain = set()

    count = 0
    for serial_ref in backbone_ids:
        ref_atom = atom_by_serial[serial_ref]
        ref_name = ref_atom['name']
        chain = ref_atom['chain']
        resseq = ref_atom['resseq']

        # Decide tipo alvo
        if ref_name == "N":
            target_type = "HN"
        elif ref_name == "CA":
            target_type = "HA"
        elif ref_name == "C":
            target_type = "O"
        else:
            continue

        # Encontra o átomo alvo no resíduo
        ligado_atom = None
        for a in atoms:
            if a['chain'] == chain and a['resseq'] == resseq and a['name'] == target_type:
                ligado_atom = a
                break

        if ligado_atom is None:
            continue

        ligado_serial = ligado_atom['serial']
        ref_coord = coords[serial_ref]

        # Busca candidato mais próximo
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

        # Swap se necessário
        if best_serial != ligado_serial:
            tmp = coords[ligado_serial]
            coords[ligado_serial] = coords[best_serial]
            coords[best_serial] = tmp

        locked_sidechain.add(ligado_serial)
        count += 1

    if verbose:
        print(f"Atribuídos {count} átomos HN/HA/O")

    return coords


# ==============================================================================
# SEÇÃO 9: FASE 3 - CONEXÃO CA-CB
# ==============================================================================

def fase3_ca_cb(
    atoms: List[Dict],
    coords: Dict[int, Tuple[float, float, float]],
    backbone_ids: List[int],
    locked: Set[int],
    verbose: bool = True
) -> Tuple[Dict[int, Tuple[float, float, float]], Set[int]]:
    """
    FASE 3: Estabelece a conexão CA-CB para todos os resíduos.

    Args:
        atoms: Lista de dicionários de átomos
        coords: Dicionário serial -> (x, y, z) de coordenadas
        backbone_ids: Lista de seriais do backbone
        locked: Conjunto de seriais já bloqueados
        verbose: Se True, imprime progresso

    Returns:
        Tupla (coordenadas atualizadas, conjunto locked atualizado)
    """
    if verbose:
        print("\n" + "="*80)
        print("FASE 3: CONEXÃO CA-CB")
        print("="*80)

    atom_by_serial = {a['serial']: a for a in atoms}
    all_serials = [a['serial'] for a in atoms]
    residues = get_backbone_residues(backbone_ids, atom_by_serial)

    count = 0
    for res in residues:
        resname = res['resname']
        resseq = res['resseq']
        chain = res['chain']

        # Glicina não tem CB
        if resname == "GLY":
            continue

        N_id = res['N']
        CA_id = res['CA']
        C_id = res['C']

        # Encontra CB
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

        # Busca candidatos por distância
        candidates = []
        for s in all_serials:
            if s in locked:
                continue

            cand_coord = coords[s]
            d = dist_tuple(CA_coord, cand_coord)

            if 1.40 <= d <= 1.60:
                # Calcula ângulo N-CA-candidato
                # Precisa converter coords para dict temporário
                ca_dict = {'x': CA_coord[0], 'y': CA_coord[1], 'z': CA_coord[2]}
                n_dict = {'x': N_coord[0], 'y': N_coord[1], 'z': N_coord[2]}
                cand_dict = {'x': cand_coord[0], 'y': cand_coord[1], 'z': cand_coord[2]}

                theta_N = angle(n_dict, ca_dict, cand_dict)
                candidates.append((s, d, theta_N))

        if not candidates:
            continue

        # Refinamento por ângulo
        if len(candidates) > 1:
            # REF1: 104-116
            filtered = [(s, d, t) for (s, d, t) in candidates if 104.0 <= t <= 116.0]
            if filtered:
                candidates = filtered

            # REF3: mais próximo de 110
            candidates.sort(key=lambda x: abs(x[2] - 110.0))

        chosen_serial = candidates[0][0]

        # Swap se necessário
        if chosen_serial != CB_id:
            tmp = coords[CB_id]
            coords[CB_id] = coords[chosen_serial]
            coords[chosen_serial] = tmp

        locked.add(CB_id)
        count += 1

    if verbose:
        print(f"Conectados {count} átomos CB")

    return coords, locked


# ==============================================================================
# SEÇÃO 10: FUNÇÕES AUXILIARES PARA CADEIAS LATERAIS
# ==============================================================================

def aplicar_coords_para_atoms(atoms: List[Dict], coords: Dict[int, Tuple[float, float, float]]) -> None:
    """
    Aplica as coordenadas do dicionário coords de volta aos átomos.

    Args:
        atoms: Lista de dicionários de átomos (será modificada in-place)
        coords: Dicionário serial -> (x, y, z)
    """
    serial_to_idx = {a['serial']: i for i, a in enumerate(atoms)}

    for serial, (x, y, z) in coords.items():
        if serial in serial_to_idx:
            idx = serial_to_idx[serial]
            atoms[idx]['x'] = x
            atoms[idx]['y'] = y
            atoms[idx]['z'] = z


def extrair_coords_de_atoms(atoms: List[Dict]) -> Dict[int, Tuple[float, float, float]]:
    """
    Extrai coordenadas dos átomos para um dicionário.

    Args:
        atoms: Lista de dicionários de átomos

    Returns:
        Dicionário serial -> (x, y, z)
    """
    return {a['serial']: (a['x'], a['y'], a['z']) for a in atoms}


def bloquear_residuo(
    atoms: List[Dict],
    locked: Set[int],
    resname: str,
    resseq: int,
    chain: str
) -> None:
    """
    Bloqueia todos os átomos de um resíduo no conjunto locked.

    Args:
        atoms: Lista de dicionários de átomos
        locked: Conjunto de índices bloqueados (será modificado)
        resname: Nome do resíduo
        resseq: Número do resíduo
        chain: Identificador da cadeia
    """
    for i, a in enumerate(atoms):
        if a['resname'] == resname and a['resseq'] == resseq and a['chain'] == chain:
            locked.add(i)


def inicializar_locked_base(atoms: List[Dict]) -> Set[int]:
    """
    Inicializa o conjunto locked com átomos base (N, CA, C, HN, HA, O, CB).

    Args:
        atoms: Lista de dicionários de átomos

    Returns:
        Conjunto de índices bloqueados
    """
    locked = set()
    base_names = {"N", "CA", "C", "HN", "HA", "O", "CB"}

    for i, a in enumerate(atoms):
        if a['name'] in base_names:
            locked.add(i)

    return locked


# ==============================================================================
# SEÇÃO 11: PROCESSAMENTO DE CADEIAS LATERAIS POR RESÍDUO
# ==============================================================================

def processar_ile(
    atoms: List[Dict],
    locked: Set[int],
    resseq: int = 576,
    chain: str = "B",
    verbose: bool = True
) -> None:
    """
    Processa a cadeia lateral de Isoleucina (ILE).

    Estrutura: CB -> CG1 -> CD1
                  -> CG2 (metil)

    Args:
        atoms: Lista de dicionários de átomos (modificada in-place)
        locked: Conjunto de índices bloqueados (modificado in-place)
        resseq: Número do resíduo
        chain: Identificador da cadeia
        verbose: Se True, imprime progresso
    """
    if verbose:
        print(f"\n  Processando ILE {chain} {resseq}...")

    res_idxs = find_residue_atoms(atoms, "ILE", resseq, chain)
    if not res_idxs:
        if verbose:
            print(f"    AVISO: ILE {chain} {resseq} não encontrado")
        return

    name2idx = build_name_index(atoms, res_idxs)

    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_CG1 = name2idx.get("CG1")
    idx_CG2 = name2idx.get("CG2")
    idx_CD1 = name2idx.get("CD1")

    # CB -> CG2
    if idx_CG2 is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.40, dmax=1.65,
                ang1_min=101.0, ang1_max=119.0,
                ang2_min=104.0, ang2_max=116.0,
                alvo_angulo=110.0, label="CG2"
            )
            if not same_coords(atoms[idx_CG2], atoms[cand_idx]):
                swap_coords(atoms[idx_CG2], atoms[cand_idx])
            locked.add(idx_CG2)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CG2: {e}")

    # HG2 (3 hidrogênios)
    alvo_hg2 = [name2idx.get(n) for n in ("1HG2", "2HG2", "3HG2") if name2idx.get(n) is not None]
    if alvo_hg2 and idx_CG2 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CG2, alvo_hg2,
            dmin_init=0.990, dmax_init=1.150,
            target_count=len(alvo_hg2), label="HG2"
        )

    # CB -> CG1
    if idx_CG1 is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.40, dmax=1.65,
                ang1_min=101.0, ang1_max=119.0,
                ang2_min=104.0, ang2_max=116.0,
                alvo_angulo=110.0, label="CG1"
            )
            if not same_coords(atoms[idx_CG1], atoms[cand_idx]):
                swap_coords(atoms[idx_CG1], atoms[cand_idx])
            locked.add(idx_CG1)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CG1: {e}")

    # HG1 e CD1 (vizinhos de CG1)
    idx_1HG1 = name2idx.get("1HG1")
    idx_2HG1 = name2idx.get("2HG1")
    alvo_cg1 = [x for x in [idx_1HG1, idx_2HG1, idx_CD1] if x is not None]
    if alvo_cg1 and idx_CG1 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CG1, alvo_cg1,
            dmin_init=0.990, dmax_init=1.150,
            target_count=len(alvo_cg1), label="CG1_viz"
        )

    # HD (3 hidrogênios de CD1)
    alvo_hd = [name2idx.get(n) for n in ("HD1", "HD2", "HD3") if name2idx.get(n) is not None]
    if not alvo_hd:
        alvo_hd = [name2idx.get(n) for n in ("1HD1", "2HD1", "3HD1") if name2idx.get(n) is not None]
    if alvo_hd and idx_CD1 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CD1, alvo_hd,
            dmin_init=0.990, dmax_init=1.150,
            target_count=len(alvo_hd), label="HD"
        )

    if verbose:
        print(f"    ILE {chain} {resseq} processado")


def processar_thr(
    atoms: List[Dict],
    locked: Set[int],
    resseq: int = 577,
    chain: str = "B",
    verbose: bool = True
) -> None:
    """
    Processa a cadeia lateral de Treonina (THR).

    Estrutura: CB -> OG1 -> HG1
                  -> CG2 (metil)
    """
    if verbose:
        print(f"\n  Processando THR {chain} {resseq}...")

    res_idxs = find_residue_atoms(atoms, "THR", resseq, chain)
    if not res_idxs:
        if verbose:
            print(f"    AVISO: THR {chain} {resseq} não encontrado")
        return

    name2idx = build_name_index(atoms, res_idxs)

    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_OG1 = name2idx.get("OG1")
    idx_HG1 = name2idx.get("HG1")
    idx_HB = name2idx.get("HB")
    idx_CG2 = name2idx.get("CG2")

    # CB -> OG1
    if idx_OG1 is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.30, dmax=1.50,
                ang1_min=101.0, ang1_max=119.0,
                ang2_min=104.0, ang2_max=116.0,
                alvo_angulo=110.0, label="OG1"
            )
            if not same_coords(atoms[idx_OG1], atoms[cand_idx]):
                swap_coords(atoms[idx_OG1], atoms[cand_idx])
            locked.add(idx_OG1)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO OG1: {e}")

    # OG1 -> HG1
    if idx_HG1 is not None and idx_OG1 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_OG1, [idx_HG1],
            dmin_init=0.90, dmax_init=1.10,
            target_count=1, label="HG1"
        )

    # CB -> HB
    if idx_HB is not None and idx_CB is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CB, [idx_HB],
            dmin_init=1.00, dmax_init=1.20,
            target_count=1, label="HB"
        )

    # CB -> CG2
    if idx_CG2 is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0,
                ang2_min=104.0, ang2_max=116.0,
                alvo_angulo=110.0, label="CG2"
            )
            if not same_coords(atoms[idx_CG2], atoms[cand_idx]):
                swap_coords(atoms[idx_CG2], atoms[cand_idx])
            locked.add(idx_CG2)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CG2: {e}")

    # CG2 -> HG2 (3 hidrogênios)
    alvo_hg2 = [name2idx.get(n) for n in ("1HG2", "2HG2", "3HG2") if name2idx.get(n) is not None]
    if alvo_hg2 and idx_CG2 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CG2, alvo_hg2,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hg2), label="HG2"
        )

    if verbose:
        print(f"    THR {chain} {resseq} processado")


def processar_leu(
    atoms: List[Dict],
    locked: Set[int],
    resseq: int = 578,
    chain: str = "B",
    verbose: bool = True
) -> None:
    """
    Processa a cadeia lateral de Leucina (LEU).

    Estrutura: CB -> CG -> CD1 (metil)
                       -> CD2 (metil)
    """
    if verbose:
        print(f"\n  Processando LEU {chain} {resseq}...")

    res_idxs = find_residue_atoms(atoms, "LEU", resseq, chain)
    if not res_idxs:
        if verbose:
            print(f"    AVISO: LEU {chain} {resseq} não encontrado")
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

    # CB -> HB1/HB2
    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CB, alvo_hb,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hb), label="HB"
        )

    # CB -> CG
    if idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0,
                ang2_min=104.0, ang2_max=116.0,
                alvo_angulo=110.0, label="CG"
            )
            if not same_coords(atoms[idx_CG], atoms[cand_idx]):
                swap_coords(atoms[idx_CG], atoms[cand_idx])
            locked.add(idx_CG)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CG: {e}")

    # CG -> HG
    if idx_HG is not None and idx_CG is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CG, [idx_HG],
            dmin_init=0.995, dmax_init=1.115,
            target_count=1, label="HG"
        )

    # CG -> CD1
    if idx_CD1 is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CG, idx_CB,
                dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0,
                ang2_min=104.0, ang2_max=116.0,
                alvo_angulo=110.0, label="CD1"
            )
            if not same_coords(atoms[idx_CD1], atoms[cand_idx]):
                swap_coords(atoms[idx_CD1], atoms[cand_idx])
            locked.add(idx_CD1)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CD1: {e}")

    # CG -> CD2
    if idx_CD2 is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CG, idx_CB,
                dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0,
                ang2_min=104.0, ang2_max=116.0,
                alvo_angulo=110.0, label="CD2"
            )
            if not same_coords(atoms[idx_CD2], atoms[cand_idx]):
                swap_coords(atoms[idx_CD2], atoms[cand_idx])
            locked.add(idx_CD2)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CD2: {e}")

    # CD1 -> HD1 (3 hidrogênios)
    alvo_hd1 = [name2idx.get(n) for n in ("1HD1", "2HD1", "3HD1") if name2idx.get(n) is not None]
    if alvo_hd1 and idx_CD1 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CD1, alvo_hd1,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hd1), label="HD1"
        )

    # CD2 -> HD2 (3 hidrogênios)
    alvo_hd2 = [name2idx.get(n) for n in ("1HD2", "2HD2", "3HD2") if name2idx.get(n) is not None]
    if alvo_hd2 and idx_CD2 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CD2, alvo_hd2,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hd2), label="HD2"
        )

    if verbose:
        print(f"    LEU {chain} {resseq} processado")


def processar_tyr(
    atoms: List[Dict],
    locked: Set[int],
    resseq: int = 579,
    chain: str = "B",
    verbose: bool = True
) -> None:
    """
    Processa a cadeia lateral de Tirosina (TYR).

    Estrutura aromática: CB -> CG -> CD1 -> CE1 -> CZ -> OH
                                  -> CD2 -> CE2 ->|
    """
    if verbose:
        print(f"\n  Processando TYR {chain} {resseq}...")

    res_idxs = find_residue_atoms(atoms, "TYR", resseq, chain)
    if not res_idxs:
        if verbose:
            print(f"    AVISO: TYR {chain} {resseq} não encontrado")
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

    # HB1/HB2
    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")
    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CB, alvo_hb,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hb), label="HB"
        )

    # CB -> CG (alifático)
    if idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.40, dmax=1.60,
                ang1_min=101.0, ang1_max=119.0,
                ang2_min=104.0, ang2_max=116.0,
                alvo_angulo=110.0, label="CG"
            )
            if not same_coords(atoms[idx_CG], atoms[cand_idx]):
                swap_coords(atoms[idx_CG], atoms[cand_idx])
            locked.add(idx_CG)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CG: {e}")

    # CG -> CD1 (aromático)
    if idx_CD1 is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CG, idx_CB,
                dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0,
                ang2_min=117.0, ang2_max=123.0,
                alvo_angulo=120.0, label="CD1"
            )
            if not same_coords(atoms[idx_CD1], atoms[cand_idx]):
                swap_coords(atoms[idx_CD1], atoms[cand_idx])
            locked.add(idx_CD1)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CD1: {e}")

    # HD1
    idx_HD1 = name2idx.get("HD1")
    if idx_HD1 is not None and idx_CD1 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CD1, [idx_HD1],
            dmin_init=0.995, dmax_init=1.115,
            target_count=1, label="HD1"
        )

    # CG -> CD2 (aromático)
    if idx_CD2 is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CG, idx_CB,
                dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0,
                ang2_min=117.0, ang2_max=123.0,
                alvo_angulo=120.0, label="CD2"
            )
            if not same_coords(atoms[idx_CD2], atoms[cand_idx]):
                swap_coords(atoms[idx_CD2], atoms[cand_idx])
            locked.add(idx_CD2)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CD2: {e}")

    # HD2
    idx_HD2 = name2idx.get("HD2")
    if idx_HD2 is not None and idx_CD2 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CD2, [idx_HD2],
            dmin_init=0.995, dmax_init=1.115,
            target_count=1, label="HD2"
        )

    # CD1 -> CE1
    if idx_CE1 is not None and idx_CD1 is not None and idx_CG is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CD1, idx_CG,
                dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0,
                ang2_min=117.0, ang2_max=123.0,
                alvo_angulo=120.0, label="CE1"
            )
            if not same_coords(atoms[idx_CE1], atoms[cand_idx]):
                swap_coords(atoms[idx_CE1], atoms[cand_idx])
            locked.add(idx_CE1)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CE1: {e}")

    # HE1
    idx_HE1 = name2idx.get("HE1")
    if idx_HE1 is not None and idx_CE1 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CE1, [idx_HE1],
            dmin_init=0.995, dmax_init=1.115,
            target_count=1, label="HE1"
        )

    # CD2 -> CE2
    if idx_CE2 is not None and idx_CD2 is not None and idx_CG is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CD2, idx_CG,
                dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0,
                ang2_min=117.0, ang2_max=123.0,
                alvo_angulo=120.0, label="CE2"
            )
            if not same_coords(atoms[idx_CE2], atoms[cand_idx]):
                swap_coords(atoms[idx_CE2], atoms[cand_idx])
            locked.add(idx_CE2)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CE2: {e}")

    # HE2
    idx_HE2 = name2idx.get("HE2")
    if idx_HE2 is not None and idx_CE2 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CE2, [idx_HE2],
            dmin_init=0.995, dmax_init=1.115,
            target_count=1, label="HE2"
        )

    # CZ (entre CE1 e CE2)
    if idx_CZ is not None and idx_CE1 is not None and idx_CD1 is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CE1, idx_CD1,
                dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0,
                ang2_min=117.0, ang2_max=123.0,
                alvo_angulo=120.0, label="CZ"
            )
            if not same_coords(atoms[idx_CZ], atoms[cand_idx]):
                swap_coords(atoms[idx_CZ], atoms[cand_idx])
            locked.add(idx_CZ)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CZ: {e}")

    # CZ -> OH
    if idx_OH is not None and idx_CZ is not None and idx_CE1 is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CZ, idx_CE1,
                dmin=1.30, dmax=1.50,
                ang1_min=114.0, ang1_max=126.0,
                ang2_min=117.0, ang2_max=123.0,
                alvo_angulo=120.0, label="OH"
            )
            if not same_coords(atoms[idx_OH], atoms[cand_idx]):
                swap_coords(atoms[idx_OH], atoms[cand_idx])
            locked.add(idx_OH)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO OH: {e}")

    # OH -> HH
    if idx_HH is not None and idx_OH is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_OH, [idx_HH],
            dmin_init=0.85, dmax_init=1.05,
            target_count=1, label="HH"
        )

    if verbose:
        print(f"    TYR {chain} {resseq} processado")


def processar_cys(
    atoms: List[Dict],
    locked: Set[int],
    resseq: int = 580,
    chain: str = "B",
    verbose: bool = True
) -> None:
    """
    Processa a cadeia lateral de Cisteína (CYS).

    Estrutura: CB -> SG -> HG1
    """
    if verbose:
        print(f"\n  Processando CYS {chain} {resseq}...")

    res_idxs = find_residue_atoms(atoms, "CYS", resseq, chain)
    if not res_idxs:
        if verbose:
            print(f"    AVISO: CYS {chain} {resseq} não encontrado")
        return

    name2idx = build_name_index(atoms, res_idxs)

    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_SG = name2idx.get("SG")
    idx_HG1 = name2idx.get("HG1")
    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")

    # HB1/HB2
    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CB, alvo_hb,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hb), label="HB"
        )

    # CB -> SG (enxofre - distâncias maiores)
    if idx_SG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.60, dmax=1.90,
                ang1_min=100.0, ang1_max=125.0,
                ang2_min=104.0, ang2_max=121.0,
                alvo_angulo=112.5, label="SG"
            )
            if not same_coords(atoms[idx_SG], atoms[cand_idx]):
                swap_coords(atoms[idx_SG], atoms[cand_idx])
            locked.add(idx_SG)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO SG: {e}")

    # SG -> HG1
    if idx_HG1 is not None and idx_SG is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_SG, [idx_HG1],
            dmin_init=1.20, dmax_init=1.45,
            target_count=1, label="HG1"
        )

    if verbose:
        print(f"    CYS {chain} {resseq} processado")


def processar_lys(
    atoms: List[Dict],
    locked: Set[int],
    resseq: int = 581,
    chain: str = "B",
    verbose: bool = True
) -> None:
    """
    Processa a cadeia lateral de Lisina (LYS).

    Estrutura longa: CB -> CG -> CD -> CE -> NZ (amino)
    """
    if verbose:
        print(f"\n  Processando LYS {chain} {resseq}...")

    res_idxs = find_residue_atoms(atoms, "LYS", resseq, chain)
    if not res_idxs:
        if verbose:
            print(f"    AVISO: LYS {chain} {resseq} não encontrado")
        return

    name2idx = build_name_index(atoms, res_idxs)

    idx_CA = name2idx.get("CA")
    idx_CB = name2idx.get("CB")
    idx_CG = name2idx.get("CG")
    idx_CD = name2idx.get("CD")
    idx_CE = name2idx.get("CE")
    idx_NZ = name2idx.get("NZ")

    # HB1/HB2
    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")
    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CB, alvo_hb,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hb), label="HB"
        )

    # CB -> CG
    if idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0,
                ang2_min=103.0, ang2_max=124.0,
                alvo_angulo=113.6, label="CG"
            )
            if not same_coords(atoms[idx_CG], atoms[cand_idx]):
                swap_coords(atoms[idx_CG], atoms[cand_idx])
            locked.add(idx_CG)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CG: {e}")

    # HG1/HG2
    idx_HG1 = name2idx.get("HG1")
    idx_HG2 = name2idx.get("HG2")
    alvo_hg = [x for x in [idx_HG1, idx_HG2] if x is not None]
    if alvo_hg and idx_CG is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CG, alvo_hg,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hg), label="HG"
        )

    # CG -> CD
    if idx_CD is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CG, idx_CB,
                dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0,
                ang2_min=103.0, ang2_max=124.0,
                alvo_angulo=113.6, label="CD"
            )
            if not same_coords(atoms[idx_CD], atoms[cand_idx]):
                swap_coords(atoms[idx_CD], atoms[cand_idx])
            locked.add(idx_CD)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CD: {e}")

    # HD1/HD2
    idx_HD1 = name2idx.get("HD1")
    idx_HD2 = name2idx.get("HD2")
    alvo_hd = [x for x in [idx_HD1, idx_HD2] if x is not None]
    if alvo_hd and idx_CD is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CD, alvo_hd,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hd), label="HD"
        )

    # CD -> CE
    if idx_CE is not None and idx_CD is not None and idx_CG is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CD, idx_CG,
                dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0,
                ang2_min=103.0, ang2_max=124.0,
                alvo_angulo=113.6, label="CE"
            )
            if not same_coords(atoms[idx_CE], atoms[cand_idx]):
                swap_coords(atoms[idx_CE], atoms[cand_idx])
            locked.add(idx_CE)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CE: {e}")

    # HE1/HE2
    idx_HE1 = name2idx.get("HE1")
    idx_HE2 = name2idx.get("HE2")
    alvo_he = [x for x in [idx_HE1, idx_HE2] if x is not None]
    if alvo_he and idx_CE is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CE, alvo_he,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_he), label="HE"
        )

    # CE -> NZ
    if idx_NZ is not None and idx_CE is not None and idx_CD is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CE, idx_CD,
                dmin=1.38, dmax=1.58,
                ang1_min=95.0, ang1_max=125.0,
                ang2_min=100.0, ang2_max=120.0,
                alvo_angulo=110.0, label="NZ"
            )
            if not same_coords(atoms[idx_NZ], atoms[cand_idx]):
                swap_coords(atoms[idx_NZ], atoms[cand_idx])
            locked.add(idx_NZ)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO NZ: {e}")

    # HZ1/HZ2/HZ3
    idx_HZ1 = name2idx.get("HZ1")
    idx_HZ2 = name2idx.get("HZ2")
    idx_HZ3 = name2idx.get("HZ3")
    alvo_hz = [x for x in [idx_HZ1, idx_HZ2, idx_HZ3] if x is not None]
    if alvo_hz and idx_NZ is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_NZ, alvo_hz,
            dmin_init=0.90, dmax_init=1.10,
            target_count=len(alvo_hz), label="HZ"
        )

    if verbose:
        print(f"    LYS {chain} {resseq} processado")


def processar_arg(
    atoms: List[Dict],
    locked: Set[int],
    resseq: int = 582,
    chain: str = "B",
    verbose: bool = True
) -> None:
    """
    Processa a cadeia lateral de Arginina (ARG).

    Estrutura com guanidínio: CB -> CG -> CD -> NE -> CZ -> NH1, NH2
    """
    if verbose:
        print(f"\n  Processando ARG {chain} {resseq}...")

    res_idxs = find_residue_atoms(atoms, "ARG", resseq, chain)
    if not res_idxs:
        if verbose:
            print(f"    AVISO: ARG {chain} {resseq} não encontrado")
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

    # HB1/HB2
    idx_HB1 = name2idx.get("HB1")
    idx_HB2 = name2idx.get("HB2")
    alvo_hb = [x for x in [idx_HB1, idx_HB2] if x is not None]
    if alvo_hb and idx_CB is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CB, alvo_hb,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hb), label="HB"
        )

    # CB -> CG
    if idx_CG is not None and idx_CB is not None and idx_CA is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CB, idx_CA,
                dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0,
                ang2_min=103.0, ang2_max=124.0,
                alvo_angulo=113.6, label="CG"
            )
            if not same_coords(atoms[idx_CG], atoms[cand_idx]):
                swap_coords(atoms[idx_CG], atoms[cand_idx])
            locked.add(idx_CG)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CG: {e}")

    # HG1/HG2
    idx_HG1 = name2idx.get("HG1")
    idx_HG2 = name2idx.get("HG2")
    alvo_hg = [x for x in [idx_HG1, idx_HG2] if x is not None]
    if alvo_hg and idx_CG is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CG, alvo_hg,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hg), label="HG"
        )

    # CG -> CD
    if idx_CD is not None and idx_CG is not None and idx_CB is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CG, idx_CB,
                dmin=1.40, dmax=1.60,
                ang1_min=100.0, ang1_max=127.0,
                ang2_min=103.0, ang2_max=124.0,
                alvo_angulo=113.6, label="CD"
            )
            if not same_coords(atoms[idx_CD], atoms[cand_idx]):
                swap_coords(atoms[idx_CD], atoms[cand_idx])
            locked.add(idx_CD)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CD: {e}")

    # HD1/HD2
    idx_HD1 = name2idx.get("HD1")
    idx_HD2 = name2idx.get("HD2")
    alvo_hd = [x for x in [idx_HD1, idx_HD2] if x is not None]
    if alvo_hd and idx_CD is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_CD, alvo_hd,
            dmin_init=0.995, dmax_init=1.115,
            target_count=len(alvo_hd), label="HD"
        )

    # CD -> NE
    if idx_NE is not None and idx_CD is not None and idx_CG is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CD, idx_CG,
                dmin=1.40, dmax=1.60,
                ang1_min=103.0, ang1_max=133.0,
                ang2_min=108.0, ang2_max=128.0,
                alvo_angulo=118.0, label="NE"
            )
            if not same_coords(atoms[idx_NE], atoms[cand_idx]):
                swap_coords(atoms[idx_NE], atoms[cand_idx])
            locked.add(idx_NE)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO NE: {e}")

    # HE
    idx_HE = name2idx.get("HE")
    if idx_HE is not None and idx_NE is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_NE, [idx_HE],
            dmin_init=0.90, dmax_init=1.10,
            target_count=1, label="HE"
        )

    # NE -> CZ (planar)
    if idx_CZ is not None and idx_NE is not None and idx_CD is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_NE, idx_CD,
                dmin=1.25, dmax=1.45,
                ang1_min=100.0, ang1_max=130.0,
                ang2_min=114.0, ang2_max=126.0,
                alvo_angulo=120.0, label="CZ"
            )
            if not same_coords(atoms[idx_CZ], atoms[cand_idx]):
                swap_coords(atoms[idx_CZ], atoms[cand_idx])
            locked.add(idx_CZ)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO CZ: {e}")

    # CZ -> NH1
    if idx_NH1 is not None and idx_CZ is not None and idx_NE is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CZ, idx_NE,
                dmin=1.25, dmax=1.45,
                ang1_min=100.0, ang1_max=130.0,
                ang2_min=114.0, ang2_max=126.0,
                alvo_angulo=120.0, label="NH1"
            )
            if not same_coords(atoms[idx_NH1], atoms[cand_idx]):
                swap_coords(atoms[idx_NH1], atoms[cand_idx])
            locked.add(idx_NH1)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO NH1: {e}")

    # 1HH1/2HH1
    idx_1HH1 = name2idx.get("1HH1")
    idx_2HH1 = name2idx.get("2HH1")
    alvo_hh1 = [x for x in [idx_1HH1, idx_2HH1] if x is not None]
    if alvo_hh1 and idx_NH1 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_NH1, alvo_hh1,
            dmin_init=0.90, dmax_init=1.10,
            target_count=len(alvo_hh1), label="HH1"
        )

    # CZ -> NH2
    if idx_NH2 is not None and idx_CZ is not None and idx_NE is not None:
        try:
            cand_idx, d, ang = escolher_pesado_com_angulo(
                atoms, locked, idx_CZ, idx_NE,
                dmin=1.25, dmax=1.45,
                ang1_min=100.0, ang1_max=130.0,
                ang2_min=114.0, ang2_max=126.0,
                alvo_angulo=120.0, label="NH2"
            )
            if not same_coords(atoms[idx_NH2], atoms[cand_idx]):
                swap_coords(atoms[idx_NH2], atoms[cand_idx])
            locked.add(idx_NH2)
        except RuntimeError as e:
            if verbose:
                print(f"    AVISO NH2: {e}")

    # 1HH2/2HH2
    idx_1HH2 = name2idx.get("1HH2")
    idx_2HH2 = name2idx.get("2HH2")
    alvo_hh2 = [x for x in [idx_1HH2, idx_2HH2] if x is not None]
    if alvo_hh2 and idx_NH2 is not None:
        atribuir_coord_alvos(
            atoms, locked, idx_NH2, alvo_hh2,
            dmin_init=0.90, dmax_init=1.10,
            target_count=len(alvo_hh2), label="HH2"
        )

    if verbose:
        print(f"    ARG {chain} {resseq} processado")


# ==============================================================================
# SEÇÃO 12: PIPELINE PRINCIPAL UNIFICADO
# ==============================================================================

def executar_pipeline_completo(
    pdb_entrada: str,
    pdb_saida: str,
    start_serial: int = 8641,
    end_serial: int = 8776,
    chain: str = "B",
    verbose: bool = True
) -> None:
    """
    Executa o pipeline completo de reordenação da estrutura CAR-T.

    Este é o ponto de entrada principal que executa todas as 10 fases
    de processamento de forma automatizada, sem gerar arquivos intermediários.

    Args:
        pdb_entrada: Caminho do arquivo PDB de entrada
        pdb_saida: Caminho do arquivo PDB de saída final
        start_serial: Serial do N inicial (default: 8641)
        end_serial: Serial do C final (default: 8776)
        chain: Identificador da cadeia (default: "B")
        verbose: Se True, imprime progresso detalhado

    Returns:
        None (escreve o arquivo PDB de saída)
    """
    print("="*80)
    print("PIPELINE UNIFICADO DE REORDENAÇÃO CAR-T")
    print("="*80)
    print(f"Arquivo de entrada: {pdb_entrada}")
    print(f"Arquivo de saída: {pdb_saida}")
    print("="*80)

    # Carregar PDB
    atoms, lines = parse_pdb_file(pdb_entrada)
    print(f"\nTotal de átomos carregados: {len(atoms)}")

    # Criar dicionário de coordenadas por serial
    coords = extrair_coords_de_atoms(atoms)

    # =========================================================================
    # FASE 1: Reconstrução do Backbone
    # =========================================================================
    coords, backbone_ids = fase1_backbone(
        atoms, coords, start_serial, end_serial, verbose
    )

    # Aplicar coordenadas de volta aos átomos
    aplicar_coords_para_atoms(atoms, coords)

    # =========================================================================
    # FASE 2: Atribuição de HN, HA, O
    # =========================================================================
    coords = extrair_coords_de_atoms(atoms)
    coords = fase2_hn_ha_o(atoms, coords, backbone_ids, verbose)
    aplicar_coords_para_atoms(atoms, coords)

    # =========================================================================
    # FASE 3: Conexão CA-CB
    # =========================================================================
    coords = extrair_coords_de_atoms(atoms)
    locked_serials = set(backbone_ids)

    # Adicionar HN/HA/O ao locked
    atom_by_serial = {a['serial']: a for a in atoms}
    residues = get_backbone_residues(backbone_ids, atom_by_serial)
    for res in residues:
        for aname in ('HN', 'HA', 'O'):
            for a in atoms:
                if a['chain'] == res['chain'] and a['resseq'] == res['resseq'] and a['name'] == aname:
                    locked_serials.add(a['serial'])

    coords, locked_serials = fase3_ca_cb(atoms, coords, backbone_ids, locked_serials, verbose)
    aplicar_coords_para_atoms(atoms, coords)

    # =========================================================================
    # FASES 4-10: Processamento de Cadeias Laterais
    # =========================================================================
    if verbose:
        print("\n" + "="*80)
        print("FASES 4-10: PROCESSAMENTO DE CADEIAS LATERAIS")
        print("="*80)

    # Inicializar locked com base nos índices (não seriais)
    locked = inicializar_locked_base(atoms)

    # FASE 4: ILE 576
    processar_ile(atoms, locked, resseq=576, chain=chain, verbose=verbose)
    bloquear_residuo(atoms, locked, "ILE", 576, chain)

    # FASE 5: THR 577
    processar_thr(atoms, locked, resseq=577, chain=chain, verbose=verbose)
    bloquear_residuo(atoms, locked, "THR", 577, chain)

    # FASE 6: LEU 578
    processar_leu(atoms, locked, resseq=578, chain=chain, verbose=verbose)
    bloquear_residuo(atoms, locked, "LEU", 578, chain)

    # FASE 7: TYR 579
    processar_tyr(atoms, locked, resseq=579, chain=chain, verbose=verbose)
    bloquear_residuo(atoms, locked, "TYR", 579, chain)

    # FASE 8: CYS 580
    processar_cys(atoms, locked, resseq=580, chain=chain, verbose=verbose)
    bloquear_residuo(atoms, locked, "CYS", 580, chain)

    # FASE 9: LYS 581
    processar_lys(atoms, locked, resseq=581, chain=chain, verbose=verbose)
    bloquear_residuo(atoms, locked, "LYS", 581, chain)

    # FASE 10: ARG 582
    processar_arg(atoms, locked, resseq=582, chain=chain, verbose=verbose)
    bloquear_residuo(atoms, locked, "ARG", 582, chain)

    # =========================================================================
    # Escrita do arquivo final
    # =========================================================================
    lines = update_pdb_lines(lines, atoms)
    write_pdb_file(pdb_saida, lines)

    print("\n" + "="*80)
    print("PIPELINE CONCLUÍDO COM SUCESSO!")
    print("="*80)
    print(f"Arquivo final: {pdb_saida}")


# ==============================================================================
# SEÇÃO 13: INTERFACE DE LINHA DE COMANDO
# ==============================================================================

def main():
    """
    Função principal para execução via linha de comando.
    """
    import sys

    # Uso no Google Colab
    try:
        from google.colab import files
        print("="*60)
        print("PIPELINE UNIFICADO DE REORDENAÇÃO CAR-T")
        print("="*60)
        print("\nSelecione o arquivo PDB de entrada:")
        uploaded = files.upload()

        if not uploaded:
            raise RuntimeError("Nenhum arquivo foi enviado.")

        pdb_entrada = list(uploaded.keys())[0]
        pdb_saida = "estrutura_cart_reordenada_final.pdb"

        executar_pipeline_completo(
            pdb_entrada=pdb_entrada,
            pdb_saida=pdb_saida,
            verbose=True
        )

        # Oferecer download
        files.download(pdb_saida)

    except ImportError:
        # Uso local via linha de comando
        if len(sys.argv) < 2:
            print("Uso: python cart_pipeline_unificado.py <arquivo_entrada.pdb> [arquivo_saida.pdb]")
            print("\nExemplo:")
            print("  python cart_pipeline_unificado.py md1amd12-frame-unico-576-583.pdb estrutura_final.pdb")
            sys.exit(1)

        pdb_entrada = sys.argv[1]
        pdb_saida = sys.argv[2] if len(sys.argv) > 2 else "estrutura_cart_reordenada_final.pdb"

        executar_pipeline_completo(
            pdb_entrada=pdb_entrada,
            pdb_saida=pdb_saida,
            verbose=True
        )


if __name__ == "__main__":
    main()
