# Documentação Técnica: Pipeline de Reordenação de Estrutura CAR-T

## Versão 1.0 - Multi-Frame
**Data:** Dezembro 2025
**Arquivo:** `reordenacao-cart-v1.py`

---

## Índice

1. [Visão Geral](#1-visão-geral)
2. [Arquitetura do Sistema](#2-arquitetura-do-sistema)
3. [Estruturas de Dados](#3-estruturas-de-dados)
4. [Algoritmos Principais](#4-algoritmos-principais)
5. [Referência de Funções](#5-referência-de-funções)
6. [Parâmetros de Configuração](#6-parâmetros-de-configuração)
7. [Fluxo de Processamento](#7-fluxo-de-processamento)
8. [Tolerâncias e Critérios Geométricos](#8-tolerâncias-e-critérios-geométricos)
9. [Tratamento de Erros](#9-tratamento-de-erros)
10. [Exemplos de Uso](#10-exemplos-de-uso)
11. [Considerações de Performance](#11-considerações-de-performance)

---

## 1. Visão Geral

### 1.1 Propósito

O Pipeline de Reordenação CAR-T é um sistema automatizado para reconstrução e validação de estruturas moleculares de receptores CAR-T (Chimeric Antigen Receptor T-cell). O pipeline processa arquivos PDB contendo trajetórias de dinâmica molecular, reordenando as coordenadas atômicas para garantir a correta conectividade química.

### 1.2 Problema Resolvido

Em simulações de dinâmica molecular geradas pelo GROMACS, as coordenadas atômicas podem estar desordenadas em relação à estrutura química esperada. Este pipeline:

- Reconstrói o backbone proteico (N-CA-C)
- Atribui corretamente átomos de hidrogênio e oxigênio do backbone (HN, HA, O)
- Estabelece conexões CA-CB
- Processa cadeias laterais de 7 resíduos específicos (ILE 576 a ARG 582)

### 1.3 Capacidades

- **Multi-frame:** Processa trajetórias completas com múltiplos frames
- **Automático:** Não requer intervenção manual entre fases
- **Sem intermediários:** Não gera arquivos PDB intermediários
- **Tolerante a falhas:** Continua processamento mesmo se um frame falhar

---

## 2. Arquitetura do Sistema

### 2.1 Diagrama de Módulos

```
┌─────────────────────────────────────────────────────────────────┐
│                    PIPELINE PRINCIPAL                            │
│  executar_pipeline_multiframe()                                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    DIVISÃO DE FRAMES                             │
│  split_trajectory_into_frames()                                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│              PROCESSAMENTO POR FRAME                             │
│  processar_frame()                                               │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │ FASE 1: Backbone (DFS)          fase1_backbone()           │ │
│  │ FASE 2: HN/HA/O                 fase2_hn_ha_o()            │ │
│  │ FASE 3: CA-CB                   fase3_ca_cb()              │ │
│  │ FASES 4-10: Cadeias Laterais    processar_*()              │ │
│  └────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    ESCRITA DE SAÍDA                              │
│  Combinação de todos os frames processados                       │
└─────────────────────────────────────────────────────────────────┘
```

### 2.2 Organização do Código

| Seção | Descrição | Linhas |
|-------|-----------|--------|
| 1 | Funções Geométricas Fundamentais | `dist`, `angle`, `swap_coords` |
| 2 | Parsing e Escrita de PDB Multi-Frame | `split_trajectory_into_frames`, `parse_frame_content` |
| 3 | Busca e Indexação de Átomos | `find_residue_atoms`, `build_name_index` |
| 4 | Construção do Backbone Canônico | `build_backbone_order`, `get_backbone_residues` |
| 5 | Limites de Ligação Química | `bond_limits_backbone` |
| 6 | Algoritmos de Seleção de Átomos | `escolher_vizinhos_dinamico`, `escolher_pesado_com_angulo` |
| 7 | Fase 1 - Backbone DFS | `find_all_backbone_paths`, `fase1_backbone` |
| 8 | Fase 2 - HN/HA/O | `fase2_hn_ha_o` |
| 9 | Fase 3 - CA-CB | `fase3_ca_cb` |
| 10 | Funções Auxiliares | `aplicar_coords_para_atoms`, `bloquear_residuo` |
| 11 | Processamento de Cadeias Laterais | `processar_ile`, `processar_thr`, etc. |
| 12 | Processamento de Frame Único | `processar_frame` |
| 13 | Pipeline Multi-Frame | `executar_pipeline_multiframe` |
| 14 | Interface CLI | `main` |

---

## 3. Estruturas de Dados

### 3.1 Dicionário de Átomo

Cada átomo é representado por um dicionário Python:

```python
{
    'line_idx': int,      # Índice da linha no arquivo PDB
    'serial': int,        # Número serial do átomo (colunas 7-11)
    'name': str,          # Nome do átomo (colunas 13-16), ex: "CA", "CB"
    'altloc': str,        # Localização alternativa (coluna 17)
    'resname': str,       # Nome do resíduo (colunas 18-20), ex: "ILE", "TYR"
    'chain': str,         # Identificador da cadeia (coluna 22)
    'resseq': int,        # Número do resíduo (colunas 23-26)
    'icode': str,         # Código de inserção (coluna 27)
    'x': float,           # Coordenada X (colunas 31-38)
    'y': float,           # Coordenada Y (colunas 39-46)
    'z': float,           # Coordenada Z (colunas 47-54)
    'raw_line': str       # Linha original do PDB
}
```

### 3.2 Dicionário de Coordenadas

Mapeamento de serial para coordenadas:

```python
coords: Dict[int, Tuple[float, float, float]]
# Exemplo: {8641: (10.5, 20.3, 15.7), 8642: (11.2, 19.8, 16.1), ...}
```

### 3.3 Conjunto Locked

Conjunto de índices de átomos cujas coordenadas já foram atribuídas:

```python
locked: Set[int]
# Exemplo: {0, 1, 2, 5, 7, 12, ...}
```

### 3.4 Estado do DFS (Fase 1)

```python
{
    "step": int,                    # Passo atual (0-23)
    "coords": Dict[int, Tuple],     # Coordenadas neste caminho
    "locked": Set[int],             # Seriais bloqueados
    "total": float,                 # Distância total acumulada
    "edges": List[Dict]             # Histórico de arestas
}
```

---

## 4. Algoritmos Principais

### 4.1 Algoritmo DFS para Backbone (Fase 1)

O backbone é reconstruído usando busca em profundidade (DFS) com backtracking:

```
ALGORITMO: find_all_backbone_paths

ENTRADA:
  - atoms: lista de átomos
  - backbone_ids: sequência canônica [N576, CA576, C576, N577, ...]
  - coords: coordenadas atuais

PROCESSO:
  1. Inicializar pilha com estado inicial (N576 como vértice 0)
  2. Enquanto pilha não vazia:
     a. Desempilhar estado
     b. Se completou 23 passos:
        - Verificar se distância total está em [30, 40] Å
        - Se sim, retornar solução
     c. Para cada átomo candidato não bloqueado:
        - Calcular distância para próximo vértice canônico
        - Se dentro dos limites de ligação:
          - Criar novo estado com swap de coordenadas
          - Empilhar novo estado

CRITÉRIOS DE PODA:
  - Distância total > 40 Å antes de completar 23 passos
  - Nenhum candidato válido para o próximo passo

SAÍDA:
  - coords: coordenadas com backbone reconstruído
  - backbone_ids: lista de seriais do backbone
```

**Tolerâncias de Ligação do Backbone:**

| Ligação | Distância Mínima (Å) | Distância Máxima (Å) |
|---------|---------------------|---------------------|
| N-CA | 1.40 | 1.60 |
| CA-C | 1.40 | 1.70 |
| C-N | 1.25 | 1.45 |

### 4.2 Algoritmo de Refinamento Angular (REF1/REF2/REF3)

Usado para seleção de átomos pesados com critérios de distância e ângulo:

```
ALGORITMO: escolher_pesado_com_angulo

ENTRADA:
  - candidatos: átomos dentro da janela de distância
  - ang1_min, ang1_max: limites REF1 (amplo)
  - ang2_min, ang2_max: limites REF2 (estreito)
  - alvo_angulo: ângulo alvo para REF3

PROCESSO:
  REF1: Filtrar candidatos por [ang1_min, ang1_max]
        Se 1 candidato → retornar
        Se 0 candidatos (vindo de ≥2) → aplicar REF3 na lista anterior

  REF2: Filtrar REF1 por [ang2_min, ang2_max]
        Se 1 candidato → retornar
        Se 0 candidatos (vindo de ≥2) → aplicar REF3 na lista REF1

  REF3: Ordenar por |ângulo - alvo_angulo|
        Retornar primeiro (mais próximo do alvo)

SAÍDA:
  - (índice, distância, ângulo) do candidato escolhido
```

### 4.3 Algoritmo de Janela Dinâmica para Hidrogênios

```
ALGORITMO: escolher_vizinhos_dinamico

ENTRADA:
  - idx_ref: átomo de referência
  - target_count: número de vizinhos desejados
  - dmin_init, dmax_init: janela inicial
  - delta: incremento (default: 0.005 Å)

PROCESSO:
  Repetir até max_iter:
    1. Contar candidatos em [dmin, dmax]
    2. Se count == target_count → retornar
    3. Se count > target_count:
       dmin += delta, dmax -= delta
    4. Se count < target_count:
       dmin -= delta, dmax += delta

  Fallback: retornar os target_count mais próximos

SAÍDA:
  - lista de índices candidatos
  - janela final usada
```

---

## 5. Referência de Funções

### 5.1 Funções Geométricas

#### `dist(a, b) → float`
Calcula distância euclidiana 3D entre dois átomos.

#### `dist_tuple(p1, p2) → float`
Calcula distância entre duas tuplas (x, y, z).

#### `angle(a, b, c) → float`
Calcula ângulo A-B-C em graus (vértice em B).

#### `same_coords(a, b, tol=1e-3) → bool`
Verifica se dois átomos têm mesmas coordenadas.

#### `swap_coords(a, b) → None`
Troca coordenadas entre dois átomos (in-place).

### 5.2 Funções de Parsing

#### `split_trajectory_into_frames(content) → List[str]`
Divide conteúdo de trajetória em frames usando `TER`/`ENDMDL`.

**Parâmetros:**
- `content`: String com conteúdo completo do arquivo

**Retorno:**
- Lista de strings, cada uma contendo um frame completo

#### `parse_frame_content(content) → Tuple[List[Dict], List[str]]`
Faz parsing de um único frame.

**Retorno:**
- `atoms`: Lista de dicionários de átomos
- `lines`: Lista de linhas originais

#### `update_pdb_lines(lines, atoms) → List[str]`
Atualiza coordenadas nas linhas PDB.

### 5.3 Funções de Busca

#### `find_residue_atoms(atoms, resname, resseq, chain=None) → List[int]`
Encontra índices de átomos de um resíduo específico.

#### `build_name_index(atoms, residue_idxs) → Dict[str, int]`
Cria mapa nome→índice para átomos de um resíduo.

### 5.4 Funções de Processamento

#### `processar_frame(frame_content, frame_num, ...) → str`
Processa um único frame da trajetória.

**Parâmetros:**
- `frame_content`: Conteúdo do frame
- `frame_num`: Número do frame (para logs)
- `start_serial`: Serial do N inicial (default: 8641)
- `end_serial`: Serial do C final (default: 8776)
- `chain`: Cadeia (default: "B")

**Retorno:**
- String com frame processado

#### `executar_pipeline_multiframe(pdb_entrada, pdb_saida, ...) → None`
Função principal que processa toda a trajetória.

### 5.5 Funções de Processamento de Resíduos

| Função | Resíduo | Estrutura |
|--------|---------|-----------|
| `processar_ile()` | ILE 576 | CB→CG1→CD, CB→CG2 |
| `processar_thr()` | THR 577 | CB→OG1→HG1, CB→CG2 |
| `processar_leu()` | LEU 578 | CB→CG→CD1, CG→CD2 |
| `processar_tyr()` | TYR 579 | Anel aromático + OH |
| `processar_cys()` | CYS 580 | CB→SG→HG1 |
| `processar_lys()` | LYS 581 | CB→CG→CD→CE→NZ |
| `processar_arg()` | ARG 582 | CB→CG→CD→NE→CZ→NH1/NH2 |

---

## 6. Parâmetros de Configuração

### 6.1 Parâmetros Globais

| Parâmetro | Valor Default | Descrição |
|-----------|---------------|-----------|
| `start_serial` | 8641 | Serial do N do primeiro resíduo (ILE 576) |
| `end_serial` | 8776 | Serial do C do último resíduo (GLY 583) |
| `chain` | "B" | Identificador da cadeia |
| `delta` | 0.005 Å | Incremento para janela dinâmica |
| `max_iter` | 200 | Máximo de iterações para busca |

### 6.2 Parâmetros do Backbone

| Parâmetro | Valor |
|-----------|-------|
| Número de passos | 23 (24 vértices) |
| Distância total mínima | 30 Å |
| Distância total máxima | 40 Å |

---

## 8. Tolerâncias e Critérios Geométricos

### 8.1 Distâncias de Ligação

#### Backbone
| Ligação | Min (Å) | Max (Å) |
|---------|---------|---------|
| N-CA | 1.40 | 1.60 |
| CA-C | 1.40 | 1.70 |
| C-N (peptídica) | 1.25 | 1.45 |

#### Carbono-Hidrogênio (sp3)
| Ligação | Min (Å) | Max (Å) |
|---------|---------|---------|
| C-H padrão | 0.995 | 1.115 |
| C-H (ILE/THR) | 0.990 | 1.150 |

#### Nitrogênio-Hidrogênio
| Ligação | Min (Å) | Max (Å) |
|---------|---------|---------|
| N-H (backbone HN) | 0.90 | 1.10 |
| N-H (amino LYS/ARG) | 0.90 | 1.10 |

#### Oxigênio-Hidrogênio
| Ligação | Min (Å) | Max (Å) |
|---------|---------|---------|
| O-H (THR, TYR) | 0.85 | 1.05 |

#### Enxofre-Hidrogênio
| Ligação | Min (Å) | Max (Å) |
|---------|---------|---------|
| S-H (CYS) | 1.20 | 1.45 |

#### Carbono-Carbono
| Ligação | Min (Å) | Max (Å) | Tipo |
|---------|---------|---------|------|
| CA-CB | 1.40 | 1.60 | sp3 |
| CB-CG | 1.40 | 1.60 | sp3 |
| CG-CD | 1.40 | 1.60 | sp3 |
| C-C aromático | 1.30 | 1.50 | sp2 |

#### Carbono-Heteroátomo
| Ligação | Min (Å) | Max (Å) | Resíduo |
|---------|---------|---------|---------|
| CB-OG1 | 1.30 | 1.50 | THR |
| CB-SG | 1.60 | 1.90 | CYS |
| CE-NZ | 1.38 | 1.58 | LYS |
| CD-NE | 1.40 | 1.60 | ARG |
| NE-CZ | 1.25 | 1.45 | ARG |
| CZ-NH | 1.25 | 1.45 | ARG |
| CZ-OH | 1.30 | 1.50 | TYR |

### 8.2 Critérios Angulares

#### Carbonos Alifáticos (sp3) - Padrão
| Refinamento | Min (°) | Max (°) | Alvo (°) |
|-------------|---------|---------|----------|
| REF1 | 101 | 119 | - |
| REF2 | 104 | 116 | - |
| REF3 | - | - | 110 |

#### Carbonos Aromáticos (sp2) - TYR
| Refinamento | Min (°) | Max (°) | Alvo (°) |
|-------------|---------|---------|----------|
| REF1 | 114 | 126 | - |
| REF2 | 117 | 123 | - |
| REF3 | - | - | 120 |

#### CYS (CB-SG)
| Refinamento | Min (°) | Max (°) | Alvo (°) |
|-------------|---------|---------|----------|
| REF1 | 100 | 125 | - |
| REF2 | 104 | 121 | - |
| REF3 | - | - | 112.5 |

#### LYS (cadeia lateral)
| Refinamento | Min (°) | Max (°) | Alvo (°) |
|-------------|---------|---------|----------|
| REF1 | 100 | 127 | - |
| REF2 | 103 | 124 | - |
| REF3 | - | - | 113.6 |

#### LYS (CE-NZ)
| Refinamento | Min (°) | Max (°) | Alvo (°) |
|-------------|---------|---------|----------|
| REF1 | 95 | 125 | - |
| REF2 | 100 | 120 | - |
| REF3 | - | - | 110 |

#### ARG (CG, CD) - sp3
| Refinamento | Min (°) | Max (°) | Alvo (°) |
|-------------|---------|---------|----------|
| REF1 | 100 | 127 | - |
| REF2 | 103 | 124 | - |
| REF3 | - | - | 113.6 |

#### ARG (CD-NE)
| Refinamento | Min (°) | Max (°) | Alvo (°) |
|-------------|---------|---------|----------|
| REF1 | 103 | 133 | - |
| REF2 | 108 | 128 | - |
| REF3 | - | - | 118 |

#### ARG (NE-CZ, CZ-NH) - planar
| Refinamento | Min (°) | Max (°) | Alvo (°) |
|-------------|---------|---------|----------|
| REF1 | 100 | 130 | - |
| REF2 | 114 | 126 | - |
| REF3 | - | - | 120 |

#### CA-CB (Fase 3)
| Ângulo | Min (°) | Max (°) | Alvo (°) |
|--------|---------|---------|----------|
| N-CA-CB | 104 | 116 | 110 |
| C-CA-CB | 105 | 118 | 111 |

---

## 9. Tratamento de Erros

### 9.1 Erros Recuperáveis

| Situação | Comportamento |
|----------|---------------|
| Frame sem átomos | Retorna frame original |
| Resíduo não encontrado | Pula processamento do resíduo |
| Nenhum candidato para átomo | Pula atribuição, continua |
| Falha em fase específica | Continua com próxima fase |

### 9.2 Erros Fatais

| Situação | Mensagem |
|----------|----------|
| Arquivo não encontrado | `FileNotFoundError` |
| Nenhum frame detectado | "Nenhum frame encontrado" |
| Backbone impossível | "Nenhum caminho válido encontrado" |

### 9.3 Relatório de Falhas

Ao final do processamento, o sistema exibe um relatório completo:

```
================================================================================
RELATÓRIO DE PROCESSAMENTO
================================================================================
Total de frames:      1000
Frames com sucesso:   997
Frames com falha:     3

--------------------------------------------------------------------------------
FRAMES COM FALHA (mantidos com coordenadas originais):
--------------------------------------------------------------------------------
  Frame 42: Nenhum caminho válido encontrado para o backbone
  Frame 156: Nenhum candidato encontrado para CG
  Frame 891: Frame vazio ou sem átomos válidos

================================================================================
Arquivo de saída: trajetoria_reordenada.pdb
================================================================================
PIPELINE CONCLUÍDO!
================================================================================
```

### 9.4 Logging de Progresso

O pipeline usa `print()` para feedback em tempo real:

```python
# Progresso durante processamento
"Processando frame 50/1000..."

# Conclusão
"Processamento concluído!"
```

---

## 10. Exemplos de Uso

### 10.1 Linha de Comando

```bash
# Uso básico
python reordenacao-cart-v1.py trajetoria.pdb saida.pdb

# Com nome de saída automático
python reordenacao-cart-v1.py trajetoria.pdb
# Saída: trajetoria_cart_reordenada.pdb
```

### 10.2 Google Colab

```python
# Executar diretamente - fará upload/download automaticamente
!python reordenacao-cart-v1.py
```

### 10.3 Como Módulo Python

```python
from reordenacao_cart_v1 import executar_pipeline_multiframe

# Processar trajetória
executar_pipeline_multiframe(
    pdb_entrada="simulacao.pdb",
    pdb_saida="simulacao_reordenada.pdb",
    start_serial=8641,
    end_serial=8776,
    chain="B",
    verbose=True
)
```

### 10.4 Processamento de Frame Único

```python
from reordenacao_cart_v1 import processar_frame

with open("frame_unico.pdb", "r") as f:
    content = f.read()

resultado = processar_frame(
    frame_content=content,
    frame_num=1,
    start_serial=8641,
    end_serial=8776,
    chain="B"
)

with open("frame_processado.pdb", "w") as f:
    f.write(resultado)
```

---

## 11. Considerações de Performance

### 11.1 Complexidade Computacional

| Fase | Complexidade | Observação |
|------|--------------|------------|
| Parsing | O(n) | n = número de linhas |
| Fase 1 (DFS) | O(m²) no pior caso | m = número de átomos; poda agressiva |
| Fase 2 (HN/HA/O) | O(r × m) | r = número de resíduos |
| Fase 3 (CA-CB) | O(r × m) | Com cálculo de ângulos |
| Fases 4-10 | O(r × m) | Por resíduo específico |

### 11.2 Uso de Memória

- **Por frame:** ~2-3 MB para estrutura com 11.200 átomos
- **Pico:** Durante DFS, múltiplos estados na pilha
- **Otimização:** Frames são processados sequencialmente, não carregados todos em memória

### 11.3 Tempo de Execução Esperado

| Tamanho | Tempo Estimado |
|---------|----------------|
| 1 frame | ~0.5-1 segundo |
| 100 frames | ~1-2 minutos |
| 1000 frames | ~10-15 minutos |

### 11.4 Otimizações Implementadas

1. **Retorno antecipado no DFS:** Retorna primeira solução válida
2. **Poda por distância:** Descarta caminhos > 40 Å antes de completar
3. **Processamento sequencial:** Não carrega toda trajetória em memória
4. **Reutilização de índices:** `build_name_index` evita buscas repetidas

---

## Apêndice A: Formato do Arquivo PDB

### Estrutura de Linha ATOM

```
ATOM    123  CA  ILE B 576      10.500  20.300  15.700  1.00  0.00           C
|       |   |   |   | |        |       |       |       |     |             |
|       |   |   |   | |        |       |       |       |     |             Elemento
|       |   |   |   | |        |       |       |       |     B-factor
|       |   |   |   | |        |       |       |       Occupancy
|       |   |   |   | |        |       |       Z
|       |   |   |   | |        |       Y
|       |   |   |   | |        X
|       |   |   |   | Número do resíduo (22-26)
|       |   |   |   Cadeia (22)
|       |   |   Nome do resíduo (18-20)
|       |   Nome do átomo (13-16)
|       Serial (7-11)
Record type (1-6)
```

### Delimitadores de Frame

```
...
ATOM  11200  OT2 ARG B 734      66.500  70.110 206.180  1.00  0.00           O
TER
ENDMDL
TITLE     Frame 2
MODEL     2
ATOM      1  N   ...
```

---

## Apêndice B: Resíduos Processados

### Sequência CAR-T (576-582)

| # | Código | Nome | Átomos Típicos |
|---|--------|------|----------------|
| 576 | ILE | Isoleucina | 19 |
| 577 | THR | Treonina | 14 |
| 578 | LEU | Leucina | 19 |
| 579 | TYR | Tirosina | 21 |
| 580 | CYS | Cisteína | 11 |
| 581 | LYS | Lisina | 22 |
| 582 | ARG | Arginina | 24 |

### Estruturas das Cadeias Laterais

```
ILE 576:           THR 577:           LEU 578:
    CB                 CB                 CB
   /  \               /  \               |
 CG1  CG2           OG1  CG2            CG
  |                  |                 /  \
 CD                 HG1              CD1  CD2

TYR 579:                    CYS 580:
    CB                          CB
    |                           |
    CG                          SG
   /  \                         |
 CD1  CD2                      HG1
  |    |
CE1  CE2
   \  /
    CZ
    |
    OH

LYS 581:                    ARG 582:
    CB                          CB
    |                           |
    CG                          CG
    |                           |
    CD                          CD
    |                           |
    CE                          NE
    |                           |
    NZ                          CZ
   /|\                         /  \
HZ1 HZ2 HZ3                  NH1  NH2
```

---

## Apêndice C: Histórico de Versões

| Versão | Data | Alterações |
|--------|------|------------|
| 1.0 | Dez/2025 | Versão inicial multi-frame |

---

## Apêndice D: Contato e Suporte

**Repositório:** doutorado
**Branch:** claude/analyze-cart-scripts-*
**Autor:** Pipeline consolidado automaticamente

---

*Documentação gerada em Dezembro de 2025*
