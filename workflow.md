# Rachel — Program Workflow Documentation

## Table of Contents

1. [Overview](#overview)
2. [System Architecture](#system-architecture)
3. [Core Workflow](#core-workflow)
4. [Data Flow](#data-flow)
5. [Key Components](#key-components)
6. [Decision Context Layering](#decision-context-layering)
7. [Chemical Analysis Pipeline](#chemical-analysis-pipeline)
8. [Complete Example: Paracetamol Synthesis](#complete-example-paracetamol-synthesis)
9. [API Reference](#api-reference)

---

## Overview

### What is Rachel?

**Rachel** is an AI-driven computer-assisted retrosynthetic analysis system that uses Large Language Models (LLMs) as the decision-making engine to break down complex organic molecules into simpler, commercially available starting materials.

### Core Purpose

- **Automated Retrosynthetic Planning**: Given a target molecule (SMILES string), the system generates a complete synthetic route by progressively disconnecting bonds and transforming functional groups
- **LLM-Driven Decision Making**: LLMs make strategic decisions about bond disconnections, protecting group strategies, and route optimization
- **Chemical Knowledge Integration**: Combines AI decision-making with rigorous chemical validation using RDKit and custom chemical tools

### Key Features

| Feature | Description |
|---------|-------------|
| **BFS Orchestration** | Breadth-first search explores synthesis tree systematically |
| **Sandbox Mechanism** | Test proposals without committing to synthesis tree |
| **Session Persistence** | JSON-based state allows interrupted sessions to resume |
| **Convergent Design** | Encourages convergent synthetic routes |
| **Forward Validation** | 5-step validation ensures chemically feasible reactions |
| **Layered Context** | Compact/full/status/tree views control token usage |
| **Smart Capping** | Automatic capping group suggestions for bond disconnections |

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              Rachel System Architecture                     │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                        LLM (Decision Maker)                         │    │
│  │                   Makes strategic chemical decisions                │    │
│  └──────────────────────────────┬──────────────────────────────────────┘    │
│                                 │ JSON-in / JSON-out                        │
│                                 ↓                                           │
│  ┌─────────────────────────────────────────────────────────────────────┐    │ 
│  │                        Main Layer (Orchestration)                   │    │
│  │  ┌──────────────┐  ┌──────────────────┐  ┌──────────────────────┐   │    │
│  │  │  RetroCmd    │  │ RetroOrchestrator│  │  RetrosynthesisTree  │   │    │
│  │  │  (23 cmds)   │  │   (BFS Engine)   │  │  (Data Model)        │   │    │
│  │  └──────────────┘  └──────────────────┘  └──────────────────────┘   │    │
│  │  ┌──────────────┐  ┌──────────────────┐  ┌──────────────────────┐   │    │
│  │  │RetroSession  │  │  RetroReport     │  │  RetroVisualizer     │   │    │
│  │  │(Persistence) │  │  (Forward Syn)   │  │  (HTML/MD Reports)   │   │    │
│  │  └──────────────┘  └──────────────────┘  └──────────────────────┘   │    │
│  └──────────────────────────────┬──────────────────────────────────────┘    │
│                                 │                                           │
│                                 ↓                                           │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                      Chemical Tools Layer (26 functions)            │    │
│  │  ┌────┐ ┌────┐ ┌────┐ ┌────┐ ┌────┐ ┌────┐ ┌────┐ ┌────┐            │    │
│  │  │ M0 │ │ M1 │ │ M2 │ │ M3 │ │ M4 │ │ M5 │ │ M6 │ │ M7 │ │ M8 │     │    │
│  │  │Utls│ │Info│ │ FG │ │Tmpl│ │Brk │ │Val │ │ CS │ │Warn│ │Cap │     │    │
│  │  └────┘ └────┘ └────┘ └────┘ └────┘ └────┘ └────┘ └────┘            │    │
│  └──────────────────────────────┬──────────────────────────────────────┘    │
│                                 │                                           │
│                                 ↓                                           │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                           RDKit                                     │    │
│  │                    (Chemistry Informatics)                          │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Layer Descriptions

| Layer | Description | Files |
|-------|-------------|-------|
| **LLM Layer** | Decision-making engine (external) | GPT-4, Claude, etc. |
| **Main Layer** | Orchestration and state management | `main/*.py` (7 files) |
| **Chemical Tools** | 26 functions across 8 modules | `chem_tools/*.py` |
| **RDKit** | Chemical informatics foundation | External library |

---

## Core Workflow

### High-Level Flow Diagram

```
┌─────────────┐
│  Start      │
└──────┬──────┘
       │
       ↓
┌─────────────────────────────────────────────────────────────────┐
│                    INITIALIZATION PHASE                         │
│  1. Validate target SMILES                                      │
│  2. Create RetrosynthesisTree with target molecule              │
│  3. Initialize BFS queue                                        │
│  4. Compute CS score for target                                 │
│  5. Create session file (JSON persistence)                      │
└──────────────────────────────┬──────────────────────────────────┘
                               │
                               ↓
┌─────────────────────────────────────────────────────────────────┐
│                     PLANNING LOOP (BFS)                         │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │  WHILE queue NOT empty AND steps < max_steps:           │    │
│  │                                                         │    │
│  │  1. prepare_next() - Get next molecule from queue       │    │
│  │     → Check termination conditions                      │    │
│  │     → Generate decision context                         │    │
│  │                                                         │    │
│  │  2. LLM Decision Process                                │    │
│  │     a) Explore options (explore_bond, explore_fgi)      │    │
│  │     b) Sandbox testing (try_bond, try_fgi,              │    │
│  │        try_precursors)                                  │    │
│  │     c) Select best option                               │    │
│  │                                                         │    │
│  │  3. commit_decision() - Commit to tree                  │    │
│  │     → Execute disconnection/FGI                         │    │
│  │     → Forward validation                                │    │
│  │     → Add precursors to queue if not terminal           │    │
│  │                                                         │    │
│  │  4. Update tree and audit state                         │    │
│  │                                                         │    │
│  └─────────────────────────────────────────────────────────┘    │
└──────────────────────────────┬──────────────────────────────────┘
                               │
                               ↓
┌─────────────────────────────────────────────────────────────────┐
│                      TERMINATION PHASE                          │
│  1. Check if all molecules are terminal                         │
│  2. OR queue is empty                                           │
│  3. OR max_steps reached                                        │
└──────────────────────────────┬──────────────────────────────────┘
                               │
                               ↓
┌─────────────────────────────────────────────────────────────────┐
│                      REPORT GENERATION                          │
│  1. finalize() - Mark tree as complete                          │
│  2. Generate forward synthesis report                           │
│  3. Create visualization (HTML/MD)                              │
│  4. Export results to output directory                          │
└──────────────────────────────┬──────────────────────────────────┘
                               │
                               ↓
                        ┌─────────────┐
                        │    End      │
                        └─────────────┘
```

### Detailed Step-by-Step Workflow

#### Phase 1: Initialization

```python
# 1. Create session
cmd = RetroCmd("my_session.json")
cmd.execute("init", {
    "target": "CC(=O)Nc1ccc(O)cc1",  # Paracetamol
    "name": "Paracetamol",
    "max_depth": 15,
    "max_steps": 50,
    "terminal_cs_threshold": 2.5,
})
```

**Internal Operations:**
1. Validate SMILES using RDKit
2. Convert to canonical SMILES
3. Create `RetrosynthesisTree` with target as root
4. Initialize BFS queue: `[(target_smiles, 0)]`
5. Initialize `seen` set for deduplication
6. Compute initial CS score
7. Create JSON session file

#### Phase 2: Planning Loop

##### Step 2.1: prepare_next()

```python
ctx = cmd.execute("next")
```

**Returns:**
- `action`: "awaiting_decision" or "queue_empty" or "auto_terminal"
- `smiles`: Current molecule to analyze
- `depth`: Current depth in tree
- `cs_score`: Complexity score
- `bond_summary`: Default compact window of disconnectable bonds
- `fgi_summary`: Default compact window of available FGI options
- `bond_summary_offset` / `bond_summary_limit`: Current bond window metadata
- `fgi_summary_offset` / `fgi_summary_limit`: Current FGI window metadata
- `bonds_omitted` / `fgi_omitted`: Remaining undisplayed items after the current window
- `warnings`: Functional group conflicts

**Termination Check (3 levels):**
1. **Small molecule**: Heavy atoms ≤ 6 → auto terminal
2. **CS threshold**: CS score ≤ terminal_cs_threshold → auto terminal
3. **No disconnectable bonds**: No viable disconnections → terminal

##### Step 2.2: Exploration (LLM-driven)

```python
# Option A: Explore bond details
cmd.execute("explore", {"bond_idx": 0})
# Returns: Full precursor SMILES for bond alternatives

# Option B: Explore FGI options
cmd.execute("explore_fgi", {})
# Returns: All functional group interconversion options
```

##### Step 2.3: Sandbox Testing

```python
# Test disconnection without committing
cmd.execute("try_bond", {
    "bond_idx": 0,
    "alt_idx": 0,
})
# Returns: precursors, forward_validation, atom_balance

# Test LLM-proposed precursors
cmd.execute("try_precursors", {
    "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
    "reaction_type": "Schotten-Baumann acylation",
})
```

**Sandbox Results Include:**
- Precursor SMILES (canonicalized)
- CS scores for each precursor
- Terminal classification
- Forward validation pass/fail
- Atom balance check
- Cycle warnings

##### Step 2.4: Commit Decision

```python
cmd.execute("commit", {
    "idx": 0,  # Selected sandbox attempt
    "reasoning": "Acyl chloride method, mild conditions",
    "confidence": "high",
})
```

**Commit Operations:**
1. Execute disconnection with selected parameters
2. Run forward validation
3. Check for cycles (duplicate precursors)
4. Create `ReactionNode` in tree
5. Create `MoleculeNode` for each precursor
6. Classify precursors as terminal or intermediate
7. Add intermediates to BFS queue
8. Update audit state

#### Phase 3: Termination

```python
# When queue is empty
cmd.execute("finalize", {
    "summary": "Retrosynthesis complete. 2 steps to commercial materials.",
})
```

#### Phase 4: Report Generation

```python
# Generate forward synthesis report
report = cmd.execute("report", {})

# Export visualizations
cmd.execute("export", {
    "name": "Paracetamol",
    "output_dir": "output/",
})
```

**Outputs:**
- `forward_report.md`: Step-by-step synthesis directions
- `tree_visualization.html`: Interactive tree with embedded molecular images
- `session_data.json`: Complete session export

---

## Data Flow

### JSON Session Persistence

```
┌─────────────────────────────────────────────────────────────┐
│                     session.json                            │
├─────────────────────────────────────────────────────────────┤
│ {                                                           │
│   "session_id": "abc123",                                   │
│   "target": "CC(=O)Nc1ccc(O)cc1",                           │ 
│   "target_name": "Paracetamol",                             │
│   "created_at": "2024-01-01T00:00:00Z",                     │
│   "tree": {                                                 │
│     "molecule_nodes": {...},                                │
│     "reaction_nodes": [...],                                │
│     "status": "in_progress"                                 │
│   },                                                        │
│   "orchestrator": {                                         │
│     "queue": [...],                                         │
│     "seen": [...],                                          │
│     "steps_executed": 2                                     │
│   },                                                        │
│   "sandbox": {                                              │
│     "attempts": [...],                                      │
│     "selected": 0                                           │
│   },                                                        │
│   "audit_state": {...}                                      │
│ }                                                           │
└─────────────────────────────────────────────────────────────┘
```

### Tree Construction Flow

```
Target Molecule (MoleculeNode: role=TARGET, depth=0)
        │
        │ commit_decision()
        ↓
   ReactionNode (step_1)
        │
        ├──→ Precursor 1 (MoleculeNode: role=INTERMEDIATE, depth=1)
        │        │
        │        │ If CS > threshold and has bonds
        │        ↓
        │   Add to BFS queue
        │
        └──→ Precursor 2 (MoleculeNode: role=TERMINAL, depth=1)
                 │
                 │ If CS ≤ threshold or heavy_atoms ≤ 6
                 ↓
            Mark as terminal (no further expansion)
```

### Sandbox vs. Tree

```
┌─────────────────────────────────────────────────────────────┐
│                        SANDBOX                              │
│  (Temporary testing, NOT saved to tree)                     │
│                                                             │
│  try_bond() → SandboxResult:                                │
│    - precursors                                             │
│    - forward_validation                                     │
│    - atom_balance                                           │
│    - cycle_warnings                                         │
│                                                             │
│  Multiple attempts stored in _sandbox_attempts              │
└─────────────────────────────────────────────────────────────┘
                          │
                          │ LLM selects best option
                          ↓
┌─────────────────────────────────────────────────────────────┐
│                          TREE                               │
│  (Persistent storage, saved to session.json)                │
│                                                             │
│  commit() → ReactionNode + MoleculeNodes                    │
│    - Written to tree.tree                                   │
│    - Saved to session.json                                  │
│    - Cannot be undone (except by reload)                    │
└─────────────────────────────────────────────────────────────┘
```

---

## Key Components

### RetroCmd — Command Interface

**File**: [main/retro_cmd.py](main/retro_cmd.py)

23 commands organized by function:

| Category | Commands | Description |
|----------|----------|-------------|
| **Session** | `init` | Create new retrosynthesis session |
| **Navigation** | `next`, `context` | Get next molecule, view context |
| **Exploration** | `explore`, `explore_fgi` | Examine bond/FGI details |
| **Sandbox** | `try_bond`, `try_fgi`, `try_precursors` | Test without committing |
| **Sandbox Mgmt** | `sandbox_list`, `sandbox_clear`, `select` | Manage sandbox state |
| **Decisions** | `commit`, `accept`, `skip` | Commit, mark terminal, skip |
| **Viewing** | `tree`, `status` | Display tree, status |
| **Completion** | `finalize`, `report`, `export` | Complete, generate reports |
| **Advanced** | `smart_cap`, `custom_cap` | Capping group assistance |

### RetrosynthesisOrchestrator — BFS Engine

**File**: [main/retro_orchestrator.py](main/retro_orchestrator.py)

**Key Methods:**

| Method | Purpose |
|--------|---------|
| `prepare_next()` | Pop from queue, analyze, return context |
| `explore_bond(idx)` | Expand bond details on demand |
| `explore_fgi()` | Expand FGI options |
| `try_disconnection()` | Sandbox test bond break |
| `try_fgi()` | Sandbox test FGI |
| `try_precursors()` | Sandbox test LLM-proposed precursors |
| `commit_decision()` | Execute disconnection, write to tree |
| `accept_terminal()` | Mark molecule as terminal |
| `finalize()` | Complete planning, return final state |

**Safety Parameters:**
- `max_depth`: Maximum tree depth (default: 15)
- `max_steps`: Maximum reaction steps (default: 50)
- `max_queue_size`: Maximum pending molecules (default: 200)
- `terminal_cs_threshold`: CS score for auto-terminal (default: 2.5)

### RetrosynthesisTree — Data Model

**File**: [main/retro_tree.py](main/retro_tree.py)

**Dual-Type Node Graph:**

```
MoleculeNode                 ReactionNode
├── smiles                   ├── step_id
├── node_id                  ├── depth
├── role (target/intermediate/terminal)  ├── reaction_smiles
├── depth                    ├── product_node
├── complexity (CS score)    ├── reactant_nodes
├── decision_context         ├── reaction_type
└── llm_analysis             ├── template_evidence
                             ├── llm_decision
                             └── forward_validation
```

**Key Features:**
- Canonical SMILES deduplication
- Depth tracking
- Role classification
- Complete JSON serialization
- Audit trail support

### Chemical Tools — 8 Modules, 26 Functions

**Directory**: [chem_tools/](chem_tools/)

| Module | ID | Functions | Purpose |
|--------|-----|-----------|---------|
| `_rdkit_utils` | M0 | 7 | RDKit utilities (internal) |
| `mol_info` | M1 | 2 | Molecular properties, scaffold matching |
| `fg_detect` | M2 | 4 | Functional group detection |
| `template_scan` | M3 | 2 | Reaction template scanning |
| `bond_break` | M4 | 6 | Bond disconnection execution |
| `forward_validate` | M5 | 2 | Forward reaction validation |
| `cs_score` | M6 | 4 | Complexity scoring |
| `fg_warnings` | M7 | 4 | Functional group conflict detection |
| `smart_cap` | M8 | 2 | Intelligent capping suggestions |

---

## Decision Context Layering

### Context Levels

```
┌─────────────────────────────────────────────────────────────┐
│                    FULL (detail="full")                     │
│  ─────────────────────────────────────────────────────────  │
│  - All disconnectable bonds with ALL precursor SMILES       │
│  - All FGI options with full details                        │
│  - Complete functional groups list                          │
│  - Full molecule properties                                 │
│  - Use for: Final decision making                           │
│  - Token cost: HIGH                                         │
└─────────────────────────────────────────────────────────────┘
                          │ (most detail)
                          ↓
┌─────────────────────────────────────────────────────────────┐
│                   COMPACT (default)                         │
│  ─────────────────────────────────────────────────────────  │
│  - Bond summary (atoms, types, reaction types, counts)      │
│  - NO precursor SMILES (use explore_bond to get them)       │
│  - FGI summary (template names only)                        │
│  - Default first window only (bond_limit=5, fgi_limit=5)    │
│  - Omitted counts + window metadata for later paging        │
│  - Use for: Initial LLM review, token efficiency            │
│  - Token cost: LOW                                          │
└─────────────────────────────────────────────────────────────┘
                          │
                          ↓
┌─────────────────────────────────────────────────────────────┐
│                    STATUS (detail="status")                 │
│  ─────────────────────────────────────────────────────────  │
│  - Queue state, step count, tree depth                      │
│  - Molecule counts (pending/terminal/intermediate)          │
│  - Elapsed time                                             │
│  - Use for: Progress monitoring                             │
│  - Token cost: MINIMAL                                      │
└─────────────────────────────────────────────────────────────┘
                          │
                          ↓
┌─────────────────────────────────────────────────────────────┐
│                     TREE (detail="tree")                    │
│  ─────────────────────────────────────────────────────────  │
│  - ASCII tree visualization                                 │
│  - Terminal/pending molecule lists                          │
│  - Use for: Human-readable overview                         │
│  - Token cost: MEDIUM                                       │
└─────────────────────────────────────────────────────────────┘
```

### On-Demand Expansion Pattern

```
1. Initial: next() → compact context (bond summaries)
                    ↓
2. LLM interested in bond 2 → explore_bond(2)
                    ↓
3. Get: Full precursor SMILES for bond 2 alternatives
                    ↓
4. LLM tries: try_bond(bond_idx=2, alt_idx=0)
                    ↓
5. Get: Sandbox result with validation
                    ↓
6. LLM satisfied: commit(idx=0)
                    ↓
7. Written to tree
```

---

## Chemical Analysis Pipeline

### Module Flow Diagram

```
                    ┌──────────────┐
                    │   Input:     │
                    │  SMILES      │
                    └──────┬───────┘
                           │
                           ↓
                    ┌──────────────┐
                    │  M0: Utils   │
                    │ parse_mol()  │
                    │ canonical()  │
                    └──────┬───────┘
                           │
         ┌─────────────────┼─────────────────┐
         │                 │                 │
         ↓                 ↓                 ↓
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│  M1: Info    │  │  M2: FG Det. │  │  M6: CS Score│
│analyze_mol() │  │detect_fgs()  │  │compute_cs()  │
└──────┬───────┘  └──────┬───────┘  └──────┬───────┘
       │                 │                 │
       └─────────────────┼─────────────────┘
                         │
                         ↓
                  ┌──────────────┐
                  │  M3: Templ.  │
                  │ scan_rxns()  │
                  └──────┬───────┘
                         │
              ┌──────────┴──────────┐
              │                     │
              ↓                     ↓
     ┌──────────────┐      ┌──────────────┐
     │  M4: Break   │      │  M5: Valid.  │
     │disconnect()  │      │validate_fwd()│
     └──────┬───────┘      └──────┬───────┘
            │                     │
            └──────────┬──────────┘
                       │
                       ↓
                ┌──────────────┐
                │  M7: Warn    │
                │check_conf.() │
                └──────┬───────┘
                       │
                       ↓
                ┌──────────────┐
                │  M8: SmartCap│
                │suggest_cap() │
                └──────────────┘
```

### M0: RDKit Utilities (Internal)

| Function | Description |
|----------|-------------|
| `parse_mol(smiles)` | Parse SMILES to RDKit Mol object |
| `canonical(smiles)` | Convert to canonical SMILES |
| `validate_smiles(smiles)` | 5-step validation |
| `smarts_match(mol, smarts)` | SMARTS pattern matching |
| `tanimoto(mol_a, mol_b)` | Fingerprint similarity |

### M1: Molecular Information

```python
from Rachel.chem_tools import analyze_molecule

result = analyze_molecule("CC(=O)Nc1ccc(O)cc1")
# Returns:
{
    "smiles": "CC(=O)Nc1ccc(O)cc1",
    "formula": "C8H9NO2",
    "mw": 151.16,
    "heavy_atoms": 10,
    "rings": [{"aromatic": True, "size": 6}],
    "stereo": {"chiral_centers": 0},
    "scaffold": "c1ccc(cc1)"  # Murcko scaffold
}
```

### M2: Functional Group Detection

```python
from Rachel.chem_tools import detect_functional_groups

result = detect_functional_groups("CC(=O)Nc1ccc(O)cc1")
# Returns:
[
    {"name": "amide", "smarts": "[NX3][CX3](=[OX1])", "count": 1, "atoms": [2,3,4]},
    {"name": "phenol", "smarts": "[OX2H][cX3]:[c]", "count": 1, "atoms": [8,9]}
]
```

### M3: Template Scanning

```python
from Rachel.chem_tools import find_disconnectable_bonds

result = find_disconnectable_bonds("CC(=O)Nc1ccc(O)cc1")
# Returns:
{
    "bonds": [
        {
            "atoms": [2, 5],
            "bond_type": "SINGLE",
            "in_ring": False,
            "heuristic_score": 0.85,
            "alternatives": [
                {"template_id": "amidation_retro", "template": "Amide bond formation"}
            ]
        }
    ]
}
```

### M4: Bond Disconnection

```python
from Rachel.chem_tools import execute_disconnection

result = execute_disconnection(
    smiles="CC(=O)Nc1ccc(O)cc1",
    bond=(2, 5),
    reaction_type="amidation"
)
# Returns:
{
    "success": True,
    "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
    "operation_log": [...]
}
```

### M5: Forward Validation

5-step validation process:

```
┌─────────────────────────────────────────────────────────────┐
│              Forward Validation Pipeline                    │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  1. Atom Balance Check (25% weight)                         │
│     → Hard fail: skeleton_imbalance, severe_imbalance       │
│     → Score: balance_score (0-1)                            │
│                                                             │
│  2. Template Forward Execution (25% weight)                 │
│     → Execute forward template if available                 │
│     → Check if products match target                        │
│                                                             │
│  3. MCS Skeleton Alignment (20% weight)                     │
│     → Maximum common substructure between                   │
│       precursors+byproducts and target                      │
│     → Hard fail: aligned=False                              │
│                                                             │
│  4. Bond Change Topology (15% weight)                       │
│     → Check bond order changes make sense                   │
│     → Hard fail: has_hard_fail                              │
│                                                             │
│  5. Functional Group Compatibility (15% weight)             │
│     → Check FG compatibility with reaction conditions       │
│     → Hard fail: forbidden FG present                       │
│                                                             │
│  → Overall feasibility_score (0-1)                          │
│  → pass/fail based on threshold                             │
└─────────────────────────────────────────────────────────────┘
```

### M6: CS Complexity Score

6-dimensional scoring:

| Dimension | Weight | Description |
|-----------|--------|-------------|
| size | 0.55 | Heavy atom count (normalized) |
| ring | 0.65 | Ring complexity (bridged, fused, etc.) |
| stereo | 0.55 | Stereochemical complexity |
| hetero | 0.40 | Heteroatom fraction |
| symmetry | -0.20 | Symmetry reduces difficulty |
| fg_density | 0.35 | Functional group density |

**Range**: 1.0 (trivial) to 10.0 (very complex)

**Classification**:
- ≤ 2.5: trivial (commercial materials)
- ≤ 6.0: moderate (synthetically accessible)
- > 6.0: complex (challenging synthesis)

### M7: Functional Group Warnings

```python
from Rachel.chem_tools import check_fg_conflicts

result = check_fg_conflicts("CC(=O)Nc1ccc(O)cc1")
# Returns:
{
    "selectivity_conflicts": [],
    "competitive_groups": [],
    "dangerous_combos": [],
    "protection_suggestions": [
        {"fg": "phenol", "suggest": "protect as acetate if needed"}
    ]
}
```

### M8: Smart Capping

```python
from Rachel.chem_tools import suggest_capping

result = suggest_capping("CC(=O)Nc1ccc(O)cc1", bond=(2, 5))
# Returns:
{
    "ok": True,
    "proposals": [
        {
            "reaction_type": "amidation",
            "cap_i": "Cl",
            "cap_j": "H",
            "fragments": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
            "confidence": 0.95,
            "description": "Standard amidation with acyl chloride"
        }
    ]
}
```

---

## Complete Example: Paracetamol Synthesis

### Target

```
SMILES: CC(=O)Nc1ccc(O)cc1
Name: Paracetamol (Acetaminophen)
```

### Step-by-Step Execution

#### Step 1: Initialize Session

```python
from Rachel.main import RetroCmd

cmd = RetroCmd("paracetamol_session.json")

result = cmd.execute("init", {
    "target": "CC(=O)Nc1ccc(O)cc1",
    "name": "Paracetamol",
    "terminal_cs_threshold": 2.5,
})

# Response:
{
    "ok": True,
    "session_id": "paracetamol_001",
    "target": "CC(=O)Nc1ccc(O)cc1",
    "hint": "调用 next 获取第一个分子的上下文"
}
```

#### Step 2: Get First Molecule Context

```python
ctx = cmd.execute("next")

# Response:
{
    "action": "awaiting_decision",
    "decision_tier": "standard",
    "smiles": "CC(=O)Nc1ccc(O)cc1",
    "depth": 0,
    "cs_score": 3.82,
    "classification": "moderate",
    "is_terminal": False,
    "n_bonds": 9,
    "n_fgi": 3,
    "bond_summary": [
        {
            "bond_idx": 0,
            "actual_bond_idx": 7,
            "atoms": [2, 5],
            "role_pair": ["aryl_sp2_c", "carbonyl_c"],
            "bond_type": "SINGLE",
            "heuristic_score": 0.92,
            "n_alternatives": 2,
            "reaction_types": ["Amide bond formation", "Amide bond formation (alt)"]
        }
    ],
    "bond_summary_offset": 0,
    "bond_summary_limit": 5,
    "bonds_omitted": 4,
    "fgi_summary": [],
    "fgi_summary_offset": 0,
    "fgi_summary_limit": 5,
    "hint": "这是精简视图。用 explore_bond(idx) 查看某键位的完整前体方案"
}
```

如果 `bonds_omitted > 0` 或 `fgi_omitted > 0`，可以显式取下一窗 summary：

```python
ctx_page_2 = cmd.execute(
    "context",
    {
        "detail": "compact",
        "bond_offset": 5,
        "bond_limit": 5,
        "fgi_offset": 5,
        "fgi_limit": 5,
    },
)
```

注意：
- `bond_summary` / `fgi_summary` 里的索引始终是全局索引
- `explore(bond_idx)` 与 `try_bond(bond_idx, alt_idx)` 不受分页影响
- `session.json` 中 `current` 复用同一 default compact 生成逻辑，与 live compact 保持一致

#### Step 3: Explore Bond Options

```python
details = cmd.execute("explore", {"bond_idx": 0})

# Response:
{
    "bond_idx": 0,
    "atoms": [2, 5],
    "bond_type": "SINGLE",
    "heuristic_score": 0.92,
    "alternatives": [
        {
            "template_id": "amidation_retro",
            "template": "Amide bond formation",
            "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
            "confidence": 0.95
        },
        {
            "template_id": "amidation_anhydride_retro",
            "template": "Amide bond from anhydride",
            "precursors": ["CC(=O)OC(=O)C", "Nc1ccc(O)cc1"],
            "confidence": 0.85
        }
    ],
    "smart_capping": [...]
}
```

#### Step 4: Sandbox Test

```python
test = cmd.execute("try_bond", {"bond_idx": 0, "alt_idx": 0})

# Response:
{
    "success": True,
    "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
    "precursor_details": [
        {"smiles": "CC(=O)Cl", "cs_score": 1.8, "is_terminal": True, "heavy_atoms": 3},
        {"smiles": "Nc1ccc(O)cc1", "cs_score": 2.1, "is_terminal": True, "heavy_atoms": 7}
    ],
    "forward_validation": {
        "ok": True,
        "pass": True,
        "feasibility_score": 0.92
    },
    "atom_balance": {"balanced": True, "balance_score": 1.0},
    "hint": "这是沙盒结果，未写入树。满意请调 commit_decision()。"
}
```

#### Step 5: Commit Decision

```python
commit = cmd.execute("commit", {
    "idx": 0,
    "reasoning": "Acyl chloride amidation is efficient. Both precursors are commercially available.",
    "confidence": "high",
})

# Response:
{
    "success": True,
    "step_id": "rxn_1",
    "new_pending": [],
    "new_terminal": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
    "tree_complete": True
}
```

#### Step 6: Check Tree Status

```python
status = cmd.execute("tree")

# Response:
{
    "tree": """
🎯 CC(=O)Nc1ccc(O)cc1 CS=3.8 [mol_0]
  ↓ Amide bond formation
    ✓ CC(=O)Cl CS=1.8 [mol_1]
    ✓ Nc1ccc(O)cc1 CS=2.1 [mol_2]
    """,
    "terminal_count": 2,
    "pending_count": 0
}
```

#### Step 7: Finalize and Export

```python
# Finalize
cmd.execute("finalize", {
    "summary": "One-step retrosynthesis complete. Both precursors are commercially available."
})

# Generate report
report = cmd.execute("report", {})

# Export visualization
cmd.execute("export", {"name": "Paracetamol"})

# Creates:
# - output/Paracetamol_forward_report.md
# - output/Paracetamol_tree_visualization.html
# - output/Paracetamol_session.json
```

### Final Forward Synthesis Report

```markdown
# Forward Synthesis Report: Paracetamol

## Target
- **SMILES**: CC(=O)Nc1ccc(O)cc1
- **Name**: Paracetamol
- **Complexity**: CS=3.82 (moderate)

## Synthetic Route (1 step)

### Step 1: Amide bond formation
- **Reactants**: Acetyl chloride + 4-aminophenol
- **Product**: Paracetamol
- **Confidence**: high
- **Reasoning**: Acyl chloride amidation is efficient. Both precursors are commercially available.

## Starting Materials (Commercial)
1. Acetyl chloride (CC(=O)Cl) - CS=1.8, trivial
2. 4-Aminophenol (Nc1ccc(O)cc1) - CS=2.1, trivial
```

---

## API Reference

### Command Quick Reference

| Command | Arguments | Returns | Description |
|---------|-----------|---------|-------------|
| `init` | target, name, max_depth, max_steps, terminal_cs_threshold | session info | Create new session |
| `next` | none | context or queue_empty | Get next molecule |
| `context` | detail (compact/full/status/tree), bond_offset, bond_limit, fgi_offset, fgi_limit | context | Get current context or later compact windows |
| `explore` | bond_idx | bond details | Expand bond info |
| `explore_fgi` | none | FGI options | Expand FGI options |
| `try_bond` | bond_idx, alt_idx | sandbox result | Test disconnection |
| `try_fgi` | fgi_idx | sandbox result | Test FGI |
| `try_precursors` | precursors (list), reaction_type | sandbox result | Test custom precursors |
| `sandbox_list` | none | attempt history | View sandbox |
| `sandbox_clear` | none | ok | Clear sandbox |
| `select` | idx | selection info | Select attempt |
| `commit` | idx, reasoning, confidence | commit result | Write to tree |
| `accept` | reason | ok | Mark as terminal |
| `skip` | reason | ok | Skip molecule |
| `tree` | none | tree text | Display tree |
| `status` | none | status info | Show status |
| `finalize` | summary | final state | Complete planning |
| `report` | none | report text | Generate report |
| `export` | name, output_dir | file paths | Export results |
| `smart_cap` | bond_idx or (smiles, bond) | capping proposals | Get capping suggestions |
| `custom_cap` | cap_i, cap_j, bond_idx or (smiles, bond) | custom result | Custom capping |

### Chemical Tools Quick Reference

| Module | Function | Signature | Returns |
|--------|----------|-----------|---------|
| M1 | `analyze_molecule` | `(smiles: str)` | Dict with molecular properties |
| M1 | `match_known_scaffolds` | `(smiles: str, max_hits: int)` | List of scaffold matches |
| M2 | `detect_functional_groups` | `(smiles: str)` | List of FG matches |
| M2 | `detect_reactive_sites` | `(smiles: str)` | List of reactive sites |
| M2 | `get_fg_reaction_mapping` | `(groups: List)` | FG→reaction mapping |
| M3 | `find_disconnectable_bonds` | `(smiles: str)` | Disconnectable bonds with alternatives |
| M3 | `scan_applicable_reactions` | `(smiles: str, mode: str)` | Applicable reaction templates |
| M4 | `execute_disconnection` | `(smiles, bond, reaction_type)` | BreakResult with precursors |
| M4 | `execute_fgi` | `(smiles, template_id)` | BreakResult for FGI |
| M4 | `preview_disconnections` | `(smiles: str)` | All possible disconnections |
| M5 | `validate_forward` | `(precursors, target, ...)` | Validation result |
| M5 | `check_atom_balance` | `(precursors, target)` | Balance check |
| M6 | `compute_cs_score` | `(smiles: str)` | CS score (1-10) |
| M6 | `classify_complexity` | `(smiles: str)` | Classification string |
| M7 | `check_fg_conflicts` | `(smiles: str)` | Conflict warnings |
| M7 | `suggest_protection_needs` | `(smiles, planned_reaction)` | Protection suggestions |
| M8 | `suggest_capping` | `(smiles, bond, ...)` | Capping proposals |
| M8 | `custom_cap` | `(smiles, bond, cap_i, cap_j)` | Custom capping result |

---

## References

- **[skill.md](skill.md)** — LLM operation manual
- **[readme.md](readme.md)** — Project overview
- **[refs.md](refs.md)** — Technical reference
- **[chem_tools/README.md](chem_tools/README.md)** — Chemical tools documentation

---

*Document generated for Rachel v1.0 - Multi-step Retrosynthetic Planning Engine*
