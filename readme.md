<div align="center">

# Rachel

**Formalized chemical reasoning for multi-step retrosynthesis**

<img alt="Python 3.10+" src="https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white">
<img alt="Active Research" src="https://img.shields.io/badge/Status-Active%20Research-2D6A4F">
<img alt="Multi-Step Retrosynthesis" src="https://img.shields.io/badge/Domain-Multi--Step%20Retrosynthesis-8C564B">
<img alt="Workflow" src="https://img.shields.io/badge/Workflow-State--Action--Validation--Commit-7B61FF">
<img alt="LLM Strategy Layer" src="https://img.shields.io/badge/LLM-Strategy%20Layer-6F42C1">

<p>
  <a href="#trace-demo">Trace Demo</a> |
  <a href="#end-to-end-example">End-to-End Example</a> |
  <a href="#highlights">Highlights</a> |
  <a href="#selected-molecules">Selected Molecules</a> |
  <a href="#minimal-quickstart">Quickstart</a>
</p>

<a href="assets/media/rachel-trace.mp4">
  <img src="assets/media/rachel-trace-poster.png" alt="Rachel trace walkthrough" width="920">
</a>

</div>

Rachel treats retrosynthetic planning as a structured `state -> action -> validation -> commit` process rather than a one-shot text generation task. The current repository is an active research codebase being cleaned up for arXiv-facing presentation while remaining in day-to-day use.

## Trace Demo

The trace above is a visual walkthrough of the Rachel planning workflow from structured state to a committed retrosynthetic route.

- Click the hero image to open the [MP4 trace](assets/media/rachel-trace.mp4)
- The trace is meant to show planning behavior, not only final route output
- It is the fastest way to see how Rachel moves from context to sandboxed candidates to validated commitment

## End-to-End Example

The figure below shows a full route-level comparison between a PaRoutes ground-truth plan and Rachel's generated result on case `n1_366`.

![Ground truth vs Rachel on PaRoutes n1_366](assets/cases/paroutes-n1-366-groundtruth-vs-rachel.png)

This example is included as a qualitative systems-level reference point. The point is not only whether a single step looks chemically plausible, but whether the route remains interpretable and structurally coherent at the full planning level.

## Highlights

| Capability | What it means in Rachel |
| --- | --- |
| Stateful planning | Rachel reasons over persistent session state instead of isolated one-shot answers. |
| Grounded operator space | Bond disconnection and FGI are treated as complementary planning operators. |
| Sandbox before commitment | Candidate steps are tried locally before they affect the main route tree. |
| Validation-gated execution | Feasibility checks, atom balance, and related validators help control commitment. |
| Structured route memory | Accepted steps become explicit route-tree objects rather than only free-form text. |
| LLM as strategy layer | The LLM helps organize search and choice, rather than acting as an unchecked chemistry oracle. |

## What Rachel Is

Rachel combines:

- a session-driven planning workflow
- chemistry-grounded operators such as bond disconnection and FGI
- sandboxed local trials before route commitment
- validator-gated route construction
- route-tree memory and audit state

The goal is not just to propose a route, but to make route construction inspectable, recoverable, and chemically checkable.

## Core Workflow

```mermaid
flowchart LR
    A["Current State"] --> B["Candidate Action"]
    B --> C["Sandbox Trial"]
    C --> D["Validation"]
    D --> E["Commit or Reject"]
    E --> F["Updated Route Tree"]
```

Candidate actions are explored before they are written into the route. Validated steps move forward into the main tree; rejected attempts remain informative planning artifacts rather than disappearing into free-form text.

## Selected Molecules

Rachel is currently showcased with three qualitative examples chosen to cover complementary strengths.

<table>
  <tr>
    <td align="center" width="33%">
      <img src="assets/showcases/qntr/mol_0.png" alt="QNTR" width="220"><br>
      <strong>QNTR</strong><br>
      <sub>Experimentally grounded example</sub>
    </td>
    <td align="center" width="33%">
      <img src="assets/showcases/losartan/mol_0.png" alt="Losartan" width="220"><br>
      <strong>Losartan</strong><br>
      <sub>Canonical medicinal chemistry target</sub>
    </td>
    <td align="center" width="33%">
      <img src="assets/showcases/rivaroxaban/mol_0.png" alt="Rivaroxaban" width="220"><br>
      <strong>Rivaroxaban</strong><br>
      <sub>Deeper drug-like example</sub>
    </td>
  </tr>
</table>

| Molecule | Role | Route depth | What it highlights |
| --- | --- | ---: | --- |
| `QNTR` | Experimentally grounded example | 2 steps | A short, interpretable route tied to real synthesis experience |
| `Losartan` | Canonical medicinal chemistry target | 4 steps | Convergent route logic with recognizable medicinal chemistry disconnections |
| `Rivaroxaban` | Deeper drug-like example | 5 steps | Longer-horizon planning with a broader transformation mix |

### QNTR

An experimentally grounded molecule connected to your own synthesis experience.

- Short 2-step route from 3 starting materials
- Useful as a compact demonstration of interpretable planning
- Especially valuable because the accepted route reflects strategy beyond simple template satisfaction

### Losartan

A canonical medicinal chemistry target with a recognizable convergent route.

- 4-step route from 4 starting materials
- Highlights tetrazole formation, N-alkylation, and Suzuki coupling logic
- Useful as a benchmark-like example that many readers can immediately interpret

### Rivaroxaban

A deeper drug-like example with a richer transformation mix.

- 5-step route from 4 starting materials
- Highlights Buchwald-Hartwig amination, FGI, cyclization, and amide formation
- Useful for showing that Rachel is not limited to short or purely toy routes

### Dual Drug-Case Comparison

The figure below places the Losartan and Rivaroxaban examples into one annotated comparison view. It gives the README a more complete qualitative picture of how Rachel behaves on two recognizable drug-like targets with different route depths and transformation profiles.

![Dual annotated comparison for Losartan and Rivaroxaban](assets/cases/rivaroxaban-losartan-dual-annotated-en.png)

Taken together, these two cases make the contrast especially clear:

- `Losartan` emphasizes classical convergent medicinal chemistry logic
- `Rivaroxaban` emphasizes deeper route depth and operator diversity
- The pair helps readers compare route style, not just isolated outcomes

See the full qualitative page in [showcases](docs/showcases.md).

## Minimal Quickstart

Current local runs assume a Python environment with the main research dependencies already available, including Python 3.10+, RDKit, `numpy`, and `Pillow`.

```python
from Rachel.main import RetroCmd

cmd = RetroCmd("my_session.json")

cmd.execute(
    "init",
    {
        "target": "CC(=O)Nc1ccc(O)cc1",
        "name": "Paracetamol",
        "terminal_cs_threshold": 1.5,
    },
)

ctx = cmd.execute("next")

cmd.execute(
    "try_precursors",
    {
        "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
        "reaction_type": "Schotten-Baumann acylation",
    },
)

cmd.execute(
    "commit",
    {
        "idx": 0,
        "reasoning": "Acylation with simple, accessible precursors.",
        "confidence": "high",
    },
)
```

This is a protocol-level example, not a full benchmark workflow. More technical notes are preserved in [usage notes](docs/usage-notes.md).

<details>
<summary><strong>Repository Map</strong></summary>

- [main](main): orchestration, session logic, route tree, reports, and command interface
- [chem_tools](chem_tools): chemistry-grounded operators and validation utilities
- [tools](tools): helper scripts for runs, analysis, visualization, and related research workflows
- [tests](tests): current validation and experiment-support material
- [plan](plan): manuscript drafts, writing materials, and paper-preparation assets

</details>

## Project Status

- Active research codebase
- Currently being prepared for arXiv-facing presentation
- Documentation is being cleaned up, but the repository remains under active use
- Core workflow is already in use
- Not yet a fully hardened OSS release
