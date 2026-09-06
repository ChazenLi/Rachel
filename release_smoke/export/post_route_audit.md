# Rachel Post-Finalize Audit

- audit id: `audit_20260906T073900.445347+0000_1edc3e0d`
- route digest: `sha256:5d24616dabb6782adecbd7b7e3eeb5906ee8c2567a1b92804599827e1edc3e0d`
- generated at: 2026-09-06T07:39:00.445347+00:00
- status: **incomplete**
- incompleteness reasons: substance_safety_lookup_incomplete, offline_mode

## Buyability

PubChem B1/B2 evidence is identity/catalog evidence only; it is not live stock, price, lead time, or procurement confirmation.

| node | SMILES | grade | CID | vendor catalog | example vendor records |
|---|---|---:|---:|---|---|
| mol_0 | `CNC(C)=O` | B0 | - | not_queried | - |

## Substance Safety Evidence

Each GHS classification below remains attached to its PubChem depositor/authority source and license; conflicting sources are not collapsed into a single verdict.

| material | roles | CID | grade | checked-source state | local structural alerts | source-attributed GHS evidence |
|---|---|---:|---:|---|---|---|
| `CNC(C)=O` | terminal | - | S0 | incomplete | - | - |

## Step Safety Screening

The explicit assessment is LLM-authored planning judgment; its semantic correctness is not independently validated.

| step | reaction | explicit assessment | local alerts | condition gaps | process status |
|---|---|---|---|---|---|

## Limitations

- This audit does not change route chemistry, terminal decisions, validation gates, or route ranking.
- No result in this payload is a generic safe/unsafe verdict.
- B1 proves a PubChem identity record, not commercial availability.
- B2 proves depositor/vendor catalog rows, not live stock, price, lead time, region, purity, pack size, or procurement.
- No hit means no evidence in the checked PubChem sources; it does not prove that a material cannot be purchased.
- GHS classifications are source-specific substance evidence and may differ across authorities or depositors.
- Supplier- and product-specific SDS records were not retrieved in this version.
- Structured condition coverage checks declared field presence and shape only; it does not validate chemical correctness or establish process safety.
- Even complete declared condition coverage supports only preliminary review and does not replace SDS, SOP, calorimetry, exposure, scale-up, or equipment assessment.
- S0 structural and combination alerts are screening hypotheses, not process hazard conclusions.
- Absence of a checked-source classification or alert does not mean safe.
