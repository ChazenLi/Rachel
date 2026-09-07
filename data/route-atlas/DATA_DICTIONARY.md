# Route Browser Data Dictionary

`data/PaRoutes120.json` and `data/RF25.json` use schema `rachel.route-browser.v1`.
IDs are local to a dataset or a route, as specified below. Source SMILES retain
stereochemical notation; drawing generation is not chemical validation.

| Field | Meaning |
| --- | --- |
| `dataset`, `generated`, `defaultVariant` | Dataset, assembly date, default route version (`final`) |
| `methods[]` | Stable comparison IDs and display names; `kind: reference` distinguishes the PaRoutes comparator from the eight entries with `kind: method` |
| `targets[]` | Complete frozen target list: `id`, `name`, `split`, `smiles`, depiction key and metadata |
| `targets[].metadata` | RF25 publication registry and molecular descriptors; empty for PaRoutes120 |
| `routes[target_id][method]` | One target-method position, retained even without a graph |
| `endpoint` | Adopted endpoint row, with local paths removed |
| `strict_closed` | Adopted manuscript target-method endpoint; null for Reference, which has no assigned method endpoint |
| `status` | `closed`, `open`, `no_route` or `excluded` for method outputs; `reference` for the PaRoutes comparator |
| `record_kind` | `reference` on reference entries; these are excluded from method-outcome totals |
| `selection`, `candidate_index` | Recorded selection policy and original candidate index when applicable |
| `sources[]` | Upstream file name and actual byte SHA-256; optional recorded hash and chemical reconciliation basis; reference sources also include `package_path` to the original subset/manifest files |
| `graph` | Selected route graph; null when no adopted route is available |
| `graph.roots` | Target molecule node ID |
| `graph.nodes[]` | Route-local ID, source SMILES, original node/occurrence identity, terminal/target flags, depiction key, original metadata and audit indices |
| `graph.reactions[]` | Route-local ID, original step identity, product node, reactant node list, forward reaction SMILES, recorded label/rationale/gate, original metadata and matching evaluation records |
| `graph.construction` | Whether connectivity preserves native identity, tree occurrences or joins a flat step list by canonical molecular identity |
| `graph.target_matches_input` | Stereochemistry-aware canonical root/input agreement, checked for every route |
| `audits[]` | Preserved method/target terminal audit rows; node audit indices address this array |
| `reviews[]` | Whole-route reviewer, case ID, route hash, four original scores/rationales, decision and source-file references |
| `reviewed_variant` | Separately preserved review input graph for the two chemically different Rachel versions; includes its own matching reviews and route hash |

When the evaluated variant is active, single-route and comparison exports carry
`activeVariant: reviewed`, its graph, its source references and its reviews.
They have an empty molecule-audit array. The `endpoint` remains the separately
reported manuscript endpoint. Full dataset exports retain both versions.

`terminal` means that no accepted/displayed reaction expands that molecule in the
route graph. It does not by itself mean buyable, audited or chemically validated.
`closed` and `open` come from the frozen endpoint table rather than native engine
success flags or visually inferred leaf status.

Reaction evaluations preserve the recorded model fields, including failures,
missing conditions, unresolved reagent roles, predictions and target ranks.
Duplicate source evaluation rows are retained when present; no averages or new
scores are calculated by the viewer. The whole-route review material is the
adopted PaRoutes120 common58 comparison. The viewer does not invent RF25 reviews
or extend common58 reviews to other targets.

PaRoutes120 has 1,080 display positions: 960 method outputs and 120 references.
The combined companion has 1,180 positions, including RF25's 100 method outputs.
Reference has 426 reaction occurrences and preserves original reaction `ID`,
`rsmi` and other metadata. Source-tree `in_stock` values remain original metadata;
they are not converted to independent terminal-audit conclusions.

`SOURCE_CATALOG.csv` maps a SHA-256 to the original file name and byte length.
Its hashes identify upstream bytes; they are different from the package checksums
of normalized exports. Original filesystem paths remain in the internal build
ledger only. A source name alone is not a unique key. Use the full SHA-256.

`CHECKSUMS.csv` covers the released companion files, excluding itself. The sibling
ZIP checksum identifies the assembled archive. Neither hash set claims recovery
of upstream files whose recorded hashes differ from the available bytes.
