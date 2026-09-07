# Rachel Route Atlas

Offline route browsers for the manuscript's main method comparisons.

Open `PaRoutes120.html` or `RF25.html` directly in a current desktop browser.
No installation, server, network connection or API key is needed. Keep the two
HTML files together to use the dataset navigation links. Each file also opens
independently. Extract the ZIP before opening an HTML file.

## Scope

| Dataset | Targets | Main methods and reference | Display positions |
| --- | ---: | ---: | ---: |
| PaRoutes120 | 120 | 8 methods + Reference | 1,080 |
| RF25 | 25 | 4 | 100 |

PaRoutes120 methods: Rachel, Direct LLM, AOT*, SyntheLite,
PaRoutesModel + Retro*, RootAligned + Retro*, RootAligned + MCTS and
LocalRetro + Retro*. PaRoutes120 also includes its 120 dataset reference routes
as an independent `Reference` comparison entry. RF25 includes the first four methods.

The pages retain unsuccessful and unavailable outcomes. They contain 827
PaRoutes120 method route graphs, 120 reference graphs and 77 RF25 route graphs, plus two separately labelled
evaluated variants. A position without recorded reaction steps remains visible;
it is not replaced by another run or another method's route. PaRoutes120 Direct
LLM `n5_7593` remains excluded under the adopted duplicate-generation rule.
Mechanical controls, F0 and F1-F4 are outside this browser's scope.

## Viewing and Export

- Search by target identifier, name or SMILES; PaRoutes120 also has n1/n5 filters.
- Select methods for side-by-side comparison, with up to four panels.
- Switch between route graphs and forward-reading reaction steps.
- Pan, zoom, fit or expand a route. Click structures or numbered reactions for details.
- Route records provide selection provenance, file hashes and matched whole-route reviews.
- Molecule and reaction details retain recorded terminal audits, decision rationales and
  matched forward-model outputs where these exist in the adopted evidence.
- Export one route, the displayed comparison or the complete normalized dataset as JSON.
- Switch Chinese/English interface labels. Original model rationales are preserved in
  their recorded language.

PaRoutes120 opens with Rachel and Reference side by side; RF25 opens with Rachel
and Direct LLM. The default method version is the selected final route or frozen process tree.
Aggregate performance tables and percentages are intentionally absent; the percentage
beside the graph controls is the drawing zoom level.

## Route Selection and Evidence

### PaRoutes Reference

All 120 reference trees come from the adopted PaRoutes subsets of Zenodo record
[7341155](https://doi.org/10.5281/zenodo.7341155): 51 n1 targets and 69 n5 targets.
Both subset files and both target manifests exactly match the SHA-256 values in
the frozen reference evaluator inputs. Their 426 reaction occurrences match those
inputs, and all target roots match the current manuscript cohort.

The original subset and manifest JSON files are included byte-for-byte in
`data/reference/`. Reference graphs preserve the source molecule occurrences,
reaction IDs, atom-mapped reaction SMILES and recorded reagent/condition fields
in their original metadata. Their matched forward-model evidence is available in
reaction details. Evaluation used the recorded condition-free inputs; original
condition metadata does not imply that it was provided to the forward models.

Reference is a dataset comparator, not a ninth planning method. Its
`record_kind` and `status` are `reference`; `strict_closed` is null because this
entry is not assigned a manuscript method-closure outcome. It is not included in
the eight-method closure totals below and does not inherit a method's terminal
audit or whole-route review.

### Main Method Outputs

Selection follows the frozen Fig. 3 endpoint tables and their method-specific route
inventories. The endpoint totals reproduced by the browser are Rachel 111,
Direct LLM 95, AOT* 95, SyntheLite 83, PaRoutesModel + Retro* 63,
RootAligned + Retro* 73, RootAligned + MCTS 63 and LocalRetro + Retro* 50
for PaRoutes120; RF25 totals are 24, 15, 14 and 6, respectively.
These are the manuscript's planning endpoints, not new evaluations or synthesis results.

Native trees retain their node identities. Nested baseline trees retain molecule
occurrences. Flat reaction lists are joined by stereochemistry-aware canonical
structure identity. An explicitly recorded transformation between distinct SMILES
that canonicalize to the same structure is retained as separate source states;
it is not silently removed. This applies to Direct LLM `n5_4328`, step 4.
Original reaction SMILES remain available for inspection. No chemistry is repaired.

Forward-model records are attached only when the product and multiset of reactants
match the displayed step after stereochemistry-aware canonicalization. Whole-route
reviews require the entire reaction multiset to match. Missing evaluations or audits
are displayed as missing; new model outputs and vendor queries were not generated.

### Two Distinct Rachel Route Versions

For PaRoutes120 `n1_2664` and `n5_2779`, the frozen process trees contain 6 and 11
reactions, whereas the adopted whole-route review inputs contain 5 and 10.
Their chemical contents differ. Both versions are included with a selector.
Review rationales appear only on the evaluated variant. Its graph is not given the
frozen tree's molecule-level terminal audit or automatically labelled strictly closed.
The route record separately retains the manuscript's target-method endpoint.
`ROUTE_VERSION_DIFFERENCES.json` lists the exact reaction differences.

### Source Reconciliation

The current bytes of 112 SyntheLite candidate bundles and 116 AOT result files do
not match their historical recorded hashes. For SyntheLite, all 83 adopted strict
route reaction multisets match the frozen evaluator inputs; 29 open candidates
match their preserved candidate views. AOT's 95 adopted strict routes match the
frozen evaluator inputs; other differing records match their preserved partial
views. The actual and recorded hashes are both retained. Chemical agreement is
not represented as byte-identical recovery. `SOURCE_RECONCILIATION.json` records
the checks. Original local files were not overwritten.

RF25 Rachel `rf25_012` (TMN-CPG) uses the adopted **eight-reaction extracted route**,
with the recorded original tree hash. The original session/tree/terminal files
remain unavailable. The five-step run is not substituted.

RF25 target records also include the frozen cohort's publication DOI and molecular
descriptors. All 25 registry structures were checked against the cohort inputs.

## Files

| File | Contents |
| --- | --- |
| `PaRoutes120.html`, `RF25.html` | Complete offline browsers with embedded molecular images and libraries |
| `data/*.json` | Complete normalized targets, selected routes, endpoints, matched evidence and source references |
| `data/reference/*.json` | Original PaRoutes n1/n5 reference subsets and target manifests, with matching source hashes |
| `ROUTE_INDEX.csv` | One row per target-method position, selection and display status |
| `SOURCE_CATALOG.csv` | Hash identities of upstream files read by the builder; raw upstream files are not all duplicated in this companion |
| `SOURCE_RECONCILIATION.json` | Explicit byte-hash differences and chemical comparison bases |
| `ROUTE_VERSION_DIFFERENCES.json` | Two process-tree versus review-input differences |
| `DATA_DICTIONARY.md` | Normalized schema and interpretation of fields |
| `DATA_VALIDATION.json`, `BROWSER_VALIDATION.json` | Data checks and real-browser interaction checks |
| `CHECKSUMS.csv`, `verify_package.py` | Portable integrity and structural checks |
| `THIRD_PARTY_NOTICES.txt` | Notices for bundled drawing-layout and icon libraries |

Run `python verify_package.py` after extraction to verify the distributed files.
This check uses only the Python standard library. It verifies the companion itself;
it does not rerun planning, forward models, PubChem audits or the whole manuscript's
source-data analysis. The existing Source Data portable package remains the
figure-level evidence release candidate; this atlas is a companion to it.

## Reproduction and Rights

Maintainer sources are in `docs/MANU/manuscript/tools/route_browser/` in the project.
The builder requires the manuscript's original source manifests, adopted evidence,
Python with RDKit/Pillow/BeautifulSoup, and the pinned Dagre/Lucide JavaScript packages.
The JSON datasets preserve the route chemistry for independent rendering without
that local directory layout. PNG depiction keys resolve inside the HTML files.

This directory distributes the manuscript's adopted route results in the
[Rachel repository](https://github.com/ChazenLi/Rachel). Rachel-authored material
is covered by the repository's [CC BY-NC 4.0 license](../../LICENSE).
Third-party libraries and PaRoutes reference material retain their respective
upstream terms and attribution; see `THIRD_PARTY_NOTICES.txt` and the
[PaRoutes source record](https://doi.org/10.5281/zenodo.7341155).
The local frozen package is preserved separately. Only this README and its
checksum entry were updated for the GitHub distribution; route data, offline
viewers and recorded validation results are unchanged.
