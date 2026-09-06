---
name: rachel
description: >
  Use Rachel for information-grounded, LLM-directed, stateful multi-step
  retrosynthesis. Rachel supplies chemical facts, candidate scaffolds,
  validation findings, persistent state, and auditable knowledge provenance;
  the LLM or chemist designs and decides the chemistry. Use for de novo route
  construction, route review, validation, reporting, export, controlled outcome
  recording, reusable-knowledge drafting, expert-governed pack review and
  publication, and debugging or maintenance inside Rachel-beta.
---
# Rachel

## 1. Role And Authority

Rachel is a stateful retrosynthesis workflow:
```text
state -> action -> validation -> commit
```
Rachel supplies:

- molecular, functional-group, site, topology, complexity, and atom-source facts;
- listed reaction and Smart CAP hints, sandbox execution, and canonical validation;
- compatibility, site-fidelity, mapping, risk, and experience observations;
- persistent route/audit state, reports, exports, and pinned knowledge.

The LLM or chemist owns:

- the route thesis and its revision;
- listed versus custom reaction choice and complete precursors/reagents;
- mechanism, selectivity, stereochemistry, and route-quality reasoning;
- reconciliation of validation with evidence;
- commit, override, continuation, terminal, and knowledge-draft decisions.

Rachel does not replace chemistry judgment. A template, family name, score,
Smart CAP hint, experience card, or `clear` gate is evidence, never proof or a
command to choose that reaction.
## 2. Public Surfaces And State Boundaries

Use `LoggedRetroCmd.execute(command, args)` for formal runs. It delegates to
the public `RetroCmd` command layer while adding host-side command provenance.
Operate through these public surfaces:

| Public surface | What to do with it |
| --- | --- |
| State and prompts | Keep one `session.json` per run. Advance/inspect with `next`, `context`, `status`, and `tree`; read every returned `prompt_brief`. |
| Strategy | Record direction with `guide`, `route_plan`, and `route_sketch`; these do not assert chemistry occurred. |
| Action space | Inspect `reaction_sites` and `explore_site`; register a complete unlisted hypothesis with `propose_action`. |
| Sandbox and validation | Run `try_action`, compare `sandbox_list`, and audit canonical `validation` before commit. |
| Route and review | Use `commit` for one reaction, `accept` for a justified terminal, and the finalize/report/export/review commands for completion or re-analysis. |
| Learning and knowledge | Record only supplied outcomes; create inactive drafts after review; select published packs only at `init`. |

State boundaries are strict:

- only a successful `commit` makes a reaction attempt a route fact;
- before finalization, only `accept` or automatic terminal handling makes a molecule terminal; `queue_empty` is necessary but not sufficient, so let `finalize` verify that no completion blocker remains and never use it to close unresolved leaves;
- `guide`, `route_plan`, and `route_sketch` change strategy state, not chemistry facts;
- `record_outcome` adds an experimental fact but never rewrites `committed`;
- `propose_knowledge` adds only an inactive draft;
- only a later `init` selecting a published pack can activate new reusable knowledge.

Use only the 34 public commands and canonical `validation` documented below;
do not depend on legacy aliases or internal session fields.

## 3. Non-Negotiable Chemistry Contract

1. Chemistry first. Chemical truth and route quality outrank templates, scores, confidence labels, convenience, action ordering, and route depth.
2. Design boldly, prove strictly. Unconventional and de novo chemistry is welcome when it is the better route hypothesis, but every commit must survive mechanism, atom-source, site, topology, compatibility, selectivity, and stereochemistry review.
3. Commit exactly one real chemical event. A one-pot operation may be one event only when its chemistry is genuinely concerted or operationally inseparable; do not compress an unrecorded multi-step sequence into one tree edge.
4. Treat listed and complete custom actions as peer hypotheses. Rachel's action space does not bound retrosynthesis.
5. Distinguish evidence classes. Route planning, validator observations, experimental outcomes, literature, and internal SOP evidence are not interchangeable.
6. True chemical contradictions block. A validator or system failure stops execution until repaired, but is not chemical disproof.
7. Continue toward simpler credible precursors. Prefer an honest advanced terminal over a fake deep disconnection or atom-source fiction.

Before every commit, explicitly ask:

- Is this a real mechanism under plausible conditions?
- Where does every new key C, N, S, halogen, carbonyl, protecting-group, or ring atom come from?
- Is the changed bond at the intended site, with the required ring size, fusion, bridge, spiro relation, and retained substituent positions?
- Are free acids, amines, alcohols, phenols, thiols, aldehydes, halides, organometallics, and temporary reactive handles compatible with the step?
- What controls chemo-, regio-, site-, and stereoselectivity?
- Does the step support the route thesis and improve route quality?
## 4. Document Boundary

Normally use only this Skill, current command output, `prompt_brief`, and the
smallest decision payload. Do not load full source, session JSON, historical
walkthroughs, `experience.md`, or long project documents during route work.

Read extra material only when required:

- `workflow.md`: state-machine and design rationale;
- `refs.md`: exact command arguments, payload fields, and exports;
- `knowledge/AUTHORING_GUIDE.md`: required before actually adding, replacing, disabling, reviewing, or publishing reusable knowledge;
- `knowledge/README.md`: profile architecture and pack CLI overview;
- source and tests: only for maintenance, debugging, or missing protocol facts.

If a field is unclear, inspect the smallest relevant reference fragment; never
guess a public argument.
## 5. Setup And Independent Run Directory

Run from the project/release root with Python 3.10+, RDKit, NumPy, and Pillow.
Use the interpreter supplied by the user/deployment; never guess an environment
name or hard-code a developer path. Give every formal run its own session:
```text
<workspace>\walkthrough_runs\YYYYMMDD_HHMMSS_<target_slug>\session.json
```
If `session.json` already exists and the user did not request resume, create a new timestamped directory. Do not mix two targets or two knowledge profiles in one run directory.
Preferred Python invocation:
```python
from pathlib import Path
from Rachel.tools.logged_runner import LoggedRetroCmd
workspace = Path.cwd()
run = workspace / "walkthrough_runs" / "YYYYMMDD_HHMMSS_target_slug"
run.mkdir(parents=True, exist_ok=True)
cmd = LoggedRetroCmd(str(run / "session.json"))
cmd.execute("init", {
    "target": "TARGET_SMILES",
    "name": "target_slug",
    "max_depth": 15,
    "max_steps": 50,
    "terminal_cs_threshold": 2.25,
})
```
Use `LoggedRetroCmd` for every formal new, resumed, terminal-review, and
node-review run. It writes append-only `NNN_<command>.input.json` and
`.output.json` files beside `session.json`; the session remains authoritative.
Give every `review_node` variant its own directory, for example
`runs/route_variant/session.json`. On success its review command and all later
logs belong to that directory; the finalized parent receives no copied logs.
Avoid inline JSON quoting in PowerShell. For multi-line Python, use a PowerShell here-string pipe, not Bash heredoc syntax:
```powershell
$env:PYTHONIOENCODING='utf-8'
$OutputEncoding=[System.Text.UTF8Encoding]::new($false)
[Console]::OutputEncoding=[System.Text.UTF8Encoding]::new($false)
@'
from Rachel.tools.logged_runner import LoggedRetroCmd
print("ok")
'@ | python -
```
Do not use `python - <<'PY'` in PowerShell.
## 6. Complete Public Command Map

Rachel exposes 34 public commands. Use these names, not hidden legacy aliases.

| Category | Commands |
| --- | --- |
| Session and context | `init`, `next`, `context`, `inspect_structures`, `status`, `tree` |
| Direction and strategy | `guide`, `route_plan`, `route_sketch` |
| Controlled learning | `record_outcome`, `learning_review`, `propose_knowledge` |
| Site and action space | `reaction_sites`, `explore_site`, `propose_action` |
| Optional node scouting | `scout_node`, `scout_record` |
| Sandbox | `try_action`, `sandbox_list`, `sandbox_clear`, `select` |
| Route mutation | `commit`, `accept`, `review_terminal`, `review_node`, `skip` |
| Continuation | `continuation_status`, `continuation_abort` |
| Completion and output | `finalize`, `audit`, `report`, `export` |
| Expert structure aids | `smart_cap`, `custom_cap` |

Every public command has this operational contract:

| Command | Use and required/important arguments | State effect and normal next step |
| --- | --- | --- |
| `init` | Start a new session. Required: `target`; optional `name`, budgets, `knowledge_profile`, `knowledge_roots`. | Creates and pins the session/profile; then `next`. |
| `next` | Activate the next pending molecule; optional `additional_steps` only after `step budget exhausted`. | May auto-close quick terminals. Budget exhaustion retains work; retry with a positive extension. On `queue_empty`, call `finalize` to verify all completion blockers. |
| `context` | Re-read current context. `detail=compact|structure|full|diagnostic|status|tree`; offsets are diagnostic-only. | Read-only; continue current state. |
| `inspect_structures` | Compare labeled SMILES with optional current molecule. `molecules`, `include_current`. | Read-only transient structure facts; then refine action or audit. |
| `guide` | Record expert direction after `next`. Required `text`; optional `intent`, site/reaction/precursor/constraint/terminal hints and `summary`. | Persists guidance, not a route fact; then normal site/action loop. |
| `record_outcome` | Record a real result for committed `step_id`. Required `outcome=success|partial|failure`; optional matching `action_id`, conditions, yield/conversion, observations, evidence. | Appends experimental fact without changing commit; normally after evidence exists. |
| `learning_review` | Review a route only after explicit complete `finalize`; no arguments. | Read-only deterministic synthesis of route and outcome evidence; then optional knowledge draft. |
| `propose_knowledge` | Create inactive draft. Required `target_pack_id`, `resource`, `entry`, `rationale`, non-empty `source_refs`; optional `evidence`. | Persists `active=false` draft; never changes runtime profile. |
| `route_plan` | Set/replace global thesis. Required `route_thesis`; use mode, evidence, risks, triggers, disconnections, precursor logic, preserve policy, terminal policy, revision reason. | Persists complete plan revision; then gather/execute local evidence. |
| `route_sketch` | Translate route-level or multi-event idea. Required `problem`, `macro_strategy`; optional disconnections, comparison reason, next executable step, terminal flag, summary, `continuation_steps`. | Persists strategy only; then exactly one listed/custom executable action. |
| `reaction_sites` | Build complete current site/reaction map; no arguments, active context required. | Read-only disclosure; choose a current `site_id`. |
| `scout_node` | Build an optional read-only scouting batch in an active standard context. Required `tasks`; optional `expansion_reason` and `frontier_reason`. | Returns immutable blind/informed task packets without saving state; dispatch them externally, then use `scout_record`. |
| `scout_record` | Validate and persist one complete scouting round. Required `round_binding`, `results`, `review_summary`; optional shortlist IDs and deferred seeds. | Records bound, idempotent scouting evidence without changing route chemistry; then return to the normal site/action loop. |
| `explore_site` | Expand same-site competing actions. Required current `site_id`. | Read-only disclosure; choose listed action or design a peer custom action. |
| `try_action` | Sandbox one listed or registered custom action. Required current `action_id`; optional `scouting_source={round_id,task_id,adoption_reason?}`. | Adds an attempt and canonical validation; source is saved only if an attempt is actually created; then compare with `sandbox_list`. |
| `propose_action` | Register a complete custom one-step hypothesis. Required complete `precursors`; provide reaction identity and positive chemical rationale; add reagents and high-risk evidence as needed. | Adds custom action only; then `try_action(returned_action_id)`. |
| `sandbox_list` | Compare all current attempts; no arguments. | Read-only view; retry, clear, select, or commit. |
| `sandbox_clear` | Remove current uncommitted attempts; no arguments. | Clears sandbox only; rerun intended action. |
| `select` | Mark one attempt. Required `idx`. | Selects sandbox row but does not mutate route; then audit and `commit`. |
| `commit` | Write one sandbox event. Required `idx`, matching `expected_action_id`, reasoning; include reaction-specific `risk_assessment` and structured `process_conditions`; optional confidence, real rejections, route-plan alignment/note, and justified override. | Adds reaction and precursor nodes; risk documentation does not rank or select chemistry; then `next`. |
| `accept` | Mark active molecule terminal. Required chemical `reason`; advanced no-actionable rescue also needs the explicit force/reason fields. | Adds terminal fact; then `next`. |
| `review_terminal` | Reopen a terminal leaf. Required `smiles`, `reason`; optional positive `additional_steps`. | Revokes route completion and requeues same node; then `next`. |
| `review_node` | Create independent variant from finalized route. Required new `variant_session_file`, `reason`, and `node_id` or `smiles`; optional directed guidance, `review_seed_id`, `adaptation_reason`, and budget. | Parent stays unchanged; a chosen seed becomes guidance only, affected complete rounds are removed from the variant, and the command rebinds to variant; then `next`. |
| `skip` | Skip current molecule only with explicit `reason`. | Records an unresolved planning stop; only a later real `commit` or reasoned `accept` for that molecule resolves it. |
| `tree` | Inspect route tree and terminals; no arguments. | Read-only. |
| `status` | Inspect session, budgets, pending nodes/continuations, and profile; no arguments. | Read-only. |
| `continuation_status` | Inspect active multi-event strategy continuations; no arguments. | Read-only; realize next event or abort. |
| `continuation_abort` | Close non-actionable continuation. Required `continuation_id`, chemical `reason`. | Closes only the continuation; it never resolves an existing skip. Then `next` or `finalize` as state permits. |
| `finalize` | Call after `next` reports `queue_empty`; optional `summary`. | Succeeds only when active context, queue, incomplete leaves, unresolved skips, and continuations are all closed; then enables learning review/report/export. |
| `audit` | Independently audit an explicitly finalized route. Optional `scopes=[buyability,safety]`, `offline`, `no_vendors`, `query_reagents`, `timeout`, and `pause`. | Stores a route-digest-bound evidence snapshot; never changes chemistry, validation, ranking, or terminal choices; then `report`/`export`. |
| `report` | Generate forward-style text and starting-material list; no arguments. | Read-only output after route construction. |
| `export` | Write artifacts. Important `output_dir`, optional `name`. | Writes report/tree/session/profile/visual assets and current audit sidecars; does not change chemistry. |
| `smart_cap` | Ask for capping suggestions using current `bond_idx` or explicit `smiles` + atom-pair `bond`; optional `max`. | Read-only heuristic; complete useful result through `propose_action`. |
| `custom_cap` | Apply expert caps. Required `cap_i`, `cap_j`, plus current `bond_idx` or explicit `smiles` + `bond`; optional `reaction_type`. | Read-only structural proposal; complete through `propose_action`. |

Do not substitute legacy or diagnostic aliases for these public commands.
## 7. Complete Workflow Map And Default Route Loop
```text
init
-> next / context(compact)
-> for a complex target, optionally record a provisional route_plan seed
   or inspect reaction_sites first when site evidence is needed
-> optional guide
-> reaction_sites
-> optional current-node scouting:
   scout_node -> host runs independent scouts -> scout_record
-> explore_site(site_id)
-> write or enrich route_plan from actual evidence
-> choose exactly one next event:
   listed action -> try_action(action_id)
   complete custom one-step action -> propose_action -> try_action(custom_id)
   route-level or multi-event idea -> route_sketch -> one executable action
-> sandbox_list
-> commit or accept
-> next
-> resolve any strategy continuation
-> finalize
-> audit
-> report
-> export
```
Commands are stateful. Do not run `next`, `reaction_sites`, `explore_site`, and `try_action` in parallel.
State rules:

- `next` selects the next active molecule from the queue.
- `reaction_sites` requires an active context.
- `explore_site` must use a current `site_id`.
- `try_action` must use a listed or registered custom `action_id`.
- `propose_action` registers a hypothesis; it does not validate it.
- `commit` accepts only a sandbox attempt.
- after `commit`, `accept`, or `skip`, call `next` before another local action.
- stop immediately on any `{"error": ...}` response.

Minimal progression:
```python
node = cmd.execute("next", {})
sites = cmd.execute("reaction_sites", {})
detail = cmd.execute("explore_site", {"site_id": "site:bond:2"})
attempt = cmd.execute("try_action", {"action_id": "site:bond:2:alt:0"})
comparison = cmd.execute("sandbox_list", {})
```

### Optional Current-Node Scouting

Scouting is opt-in for a directionally difficult active node. Difficulty here
means unresolved uncertainty about what bond disconnection or FGI to pursue at
the current structure; molecular size, atom count, route depth, or a complexity
score alone is not a trigger. Rachel neither auto-triggers Scouting nor creates
a task list. After reading the current structure and site/action evidence, the
main LLM calls `scout_node` only when the expected local information gain
justifies it: the listed directions are materially narrow or repetitive,
several distinct disconnections or FGI strategies remain genuinely competitive,
prior blocked/inconclusive attempts have not resolved direction, or the user
requests an independent challenge. Continue directly through `explore_site`
when one well-supported local direction is sufficient; routine nodes do not
need a scouting round.

Treat results as advisory current-node hypotheses. Require each scout to name
one concrete reaction idea and its indexed target site; for example,
"disconnect `bond_index=k` at `atoms=[i,j]` with named forward reaction logic"
or "apply a named FGI at `atom_indices=[...]`". The main LLM collapses duplicate
directions and independently adopts, revises, or ignores each through Rachel's
normal action and sandbox flow. Scouting is not a route planner, validator, or
vote.

When Scouting is chosen, the main LLM supplies the whole batch before launch;
`tasks` has no code-generated default. The routine batch is exactly 5 tasks
with at least 2 `blind` and 3 `informed`; 6-10 tasks require a non-empty
round-level `expansion_reason`. Every task needs a distinct one-line `focus`.

Optionally give a task one supported `methodology_lens` when a concrete
methodological angle is useful; the six stable lens IDs and their exact packet
prompts are in [refs.md](refs.md). A lens is a lightweight focus, not a persona,
reaction whitelist, or automatic trigger.

Use `frontier_hypothesis=true` only for a deliberately unprecedented
mechanistic proposal and always pair it with a methodology lens. Normally use
at most one such task in a round. A high-value node may use exactly two only
with a non-empty round-level `frontier_reason`; more than two are rejected.
Its direction must say `hypothesis_status=unprecedented_hypothesis` and give a
compact `hypothesis_basis`. Keep the full mechanism, proposed operation,
workup/purification, and falsification analysis on the host; Rachel stores only
the compact structured direction and hypothesis label.
A deferred review seed sourced from a frontier direction automatically retains
that label and basis when copied into a later `review_node` variant.
Related-node prompts retain a frontier seed's compact label and basis. Human
reports prominently disclose only `unprecedented_hypothesis`, not other statuses.

```python
prepared = cmd.execute("scout_node", {
    "tasks": [
        {"visibility": "blind", "focus": "independent topology challenge"},
        {"visibility": "blind", "focus": "alternative bond logic"},
        {"visibility": "informed", "focus": "challenge the listed action space"},
        {"visibility": "informed", "focus": "late FGI opportunity"},
        {"visibility": "informed", "focus": "selectivity-aware local direction"},
    ],
})
```

`scout_node` is read-only. Dispatch its packets independently outside Rachel and
retain full free-text analyses on the host. Submit exactly one complete batch to
`scout_record`; Rachel stores only each task's structured direction,
abstain/unavailable state, the binding, the main LLM review, 0-3 shortlist IDs,
and 0-3 deferred review seeds. Shortlist is advisory and no consensus score is
computed. Keep `review_summary` and shortlisted structured cores concise and
site-specific: the current schema normalizes these strings but does not impose
a hard character limit. Full reasoning belongs in the host-side scout analyses.
An identical recorded-round replay remains idempotent after context advancement;
a first submission still requires the active-node binding, and conflicts reject.

A recorded direction still follows the existing chemistry path:

```text
listed action -> try_action(action_id, scouting_source)
unlisted one-step direction -> propose_action -> try_action(custom_id, scouting_source)
multi-event direction -> route_sketch -> materialize one real event -> try_action(..., scouting_source)
```

`scouting_source` identifies one directly originating task. A non-shortlist task
requires `adoption_reason`; matching is checked against available structured
action site fields. The source is written only when `try_action` creates a real
sandbox attempt and then propagates through commit/tree/report/export without
changing validation policy.
System FGI actions expose stable `fgi_idx` identity but no atom map, so Rachel
checks functional-group kind plus that identity; the main LLM must verify the
recorded target site. Custom actions require matching `changed_bonds.product_atoms`.

All rounds remain in session/export. Ordinary `prompt_brief` exposes only the
active node's latest `scouting_review`, its shortlist cores,
`scouting_round_count`, and one short `scouting_guidance`. That guidance remains
present throughout the related node's prompt lifecycle and is absent on nodes
without a recorded round. A prepared but unrecorded `scout_node` batch is not
mounted into ordinary prompts. Treat `scouting_source` as provenance only;
validation comes from the realized action. Reports also show
reviewed-but-not-adopted rounds and deferred seed cores without writing those
uncommitted facts into the route tree.

Use `review_node(review_seed_id=...)` only to copy one deferred seed into a
finalized-route variant as self-contained guidance. It is not an action and must
return through `next` and the normal sandbox flow. Changes to seed-derived
`instruction`, `site_hint`, `reaction_hint`, `precursors`, or `constraints`
require `adaptation_reason`; Rachel records `adapted_fields`. Other semantic
changes remain the main LLM's audit responsibility. The variant removes complete
rounds for the reviewed node and pruned subtree, while the parent session stays
unchanged.

Use this branch map; later sections provide complete payloads and gates:

| Situation | Required path or constraint |
| --- | --- |
| New run | New directory -> `LoggedRetroCmd(new session)` -> `init` -> `next`. |
| Resume | `LoggedRetroCmd(existing session)` -> `status` -> `context(compact)` or `next`; never call `init`. Missing/changed pinned packs are hard failures, not permission to fall back to base. |
| Strategy grounding | With useful compact evidence, write provisional `route_plan` revision 0, inspect sites, then replace it with evidence-enriched revision 1. Otherwise inspect sites first and write complete revision 0. Revise when later evidence changes the thesis. |
| Action type | Listed -> `try_action`; complete custom -> `propose_action -> try_action`; wildcard/Smart CAP -> complete real structures, then the custom path; multi-event idea -> `route_sketch`, then exactly one executable event. |
| Validation | `try_action -> sandbox_list`; repair/change `blocked`, prove or revise `proof_required`, separate gaps from tool limits for `inconclusive`, address `warning`, and still audit `clear`; then optional `select` -> `commit(idx + expected_action_id)` -> `next`. |
| Terminal/continuation | Simple credible terminal -> `accept(reason)`; advanced terminal -> bounded mini-route review; realize a credible event normally or give the explicit no-actionable rationale. `continuation_abort` closes only the continuation. Finalize only with no active context, queue, incomplete leaf, unresolved skip, or continuation. |
| Reopen/review | `review_terminal` requeues the same leaf and changes the original tree; `review_node` preserves the finalized parent, rebinds to a new variant, and runs the normal loop there. Finalize the resulting route again. |
| Experimental learning | Supplied experiment -> `record_outcome`; after explicit `finalize`, run `learning_review -> propose_knowledge(active=false) -> import-drafts -> complete metadata -> expert approve/reject -> conflicts/validate -> publish new version`; only a later `init` can select it. |

All action paths converge on the same sandbox, validation, comparison, and
commit gate. Never infer an experiment, and never let a draft or approval alter
the active session. Never let local template convenience silently replace the
route thesis.

## 8. Context And Payload Reading

Read `next()["current"]` for the active SMILES, compact cognition, and `prompt_brief`.
Use context layers deliberately:

- `context(detail="compact")`: molecular summary and opportunity brief;
- `reaction_sites()`: complete first-layer site/reaction menu;
- `explore_site(site_id)`: selected-site actions and local bond/FG facts;
- `context(detail="structure")`: indexed atoms, bonds, rings, stereo, symmetry, and scaffold topology for the active molecule;
- `inspect_structures(...)`: read-only comparison of current product and up to five labeled candidate structures.

```python
compact = cmd.execute("context", {"detail": "compact"})
structure = cmd.execute("context", {"detail": "structure"})
```

Use `detail="status"` for the smallest state-only view and `detail="tree"` when
the printed tree is specifically needed. Use `detail="full"` or
`detail="diagnostic"` only for debugging, maintenance, or a missing protocol
fact; offsets and limits are diagnostic controls, not the normal LLM workflow.

Example structural inspection:
```python
cmd.execute("inspect_structures", {
    "include_current": True,
    "molecules": [
        {"label": "precursor_A", "smiles": "..."},
        {"label": "precursor_B", "smiles": "..."},
    ],
})
```
Labels must be unique. Invalid SMILES fail explicitly. The command is read-only and does not change the session or route tree.
Important public fields:

- `reaction_sites`: read `site_reaction_map`, not guessed `sites`/`items`; rows include `site_id`, `site_type`, `site_hint`, `action_count`, `reaction_count`, `competition_hint`, `risk_hint`, and `reactions`.
- `explore_site`: read `actions`; rows include `action_id`, `site_id`, `reaction.id`, `reaction.name`, `source`, `precursors_preview`, `risk_tags`, and relevant bond/FG identifiers.
- `sandbox_list`: read `attempts`; `by_site`/`by_reaction` contain indexes. Each attempt exposes canonical `rachel.validation.v2` `validation`.
- absent `mechanism_interpretation` means no registered interpretation added useful evidence; an unregistered reaction name is not itself a warning.

Read every `prompt_brief`: it may contain `route_plan_brief`,
`route_strategy_brief`, `strategy_continuation_brief`, `chemist_guidance`,
`experience_prompts`, `required_audit_fields`, and `next_actions`.

Experience cards are reminders, not verdicts and not a quota. Discovery may mount all relevant cards or none. Audit mounting remains deliberately capped. Do not load the full experience files during normal execution.
## 9. Route Strategy Commands
### `route_plan`

Use `route_plan` for the persistent global synthesis thesis. It does not mutate the route tree.
Follow the provisional-first or evidence-first path in the workflow map. A seed
is a hypothesis, not proof. Every revision replaces the prior plan, so restate
all applicable fields instead of sending a partial patch.

Use these command-enforced writing budgets: `route_thesis<=360` characters;
`route_mode<=100`; `key_disconnections<=5` items of `<=160` each;
`preferred_precursor_logic<=4x180`; `protect_or_preserve<=4x160`;
`mode_evidence`, `strategic_risks`, and `revision_triggers<=5x180` each;
`terminal_rescue_policy<=220`; and `revision_reason<=160`. If any field is
over limit or has the wrong type, Rachel rejects the whole update and identifies
the field; shorten it and retry. Rachel never truncates an accepted plan.
`prompt_brief.route_plan_brief` contains the complete bounded current plan;
revision history remains audit-only.

```python
cmd.execute("route_plan", {
    "route_thesis": "short current synthesis thesis",
    "route_mode": "late_fgi",
    "mode_evidence": ["why this route mode fits"],
    "strategic_risks": ["what may falsify it"],
    "revision_triggers": ["evidence that requires revision"],
    "key_disconnections": ["key bond or FG logic"],
    "preferred_precursor_logic": ["precursor family or handle"],
    "protect_or_preserve": ["motifs and scaffolds to preserve"],
    "terminal_rescue_policy": "Before accepting an advanced terminal, test a bounded mechanistic rollback.",
    "revision_reason": "initial",
})
```
Compare route modes rather than jumping to a local template: `late_fgi`
preserves a mature core for late handle edits; `scaffold_assembly` builds the
core/rings from simpler pieces; `electronic_state_strategy` builds through a
non-final electronic state; `hybrid` combines modes explicitly.

Revise when evidence changes a substantive claim: key disconnection, precursor family, event sequence, handle timing, preserve/build choice, selectivity strategy, strategic risk, or revision trigger. A reagent or solvent refinement inside the same planned event normally does not require a plan revision.
### `guide`

Use `guide` to record high-priority chemist direction for the current active node. Guidance affects prompts and audit but is not validation evidence and does not directly mutate the route tree.
```python
cmd.execute("guide", {
    "text": "Preserve the fused core and test late N-functionalization.",
    "intent": "reaction_hint",
    "reaction_hint": "late N-functionalization",
    "constraints": ["do not break the fused core"],
    "summary": "Preserve core; test late N-functionalization.",
})
```
Convert guidance into site, reaction, precursor, or terminal constraints. Do not skip candidate comparison, sandboxing, or commit audit because an expert specified a direction.
### `route_sketch`

Use `route_sketch` for strategy-to-action translation, a multi-event idea, a
route-thesis change, or bounded terminal review. It is not required for every
custom action. The strict multi-event form is:

```python
sketch = cmd.execute("route_sketch", {
    "problem": "A terminal needs a two-event mechanistic rollback.",
    "macro_strategy": "Execute one FG rollback, then expose the simpler precursor.",
    "key_disconnections": ["first rollback", "follow-up disconnection"],
    "rejected_action_space_reason": "Only when actual listed actions were rejected.",
    "next_executable_step": "propose_action",
    "terminal_review": True,
    "summary": "Execute and validate one rollback event at a time.",
    "continuation_steps": [
        {
            "step_idx": 0,
            "target_smiles": "CURRENT_SMILES",
            "reaction_name": "first real event",
            "expected_precursors": ["PRECURSOR_A", "CO_REAGENT"],
            "continuation_precursor": "PRECURSOR_A",
        },
        {
            "step_idx": 1,
            "target_smiles": "PRECURSOR_A",
            "reaction_name": "follow-up real event",
            "expected_precursors": ["SIMPLER_PRECURSOR"],
        },
    ],
})
```

For one next event, use the same payload without `continuation_steps` and set
`terminal_review=False`. Convert and commit only one real event at a time.

Allowed `next_executable_step` values are `propose_action`, `explore_site`,
`try_action`, and `reaction_sites`. Never put `accept` there.

Bind the first custom event to the returned sketch:

```python
proposal = cmd.execute("propose_action", {
    "route_sketch_id": sketch["route_sketch_id"],
    "precursors": ["PRECURSOR_1", "PRECURSOR_2"],
    "reaction_name": "first real event",
    "rationale_summary": "mechanism, atom source, site, selectivity, and route fit",
})
```

After that event is committed, read `prompt_brief.strategy_continuation_brief`.
When it requests the next event, pass its public identifiers back with the next
custom action:

```python
cmd.execute("propose_action", {
    "continuation_id": continuation_id,
    "continuation_step_idx": continuation_step_idx,
    "continuation_precursor": "PRECURSOR_BEING_EXPANDED",
    "precursors": ["NEXT_PRECURSOR"],
    "reaction_name": "next real event",
    "rationale_summary": "why this is exactly the requested continuation event",
})
```

Copy these identifiers from current Rachel output; never invent them. Every
continuation event still requires `try_action`, validation, and its own `commit`.
## 10. Site, Listed Action, And Custom Action Workflow

First inspect the real site menu, then the actions at one selected site. Do not infer the full action space from compact context.

Listed and complete custom actions are peers. A custom action may be a stronger
same-site reaction, an unlisted disconnection, a completed wildcard, or the
first event from a route sketch.

Lead with the positive chemical case. Do not invent a rejection story merely to gain access to `propose_action`.
```python
proposal = cmd.execute("propose_action", {
    "precursors": ["PRECURSOR_1", "PRECURSOR_2"],
    "reaction_name": "specific one-step reaction",
    "action_label": "short label",
    "rationale_summary": "mechanism, site, atom source, selectivity, route fit",
    "reagents": ["reagent or catalyst"],
    "risk_tags": ["chemoselectivity"],
})
custom_id = proposal["action_id"]
attempt = cmd.execute("try_action", {"action_id": custom_id})
```
Record rejected alternatives only when they are actually rejected after real
comparison; keep their action IDs and reasons auditable.

Put route-bearing tree nodes in `precursors`; put current-event catalysts,
metals, donors, and non-leaf components in `reagents`. Both participate in atom
accounting, while scaffold/topology/site audits treat `precursors` as
route-bearing. Repeat identical SMILES when two copies supply two product
sites; prose such as "two equivalents" supplies no atoms to validation.

Additional selection rules:

- preserve a mature rigid/heteroaryl scaffold when routes are comparably
  credible; permit scaffold assembly when mechanism and precursors make it better;
- install highly reactive temporary handles late when practical;
- record protection and deprotection as real events, never as hidden reasoning;
- when chemistry is comparable, prefer the action that moves route-bearing
  nodes toward simpler, stabler, and more purchasable precursors.
### Wildcard previews

A `*` in `precursors_preview` is an attachment cue, not a reagent or proof of a
broken template. Complete that same disconnection with closed precursors and
required leaving groups/reagents through `propose_action -> try_action`; do not
count preview and completion as different chemistry alternatives.
### Smart CAP and custom CAP

`smart_cap` and `custom_cap` are expert structure aids. They propose structural capping hypotheses and cannot be committed directly.
```python
cmd.execute("smart_cap", {"bond_idx": 0})
cmd.execute("custom_cap", {
    "bond_idx": 0,
    "cap_i": "Br",
    "cap_j": "B(O)O",
    "reaction_type": "llm_custom",
})
```
Complete a useful cap as a real one-step `propose_action`, then `try_action`.
Trying a raw hint may return `llm_completion_required` without an attempt.
### Sandbox discipline

Use `sandbox_list` to compare attempts. If the sandbox is polluted by obsolete or ambiguous trials, call:
```python
cmd.execute("sandbox_clear", {})
```
Then rerun the intended action. `select(idx)` marks a comparison choice;
`commit` remains the only reaction-attempt route mutation.

```python
selected = cmd.execute("select", {"idx": 1})
```

Read the selected row's canonical `validation`, then use the same `idx` and
displayed `action_id` as `expected_action_id` in `commit`.
## 11. Validation V2 And Commit

Sandbox `success` means execution, not permission to commit. Read `validation.decision_gate.state` first:

- `blocked`: do not commit; repair a `system_error` and rerun, or revise the chemistry/precursor.
- `proof_required`: add the named atom-source, tether, anchor, mechanism, site, or selectivity proof, or revise the action.
- `inconclusive`: separate `evidence_gaps` from `tool_limits`; template absence is neither disproof nor support.
- `warning`: proceed only after addressing the risk in reasoning.
- `clear`: no public objection, but still perform the full chemistry review.

Inspect the named evidence, not only the state: `contradictions`, `proof_obligations`, `evidence_gaps`, `tool_limits`, `warnings`, and `system_errors`; reconcile them with `execution`, `observations`, `declared_action`, and `mechanism_interpretation`.

Before commit, repeat this concrete chemistry audit:

1. Mechanism: name the productive event and why substrate and conditions support it.
2. Atom source: account for every new key atom, carbonyl, protecting group, partner, and ring bond.
3. Site/topology: identify changed bonds, anchors, ring size, fusion/bridge/spiro relation, and retained substituent positions.
4. Compatibility: address exposed groups, handle timing, redox state, organometallic stability, and protection/deprotection.
5. Selectivity/stereochemistry: explain chemo-, regio-, facial-, and stereocontrol; never rewrite unspecified stereo silently.
6. Route quality: state `route_plan` alignment/revision, convergence, and actual rejections.
7. Risk: bind hazards to named materials, byproducts, or mechanistic events; separate sourced facts, chemical inference, and unresolved evidence.

Write `risk_assessment` as a concise, reaction-specific planning judgment. Cover applicable material hazards; energetic, gas, pressure, unstable-intermediate, or incompatibility hazards implied by the reaction; credible exposure routes; required controls; and missing SDS/SOP evidence. State uncertainty explicitly. A generic safe/unsafe label, an unbound hazard list, or invented operating detail does not satisfy this audit.

Construct a coherent provisional operating plan for the selected attempt and record it in `process_conditions`. The seven core fields are `solvent`, `equivalents`, `addition_order`, `temperature`, `atmosphere`, `reaction_time`, and `workup`. Add `scale`, `concentration`, or `pressure` only when supplied or materially relevant. Each declared value is one of:

```python
{"status": "specified", "value": "value or range", "basis_kind": "user_provided | literature | experimental | sop | planning_hypothesis", "basis": "source or explicit chemical planning basis"}
{"status": "unknown", "reason": "why no responsible planning proposal can be made"}
```

Every new `specified` value should carry one `basis_kind`. `user_provided`, `literature`, `experimental`, and `sop` identify evidence-bearing sources. When those sources are absent, turn useful reaction- and substrate-specific judgment into an explicit `planning_hypothesis`; absence of literature or SOP alone does not make a field unknown. Use `unknown` only when neither evidence nor a responsible planning suggestion exists. Use `addition_order` for the candidate charge, cooling, addition, and hold sequence; let `workup` contain candidate quench, isolation, and purification. Recheck every planning hypothesis against the exact substrate, reaction-specific risk, literature, SDS, and SOP before execution, and never present it as an established operating condition.

Legacy sessions without `basis_kind` remain valid but reports identify their specified values as source-unclassified. The deterministic audit checks structure and coverage, not chemical correctness. Seven valid non-placeholder core values yield `conditions_available_for_preliminary_process_review`; missing optional fields do not block that planning-review status, and the status never means the reaction or process is safe.

Protection and deprotection are real events. Resolve `functional_group_condition_conflict` with explicit protection or compatible conditions; "assume protected" in reasoning is not a solution.

Commit with both sandbox index and action identity:
```python
cmd.execute("commit", {
    "idx": idx,
    "expected_action_id": "action_id shown at this sandbox idx",
    "reasoning": "mechanism, atom source, site/topology, compatibility, selectivity, stereo, route fit",
    "confidence": "high",
    "rejected": [
        {"action_id": "...", "reason": "only a real compared rejection"}
    ],
    "risk_assessment": "named material/mechanistic hazards, exposure, controls, evidence basis, and unresolved SDS/SOP questions",
    "process_conditions": {
        "solvent": {"status": "specified", "value": "candidate solvent or solvent pair after a solubility screen", "basis_kind": "planning_hypothesis", "basis": "reaction class and current substrate polarity"},
        "equivalents": {"status": "specified", "value": {"substrate": 1.0, "reagent": "1.1-1.5"}, "basis_kind": "planning_hypothesis", "basis": "bounded screening range for the selected event"},
        "addition_order": {"status": "specified", "value": ["charge substrate and solvent", "cool if exothermic", "add reagent gradually", "hold at reaction temperature"], "basis_kind": "planning_hypothesis", "basis": "limit local concentration and thermal excursion"},
        "temperature": {"status": "specified", "value": "candidate range for the selected reaction", "basis_kind": "planning_hypothesis", "basis": "mechanism and functional-group compatibility"},
        "atmosphere": {"status": "specified", "value": "air or inert atmosphere selected from reagent sensitivity", "basis_kind": "planning_hypothesis", "basis": "oxidation and moisture sensitivity review"},
        "reaction_time": {"status": "specified", "value": "screening window with TLC or LC-MS endpoint", "basis_kind": "planning_hypothesis", "basis": "reaction-class kinetics remain unverified"},
        "workup": {"status": "specified", "value": {"quench": "candidate controlled quench", "isolation": "candidate phase handling or filtration", "purification": "recrystallization or chromatography screen"}, "basis_kind": "planning_hypothesis", "basis": "expected reagent residues, byproducts, and product properties"},
    },
    "route_plan_alignment": "supports",
})
```
Always send `expected_action_id`; it prevents a stale/reordered index from committing an action different from the reasoning. Set `route_plan_alignment` to `supports`, `revises`, `conflicts`, or `not_applicable`; add `route_plan_note` only when the label is insufficient.

Use `validation_override` only for a documented `proof_required` false positive or allowed independent-evidence case:
```python
cmd.execute("commit", {
    "idx": idx,
    "expected_action_id": "action_id shown at this sandbox idx",
    "reasoning": "complete chemistry audit plus why independent evidence resolves the finding",
    "validation_override": {
        "allowed": True,
        "reason": "why the named proof obligation is resolved",
        "evidence": "independent structure, mechanism, literature, experiment, or SOP evidence",
    },
})
```
Never override unavailable validation, true atom imbalance, wrong-site chemistry, impossible topology, a missing protecting-group donor, or unsupported selectivity.
## 12. High-Risk Topology, Atom Source, And EAS

A custom action that builds, opens, fuses, bridges, spiro-connects, or rearranges a ring/scaffold needs structured proof; a family label is insufficient. Before sandboxing, provide when applicable:

- `intended_deltas`;
- `expected_ring_change`;
- exact `changed_bonds`;
- at least two independent `preserved_anchors`;
- `mechanistic_evidence`;
- atom-source and tether evidence in `family_evidence`;
- focused `risk_tags`.

```python
cmd.execute("propose_action", {
    "precursors": ["..."],
    "reaction_name": "one custom ring event",
    "action_label": "ring closure",
    "rationale_summary": "single event with explicit atom source",
    "intended_deltas": ["ring_closure"],
    "expected_ring_change": "fused_ring",
    "changed_bonds": [
        {"product_atoms": [0, 1], "precursor_atoms": [0, 1], "event": "formed"}
    ],
    "preserved_anchors": ["anchor one", "anchor two"],
    "mechanistic_evidence": ["tether and closure mechanism"],
    "family_evidence": {
        "same_precursor_tether": "mapped tether",
        "new_ring_bond_atom_source": "named source atoms"
    },
    "risk_tags": ["ring_construction", "scaffold_edit"],
})
```
`validation.observations.atom_mapping` is an MCS heuristic, not a reaction oracle. Use it to challenge, then reconcile it with the declared mechanism, anchors, atom sources, and independent structure evidence.

Carbene, radical, or radical-ion precursors normally create `proof_required`, not automatic contradiction. Correct an unintended placeholder to a closed-shell precursor, or prove generation, lifetime, atom source, mechanism, and selectivity. Independent atom-balance, site, or topology contradictions still block. Standalone Li, Mg, Zn, or Cu may be `elemental_metal_reagent`; do not mistake representation electrons for an unsupported molecular radical.

For organometallic normalization, keep `current_reagent` in this event; `upstream_source_precursors` belong to a possible separate preparation and do not replace the validated reactant.

The EAS audit is deliberately narrow: declared nitration with one unambiguous new aryl nitro group, one uniquely locatable six-membered carbocyclic aromatic ring, a same-core site map, and complete preserved-substituent rule coverage. When applicable it exposes `validation.observations.eas_site_selectivity`:

- `supported`: heuristic support, not measured regioselectivity;
- `competing`: warning that another site is plausible;
- `unsupported`: revise the site or supply real substrate-specific evidence;
- `inconclusive`: ring, site, or rule coverage is insufficient, not chemical disproof;
- absent: the audit did not apply.

Inspect `reason_code`, `candidate_sites`, optional `required_evidence`, and knowledge references. Do not generalize to non-nitration EAS, denitration, fused/heteroaromatic rings, multiple ambiguous nitro changes, or partial rule coverage; a risk tag alone does not activate the audit. Atom mapping and a directing-rule hit never prove selectivity.
## 13. Terminal, Continuation, And Variant Review

Before accepting an advanced terminal, test whether a bounded 1-3-step mini-route reaches simpler, stable, purchasable, or more rational precursors.

Low complexity does not justify automatic acceptance when stereochemistry is assigned; review its source and a credible rollback. Consider:

- FG rollback and oxidation-state adjustment;
- explicit protection/deprotection;
- regio-, chemo-, or stereoselectivity sources;
- physical-organic reactivity differences;
- scaffold assembly or non-final electronic states;
- one credible next event from a route sketch.

Executable rescue steps return to the normal validation and commit path. For
advanced-terminal review, follow:
`route_sketch(terminal_review=True) -> propose_action(route_sketch_id=...) -> try_action -> sandbox_list`.
A listed action alone, or registration without `try_action`, does not satisfy the gate. Commit a credible custom event, then `next`; if rejected, accept only with a chemical and route-quality reason.

With multiple `continuation_steps`, sandbox alone cannot close the mini-route. Commit the credible first event to create the continuation; use force only if no event can be defined or repaired. Tool/template failure is insufficient.

Accept an ordinary simple, credible terminal with:

```python
cmd.execute("accept", {"reason": "specific chemical reason this is a credible starting material"})
```

Forced advanced form:
```python
cmd.execute("accept", {
    "reason": "specific advanced-terminal chemistry rationale",
    "force_accept_without_rescue": True,
    "rescue_not_actionable_reason": "why no real mini-route event is actionable",
})
```
`skip` is not chemistry success or an accepted starting material. Use it only
with a non-empty explicit reason when the item cannot be processed in scope:
```python
cmd.execute("skip", {"reason": "specific reason"})
```
Rachel records `outcome="skipped"`. The planning stop remains unresolved and
blocks `finalize`, `audit`, `learning_review`, and `review_node`. Reopen it with
`review_terminal`, then resolve it through a real commit or a reasoned `accept`.

When `strategy_continuation_brief` is active, copy its public ID and inspect or close it with:
```python
cmd.execute("continuation_status", {})
cmd.execute("continuation_abort", {
    "continuation_id": continuation_id,
    "reason": "specific chemical reason",
})
```
Resolve the next event through the normal loop. Never finalize with a pending continuation. If abort restores an auto-terminal but a better strategy exists, use `review_terminal`.

Use `review_terminal` to decompose the same unexpanded terminal leaf further:
```python
cmd.execute("review_terminal", {
    "smiles": "TERMINAL_SMILES",
    "reason": "why the leaf should be reopened",
    "additional_steps": 5,
})
```
`reason` is required; `additional_steps` is a non-negative integer and must be positive when the budget is exhausted. Reopening revokes completion; resolve all leaves and `finalize` again.

Use `review_node` only on an explicitly finalized parent with no pending continuation. It creates an independent variant and leaves the parent unchanged:
```python
variant = cmd.execute("review_node", {
    "node_id": "mol_3",
    "reason": "test a better route hypothesis",
    "instruction": "Use a mild alternative while preserving the fused core.",
    "intent": "directive",
    "site_hint": "edit only the selected functional group",
    "constraints": ["preserve the fused core"],
    "variant_session_file": r"runs\variant\session.json",
    "additional_steps": 10,
})
```
Require `reason`, a new `variant_session_file`, and `node_id` or `smiles`; both selectors must identify the same node when supplied together. Guidance fields match `guide`: `intent`, `site_hint`, `reaction_hint`, `precursors`, `constraints`, and `terminal_hint`. Optional `review_seed_id` must belong to that node; the seed is copied before round pruning and remains guidance rather than executable chemistry. Supply `adaptation_reason` when a structured override materially changes the selected seed.

There is no public `mode`. Rachel infers `extend_terminal` for a terminal leaf and `replace_subtree` for the target/expanded intermediate. The current command instance rebinds after successful reload; a new process must construct `LoggedRetroCmd` from the returned `variant_session_file`. Variants inherit the pinned profile and fail on missing/mismatched packs.

Before terminal acceptance or route finalization, repeat the chemistry quality boundary: every committed edge must remain one real, atom-sourced, correctly located, topology-faithful, compatible, selective, stereochemically honest event. A shorter truthful route is better than deeper invented chemistry.
## 14. Pinned Knowledge Profiles

The base pack loads implicitly; external references compose only in this order:
```text
rachel.base@1.0.0 -> team pack -> project pack
```
Initialize only immutable `id@version` references:
```python
cmd.execute("init", {
    "target": "TARGET_SMILES",
    "knowledge_profile": ["team.acme@1.2.0", "project.alpha@3.0.0"],
    "knowledge_roots": ["knowledge_packs"],
})
```
The session pins exact pack IDs, versions, and content hashes. Reload, export, `review_terminal`, and `review_node` reproduce that selection; missing or changed packs fail loudly, and active sessions never hot-reload.

Packs are JSON only and cannot execute Python. External packs may add advisory or stricter knowledge but cannot weaken locked state-machine, atom-accounting, topology, contradiction, family, risk, or validation gates.

Public audit exposes only IDs, versions, hashes, and hit entry IDs. Never copy proprietary entry text or local pack paths into public reasoning or reports.
## 15. Controlled Outcomes And Experience Learning

`committed` means only that an action entered the planning tree. It does not mean the reaction succeeded experimentally.

Use `record_outcome` only for facts supplied by the user, an experiment or literature record, or an authorized expert. `outcome` must be `success`, `partial`, or `failure`; any yield or conversion must be numeric and within 0 to 100:
```python
cmd.execute("record_outcome", {
    "step_id": "rxn_1",
    "action_id": "custom:amide:1",
    "outcome": "partial",
    "conditions": {"solvent": "DMF", "temperature_c": 25, "time_h": 4},
    "yield_percent": 61.5,
    "conversion_percent": 78,
    "observations": "Incomplete consumption.",
    "evidence": [{"kind": "internal_run", "reference": "ELN-2026-0815-07"}],
})
```
Never infer experimental success from `clear`, template execution, commit, route finalization, model confidence, or a plausible mechanism.

Only after explicit `finalize`, build the reusable-evidence summary:
```python
review = cmd.execute("learning_review", {})
```
Inspect guidance, route deviations, selected and rejected actions, validation, experience-card hits, experimental outcomes, and `steps_without_experimental_outcome`.

An LLM may propose a draft only when the user requested reusable-knowledge work and the finalized review supports it:
```python
cmd.execute("propose_knowledge", {
    "target_pack_id": "team.acme",
    "resource": "prompt.experience_cards",
    "entry": {
        "id": "team.acme.exp.mild_amidation",
        "tags": ["amide_formation"],
        "activation": {"any_events": ["stage.commit"]},
        "one_line": "Check conversion before increasing temperature.",
        "action_prompt": "Compare time and conversion before raising heat.",
        "avoid": "Do not infer a temperature benefit from one partial run.",
    },
    "rationale": "Bounded conclusion from the reviewed route evidence.",
    "source_refs": [{"step_id": "rxn_1", "action_id": "custom:amide:1"}],
    "evidence": [
        {"kind": "route_review", "reference": "rxn_1"},
        {"kind": "experimental_outcome", "outcome_id": "outcome_0001"},
    ],
})
```
The command verifies the resource, real source steps/actions, target-pack namespace, and absence of `_knowledge`. Its result is `status=staging, active=false` and cannot change runtime behavior. Never add executable code or fabricate evidence.

The LLM may draft; only an authorized expert may approve or reject. Approval remains inactive. Only a newly published immutable pack version can be selected by a later `init`.

Use the narrowest resource:

| Need | Resource |
| --- | --- |
| Contextual check and avoid reminder | `prompt.experience_cards` |
| Stage/event/tag, guardrail, self-prompt, or command policy | `prompt.runtime_policy` |
| Executable reaction or retro transform | `chem.reactions` |
| Functional-group or site match | `chem.functional_groups` / `chem.reactive_sites` |
| Byproduct lookup | `chem.byproducts` |
| Known scaffold identification | `chem.known_scaffolds` |
| Protecting-group definition or orthogonality | `chem.protecting_groups` |
| Reaction-role requirement | `chem.reaction_roles` |
| Local capping suggestion | `chem.smart_cap_rules` |
| Family delta and evidence interpretation | `chem.reaction_families` |
| Narrow aromatic-nitration directing heuristic | `chem.eas_directing_rules` |
| Complexity/progress SMARTS | `chem.cs_smarts` |
| Compatibility, selectivity, danger, or structural alert | `chem.fg_compatibility`, `chem.selectivity_conflicts`, `chem.selectivity_reactivity`, `chem.dangerous_combos`, or `chem.structural_alerts` |

Before authoring, replacing, disabling, reviewing, or publishing any prompt, experience, reaction, SMARTS, CAP, family, or risk entry, stop route work and read `knowledge/AUTHORING_GUIDE.md` completely for exact schemas, evidence kinds, namespaces, scopes, conflict relations, and gates.
## 16. Pack Publication Boundary

The standalone CLI provides seven commands:
```powershell
python -m Rachel.knowledge.cli validate --pack "PATH_TO_PACK"
python -m Rachel.knowledge.cli diff --left "PATH_TO_OLD_PACK" --right "PATH_TO_NEW_PACK"
python -m Rachel.knowledge.cli conflicts --pack "PATH_TO_TEAM_PACK" --pack "PATH_TO_PROJECT_PACK"
python -m Rachel.knowledge.cli import-drafts --session "PATH_TO_SESSION_JSON" --workspace "PATH_TO_STAGING_JSON" --pack-id "team.id"
python -m Rachel.knowledge.cli approve --workspace "PATH_TO_STAGING_JSON" --draft-id "DRAFT_ID" --reviewer "EXPERT_ID" --reason "Scoped evidence accepted."
python -m Rachel.knowledge.cli reject --workspace "PATH_TO_STAGING_JSON" --draft-id "DRAFT_ID" --reviewer "EXPERT_ID" --reason "Evidence or scope is insufficient."
python -m Rachel.knowledge.cli publish --workspace "PATH_TO_STAGING_JSON" --pack-root "knowledge_packs" --version "NEW_VERSION" --kind team
```
Publication creates a new version and never overwrites an existing one. Append `--dependency team.id@version` once per dependency; the option may repeat.

Every entry requires real route evidence, source step/action references, and expert approval. Reactions, functional/reactive SMARTS, protecting groups, Smart CAP, families, EAS rules, CS SMARTS, and locked risk/compatibility resources additionally require a parseable structured reaction instance, structured validation result, and experimental, literature, or internal SOP evidence.

Publication blocks duplicate IDs, implicit overrides, locked modifications, opposite prompt effects in the same scope, conflicting risk severity, and family allowed/forbidden delta conflicts. Competing reactions on one substrate remain parallel actions. Use `coexist` or `supersedes` only when the semantic relationship is true; otherwise narrow scope or reject the draft.
## 17. Completion And Export

Finish every molecule and continuation before finalization. Do not use
`finalize` to turn unresolved molecules into terminals. Active context, queued
work, incomplete leaves, unresolved skips, or pending continuations make it
return `route_not_ready_for_finalize` with `completion_blockers`; the failed
call does not drain or rewrite route state. `queue_empty` is necessary but not
sufficient because a prior skip can remain unresolved. Require `next` to report
`queue_empty`, `continuation_status` to report no pending continuation, and the
final `finalize` call to accept the complete state:
```python
end = cmd.execute("next", {})
status = cmd.execute("status", {})
continuations = cmd.execute("continuation_status", {})
if any("error" in result for result in (end, status, continuations)):
    raise RuntimeError({"end": end, "status": status, "continuations": continuations})
if end.get("action") != "queue_empty" or status.get("pending_count") or continuations.get("pending_continuation_count"):
    raise RuntimeError({"end": end, "status": status, "continuations": continuations})
tree = cmd.execute("tree", {})
if "error" in tree:
    raise RuntimeError(tree)
finalized = cmd.execute("finalize", {"summary": "route-level conclusion and uncertainty"})
if "error" in finalized or finalized.get("tree", {}).get("status") != "complete":
    raise RuntimeError(finalized)
audit = cmd.execute("audit", {
    "scopes": ["buyability", "safety"],
})
if "error" in audit:
    raise RuntimeError(audit)
report = cmd.execute("report", {})
exported = cmd.execute("export", {
    "name": "target_slug",
    "output_dir": str(run / "export"),
})
if "error" in report or "error" in exported or exported.get("visualization_ok") is False:
    raise RuntimeError({"report": report, "export": exported})
if exported.get("post_route_audit_status", {}).get("state") != "current":
    raise RuntimeError(exported)
```
Typical exports include session and route-tree JSON/text, terminal and knowledge-profile JSON, Markdown/HTML/text reports, visualization JSON, molecule/reaction images, a vector `images/synthesis_tree.svg` route overview with PNG fallback, and current `post_route_audit.json`/`.md` sidecars. Text, Markdown, and HTML share one forward-route projection and show the final thesis once, complete step chemistry/risk, compact source-aware conditions, gate findings, outcomes, scouting/frontier summary, and a compact audit summary with sidecar links. Accepted terminals appear as starting materials; unresolved skips are excluded from `terminals.json` and starting-material counts and instead appear under `Unresolved Planning Stops` with their reasons. They omit sandbox, prompt events, experience-card IDs, task/round source IDs, repeated plan snapshots, and full provenance retained in `session.json`/`tree.json`. Complete audit evidence remains in the sidecars; safety steps use forward `Step N (rxn_id)`. HTML reuses `images/` assets rather than embedding repeated base64 copies. Its vertical route workspace defaults to fit width and supports scroll, zoom, pan, fit route, and fullscreen without changing `visualization.json` or route state. Exported `tree.json` still carries the complete export-only current audit snapshot without changing `RetrosynthesisTree` or `route_digest`.
When reusing an export directory, Rachel removes only its two audit sidecars
before writing; if no current audit exists, they remain absent and other files stay.
`output_dir` must not make the exported `session.json` resolve to the active
authoritative session file; Rachel rejects that exact collision before creating
or deleting any export artifact.
`B1` means PubChem identity evidence; `B2` adds vendor/depositor catalog rows but not live stock, price, lead time, region, or procurement confirmation. Safety `S0` is local screening and `S1` is source-attributed PubChem GHS evidence. The LLM `risk_assessment` remains attributed planning judgment and is not independently semantically validated.

For each step, the audit reports declared conditions, source classes, and core coverage gaps. `conditions_insufficient_for_process_assessment` means at least one of the seven core fields is missing, invalid, or explicitly unknown. `conditions_available_for_preliminary_process_review` means all seven core fields have valid non-placeholder `specified` values; optional `scale`, `concentration`, and `pressure` remain visible when supplied. Planning suggestions are still unverified, and the status never replaces exact SDS, laboratory SOP, calorimetry, exposure, scale-up, pressure, concentration, or equipment assessment.

## 18. Error Recovery

- A failed session save is reported as a persistence error. The command runner
  reloads the last complete authoritative `session.json` (or resets an
  unpersisted first `init`) before serving another command; if reload itself
  fails, inspect `session_recovery_error`. Validation and other read-only errors
  do not trigger a reload.
- `no active context`: call `next`.
- PowerShell quoting, redirection, or `????` output: use Python dictionaries and the UTF-8 here-string pattern from Section 5, never Bash heredoc syntax.
- summary code cannot find `sites` or `items`: use `site_reaction_map`.
- `precursors_preview` contains `*` or `source="smart_capping"`: complete and correct the structural suggestion as a real `propose_action`, then validate it; neither form is commit-ready evidence. Read the original Smart CAP rule label from `execution.heuristic_reaction_hint`, and treat both the label and fragments only as structural hints.
- `protecting_group_source_required`: add the actual donor reagent and revalidate; do not override missing atoms.
- detected carbene, radical, or radical-ion precursor: treat as `proof_required` unless independent contradictions exist; provide generation, lifetime, atom source, mechanism, and selectivity evidence or correct an unintended placeholder.
- `elemental_metal_reagent`: RDKit radical electrons on isolated Li/Mg/Zn/Cu are a representation observation, not automatically an unsupported radical.
- `inconclusive`: separate chemistry evidence gaps from tool limits.
- `step budget exhausted`: unresolved queue/tree state was retained; inspect `status`, then retry `next(additional_steps=N)` with a justified positive extension.
- `queue capacity exceeded`: the commit was rejected before route-tree mutation; reduce or revise the precursor set or raise the configured capacity before retrying.
- an aborted continuation does not resolve a prior skip for the same molecule; reopen it and finish with a real `commit` or reasoned `accept`.
- successful `commit`: read `step_id`; do not expect removed legacy `success`.
- any `{"error": ...}`: stop, fix the state or payload, rerun the command, and never claim completion from partial output.
## 19. Final Checklist

Before claiming a route is complete, verify:

- every commit is one real chemical event;
- every new key atom and bond has an explicit precursor or reagent source, and the reaction occurs at the declared site;
- ring size, fusion, bridge, spiro relation, and retained topology are correct;
- functional groups and temporary handles survive the conditions, with protection and deprotection represented as real steps when required;
- chemo-, regio-, site-, and stereoselectivity are explained without invented certainty;
- actual alternative actions and rejections remain auditable;
- terminal choices are chemically honest and not caused only by tool failure, and all continuations are resolved;
- route commits are not mislabeled as experimental outcomes;
- knowledge drafts remain inactive until expert review and immutable publish;
- profile ID, version, and hash remain pinned through reload and variants;
- `report` and `export` complete without hidden errors.
- the latest post-finalize audit is `current`, B1/B2 and S0/S1 are interpreted within their evidence limits, and missing process conditions remain explicit;

Rachel informs, challenges, validates, and records. The LLM or chemist designs and decides. If mechanism, atom source, site/topology, compatibility, selectivity, stereochemistry, or evidence does not close, revise the action instead of obeying a template or forcing route depth.
