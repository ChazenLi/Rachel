# Site Anchor Audit Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Block hand-written same-core precursor commits when substituent anchor positions drift across a preserved scaffold.

**Architecture:** Keep full structural comparison in a new read-only chemistry helper, keep sandbox/session payload compact via a micro-summary, and enforce the hard gate in `RetroSession.commit()` so failed attempts remain in sandbox for correction. Preserve existing flow and exports by only adding fields and checks.

**Tech Stack:** Python, RDKit, pytest

---

### Task 1: Add a read-only site-anchor audit helper

**Files:**
- Create: `E:\Python\skills\Rachel\chem_tools\site_audit.py`
- Reference: `E:\Python\skills\Rachel\chem_tools\mol_info.py`

**Steps:**
1. Add a helper that identifies the major scaffold-retaining precursor.
2. Compare scaffold attachment sites between current molecule and the major precursor.
3. Return a compact audit dict with `site_retentive`, `changed_site_count`, `pass`, and `summary`.

### Task 2: Preserve minimal validation signals in sandbox results

**Files:**
- Modify: `E:\Python\skills\Rachel\main\retro_orchestrator.py`

**Steps:**
1. Add a compact `validation_micro` extractor.
2. Extend `SandboxResult` with additive `validation_micro` and `site_audit` fields.
3. Populate those fields from `try_precursors()` without changing the existing flattened `forward_validation` behavior.

### Task 3: Enforce commit-time hard gate without clearing sandbox

**Files:**
- Modify: `E:\Python\skills\Rachel\main\retro_session.py`
- Modify: `E:\Python\skills\Rachel\main\retro_tree.py`

**Steps:**
1. Add a minimal helper in `RetroSession` to detect site-retention claims and evaluate the new gate.
2. Block `llm_proposed` site-retentive commits when the audit fails or when `template_attempted=false` with `scaffold_alignment < 1.0`.
3. Record `gate_failed` audit entries and keep sandbox state intact.
4. Add a new `site_audit_summary` field to `LLMDecision` for export visibility.

### Task 4: Add regression tests for the 095 failure mode

**Files:**
- Create: `E:\Python\skills\Rachel\tests\test_site_anchor_audit.py`

**Steps:**
1. Verify the new audit rejects the swapped-site indole precursor.
2. Verify the audit accepts the correct single-site rollback precursor.
3. Verify `RetroSession.commit()` blocks the bad hand-written precursor and preserves sandbox state.
