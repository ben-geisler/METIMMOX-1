# Session summary — adversarial review, fixes, and issue triage

**Date:** 2026-07-19 to 2026-07-29
**Branch:** `refactor/simplify-legacy-scaffolding` → merged to `main` via PR #148
**Model:** Claude Fable 5

This session began as an adversarial review of the R + QMD codebase (building on the
prior Opus/ValidateHE report in `report_Opus/`) and grew into verifying an existing
refactor, fixing three newly found correctness defects, regenerating all caches,
opening and merging PR #148, and triaging the open-issue backlog.

---

## What we did

### 1. Verified the existing refactor (issue #146)

Re-derived the prior refactor's fixes (labelled `147-A1`…`A13`, `B1`…`B10`) against
the actual source rather than trusting commit messages. Confirmed the Opus report's
#1 "Fail" finding (OS/PFS independent-AIC fitting with no ordering constraint) was
genuinely resolved by ordering-constrained joint distribution selection. Confirmed
the ~2,989-line reduction (20,202 → 17,213 R + QMD source, excl. archive) and, via
snapshot comparison, that it changed **no** result.

### 2. Fixed three correctness defects the review found (issue #147)

These were missed by both the Opus review and static inspection; they surfaced only
when caches were regenerated and reports rendered.

- **A14 (`fc4ae5b`) — interaction-coefficient extraction returned 0%. CHANGES RESULTS.**
  `extract_interaction_coefficients()` matched `^crp.*:Rx`, but R names the term
  `Rx:crp` (`RxExperimental arm:crp1`). It matched nothing, warned, and continued —
  so the four `b_*_rx_*` rows were *absent* from every EVPPI table, not zero. Fixed to
  match structurally (order-agnostic). Result: EVPPI table went **20 → 26 rows**, and
  **`b_tmb_braf_rx_pfs` now carries 41.2% of total EVPI** — the second-largest
  parameter after `u_np`, previously invisible. Pre-existing defect, not an A13
  regression.
- **A15 (`95cdaf0`) — `scenario_effect.qmd` could not render.** It compared the cache
  `failed_draw_policy` against a stale string literal instead of calling
  `psa_failed_draw_policy()`; it rejected every fresh cache. Now calls the accessor.
- **A7 (`1bb4d53`) — enriched-population path silently clamped OS ≥ PFS.** Replaced the
  silent `pmin(pfs, os)` with a labelled assertion matching `model_fun`. Latent
  (0 violations under gamma/gamma), **no result change**.

### 3. Fixed the discount-cycle-length hardcode (issue #82)

- **A17 (`d2663a4`)** — `discount_weights()` computed years as `/ 52`, silently
  assuming weekly cycles. Now takes `cl` explicitly and uses `* cl`. Verified the
  deterministic base case is **exactly** unchanged (`max |dCost| = 0`, `max |dEff| = 0`);
  `x/52` vs `x*(1/52)` differ at ~1e-16 but cancel in the aggregated sums, so no cache
  regeneration was needed.

### 4. Added regression tests (12 → 14 test files)

- `test_interaction_coefficient_extraction.R`, `test_enriched_survival_ordering.R`,
  `test_report_contracts.R` (A16, 163 checks covering the report/cache contracts that
  A14/A15 slipped through), and `test_discount_weights_cycle_length.R`.
- Every new test was **mutation-verified**: reverting the corresponding fix makes it
  fail with a message naming the fix.

### 5. Regenerated all caches from corrected code

Sampling (5000 models, 5 benign singular-matrix failures now correctly substituted),
PSA, EVPPI, and scenario EVPPI. PSA diagnostics confirmed the A1–A4 robustness work is
live: **fallback 5/5000, all replaced, 0 dropped**. Deleted 537 MB of orphaned
Model-A/B/C caches (`*_focused`/`*_joint`/`*_separate`).

### 6. Snapshots and the bug-fix impact report

- Corrected the `#147` baseline snapshot's stale `issue_number` metadata (146 → 147),
  preserving its honest `git_commit` provenance and byte-identical results.
- Rewrote `bug_fix_impact.qmd` to compare every baseline/fixed snapshot pair plus a
  cumulative row; rendered to PDF **and** HTML and tracked both (with a scoped
  `.gitignore` negation for the PDF, since the repo otherwise ignores rendered reports).

### 7. PR #148 — created, structured, merged

One PR splitting **refactor (#146, no result change)** from **correctness fixes
(#147, one changes results)**, with a commit-label key decoding the `147-*` tags, the
snapshot-delta tables, and the poster flagged as knowingly-failing. Rewrote both issue
bodies to match the split. Merged to `main`.

### 8. Issue triage — 9 issues closed

- **Closed by PR #148:** #146, #147, plus #55, #50, #59 (verified already fixed by the
  branch, `main` vs `HEAD` diff), and #82 (fixed in A17).
- **#101** (P0, "3 models in para_models") closed as **superseded** — it asks for the
  Model A/B/C design #146 retired.
- **8 more closed as already-resolved-in-main**, each with a code-evidence comment:
  #15, #51, #52, #53, #54, #57, #61, #63.

---

## What we did NOT do

### Left failing / open by deliberate decision

- **`poster_biomarker_correlation.qmd` — left intentionally failing.** It hard-codes a
  weibull-OS/lognormal-PFS pair that violates OS ≥ PFS at **251 of 521 time points**
  (the original WB-001 defect). `model_fun` now refuses it rather than clamping.
  Fixing it means porting `curve_draw()` to the gamma/gamma pair, which changes the
  poster's published numbers — a stats decision left to the author.
- **Issue #145 — left open.** Its code fixes are already in `main`; its bug #1 (OS/PFS
  clamp) still has the poster instance above. Deliberately excluded from PR #148.
- **Issue #124 — left open, needs a spec rewrite, not implementation.** Its figure spec
  references retired Model A/B/C panels; the underlying figures exist but the spec is
  stale against the single-model design.

### Flagged for the author's judgement, not changed

- **Treatment-schedule week indexing** (issue #49 / "A12"): schedules are hard-coded to
  a 520-week horizon, and the week-vs-cycle indexing (weeks 4,6,12… vs the documented
  5,7,13…) needs resolving against the trial protocol.
- **12-vs-13-week schedules**: CT every 12 weeks vs follow-up every 13 — may be two
  legitimate schedules; not assumed a bug.
- **Monitoring without visits** after week 39; and the manuscript EVPPI section, which
  now needs reading against the new `b_tmb_braf_rx_pfs` = 41.2% result.

### Open issues not addressed (14 remain)

- **Code:** #36 (P1 — visit costs from Eline's DRG weights; would change results;
  blocked on the DRG numbers), #60 (OS[1]==1 not asserted), #56 (CV values
  undocumented), #49, #58 (identity comparisons remain), #66, #16, #6, #137 (kNN EVPPI
  unvalidated against `voi`/`BCEA` — more relevant now that A14 made interaction EVPPI
  real).
- **Clinical/research:** #29, #31 (Norwegian EQ-5D tariffs), #33, #34.

---

## Key results and caveats

- **Deterministic base case unchanged throughout:** control €21,055.79 / 1.3498 QALY,
  CRP €53,790.65 / 1.3886, TMB/BRAF €63,235.05 / 1.3514. CRP ICER €843,973/QALY.
- **The one result that moved (A14):** the EVPPI/value-of-information table now includes
  the treatment-by-biomarker interactions; `b_tmb_braf_rx_pfs` at 41.2% of EVPI was
  absent from all prior results and **should be reviewed before entering the manuscript.**
- **PSA means shifted slightly** (control −€99/−0.0047 QALY, CRP −€65/−0.0027,
  TMB/BRAF −€155/−0.0072) — the expected fingerprint of the sampling/seeding/replacement
  fixes touching only the stochastic layer.
- **14/14 focused tests pass.**
