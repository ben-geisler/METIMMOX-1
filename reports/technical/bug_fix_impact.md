# Bug Fix Impact Analysis
Ben Geisler
2026-09-08

- [Snapshot Inventory](#snapshot-inventory)
- [Impact Comparisons](#impact-comparisons)
  - [Issue 145: baseline vs fixed](#issue-145-baseline-vs-fixed)
  - [Issue 146: baseline vs fixed](#issue-146-baseline-vs-fixed)
  - [Issue 147: baseline vs fixed](#issue-147-baseline-vs-fixed)
  - [Issue 149: baseline vs fixed](#issue-149-baseline-vs-fixed)
  - [Issue 151: baseline vs fixed](#issue-151-baseline-vs-fixed)
  - [Issue 152: baseline vs fixed](#issue-152-baseline-vs-fixed)
  - [Issue 153: baseline vs fixed](#issue-153-baseline-vs-fixed)
  - [Issue 154: baseline vs fixed](#issue-154-baseline-vs-fixed)
  - [Issue 155: baseline vs fixed](#issue-155-baseline-vs-fixed)
  - [Issue 156: baseline vs fixed](#issue-156-baseline-vs-fixed)
  - [Issue 157: baseline vs fixed](#issue-157-baseline-vs-fixed)
  - [Cumulative: issue 145 baseline vs issue 157
    fixed](#cumulative-issue-145-baseline-vs-issue-157-fixed)
- [Interpretation](#interpretation)

# Snapshot Inventory

| filename                          | issue | status   | commit  | mtime               |
|:----------------------------------|:------|:---------|:--------|:--------------------|
| snapshot_145_baseline_ce8ec08.rds | 145   | baseline | ce8ec08 | 2026-07-12 11:55:14 |
| snapshot_145_fixed_ea4431f.rds    | 145   | fixed    | ea4431f | 2026-07-12 19:54:11 |
| snapshot_146_baseline_a3befd2.rds | 146   | baseline | a3befd2 | 2026-07-27 13:57:12 |
| snapshot_146_fixed_7f76f9d.rds    | 146   | fixed    | 7f76f9d | 2026-07-27 13:57:12 |
| snapshot_147_baseline_7f76f9d.rds | 147   | baseline | 7f76f9d | 2026-07-27 13:57:12 |
| snapshot_147_fixed_d1f203f.rds    | 147   | fixed    | d1f203f | 2026-07-27 13:57:12 |
| snapshot_149_baseline_8f21e24.rds | 149   | baseline | 8f21e24 | 2026-09-08 08:51:07 |
| snapshot_149_fixed_4140474.rds    | 149   | fixed    | 4140474 | 2026-09-08 08:51:07 |
| snapshot_151_baseline_89d7a28.rds | 151   | baseline | 89d7a28 | 2026-09-08 08:51:07 |
| snapshot_151_fixed_70d672d.rds    | 151   | fixed    | 70d672d | 2026-09-08 08:51:07 |
| snapshot_152_baseline_e27b60f.rds | 152   | baseline | e27b60f | 2026-09-08 08:51:07 |
| snapshot_152_fixed_9cb460b.rds    | 152   | fixed    | 9cb460b | 2026-09-08 08:51:07 |
| snapshot_153_baseline_66a995c.rds | 153   | baseline | 66a995c | 2026-09-08 08:51:07 |
| snapshot_153_fixed_3b291ff.rds    | 153   | fixed    | 3b291ff | 2026-09-08 08:51:07 |
| snapshot_154_baseline_3d5938d.rds | 154   | baseline | 3d5938d | 2026-09-08 08:51:07 |
| snapshot_154_fixed_1f7202c.rds    | 154   | fixed    | 1f7202c | 2026-09-08 08:51:07 |
| snapshot_155_baseline_331a1a2.rds | 155   | baseline | 331a1a2 | 2026-09-08 08:51:07 |
| snapshot_155_fixed_9e29ea6.rds    | 155   | fixed    | 9e29ea6 | 2026-09-08 08:51:07 |
| snapshot_156_baseline_2806457.rds | 156   | baseline | 2806457 | 2026-09-08 08:51:07 |
| snapshot_156_fixed_8a89056.rds    | 156   | fixed    | 8a89056 | 2026-09-08 08:51:07 |
| snapshot_157_baseline_efcd775.rds | 157   | baseline | efcd775 | 2026-09-08 08:51:07 |
| snapshot_157_fixed_cdf852f.rds    | 157   | fixed    | cdf852f | 2026-09-08 08:51:07 |

Available single-model snapshots

# Impact Comparisons

| Comparison | Before | After |
|:---|:---|:---|
| Issue 145: baseline vs fixed | snapshot_145_baseline_ce8ec08.rds | snapshot_145_fixed_ea4431f.rds |
| Issue 146: baseline vs fixed | snapshot_146_baseline_a3befd2.rds | snapshot_146_fixed_7f76f9d.rds |
| Issue 147: baseline vs fixed | snapshot_147_baseline_7f76f9d.rds | snapshot_147_fixed_d1f203f.rds |
| Issue 149: baseline vs fixed | snapshot_149_baseline_8f21e24.rds | snapshot_149_fixed_4140474.rds |
| Issue 151: baseline vs fixed | snapshot_151_baseline_89d7a28.rds | snapshot_151_fixed_70d672d.rds |
| Issue 152: baseline vs fixed | snapshot_152_baseline_e27b60f.rds | snapshot_152_fixed_9cb460b.rds |
| Issue 153: baseline vs fixed | snapshot_153_baseline_66a995c.rds | snapshot_153_fixed_3b291ff.rds |
| Issue 154: baseline vs fixed | snapshot_154_baseline_3d5938d.rds | snapshot_154_fixed_1f7202c.rds |
| Issue 155: baseline vs fixed | snapshot_155_baseline_331a1a2.rds | snapshot_155_fixed_9e29ea6.rds |
| Issue 156: baseline vs fixed | snapshot_156_baseline_2806457.rds | snapshot_156_fixed_8a89056.rds |
| Issue 157: baseline vs fixed | snapshot_157_baseline_efcd775.rds | snapshot_157_fixed_cdf852f.rds |
| Cumulative: issue 145 baseline vs issue 157 fixed | snapshot_145_baseline_ce8ec08.rds | snapshot_157_fixed_cdf852f.rds |

Impact comparisons included in this report

Each comparison below reports the change in base case costs and QALYs,
and in net monetary benefit at the willingness-to-pay threshold, from
the *before* snapshot to the *after* snapshot. Where both snapshots
carry the EVPPI cache (snapshots taken from issue \#152 onwards), the
per-patient EVPI and the EVPPI table are compared as well; the *after*
column includes the Monte Carlo standard error of the regression
estimator. From issue \#156 onwards each snapshot records the md5 and
modification time of the PSA cache it copied, and a *fixed* snapshot
regenerates the PSA and EVPPI caches when they predate the fix commit;
the provenance line under each comparison states whether the two PSA
files are the same cache.

## Issue 145: baseline vs fixed

Table: Base case impact – Issue 145: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 20,867\| 1.4141\| EUR 21,056\| 1.3498\| EUR
189\| -0.0644\| \|CRP-guided \| EUR 55,026\| 1.3842\| EUR 53,791\|
1.3886\| -EUR 1,235\| 0.0043\| \|TMB/BRAF-guided \| EUR 60,798\|
1.3530\| EUR 63,235\| 1.3514\| EUR 2,438\| -0.0017\|

Table: Net monetary benefit impact – Issue 145: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 51,254\| EUR
47,783\| -EUR 3,471\| \|CRP-guided \| EUR 15,570\| EUR 17,026\| EUR
1,456\| \|TMB/BRAF-guided \| EUR 8,207\| EUR 5,684\| -EUR 2,523\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 146: baseline vs fixed

Table: Base case impact – Issue 146: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,056\| 1.3498\| EUR 21,056\| 1.3498\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,791\| 1.3886\| EUR 53,791\| 1.3886\|
EUR 0\| 0.0000\| \|TMB/BRAF-guided \| EUR 63,235\| 1.3514\| EUR 63,235\|
1.3514\| EUR 0\| 0.0000\|

Table: Net monetary benefit impact – Issue 146: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 47,783\| EUR
47,783\| EUR 0\| \|CRP-guided \| EUR 17,026\| EUR 17,026\| EUR 0\|
\|TMB/BRAF-guided \| EUR 5,684\| EUR 5,684\| EUR 0\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 147: baseline vs fixed

Table: Base case impact – Issue 147: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,056\| 1.3498\| EUR 21,056\| 1.3498\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,791\| 1.3886\| EUR 53,791\| 1.3886\|
EUR 0\| -0.0000\| \|TMB/BRAF-guided \| EUR 63,235\| 1.3514\| EUR
63,235\| 1.3514\| EUR 0\| -0.0000\|

Table: Net monetary benefit impact – Issue 147: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 47,783\| EUR
47,783\| EUR 0\| \|CRP-guided \| EUR 17,026\| EUR 17,026\| EUR 0\|
\|TMB/BRAF-guided \| EUR 5,684\| EUR 5,684\| EUR 0\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 149: baseline vs fixed

Table: Base case impact – Issue 149: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,056\| 1.3498\| EUR 21,523\| 1.3714\| EUR
468\| 0.0216\| \|CRP-guided \| EUR 53,791\| 1.3886\| EUR 53,928\|
1.3918\| EUR 137\| 0.0033\| \|TMB/BRAF-guided \| EUR 63,235\| 1.3514\|
EUR 63,399\| 1.3609\| EUR 164\| 0.0095\|

Table: Net monetary benefit impact – Issue 149: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 47,783\| EUR
48,419\| EUR 636\| \|CRP-guided \| EUR 17,026\| EUR 17,056\| EUR 30\|
\|TMB/BRAF-guided \| EUR 5,684\| EUR 6,005\| EUR 320\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 151: baseline vs fixed

Table: Base case impact – Issue 151: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,523\| 1.3714\| EUR 21,523\| 1.3714\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,928\| 1.3918\| EUR 53,928\| 1.3918\|
EUR 0\| 0.0000\| \|TMB/BRAF-guided \| EUR 63,399\| 1.3609\| EUR 63,399\|
1.3609\| EUR 0\| 0.0000\|

Table: Net monetary benefit impact – Issue 151: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 48,419\| EUR
48,419\| EUR 0\| \|CRP-guided \| EUR 17,056\| EUR 17,056\| EUR 0\|
\|TMB/BRAF-guided \| EUR 6,005\| EUR 6,005\| EUR 0\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

## Issue 152: baseline vs fixed

Table: Base case impact – Issue 152: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,523\| 1.3714\| EUR 21,523\| 1.3714\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,928\| 1.3918\| EUR 53,928\| 1.3918\|
EUR 0\| 0.0000\| \|TMB/BRAF-guided \| EUR 63,399\| 1.3609\| EUR 63,399\|
1.3609\| EUR 0\| 0.0000\|

Table: Net monetary benefit impact – Issue 152: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 48,419\| EUR
48,419\| EUR 0\| \|CRP-guided \| EUR 17,056\| EUR 17,056\| EUR 0\|
\|TMB/BRAF-guided \| EUR 6,005\| EUR 6,005\| EUR 0\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 45.20 before, EUR 45.20 after.

Table: EVPPI impact – Issue 152: baseline vs fixed \|Parameter / group
\| EVPPI Before\| % EVPI Before\| EVPPI After\| SE After\| % EVPI
After\| \|:—————————————\|————:\|————-:\|———–:\|——–:\|————:\| \|All
costs (group) \| EUR 45.20\| 100.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|Drug costs (group) \| EUR 45.20\| 100.0%\| EUR 0.00\| EUR 0.00\|
0.0%\| \|TMB/BRAF-treatment interaction (group) \| EUR 45.20\| 100.0%\|
EUR 0.00\| EUR 0.00\| 0.0%\| \|Resource-use costs (group) \| EUR 45.20\|
100.0%\| EUR 0.00\| EUR 0.04\| 0.0%\| \|b_tmb_braf_rx_os \| EUR 45.20\|
100.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|b_tmb_braf_rx_pfs \| EUR 45.20\|
100.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|c_drug_FLOX \| EUR 45.20\|
100.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|c_drug_nivo \| EUR 45.20\|
100.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|u_np \| EUR 45.20\| 100.0%\| EUR
0.00\| EUR 0.00\| 0.0%\| \|c_test_CT \| EUR 43.76\| 96.8%\| EUR 0.00\|
EUR 0.00\| 0.0%\| \|c_other_last \| EUR 41.57\| 92.0%\| EUR 0.00\| EUR
0.00\| 0.0%\| \|p_tmb_braf \| EUR 40.09\| 88.7%\| EUR 0.00\| EUR 0.00\|
0.0%\| \|CRP-treatment interaction (group) \| EUR 0.00\| 0.0%\| EUR
0.00\| EUR 0.74\| 0.0%\| \|Biomarker prevalence (group) \| EUR 0.00\|
0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|Test costs (group) \| EUR 0.00\|
0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|Utilities (group) \| EUR 0.00\|
0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|b_crp_rx_os \| EUR 0.00\| 0.0%\|
EUR 0.10\| EUR 1.60\| 0.2%\| \|b_crp_rx_pfs \| EUR 0.00\| 0.0%\| EUR
0.00\| EUR 0.12\| 0.0%\| \|c_other_baseline \| EUR 0.00\| 0.0%\| EUR
0.00\| EUR 0.00\| 0.0%\| \|c_other_follow \| EUR 0.00\| 0.0%\| EUR
0.00\| EUR 0.00\| 0.0%\| \|c_other_visit \| EUR 0.00\| 0.0%\| EUR 0.00\|
EUR 0.00\| 0.0%\| \|c_test_CRP \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\| \|c_test_NGS \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\|
0.0%\| \|c_test_blood \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|p_crp \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|u_p \| EUR
0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|Biomarker-treatment
interaction (group) \| –\| –\| EUR 0.00\| EUR 1.61\| 0.0%\|

## Issue 153: baseline vs fixed

Table: Base case impact – Issue 153: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,523\| 1.3714\| EUR 21,523\| 1.3714\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,928\| 1.3918\| EUR 53,928\| 1.3918\|
EUR 0\| 0.0000\| \|TMB/BRAF-guided \| EUR 63,399\| 1.3609\| EUR 63,399\|
1.3609\| EUR 0\| 0.0000\|

Table: Net monetary benefit impact – Issue 153: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 48,419\| EUR
48,419\| EUR 0\| \|CRP-guided \| EUR 17,056\| EUR 17,056\| EUR 0\|
\|TMB/BRAF-guided \| EUR 6,005\| EUR 6,005\| EUR 0\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 45.20 before, EUR 45.20 after.

Table: EVPPI impact – Issue 153: baseline vs fixed \|Parameter / group
\| EVPPI Before\| % EVPI Before\| EVPPI After\| SE After\| % EVPI
After\| \|:—————————————\|————:\|————-:\|———–:\|——–:\|————:\|
\|b_crp_rx_os \| EUR 0.10\| 0.2%\| EUR 0.10\| EUR 1.60\| 0.2%\| \|All
costs (group) \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|Drug
costs (group) \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|Biomarker-treatment interaction (group) \| EUR 0.00\| 0.0%\| EUR
0.00\| EUR 1.61\| 0.0%\| \|CRP-treatment interaction (group) \| EUR
0.00\| 0.0%\| EUR 0.00\| EUR 0.74\| 0.0%\| \|TMB/BRAF-treatment
interaction (group) \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|Resource-use costs (group) \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.04\|
0.0%\| \|Biomarker prevalence (group) \| EUR 0.00\| 0.0%\| EUR 0.00\|
EUR 0.00\| 0.0%\| \|Test costs (group) \| EUR 0.00\| 0.0%\| EUR 0.00\|
EUR 0.00\| 0.0%\| \|Utilities (group) \| EUR 0.00\| 0.0%\| EUR 0.00\|
EUR 0.00\| 0.0%\| \|b_crp_rx_pfs \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.12\| 0.0%\| \|b_tmb_braf_rx_os \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\| \|b_tmb_braf_rx_pfs \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\| \|c_drug_FLOX \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\|
0.0%\| \|c_drug_nivo \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_other_baseline \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_other_follow \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_other_last \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_other_visit \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_test_CRP \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_test_CT \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_test_NGS \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_test_blood \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|p_crp
\| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|p_tmb_braf \| EUR
0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|u_np \| EUR 0.00\| 0.0%\|
EUR 0.00\| EUR 0.00\| 0.0%\| \|u_p \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\|

## Issue 154: baseline vs fixed

Table: Base case impact – Issue 154: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,523\| 1.3714\| EUR 21,523\| 1.3714\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,928\| 1.3918\| EUR 53,928\| 1.3918\|
EUR 0\| 0.0000\| \|TMB/BRAF-guided \| EUR 63,399\| 1.3609\| EUR 63,399\|
1.3609\| EUR 0\| 0.0000\|

Table: Net monetary benefit impact – Issue 154: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 48,419\| EUR
48,419\| EUR 0\| \|CRP-guided \| EUR 17,056\| EUR 17,056\| EUR 0\|
\|TMB/BRAF-guided \| EUR 6,005\| EUR 6,005\| EUR 0\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 45.20 before, EUR 22.13 after.

Table: EVPPI impact – Issue 154: baseline vs fixed \|Parameter / group
\| EVPPI Before\| % EVPI Before\| EVPPI After\| SE After\| % EVPI
After\| \|:—————————————\|————:\|————-:\|———–:\|——–:\|————:\|
\|b_crp_rx_os \| EUR 0.10\| 0.2%\| EUR 0.09\| EUR 0.95\| 0.4%\| \|All
costs (group) \| EUR 0.00\| 0.0%\| –\| –\| –\| \|Drug costs (group) \|
EUR 0.00\| 0.0%\| –\| –\| –\| \|Biomarker-treatment interaction (group)
\| EUR 0.00\| 0.0%\| EUR 1.18\| EUR 1.94\| 5.3%\| \|CRP-treatment
interaction (group) \| EUR 0.00\| 0.0%\| EUR 0.04\| EUR 0.41\| 0.2%\|
\|TMB/BRAF-treatment interaction (group) \| EUR 0.00\| 0.0%\| EUR 0.00\|
EUR 0.00\| 0.0%\| \|Resource-use costs (group) \| EUR 0.00\| 0.0%\| EUR
0.00\| EUR 0.00\| 0.0%\| \|Biomarker prevalence (group) \| EUR 0.00\|
0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|Test costs (group) \| EUR 0.00\|
0.0%\| –\| –\| –\| \|Utilities (group) \| EUR 0.00\| 0.0%\| EUR 0.00\|
EUR 0.00\| 0.0%\| \|b_crp_rx_pfs \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.08\| 0.0%\| \|b_tmb_braf_rx_os \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\| \|b_tmb_braf_rx_pfs \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\| \|c_drug_FLOX \| EUR 0.00\| 0.0%\| –\| –\| –\|
\|c_drug_nivo \| EUR 0.00\| 0.0%\| –\| –\| –\| \|c_other_baseline \| EUR
0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|c_other_follow \| EUR
0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|c_other_last \| EUR 0.00\|
0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|c_other_visit \| EUR 0.00\| 0.0%\|
EUR 0.00\| EUR 0.00\| 0.0%\| \|c_test_CRP \| EUR 0.00\| 0.0%\| –\| –\|
–\| \|c_test_CT \| EUR 0.00\| 0.0%\| –\| –\| –\| \|c_test_NGS \| EUR
0.00\| 0.0%\| –\| –\| –\| \|c_test_blood \| EUR 0.00\| 0.0%\| –\| –\|
–\| \|p_crp \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|p_tmb_braf \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|u_np \|
EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|u_p \| EUR 0.00\|
0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|u_decrement \| –\| –\| EUR 0.00\|
EUR 0.00\| 0.0%\|

## Issue 155: baseline vs fixed

Table: Base case impact – Issue 155: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,523\| 1.3714\| EUR 21,523\| 1.3714\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,928\| 1.3918\| EUR 53,928\| 1.3918\|
EUR 0\| 0.0000\| \|TMB/BRAF-guided \| EUR 63,399\| 1.3609\| EUR 63,399\|
1.3609\| EUR 0\| 0.0000\|

Table: Net monetary benefit impact – Issue 155: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 48,419\| EUR
48,419\| EUR 0\| \|CRP-guided \| EUR 17,056\| EUR 17,056\| EUR 0\|
\|TMB/BRAF-guided \| EUR 6,005\| EUR 6,005\| EUR 0\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after snapshot predates issue \#156 and did not
record the PSA cache md5. Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 22.13 before, EUR 22.13 after.

Table: EVPPI impact – Issue 155: baseline vs fixed \|Parameter / group
\| EVPPI Before\| % EVPI Before\| EVPPI After\| SE After\| % EVPI
After\| \|:—————————————\|————:\|————-:\|———–:\|——–:\|————:\|
\|Biomarker-treatment interaction (group) \| EUR 1.18\| 5.3%\| EUR
1.18\| EUR 1.94\| 5.3%\| \|b_crp_rx_os \| EUR 0.09\| 0.4%\| EUR 0.09\|
EUR 0.95\| 0.4%\| \|CRP-treatment interaction (group) \| EUR 0.04\|
0.2%\| EUR 0.04\| EUR 0.41\| 0.2%\| \|TMB/BRAF-treatment interaction
(group) \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|Resource-use
costs (group) \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|Biomarker prevalence (group) \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\| \|Utilities (group) \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\| \|b_crp_rx_pfs \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.08\|
0.0%\| \|b_tmb_braf_rx_os \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\|
0.0%\| \|b_tmb_braf_rx_pfs \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\|
0.0%\| \|c_other_baseline \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\|
0.0%\| \|c_other_follow \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\|
0.0%\| \|c_other_last \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|c_other_visit \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\|
\|p_crp \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|p_tmb_braf
\| EUR 0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|u_decrement \| EUR
0.00\| 0.0%\| EUR 0.00\| EUR 0.00\| 0.0%\| \|u_np \| EUR 0.00\| 0.0%\|
EUR 0.00\| EUR 0.00\| 0.0%\| \|u_p \| EUR 0.00\| 0.0%\| EUR 0.00\| EUR
0.00\| 0.0%\|

## Issue 156: baseline vs fixed

Table: Base case impact – Issue 156: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,523\| 1.3714\| EUR 21,523\| 1.3714\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,928\| 1.3918\| EUR 53,928\| 1.3918\|
EUR 0\| 0.0000\| \|TMB/BRAF-guided \| EUR 63,399\| 1.3609\| EUR 63,399\|
1.3609\| EUR 0\| 0.0000\|

Table: Net monetary benefit impact – Issue 156: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 48,419\| EUR
48,419\| EUR 0\| \|CRP-guided \| EUR 17,056\| EUR 17,056\| EUR 0\|
\|TMB/BRAF-guided \| EUR 6,005\| EUR 6,005\| EUR 0\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after PSA cache md5 fa17a8dab004ee9fc4c850ff9897b656
(modified 2026-09-06 23:39). Identity of the two PSA files cannot be
established from metadata.

Per-patient EVPI: EUR 22.13 before, EUR 0.00 after.

## Issue 157: baseline vs fixed

Table: Base case impact – Issue 157: baseline vs fixed \|Strategy \|
Cost Before\| QALYs Before\| Cost After\| QALYs After\| Cost Change\|
QALY Change\| \|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\|
\|Standard of Care \| EUR 21,523\| 1.3714\| EUR 21,523\| 1.3714\| EUR
0\| 0.0000\| \|CRP-guided \| EUR 53,928\| 1.3918\| EUR 53,928\| 1.3918\|
EUR 0\| 0.0000\| \|TMB/BRAF-guided \| EUR 63,399\| 1.3609\| EUR 63,399\|
1.3609\| EUR 0\| 0.0000\|

Table: Net monetary benefit impact – Issue 157: baseline vs fixed
\|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 48,419\| EUR
48,419\| EUR 0\| \|CRP-guided \| EUR 17,056\| EUR 17,056\| EUR 0\|
\|TMB/BRAF-guided \| EUR 6,005\| EUR 6,005\| EUR 0\|

PSA cache md5 before fa17a8dab004ee9fc4c850ff9897b656 (modified
2026-09-06 23:39), after a057d7d7856e7c3bc7602058118d9a92 (modified
2026-09-07 09:28). The after PSA was regenerated inside the snapshot run
because the cache predated the fix commit.

Per-patient EVPI: EUR 0.00 before, EUR 0.00 after.

## Cumulative: issue 145 baseline vs issue 157 fixed

Table: Base case impact – Cumulative: issue 145 baseline vs issue 157
fixed \|Strategy \| Cost Before\| QALYs Before\| Cost After\| QALYs
After\| Cost Change\| QALY Change\|
\|:—————-\|———–:\|————:\|———-:\|———–:\|———–:\|———–:\| \|Standard of Care
\| EUR 20,867\| 1.4141\| EUR 21,523\| 1.3714\| EUR 656\| -0.0427\|
\|CRP-guided \| EUR 55,026\| 1.3842\| EUR 53,928\| 1.3918\| -EUR 1,098\|
0.0076\| \|TMB/BRAF-guided \| EUR 60,798\| 1.3530\| EUR 63,399\|
1.3609\| EUR 2,602\| 0.0078\|

Table: Net monetary benefit impact – Cumulative: issue 145 baseline vs
issue 157 fixed \|Strategy \| NMB Before\| NMB After\| NMB Change\|
\|:—————-\|———-:\|———-:\|———-:\| \|Standard of Care \| EUR 51,254\| EUR
48,419\| -EUR 2,835\| \|CRP-guided \| EUR 15,570\| EUR 17,056\| EUR
1,486\| \|TMB/BRAF-guided \| EUR 8,207\| EUR 6,005\| -EUR 2,203\|

PSA provenance: before snapshot predates issue \#156 and did not record
the PSA cache md5; after PSA cache md5 a057d7d7856e7c3bc7602058118d9a92
(modified 2026-09-07 09:28, regenerated inside the snapshot run).
Identity of the two PSA files cannot be established from metadata.

# Interpretation

This report compares baseline and fixed single-model snapshots for each
bug-fix issue that has both, plus a cumulative comparison from the
earliest baseline to the latest fixed snapshot. Snapshots are produced
by the `13_save_snapshot.R` analysis script (in `analysis/`). If older
snapshots contain obsolete multi-model components, `load_snapshot()`
warns during render and the tables use the top-level single-model
results.

------------------------------------------------------------------------

**Report completed on:** 2026-09-08  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.2
