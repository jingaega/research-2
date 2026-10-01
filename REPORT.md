# REPORT — Classifying psychiatric disorders from pooled postmortem brain microarrays

Every number below comes from a JSON in `results/` (paths given). The full decision record is in
`LOG.md`, the plan in `PLAN.md`, the audit in `DATA_AUDIT.json` and the sources in `LITERATURE.md`.
The generated tables are in `results/report_tables.md`.

## 1. Answers in brief

1. **Is there signal? Modest, for schizophrenia only, and not cleanly separable from study
   composition.**
   - SCZ vs control, locked test: macro-F1 **0.653 [95% CI 0.543–0.750]** (chance 0.492), AUC 0.703.
   - The CV estimate agrees: 0.669 ± 0.065.
   - On the same test donors a model that only knows *which cohort a donor came from* scores
     **0.618 [0.505–0.714]**. In CV the gap is not significant under the corrected test.
   - What distinguishes the expression models is **transfer**: they keep SCZ signal in 3 of 5 unseen
     cohorts (macro-F1 0.66–0.71 for V6b, 0.63–0.74 for V6a), while the study-only model collapses there (0.22–0.44).
2. **BD vs control** reaches test macro-F1 0.636 [0.507–0.753], AUC 0.758, but **the signal depends on one
   cohort**. Without the Stanley donors, dev CV drops to chance (0.502 vs 0.486). In unseen cohorts it is
   at chance outside Pittsburgh.
3. **MDD vs control** reaches test macro-F1 0.623 [0.494–0.749]; the CI touches chance. In unseen cohorts it
   is at or below chance (0.41–0.53). **No transferable MDD signal.**
4. **Task B (SCZ vs BD vs MDD) is at chance**: test macro-F1 0.367 [0.250–0.476] (chance 0.331). The
   study-only model gets 0.374. **These data cannot support discrimination between disorders.**
5. **Graph structure does not help, for strong or weak contrasts.**
   - STRING interaction smoothing: −0.000 / −0.004 / −0.009 / −0.003 (SCZ / BD / MDD / Task B) in CV, and
     never better than a random graph with the same edge count.
   - Co-expression graphs: the same.
   - Curated Reactome pathways are *worse* than genes, and random gene sets of equal size beat them on
     3 of 4 tasks.
   - On the test set the graph effect (best model minus its no-graph twin) is +0.026, +0.013, +0.013 and
     −0.031, and every 95% CI includes 0.
   - The graph neural network was not run, because its pre-registered trigger (a cheaper graph variant
     beating its random control) was never met.
6. **Most of the SCZ signal is cell-type composition.** Ten cell-type indices recover 88% of the
   baseline's above-chance SCZ performance in CV. This could reflect disease, medication, agonal state
   or dissection; the data cannot tell them apart.

The plausible "modest SCZ signal, Task B at chance, no graph comparison survives" outcome is what happened.
It is a finding about these data, not a failed run.

## 2. Data (`DATA_AUDIT.json`)

- 5 cohorts after auditing 16 GEO series:
  - Stanley: GSE35978 (GSE12649 used for the audit only).
  - Pittsburgh: GSE53987 + GSE54567/68 (BA9), with GSE54571/72 used for the audit.
  - Charing Cross: GSE17612.
  - "Victoria": GSE21138 (provenance unverified).
  - Pritzker: GSE92538, which I added.
- Preprocessing: raw CEL files, single-array frozen RMA, and 11,955 genes shared by all platforms.
- Locked classes: **SCZ, BD, MDD**. Locked test set: 135 donors.
- Audit findings (all reported in DATA_AUDIT.json):
  - GSE11223 in the candidate list is colon tissue.
  - 129 of 569 donors appear in more than one series, including at least 20 Pittsburgh donors shared
    between two labs.
  - The GSE53987 hippocampal arrays are mislabelled relative to their donors.
  - Sex mismatches, mixed samples and identity conflicts give 7 "affected" donors, excluded from the
    primary analysis, plus 1 donor with an unresolvable label.

## 3. Final comparison (both final configurations, every task, all three evaluations)

Chance macro-F1 on test: uniform random / class-prior random / majority class.
- SCZ: 0.492 / 0.499 / 0.378.
- BD: 0.488 / 0.497 / 0.390.
- MDD: 0.487 / 0.498 / 0.394.
- Task B: 0.331 / 0.329 / 0.167.
- MCC chance is 0 and AUC chance is 0.5 throughout.

**Best model = V6b** (in-fold ANOVA top-1,000 genes → STRING graph smoothing → L2 logistic regression).
**Best simple baseline = V6a** (the same without the graph). C1 = study-only reference control. CV
values are mean ± fold SD; logistic regression is deterministic, so seed SD = 0. Test CIs are percentile
bootstraps over test donors (2,000 resamples).

| Task | Config | (1) Locked test macro-F1 [95% CI] | test AUC | (2) Dev 5×5 CV macro-F1 | (3) Leave-one-cohort-out macro-F1 |
|---|---|---|---|---|---|
| SCZ vs CTL | V6b | **0.653** [0.543, 0.750] | 0.703 | 0.669 ± 0.065 | CharingX 0.656, Pitt 0.684, Pritzker 0.516, Stanley 0.543, Victoria 0.712 |
| | V6a | **0.627** [0.518, 0.724] | 0.719 | 0.668 ± 0.072 | 0.631, 0.741, 0.516, 0.621, 0.664 |
| | C1 | 0.618 [0.505, 0.714] | 0.654 | 0.626 ± 0.013 | 0.309, 0.224, 0.439, 0.339, 0.344 |
| BD vs CTL | V6b | **0.636** [0.507, 0.753] | 0.758 | 0.568 ± 0.067 | Pitt 0.664, Pritzker 0.500, Stanley 0.522 |
| | V6a | **0.623** [0.490, 0.739] | 0.779 | 0.561 ± 0.050 | 0.622, 0.488, 0.496 |
| | C1 | 0.613 [0.484, 0.729] | 0.624 | 0.596 ± 0.030 | 0.403, 0.439, 0.323 |
| MDD vs CTL | V6b | **0.623** [0.494, 0.749] | 0.647 | 0.576 ± 0.053 | Pitt 0.532, Pritzker 0.453, Stanley 0.487 |
| | V6a | **0.610** [0.476, 0.733] | 0.645 | 0.592 ± 0.065 | 0.529, 0.438, 0.406 |
| | C1 | 0.597 [0.472, 0.721] | 0.628 | 0.595 ± 0.018 | 0.329, 0.401, 0.193 |
| Task B | V6b | **0.367** [0.250, 0.476] | 0.587 | 0.464 ± 0.075 | Pitt 0.443, Pritzker 0.386, Stanley 0.286 |
| | V6a | **0.398** [0.280, 0.509] | 0.600 | 0.467 ± 0.080 | 0.441, 0.406, 0.311 |
| | C1 | 0.374 [0.286, 0.456] | 0.646 | 0.379 ± 0.013 | 0.133, 0.140, 0.081 |

- **Width of the test CIs.** The test set has 35 SCZ / 54 CTL, 23 BD / 41 CTL, 22 MDD / 41 CTL and
  21 / 23 / 22 donors for Task B. Every test macro-F1 CI is about ±0.09–0.13 wide. The CIs of V6b, V6a
  and the study-only control overlap almost entirely on every task, so **the test set cannot rank these
  three**. It can say that SCZ and BD performance is above chance (the lower bound exceeds uniform
  chance), and that MDD and Task B are not clearly above chance.
- **Task B on test is worse than in CV** (0.37–0.40 vs 0.46). Task B's CV score is inflated by
  single-diagnosis batches (§5.3); without them CV is 0.508 for C0. Neither number survives to test.
- **Per-class results and confusion matrices** for test and CV (summed over 25 folds) are in
  `results/report_tables.md` and `results/final_comparison.json`.
  - SCZ test: V6b recalls 19/35 SCZ and 41/54 controls.
  - Task B test: BD recall is 0.17–0.22, and the errors are spread over all three classes, with no
    structure.
- **Seed variance vs fold variance.** For the stochastic models: random forest has fold SD 0.056–0.074 and
  seed SD 0.018–0.050. The random-graph and random-gene-set controls have seed SD 0.006–0.050. Logistic regression is
  deterministic (seed SD 0). Fold variance dominates everywhere.
- **Selection optimism.** The bias-corrected CV estimate of the selected configuration is 0.571 (mean over
  tasks), against a raw 0.572 (`results/final_selection.json`). Choosing among nearly equal candidates cost
  almost nothing.

## 4. Did graph structure help, and does that depend on signal strength? (`results/summary/V1_*, V2_*, V3_*`)

Development 5×5 CV; Δ macro-F1 vs the no-graph baseline C0; Nadeau–Bengio corrected p on the 25 fold
differences. The test assumes the correlation between folds is about n_test/n_train.

| Task (signal) | STRING Δ (p) | vs random graph | Co-expression Δ (p) | vs random | Reactome Δ (p) | vs random sets |
|---|---|---|---|---|---|---|
| SCZ (strong) | −0.000 (0.99) | −0.003 (0.66) | −0.004 (0.64) | −0.008 (0.35) | −0.052 (0.16) | −0.049 (0.11) |
| BD | −0.004 (0.87) | +0.000 (1.00) | −0.007 (0.76) | −0.000 (0.99) | −0.023 (0.50) | +0.015 (0.65) |
| MDD | −0.009 (0.52) | −0.004 (0.83) | −0.023 (0.27) | −0.015 (0.53) | −0.035 (0.32) | −0.007 (0.83) |
| Task B (weakest) | −0.003 (0.86) | −0.004 (0.84) | −0.009 (0.44) | −0.009 (0.62) | −0.068 (0.07) | −0.031 (0.35) |

No comparison survives correction. The graph contribution is ≈ 0 for the strong contrast and ≈ 0 for
the weak ones. If anything, pathway aggregation loses most where the signal is weakest. That is the
opposite of the "graphs help weak signals" hypothesis, but it is not significant. This agrees with the
skeptical benchmarks (Brouard 2024; Staiger 2012) and not with the uncontrolled positive reports
(Chereda 2019; Ramirez 2020).

## 5. Controls, and what they show (`results/summary/C*`)

1. **Study-only baseline (C1).** In CV it equals or beats the expression baseline on BD (0.596 vs 0.547)
   and MDD (0.595 vs 0.574), and is close on SCZ (0.626 vs 0.649; corrected p 0.55). The cohorts differ in
   case/control proportions. **By the protocol's criterion, no expression model clearly beats study
   composition in CV or on the test set.** Only the leave-one-cohort-out column separates them.
2. **Study-identity probe (C2).** Cohort is perfectly predictable from raw features (macro-F1 1.00). After
   in-fold per-batch standardisation it falls to 0.10 (chance 0.19); after ComBat it is 0.21.
3. **Batch correction can leak through composition (V5, C8, C9).** This is the main methodological finding.
   - Add-on ComBat with diagnosis as a protected covariate, fitted strictly in-fold, *raised* CV macro-F1
     by up to +0.11.
   - Under a **within-batch label permutation** it still scored above chance (+0.05 SCZ, +0.06 BD, +0.08
     Task B), while the baseline sat at chance.
   - ComBat shifts each held-out donor by its batch's training class composition, which transfers the
     cohort prior into the features without ever seeing a held-out label.
   - Batches containing only one diagnosis (the Sibille series in the SCZ, BD and Task B tasks) make it
     worst. The MDD task has none, and showed no gain.
   - The standard "fit on all data vs in-fold" demonstration **cannot** detect this, because both arms share
     it (the all-data arm scored *lower*).
   - V5 was disqualified by a rule fixed before the permutation result was seen.
4. **Leakage demos.**
   - Sample-level vs donor-grouped splitting on 530 arrays from 262 donors: 0.711 vs 0.702 macro-F1. The
     inflation is small here because a donor's extra arrays come from other regions or platforms, which
     are standardised as separate batches.
   - ComBat fitted on all data vs in-fold: see item 3.
5. **Sex-mismatch sensitivity (C6).** Removing the affected donors changes CV macro-F1 by +0.013 (SCZ),
   +0.001 (BD), +0.044 (MDD) and −0.013 (Task B). Matched random removals give +0.002, −0.009, +0.023 and
   +0.004, and the observed change is matched or exceeded by 15%, 25%, 5% and 85% of random draws. Only
   MDD is borderline. Removing just 2 donors moves CV by ±0.02 through fold reshuffling alone.
6. **Single-study dependence (C7).** Stanley supplies 52% of dev BD donors. Without Stanley, BD vs CTL is
   at chance for C0 (0.497), V6a (0.501) and V6b (0.502).
7. **Cell-type composition (V9).** Ten indices give SCZ 0.630 vs C0 0.649, recovering 88% of the excess
   over chance. For BD and MDD they recover about 19%, and for Task B 51%.

## 6. Variants tried (budget: 8 of 12 used; `results/variants.json`)

| # | Variant | Mean CV (4 tasks) | Verdict |
|---|---|---|---|
| 1 | V5 add-on ComBat | 0.619 | **invalid**: composition leakage |
| 2 | V1 STRING smoothing | 0.550 | failed |
| 3 | V2 co-expression kNN graph | 0.543 | failed |
| 4 | V3 Reactome pathway scores | 0.509 | failed (random sets better) |
| 5 | V6 ANOVA top-1,000 filter (V6a; V6b = graph twin) | 0.572 / 0.569 | no gain beyond noise; selected as final configurations |
| 6 | V7 random forest | 0.556 | failed |
| 7 | V8 ensemble C0+V1 | 0.551 | failed |
| 8 | V9 cell-type composition | 0.506 | failed as a model; informative |
| – | V4 GCN | – | not run (trigger not met) |
| ref | C0 baseline / C1 study-only | 0.554 / 0.549 | controls |

## 7. Limitations and deviations (all logged with reasons in LOG.md)

- **Task definitions admitted single-class batches** (my error). It is reported, not silently fixed: the
  primary tasks follow the approved plan, and the complete-batch versions (C9) are reported as a
  sensitivity analysis.
- **Duplicate detection is incomplete** (validation sensitivity 0.42). The mitigation is same-bank cohort
  merging plus liberal split groups. Cross-cohort leakage is not possible, but within-cohort CV may be
  slightly optimistic.
- **Researcher degrees of freedom:**
  - GSE92538 was added before any expression data was seen, but knowing it would make BD eligible.
  - The cohort-merge rule was amended after measuring detection sensitivity.
  - The C8 null was changed from within-cohort to within-batch permutation **before** any C8 result was
    read, because the first version would have hidden the leak.
- **Effects that cannot be separated:** medication, illness stage, agonal factors and cell composition.
  Isoform-level effects, where the largest disorder effects lie, are invisible on these gene-level arrays.
- **GSE21138 provenance is unverified.** Its "Victoria" label comes from my background knowledge.
- **The test set is small** (63–89 donors per task). It can confirm above-chance SCZ and BD performance,
  but it cannot rank models whose CV scores differ by < 0.05.

## 8. Reproducibility

Every run JSON records the Python, scikit-learn, numpy, pandas, scipy, R, limma and frma versions,
the platform, and the fold sizes (`results/versions.json`). The scripts are numbered by stage in
`scripts/`. Raw data is regenerated by `01`–`10`; modeling uses `pipeline.py` and the `run_*.py` scripts;
the single test evaluation is `final_test.py`. Library versions affect fold splits and gene selection,
so results are reproducible only with the recorded versions.
