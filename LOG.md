# LOG

Chronological record of decisions, pre-registrations and results. Newest entries at the bottom.

## 2026-09-30 — Environment

- Python 3.11.15, scikit-learn 1.9.1, numpy 2.4.6, pandas 3.0.6, scipy 1.17.1 (pip).
- R was NOT installed in the container at session start. Installed from Ubuntu apt:
  R 4.3.3, limma 3.58.1, affy 1.80.0, GEOquery 2.70.0, Biobase 2.62.0 (r-bioc-* packages).
  Verified with `packageVersion()` after install.
- 4 CPUs, 15 GB RAM, no GPU.

## 2026-09-30 — Pre-declared inclusion rules (written BEFORE any expression data was downloaded)

At the time of writing, only GEO metadata (series/sample annotation text, no expression values)
had been read.

1. **Minimum donors per contributing study:** a study contributes a disorder only if it
   provides at least **8 unique donors** with that diagnosis, counted after cross-study donor
   deduplication.
2. **Independence of studies.** The unit of "study" for the class-inclusion rule, for
   leave-one-study-out and for test-set stratification is the **cohort**, not the GEO accession:
   - GEO series that are artificial splits of one experiment (same lab, same brain bank, same
     design, split by sex or region) are one cohort.
   - Two series are merged into one cohort when at least 50% of the donors of the smaller one
     are found (metadata or expression identity) in the other. Otherwise a shared donor is
     assigned to one cohort only (the one where it has the most samples; ties → earlier GEO
     accession) and counts toward that cohort alone.
   - A disorder enters Task A/B only if at least **3 distinct cohorts** each contribute ≥ 8
     unique donors of it after this assignment.
3. Only psychiatric diagnoses (DSM-defined: SCZ, BD, MDD, ASD, alcohol use disorder, etc.)
   and unaffected controls are admissible. Non-brain tissue is excluded outright.

## 2026-09-30 — Metadata screening of the candidate list and additions (before any expression data)

- GSE11223 is "Colon biopsies from UC patients and healthy controls" (GPL1708, 202 samples) per GEO:
  not brain, not psychiatric. **Excluded.** It appears in the Gandal et al. (2018) accession list,
  so the candidate list contains at least one wrong accession.
- GSE28475 is a methods study (frozen vs fixed, DASL vs IVT, reference RNA dilutions, technical
  replicates) with ASD/control donors and no donor IDs linking to GSE28521. ASD comes from at most
  2 series (GSE28521, GSE28475), both Illumina, so ASD cannot reach 3 cohorts. Both are **excluded from
  modeling** (audit only).
- GSE29555 (alcohol use disorder, 4 amygdala/cortex regions, Illumina): the only AUD series, so
  **excluded from modeling**.
- GSE54567/54568/54571/54572 are one Sibille-lab (Pittsburgh) experiment split by sex (M/F) and region
  (BA9/BA25), so they form one cohort.
- GSE35978 is a SuperSeries (GSE35974 cerebellum + GSE35977 parietal). GEO text: "Stanley Medical
  Research Institute's Neuropathology Consortium and Array collections". GSE12649 GEO text: "102
  postmortem brains obtained from the Stanley Medical Research Institute". Expected heavy overlap.
- Additional series screened (metadata only):
  - **Added: GSE92538** (Pritzker; DLPFC; U133A + U133 Plus 2 CEL; SCZ/BD/MDD/CTL; subject IDs in GEO).
    It is an independent brain bank (UC Davis / UC Irvine / U Michigan). **Honesty note:** I noticed
    that without it BD would have only 2 candidate cohorts (Stanley, Pittsburgh). Its inclusion is a
    researcher degree of freedom exercised before any expression data was seen. It adds a third
    independent cohort for SCZ, BD and MDD alike.
  - Not added: GSE5388 (BD, U133A). Its characteristics fields (side of brain, fluphenazine
    equivalents, lithium, valproate, suicide, drug/alcohol ratings) are Stanley Array Collection
    fields, and 10 exact sex/age/PMI/pH matches were found with GSE35978 (chance ≈ 1.9). It adds no
    new cohort.
  - Not added: GSE21935 (SCZ BA22). It shares donor-ID strings with GSE17612 (36 shared IDs; sex and
    age agree for 35), so it is the same Charing Cross cohort in a second region. No new cohort.
  - Not added: GSE12654 (Stanley Consortium, U95Av2): Stanley again, and an older platform that would
    shrink the gene intersection. GSE62191 (Agilent two-colour ratios vs a cell-line pool; brain bank
    not stated in GEO; 30/29/29 CTL/BD/SCZ, suspiciously close to Stanley Array proportions; not
    verifiable). GSE87610 and GSE12679 (laser-capture cell-type profiles, not bulk tissue).
    GSE45642 (Pritzker controls only; GEO states they overlap GSE92538). GSE54562–54566 (further
    Sibille MDD cohorts; same lab and bank, no new cohort, three on Illumina without CELs).

## 2026-09-30 — Data starting point and preprocessing (decided before modeling)

- **Raw CEL files for every retained series** (all retained series are Affymetrix). Normalisation is
  **frozen RMA (fRMA; frma 1.54.0)** with the platform's frozen vectors (hgu133afrmavecs 1.5.0,
  hgu133plus2frmavecs 1.5.0, hugene.1.0.st.v1frmavecs 1.1.0; oligo 1.66.0 for HuGene). fRMA
  normalises and summarises each array against externally estimated frozen parameters, so no
  information passes between arrays, folds or the test set. Standard RMA would quantile-normalise
  test arrays together with training arrays, which rule 4 forbids.
- The deposited series matrices were used **only** for a first QC pass and were then abandoned:
  - The GSE12649 deposit is per-gene median-normalised (every gene's median is exactly 1.0).
  - The Sibille deposits give meaningless sex signals (Y/XIST not bimodal), whereas fRMA on their CELs
    gives clean single-sex series (0 mismatches in 108 arrays). The deposited values appear to be
    residualised or pair-normalised.
- Probe → gene: Bioconductor annotation packages (hgu133a.db 3.13.0, hgu133plus2.db 3.13.0,
  hugene10sttranscriptcluster.db 8.8.0, org.Hs.eg.db 3.18.0). Probesets that map to >1 Entrez gene are
  dropped. Multiple probesets per gene are averaged (fixed rule, no data-driven probe choice). For
  GSE92538 U133 Plus 2 arrays only U133A probesets are used, so both of its chip types summarise the
  same probes.
- Genes on all modeling platforms (U133A ∩ U133 Plus 2 ∩ HuGene 1.0 ST): 11,972 Entrez genes (11,480
  autosomal).

## 2026-09-30 — Sex check (fRMA data; label-free)

- Toker et al. (2016) style: Y score = mean of RPS4Y1, DDX3Y, KDM5D, USP9Y, EIF1AY, UTY, ZFY, NLGN4Y;
  XIST separately. Two-means clustering per chip type, pooled across series. A per-series split fails
  for single-sex series: the Sibille M/F series and the mostly-male Pritzker set produced spurious calls
  in my first attempt. XIST has no core transcript cluster on HuGene 1.0 ST, so GSE35978 is called on
  Y genes only.
- Categories: M, F, MIXED (Y high **and** XIST high: a two-donor mixture or XXY), LOW_BOTH (both low:
  technical).

## 2026-09-30 — Donor identity audit (expression; label-free)

- Method: residualise each (series, chip, region) unit by z-scoring genes and removing the top 10 PCs.
  Select 300 autosomal identity genes whose residuals reproduce across known same-donor pairs. Match by
  mutual best correlation, scored as min(row z, column z).
- The first attempt without PC removal picked agonal/stress genes (APOLD1, ADM, DDIT4, ZFP36) and
  recovered only 5/123 validation pairs at z ≥ 4. Removing 10 PCs (chosen on calibration pairs, tied
  with 5; 10 kept) gave genotype-driven genes (GSTM1/3/5, GSTT1, RPS26, PEX6, HLA-DPB1, NQO2).
- Thresholds from **negative controls** (1,927 mutual-best hits between series from different brain
  banks, which cannot share donors): max z = 5.08, 99th percentile 3.81 → **match** z ≥ 5.2,
  **probable** z ≥ 3.9.
- Validation on held-out known pairs (GSE35978 cerebellum ↔ parietal, 123 pairs, a hard cross-tissue
  case): 52/58 calls correct at z ≥ 3.9 (sensitivity 0.42); 10/10 correct at z ≥ 5.2 (sensitivity 0.08).
  **Detection sensitivity is limited, so undetected duplicates must be assumed to exist.** This drives
  the amendment below.

## 2026-09-30 — AMENDMENT to my pre-declared cohort rule (before modeling; no labels or results used)

Original rule: merge two series into one cohort when ≥ 50% of the smaller series' donors are *found*
in the other. Problem: expression detection of shared donors has limited sensitivity (0.42 at the
probable tier, cross-region), so "found" understates true overlap. Pittsburgh: GSE53987 PFC/STR vs
the Sibille series shows ≈ 20–23 detected shared donors out of 58 Sibille donors (34–40%). The true
overlap is plausibly ≥ 50%, and leaving them as separate cohorts would leak undetected duplicates
across leave-one-study-out splits.
**Amended rule:** series from the same brain bank with *any* confirmed shared donors are one cohort
(Stanley = GSE12649 + GSE35978; Pittsburgh = GSE53987 + GSE54567/8/71/72). This does not change the
class list: the cohort counts per disorder are unchanged or lower. It only makes splits more
conservative.
For splitting (CV folds, test set, LOSO), donors are grouped **liberally**: any same-bank mutual-best
hit with z ≥ 3.0 whose recorded diagnosis and expression sex agree joins the same split group. A
false grouping only makes two different donors share a fold, which is harmless. A missed duplicate
leaks.
**Stanley:** 80 of 102 GSE12649 samples were matched to GSE35978 donors (41 at match tier). The
remaining 22 are most likely undetected duplicates, since GEO says all 102 are Stanley brains. GSE12649 is therefore **used
for the audit only, not for modeling**. Stanley donors are modelled from GSE35978 parietal cortex, one
platform and one region.

## 2026-09-30 — Donor table, locked class list, locked test set, runtime (no modeling)

- `scripts/14_build_audit.py` → `DATA_AUDIT.json`, `data/donors.tsv`.
  - 569 donors across the modeling series; 129 appear in more than one series.
  - 7 eligible donors are "affected" (sex mismatch / both sexes / MIXED / identity-matched to a different
    diagnosis). 1 Pritzker donor has an unresolvable label (two subject IDs, one person, CTL vs BD) and is
    excluded.
  - Label rule: a donor's label is the recorded diagnosis of its modeling sample. A different label on
    another sample of the same donor makes it "affected".
- **LOCKED CLASS LIST: SCZ, BD, MDD** (5 / 3 / 3 cohorts with ≥ 8 donors; unchanged after removing
  affected donors). ASD and AUD fail (≤ 2 series, 1 lab each).
- Decisions written before modeling: primary analysis excludes affected donors. Task A uses only the
  cohorts contributing the disorder. Task B uses the cohorts contributing ≥ 2 disorders (Stanley,
  Pittsburgh, Pritzker).
- `scripts/15_lock_test_split.py`: StratifiedGroupKFold(4, shuffle, seed 20260930), first fold = test.
  135 test / 405 dev donors. scikit-learn 1.9.1. Saved to `data/split_locked.tsv`,
  `results/test_split.json`.
- `scripts/16_features.py`: 11,955 genes (chrY and XIST/TSIX removed). Dev and test matrices are in
  separate files.
- `scripts/17_string_graph.py`: STRING v12 ≥ 700 → 153,274 edges; 1,463 isolated genes.
- `scripts/18_timing.py`: **permuted labels, dev only, no metric computed.** LR + inner CV takes 62 s
  per 25 folds on the largest task; with graph smoothing 62.5 s; one GCN fold (150 epochs) 89 s.
- PLAN.md written. **Waiting for approval before any training or evaluation.**

## 2026-09-30 — Literature review expanded at the user's request (≥ 20 papers)

- 19 further searches (PubMed E-utilities, arXiv, Europe PMC). The queries and returned PMIDs are logged in
  `results/literature_queries.json`, and the abstracts are saved in `data/lit/`.
- Two searches returned the wrong paper (Iwamoto → Munakata 2005; Varma & Simon → Tsamardinos 2018). Both
  were logged and retried by PMID.
- Total now **39 sources read** (15 with full text or quoted passages, 24 abstracts). Pro-graph: 5 (+1
  design); skeptical: 6.
- **Plan changes (PLAN.md §12):** BBC-CV bias correction for the selected configurations; reserve V9
  cell-type-composition classifier; stated limitations (medication, illness stage, isoforms); the
  GSE21138 provenance caveat (the cohort key "Victoria" is unverified); three paper-vs-GEO conflicts
  recorded.
- The class list, test split and audit are **unchanged**; nothing in the literature alters them.
- Still no modeling. Awaiting approval.

## 2026-09-30 — PLAN APPROVED by the user. Implementation notes (no results yet)

- `scripts/pipeline.py` implements PLAN §5–7: identical folds for every configuration, in-fold batch
  adjustment, in-fold transforms, and one JSON per (configuration, task, scheme, fold, seed) under
  `results/runs/`, with skip-if-exists.
- **Computational change, not a modeling change:** L2 logistic regression is fitted in the n-dimensional
  row space of the training matrix (X·Vᵀ from its SVD). By the representer theorem this is exactly the same
  model. Checked on one dev fold: max |Δ probability| = 3e-16 at C = 1e-4 and 7e-6 at C = 1 (optimiser
  tolerance), for a 40× speed-up (43 s → 1 s per fold). No metric was computed in this check.
- ComBat: neuroCombat 0.2.12 fit on training donors (diagnosis as protected covariate). The apply-to-new
  step is re-implemented (same formula as neuroCombatFromTraining, which is marked "under development"):
  test donors get the training batch's γ*/δ* and the training-average covariate term, so no test labels
  are used.
- Graph smoothing re-standardises the smoothed features on the training part. This is a no-op for the
  no-graph baseline, and it removes the scale shrinkage of isolated genes (1,463 isolated in STRING).

## 2026-09-30 — Controls C0 (no-structure baseline) and C1 (study-only): expectations stated before running

- C0: per-batch standardisation → L2-LR (inner-CV C). Expectation: SCZ vs CTL modestly above chance
  (macro-F1 roughly 0.55–0.65, from small effect sizes, Mistry 2013: ~15% changes). BD/MDD vs CTL and
  Task B near chance.
- C1 study-only (one-hot cohort → balanced LR). Expectation: close to chance in CV, because every task
  cohort contributes both/all classes and folds are stratified by cohort. The cohort-conditional class
  priors differ across cohorts, so it may beat uniform chance somewhat (for example Pritzker is 75%
  CTL in the BD task). LOSO: an unseen cohort gets the global prior (≈ majority/balanced prediction).

## 2026-09-30 — RESULT: C0 and C1 (controls)

`results/summary/C0_baseline_LR.json`, `results/summary/C1_study_only.json`. CV = 5×5 donor-grouped
folds; values are mean ± fold SD; LR is deterministic, so seed SD = 0.

| Task | C0 macro-F1 | C1 study-only | Δ (C0−C1) | NB p (corr) | p (uncorr) | chance (uniform) |
|---|---|---|---|---|---|---|
| SCZ vs CTL | 0.649 ± 0.064 | 0.626 ± 0.013 | +0.023 | 0.548 | 0.114 | 0.494 |
| BD vs CTL | 0.547 ± 0.068 | 0.596 ± 0.030 | −0.049 | 0.224 | 0.003 | 0.486 |
| MDD vs CTL | 0.574 ± 0.080 | 0.595 ± 0.018 | −0.021 | 0.640 | 0.215 | 0.490 |
| Task B | 0.445 ± 0.065 | 0.379 ± 0.013 | +0.066 | 0.078 | <0.001 | 0.333 |

AUC Δ (C0−C1): SCZ +0.068 (corr p = 0.053), BD −0.057, MDD −0.010, Task B +0.011.

**Interpretation.** The study-only model exploits differing class proportions across cohorts (e.g.
Stanley ≈ 50% SCZ, Pritzker ≈ 22% SCZ). In CV, C0 does **not** significantly beat it on any task. This is
the protocol's failure condition for "learned something beyond study composition". Leave-one-cohort-out
separates them: C1 collapses on unseen cohorts (macro-F1 0.08–0.44), while C0 keeps SCZ signal in 3/5
unseen cohorts (Pittsburgh 0.79, CharingCross 0.66, Victoria 0.66; Pritzker 0.52, Stanley 0.54). C0's
per-batch standardisation deliberately removes cohort means, so it *cannot* use cohort priors. Its CV
score is expression signal, but not more than the priors alone give.

## 2026-09-30 — PRE-REGISTRATION V5 (variant 1/12): add-on ComBat instead of per-batch standardisation

- Change: in-fold neuroCombat (parametric EB; diagnosis as protected covariate), applied to C0 now and to
  the best graph model later (same variant, equal treatment).
- Mechanism: EB shrinkage of batch location/scale (Johnson 2007) stabilises small batches.
- Expectation: Δ macro-F1 vs C0 within ±0.02. Nygaard (2016) predicts that protected-covariate ComBat
  over-separates classes *in-sample*. In-fold fitting should not carry that into held-out donors, but
  preserved diagnosis effects in training could make training separation look larger than it transfers.
- Failure: Δ ≤ 0 vs C0 (Nadeau–Bengio), or no gain beyond noise (|Δ| < 0.03).
- The same runs double as the in-fold arm of leakage demo C4(ii); the all-data arm (ComBat fitted with
  held-out labels) is the control.

## 2026-10-01 — PRE-REGISTRATIONS V1, V2, V3 (variants 2, 3, 4 of 12), written before any graph run

All three: same folds, per-batch standardisation and inner-CV L2-LR as C0. The no-graph ablation **is**
C0. Each has a density-matched random control rebuilt per fold with 3 seeds (control, not counted).
Comparisons use Nadeau–Bengio on the 25 fold-level differences, per task.

**V1 STRING-smoothed features (curated graph).** X' = ½X + ½·D^-½AD^-½X over STRING v12 ≥ 700 (153,274
edges), then training-fit re-standardisation and LR (SGC form, Wu 2019).
- Mechanism: neighbourhood averaging reduces per-gene noise if disease effects are coordinated within
  interaction modules (Chuang 2007).
- Expectation: Δ vs C0 within ±0.02 on SCZ; if any gain, larger on BD/MDD/Task B (diffuse signal). Not
  better than the random graph (Staiger 2012; Brouard 2024).
- Failure: Δ vs C0 ≤ 0, or Δ vs random graph not significant (corr p ≥ 0.05), on a task.

**V2 co-expression kNN graph (expression-derived).** k = 10 neighbours by |Pearson r| among training
donors after batch standardisation, symmetrised; same smoothing. The random control matches V2's own
edge count per fold.
- Mechanism: data-driven modules (Chen 2013 found SCZ/BD signal in co-expression modules).
- Expectation: ≈ C0. Its graph partly encodes cell-type composition axes (Hagenauer 2018), which may
  help BD/SCZ slightly.
- Failure: as V1.

**V3 Reactome pathway scores (curated gene sets).** Mean of standardised genes per Reactome pathway
(10–300 member genes in our feature set; 1,098 sets) replaces genes, then LR. Control: random gene sets
with identical sizes, rebuilt per fold (3 seeds).
- Mechanism: pathway-activity aggregation (Lee-type, evaluated in Staiger 2012).
- Expectation: ≤ C0 on SCZ (information loss); possible small gain on weak tasks.
- Failure: Δ vs C0 ≤ 0, or not better than random gene sets (corr p ≥ 0.05).

**Conditional V4 (GCN)** runs only if V1 or V3 beats its random control with corrected p < 0.10 on any task
(PLAN §7.2).

## 2026-10-01 — RESULTS: V5 (ComBat), C4(ii) ComBat leakage demo, C7, C2

**V5 add-on ComBat vs C0** (CV macro-F1; Δ; Nadeau–Bengio corrected p / uncorrected p):
- SCZ 0.720 ± 0.064 vs 0.649: Δ +0.071, p 0.001 / <1e-4
- BD 0.623 vs 0.547: Δ +0.076, p 0.021
- MDD 0.577 vs 0.574: Δ +0.003, p 0.92
- Task B 0.554 vs 0.445: Δ +0.109, p 0.003
- LOSO: ComBat cannot adjust an unseen cohort (it falls back to self-standardisation), so LOSO is ≈ C0's
  (SCZ: Pittsburgh 0.79, Victoria 0.63, CharingCross 0.55, Pritzker 0.52, Stanley 0.54).
- **Against my pre-registered expectation (±0.02).**

**C4(ii) ComBat fitted on ALL dev data with held-out labels** (the "leaky" arm) vs in-fold V5: SCZ 0.667
(Δ −0.054, p 0.053), BD 0.554 (−0.069, p 0.017), MDD 0.522 (−0.055, p 0.057), Task B 0.485 (−0.069, p 0.045).
**The leaky arm is LOWER than in-fold, the opposite of what a leakage demo should show.** Not
understood yet; see the diagnostic below.

**C7 single-study dependence.** Stanley supplies 51.6% of dev BD donors. Without Stanley (C0):
- BD vs CTL 0.497 ± 0.083 (chance 0.486), versus 0.547 with Stanley. **The BD signal depends on the Stanley
  cohort.**
- Task B without Stanley: 0.511 ± 0.096.

**C2 study-identity probe** (multinomial LR, 5×5 grouped CV, primary dev donors):
- Cohort (5 classes): raw 1.000; after per-batch standardisation 0.101 (uniform chance 0.189, majority
  0.095); after ComBat 0.212.
- Batch (8 classes): raw 0.873; after per-batch standardisation 0.053 (chance 0.113); after ComBat 0.109.
- Per-batch standardisation removes study identity entirely. ComBat leaves a little cohort information
  (0.21 vs 0.19 chance).

## 2026-10-01 — ADDED CONTROL C8 (not a variant): within-cohort label permutation

Why: V5's large gain and the inverted leakage demo need a check that V5 is not exploiting cohort
composition. ComBat preserves diagnosis effects when estimating batch effects, so a batch's residual
mean keeps its case fraction × disease effect, and C2 shows cohort is weakly identifiable after ComBat.

Test: permute diagnosis labels **within each cohort** (cohort × class counts preserved, expression–label
link destroyed), 5 permutations × 5×5 CV, for C0 and V5. The study-only model keeps its full CV score under
this permutation. A model that learns biology must drop to chance.

Prediction if V5 is clean: both ≈ chance (≈ 0.49 binary, 0.33 Task B). If V5 stays clearly above chance,
its CV gain is composition leakage and V5 is rejected as an evaluation artefact.

## 2026-10-01 — Diagnosis of V5 / inverted leakage demo; C8 REVISED before any C8 result was read

**One-fold diagnostic** (SCZ, folds r0f0–r0f4):
- The two Sibille BA9 batches (GSE54567, GSE54568) contain **only controls** in the SCZ task. They enter
  via the Pittsburgh cohort, which contributes SCZ through GSE53987.
- With diagnosis as a protected covariate, ComBat aligns such a batch to the control mean, so held-out
  donors of that batch get P(SCZ) ≈ 0.18–0.20 versus ≈ 0.40 for controls elsewhere.
- More generally, ComBat's add-on step shifts every held-out donor by its batch's *training* class
  composition × the estimated disease effect. This injects the cohort/batch prior into the features.
  That is study-composition leakage: no held-out label is used, but batch priors learned from training
  labels transfer.
- The "all-data" leakage arm does the same, which is why it could not reveal the problem. Its lower
  score is consistent with this (its batch estimates net out true test labels, so it injects less
  composition shift into test donors).

**Batches lacking a class, per task** (`results/task_batch_composition.json`):
- SCZ and BD tasks: GSE54567 (6 CTL only) and GSE54568 (7 CTL only).
- Task B: GSE54567 (10 MDD only) and GSE54568 (6 MDD only).
- MDD task: none.
- This matches the V5 pattern: big gains on SCZ, BD and Task B; none on MDD (+0.003).

**My error (escalation, rule 13):** my task definition ("all CTL donors of cohorts contributing D")
lets single-class *batches* into Tasks A-SCZ, A-BD and B, even though every cohort contributes both
classes. The plan is kept as approved, and the issue is reported, not silently fixed.

**C8 revised:** the first C8 launch (within-*cohort* permutation) was stopped and deleted after 0
results were read. It would have mixed labels into the single-class batches and so destroyed exactly the
shortcut under test, a false "clean" verdict. **C8 now permutes labels within BATCH**: batch × class
counts are preserved, and the expression–label link is destroyed.
- Prediction if V5 exploits composition: V5 stays near the study-only level (≈ 0.6 binary) under the
  permutation, while C0 ≈ chance.
- **Decision rule (fixed now):** if V5's permuted mean macro-F1 exceeds the uniform-chance level by
  ≥ 0.03 on a task, V5's CV gain on that task is declared composition leakage, and V5 is not eligible as
  a final model.

**Added control C9:** C0 and V5 re-run on each task restricted to batches containing every class of the
task. This is a sensitivity analysis, not a change to the primary tasks.

## 2026-10-01 — RESULT C8 (within-batch label permutation null; 5 permutations × 5×5 CV)

| Task | uniform chance | C0 permuted | V5 permuted | V5 − chance | C0 real | V5 real |
|---|---|---|---|---|---|---|
| SCZ vs CTL | 0.494 | 0.517 | 0.546 | +0.052 | 0.649 | 0.720 |
| BD vs CTL | 0.486 | 0.497 | 0.547 | +0.061 | 0.547 | 0.623 |
| MDD vs CTL | 0.490 | 0.500 | 0.514 | +0.024 | 0.574 | 0.577 |
| Task B | 0.333 | 0.348 | 0.409 | +0.076 | 0.445 | 0.554 |

- C0 is at chance under the null (largest excess +0.023 on SCZ). Its pipeline does not exploit batch
  composition.
- **V5 stays above chance by ≥ 0.03 on SCZ, BD and Task B → by the pre-fixed rule, V5's CV gain is
  declared composition leakage; V5 is NOT eligible as a final model** (variant 1/12 counted, outcome:
  failed / invalid).
- Nuance (reported, not used to rescue V5): on real labels V5 exceeds its own null by more than C0 exceeds
  its own (SCZ +0.174 vs +0.132; Task B +0.145 vs +0.097). Part of ComBat's gain may be the legitimate
  mechanism: it does not remove disease signal from batches with unequal case fractions, which per-batch
  standardisation does. The two parts cannot be separated with these data.
- Methodological finding: **add-on ComBat with diagnosis as a protected covariate transfers batch class
  composition from training labels into held-out donors.** The protocol's "in-fold" requirement is met
  and no held-out label is used, yet the CV estimate is still inflated. The usual leakage demo
  (fit on all data vs in-fold) cannot detect this, because both arms share it.

## 2026-10-01 — RESULT C9 (tasks restricted to batches containing every class)

Dropped: 13 Sibille CTL-only donors (SCZ and BD tasks) and 16 Sibille MDD-only donors (Task B). MDD task
unchanged.

| Task | C0 primary | C0 complete-batches | V5 complete-batches | Δ V5−C0 (corr p) |
|---|---|---|---|---|
| SCZ | 0.649 | 0.659 ± 0.066 | 0.706 ± 0.071 | +0.047 (0.049) |
| BD | 0.547 | 0.540 ± 0.102 | 0.594 ± 0.095 | +0.054 (0.053) |
| MDD | 0.574 | 0.574 | 0.577 | +0.003 (0.92) |
| Task B | 0.445 | 0.508 ± 0.061 | 0.540 ± 0.061 | +0.032 (0.139) |

- Single-class batches **hurt C0** (per-batch centring of an MDD-only batch makes those donors look
  average; Task B C0 0.445 → 0.508 without them) and **inflate V5**.
- Removing them roughly halves V5's advantage. The remainder is still partly the mixed-batch composition
  shift shown by C8.
- **Escalation (rule 13):** the approved task definition admits single-class batches. I am following the
  approved plan (primary = as defined); C9 is reported alongside every primary number as the
  sensitivity analysis. In my judgement, the complete-batch tasks are the better-posed versions of Tasks
  A-SCZ, A-BD and B.

## 2026-10-01 — RESULTS V1, V2, V3 (graph variants; base = per-batch standardisation; no-graph = C0)

CV macro-F1 (mean ± fold SD; seed SD for random controls):

| Task | C0 no-graph | V1 STRING | V1 random (seed SD) | V2 coexpr | V2 random | V3 Reactome | V3 random sets |
|---|---|---|---|---|---|---|---|
| SCZ | 0.649 | 0.649 ± 0.067 | 0.652 (0.006) | 0.645 | 0.653 (0.009) | 0.597 | 0.645 (0.026) |
| BD | 0.547 | 0.543 | 0.543 (0.018) | 0.540 | 0.541 (0.016) | 0.524 | 0.509 (0.050) |
| MDD | 0.574 | 0.565 | 0.569 (0.024) | 0.551 | 0.566 (0.029) | 0.539 | 0.547 (0.041) |
| Task B | 0.445 | 0.442 | 0.446 (0.024) | 0.436 | 0.445 (0.026) | 0.377 | 0.408 (0.049) |

Nadeau–Bengio (25 fold differences; corrected p / uncorrected p):
- V1 vs no-graph: Δ −0.000 / −0.004 / −0.009 / −0.003 (SCZ/BD/MDD/B), all p_corr ≥ 0.52. V1 vs random:
  |Δ| ≤ 0.004, all p_corr ≥ 0.66.
- V2 vs no-graph: Δ −0.004 / −0.007 / −0.023 / −0.009, p_corr ≥ 0.27. V2 vs random: Δ −0.008 / −0.000 /
  −0.015 / −0.009.
- V3 vs no-graph: Δ −0.052 (p 0.16; unc 0.001) / −0.023 / −0.035 / −0.068 (p 0.070; unc <0.001). V3 vs
  random sets: −0.049 / +0.015 / −0.007 / −0.031, all p_corr ≥ 0.11.

**Verdicts (pre-registered failure criteria):**
- **V1 FAILED, V2 FAILED, V3 FAILED** on every task.
- No graph beats its density-matched random control; random gene sets beat Reactome on 3/4 tasks.
- **V4 (GCN) is NOT triggered** (condition: V1 or V3 beats its random control with corr p < 0.10 on any
  task; none did).
- Central question: graph − no-graph is ≈ 0 on the strong contrast (SCZ: −0.000) and ≈ 0 on the weak ones
  (Task B: −0.003 STRING). No difference by signal strength for STRING or co-expression. Pathway
  aggregation loses more on the weakest task (Task B −0.068) than on SCZ (−0.052). That is the opposite of
  the "graphs help weak signal" hypothesis, but it is not significant under the corrected test.
- Running variant count: **4 / 12** (V5, V1, V2, V3).

## 2026-10-01 — PRE-REGISTRATIONS V6, V7, V8 (variants 5, 6, 7 of 12)

- **V6 in-fold ANOVA-F top-1,000 filter** before LR. Applied to C0 (V6a) and to V1 (V6b; equal treatment
  of the best graph model by mean CV, V1 = 0.550 vs V2 0.543, V3 0.509).
  - Mechanism: fewer noise dimensions for n ≪ p.
  - Expectation: ≈ C0; possibly worse, because a univariate filter is unstable at n ≈ 200 (Varma & Simon
    2006 warn about selection-induced variance).
  - Failure: Δ vs C0 ≤ 0.
- **V7 random forest** (500 trees, √p features, balanced, 3 seeds) on per-batch-standardised genes.
  - Expectation: below C0 (Hornung 2017; Kapoor & Narayanan 2023: LR is hard to beat).
  - Failure: Δ vs C0 ≤ 0.
- **V8 ensemble**: mean probability of C0 and V1, computed from the stored out-of-fold probabilities on
  identical folds (no refit).
  - Expectation: ≈ C0 (C0 and V1 are nearly identical models).
  - Failure: Δ vs the better member ≤ 0 or < 0.03.

## 2026-10-01 — PRE-REGISTRATION V9 (variant 8 of 12): cell-type composition classifier

- Features: 10 BrainInABlender cell-type indices (Hagenauer et al. 2018 marker database, repo commit
  015cc35). Each is the mean of in-fold per-batch-standardised marker genes; markers listed for > 1
  primary cell type are dropped. Set sizes 11–285. Then the same inner-CV LR. Control: random gene sets
  of identical sizes (3 seeds).
- Purpose and mechanism: Ramaker 2017 and Hagenauer 2018 report neuron↓/astrocyte↑ in SCZ/BD and
  composition as the main axis of variance. This measures how much of C0's above-chance signal a
  10-number composition summary carries.
- Expectation: SCZ and BD recover roughly half of C0's excess over chance; MDD and Task B ≈ chance.
- As a candidate final model it fails if Δ vs C0 ≤ 0 (expected). The informative quantity is
  (V9 − chance)/(C0 − chance) per task, reported regardless.

## 2026-10-01 — RESULTS V6, V7, V8

| Task | C0 | V6a filter+LR | V6b filter+STRING | V7 RF (seed SD) | V8 ensemble |
|---|---|---|---|---|---|
| SCZ | 0.649 | 0.668 ± 0.072 | 0.669 ± 0.065 | 0.680 ± 0.056 (0.018) | 0.650 |
| BD | 0.547 | 0.561 | 0.568 | 0.520 (0.050) | 0.539 |
| MDD | 0.574 | 0.592 | 0.576 | 0.543 (0.049) | 0.569 |
| Task B | 0.445 | 0.467 | 0.464 | 0.480 (0.038) | 0.445 |
| mean | 0.554 | **0.572** | 0.569 | 0.556 | 0.551 |

Nadeau–Bengio, corrected p (uncorrected):
- V6a vs C0: +0.019 (0.37; 0.021), +0.014 (0.60), +0.018 (0.67), +0.023 (0.47; 0.058).
- V6b vs V1: +0.020 (0.28), +0.025 (0.30), +0.011 (0.75), +0.022 (0.55).
- V7 vs C0: +0.032 (0.25; 0.004), −0.027 (0.58), −0.031 (0.41), +0.035 (0.26; 0.005).
- V8 vs C0: +0.001, −0.008, −0.006, −0.000.

Verdicts:
- **V6: no gain beyond noise** (all |Δ| < 0.03, none significant under correction). Not a strict
  pre-registered "failure" (Δ > 0), but it does not hold across folds.
- **V7: FAILED** (Δ ≤ 0 on BD and MDD; the SCZ and Task B gains are not significant).
- **V8: FAILED.**
- Graph inside the filtered model (V6b − V6a): −0.000 / +0.007 / −0.016 / −0.003, again ≈ 0.

Running variant count: **8 / 12** (V5, V1, V2, V3, V6, V7, V8; V9 pre-registered and queued).
(Note: an old shell "wait" helper was killed by its time limit; it ran no computation.)

## 2026-10-01 — RESULT V9 (cell-type composition) and FINAL SELECTION (frozen before the test set is opened)

**V9** (10 BrainInABlender indices; CV macro-F1; Δ vs C0 corr p; Δ vs random sets corr p; fraction of C0's excess
over chance recovered):
- SCZ 0.630 (−0.019, 0.60; +0.043, 0.14; **0.88**)
- BD 0.498 (−0.049; +0.050; 0.19)
- MDD 0.505 (−0.069, 0.085; +0.009; 0.18)
- Task B 0.390 (−0.055; +0.032; 0.51)
- **V9 FAILED** as a model.
- Finding: ten composition numbers carry about 88% of the detectable SCZ signal and about half of Task B's.
  This is consistent with neuron/glia shifts (Ramaker 2017; Hagenauer 2018), which could be disease biology,
  medication, agonal state or dissection; these data cannot separate them.

**Final variant count: 8 / 12** (V5, V1, V2, V3, V6, V7, V8, V9; V4 not triggered). Table:
`results/variants.json`.

**Selection (PLAN §9 rule; `results/final_selection.json`):**
- Eligible pool excludes V5 (composition leakage) and C1 (control).
- **Best model = V6b** (in-fold ANOVA top-1,000 → STRING smoothing α = 0.5 → inner-CV L2-LR). Mean CV
  0.569; the best non-baseline (V6b 0.569 > V8 0.550 ≈ V1 0.550 > V2 0.543 > V3 0.509 > V9 0.506).
- **Best simple baseline = V6a** (in-fold ANOVA top-1,000 → inner-CV L2-LR). Mean CV 0.572 (> V7 0.556 >
  C0 0.554).
- BBC-CV (bias-corrected estimate of the selected configuration; 5 repeats × 200 bootstraps): 0.571 mean
  (SCZ 0.679, BD 0.558, MDD 0.580, B 0.469). Selection optimism is negligible because candidates are nearly
  equal.

**Frozen test protocol (declared now, before any test feature is read):**
1. Per task: train on all primary dev donors of the task (affected excluded); test on the task's
   primary test donors. Per-batch standardisation uses training statistics of the same batch. The
   filter and inner-CV C (seed 0) are fitted on training only. One fit per configuration per task.
2. Configurations evaluated: **V6b (best model), V6a (best simple baseline)**, and **C1 study-only as a
   reference control** (untuned; needed to judge "beats study composition" on test).
3. Reported per task: macro-F1, MCC, AUC, per-class precision and recall, confusion matrix, with
   percentile bootstrap 95% CIs (2,000 resamples of test donors), test chance levels (uniform, prior,
   majority) and per-class test n.
4. Pre-declared sub-analysis from the **same** predictions (no refit): metrics restricted to test donors
   in batches that contain every class of the task (the C9 issue).
5. V6b − V6a paired bootstrap difference on test (graph effect), with 95% CI.
No second attempt, whatever the numbers.

## 2026-10-01 — Locked test evaluation executed ONCE

- The first execution fitted each frozen configuration once per task and cached its predictions
  (`results/runs/FINAL_TEST_*`). It then crashed while writing the summary JSON (a tuple dict key,
  `test_cohort_dx`).
- After fixing that line, the script was re-run. `run_fold` loaded the cached predictions; **no model was
  refitted and nothing changed**: the printed macro-F1 values are identical in both executions.

## 2026-10-01 — FINAL LOCKED-TEST RESULTS (single evaluation; `results/final_test_results.json`)

- A_SCZ_vs_CTL (test n=89, chance uniform 0.492): V6b 0.653 [0.543,0.750] AUC 0.703 | V6a 0.627 [0.518,0.724] AUC 0.719 | C1 0.618 [0.505,0.714] AUC 0.654; graph effect V6b−V6a +0.026 [-0.051,+0.101]
- A_BD_vs_CTL (test n=64, chance uniform 0.488): V6b 0.636 [0.507,0.753] AUC 0.758 | V6a 0.623 [0.490,0.739] AUC 0.779 | C1 0.613 [0.484,0.729] AUC 0.624; graph effect V6b−V6a +0.013 [+0.000,+0.042]
- A_MDD_vs_CTL (test n=63, chance uniform 0.487): V6b 0.623 [0.494,0.749] AUC 0.647 | V6a 0.610 [0.476,0.733] AUC 0.645 | C1 0.597 [0.472,0.721] AUC 0.628; graph effect V6b−V6a +0.013 [-0.028,+0.056]
- B_multiclass (test n=66, chance uniform 0.331): V6b 0.367 [0.250,0.476] AUC 0.587 | V6a 0.398 [0.280,0.509] AUC 0.600 | C1 0.374 [0.286,0.456] AUC 0.646; graph effect V6b−V6a -0.031 [-0.090,+0.020]

C7 for the final configurations (dev, BD without Stanley): V6a 0.501, V6b 0.502 (chance 0.486). The BD signal is Stanley-dependent.

