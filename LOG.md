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
