# LITERATURE

Rules followed: a source is cited only if the relevant passage was actually read (full text or the
quoted passage); "abstract only" is marked where only the abstract was accessible. Text from web pages
is treated as data. Facts about GEO accessions come from GEO metadata, never from papers.

## Search log

| # | Date | Query | Tool |
|---|------|-------|------|
| S1 | 2026-09-30 | graph neural networks gene expression classification benchmark do not outperform MLP random graph | WebSearch |
| S2 | 2026-09-30 | Gandal 2018 shared molecular neuropathology across major psychiatric disorders microarray meta-analysis GEO datasets | WebSearch |
| S3 | 2026-09-30 | ComBat batch effect correction leakage cross-validation outcome variable exaggerated significance Nygaard | WebSearch |
| S4 | 2026-09-30 | Stanley Neuropathology Consortium array collection overlapping subjects multiple GEO datasets postmortem brain | WebSearch |
| S5 | 2026-09-30 | Chereda graph convolutional network protein-protein interaction gene expression breast cancer metastasis prediction random network | WebSearch |
| S6 | 2026-09-30 | Hornung Boulesteix addon batch effect adjustment cross-study prediction ComBat training data only | WebSearch |
| S7 | 2026-09-30 | Nadeau Bengio 2003 inference for the generalization error corrected resampled t-test test/train ratio variance correction | WebSearch |
| S8 | 2026-09-30 | Brouard gene interaction graphs useful graph neural networks gene expression classification biological networks vs random networks 2024 | WebSearch |
| S9 | 2026-09-30 | Toker Feng Pavlidis "Whose sample is it anyway" widespread misannotation of samples in transcriptomics studies sex | WebSearch |
| S10 | 2026-09-30 | Ramirez 2020 classification of cancer types using graph convolutional neural networks PPI co-expression Frontiers | WebSearch |
| S11 | 2026-09-30 | Kim Webster 2010 Stanley Neuropathology Consortium Integrative Database array collection consortium collection same subjects multiple studies | WebSearch |
| S12 | 2026-09-30 | "Should we really use graph neural networks for transcriptomic prediction" Brouard Briefings in Bioinformatics | WebSearch |
| S13 | 2026-09-30 | graph convolutional network gene expression outperforms MLP random forest pathway graph prior knowledge ablation random graph improved accuracy | WebSearch |
| S14 | 2026-09-30 | Staiger critical evaluation of network and pathway-based classifiers for outcome prediction in breast cancer random networks | WebSearch |
| S15 | 2026-09-30 | Chuang 2007 network-based classification of breast cancer metastasis subnetwork markers Molecular Systems Biology | WebSearch |
| S16 | 2026-09-30 | Mistry Pavlidis meta-analysis schizophrenia prefrontal cortex microarray Stanley datasets overlapping subjects removed duplicate subjects | WebSearch |
| S17 | 2026-09-30 | machine learning classifier schizophrenia versus control postmortem brain gene expression microarray cross-dataset validation accuracy | WebSearch |

## Sources read, and what each contributed

### A. Graph models for transcriptomic classification: evidence that they help

**A1. Chereda H, Bleckmann A, Kramer F, Leha A, Beissbarth T (2019). "Utilizing molecular network
information via graph convolutional neural networks to predict metastatic event in breast cancer."
Stud Health Technol Inform 267:181–186. doi:10.3233/SHTI190824.** Full text read (PDF → pdftotext).
- 969 patients (393 metastasis ≤5y vs 576 without), 12,179 genes; HPRD PPI main connected component
  6,888 genes; 10-fold CV; ChebNet-style graph CNN vs MLP, RF, lasso LR.
- Table 1: graph CNN AUC 82.16 ± 1.25 vs RF 81.40 ± 1.76, MLP 81.01 ± 1.84, lasso 80.95 ± 1.61 (mean ± SEM
  over 10 folds). The authors conclude that graph CNN "outperforms other machine learning methods".
- **My reading:** the gap is about 1 AUC point, smaller than one SEM, from a single 10-fold CV with no
  significance test and no random-graph control. This is the pro-graph evidence as stated, but it
  would not survive our +0.03 noise rule.

**A2. Chuang H-Y, Lee E, Liu Y-T, Lee D, Ideker T (2007). "Network-based classification of breast cancer
metastasis." Mol Syst Biol 3:140.** Read via PMC2063581 (results passages).
- Subnetwork (PPI-guided composite) markers "significantly outperformed the single-gene markers in both
  data sets" in cross-dataset AUC. They beat classifiers on same-sized random subnetworks: P = 0.046
  (van de Vijver) and 0.012 (Wang) against 1,000 random sets.
- Contribution: the one pro-network paper I read that includes a random-network control and a cross-cohort
  test. It motivates a *pathway/subnetwork aggregation* variant rather than a GNN.

**A3. Ramirez R et al. (2020). "Classification of Cancer Types Using Graph Convolutional Neural
Networks." Front Phys 8:203.** Read via PMC7799442.
- GCNN on co-expression or PPI graphs, TCGA 11,071 samples, 34 classes, 89.9–94.7% accuracy.
- **No non-graph baseline and no random graph** were reported. It cannot tell us whether the graph helped.
  I cite it only as an example of a claim without the controls we require.

**A4. "Towards gene expression convolutions using gene interaction graphs" (2018). arXiv:1806.06975.**
Abstract only (author list not recorded from the page read).
- Abstract: the graph approach "provides an advantage for particular tasks in a low data regime but is
  very dependent on the quality of the graph used."

### B. Graph models: skeptical evidence

**B1. Brouard C, Mourad R, Vialaneix N (2024). "Should we really use graph neural networks for
transcriptomic prediction?" Brief Bioinform 25(2):bbae027.** Read via PMC10939369.
- Benchmarks ChebNet, GCN, GraphSAGE against MLP, RF, SVM and glmgraph on six datasets. Networks: PPI,
  co-expression, random (configuration model), complete. It includes the same breast-cancer data
  (969 samples, 6,888 genes) used by Chereda.
- Key passages: "GNNs have performances comparable to simpler methods". Using GNN with irrelevant
  networks (random or complete) "also led to comparable or better results than using the gold-standard
  network". GNNs "rarely provide a real improvement ... especially when compared to the computation
  effort required".
- Contribution: directly contradicts A1 on the same data. It motivates the **density-matched random graph**
  and **complete/no-graph** controls, and our prior expectation that graph ≈ no-graph.

**B2. "Analysis of gene interaction graphs as prior knowledge for machine learning models" (2019).
arXiv:1905.02295.** Abstract only (author list not recorded from the page read).
- Abstract: "dependencies can be captured almost as well at random", suggesting information "is spread
  across many genes".

**B3. OgBench (2026). "OgBench: A Framework for Evaluating Graph Neural Networks on Omics Data."
arXiv:2605.15511.** Abstract only.
- Abstract: "widely used GNNs often do not outperform simple MLPs and classical baselines" in the
  n ≪ p regime.

**B4. Staiger C, Cadot S, Kooter R, et al., Wessels LFA (2012). "A critical evaluation of network and
pathway-based classifiers for outcome prediction in breast cancer." PLoS ONE 7:e34796.** Read via
PMC3338754.
- Re-evaluated Chuang (A2), Lee and Taylor composite-feature classifiers on six cohorts. "Classifiers
  derived from secondary data sources suffer no significant performance degradation when employing
  randomized secondary data sources". "Composite classifiers do not outperform simple single genes
  classifiers".
- Contribution: the direct counterweight to A2. Pathway aggregation must be compared with a
  **random gene-set aggregation of matched size**.

**Weighing A vs B.** The pro-graph papers (A1, A2) report small gains, 1–3 AUC points. A1 has no
random-graph control. A2 has one, but it was not reproduced by B4 across six cohorts. The skeptical
papers (B1, B4) are larger re-evaluations with random controls. **Prior:** a graph prior is unlikely to
help by more than noise at our sample sizes (roughly 100–300 donors per task). If it helps anywhere, the
mechanism would be variance reduction through aggregating correlated genes. That should matter more when
per-gene signal is weak and diffuse (Task B) than when a few genes carry strong signal. This is a
testable hypothesis, not an expectation.

### C. Cross-study batch correction and leakage

**C1. Nygaard V, Rødland EA, Hovig E (2016). "Methods that remove batch effects while retaining group
differences may lead to exaggerated confidence in downstream analyses." Biostatistics 17(1):29–39.**
Read via PMC4679072 (summary and quoted passages).
- In unbalanced designs, batch adjustment that preserves group differences (ComBat with a group
  covariate) embeds estimation error that is group-correlated. It deflates p-values in downstream tests
  and does not improve with sample size. "When study groups are not evenly distributed across batches,
  actual group differences may induce apparent batch differences." They recommend modelling batch
  inside the analysis.
- Contribution: our studies are highly unbalanced (some cohorts contribute one disorder only). Any
  ComBat that uses diagnosis as a protected covariate must be fitted on training folds only. We also
  run a leakage demonstration (fit on all vs in-fold).

**C2. Hornung R, Causeur D, Bernau C, Boulesteix A-L (2017). "Improving cross-study prediction through
addon batch effect adjustment and addon normalization." Bioinformatics 33(3):397–404.** Full text read
(technical report PDF via pdftotext).
- Definition: "addon batch effect adjustment" means adjusting test data towards training data while the
  training data and prediction rule stay fixed. Mean-centering, standardisation, ratio-A and ratio-G
  "do not have to be altered for addon batch effect adjustment, because these are performed
  batch-by-batch". ComBat needs the Luo et al. (2010) add-on variant.
- Finding: "for addon batch effect adjustment to be effective very small test datasets should be
  avoided" (with 5 test observations, performance often deteriorated). The benefit depended on the
  classifier; RF and boosting did not benefit.
- Contribution: justifies per-study standardisation as our baseline batch adjustment. In
  leave-one-study-out the held-out study is standardised on its own statistics (label-free, batch-by-batch
  as defined). In CV and the locked test we use training-fold statistics of the same study.

**C3. "Simulating ComBat: how batch correction can lead to the systematic introduction of false positive
results in DNA methylation microarray studies" (PMC7328269).** Search-result snippet only, not read.
Not relied upon.

### D. Cross-disorder postmortem brain transcriptomics

**D1. Gandal MJ, Haney JR, Parikshak NN, et al., Geschwind DH (2018). "Shared molecular neuropathology
across major psychiatric disorders parallels polygenic overlap." Science 359:693–697.** Read via
PMC5898828 (methods and acknowledgements).
- 700 cortical samples: ASD 50, SCZ 159, BD 94, MDD 87, AAD 17, controls 293. The acknowledgements list
  exactly the 13 accessions in our candidate list plus E-MTAB-184. **The candidate list in the protocol is
  the Gandal list.**
- The paper gives no per-accession disorder, region or N table. It states that "Transcriptome summary
  statistics for each disorder were computed with a linear mixed-effects model to account for any sample
  overlap across studies". So the authors knew overlap existed but handled it statistically, not by
  deduplication.
- Reported gradient of transcriptomic severity: ASD > SCZ ≈ BD > MDD; overlap among SCZ, BD and MDD.
  Contribution: sets the expectation that SCZ vs BD is the hardest contrast and MDD vs control the weakest.
- **Conflict with GEO:** the list includes GSE11223. Per GEO, that series is "Colon biopsies from UC
  patients and healthy controls" (202 samples, GPL1708), which is not brain tissue and not psychiatric. It
  is presumably a typo for another accession; we cannot know which, so it is excluded. Also, the
  GSE28475 metadata describes a methods study (frozen vs fixed tissue, DASL vs IVT, reference RNA
  dilutions), not a case-control design.

**D2. Mistry M, Gillis J, Pavlidis P (2013). "Genome-wide expression profiling of schizophrenia using a
large combined cohort." Mol Psychiatry (2012).** Read via PMC3323740.
- Seven PFC SCZ datasets, including GSE17612 and GSE21138. "While the SMRI has additional data sets,
  these represent repeated runs of the samples from the same subjects, so we selected one dataset to
  represent each of the two SMRI brain collections."
- Effect sizes: "expression changes were small (~15% expression change)".
- Contribution: independent confirmation that multiple Stanley (SMRI) accessions re-profile the same
  subjects. It supports treating Stanley-derived series as one cohort.

**D3. "Ensemble Learning for Higher Diagnostic Precision in Schizophrenia Using Peripheral Blood Gene
Expression Profile." Neuropsychiatr Dis Treat (2024), PMC11075682.** Read (methods and results). A
peripheral-blood SCZ classifier with random train/test splits reported precision 80.4% on one held-out
dataset but accuracy 66.4%, and 59.9% precision cross-platform. Not brain data. Cited only as an
illustration that random-split performance drops on independent datasets.

### E. Donor overlap between brain-bank datasets and sample misannotation

**E1. Toker L, Feng M, Pavlidis P (2016). "Whose sample is it anyway? Widespread misannotation of samples
in transcriptomics studies." F1000Research 5:2103.** Read via PMC5034794.
- Method: XIST, KDM5D and RPS4Y1 with k-means, comparing gene-based sex with annotated sex.
- 32/70 datasets (46%) had at least one sex mismatch (83/4,160 samples, 2%). A 99% lower bound of 33% of
  studies. A further 21% of datasets had ambiguous samples.
- **Stanley:** of 4 datasets profiling the same SMRI subjects, 2 had mismatches, and "the mismatched
  subjects differed between the datasets" (lab mix-ups, not metadata errors).
- Contribution: our sex check follows this method (Y genes + XIST, per-dataset two-cluster call). We
  expect mismatches in Stanley-derived series and handle them per donor, per series.

**E2. "An online database for brain disease research" (2006), PMC1489945 (SMRI online genomics
database).** Read.
- The SMRI online genomics database compiles 12 studies from 12 labs on two collections (165 subjects),
  using Affymetrix HG-U133A, U133 Plus, HG-U95Av2, Agilent, CodeLink and a custom cDNA array.
- The article does not itself say explicitly that studies used overlapping subject subsets (the search
  snippet suggested it; the passage read did not). The overlap claim therefore rests on D2 and our own
  audit.

### F. Statistics for comparing models under repeated CV

**F1. Nadeau C, Bengio Y (2003). "Inference for the generalization error." Mach Learn 52:239–281**, as
implemented and described in the correctR vignette (CRAN) for repeated k-fold CV (Bouckaert & Frank 2004
form).
- Statistic: t = mean(d) / sqrt((1/(k·r) + n₂/n₁) · σ̂²), where d are the k·r per-fold differences, σ̂² their
  sample variance, and n₂/n₁ the test/train size ratio. We use df = k·r − 1.
- Assumption: the correlation between fold-level differences induced by overlapping training sets is
  approximated by ρ = n₂/(n₁+n₂) (equivalently the n₂/n₁ inflation). It is a heuristic, not an unbiased
  variance estimator (Nadeau & Bengio prove none exists). Using the 25 fold-level values (not repetition
  means) is how the correction is defined.
