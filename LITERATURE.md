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
| S18 | 2026-09-30 | PubMed esearch: "Adjusting batch effects in microarray expression data using empirical Bayes methods[Title]" | E-utilities |
| S19 | 2026-09-30 | PubMed: Luo 2010 comparison of batch effect removal methods ... MAQC-II | E-utilities |
| S20 | 2026-09-30 | PubMed: McCall Bolstad Irizarry 2010 Frozen robust multiarray analysis fRMA | E-utilities |
| S21 | 2026-09-30 | PubMed: Leek 2010 Tackling the widespread and critical impact of batch effects | E-utilities |
| S22 | 2026-09-30 | PubMed: "MixupMapper[Title/Abstract]" (first query with year/author terms returned nothing) | E-utilities |
| S23 | 2026-09-30 | PubMed: Iwamoto Bundo Kato mitochondria-related genes ... 2005 → returned the WRONG paper (Munakata 2005, LARS2); retried as 15563509[uid] (the PMID GEO lists for GSE12649) | E-utilities |
| S24 | 2026-09-30 | PubMed: Chen C 2013 two gene co-expression modules differentiate psychotics and controls | E-utilities |
| S25 | 2026-09-30 | PubMed: Maycox 2009 analysis of gene expression in two large schizophrenia cohorts | E-utilities |
| S26 | 2026-09-30 | PubMed: Narayan 2008 molecular profiles of schizophrenia in the CNS at different stages of illness | E-utilities |
| S27 | 2026-09-30 | PubMed: Hagenauer 2018 inference of cell type content from human brain transcriptomic datasets | E-utilities |
| S28 | 2026-09-30 | PubMed: Ramaker 2017 post-mortem molecular profiling of three psychiatric disorders | E-utilities |
| S29 | 2026-09-30 | PubMed: Gandal 2018 transcriptome-wide isoform-level dysregulation in ASD, schizophrenia, and bipolar disorder | E-utilities |
| S30 | 2026-09-30 | PubMed: Varma Simon 2006 bias in error estimation ... → returned Tsamardinos 2018 (kept, relevant); retried as 16504092[uid] | E-utilities |
| S31 | 2026-09-30 | PubMed: Dietterich 1998 approximate statistical tests for comparing supervised classification learning algorithms | E-utilities |
| S32 | 2026-09-30 | PubMed: Kapoor Narayanan 2023 leakage and the reproducibility crisis in machine-learning-based science | E-utilities |
| S33 | 2026-09-30 | PubMed: Wang 2021 MOGONET integrates multi-omics data using graph convolutional networks | E-utilities |
| S34 | 2026-09-30 | PubMed by PMID: 31123247 (Lanz 2019), 24608543 (Chang/Sibille 2014) | E-utilities |
| S35 | 2026-09-30 | arXiv abstract pages: 1902.07153 (SGC), 1912.09893 (fair GNN comparison), 1711.05859 (Rhee), 1811.05868 (GNN evaluation pitfalls) | WebFetch |
| S36 | 2026-09-30 | Europe PMC full-text XML: PMC6533277 (Lanz 2019), PMC3946570 (Chang 2014); PMC2783475 (Narayan 2008) not retrievable (PMC CAPTCHA; Europe PMC returned no text/PDF) | curl |

Abstracts retrieved via E-utilities are saved verbatim in `data/lit/pmid_<id>.txt` (git-ignored; regenerate with
`scripts/19_fetch_literature.py`), and every query with its returned PMIDs is logged in `results/literature_queries.json`.

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

### A (continued). Graph models: evidence and design sources

**A5. Rhee S, Seo S, Kim S (2018). "Hybrid approach of relation network and localized graph convolutional
filtering for breast cancer subtype classification." IJCAI 2018 (arXiv:1711.05859).** Abstract only.
- A graph CNN on a PPI network plus a relation network for PAM50 subtypes; claims "significantly better
  performances than existing methods". The abstract names no baselines and no random-graph control.
  Pro-graph claim, weak controls (same pattern as A1 and A3).

**A6. Wang T, Shao W, Huang Z, et al., Huang K (2021). "MOGONET integrates multi-omics data using graph
convolutional networks allowing patient classification and biomarker identification." Nat Commun
12:3445.** Abstract read (PMID 34103512).
- Claims to outperform other multi-omics integrative methods. **Its graphs are patient-similarity graphs,
  not gene graphs.** That is a different use of structure from ours (a sample graph would couple test
  donors to training donors and needs care to avoid transductive leakage). Not adopted.

**A7. Wu F, Zhang T, de Souza AH Jr, Fifty C, Yu T, Weinberger KQ (2019). "Simplifying Graph
Convolutional Networks." ICML 2019 (arXiv:1902.07153).** Abstract read.
- A GCN reduces to "a fixed low-pass filter followed by a linear classifier" with competitive accuracy
  and large speedups.
- **Design contribution:** variants V1/V2 (graph-smoothed features → logistic regression) are exactly
  this SGC form, applied over the gene graph. That makes the graph the *only* difference from the
  no-structure baseline, and makes the density-matched random-graph control directly comparable.

### B (continued). Skeptical / methodological evidence on graph models

**B5. Errica F, Podda M, Bacciu D, Micheli A (2020). "A Fair Comparison of Graph Neural Networks for
Graph Classification." ICLR 2020 (arXiv:1912.09893).** Abstract read.
- Over 47,000 experiments: under a rigorous, identical protocol, structure-agnostic baselines were
  competitive, and on some datasets "structural information has not been exploited yet". Supports our
  equal-treatment protocol and the no-graph control.

**B6. Shchur O, Mumme M, Bojchevski A, Günnemann S (2018). "Pitfalls of Graph Neural Network
Evaluation." arXiv:1811.05868.** Abstract read.
- A single fixed train/test split gives unstable model rankings ("different splits of the data leads to
  dramatically different rankings"). Simpler GNNs beat complex ones when all are tuned fairly. Supports
  repeated CV with identical folds for every model and identical tuning budgets.

### C (continued). Batch correction, normalisation and prediction

**C4. Johnson WE, Li C, Rabinovic A (2007). "Adjusting batch effects in microarray expression data using
empirical Bayes methods." Biostatistics 8:118–127.** Abstract read (PMID 16632515).
- Earlier methods "require large batch sizes (> 25)". Empirical-Bayes location/scale adjustment is
  "robust to outliers in small sample sizes". This is the mechanism behind V5: our smallest batches
  (e.g. GSE54567, about 20 dev donors) are where shrinkage could help.

**C5. Luo J, Schumacher M, Scherer A, et al., Zhang J (2010). "A comparison of batch effect removal
methods for enhancement of prediction performance using MAQC-II microarray gene expression data."
Pharmacogenomics J 10:278–291.** Abstract read (PMID 20676067).
- Over 120 cross-batch prediction cases (SVM, KNN, MCC), Ratio-G, Ratio-A, EJLR, **mean-centering and
  standardization** performed "better or equivalent to no batch effect removal" in 89/85/83/79/75% of
  cases. This supports per-batch standardisation as the baseline, and shows that it can also hurt
  (about 25% of cases), which is why V5 compares alternatives.

**C6. McCall MN, Bolstad BM, Irizarry RA (2010). "Frozen robust multiarray analysis (fRMA)."
Biostatistics 11:242–253.** Abstract read (PMID 20097884).
- RMA's normalisation and summarisation "require multiple arrays to be analyzed simultaneously". fRMA
  uses precomputed frozen probe effects to process arrays "individually or in small batches". It is
  "comparable to RMA when the data are analyzed as a single batch and outperforms RMA when analyzing
  multiple batches". This is **the justification for our preprocessing choice**: per-array
  normalisation cannot leak across folds.

**C7. Leek JT, Scharpf RB, Bravo HC, et al., Irizarry RA (2010). "Tackling the widespread and critical
impact of batch effects in high-throughput data." Nat Rev Genet 11:733–739.** Abstract read
(PMID 20838408).
- Batch effects become "a major problem when batch effects are correlated with an outcome of interest".
  In our pooled design batch (study/cohort) *is* correlated with diagnosis. This is the rationale for
  the study-only baseline and the study-identity probe.

### D (continued). The datasets' own papers and cross-disorder transcriptomics

**D4. Iwamoto K, Bundo M, Kato T (2005). "Altered expression of mitochondria-related genes in postmortem
brains of patients with bipolar disorder or schizophrenia..." Hum Mol Genet 14:241–253** (PMID 15563509,
the PMID GEO lists for **GSE12649**). Abstract read.
- Global down-regulation of mitochondrial genes in BD and SZ even after controlling for pH, which the
  authors say was "likely due to the effects of medication". **Contribution:** medication is a
  plausible confounder of any classifier signal (treated patients vs untreated controls). The data
  cannot separate it from disease; this is noted as a limitation.
- (S23 initially retrieved Munakata et al. 2005, Biol Psychiatry, a LARS2 study from the same lab. It
  states that its prefrontal cortices came from "the Stanley Foundation Brain Collection", which is
  consistent with GEO's GSE12649 text.)

**D5. Chen C, Cheng L, Grennan K, et al., Liu C (2013). "Two gene co-expression modules differentiate
psychotics and controls." Mol Psychiatry 18:1308–1314** (PMID 23147385, the GEO PMID for **GSE35978**).
Abstract read.
- WGCNA modules (neuron differentiation and development; neuron protection) differed between patients
  and controls in cortex and cerebellum, "preserved in five expression data sets". Suggests SCZ/BD
  signal is spread over co-expression modules, which motivates the expression-derived graph (V2).

**D6. Maycox PR, Kelly F, Taylor A, et al., de Belleroche J (2009). Mol Psychiatry 14:1083–1094**
(PMID 19255580, GEO PMID for **GSE17612**). Abstract read.
- Tissue source: "Charing Cross Hospital prospective collection", which confirms our cohort name.
  It compared this cohort with a Harvard Brain Bank PFC dataset: 51 common changes, 49 in the same
  direction, "multiple, small but synergistic changes". This is the pattern under which aggregation
  (graph or pathway) *could* help.
- (Note: GEO does not state the brain bank for GSE17612; the name comes from this paper.)

**D7. Narayan S, Tang B, Head SR, et al., Dean B, Thomas EA (2008). "Molecular profiles of schizophrenia
in the CNS at different stages of illness." Brain Res 1239:235–248** (PMID 18778695, GEO PMID for
**GSE21138**). Abstract read; full text not retrievable.
- Expression differences are largest in short-duration illness, and the stages differ in the systems
  affected. **Contribution:** SCZ is heterogeneous by illness stage within one cohort, so any SCZ
  classifier's performance will depend on the stage mix per cohort.
- **The brain bank is not stated in GEO or in the abstract.** Our internal cohort key "Victoria" came from
  background knowledge (co-author B. Dean), **not from a source read here**. The key is kept (renaming
  would change the locked split's stratum encoding), but the cohort's provenance is unverified.

**D8. Lanz TA, Reinhart V, Sheehan MJ, et al., Lewis DA, Kleiman RJ (2019). Transl Psychiatry 9:151**
(PMID 31123247, GEO PMID for **GSE53987**). Abstract and methods (Europe PMC full text) read.
- Methods: "Brains from tetrads of subjects with SCZ, MDD, or BD and unaffected comparison subjects were
  obtained at the University of Pittsburgh". All samples were collected "at the coroner's office". 19
  tetrads matched for sex and age.
- The largest SCZ effects were in hippocampus. **Conflict with our data:** our identity audit shows the
  GSE53987 hippocampal arrays are mislabelled relative to PFC/STR (hippocampal arrays match other
  recorded donors' PFC/STR arrays at z up to 8.7). The paper's hippocampus-specific findings may be
  affected. We report this conflict and do not use the hippocampus.

**D9. Chang L-C, Jamain S, Lin C-W, Rujescu D, Tseng GC, Sibille E (2014). PLoS ONE 9:e90980**
(PMID 24608543, GEO PMID for **GSE54567/68/71/72**). Abstract and methods table (Europe PMC) read.
- Meta-analysis of 11 MDD cohorts, all from the Sibille lab (University of Pittsburgh). Table 1 lists
  BA25_F and BA25_M (26 subjects each), BA9_F and BA9_M, and the MD1–MD3 cohorts.
- Its table gives BA9_F as 32 subjects, while GEO GSE54568 has 30 samples, and BA25_M as 26 while GEO
  GSE54572 has 24 samples. **Conflict: paper vs GEO sample counts; GEO is used.**
- The passages read do not state whether subjects overlap with other Pittsburgh studies. Our audit
  finds ≥ 20 donors shared with GSE53987.

**D10. Hagenauer MH, Schulmann A, Li JZ, et al., Akil H (2018). PLoS ONE 13:e0200003** (PMID 30016334,
GEO PMID for **GSE92538**). Abstract read.
- Principal components of variation in GSE92538, GSE53987, GSE21935, GSE21138 and CommonMind "strongly
  correlated with the predicted neuronal/glial content". "Prolonged hypoxia around the time of death
  predicted increased astrocytic and endothelial gene expression".
- **Contribution:** (i) independently explains why our first identity attempt was dominated by
  agonal/hypoxia genes until top PCs were removed. (ii) Cell-type composition is a major axis that can
  differ by diagnosis *and* by cohort and agonal state, so a classifier may exploit composition. This is
  a reserve variant (cell-type marker scores only) to test how much signal is just composition.

**D11. Ramaker RC, Bowling KM, Lasseigne BN, et al., Myers RM (2017). "Post-mortem molecular profiling of
three psychiatric disorders." Genome Med 9:72.** Abstract read (PMID 28754123). RNA-seq, so used as
literature only, not data.
- 24 subjects per group (SCZ, BD, MDD, CTL) in 3 regions. The strongest effects were in SCZ; "broad
  down-regulation of genes specific to neurons and concordant up-regulation of genes specific to
  astrocytes" in SCZ and BD.
- Sets the expectation that MDD vs CTL is the weakest contrast, and that SCZ/BD signal is partly a
  neuron/astrocyte composition shift.

**D12. Gandal MJ, Zhang P, Hadjimichael E, et al.; PsychENCODE Consortium; Geschwind DH (2018).
"Transcriptome-wide isoform-level dysregulation in ASD, schizophrenia, and bipolar disorder." Science
362:eaat8127.** Abstract read (PMID 30545856). RNA-seq, literature only.
- 1,695 individuals. "Isoform-level changes capturing the largest disease effects", with microglial,
  astrocyte and interferon modules. **Contribution:** gene-level microarray summaries (our only option
  across these platforms) discard the isoform level where effects are largest. This is a stated
  ceiling on what our features can capture.

### E (continued). Sample mix-ups

**E3. Westra H-J, Jansen RC, Fehrmann RSN, et al., Franke L (2011). "MixupMapper: correcting sample
mix-ups in genome-wide datasets increases power to detect small genetic effects." Bioinformatics
27:2104–2111.** Abstract read (PMID 21653519).
- Uses cis-eQTL genotype–expression links. On average 3% of samples had incorrect expression
  phenotypes, and 23% in one dataset.
- **Contribution:** the principle behind our identity audit (genotype-driven genes identify
  individuals). Our rates (for example 5 of 68 GSE53987 hippocampal arrays with sex conflicts, plus
  identity mismatches) are in the range they report.

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

### F (continued). Evaluation statistics, tuning bias, leakage

**F2. Dietterich TG (1998). "Approximate statistical tests for comparing supervised classification
learning algorithms." Neural Comput 10:1895–1923.** Abstract read (PMID 9744903).
- A paired t-test over repeated random train/test splits has "high probability of type I error ... and
  should never be used". The 10-fold CV t-test has "somewhat elevated" type I error. This is the
  problem the Nadeau–Bengio correction (F1) addresses, and why we do not use an uncorrected t-test on
  our 25 folds as the primary test.

**F3. Varma S, Simon R (2006). "Bias in error estimation when using cross-validation for model
selection." BMC Bioinformatics 7:91.** Abstract read (PMID 16504092).
- On null data, the CV error of the tuned classifier was < 30% on 18.5% (shrunken centroids) and 38%
  (SVM) of datasets, while "performance ... on the independent test set was no better than chance".
  Nested CV "reduces the bias considerably".
- **Contribution:** hyperparameters are tuned by an inner CV nested in each outer fold. The dev-CV score
  of whichever configuration we select as "best" among ≤ 12 is itself optimistically biased, which is
  why only the single locked-test evaluation is unbiased. PLAN.md says this explicitly.

**F4. Tsamardinos I, Greasidou E, Borboudakis G (2018). "Bootstrapping the out-of-sample predictions for
efficient and accurate cross-validation." Mach Learn 107:1895–1922.** Abstract read (PMID 30393425).
- "The cross-validated performance of the best configuration is optimistically biased". It proposes
  bootstrap bias correction (BBC-CV). **Contribution:** we report the BBC-CV-corrected dev estimate for
  the selected configurations as a supplementary number. It is computed only from stored out-of-fold
  predictions and needs no extra model fits.

**F5. Kapoor S, Narayanan A (2023). "Leakage and the reproducibility crisis in machine-learning-based
science." Patterns 4:100804.** Abstract read (PMID 37720327).
- Leakage found in 17 fields and 294 papers, with a taxonomy of eight leakage types. In their
  reproduction, "when the errors are corrected, complex ML models do not perform substantively better
  than decades-old LR models". **Contribution:** our leakage checklist (preprocessing, duplicates, test
  set, feature selection, temporal/batch) follows their taxonomy, and it supports logistic regression as
  the no-structure baseline.

## Source count

Read (full text or the quoted passage/abstract; "abstract only" as marked): A1–A7, B1–B6, C1, C2, C4–C7,
D1–D12, E1–E3, F1–F5 = **39 sources read** (15 with full text or quoted passages, 24 abstract-level). Snippet-only, not relied on: C3.
Graph-model evidence is weighed as **pro: A1, A2, A3, A5, A6 (A7 is neutral design)** vs
**skeptical: B1–B6**. The pro papers either lack random-graph and no-graph controls (A1, A3, A5) or use
patient graphs (A6). The only pro paper with a random-network control (A2) failed to replicate under
re-evaluation (B4). Prior unchanged: a gene-graph gain beyond a density-matched random graph is
unlikely at our sample sizes.

## Conflicts between papers and GEO / our data (reported per protocol)

1. Gandal 2018 (D1) accession list includes **GSE11223**, which GEO describes as colon biopsies (UC).
2. Chang 2014 (D9) Table 1 subject counts (BA9_F 32, BA25_M 26) ≠ GEO sample counts (GSE54568: 30;
   GSE54572: 24). GEO used.
3. Lanz 2019 (D8) reports hippocampus as the region with the largest SCZ effects. Our expression
   identity audit shows the deposited GSE53987 hippocampal arrays do not belong to the donors their
   labels say.
4. Mistry 2013 (D2) and Gandal 2018 (D1) treat Stanley overlap as known. Neither the Lanz nor the Sibille
   papers read mention shared Pittsburgh donors, but our audit finds ≥ 20.
