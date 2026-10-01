# When do graph/GNN models beat both standard ML and random graphs? A targeted review

Question (user request, 2026-10-01): across the literature and our own references, which graph-based
models for molecular-profile prediction succeeded, beat general ML, **and** beat a random-graph
control, and what made them succeed?

Only sources whose relevant passages I actually read are used (full text via Europe PMC or arXiv;
see LITERATURE.md §G for the search log). "Beat random" means a structure-matched random or permuted
graph control was reported.

## 1. Scorecard

| # | Study | Task / data | n | Beat general ML? | Beat random graph? | Notes |
|---|---|---|---|---|---|---|
| S1 | Erion et al., attribution priors (arXiv 1906.10670) | AML drug response (IC50) from expression, regression | moderate (one AML cohort) | **Yes**, significant (R²) | **Yes**: a randomized graph gave "no better than … a neural network trained without any attribution prior" | Graph = **tissue-specific HumanBase** functional network; used as a *smoothness prior on attributions*, not as an architecture |
| S2 | Proteomics GNN, UK Biobank (PMC12676398) | HbA1c and other traits from plasma proteomics | ~23,400 (subsampled to 11,715 and 5,858) | **Yes** (vs linear and other deep models) | **Yes**: the degree-preserving permuted GNN was worse at every training size | Graph = **curated GO Molecular Function** sets; large n; continuous, molecularly proximal phenotype |
| S3 | KPNN, Fortelny & Bock 2020 (PMC7397672) | T-cell receptor stimulation from single-cell RNA-seq | thousands of cells | Comparable (the goal was interpretability) | **Yes, "generally"**: edge-shuffled controls had reduced performance | Graph = **the TCR signalling network itself**; the label is the activation of that very pathway. The authors note shuffled controls still retained some biology |
| S4 | GINCCo (PMC8826027) | METABRIC breast cancer subtype, grade and relapse | 1,980 | Mostly **not significant**; significant only for IC10 subtype | **Yes**: beat randomly constructed computational graphs | Gain explained as sparsity (< 0.05% of MLP parameters) **plus** complex structure |
| S5 | Brouard et al. 2024 simulated data (PMC10939369) | Phenotype simulated *from a known gene network* | 100 | GNN > MLP and RF | **Yes**, only when the **true generating network** is used | On real data (PPI, co-expression) random or complete graphs were as good or better |
| S6 | Chuang et al. 2007 (our ref A2) | Breast cancer metastasis, cross-cohort | 286–295 per cohort | Yes (vs single-gene markers) | **Yes** (P = 0.046 / 0.012 vs random subnetworks) | Not reproduced across 6 cohorts by Staiger 2012 (B4) |
| S7 | L1000 GNN (PMC12161499) | Infer full transcriptome from 970 landmark genes | 2,500 train | **Yes** (vs LR, MLP) | Not tested | Task = predicting genes from genes, so gene–gene edges *are* the target relationship |
| — | P-NET, Elmarakeby 2021 (PMC8514339) | Prostate cancer state from mutations/CNV | 1,013 | Yes (vs SVM, LR, trees, dense NN; strongest at n ≤ 500) | **Not tested in the original** (compared with dense nets only) | A re-evaluation found randomized sparse versions reached the same AUC (0.896 vs 0.899) |
| — | "Sparsity is all you need" (PMC13446513) | 29 pathway-informed networks | many | — | **No**: structure-matched randomized versions "consistently match or outperform" | Gains "arise predominantly from sparsity-induced regularization rather than from biological knowledge" |
| — | Chereda 2019, Ramirez 2020, Rhee 2018, MOGONET 2021 (our refs) | Cancer classification | 969 – 11k | Claimed | **Not tested** | Positive claims without random controls |
| — | Brouard 2024 (real data), Staiger 2012, **this study** | Various, incl. our psychiatric tasks | 185 – 156k | Rarely | **No** | Random/complete graphs ≈ or > biological graphs |

## 2. What the successes have in common

1. **The graph encodes the mechanism that produces the label (task–graph alignment).** This is the
   strongest common factor.
   - KPNN predicts TCR stimulation using the TCR signalling network.
   - The AML model uses a network for the tissue of the disease.
   - HbA1c (a glycaemic marker) uses molecular-function sets of circulating proteins.
   - Brouard's simulation only helped when the graph *was* the generative network.
   - L1000 predicts genes from genes.
   - When the graph is a generic interaction map with no specific link to the label (PPI or
     co-expression for a complex clinical outcome), random graphs do as well (Brouard real data;
     Staiger; 29-model re-evaluation; us).
2. **Context-specific, functional networks rather than generic PPI.** S1's gain came from a
   tissue-specific HumanBase network, and its random version gave nothing. Generic STRING/HPRD PPI
   (Chereda, Brouard, us) never beat random.
3. **The graph acts as a regulariser, and much of the "beat ML" gain is sparsity.**
   - P-NET's advantage over dense networks was significant only at small training sizes (n ≤ 500;
     e.g. P = 0.004 at n = 155).
   - GINCCo uses < 0.05% of an MLP's parameters.
   - The 29-model study shows random sparsity reproduces most of these gains.
   - **Beating general ML ≠ beating random graphs.** Most papers show only the first.
4. **Enough samples relative to signal, with a strong, molecularly proximal signal.**
   - The successes with clean random controls had large n (S2 ≈ 23k; S3 thousands of cells) or a
     label that is a direct biochemical readout (drug IC50, HbA1c, pathway stimulation).
   - S2 showed the prior's benefit persisted as data shrank, but "limited" there still meant
     5,858 samples.
5. **Proper controls existed.** The convincing successes used *degree- or structure-preserving*
   permutations (S2, S3) or density-matched random graphs (S1, S4). Positive claims without such
   controls (Chereda, Ramirez, Rhee, MOGONET, original P-NET) cannot separate biology from
   sparsity, capacity or regularisation.
6. **Single-cohort evaluation.** None of the successes faced cross-study batch confounding.
   Their CV did not have to transfer across brain banks or platforms, as ours did.

## 3. Why our setting missed every success condition

| Success condition | Our study |
|---|---|
| Graph aligned with the label mechanism | Generic STRING ≥ 700 PPI, Reactome, in-fold co-expression; no known mechanism for SCZ/BD/MDD in bulk cortex |
| Context-specific network | Not brain- or cell-type-specific |
| Strong, proximal signal | Weak, diffuse; ~88% of the SCZ signal is cell-type composition (V9); BD/MDD near chance |
| Enough samples | 185–268 dev donors per task |
| Model where the prior can act as regulariser | Already a strongly L2-regularised linear model (C ≈ 1e-4), so extra smoothing has little to add |
| Single homogeneous cohort | 5 cohorts, 3 platforms, batch-confounded |

Our null result for graphs (Δ ≈ 0 vs no-graph and vs random) is therefore what this literature
predicts, not an anomaly.

## 4. If one wanted to give graphs a fair chance here (not run; would need new pre-registration)

1. **Brain- and cell-type-specific networks** built from independent data, never from our training
   folds: for example HumanBase brain/cortex functional networks, or PsychENCODE brain regulatory
   networks.
2. **Graph as an attribution or smoothness prior** on a regularised model (the S1 design), rather than
   graph convolution.
3. **Sparse random-connectivity and degree-preserving permuted controls** as the primary comparison,
   not only Erdős–Rényi. Report "beats random sparsity" separately from "beats dense ML".
4. **Expected ceiling.** With n ≈ 200 per task and signal dominated by composition, any gain would
   likely be < 0.03 macro-F1, below this study's noise floor. The locked test set has already been
   used, so any such analysis would be exploratory and evaluated on CV and leave-one-cohort-out
   only.
