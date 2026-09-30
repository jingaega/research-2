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
