# research-2: pooled psychiatric postmortem brain transcriptomics

Start with **REPORT.md** (findings). Then:

- **PLAN.md**: the approved plan.
- **LOG.md**: every decision, pre-registration and result, in order.
- **DATA_AUDIT.json**: per-series audit, donor deduplication, sex checks, locked class list.
- **LITERATURE.md**: searches and the 39 sources read.
- `results/report_tables.md`, `results/final_comparison.json`, `results/final_test_results.json`:
  final tables.
- `results/variants.json`: every variant, with the running budget count.
- `results/summary/*.json`: per-configuration CV and leave-one-cohort-out summaries.
- `results/runs/`: one JSON per configuration, task, fold and seed.
- `scripts/`: the numbered pipeline. Raw data is not committed; regenerate it with scripts 01–10.
