# Biological Context

This project is a computational data-analysis exercise. It does not include any wet-lab work — no compounds were synthesized or tested, and no biophysical assays were run. This document explains where a model like this one fits in a real drug discovery workflow, and is explicit about what this project does and does not demonstrate.

## Why Binding Affinity Prediction Matters

Binding affinity (how strongly a compound binds a target protein) is one of the earliest filters in drug discovery. Before a compound reaches synthesis or wet-lab testing, computational estimates of binding affinity are used to narrow a large candidate pool down to a smaller set worth the time and cost of experimental follow-up. A useful computational model does not need to be perfectly accurate — it needs to reliably separate likely-poor candidates from those worth testing further.

## Computational Screening vs Experimental Validation

```text
Computational screening
        ↓
Candidate prioritization
        ↓
Experimental screening
        ↓
Biophysical validation
```

Computational prediction is a prioritization step, not a replacement for experimental validation. In practice, a model like the one in this repository would be used to reduce a large candidate list to a shorter one, which is then still subject to actual experimental screening and biophysical confirmation. It does not, and should not, replace those downstream steps.

## What This Project Does and Does Not Show

- **Does show**: end-to-end handling of a compound-protein binding dataset — data quality checks, target-leakage detection, feature selection under multicollinearity, model comparison with proper validation, and an honest assessment of where the model fails.
- **Does not show**: performance on real experimental affinity measurements, validation against any wet-lab assay, or evidence that the underlying molecular descriptors used here are sufficient for production-grade candidate ranking. The dataset is a virtual-screening / simulated dataset (see the main [README](../README.md) for the source), not a set of experimentally measured affinities.

Given the residual analysis findings documented in the main README (systematic underprediction of high-affinity compounds), the model in this repository is better positioned as a preliminary filtering tool — useful for deprioritizing candidates the model is confident are weak binders — rather than a final candidate-ranking model. Any real-world use would require validation against experimental data before being trusted for ranking decisions.
