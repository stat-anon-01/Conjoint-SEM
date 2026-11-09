# Understanding Conjoint Responses Through Structural Equation Modeling

Companion R code for *Understanding Conjoint Responses Through Structural Equation Modeling: A Tutorial on Assumptions and Model Comparisons*.

This repository has two purposes: to show how rating, binary, and paired-comparison conjoint models are specified in `lavaan`, and to reproduce the three demonstrations in the article. All examples use simulated data; no external dataset is required.

## Where to start

| Goal | Folder | What you will find |
|---|---|---|
| Learn the SEM syntax for each response format | [`01_basic_models/`](01_basic_models/) | Six standalone examples: rating, binary, and paired comparisons, each with homogeneous (HOMO) or heterogeneous (HETERO) residual variances |
| Reproduce the article's demonstrations | [`02_demonstrations/`](02_demonstrations/) | Separate data generation, model estimation, and table scripts for each demonstration |
| Examine numerical identification diagnostics | [`03_technical_checks/`](03_technical_checks/) | DEMO 1 identification/scale check and its saved summary |

## Requirements

Install R and the packages `lavaan` and `MASS`:

```r
install.packages(c("lavaan", "MASS"))
```

Run scripts with the working directory set to the folder containing them. For example, from the repository root:

```r
setwd("01_basic_models")
source("01_RATING_HOMO.R")
```

Each basic-model script automatically sources `00_DGP.R`; you can run the six examples independently. See the [basic-model guide](01_basic_models/README.md) for the complete list.

For the manuscript demonstrations, set the working directory to the relevant demonstration folder and follow the filenames below:

| Demonstration | Question | Execution order |
|---|---|---|
| [DEMO 1](02_demonstrations/demo01/) | Binary-response residual variance and threshold restrictions | `Demo01_00_DGP.R` → `Demo01_01_ANALYSIS.R` → `Demo01_02_TABLE.R` |
| [DEMO 2](02_demonstrations/demo02/) | Residual dependence in paired comparisons | `01_ANALYSIS.R` → `02_TABLE.R` (`01_ANALYSIS.R` sources `00_DGP.R`) |
| [DEMO 3](02_demonstrations/demo03/) | Comparing rating and binary preference parameters | `01_ANALYSIS.R` → `02_TABLE.R` (`01_ANALYSIS.R` sources `00_DGP.R`) |

For example, from `02_demonstrations/demo01`:

```r
source("Demo01_00_DGP.R")
source("Demo01_01_ANALYSIS.R")
source("Demo01_02_TABLE.R")
```

The [DEMO 1 technical check](03_technical_checks/Demo01_Identification_Lavaan_Check.R) can be run from the `03_technical_checks` directory. It loads the original DEMO 1 scripts and writes `Demo01_Identification_Lavaan_Results.csv`.

## Interpretation notes

- `01_basic_models` is instructional. Its data and analyses are **not** the results reported for DEMO 1–3.
- In the basic examples, HOMO and HETERO are alternative *fitted models* for the same generated data. The shared DGP has independent, equal-variance residuals.
- Binary and paired-comparison HOMO/HETERO examples use different latent-response scale anchors. Do not compare raw, unstandardized estimates across those models without adjusting for their normalizations.
- Numerical full rank at an estimated parameter point supports **local**, not global, identification. The technical check is a numerical diagnostic rather than a general proof.

The author has run the basic examples and the three demonstrations, and has run the DEMO 1 identification check. A separate, independent reproduction in a clean environment has not been established here. The provided scripts retain their original simulation seeds and model specifications.

## Citation

Please cite the accompanying article when using these materials. Add its final bibliographic details and the public repository URL when available.
