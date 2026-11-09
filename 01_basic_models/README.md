# Basic conjoint SEM examples

These six R/`lavaan` scripts illustrate how the same three-attribute, eight-profile design can be modeled for three response formats. They are **teaching examples**, separate from the article's three demonstrations.

| Response format | Equal residual variances (HOMO) | Freely varying residual variances (HETERO) |
|---|---|---|
| Ratings | [`01_RATING_HOMO.R`](01_RATING_HOMO.R) | [`02_RATING_HETERO.R`](02_RATING_HETERO.R) |
| Binary evaluations | [`03_BINARY_HOMO.R`](03_BINARY_HOMO.R) | [`04_BINARY_HETERO.R`](04_BINARY_HETERO.R) |
| Paired comparisons | [`05_PAIRED_HOMO.R`](05_PAIRED_HOMO.R) | [`06_PAIRED_HETERO.R`](06_PAIRED_HETERO.R) |

## Run an example

Install `lavaan` and `MASS` if needed. Set your working directory to `01_basic_models`, then run **any one** of the six scripts:

```r
source("01_RATING_HOMO.R")
```

Each script automatically sources `00_DGP.R`, which generates the common simulated sample (N = 1,000; seed = 12345). There is no need to run the scripts in numerical order.

## What HOMO and HETERO mean

HOMO constrains the fitted residual variances to be equal; HETERO relaxes that constraint. The data are the **same** within each response format, and the DGP itself uses independent, equal-variance residuals. HETERO is not a second data-generating condition.

For binary evaluations and paired comparisons, the scripts employ different latent-response scale normalizations. Consequently, the unstandardized HOMO and HETERO estimates cannot be compared directly without rescaling. In paired comparisons, the common baseline factor cancels when one profile's utility is subtracted from another's.

For the analyses and tables reported in the article, see [`02_demonstrations`](../02_demonstrations/). For the numerical identification check of DEMO 1, see [`03_technical_checks`](../03_technical_checks/).
