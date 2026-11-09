# Demonstration 2: Tables from fitted paired-comparison models
# Run 01_ANALYSIS.R first. Requires lavaan.

library(lavaan)
obj <- readRDS("Demo02_results.rds")
fit_A <- obj$results$A$fit
fit_B <- obj$results$B$fit
fit_C <- obj$results$C$fit

# Model fit: Table 3 in the manuscript.
indices <- c("chisq.scaled", "df.scaled", "cfi.scaled", "rmsea.scaled", "srmr")
model_fit <- rbind(
  A = fitMeasures(fit_A, indices),
  B = fitMeasures(fit_B, indices),
  C = fitMeasures(fit_C, indices)
)
print(round(model_fit, 3))

# Parameter estimates for the best-fitting combined residual model.
parameters_C <- parameterEstimates(fit_C, ci = TRUE)
preference_C <- subset(parameters_C,
  lhs %in% paste0("F", 1:3) & op %in% c("~1", "~~"))
residual_C <- subset(parameters_C,
  op == "~~" & lhs == rhs &
    lhs %in% c(paste0("U", 1:8), obj$pn))
print(preference_C)
print(residual_C)

# Save transparent, machine-readable results.
dir.create("Demo02_tables", showWarnings = FALSE)
write.csv(model_fit, "Demo02_tables/model_fit.csv")
write.csv(preference_C, "Demo02_tables/preference_parameters.csv", row.names = FALSE)
write.csv(residual_C, "Demo02_tables/residual_variances.csv", row.names = FALSE)
