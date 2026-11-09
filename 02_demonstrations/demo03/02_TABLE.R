# Demonstration 3: Tables from fitted cross-format models
# Run 01_ANALYSIS.R first. Requires lavaan.

library(lavaan)
obj <- readRDS("Demo03_results.rds")
fit_M0 <- obj$fits$M0
fit_M1 <- obj$fits$M1
fit_M2 <- obj$fits$M2
fit_M3 <- obj$fits$M3

# Model fit: Table 4 in the manuscript.
indices <- c("chisq.scaled", "df.scaled", "cfi.scaled", "rmsea.scaled", "srmr")
model_fit <- rbind(
  M0 = fitMeasures(fit_M0, indices),
  M1 = fitMeasures(fit_M1, indices),
  M2 = fitMeasures(fit_M2, indices),
  M3 = fitMeasures(fit_M3, indices)
)
print(round(model_fit, 3))

# WLSMV scaled difference tests for the nested population-level restrictions.
# M2 vs. M3 is not tested using a conventional chi-square difference test.
test_M0_M1 <- lavTestLRT(fit_M0, fit_M1, method = "satorra.2000")
test_M1_M2 <- lavTestLRT(fit_M1, fit_M2, method = "satorra.2000")
print(test_M0_M1)
print(test_M1_M2)

# M2 population-level preference estimates and cross-format correlations.
parameters_M2 <- parameterEstimates(fit_M2, ci = TRUE)
preference_M2 <- subset(parameters_M2,
  grepl("^[RB]_F[0-3]$", lhs) & op %in% c("~1", "~~"))
standardized_M2 <- standardizedSolution(fit_M2, type = "std.all", ci = TRUE)
cross_format <- subset(standardized_M2,
  op == "~~" & grepl("^R_F[0-3]$", lhs) & grepl("^B_F[0-3]$", rhs))
print(preference_M2)
print(cross_format)

# Save results for inspection and reuse.
dir.create("Demo03_tables", showWarnings = FALSE)
write.csv(model_fit, "Demo03_tables/model_fit.csv")
write.csv(as.data.frame(test_M0_M1), "Demo03_tables/M0_vs_M1.csv")
write.csv(as.data.frame(test_M1_M2), "Demo03_tables/M1_vs_M2.csv")
write.csv(preference_M2, "Demo03_tables/preference_parameters.csv", row.names = FALSE)
write.csv(cross_format, "Demo03_tables/cross_format_correlations.csv", row.names = FALSE)
