# Insurance Claims Stochastic Modeling Analysis
# Author: Francesco
# Date: 28/08/2025
# 
# This script applies stochastic modeling techniques to insurance claims data:
# - Frequency Modeling with discrete distributions
# - Severity Modeling with continuous distributions  
# - Monte Carlo Simulation with variance reduction
# - Risk Assessment (VaR and Combined Ratios)

# Load required libraries
library(fitdistrplus)
library(ggplot2)
library(knitr)
library(dplyr)

# ==============================================================================
# DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load and examine the dataset
data <- read.csv("DATA_SET_4.csv", sep = ";")
cat("Original dataset dimensions:", dim(data), "\n")

# Data cleaning pipeline
data$X <- NULL  # Remove unnecessary column

# Remove policies with zero claims (no claims experience)
data <- data[data$CLM_FREQ != 0, ]

# Clean premium column - remove spaces and convert to numeric
data$PREMIUM <- as.numeric(gsub(" ", "", data$PREMIUM))

# Clean claim amount columns (CLM_AMT_1 through CLM_AMT_8)
for (i in 1:8) {
  col_name <- paste0("CLM_AMT_", i)
  data[[col_name]] <- gsub(" ", "", data[[col_name]])  # Remove spaces
  data[[col_name]][data[[col_name]] %in% c("-", "")] <- "0"  # Replace missing with 0
  data[[col_name]] <- as.integer(data[[col_name]])  # Convert to integer
  data <- data[data[[col_name]] >= 0, ]  # Remove negative values
}

# Calculate total claim amounts per policy
data$TOTAL_CLM_AMT <- rowSums(data[, paste0("CLM_AMT_", 1:8)], na.rm = TRUE)

cat("Cleaned dataset dimensions:", dim(data), "\n")

# ==============================================================================
# CLAIM FREQUENCY ANALYSIS
# ==============================================================================

# Visualize claim frequency distribution
frequency_plot <- ggplot(data, aes(x = CLM_FREQ)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Claim Frequencies",
       x = "Number of Claims", y = "Frequency") +
  theme_minimal()
print(frequency_plot)

# Fit discrete distributions to claim frequency
fit_poisson <- fitdist(data$CLM_FREQ, "pois", method = "mle")
fit_nbinom <- fitdist(data$CLM_FREQ, "nbinom", method = "mle")
fit_geom <- fitdist(data$CLM_FREQ, "geom", method = "mle")

# Compare model fits using AIC and Chi-squared
fit_list <- list(Poisson = fit_poisson, 
                 NegBinomial = fit_nbinom,
                 Geometric = fit_geom)

gof <- gofstat(fit_list)

# Create comparison table
comparison_table <- data.frame(
  AIC = round(gof$aic, 2),
  Chi_Squared = round(gof$chisq, 2)
) %>% arrange(AIC)

print("Model Comparison for Claim Frequency:")
print(comparison_table)
# Result: Poisson distribution provides the best fit based on AIC criterion

# Visualize fitted vs empirical distributions
obs_counts <- table(data$CLM_FREQ)
obs_probs <- obs_counts / sum(obs_counts)
x_vals <- as.integer(names(obs_counts))

# Calculate fitted probabilities for each distribution
pois_probs <- dpois(x_vals, lambda = fit_poisson$estimate)
nbinom_probs <- dnbinom(x_vals, size = fit_nbinom$estimate["size"],
                        mu = fit_nbinom$estimate["mu"])
geom_probs <- dgeom(x_vals - 1, prob = fit_geom$estimate)

# Create visualization comparing all distributions
plot_data <- data.frame(
  x = rep(x_vals, 4),
  probability = c(obs_probs, pois_probs, nbinom_probs, geom_probs),
  distribution = rep(c("Empirical", "Poisson", "Neg. Binomial", "Geometric"), 
                     each = length(x_vals))
)

dist_comparison_plot <- ggplot(plot_data, aes(x = x, y = probability, color = distribution)) +
  geom_point(size = 3) + geom_line() +
  labs(title = "Empirical vs Fitted Discrete Distributions",
       x = "Claim Frequency", y = "Probability") +
  theme_minimal() +
  scale_color_manual(values = c("black", "blue", "red", "darkgreen"))
print(dist_comparison_plot)

# ==============================================================================
# MONTE CARLO ESTIMATION WITH VARIANCE REDUCTION
# ==============================================================================

set.seed(2025)
observed_mean <- mean(data$CLM_FREQ)

# Standard Monte Carlo estimation function
plot_monte_carlo <- function(max_simulations = 50000, step = 1000, sample_size = 100) {
  results <- data.frame(Simulation = integer(), Estimator = numeric(), 
                        Lower_CI = numeric(), Upper_CI = numeric())
  
  for (num_sims in seq(100, max_simulations, by = step)) {
    # Generate simulations and calculate estimator
    sims <- replicate(num_sims, mean(rpois(sample_size, lambda = observed_mean)))
    
    # Calculate confidence intervals
    mean_est <- mean(sims)
    sd_est <- sd(sims)
    alpha <- 0.05
    z_value <- qnorm(1 - alpha/2)
    lower_ci <- mean_est - z_value * (sd_est / sqrt(num_sims))
    upper_ci <- mean_est + z_value * (sd_est / sqrt(num_sims))
    
    results <- rbind(results, data.frame(
      Simulation = num_sims, Estimator = mean_est,
      Lower_CI = lower_ci, Upper_CI = upper_ci
    ))
  }
  
  return(list(results = results, final_variance = var(sims), final_estimate = mean(sims)))
}

# Run standard Monte Carlo
mc_standard <- plot_monte_carlo()

# Visualize convergence
mc_plot <- ggplot(mc_standard$results, aes(x = Simulation)) +
  geom_line(aes(y = Estimator), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, fill = "red") +
  geom_hline(yintercept = observed_mean, linetype = "dashed", 
             color = "darkgreen", size = 1) +
  labs(title = "Monte Carlo Convergence: Claim Frequency Estimation",
       subtitle = paste("True mean =", round(observed_mean, 3)),
       x = "Number of Simulations", y = "Estimated Mean") +
  theme_minimal()
print(mc_plot)

# Antithetic variates for variance reduction
set.seed(2025)
sample_size <- 100
num_simulations <- 50000

# Generate antithetic variates using inverse CDF method
anthitetic_sim <- replicate(num_simulations, {
  u <- runif(sample_size)
  x1 <- qpois(u, lambda = observed_mean)        # Original sample
  x2 <- qpois(1 - u, lambda = observed_mean)    # Antithetic sample
  mean(c(x1, x2))
})

# Compare variance reduction effectiveness
variance_comparison <- data.frame(
  Method = c("Standard MC", "Antithetic Variates"),
  Estimate = c(mc_standard$final_estimate, mean(anthitetic_sim)),
  Variance = c(mc_standard$final_variance, var(anthitetic_sim))
)

print("Variance Reduction Comparison:")
print(variance_comparison)

# ==============================================================================
# CLAIM SEVERITY ANALYSIS
# ==============================================================================

# Extract individual claim amounts (severity data)
selected_columns <- data[, paste0("CLM_AMT_", 1:8)]
severity_matrix <- as.matrix(selected_columns)
severity_vector <- as.vector(severity_matrix)
severity_vector <- severity_vector[severity_vector > 0]  # Remove zero claims

# Visualize severity distribution
severity_hist <- ggplot(data.frame(severity = severity_vector), aes(x = severity)) +
  geom_histogram(bins = 50, fill = "lightcoral", alpha = 0.7) +
  labs(title = "Distribution of Claim Severities",
       x = "Claim Amount", y = "Frequency") +
  theme_minimal()
print(severity_hist)

# Fit continuous distributions to claim severity
fit_exp <- fitdist(severity_vector, "exp", method = "mle")
fit_gamma <- fitdist(severity_vector, "gamma", method = "mle")
fit_lnorm <- fitdist(severity_vector, "lnorm", method = "mle")
fit_norm <- fitdist(severity_vector, "norm", method = "mle")

# Density comparison plot
denscomp(
  list(fit_exp, fit_gamma, fit_lnorm, fit_norm),
  legendtext = c("Exponential", "Gamma", "Lognormal", "Normal"),
  main = "Density Comparison: Claim Severity Models",
  fitlwd = 2
)

# Goodness-of-fit comparison using AIC and Kolmogorov-Smirnov test
gof_severity <- gofstat(
  list(fit_exp, fit_gamma, fit_lnorm, fit_norm),
  fitnames = c("Exponential", "Gamma", "Lognormal", "Normal")
)

severity_results <- data.frame(
  AIC = c(fit_exp$aic, fit_gamma$aic, fit_lnorm$aic, fit_norm$aic),
  KS_Statistic = gof_severity$ks
) %>% arrange(KS_Statistic)

print("Model Comparison for Claim Severity:")
print(severity_results)
# Result: Normal distribution provides the best fit based on KS test

# Stratified sampling for variance reduction
set.seed(2025)
mean_norm <- fit_norm$estimate["mean"]
sd_norm <- fit_norm$estimate["sd"]

sample_size <- 100
num_simulations <- 10000

# Standard Monte Carlo for severity
sims_standard <- replicate(num_simulations, {
  x <- rnorm(sample_size, mean = mean_norm, sd = sd_norm)
  mean(x)
})

# Stratified sampling - divide [0,1] into equal intervals
sims_stratified <- replicate(num_simulations, {
  u <- (1:sample_size - runif(sample_size)) / sample_size
  x <- qnorm(u, mean = mean_norm, sd = sd_norm)
  mean(x)
})

# Compare results
severity_mc_results <- data.frame(
  Method = c("Standard MC", "Stratified Sampling"),
  Estimate = c(mean(sims_standard), mean(sims_stratified)),
  Variance = c(var(sims_standard), var(sims_stratified))
)

print("Severity Estimation: Variance Reduction Results:")
print(severity_mc_results)

# ==============================================================================
# RISK PREMIUM CALCULATION
# ==============================================================================

# Compound Poisson model: S = X₁ + X₂ + ... + Xₙ where N ~ Poisson(λ)
# N represents claim frequency, X_i represents claim severity
set.seed(2025)
n_simulations <- 10000
risk_premiums <- numeric(n_simulations)

# Simulate aggregate claims for each policy
for (i in 1:n_simulations) {
  n_claims <- rpois(1, lambda = observed_mean)  # Draw number of claims
  if (n_claims > 0) {
    # Draw claim amounts and sum them
    claim_amounts <- rnorm(n_claims, mean = mean_norm, sd = sd_norm)
    risk_premiums[i] <- sum(claim_amounts)
  } else {
    risk_premiums[i] <- 0  # No claims
  }
}

# Visualize risk premium distribution
risk_premium_plot <- ggplot(data.frame(premium = risk_premiums), aes(x = premium)) +
  geom_histogram(bins = 50, fill = "lightgreen", alpha = 0.7) +
  labs(title = "Simulated Risk Premium Distribution",
       x = "Total Claim Amount", y = "Frequency") +
  theme_minimal()
print(risk_premium_plot)

# Compare simulated vs empirical risk premiums
risk_comparison <- data.frame(
  Metric = c("Simulated Risk Premium", "Empirical Risk Premium", "Average Premium Paid"),
  Value = c(mean(risk_premiums), mean(data$TOTAL_CLM_AMT), mean(data$PREMIUM))
)

print("Risk Premium Comparison:")
print(risk_comparison)

# ==============================================================================
# RISK MANAGEMENT METRICS
# ==============================================================================

# Value-at-Risk (VaR) at 99.5% confidence level
confidence_level <- 0.995
var_99_5 <- quantile(risk_premiums, confidence_level)

cat("99.5% Value-at-Risk:", round(var_99_5, 2), "\n")
cat("This represents the economic capital needed to cover claims 99.5% of the time.\n")
cat("Probability of ruin:", (1 - confidence_level) * 100, "%\n")

# Visualize VaR
var_plot <- ggplot(data.frame(premium = risk_premiums), aes(x = premium)) +
  geom_histogram(bins = 50, fill = "lightblue", alpha = 0.7) +
  geom_vline(xintercept = var_99_5, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Risk Premium Distribution with 99.5% VaR",
       x = "Total Claim Amount", y = "Frequency",
       subtitle = paste("VaR₉₉.₅% =", round(var_99_5, 2))) +
  theme_minimal()
print(var_plot)

# Combined Ratio Analysis
combined_ratio <- mean(data$TOTAL_CLM_AMT) / mean(data$PREMIUM)

cat("Combined Ratio:", round(combined_ratio, 3), "\n")

if (combined_ratio > 1) {
  cat("Result: Reinsurance is recommended (Combined Ratio > 1)\n")
  cat("Claims exceed premiums by", round((combined_ratio - 1) * 100, 1), "%\n")
} else {
  cat("Result: No reinsurance needed (Combined Ratio ≤ 1)\n")
}

# Reinsurance threshold analysis including administrative costs
admin_cost_rate <- 0.10  # 10% administrative costs
admin_costs <- admin_cost_rate * mean(data$PREMIUM)

sorted_claims <- sort(data$TOTAL_CLM_AMT)
combined_ratios <- (sorted_claims + admin_costs) / mean(data$PREMIUM)

# Find the first claim amount where combined ratio exceeds 1
reinsurance_threshold <- sorted_claims[which(combined_ratios > 1)[1]]

cat("\nRecommended reinsurance threshold (including 10% admin costs):", 
    round(reinsurance_threshold, 2), "\n")
