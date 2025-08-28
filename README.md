# Insurance Claims Stochastic Modeling

## Repository Structure

This repository contains a comprehensive stochastic modeling analysis of insurance claims data. The project is organized as follows:

- **`claims_sim.R`** - Complete R script containing all analysis functions and methods
- **`doc/claims_simulation.html`** - Interactive HTML report with code, visualizations, and detailed commentary ([View Report](https://cesco2002.github.io/Insurance-Claims-Simulation/claims_simulation.html))
- **`data/DATA_SET_4.csv`** - Source dataset (motor insurance portfolio data)
- **`README.md`** - This file

**For optimal viewing experience:** While the R code demonstrates the technical implementation, I recommend viewing the [HTML report](https://cesco2002.github.io/Insurance-Claims-Simulation/claims_simulation.html) for a complete analysis with integrated results, visualizations, and interpretations.

## Project Overview

This analysis addresses a real-world actuarial problem: determining whether an insurance company can reduce premium rates while maintaining profitability and regulatory compliance. The project combines advanced statistical modeling with Monte Carlo simulation techniques to assess insurance risk.

### Business Context

An insurance company's marketing department requested a tariff reduction analysis for their motor insurance portfolio. As the actuarial analyst, I needed to evaluate whether current pricing strategies align with the underlying risk profile of the portfolio, considering regulatory solvency requirements and profitability constraints.

## Methodology & Implementation

### 1. Data Preparation & Exploration
- **Dataset**: 1,002 motor insurance policies with claim frequency, severity, and premium information
- **Data Cleaning**: standardized numerical formats, handled outliers

### 2. Frequency Modeling
**Objective**: Model the number of claims per policyholder

- **Distribution Candidates**: Poisson, Negative Binomial, Geometric
- **Model Selection**: Maximum Likelihood Estimation (MLE) with AIC and Chi-squared criteria
- **Result**: Poisson distribution (λ = 3.094) provided optimal fit
- **Validation**: Monte Carlo simulation with 50,000 iterations

### 3. Severity Modeling  
**Objective**: Model individual claim amounts

- **Distribution Candidates**: Exponential, Gamma, Log-Normal, Normal
- **Model Selection**: MLE with AIC and Kolmogorov-Smirnov statistics
- **Result**: Normal distribution (μ = 598.72, σ = 47.33) best captured severity patterns
- **Validation**: Monte Carlo estimation with empirical comparison

### 4. Variance Reduction Techniques
**Advanced Simulation Methods**:

- **Antithetic Variates Method**: Applied to frequency modeling
  - Reduced variance from 0.0311 to 0.00098 (96.8% reduction)
  - Uses complementary uniform random variables (U, 1-U) for Poisson generation

- **Stratified Sampling Method**: Applied to severity modeling  
  - Partitions probability space for improved coverage
  - Significant variance reduction while maintaining unbiased estimates

### 5. Risk Premium Analysis
**Compound Poisson Model**: 
```
S = Σ(i=1 to N) X_i
```
Where:
- N ~ Poisson(3.094) represents claim frequency
- X_i ~ Normal(598.72, 47.33²) represents claim severity

**Key Findings**:
- Simulated risk premium: €1,857.53
- Empirical risk premium: €1,852.50  
- Average premium paid: €1,688.25
- **Gap**: Customers paying 9.1% below expected risk cost

### 6. Risk Management Metrics

**Value-at-Risk (99.5% confidence)**:
- VaR₉₉.₅% = €5,317.72 per policy
- Economic capital requirement for 0.5% probability of ruin
- Compliant with Solvency II regulatory standards

**Combined Ratio Analysis**:
- Combined ratio: 1.097 (claims exceed premiums by 9.7%)
- Threshold: Combined ratio > 1.0 indicates need for reinsurance
- Recommended reinsurance threshold: €1,578 (including 10% administrative costs)


---
