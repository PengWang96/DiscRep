# README

DiscRep provides tools for assessing replicability using the distinguishability criterion and posterior-PRP. It supports random and fixed effect models, tests for heterogeneity and publication bias, and tissue-specific eQTL analysis.

## System Requirements

- Software Dependencies:
  - R (>= v3.5.0).
  - R packages: stats (>= v4.4.3, base R), devtools (>= v2.4.5, for installation only).

- Tested Environments:
  - Windows 11 (22H2), 13th Gen Intel Core i9-13900HX, 16GB RAM.

- Hardware: Standard desktop/laptop (no specialized hardware required).


## Installation

To install this R package from GitHub, you'll need to use the `devtools` package. If you haven't already installed `devtools`, you can do so by running:

```r
install.packages("devtools")
```

Once `devtools` is installed, you can install the `DiscRep` package directly from the GitHub repository:

```r
devtools::install_github("PengWang96/DiscRep")
```

The installation time is typically within 5 seconds (depending on your network speed). After installation, you can load the package into your R environment using:

```r
library(DiscRep)
```

## Data

The package provides three datasets for illustrative purposes: `dat.ha`, `dat.lcj`, and `dat.slf`. Each dataset represents meta-analytic data.

- **`dat.ha`**: A data frame of 109 studies investigating the effect of placebo interventions on patient-reported outcomes across clinical conditions.
  - `y`: Observed effect size for each study.
  - `s2`: Within-study variance for each study.
  - **Source**: Hróbjartsson A, Gøtzsche PC (2010). *Cochrane Database of Systematic Reviews*, **1**. Art. No.: CD003974. [doi:10.1002/14651858.CD003974.pub3](https://doi.org/10.1002/14651858.CD003974.pub3).

- **`dat.lcj`**: A data frame of 33 studies examining the effect of progressive resistance strength training on physical function in older adults.
  - `y`: Observed effect size for each study.
  - `s2`: Within-study variance for each study.
  - **Source**: Liu CJ, Latham NK (2009). *Cochrane Database of Systematic Reviews*, **3**. Art. No.: CD002759. [doi:10.1002/14651858.CD002759.pub2](https://doi.org/10.1002/14651858.CD002759.pub2).

- **`dat.slf`**: A data frame of 56 studies exploring the effect of nicotine gum on smoking cessation.
  - `y`: Observed effect size for each study.
  - `s2`: Within-study variance for each study.
  - **Source**: Stead LF et al. (2012). *Cochrane Database of Systematic Reviews*, **11**. Art. No.: CD000146. [doi:10.1002/14651858.CD000146.pub4](https://doi.org/10.1002/14651858.CD000146.pub4).

### Usage

Load these datasets as follows:

```r
data("dat.ha")
data("dat.lcj")
data("dat.slf")
```



## Functions

### 1. Get posterior-predictive replication p-values (posterior-PRPs) under Random Effect Model: `metropolis_hastings`

The `metropolis_hastings` function fits a random effect model using the Metropolis-Hastings MCMC algorithm. You can specify different test statistics and heterogeneity levels.

#### Usage Example:

```r
# Load necessary data
data("dat.slf")
m <- nrow(dat.slf)
hat_beta <- dat.slf$y
hat_sigma_sq <- dat.slf$s2

# Optional: precompute k_vec once and reuse
# (recommended for large simulations to avoid repeated setup cost).
pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05)
k_vec <- sapply(pvec, inverse_P_mis)

# Fit model using Metropolis-Hastings algorithm (Egger test)
results_random_egger <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq, test = "Egger", k_vec = k_vec)

# Print results
print(results_random_egger)
```

Expected output is `results_random_egger$p_value`, which is posterior-PRP. Posterior-PRP < 0.05 supports irreplicability. 

The program run time for repeating this example 2000 times using parallel computation in R on the previously provided computer specifications is 20 minutes.

#### Additional Examples:
- Using the `Q` test:
  
```r
results_random_q <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq, test = "Q", k_vec = k_vec)
```

- Specifying a fixed tolerable heterogeneity level (e.g., `k = 0.2726814` corresponds to $P_{mis} = 0.05$):

```r
results_random_fixed_k <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq, k_vec = 0.2726814)
```

#### Note

In `metropolis_hastings`, `k_vec` is optional.
If omitted, the function uses default
`pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05); k_vec <- sapply(pvec, inverse_P_mis)`
and caches it in the current R process.
For large simulation loops, precompute `k_vec` once outside the loop and pass it in.
The function warns if any `k_vec > 0.2726814`, because only `k` in `[0, 0.2726814]` corresponds to `P_mis < 0.05`.


### 2. Get posterior-PRPs under Fixed Effect Model: `fixed_effect`

The `fixed_effect` function fits a fixed effect model using a Monte Carlo simulation.

#### Usage Example:

```r
# Fit the fixed effect model
results_fixed <- fixed_effect(10000, m, hat_beta, hat_sigma_sq)

# Print the results
print(results_fixed)
```

### 3. Cochran's Q Test: `frequency_pvalue`

The `frequency_pvalue` function calculates the p-value for Cochran's Q test in the fixed effect model. This is useful for assessing heterogeneity across studies.

#### Usage Example:

```r
p_value_q <- frequency_pvalue(hat_beta, hat_sigma_sq)
print(p_value_q)
```

### 4. Numerical Calculation of Posterior-PRP: `calc_PRP_FE_NumInt`

This function calculates the Posterior-PRP (Posterior Predictive Replication Probability) using numerical integration under the fixed effect model. 
It integrates over the posterior distribution of true effect sizes.

#### Usage Example:

```r
prp_result <- calc_PRP_FE_NumInt(hat_beta, hat_sigma_sq)
print(prp_result)
```

## Example Workflow

Here's a complete workflow using the provided functions:

```r
# Load necessary data
data("dat.slf")

# Extract values
m <- nrow(dat.slf)
hat_beta <- dat.slf$y
hat_sigma_sq <- dat.slf$s2

# Optional: precompute k_vec once and reuse
# (recommended for large simulations to avoid repeated setup cost).
pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05)
k_vec <- sapply(pvec, inverse_P_mis)

# Random Effect Model (Egger Test)
results_random_egger <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq, test = "Egger", k_vec = k_vec)
print(results_random_egger)

# Random Effect Model (Q Test)
results_random_q <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq, test = "Q", k_vec = k_vec)
print(results_random_q)

# Random Effect Model with fixed heterogeneity level
results_random_fixed_k <- metropolis_hastings(10000, 0.05, m, hat_beta, hat_sigma_sq, k_vec = 0.2726814)
print(results_random_fixed_k)

# Fixed Effect Model
results_fixed <- fixed_effect(10000, m, hat_beta, hat_sigma_sq)
print(results_fixed)

# Cochran's Q Test
p_value_q <- frequency_pvalue(hat_beta, hat_sigma_sq)
print(p_value_q)

# Numerical calculation of Posterior-PRP
prp_result <- calc_PRP_FE_NumInt(hat_beta, hat_sigma_sq)
print(prp_result)
```


## Data and reproducibility

All the R scripts to reproduce the numeric simulation results and real data analysis in 
the manuscript can be found in “paper_repro�?folder.


# References

The datasets used in this package are based on results from systematic reviews. Below are the sources for the data:

- Hróbjartsson A, Gøtzsche PC (2010). "Placebo interventions for all clinical conditions." *Cochrane Database of Systematic Reviews*, **1**. Art. No.: CD003974. [https://doi.org/10.1002/14651858.CD003974.pub3](https://doi.org/10.1002/14651858.CD003974.pub3)

- Liu CJ, Latham NK (2009). "Progressive resistance strength training for improving physical function in older adults." *Cochrane Database of Systematic Reviews*, **3**. Art. No.: CD002759. [https://doi.org/10.1002/14651858.CD002759.pub2](https://doi.org/10.1002/14651858.CD002759.pub2)

- Stead LF, Perera R, Bullen C, Mant D, Hartmann-Boyce J, Cahill K, Lancaster T (2012). "Nicotine replacement therapy for smoking cessation." *Cochrane Database of Systematic Reviews*, **11**. Art. No.: CD000146. [https://doi.org/10.1002/14651858.CD000146.pub4](https://doi.org/10.1002/14651858.CD000146.pub4)

---

These references can be cited when using the datasets included in this package, as they originate from these systematic reviews.


## Contributing

Feel free to contribute to the project by reporting issues, suggesting features, or submitting pull requests.

## License

This project is licensed under the GPL-3 License.

---

Enjoy using the package and explore the various statistical methods provided! For questions or feedback, please open an issue on GitHub.

---

To view the full repository and more detailed instructions, visit the [DiscRep GitHub repository](https://github.com/PengWang96/DiscRep).
