# THEORETICAL DISTRIBUTION ESTIMATION
# The proportion of most diseases caused by alcohol in
# the component cause model in a population is determined by:
#  • The distribution of the volume of exposure
#  • The relative risk associated with each level of exposure
rm(list = ls())
gc()
# List of required packages
required_packages <- c("dplyr", "readr", "haven", 
                       "ggplot2", "gridExtra", "fitdistrplus")

# Install and load required packages
sapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

data <- readRDS("enpg_full.RDS")

cd_collapsed <- data %>% 
  filter(volajohdia > 0) %>% 
  mutate(volajohdia = ifelse(volajohdia > 150, 150, volajohdia)) %>% 
  pull(volajohdia)

# Fit log-normal distribution
fit_lognorm_col <- fitdist(cd_collapsed, "lnorm")
summary(cd_collapsed)

x_vals <- seq(0, 150, length.out = 1000)
y_lnorm <- dlnorm(x_vals, meanlog = fit_lognorm_col$estimate["meanlog"],
                  sdlog = fit_lognorm_col$estimate["sdlog"])
plot_lnorm <- data.frame(x = x_vals, y = y_lnorm)

p_lnorm <- ggplot(plot_lnorm, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = "Fitted lnorm Distribution",
       x = "Alcohol Consumption (grams per day)",
       y = "Density") +
  theme_minimal()
# Fit gamma distribution
fit_gamma_col <- fitdist(cd_collapsed, "gamma")

y_gamma <- dgamma(x_vals, shape = fit_gamma_col$estimate["shape"], rate = fit_gamma_col$estimate["rate"])

# Plot the fitted gamma distribution
plot_gamma <- data.frame(x = x_vals, y = y_gamma)

p_gamma <- ggplot(plot_gamma, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = "Fitted Gamma Distribution",
       x = "Alcohol Consumption (grams per day)",
       y = "Density") +
  theme_minimal()

# Fit Weibull distribution
fit_weibull_col <- fitdist(cd_collapsed, "weibull")

y_weibull <- dweibull(x_vals, shape = fit_weibull_col$estimate["shape"], scale = fit_weibull_col$estimate["scale"])

# Plot the fitted gamma distribution
plot_weibull <- data.frame(x = x_vals, y = y_weibull)

p_weibull <- ggplot(plot_weibull, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = "Fitted Weibull Distribution",
       x = "Alcohol Consumption (grams per day)",
       y = "Density") +
  theme_minimal()

grid.arrange(p_lnorm, p_gamma, p_weibull, nrow = 3)

# Q-Q Plot function
create_qqcomp_ggplot <- function(fit_list, title) {
  ggplot_fit <- qqcomp(fit_list, plotstyle = "ggplot")
  ggplot_fit + 
    ggtitle(title) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    labs(color = "Distribution")+
    coord_cartesian(xlim = c(0, 150))
}

# Generate Q-Q plots for collapsed data
qqplot_collapsed <- create_qqcomp_ggplot(list(fit_lognorm_col, fit_gamma_col, fit_weibull_col), "Q-Q Plot Comparison for Collapsed Data")

# Generate Q-Q plots for filtered data
qqplot_filtered <- create_qqcomp_ggplot(list(fit_lognorm_fil, fit_gamma_fil, fit_weibull_fil), "Q-Q Plot Comparison for Filtered Data")

# Combine plots in a grid
grid.arrange(qqplot_collapsed, qqplot_filtered, nrow = 2)

# PAF FOR TUBERCULOSIS

calculate_paf <- function(beta_1, data, distribution) {
  # Fit the chosen distribution
  if (distribution == "lognorm") {
    fit <- fitdist(data, "lnorm")
    prevalence_function <- function(x) {
      dlnorm(x, meanlog = fit$estimate["meanlog"], sdlog = fit$estimate["sdlog"])
    }
  } else if (distribution == "gamma") {
    fit <- fitdist(data, "gamma")
    prevalence_function <- function(x) {
      dgamma(x, shape = fit$estimate["shape"], rate = fit$estimate["rate"])
    }
  } else if (distribution == "weibull") {
    fit <- fitdist(data, "weibull")
    prevalence_function <- function(x) {
      dweibull(x, shape = fit$estimate["shape"], scale = fit$estimate["scale"])
    }
  } else {
    stop("Unsupported distribution")
  }
  
  rr_function <- function(x) {
    exp(beta_1 * x)
  }
  
  integrand <- function(x) {
    prevalence_function(x) * (rr_function(x) - 1)
  }
  
  numerator <- integrate(integrand, lower = 0, upper = 150)$value
  denominator <- numerator + 1
  
  paf <- numerator / denominator
  return(paf)
}

# Calculate the PAF for TB using different distributions
beta_1 <- 0.0179695
paf_tb_lognorm <- calculate_paf(beta_1, cd_collapsed, "lognorm")
paf_tb_gamma <- calculate_paf(beta_1, cd_collapsed, "gamma")
paf_tb_weibull <- calculate_paf(beta_1, cd_collapsed, "weibull")

list(
  lognorm = paf_tb_lognorm,
  gamma = paf_tb_gamma,
  weibull = paf_tb_weibull
)


# PLOT
# Define the beta_1 for the relative risk function
beta_1 <- 0.0179695

# Define the relative risk function
rr_function <- function(x) {
  exp(beta_1 * x)
}

# Define the prevalence function using the fitted gamma distribution
prevalence_function <- function(x) {
  dgamma(x, shape = fit_gamma_col$estimate["shape"], rate = fit_gamma_col$estimate["rate"])
}

# Define the integrand function
integrand <- function(x) {
  prevalence_function(x) * (rr_function(x) - 1)
}

# Generate values for plotting the integrand function
x_vals <- seq(0, 150, length.out = 1000)
y_vals <- integrand(cd_collapsed)

# Plot the integrand function
plot_integrand <- data.frame(x = cd_collapsed, y = y_vals)
p_integrand <-  ggplot(plot_integrand, aes(x = x, y = y_vals)) +
  geom_line(color = "blue") +
  labs(title = "Integrand Function Using Gamma Distribution",
       x = "Alcohol Consumption (grams per day)",
       y = "Prevalence x (RR-1)") +
  theme_minimal() 

# PLOT THE RISK FUNCTION
rr <- rr_function(cd_collapsed)
plot_rr <- data.frame(x = cd_collapsed, y = rr)
p_rr <- ggplot(plot_rr, aes(x = x, y = rr)) +
  geom_line(color = "blue") +
  labs(title = "Relative Risk of Tuberculosis",
       x = "Alcohol Consumption (grams per day)",
       y = "Relative Risk (Tuberculosis)") +
  theme_minimal() +
  geom_hline(yintercept = 1) 

grid.arrange(p_gamma, p_rr, p_integrand, nrow = 3)
# REVISAR LA FUNCIÓN DE PREVALENCIA



integrand <- function(x) {
  prevalence_function(x) * (rr_function(x) - 1)
}

numerator <- integrate(y_vals)$value
denominator <- numerator + 1

paf <- numerator / denominator
return(paf)

