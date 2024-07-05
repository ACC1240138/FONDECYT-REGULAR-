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

cd_filter <- data %>% 
  filter(volajohdia > 0 & volajohdia <= 150) %>% 
  pull(volajohdia)

# Fit log-normal distribution
fit_lognorm_col <- fitdist(cd_collapsed, "lnorm")
fit_lognorm_fil <- fitdist(cd_filter, "lnorm")

# Fit gamma distribution
fit_gamma_col <- fitdist(cd_collapsed, "gamma")
fit_gamma_fil <- fitdist(cd_filter, "gamma")

# Fit Weibull distribution
fit_weibull_col <- fitdist(cd_collapsed, "weibull")
fit_weibull_fil <- fitdist(cd_filter, "weibull")



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
paf_tb_lognorm <- calculate_paf(beta_1, cd_filter, "lognorm")
paf_tb_gamma <- calculate_paf(beta_1, cd_filter, "gamma")
paf_tb_weibull <- calculate_paf(beta_1, cd_filter, "weibull")

list(
  lognorm = paf_tb_lognorm,
  gamma = paf_tb_gamma,
  weibull = paf_tb_weibull
)

calculate_paf(beta_1, cd_collapsed, "lognorm")
calculate_paf(beta_1, cd_collapsed, "gamma")
calculate_paf(beta_1, cd_collapsed, "weibull")



