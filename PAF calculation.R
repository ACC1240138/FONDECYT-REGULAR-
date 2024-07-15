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

data <- readRDS("enpg_full.RDS") %>% 
  mutate(aux = ifelse(oh1 == "No" & !is.na(oh2) ,1,0)) %>% 
  filter(aux == 0) %>% 
  dplyr::select(-aux)

cd_collapsed <- data %>% 
  filter(volajohdia > 0) %>% 
  mutate(volajohdia = ifelse(volajohdia > 150, 150, volajohdia)) %>% 
  pull(volajohdia)

# Fit log-normal distribution
fit_lognorm_col <- fitdist(cd_collapsed, "lnorm")
x_vals <- seq(0, 150, length.out = 1000)
y_lnorm <- dlnorm(x_vals, meanlog = fit_lognorm_col$estimate["meanlog"],
                  sdlog = fit_lognorm_col$estimate["sdlog"])

p_lnorm <- ggplot(data.frame(x = x_vals, y = y_lnorm), aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = "Fitted lnorm Distribution",
       x = "Alcohol Consumption (grams per day)",
       y = "Density") +
  theme_minimal()

# Fit gamma distribution
fit_gamma_col <- fitdist(cd_collapsed, "gamma")

y_gamma <- dgamma(x_vals, shape = fit_gamma_col$estimate["shape"], 
                  rate = fit_gamma_col$estimate["rate"])

p_gamma <- ggplot(data.frame(x = x_vals, y = y_gamma), aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = "Fitted Gamma Distribution",
       x = "Alcohol Consumption (grams per day)",
       y = "Density") +
  theme_minimal()

# Fit Weibull distribution
fit_weibull_col <- fitdist(cd_collapsed, "weibull")

y_weibull <- dweibull(x_vals, shape = fit_weibull_col$estimate["shape"], scale = fit_weibull_col$estimate["scale"])

p_weibull <- ggplot(data.frame(x = x_vals, y = y_weibull), aes(x = x, y = y)) +
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

# PAF FOR TUBERCULOSIS
p_nd <- data %>% 
  mutate(non_drinkers = ifelse(oh1 == "No", 1, 0),
         former_drinkers = ifelse(oh2 == ">30" | oh2 == ">1 año", 1 , 0)) %>% 
  summarise(p_nd = mean(non_drinkers, na.rm = T)) %>% 
  pull(p_nd)

p_fd <- data %>% 
  mutate(non_drinkers = ifelse(oh1 == "No", 1, 0),
         former_drinkers = ifelse(oh1 == "Si" & (oh2 == ">30" | oh2 == ">1 año"), 1 , 0)) %>% 
  summarise(p_fd = mean(former_drinkers, na.rm = T)) %>% 
  pull(p_fd)

# FUNCTION TO CALCULATE PAF USING CURRENT DRINKERS INFORMATION ONLY
calculate_paf_cr <- function(beta_1, data, distribution) {

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
paf_tb_lognorm_cr <- calculate_paf_cr(beta_1, cd_collapsed, "lognorm")
paf_tb_gamma_cr <- calculate_paf_cr(beta_1, cd_collapsed, "gamma")
paf_tb_weibull_cr <- calculate_paf_cr(beta_1, cd_collapsed, "weibull")

list(
  lognorm = paf_tb_lognorm_cr,
  gamma = paf_tb_gamma_cr,
  weibull = paf_tb_weibull_cr
)

# FUNCTION TO CALCULATE PAF USING PROPORTION OF FORMER DRINKERS 
 
calculate_paf_fd <- function(beta_1, var,data, distribution, p_fd, rr_fd) {
  
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
  
  integral_value <- integrate(integrand, lower = 0, upper = 150)$value
  
  numerator <- p_fd*(rr_fd-1)+integral_value
  denominator <- p_fd*(rr_fd-1)+(integral_value + 1)
  
  paf <- numerator / denominator

}

set.seed(123)
num_simulations <- 1000
simulated_pafs <- numeric(num_simulations)

for (i in 1:num_simulations) {
  beta_sim <- rnorm(1, mean = beta_1, sd = sqrt(0.0072152**2))
  simulated_pafs[i] <- calculate_paf_fd(beta_sim, var= 0.0072152**2, cd_collapsed, "lognorm", p_fd, 1)
}

# Calculate 95% Uncertainty Intervals
paf_mean <- mean(simulated_pafs)
lower_ui <- quantile(simulated_pafs, 0.025)
upper_ui <- quantile(simulated_pafs, 0.975)

list(
  mean_paf = paf_mean,
  lower_ui = lower_ui,
  upper_ui = upper_ui
)
# Calculate the PAF for TB using different distributions
paf_tb_lognorm_fd <- calculate_paf_fd(beta_1, var= 0.0072152**2,cd_collapsed,"lognorm", p_fd, 1)
paf_tb_gamma_fd <- calculate_paf_fd(beta_1, var= 0.0072152**2, cd_collapsed, "gamma", p_fd, 1)
paf_tb_weibull_fd <- calculate_paf_fd(beta_1, var= 0.0072152**2, cd_collapsed, "weibull", p_fd, 1)

list(
  lognorm = paf_tb_lognorm_nd,
  gamma = paf_tb_gamma_nd,
  weibull = paf_tb_weibull_nd
)

# Calculate the PAF for pancreatitis using different distributions
paf_pc_lognorm_nd <- calculate_paf_nd(0.0173451 , cd_collapsed, "lognorm")
paf_pc_gamma_nd <- calculate_paf_nd(0.0173451 , cd_collapsed, "gamma")
paf_pc_weibull_nd <- calculate_paf_nd(0.0173451 , cd_collapsed, "weibull")

paf_pc_lognorm_cr <- calculate_paf_cr(0.0173451 , cd_collapsed, "lognorm")
paf_pc_gamma_cr <- calculate_paf_cr(0.0173451 , cd_collapsed, "gamma")
paf_pc_weibull_cr <- calculate_paf_cr(0.0173451 , cd_collapsed, "weibull")


















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
y_vals <- integrate(integrand, 0, 150)$value
(p_nd*p_fd*y_vals)/(p_nd*p_fd*(y_vals+1))

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


