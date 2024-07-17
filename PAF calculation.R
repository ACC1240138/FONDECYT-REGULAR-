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

cd_vector <- data %>% 
  filter(volajohdia > 0) %>% 
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



# FUNCTION TO CALCULATE PAF USING PROPORTION OF FORMER DRINKERS 
# THIS FUNCTION APPLIES FOR ANY DISEASE WITH A RR GIVEN BY exp(B1*X)
# TUBERCULOSIS, PANCREATITIS, COLON AND RECTUM CANCER, LIVER CANCER, 
# BREAST CANCER
 
calculate_paf_fd <- function(beta_1, data, distribution, p_fd, rr_fd) {
  
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

# BETAS OF RELATIVE RISK

# TUBERCULOSIS
b_tb <- 0.0179695
var_tb <- 0.0072152^2

# HIV/AIDS (Male)
b_hiv_male <- ln(1.54)
# for consumption > 61 grams/day
var_hiv_male <- 0.0782107722^2

# HIV/AIDS (Female)
b_hiv_fem <- ln(1.54)
# for consumption > 49 grams/day
var_hiv_fem <- 0.0782107722^2

# LOWER RESPORATORY INFECTIONS
b_lri <- 0.4764038
var_lri <- 0.19220552^2

# LIP AND ORAL CAVITY CANCER
b1_locan <- 0.02474
b2_locan <- -0.00004
# Variance-covariance matrix:
# Variance (b1): 0.000002953
# Variance (b2): 0.000000000102
# Covariance: -0.0000000127

# Other Pharyngeal Cancers
b1_opcan <- 0.02474
b2_opcan <- -0.00004
# Variance-covariance matrix:
# Variance (b1): 0.000002953
# Variance (b2): 0.000000000102
# Covariance: -0.0000000127

# Oesophagus Cancer
b1_oescan <- 0.0132063596418668
b2_oescan <- -4.14801974664481*10^-08
# Variance-covariance matrix:
# Variance (b1): 1.525706*10^-07
# Variance (b2): 0.000002953
# Covariance: -6.885205*10^-13

# Colon and Rectum Cancers (Male)
b1_crcan_male <- 0.006279
var_crcan <- 0.000000907
rr_crcan_fd_male <- 2.19
var_rr_crcan_fd_male <- 0.04651062

# Colon and Rectum Cancers (Female)
b1_crcan_fem <- 0.006279
rr_crcan_fd_fem <- 1.05
var_rr_crcan_fd_fem <- 0.1459680025873172

# Liver Cancer (Male)
b1_lican_male <- 0.005041
var_lican <- 0.000003097
rr_lican_fd_male <- 2.23
var_rr_lican_fd_male <- 0.2590977572

# Liver Cancer (Female)
b1_lican_fem <- 0.005041
rr_lican_fd_fem <- 2.68
var_rr_lican_fd_fem <- 0.2725606092

# Breast Cancer (Female)
b1_bcan <- 0.01018
var_bcan <- 0.0000007208
rr_bcan_fd <- 1
var_rr_bcan_fd <- 0

# Larynx Cancer
b1_lxcan <- 0.01462
b2_lxcan <- -0.00002
# Variance-covariance matrix:
# Variance (b1): 0.000003585
# Variance (b2): 0.000000000126
# Covariance: -0.0000000162

# Diabetes Mellitus (Male)
b1_dm_male <- 0.1763703
var_dm_male <- -0.0728256
# Variance-covariance matrix:
# Variance (b1): 0.16812525
# Variance (b2): 0.29964479
# Covariance: -0.2240129

# Diabetes Mellitus (Female)
b1_dm_fem <- -1.3133910
b2_dm_fem <- 1.0142390
# Variance-covariance matrix:
# Variance (b1): 0.16812525
# Variance (b2): 0.29964479
# Covariance: -0.2240129

# Epilepsy
b1_epi <- 1.22861
var_epi <- 0.13919742

# Hypertensive Heart Disease (Male)
b1_hhd_male <- 0.0150537
b2_hhd_male <- -0.0156155
# Variance-covariance matrix:
# Variance (b1): 0.00241962
# Variance (b2): 0.00416022
# Covariance: -0.000009788

# Hypertensive Heart Disease (Female)
b1_hhd_fem <- -0.0154196
b2_hhd_fem <- 0.0217586
# Variance-covariance matrix:
# Variance (b1): 0.0064592
# Variance (b2): 0.00680182
# Covariance: -0.00004194

# PANCREATITIS MALE
b1_panc <- 0.0173451
var_panc <- 0.003803**2
rr_panc_fd <- 2.2
var_rr_panc_fd <- 0.213**2

# PANCREATITIS FEMALE















# PROPORTION OF NON DRINKERS
p_nd <- data %>% 
  mutate(non_drinkers = ifelse(oh1 == "No", 1, 0),
         former_drinkers = ifelse(oh2 == ">30" | oh2 == ">1 año", 1 , 0)) %>% 
  summarise(p_nd = mean(non_drinkers, na.rm = T)) %>% 
  pull(p_nd)

# PROPORTION OF FORMER DRINKERS
p_fd <- data %>% 
  mutate(non_drinkers = ifelse(oh1 == "No", 1, 0),
         former_drinkers = ifelse(oh1 == "Si" & (oh2 == ">30" | oh2 == ">1 año"), 1 , 0)) %>% 
  summarise(p_fd = mean(former_drinkers, na.rm = T)) %>% 
  pull(p_fd)

#FD FEM
cd_fem <- data %>% 
  filter(sexo == "Mujer", volajohdia > 0) %>% 
  mutate(volajohdia = ifelse(volajohdia > 150, 150, volajohdia)) %>% 
  pull(volajohdia)
p_fd_fem <- data %>% 
  filter(sexo == "Mujer") %>% 
  mutate(non_drinkers = ifelse(oh1 == "No", 1, 0),
         former_drinkers = ifelse(oh1 == "Si" & (oh2 == ">30" | oh2 == ">1 año"), 1 , 0)) %>% 
  summarise(p_fd = mean(former_drinkers, na.rm = T)) %>% 
  pull(p_fd)

#FD MALE
cd_male <- data %>% 
  filter(sexo == "Hombre", volajohdia > 0) %>% 
  mutate(volajohdia = ifelse(volajohdia > 150, 150, volajohdia)) %>% 
  pull(volajohdia)
p_fd_male <- data %>% 
  filter(sexo == "Hombre") %>% 
  mutate(non_drinkers = ifelse(oh1 == "No", 1, 0),
         former_drinkers = ifelse(oh1 == "Si" & (oh2 == ">30" | oh2 == ">1 año"), 1 , 0)) %>% 
  summarise(p_fd = mean(former_drinkers, na.rm = T)) %>% 
  pull(p_fd)

# PAF TUBERCULOSIS MALES
paf_tb_lognorm_fd <- calculate_paf_fd(b_tb,cd_male,"lognorm", p_fd_male, 1)
paf_tb_gamma_fd <- calculate_paf_fd(b_tb, cd_male, "gamma", p_fd_male, 1)
paf_tb_weibull_fd <- calculate_paf_fd(b_tb, cd_male, "weibull", p_fd_male, 1)

list(
  lognorm = paf_tb_lognorm_fd,
  gamma = paf_tb_gamma_fd,
  weibull = paf_tb_weibull_fd
)

# PAF TUBERCULOSIS FEM
paf_tb_lognorm_fem <- calculate_paf_fd(b_tb,cd_fem,"lognorm", p_fd_fem, 1)
paf_tb_gamma_fem <- calculate_paf_fd(b_tb, cd_fem, "gamma", p_fd_fem, 1)
paf_tb_weibull_fem <- calculate_paf_fd(b_tb, cd_fem, "weibull", p_fd_fem, 1)

list(
  lognorm = paf_tb_lognorm_fem,
  gamma = paf_tb_gamma_fem,
  weibull = paf_tb_weibull_fem
)

# PAF PANCREATITIS MALES
  
paf_pc_lognorm_male <- calculate_paf_fd(b1_panc , cd_male, "lognorm", p_fd_male, 2.2)
paf_pc_gamma_male <- calculate_paf_fd(b1_panc , cd_male, "gamma", p_fd_male, 2.2)
paf_pc_weibull_male <- calculate_paf_fd(b1_panc , cd_male, "weibull", p_fd_male, 2.2)

list(
  lognorm = paf_pc_lognorm_male,
  gamma = paf_pc_gamma_male,
  weibull = paf_pc_weibull_male
)

# PAF COLON AND RECTUM CANCER MALE
paf_crcan_lognorm_male <- calculate_paf_fd(b1_crcan_male , cd_male, "lognorm", p_fd_male, rr_crcan_fd_male)
paf_crcan_gamma_male <- calculate_paf_fd(b1_crcan_male , cd_male, "gamma", p_fd_male, rr_crcan_fd_male)
paf_crcan_weibull_male <- calculate_paf_fd(b1_crcan_male , cd_male, "weibull", p_fd_male, rr_crcan_fd_male)

list(
  lognorm = paf_crcan_lognorm_male,
  gamma = paf_crcan_gamma_male,
  weibull = paf_crcan_weibull_male
)

# PAF COLON AND RECTUM CANCER FEMALE

paf_crcan_lognorm_fem <- calculate_paf_fd(b1_crcan_fem , cd_fem, "lognorm", p_fd_fem, rr_crcan_fd_fem)
paf_crcan_gamma_fem <- calculate_paf_fd(b1_crcan_fem , cd_fem, "gamma", p_fd_fem, rr_crcan_fd_fem)
paf_crcan_weibull_fem <- calculate_paf_fd(b1_crcan_fem , cd_fem, "weibull", p_fd_fem, rr_crcan_fd_fem)

list(
  lognorm = paf_crcan_lognorm_fem,
  gamma = paf_crcan_gamma_fem,
  weibull = paf_crcan_weibull_fem
)

calculate_paf_fd(b1_crcan)
# PAF LIVER CANCER MALE
paf_lican_lognorm_male <- calculate_paf_fd(b1_lican_male , cd_male, "lognorm", p_fd_male, rr_lican_fd_male)
paf_lican_gamma_male <- calculate_paf_fd(b1_lican_male , cd_male, "gamma", p_fd_male, rr_lican_fd_male)
paf_lican_weibull_male <- calculate_paf_fd(b1_lican_male , cd_male, "weibull", p_fd_male, rr_lican_fd_male)

list(
  lognorm = paf_lican_lognorm_male,
  gamma = paf_lican_gamma_male,
  weibull = paf_lican_weibull_male
)

# PAF LIVER CANCER FEMALE
paf_lican_lognorm_fem <- calculate_paf_fd(b1_lican_fem , cd_fem, "lognorm", p_fd_fem, rr_lican_fd_fem)
paf_lican_gamma_fem <- calculate_paf_fd(b1_lican_fem ,  cd_fem, "gamma", p_fd_male, rr_lican_fd_fem)
paf_lican_weibull_fem <- calculate_paf_fd(b1_lican_fem , cd_fem, "weibull", p_fd_male, rr_lican_fd_fem)

list(
  lognorm = paf_lican_lognorm_fem,
  gamma = paf_lican_gamma_fem,
  weibull = paf_lican_weibull_fem
)

# PAF BREAST CANCER
paf_bcan_lognorm_fem <- calculate_paf_fd(b1_bcan , cd_fem, "lognorm", p_fd_fem, rr_bcan_fd)
paf_bcan_gamma_fem <- calculate_paf_fd(b1_bcan , cd_fem, "gamma", p_fd_male, rr_bcan_fd)
paf_bcan_weibull_fem <- calculate_paf_fd(b1_bcan , cd_fem, "weibull", p_fd_male, rr_bcan_fd)

list(
  lognorm = paf_bcan_lognorm_fem,
  gamma = paf_bcan_gamma_fem,
  weibull = paf_bcan_weibull_fem
)


# PROBAR LA VERSIÓN CATEGÓRICA
# PROBAR LA DISTRIBUCIÓN LA EMPÍRICA
# PROBAR REGRESIÓN GAMMA
# F INVERSA




fit <- fitdist(cd_vector, "gamma")

prevalence_function <- function(x) {
  dgamma(x, shape = fit$estimate["shape"], rate = fit$estimate["rate"])
}
rr_function <- function(x) {
  exp(beta_1 * x)
}
integrand <- function(x) {
  prevalence_function(x) * (rr_function(x) - 1)
}
integrate(integrand, 0, 150)$value



shape <- fit$estimate[["shape"]]
rate <- fit$estimate[["rate"]]
scale <- 1 / rate
gamma_mean <- shape * scale
gamma_variance <- shape * (scale^2)
gamma_sd <- sqrt(gamma_variance)
mean(cd_vector)

sd(cd_vector)



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


