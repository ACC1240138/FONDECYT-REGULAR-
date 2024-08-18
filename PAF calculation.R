# THEORETICAL DISTRIBUTION ESTIMATION
# The proportion of most diseases caused by alcohol in
# the component cause model in a population is determined by:
#  • The distribution of the volume of exposure
#  • The relative risk associated with each level of exposure
rm(list = ls())
gc()
# List of required packages
required_packages <- c("dplyr", "readr", "haven", 
                       "ggplot2", "gridExtra", "fitdistrplus", "MASS")

# Install and load required packages
sapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

data <- readRDS("enpg_full.RDS") %>% 
  mutate(aux = ifelse(oh1 == "No" & !is.na(oh2) ,1,0)) %>% 
  filter(aux == 0) %>% 
  dplyr::select(-aux)



# BETAS OF RELATIVE RISK

# TUBERCULOSIS
b_tb <- 0.0179695
var_tb <- 0.0072152^2

# HIV/AIDS (Male)
b_hiv_male <- log(1.54)
# for consumption > 61 grams/day
var_hiv_male <- 0.0782107722^2

# HIV/AIDS (Female)
b_hiv_fem <- log(1.54)
# for consumption > 49 grams/day
var_hiv_fem <- 0.0782107722^2

# LOWER RESPORATORY INFECTIONS
b_lri <- 0.4764038
var_lri <- 0.19220552^2

# LIP AND ORAL CAVITY CANCER
b1_locan <- 0.02474
b2_locan <- -0.00004
rr_locan_fd <- 1.2
# Variance-covariance matrix:
# Variance (b1): 0.000002953
# Variance (b2): 0.000000000102
# Covariance: -0.0000000127

# Other Pharyngeal Cancers
b1_opcan <- 0.02474
b2_opcan <- -0.00004
rr_opcan_fd <- 1.2
var_b1_opcan <- 0.000002953
var_b2_opcan <- 0.000000000102
cov_opcan <- -0.0000000127

# Oesophagus Cancer
b1_oescan <- 0.0132063596418668
b2_oescan <- -4.14801974664481*10^-08
rr_oescan_fd <- 1.16
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
b2_dm_male <- -0.0728256

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
b1_panc_fem <- -0.0272886
b2_panc_fem <- 0.0611466

# ISCHAEMIC HEART DISEASE MALE 15 TO 34
b1_ihd_15_male <- 1.111874
b2_ihd_15_male <- -0.4870068
b3_ihd_15_male <- 1.550984
b4_ihd_15_male <- 0.012

# ISCHAEMIC HEART DISEASE MALE 35 TO 64
b1_ihd_35_male <- 0.757104
b2_ihd_35_male <- -0.4870068
b3_ihd_35_male <- 1.550984
b4_ihd_35_male <- 0.012

# ISCHAEMIC HEART DISEASE FEMALE 15 TO 34
b1_ihd_15_fem <- 1.111874
b2_ihd_15_fem <- 1.832441
b3_ihd_15_fem <- 1.538557
b4_ihd_15_fem <- 0.01

# ISCHAEMIC HEART DISEASE FEMALE 35 TO 64
b1_ihd_35_fem = 1.035623
b2_ihd_35_fem = 1.832441
b3_ihd_35_fem = 1.538557
b4_ihd_35_fem = 0.01

###############################################
# CONTINOUS VERSION USING GAMMA DISTRIBUTION #
###############################################
data_input <- data %>% 
  mutate(edad_tramo = case_when(between(edad, 15, 29)~1,
                                between(edad, 30,44)~2,
                                between(edad,45,59)~3,
                                between(edad,60,65)~4),
         cvolaj = case_when(oh1 == "No" ~ "ltabs",
                            oh2 == ">30" | oh2 == ">1 año" ~ "fd",
                            sexo == "Mujer" & volajohdia > 0 & volajohdia <= 19.99 ~ "cat1",
                            sexo == "Mujer" & volajohdia >= 20 & volajohdia <= 39.99 ~ "cat2",
                            sexo == "Mujer" & volajohdia >= 40 & volajohdia <= 60  ~ "cat3",
                            sexo == "Mujer" & volajohdia > 60 ~ "cat4",
                            sexo == "Hombre" & volajohdia > 0 & volajohdia <= 39.99 ~ "cat1",
                            sexo == "Hombre" & volajohdia >= 40 & volajohdia <= 59.99 ~ "cat2",
                            sexo == "Hombre" & volajohdia >= 60 & volajohdia<= 100 ~ "cat3",
                            sexo == "Hombre" & volajohdia > 100 ~ "cat4",
                            TRUE ~ NA)) %>% 
  dplyr::select(sexo, exp, edad_tramo, volajohdia, cvolaj)

# WE NEED THE PROPORTION OF MALES AND FEMALES IN POPULATION FOR POOLED PAF
input <- data_input %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo, edad_tramo, cvolaj) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2))

props <- data_input %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo, edad_tramo) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2))



p_abs_male1 <- sum(input[6,5])
p_abs_male2 <- sum(input[12,5])
p_abs_male3 <- sum(input[18,5])
p_abs_male4 <- sum(input[24,5])
p_abs_male <- c(p_abs_male1,p_abs_male2, p_abs_male3, p_abs_male4)
p_form_male1 <- sum(input[5,5])
p_form_male2 <- sum(input[11,5])
p_form_male3 <- sum(input[17,5])
p_form_male4 <- sum(input[23,5])
p_form_male <- c(p_form_male1, p_form_male2, p_form_male3, p_form_male4)
p_abs_fem1 <- sum(input[30,5])
p_abs_fem2 <- sum(input[36,5])
p_abs_fem3 <- sum(input[42,5])
p_abs_fem4 <- sum(input[48,5])
p_abs_fem <- c(p_abs_fem1,p_abs_fem2,p_abs_fem3,p_abs_fem4)
p_form_fem1 <- sum(input[29,5])
p_form_fem2 <- sum(input[35,5])
p_form_fem3 <- sum(input[41,5])
p_form_fem4 <- sum(input[47,5])
p_form_fem <- c(p_form_fem1,p_form_fem2,p_form_fem3,p_form_fem4)

x_vals <- seq(0.1, 150, length.out = 1500)

cd_fem1 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1) %>%
  pull(volajohdia)
cd_fem2 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2) %>%
  pull(volajohdia)
cd_fem3 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3) %>%
  pull(volajohdia)
cd_fem4 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4) %>%
  pull(volajohdia)

gamma_fem1 <- fitdist(cd_fem1, "gamma")
gamma_fem2 <- fitdist(cd_fem2, "gamma")
gamma_fem3 <- fitdist(cd_fem3, "gamma")
gamma_fem4 <- fitdist(cd_fem4, "gamma")
y_gamma_fem1 <- dgamma(x_vals, shape = gamma_fem1$estimate["shape"],
                       rate = gamma_fem1$estimate["rate"])
y_gamma_fem2 <- dgamma(x_vals, shape = gamma_fem2$estimate["shape"],
                       rate = gamma_fem2$estimate["rate"])
y_gamma_fem3 <- dgamma(x_vals, shape = gamma_fem3$estimate["shape"],
                       rate = gamma_fem3$estimate["rate"])
y_gamma_fem4 <- dgamma(x_vals, shape = gamma_fem4$estimate["shape"],
                       rate = gamma_fem4$estimate["rate"])

#
cd_male1 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1) %>%
  pull(volajohdia)
cd_male2 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2) %>%
  pull(volajohdia)
cd_male3 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3) %>%
  pull(volajohdia)
cd_male4 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4) %>%
  pull(volajohdia)
gamma_male1 <- fitdist(cd_male1, "gamma")
gamma_male2 <- fitdist(cd_male2, "gamma")
gamma_male3 <- fitdist(cd_male3, "gamma")
gamma_male4 <- fitdist(cd_male4, "gamma")
y_gamma_male1 <- dgamma(x_vals, shape = gamma_male1$estimate["shape"],
                        rate = gamma_male1$estimate["rate"])
y_gamma_male2 <- dgamma(x_vals, shape = gamma_male2$estimate["shape"],
                        rate = gamma_male2$estimate["rate"])
y_gamma_male3 <- dgamma(x_vals, shape = gamma_male3$estimate["shape"],
                        rate = gamma_male3$estimate["rate"])
y_gamma_male4 <- dgamma(x_vals, shape = gamma_male4$estimate["shape"],
                        rate = gamma_male4$estimate["rate"])

# TRAPEZOIDAL INTEGRATION FUNCTION
trap_int <- function(x, y, rr, prop_abs, rr_form, prop_form) {
  
  # Interval width
  dx <- x[2] - x[1]
  ncgamma <- sum(y[-1] + y[-length(x)]) * dx / 2
  # Normalize the gamma function
  normalized_y <- (1 -  (prop_abs+prop_form)) * y/ncgamma
  
  # Calculate the excess relative risk
  excess_rr <- rr - 1
  
  # Calculate the weighted excess relative risk
  weighted_excess_rr <- normalized_y * excess_rr
  
  # Apply the trapezoidal rule to calculate the numerator
  numerator <- (rr_form-1)*prop_form+sum((weighted_excess_rr[-1] + weighted_excess_rr[-length(x)]) / 2) * dx
  
  # Calculate the denominator
  denominator <- numerator + 1
  
  # Calculate PAF
  paf <- round(numerator / denominator, 3)
  
  return(paf)
}

# LINEAR RELATIVE RISK FUNCTION
rr_linear <- function(x, b){
  exp(x*b)
}

# FUNCTION FOR POOLED BY SEX
paf_pool <- function(paf1, paf2, paf3, paf4, sexo) {
  if (sexo == "Mujer") {
    return(paf1 * 0.3 + paf2 * 0.3 + paf3 * 0.29 + paf4 * 0.11)
  } else {
    return(paf1 * 0.36 + paf2 * 0.28 + paf3 * 0.29 + paf4 * 0.09)
  }
}
# FUNCTION FOR OVERALL PAF
paf_overall <- function(paf_fem, paf_male){
  paf_fem*0.42+paf_male*0.58
}
####################################
# ESTIMATING AAF FOR BREAST CANCER #
####################################

# Calculate the relative risk for each x
rr_bcan <- rr_linear(x_vals, b1_bcan)

paf_bc1 <- trap_int(x = x_vals, y = y_gamma_fem1, rr = rr_bcan, prop_abs = p_abs_fem1,rr_form = 1, prop_form = p_form_fem1)
paf_bc2 <- trap_int(x = x_vals, y = y_gamma_fem2, rr = rr_bcan, prop_abs = p_abs_fem2,rr_form = 1,prop_form = p_form_fem2)
paf_bc3 <- trap_int(x = x_vals, y = y_gamma_fem3, rr = rr_bcan, prop_abs = p_abs_fem3,rr_form = 1,prop_form = p_form_fem3)
paf_bc4 <- trap_int(x = x_vals, y = y_gamma_fem4, rr = rr_bcan, prop_abs = p_abs_fem4,rr_form = 1,prop_form = p_form_fem4)
paf_bc <- paf_pool(paf_bc1,paf_bc2,paf_bc3,paf_bc4,"Mujer")

# Simulation parameters
set.seed(145)
# Simulate PCA, proportions of abstainers and former drinkers
# Assuming the same setup as before
n_sim <- 10000

# Pre-allocate a vector for the simulated PAFs
library(MASS)
confint_paf <- function(gamma, beta, var_beta, p_abs, p_form){
set.seed(145)
  n_sim <- 10000
  simulated_pafs <- numeric(n_sim)
  for (i in 1:n_sim) {
  pca_sim <- rgamma(1000, shape = gamma$estimate["shape"], rate = gamma$estimate["rate"])
  # Calculate the shape and rate parameters from the simulated PCA values
  mean_sim <- mean(pca_sim)
  sd_sim <- sd(pca_sim)
  
  # Calculate shape and rate for the new gamma distribution
  shape_sim <- (mean_sim / sd_sim)^2
  rate_sim <- mean_sim / (sd_sim^2)
  
  # Simulate the gamma distribution based on `k_sim` and `theta_sim`
  y_gamma_sim <- dgamma(x_vals, shape = shape_sim, rate = rate_sim)
  
  # Skip iteration if y_gamma_sim contains NaN values
  if (any(is.nan(y_gamma_sim))) next
  
  # Simulate beta coefficient
  beta_sim <- rnorm(1000, beta, sqrt(var_beta))
  
  # Define the relative risk function
  rr_function <- function(x) {
    exp(beta_sim * x)
  }
  # Simulate proportions of lifetime abstainers and former drinkers for the current age category
  prop_abs_sim <- rnorm(1000, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000))
  prop_form_sim <- rnorm(1000, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000))
  
  # Ensure the simulated proportions are positive
  prop_abs_sim <- max(prop_abs_sim, 0.001)
  prop_form_sim <- max(prop_form_sim, 0.001)
  # Calculate the PAF using the trapezoidal method
  simulated_pafs[i] <- trap_int(x = x_vals, y = y_gamma_sim, rr = rr_function(x_vals), 
                                prop_abs = prop_abs_sim, rr_form = 1, prop_form = prop_form_sim)
}

# Remove NaN values from simulated PAFs
simulated_pafs <- simulated_pafs[!is.nan(simulated_pafs)]

# Calculate the 95% confidence interval
paf_lower <- quantile(simulated_pafs, 0.025)
paf_upper <- quantile(simulated_pafs, 0.975)
paf_point_estimate <- mean(simulated_pafs)

return(list(
  Point_Estimate = round(paf_point_estimate,3),
  Lower_CI = paf_lower,
  Upper_CI = paf_upper)
)
}

confint_paf(gamma_fem1, b1_bcan, var_bcan,p_abs_fem1, p_form_fem1)
confint_paf(gamma_fem2, b1_bcan, var_bcan,p_abs_fem2, p_form_fem2)
confint_paf(gamma_fem3, b1_bcan, var_bcan,p_abs_fem3, p_form_fem3)
confint_paf(gamma_fem4, b1_bcan, var_bcan,p_abs_fem4, p_form_fem4)

####################################
# ESTIMATING AAF FOR TUBERCULOSIS  #
####################################
rr_tb <- rr_linear(x_vals, b_tb)
paf_tb_fem1 <- trap_int(x = x_vals, y = y_gamma_fem1, rr = rr_tb, prop_abs = p_abs_fem1,rr_form = 1, prop_form = p_form_fem1)
paf_tb_fem2 <- trap_int(x = x_vals, y = y_gamma_fem2, rr = rr_tb, prop_abs = p_abs_fem2,rr_form = 1, prop_form = p_form_fem2)
paf_tb_fem3 <- trap_int(x = x_vals, y = y_gamma_fem3, rr = rr_tb, prop_abs = p_abs_fem3,rr_form = 1, prop_form = p_form_fem3)
paf_tb_fem4 <- trap_int(x = x_vals, y = y_gamma_fem4, rr = rr_tb, prop_abs = p_abs_fem4,rr_form = 1, prop_form = p_form_fem4)
paf_tb_fem <- paf_pool(paf_tb_fem1,paf_tb_fem2,paf_tb_fem3,paf_tb_fem4,"Mujer")

paf_tb_male1 <- trap_int(x = x_vals, y = y_gamma_male1, rr = rr_tb, prop_abs = p_abs_male1,rr_form = 1, prop_form = p_form_male1)
paf_tb_male2 <- trap_int(x = x_vals, y = y_gamma_male2, rr = rr_tb, prop_abs = p_abs_male2,rr_form = 1, prop_form = p_form_male2)
paf_tb_male3 <- trap_int(x = x_vals, y = y_gamma_male3, rr = rr_tb, prop_abs = p_abs_male3,rr_form = 1, prop_form = p_form_male3)
paf_tb_male4 <- trap_int(x = x_vals, y = y_gamma_male4, rr = rr_tb, prop_abs = p_abs_male4,rr_form = 1, prop_form = p_form_male4)
paf_tb_male <- paf_pool(paf_tb_male1,paf_tb_male2,paf_tb_male3,paf_tb_male4,"Hombre")
paf_tb <- paf_overall(paf_tb_fem, paf_tb_male)

confint_paf(gamma_fem1, b_tb, var_tb,p_abs_fem1, p_form_fem1)
confint_paf(gamma_fem2, b_tb, var_tb,p_abs_fem2, p_form_fem2)
confint_paf(gamma_fem3, b_tb, var_tb,p_abs_fem3, p_form_fem3)
confint_paf(gamma_fem4, b_tb, var_tb,p_abs_fem4, p_form_fem4)

confint_paf(gamma_male1, b_tb, var_tb,p_abs_male1, p_form_male1)
confint_paf(gamma_male2, b_tb, var_tb,p_abs_male2, p_form_male2)
confint_paf(gamma_male3, b_tb, var_tb,p_abs_male3, p_form_male3)
confint_paf(gamma_male4, b_tb, var_tb,p_abs_male4, p_form_male4)

################################################
# ESTIMATING AAF FOR OTHER PHARINGEAL CANCER   #
################################################

confint_paf_vcov <- function(gamma, b1, b2, var_b1, var_b2, cov_b1_b2, p_abs, p_form, rr_fd) {
  set.seed(145)
  n_sim <- 10000
  simulated_pafs <- numeric(n_sim)
  
  # Construct the variance-covariance matrix
  var_cov_matrix <- matrix(c(var_b1, cov_b1_b2, cov_b1_b2, var_b2), nrow = 2)
  
  for (i in 1:n_sim) {
    # Simulate PCA mean and SD using the gamma distribution
    pca_sim <- rgamma(1000, shape = gamma$estimate["shape"], rate = gamma$estimate["rate"])
    
    # Calculate the shape and rate parameters from the simulated PCA values
    mean_sim <- mean(pca_sim)
    sd_sim <- sd(pca_sim)
    
    # Calculate shape and rate for the new gamma distribution
    shape_sim <- (mean_sim / sd_sim)^2
    rate_sim <- mean_sim / (sd_sim^2)
    
    # Simulate the gamma distribution
    y_gamma_sim <- dgamma(x_vals, shape = shape_sim, rate = rate_sim)
    
    # Skip iteration if y_gamma_sim contains NaN values
    if (any(is.nan(y_gamma_sim))) next
    
    # Simulate the beta coefficients jointly from a multivariate normal distribution
    beta_sim <- MASS::mvrnorm(1, mu = c(b1, b2), Sigma = var_cov_matrix)
    
    # Define the relative risk function
    rr_function <- function(x) {
      exp(beta_sim[1] * x + beta_sim[2] * x^2)
    }
    
    # Simulate proportions of lifetime abstainers and former drinkers
    prop_abs_sim <- max(rnorm(1, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000)), 0.001)
    prop_form_sim <- max(rnorm(1, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000)), 0.001)
    
    # Calculate the PAF using the trapezoidal method
    simulated_pafs[i] <- trap_int(x = x_vals, y = y_gamma_sim, rr = rr_function(x_vals), 
                                  prop_abs = prop_abs_sim, rr_form = rr_fd, prop_form = prop_form_sim)
  }
  
  # Remove NaN values from simulated PAFs
  simulated_pafs <- simulated_pafs[!is.nan(simulated_pafs)]
  
  # Calculate the 95% confidence interval
  paf_lower <- quantile(simulated_pafs, 0.025)
  paf_upper <- quantile(simulated_pafs, 0.975)
  paf_point_estimate <- mean(simulated_pafs)
  
  return(list(point_estimate = paf_point_estimate, lower_ci = paf_lower, upper_ci = paf_upper))
}

rr_opcan_cal <- function(b1,b2,x){
  exp(b1*x+b2*x**2)
}
rr_opcan <- rr_opcan_cal(b1_opcan, b2_opcan, x_vals)

paf_opcan_fem1 <- trap_int(x = x_vals, y = y_gamma_fem1, rr = rr_opcan, prop_abs = p_abs_fem1,rr_form = rr_opcan_fd, prop_form = p_form_fem1)
paf_opcan_fem2 <- trap_int(x = x_vals, y = y_gamma_fem2, rr = rr_opcan, prop_abs = p_abs_fem2,rr_form = rr_opcan_fd, prop_form = p_form_fem2)
paf_opcan_fem3 <- trap_int(x = x_vals, y = y_gamma_fem3, rr = rr_opcan, prop_abs = p_abs_fem3,rr_form = rr_opcan_fd, prop_form = p_form_fem3)
paf_opcan_fem4 <- trap_int(x = x_vals, y = y_gamma_fem4, rr = rr_opcan, prop_abs = p_abs_fem4,rr_form = rr_opcan_fd, prop_form = p_form_fem4)
paf_opcan_fem <- paf_pool(paf_opcan_fem1,paf_opcan_fem2,paf_opcan_fem3,paf_opcan_fem4, "Mujer")

paf_opcan_male1 <- trap_int(x = x_vals, y = y_gamma_male1, rr = rr_opcan, prop_abs = p_abs_male1,rr_form = rr_opcan_fd, prop_form = p_form_male1)
paf_opcan_male2 <- trap_int(x = x_vals, y = y_gamma_male2, rr = rr_opcan, prop_abs = p_abs_male2,rr_form = rr_opcan_fd, prop_form = p_form_male2)
paf_opcan_male3 <- trap_int(x = x_vals, y = y_gamma_male3, rr = rr_opcan, prop_abs = p_abs_male3,rr_form = rr_opcan_fd, prop_form = p_form_male3)
paf_opcan_male4 <- trap_int(x = x_vals, y = y_gamma_male4, rr = rr_opcan, prop_abs = p_abs_male4,rr_form = rr_opcan_fd, prop_form = p_form_male4)
paf_opcan_male <- paf_pool(paf_opcan_male1,paf_opcan_male2,paf_opcan_male3,paf_opcan_male4, "Hombre")
paf_opcan <- paf_overall(paf_opcan_fem, paf_opcan_male)

confint_paf_vcov(gamma_fem1, b1_opcan, b2_opcan, var_b1_opcan, var_b2_opcan, cov_opcan, p_abs_fem1, p_form_fem1, rr_opcan_fd)
confint_paf_vcov(gamma_fem2, b1_opcan, b2_opcan, var_b1_opcan, var_b2_opcan, cov_opcan, p_abs_fem2, p_form_fem2, rr_opcan_fd)
confint_paf_vcov(gamma_fem3, b1_opcan, b2_opcan, var_b1_opcan, var_b2_opcan, cov_opcan, p_abs_fem3, p_form_fem3, rr_opcan_fd)
confint_paf_vcov(gamma_fem4, b1_opcan, b2_opcan, var_b1_opcan, var_b2_opcan, cov_opcan, p_abs_fem4, p_form_fem4, rr_opcan_fd)

confint_paf_vcov(gamma_male1, b1_opcan, b2_opcan, var_b1_opcan, var_b2_opcan, cov_opcan, p_abs_male1, p_form_male1, rr_opcan_fd)
confint_paf_vcov(gamma_male2, b1_opcan, b2_opcan, var_b1_opcan, var_b2_opcan, cov_opcan, p_abs_male2, p_form_male2, rr_opcan_fd)
confint_paf_vcov(gamma_male3, b1_opcan, b2_opcan, var_b1_opcan, var_b2_opcan, cov_opcan, p_abs_male3, p_form_male3, rr_opcan_fd)
confint_paf_vcov(gamma_male4, b1_opcan, b2_opcan, var_b1_opcan, var_b2_opcan, cov_opcan, p_abs_male4, p_form_male4, rr_opcan_fd)

##########################################
# ESTIMATING AAF FOR OESOPHAFUS CANCER   #
##########################################
rr_oescan_cal <- function(b1,b2,x){
  exp(b1*x+b2*x**3)
}
rr_oescan <- rr_oescan_cal(b1_oescan, b2_oescan, x_vals)
paf_oescan_fem1 <- trap_int(x = x_vals, y = y_gamma_fem1, rr = rr_opcan, prop_abs = p_abs_fem1,rr_form = rr_oescan_fd, prop_form = p_form_fem1)
paf_oescan_fem2 <- trap_int(x = x_vals, y = y_gamma_fem2, rr = rr_opcan, prop_abs = p_abs_fem2,rr_form = rr_oescan_fd, prop_form = p_form_fem2)
paf_oescan_fem3 <- trap_int(x = x_vals, y = y_gamma_fem3, rr = rr_opcan, prop_abs = p_abs_fem3,rr_form = rr_oescan_fd, prop_form = p_form_fem3)
paf_oescan_fem4 <- trap_int(x = x_vals, y = y_gamma_fem4, rr = rr_opcan, prop_abs = p_abs_fem4,rr_form = rr_oescan_fd, prop_form = p_form_fem4)
paf_oescan_fem <- paf_pool(paf_oescan_fem1,paf_oescan_fem2,paf_oescan_fem3,paf_oescan_fem4, "Mujer")

paf_oescan_male1 <- trap_int(x = x_vals, y = y_gamma_male1, rr = rr_opcan, prop_abs = p_abs_male1,rr_form = rr_oescan_fd, prop_form = p_form_male1)
paf_oescan_male2 <- trap_int(x = x_vals, y = y_gamma_male2, rr = rr_opcan, prop_abs = p_abs_male2,rr_form = rr_oescan_fd, prop_form = p_form_male2)
paf_oescan_male3 <- trap_int(x = x_vals, y = y_gamma_male3, rr = rr_opcan, prop_abs = p_abs_male3,rr_form = rr_oescan_fd, prop_form = p_form_male3)
paf_oescan_male4 <- trap_int(x = x_vals, y = y_gamma_male4, rr = rr_opcan, prop_abs = p_abs_male4,rr_form = rr_oescan_fd, prop_form = p_form_male4)
paf_oescan_male <- paf_pool(paf_oescan_male1,paf_oescan_male2,paf_oescan_male3,paf_oescan_male4, "Hombre")

paf_oescan <- paf_overall(paf_oescan_fem,paf_oescan_male)

################################################
# ESTIMATING AAF FOR COLON AND RECTUM CANCER   #
################################################
rr_crcan_male <- rr_linear(b1_crcan_male, x_vals)
rr_crcan_fem <- rr_linear(b1_crcan_fem, x_vals)
paf_crcan_fem <- trap_int_fd(x_vals, y_gamma_fem, rr_crcan_fem, 1.05, 0.3)
paf_crcan_male <- trap_int_fd(x_vals, y_gamma_male, rr_crcan_male, 2.19, 0.39)
paf_crcan <- paf_crcan_fem * 0.42 + paf_crcan_male * 0.58

######################################
# ESTIMATING AAF FOR LIVER CANCER   #
#####################################
rr_lican_male <- rr_linear(b1_lican_male, x_vals)
rr_lican_fem <- rr_linear(b1_lican_fem, x_vals)
paf_lican_male <- trap_int_fd(x_vals, y_gamma_male, rr_lican_male, 2.23, 0.39)
paf_lican_fem <- trap_int_fd(x_vals, y_gamma_fem, rr_lican_fem, 2.68, 0.3)
paf_lican <- paf_lican_fem * 0.42 + paf_lican_male * 0.58

######################################
# ESTIMATING AAF FOR LARYNX CANCER   #
#####################################
rr_lxcan <- rr_linear(b1_lxcan, x_vals)
paf_lxcan_male <- trap_int(x_vals, y_gamma_male, rr_lxcan)
paf_lxcan_fem <- trap_int(x_vals, y_gamma_fem, rr_lxcan)
paf_lxcan <- paf_lxcan_fem * 0.42 + paf_lxcan_male * 0.58

####### CATEGORICAL VERSION
table(data$cvolaj)


# STRATIFIED CONSUMPTION LEVEL
oh_level_median <- capped_data %>% 
  filter(!is.na(volajohdia),!is.na(cvolaj)) %>% 
  group_by(sexo, edad_tramo,cvolaj) %>%
  summarise(oh_level = round(median(volajohdia),1))

oh_level_median2 <- capped_data %>% 
  filter(!is.na(volajohdia),!is.na(cvolaj)) %>% 
  group_by(sexo,cvolaj) %>%
  summarise(oh_level = round(median(volajohdia),1))

# STRATIFIED PROPORTION OF CONSUMPTION LEVEL
prop_level <- capped_data %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo, edad_tramo,cvolaj) %>%
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>%
  group_by(sexo, edad_tramo) %>%
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count)

prop_level2 <- capped_data %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo,cvolaj) %>%
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>%
  group_by(sexo) %>%
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count)

# JOIN BOTH TABLES
paf_base_median <- oh_level_median %>% 
  left_join(prop_level, by = c("sexo", "edad_tramo","cvolaj"))

paf_base_median2 <- oh_level_median2 %>% 
  left_join(prop_level2, by = c("sexo","cvolaj"))

####################################
# ESTIMATING AAF FOR BREAST CANCER #
####################################

paf_bc_median <- paf_base_median %>% 
  filter(sexo=="Mujer") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b1_bcan*oh_level),1),NA))
paf_bc_median2 <- paf_base_median2 %>% 
  filter(sexo=="Mujer") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b1_bcan*oh_level),1),NA))


cat_paf_calculator <- function(data){
  num1 = data[1,5]*(data[1,6]-1)+data[2,5]*(data[2,6]-1)+data[3,5]*(data[3,6]-1)+data[4,5]*(data[4,6]-1)
  den1 = num1+1
  num2 = data[7,5]*(data[7,6]-1)+data[8,5]*(data[8,6]-1)+data[9,5]*(data[9,6]-1)+data[10,5]*(data[10,6]-1)
  den2 = num2+1
  num3 = data[13,5]*(data[13,6]-1)+data[14,5]*(data[14,6]-1)+data[15,5]*(data[15,6]-1)+data[16,5]*(data[16,6]-1)
  den3 = num3+1
  num4 = data[19,5]*(data[19,6]-1)+data[20,5]*(data[20,6]-1)+data[21,5]*(data[21,6]-1)+data[22,5]*(data[22,6]-1)
  den4 = num4+1
  
  # Calculate PAFs
  paf_15_29 <- num1 / den1
  paf_30_44 <- num2 / den2
  paf_45_59 <- num3 / den3
  paf_60_65 <- num4 / den4
  
  # Store the results in a list
  result <- list(
    "PAF for 15-29" = round(paf_15_29,2),
    "PAF for 30-44" = round(paf_30_44,2),
    "PAF for 45-59" = round(paf_45_59,2),
    "PAF for 60-65" = round(paf_60_65,2)
  )
  
  return(result)
}
cat_paf_calculator(paf_bc_median)
paf_global <- function(data){
  num1 = data[1,4]*(data[1,5]-1)+data[2,4]*(data[2,5]-1)+data[3,4]*(data[3,5]-1)+data[4,4]*(data[4,5]-1)
  den1 = num1+1
  print(round(num1/den1,2))
}
paf_global(paf_bc_median2)

####################################
# ESTIMATING AAF FOR TUBERCULOSIS  #
####################################
paf_tb_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs", cvolaj != "fd") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b_tb*oh_level),1),NA))
paf_tb_median2 <- paf_base_median2 %>% 
  filter(cvolaj != "ltabs", cvolaj != "fd") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b_tb*oh_level),1),NA))

cat_paf <- function(data){
  num1 = data[1,5]*(data[1,6]-1)+data[2,5]*(data[2,6]-1)+data[3,5]*(data[3,6]-1)+data[4,5]*(data[4,6]-1)
  den1 = num1+1
  num2 = data[5,5]*(data[5,6]-1)+data[6,5]*(data[6,6]-1)+data[7,5]*(data[7,6]-1)+data[8,5]*(data[8,6]-1)
  den2 = num2+1
  num3 = data[9,5]*(data[9,6]-1)+data[10,5]*(data[10,6]-1)+data[11,5]*(data[11,6]-1)+data[12,5]*(data[12,6]-1)
  den3 = num3+1
  num4 = data[13,5]*(data[13,6]-1)+data[14,5]*(data[14,6]-1)+data[15,5]*(data[15,6]-1)+data[16,5]*(data[16,6]-1)
  den4 = num4+1
  
  num5 = data[17,5]*(data[17,6]-1)+data[18,5]*(data[18,6]-1)+data[19,5]*(data[19,6]-1)+data[20,5]*(data[20,6]-1)
  den5 = num5+1
  num6 = data[21,5]*(data[21,6]-1)+data[22,5]*(data[22,6]-1)+data[23,5]*(data[23,6]-1)+data[24,5]*(data[24,6]-1)
  den6 = num6+1
  num7 = data[25,5]*(data[25,6]-1)+data[26,5]*(data[26,6]-1)+data[27,5]*(data[27,6]-1)+data[28,5]*(data[28,6]-1)
  den7 = num7+1
  num8 = data[29,5]*(data[29,6]-1)+data[30,5]*(data[30,6]-1)+data[31,5]*(data[31,6]-1)+data[32,5]*(data[32,6]-1)
  den8 = num8+1 
   # Calculate PAFs
  paf_15_29_m <- num1 / den1
  paf_30_44_m <- num2 / den2
  paf_45_59_m <- num3 / den3
  paf_60_65_m <- num4 / den4
  paf_15_29_f <- num5 / den5
  paf_30_44_f <- num6 / den6
  paf_45_59_f <- num7 / den7
  paf_60_65_f <- num8 / den8
  # Store the results in a list
  result <- list(
    "PAF for 15-29 (Male)" = round(paf_15_29_m,3),
    "PAF for 30-44 (Male)" = round(paf_30_44_m,3),
    "PAF for 45-59 (Male)" = round(paf_45_59_m,3),
    "PAF for 60-65 (Male)" = round(paf_60_65_m,3),
    "PAF for 15-29 (Fem)" = round(paf_15_29_f,3),
    "PAF for 30-44 (Fem)" = round(paf_30_44_f,3),
    "PAF for 45-59 (Fem)" = round(paf_45_59_f,3),
    "PAF for 60-65 (Fem)" = round(paf_60_65_f,3)
  )
  
  return(result)
}
cat_paf(paf_tb_median)
paf_global2 <- function(data){
  num1 = data[1,4]*(data[1,5]-1)+data[2,4]*(data[2,5]-1)+data[3,4]*(data[3,5]-1)+data[4,4]*(data[4,5]-1)
  den1 = num1+1
  num2 = data[5,4]*(data[5,5]-1)+data[6,4]*(data[6,5]-1)+data[7,4]*(data[7,5]-1)+data[8,4]*(data[8,5]-1)
  den2 = num2+1
  paf1 = num1/den1
  paf2 = num2/den2
  global = paf1*0.58+paf2*0.42
  result <- list(
    "PAF Male" = round(paf1,2),
    "PAF Female" = round(paf2,2),
    "PAF Global" = round(global,2)
  )
  return(result)
}
paf_global2(paf_tb_median2)
###########################
# ESTIMATING AAF FOR HIV  #
###########################

paf_hiv_male <- capped_data %>% 
  filter(sexo == "Hombre", !is.na(volajohdia)) %>% 
  mutate(aux = ifelse(volajohdia >= 61, 1, 0)) %>% 
  group_by(edad_tramo, aux) %>%
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>%
  group_by(edad_tramo) %>%
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count) %>% 
  mutate(rr = ifelse(aux == 1, exp(b_hiv_male), 1)) %>% 
  filter(aux == 1)

paf_hiv_fem <- capped_data %>% 
  filter(sexo == "Mujer", !is.na(volajohdia)) %>% 
  mutate(aux = ifelse(volajohdia >= 49, 1, 0)) %>% 
  group_by(edad_tramo, aux) %>%
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>%
  group_by(edad_tramo) %>%
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count) %>% 
  mutate(rr = ifelse(aux == 1, exp(b_hiv_male), 1)) %>% 
  filter(aux == 1)

paf_hiv_calc <- function(data){
  num1 = data[1,3]*(data[1,4]-1)
  den1 = num1+1
  num2 = data[2,3]*(data[2,4]-1)
  den2 = num2+1
  num3 = data[3,3]*(data[3,4]-1)
  den3 = num3+1
  num4 = data[4,3]*(data[4,4]-1)
  den4 = num4+1
  
  # Calculate PAFs
  paf_15_29 <- num1 / den1
  paf_30_44 <- num2 / den2
  paf_45_59 <- num3 / den3
  paf_60_65 <- num4 / den4
  
  # Store the results in a list
  result <- list(
    "PAF for 15-29 (Male)" = round(paf_15_29,3),
    "PAF for 30-44 (Male)" = round(paf_30_44,3),
    "PAF for 45-59 (Male)" = round(paf_45_59,3),
    "PAF for 60-65 (Male)" = round(paf_60_65,3)
  )
  
  return(result)
}

paf_hiv_calc(paf_hiv_male)
paf_hiv_calc(paf_hiv_fem)

####################################################
# ESTIMATING AAF FOR LOWER RESPIRATORY INFECTIONS  #
####################################################
paf_lri_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs", cvolaj != "fd") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b_lri*((oh_level+0.0399999618530273)/100)),1),NA))

cat_paf(paf_lri_median)


####################################################
# ESTIMATING AAF FOR LIP AND ORAL CAVITY CANCER   #
####################################################
paf_locan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_locan*oh_level+b2_locan*oh_level**2),1)))

cat_paf_fd <- function(data){
  num1 = data[1,5]*(data[1,6]-1)+data[2,5]*(data[2,6]-1)+data[3,5]*(data[3,6]-1)+data[4,5]*(data[4,6]-1)+data[5,5]*(data[5,6]-1)
  den1 = num1+1
  num2 = data[6,5]*(data[6,6]-1)+data[7,5]*(data[7,6]-1)+data[8,5]*(data[8,6]-1)+data[9,5]*(data[9,6]-1)+data[10,5]*(data[10,6]-1)
  den2 = num2+1
  num3 = data[11,5]*(data[11,6]-1)+data[12,5]*(data[12,6]-1)+data[13,5]*(data[13,6]-1)+data[14,5]*(data[14,6]-1)+data[15,5]*(data[15,6]-1)
  den3 = num3+1
  num4 = data[16,5]*(data[16,6]-1)+data[17,5]*(data[17,6]-1)+data[18,5]*(data[18,6]-1)+data[19,5]*(data[19,6]-1)+data[20,5]*(data[20,6]-1)
  den4 = num4+1
  
  num5 = data[21,5]*(data[21,6]-1)+data[22,5]*(data[22,6]-1)+data[23,5]*(data[23,6]-1)+data[24,5]*(data[24,6]-1)+data[25,5]*(data[25,6]-1)
  den5 = num5+1
  num6 = data[26,5]*(data[26,6]-1)+data[27,5]*(data[27,6]-1)+data[28,5]*(data[28,6]-1)+data[29,5]*(data[29,6]-1)+data[30,5]*(data[30,6]-1)
  den6 = num6+1
  num7 = data[31,5]*(data[31,6]-1)+data[32,5]*(data[32,6]-1)+data[33,5]*(data[33,6]-1)+data[34,5]*(data[34,6]-1)+data[35,5]*(data[35,6]-1)
  den7 = num7+1
  num8 = data[36,5]*(data[36,6]-1)+data[37,5]*(data[37,6]-1)+data[38,5]*(data[38,6]-1)+data[39,5]*(data[39,6]-1)+data[40,5]*(data[40,6]-1)
  den8 = num8+1 
  # Calculate PAFs
  paf_15_29_m <- num1 / den1
  paf_30_44_m <- num2 / den2
  paf_45_59_m <- num3 / den3
  paf_60_65_m <- num4 / den4
  paf_15_29_f <- num5 / den5
  paf_30_44_f <- num6 / den6
  paf_45_59_f <- num7 / den7
  paf_60_65_f <- num8 / den8
  # Store the results in a list
  result <- list(
    "PAF for 15-29 (Male)" = round(paf_15_29_m,3),
    "PAF for 30-44 (Male)" = round(paf_30_44_m,3),
    "PAF for 45-59 (Male)" = round(paf_45_59_m,3),
    "PAF for 60-65 (Male)" = round(paf_60_65_m,3),
    "PAF for 15-29 (Fem)" = round(paf_15_29_f,3),
    "PAF for 30-44 (Fem)" = round(paf_30_44_f,3),
    "PAF for 45-59 (Fem)" = round(paf_45_59_f,3),
    "PAF for 60-65 (Fem)" = round(paf_60_65_f,3)
  )
  
  return(result)
  
}

cat_paf_fd(paf_locan_median)
cat_paf(paf_locan_median %>% filter(cvolaj != "fd"))


################################################
# ESTIMATING AAF FOR OTHER PHARINGEAL CANCER   #
################################################
paf_opcan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_locan*oh_level+b2_locan*oh_level**2),1)))

paf_opcan_median2 <- paf_base_median2 %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_locan*oh_level+b2_locan*oh_level**2),1)))

cat_paf_fd(paf_opcan_median)
cat_paf(paf_opcan_median %>% filter(cvolaj != "fd"))

paf_global_fd <- function(data){
  num1 = data[5,4]*(data[5,5]-1)+data[1,4]*(data[1,5]-1)+data[2,4]*(data[2,5]-1)+data[3,4]*(data[3,5]-1)+data[4,4]*(data[4,5]-1)
  den1 = num1+1
  num2 = data[10,4]*(data[10,5]-1)+data[6,4]*(data[6,5]-1)+data[7,4]*(data[7,5]-1)+data[8,4]*(data[8,5]-1)+data[9,4]*(data[9,5]-1)
  den2 = num2+1
  paf1 = num1/den1
  paf2 = num2/den2
  global = paf1*0.58+paf2*0.42
  result <- list(
    "PAF Male" = round(paf1,2),
    "PAF Female" = round(paf2,2),
    "PAF Global" = round(global,2)
  )
  return(result)
}
paf_global_fd(paf_opcan_median2)
##########################################
# ESTIMATING AAF FOR OESOPHAFUS CANCER   #
##########################################
paf_oescan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.16,round(exp(b1_oescan*oh_level+b2_oescan*oh_level**3),1)))

paf_oescan_median2 <- paf_base_median2 %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.16,round(exp(b1_oescan*oh_level+b2_oescan*oh_level**3),1)))

cat_paf_fd(paf_oescan_median)

cat_paf(paf_oescan_median %>% filter(cvolaj != "fd"))
paf_global_fd(paf_oescan_median2)
################################################
# ESTIMATING AAF FOR COLON AND RECTUM CANCER   #
################################################
paf_crcan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.19,
                        cvolaj == "fd" & sexo == "Mujer"~1.05,
                        sexo == "Hombre"~exp(b1_crcan_male*oh_level),
                        sexo == "Mujer"~exp(b1_crcan_fem*oh_level)))
paf_crcan_median2 <- paf_base_median2 %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.19,
                        cvolaj == "fd" & sexo == "Mujer"~1.05,
                        sexo == "Hombre"~exp(b1_crcan_male*oh_level),
                        sexo == "Mujer"~exp(b1_crcan_fem*oh_level)))
cat_paf_fd(paf_crcan_median)
cat_paf(paf_crcan_median %>% filter(cvolaj != "fd"))
paf_global_fd(paf_crcan_median2)
######################################
# ESTIMATING AAF FOR LIVER CANCER   #
#####################################
paf_lican_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.23,
                        cvolaj == "fd" & sexo == "Mujer"~2.68,
                        sexo == "Hombre"~exp(b1_lican_male*oh_level),
                        sexo == "Mujer"~exp(b1_lican_fem*oh_level)))

cat_paf_fd(paf_lican_median)
cat_paf(paf_lican_median %>% filter(cvolaj != "fd"))

paf_lican_median2 <- paf_base_median2 %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.23,
                        cvolaj == "fd" & sexo == "Mujer"~2.68,
                        sexo == "Hombre"~exp(b1_lican_male*oh_level),
                        sexo == "Mujer"~exp(b1_lican_fem*oh_level)))
paf_global2(paf_lican_median2)
######################################
# ESTIMATING AAF FOR LARYNX CANCER   #
#####################################
paf_lxcan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.18,exp(b1_lxcan*oh_level)))
                       
paf_lxcan_median2 <- paf_base_median2 %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.18,exp(b1_lxcan*oh_level)))

cat_paf_fd(paf_lxcan_median)
cat_paf(paf_lxcan_median %>% filter(cvolaj != "fd"))

paf_global2(paf_lxcan_median2)
########################################
# ESTIMATING AAF FOR DIABETES MELLITUS #
########################################
paf_dm_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~1.18,
                        cvolaj == "fd" & sexo == "Mujer"~1.14,
                        sexo == "Hombre"~exp((b1_dm_male*(oh_level/100)**2)+(b2_dm_male*(oh_level/100)**3)),
                        sexo == "Mujer"~exp(b1_dm_fem*sqrt(oh_level/100)+b2_dm_fem*(oh_level/100))))


cat_paf_fd(paf_dm_median)
cat_paf(paf_dm_median %>% filter(cvolaj != "fd"))

################################
# ESTIMATING AAF FOR EPILEPSY #
###############################
paf_epi_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1,exp(b1_epi*((oh_level+0.5)/100))))


cat_paf_fd(paf_epi_median)
cat_paf(paf_epi_median %>% filter(cvolaj != "fd"))

#################################################
# ESTIMATING AAF FOR HYPERTENSIVE HEART DISEASE #
#################################################

paf_hhd_male <- capped_data %>% 
  filter(sexo == "Hombre", !is.na(volajohdia), volajohdia > 0) %>% 
  mutate(aux = case_when(volajohdia <21~1,
                         volajohdia >=21 & volajohdia < 75~2,
                         volajohdia >= 75~3),
         sexo = "Hombre") %>% 
  group_by(sexo,edad_tramo, aux) %>%
  summarise(weighted_count = sum(exp, na.rm = TRUE),
            oh_level = median(volajohdia)) %>%
  group_by(edad_tramo) %>%
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count) %>% 
  mutate(rr = case_when(aux == 1 ~ exp(b1_hhd_male*oh_level-b2_hhd_male*(oh_level**3/75**2)),
                        aux == 2 ~ exp(b1_hhd_male*oh_level+b2_hhd_male*((oh_level**3-(oh_level-21)**3*75)/(75-21))/75**2),
                        aux == 3 ~ exp(b1_hhd_male*oh_level+b2_hhd_male*(oh_level**3-((oh_level-21)**3*75-(oh_level-75)**3*21)/(75-21))/75**2)))

paf_hhd_fem <- capped_data %>% 
  filter(sexo == "Mujer", !is.na(volajohdia), volajohdia >0) %>% 
  mutate(aux = case_when(volajohdia <18.99~1,
                         volajohdia >=18.99 & volajohdia < 75~2,
                         volajohdia >= 75~3),
         sexo = "Mujer") %>%  
  group_by(sexo,edad_tramo, aux) %>%
  summarise(weighted_count = sum(exp, na.rm = TRUE),
            oh_level = median(volajohdia)) %>% 
  group_by(edad_tramo) %>%
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count) %>% 
  mutate(rr = case_when(aux == 1 ~ 1,
                        aux == 2 ~ exp(b1_hhd_fem*oh_level+b2_hhd_fem*(oh_level**3-((oh_level-10)**3*20-(oh_level-20)**3*10)/(20-10))/20**2),
                        aux == 3 ~ exp(b1_hhd_fem*116+b2_hhd_fem*((75**3-(((75-10)**3*20-(75-20)**3*10)/(20-10)))/20**2))))
paf_hhd <- bind_rows(paf_hhd_male, paf_hhd_fem)
paf_hhd <- paf_hhd %>% 
  mutate(rxp = prop * (rr - 1)) %>%  # Calculate the product (proportion * (rr - 1))
  group_by(sexo,edad_tramo) %>% 
  summarise(numerator = sum(rxp, na.rm = TRUE),  # Sum the products within each age group
            denominator = numerator + 1) %>%     # Add 1 to the numerator to get the denominator
  mutate(paf = round(numerator / denominator,2)) %>%      # Calculate PAF
  ungroup()

#########################################
# ISCHAEMIC HEART DISEASE MALE 15 TO 34 #
#########################################

# Function to calculate y1 and y2
calculate_y <- function(x) {
  return((x + 0.0099999997764826) / 100)
}

ihd_15_male <- capped_data %>% 
  filter(edad < 34 & sexo == "Hombre" & volajohdia > 0) %>% 
  mutate(aux = case_when(volajohdia < 60 ~ 1, 
                         volajohdia >= 60 & volajohdia < 100~2,
                         volajohdia > 100~3)) %>% 
  group_by(aux) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE),
            oh_level = median(volajohdia)) %>% 
  ungroup() %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count) %>% 
  mutate(y1 = calculate_y(oh_level),
         y2 = calculate_y(60),
         rr = case_when(aux == 1 ~ exp(b1_ihd_15_male * (b2_ihd_15_male * sqrt(y1) + b3_ihd_15_male * y1^3)),
                   aux == 2 ~ 0.04571551 + exp(b1_ihd_15_male * (b2_ihd_15_male * sqrt(y2) + b3_ihd_15_male * y2^3)),
                   aux == 3 ~ exp(b4_ihd_15_male * (oh_level - 100)) - 1 + 0.04571551 + exp(b1_ihd_15_male * (b2_ihd_15_male * sqrt(y2) + b3_ihd_15_male * y2^3))),
         sexo = "Hombre",
         edad = "15 a 34") %>% 
  dplyr::select(-y1, -y2,-aux) %>% 
  dplyr::select(sexo,edad,oh_level,prop,rr)
  
          
#########################################
# ISCHAEMIC HEART DISEASE MALE 35 TO 64 #
#########################################

ihd_35_male <- capped_data %>% 
  filter(edad > 35 & sexo == "Hombre" & volajohdia > 0) %>% 
  mutate(aux = case_when(volajohdia < 60 ~ 1, 
                         volajohdia >= 60 & volajohdia < 100~2,
                         volajohdia > 100~3)) %>% 
  group_by(aux) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE),
            oh_level = median(volajohdia)) %>% 
  ungroup() %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count) %>% 
  mutate(y1 = calculate_y(oh_level),
         y2 = calculate_y(60),
         rr = case_when(aux == 1 ~ exp(b1_ihd_15_male * (b2_ihd_15_male * sqrt(y1) + b3_ihd_15_male * y1^3)),
                        aux == 2 ~ 0.04571551 + exp(b1_ihd_15_male * (b2_ihd_15_male * sqrt(y2) + b3_ihd_15_male * y2^3)),
                        aux == 3 ~ exp(b4_ihd_15_male * (oh_level - 100)) - 1 + 0.04571551 + exp(b1_ihd_15_male * (b2_ihd_15_male * sqrt(y2) + b3_ihd_15_male * y2^3))),
         sexo = "Hombre",
         edad = "35 a 65") %>% 
  dplyr::select(-y1, -y2,-aux) %>% 
  dplyr::select(sexo,edad,oh_level,prop,rr)

###########################################
# ISCHAEMIC HEART DISEASE FEMALE 15 TO 34 #
###########################################
ihd_15_fem <- capped_data %>% 
  filter(edad < 34 & sexo == "Mujer" & volajohdia > 0) %>% 
  mutate(aux = case_when(volajohdia <= 30.4 ~ 1, 
                         volajohdia > 30.4 ~2)) %>% 
  group_by(aux) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE),
            oh_level = median(volajohdia)) %>% 
  ungroup() %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count) %>% 
  mutate(y1 = calculate_y(oh_level),
         y2 = calculate_y(30.4),
         rr = case_when(aux == 1 ~ exp(b1_ihd_15_fem * (b2_ihd_15_fem * y1 + b3_ihd_15_fem * y1 * log(y1))),
                        aux == 2 ~ exp(b4_ihd_15_fem * (oh_level - 30.3814)) - 1 + exp(b1_ihd_15_fem * (b2_ihd_15_fem * y2 + b3_ihd_15_fem * y2 * log(y2)))),
         sexo = "Mujer",
         edad = "15 a 34") %>% 
  dplyr::select(-y1, -y2,-aux) %>% 
  dplyr::select(sexo,edad,oh_level,prop,rr)

###########################################
# ISCHAEMIC HEART DISEASE FEMALE 35 TO 64 #
###########################################



