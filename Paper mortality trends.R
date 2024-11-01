#---------------------------#
# MORTALITY TRENDS IN CHILE #
#---------------------------#

required_packages <- c("dplyr", "readr", "haven", 
                       "ggplot2", "gridExtra", 
                       "fitdistrplus", "MASS")

# Install and load required packages
sapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

data <- readRDS("enpg_full.RDS") %>% 
  mutate(aux = ifelse(oh1 == "No" & !is.na(oh2) ,1,0)) %>% 
  filter(aux == 0, year <= 2018) %>% 
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
var_b1_oescan <- 1.525706*10^-07
var_b2_oescan <- 0.000002953
cov_oescan <- -6.885205*10^-13

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
  dplyr::select(year, sexo, exp, edad_tramo, volajohdia, cvolaj)


input <- data_input %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo,year, edad_tramo, cvolaj) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>% 
  filter(cvolaj == "fd" | cvolaj == "ltabs") %>% 
  dplyr::select(-weighted_count)

input_male <- input %>% 
  filter(sexo == "Hombre")
input_female <- input %>% 
  filter(sexo == "Mujer")

# Function to fit gamma distribution
# Define the x_vals range (adjust as needed)
x_vals <- seq(0.1, 150, length.out = 1500)

# Function to fit gamma distribution and return both parameters and density
fit_gamma_with_density <- function(volajohdia_data, x_vals) {
  if (length(volajohdia_data) > 1) {
    fit <- fitdist(volajohdia_data, "gamma")
    shape <- fit$estimate["shape"]
    rate <- fit$estimate["rate"]
    density <- dgamma(x_vals, shape = shape, rate = rate)
    return(list(shape = shape, rate = rate, density = density))
  } else {
    return(list(shape = NA, rate = NA, density = rep(NA, length(x_vals))))
  }
}

input_gamma <- data_input %>%
  filter(!is.na(volajohdia), volajohdia > 0) %>%
  group_by(sexo, year, edad_tramo) %>%
  summarise(pca = mean(volajohdia),
    gamma_fit = list(fit_gamma_with_density(volajohdia, x_vals))) %>%
  mutate(shape = purrr::map_dbl(gamma_fit, "shape"),
         rate = purrr::map_dbl(gamma_fit, "rate")) %>%
  ungroup() %>%
  dplyr::select(-gamma_fit)

input_female <- input_female %>% 
  left_join(input_gamma, by=c("sexo", "year", "edad_tramo"))

input_male <- input_male %>% 
  left_join(input_gamma, by=c("sexo", "year", "edad_tramo"))

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

confint_paf <- function(gamma, beta, var_beta, p_abs, rr_form, p_form, rr_function){
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
    beta_sim <- rnorm(1, beta, sqrt(var_beta))
    rr_sim <- rr_function(x_vals, beta_sim)
    # Simulate proportions of lifetime abstainers and former drinkers for the current age category
    prop_abs_sim <- max(rnorm(1, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000)),0.001)
    prop_form_sim <- max(rnorm(1, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000)), 0.001)
    
    # Calculate the PAF using the trapezoidal method
    simulated_pafs[i] <- trap_int(x = x_vals, y = y_gamma_sim, rr = rr_sim, 
                                  prop_abs = prop_abs_sim, rr_form = rr_form, prop_form = prop_form_sim)
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

confint_paf_vcov <- function(gamma, betas, cov_matrix, p_abs, p_form, rr_fd, rr_function) {
  set.seed(145)
  n_sim <- 10000
  simulated_pafs <- numeric(n_sim)
  
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
    beta_sim <- MASS::mvrnorm(1, mu = betas, Sigma = cov_matrix)
    
    rr_sim <- rr_function(x_vals, beta_sim)
    
    # Simulate proportions of lifetime abstainers and former drinkers
    prop_abs_sim <- max(rnorm(1, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000)), 0.001)
    prop_form_sim <- max(rnorm(1, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000)), 0.001)
    
    # Calculate the PAF using the trapezoidal method
    simulated_pafs[i] <- trap_int(x = x_vals, y = y_gamma_sim, rr = rr_sim, 
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

cd_fem1_08 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, year == 2008) %>%
  pull(volajohdia)
cd_fem2_08 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, year == 2008) %>%
  pull(volajohdia)
cd_fem3_08 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, year == 2008) %>%
  pull(volajohdia)
cd_fem4_08 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, year == 2008) %>%
  pull(volajohdia)

cd_male1_08 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1, year == 2008) %>%
  pull(volajohdia)
cd_male2_08 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2, year == 2008) %>%
  pull(volajohdia)
cd_male3_08 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3, year == 2008) %>%
  pull(volajohdia)
cd_male4_08 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4, year == 2008) %>%
  pull(volajohdia)

g_fem1_08 <- fitdist(cd_fem1_08, "gamma")
g_fem2_08 <- fitdist(cd_fem2_08, "gamma")
g_fem3_08 <-fitdist(cd_fem3_08, "gamma")
g_fem4_08 <-fitdist(cd_fem4_08, "gamma")

g_male1_08 <- fitdist(cd_male1_08, "gamma")
g_male2_08 <- fitdist(cd_male2_08, "gamma")
g_male3_08 <-fitdist(cd_male3_08, "gamma")
g_male4_08 <-fitdist(cd_male4_08, "gamma")

cd_fem1_10 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, year == 2010) %>%
  pull(volajohdia)
cd_fem2_10 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, year == 2010) %>%
  pull(volajohdia)
cd_fem3_10 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, year == 2010) %>%
  pull(volajohdia)
cd_fem4_10 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, year == 2010) %>%
  pull(volajohdia)

cd_male1_10 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1, year == 2010) %>%
  pull(volajohdia)
cd_male2_10 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2, year == 2010) %>%
  pull(volajohdia)
cd_male3_10 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3, year == 2010) %>%
  pull(volajohdia)
cd_male4_10 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4, year == 2010) %>%
  pull(volajohdia)


g_fem1_10 <- fitdist(cd_fem1_10, "gamma")
g_fem2_10 <- fitdist(cd_fem2_10, "gamma")
g_fem3_10 <-fitdist(cd_fem3_10, "gamma")
g_fem4_10 <-fitdist(cd_fem4_10, "gamma")

g_male1_10 <- fitdist(cd_male1_10, "gamma")
g_male2_10 <- fitdist(cd_male2_10, "gamma")
g_male3_10 <-fitdist(cd_male3_10, "gamma")
g_male4_10 <-fitdist(cd_male4_10, "gamma")

cd_fem1_12 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, year == 2012) %>%
  pull(volajohdia)
cd_fem2_12 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, year == 2012) %>%
  pull(volajohdia)
cd_fem3_12 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, year == 2012) %>%
  pull(volajohdia)
cd_fem4_12 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, year == 2012) %>%
  pull(volajohdia)

cd_male1_12 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1, year == 2012) %>%
  pull(volajohdia)
cd_male2_12 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2, year == 2012) %>%
  pull(volajohdia)
cd_male3_12 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3, year == 2012) %>%
  pull(volajohdia)
cd_male4_12 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4, year == 2012) %>%
  pull(volajohdia)


g_fem1_12 <- fitdist(cd_fem1_12, "gamma")
g_fem2_12 <- fitdist(cd_fem2_12, "gamma")
g_fem3_12 <-fitdist(cd_fem3_12, "gamma")
g_fem4_12 <-fitdist(cd_fem4_12, "gamma")

g_male1_12 <- fitdist(cd_male1_12, "gamma")
g_male2_12 <- fitdist(cd_male2_12, "gamma")
g_male3_12 <-fitdist(cd_male3_12, "gamma")
g_male4_12 <-fitdist(cd_male4_12, "gamma")

cd_fem1_14 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, year == 2014) %>%
  pull(volajohdia)
cd_fem2_14 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, year == 2014) %>%
  pull(volajohdia)
cd_fem3_14 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, year == 2014) %>%
  pull(volajohdia)
cd_fem4_14 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, year == 2014) %>%
  pull(volajohdia)

cd_male1_14 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1, year == 2014) %>%
  pull(volajohdia)
cd_male2_14 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2, year == 2014) %>%
  pull(volajohdia)
cd_male3_14 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3, year == 2014) %>%
  pull(volajohdia)
cd_male4_14 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4, year == 2014) %>%
  pull(volajohdia)


g_fem1_14 <- fitdist(cd_fem1_14, "gamma")
g_fem2_14 <- fitdist(cd_fem2_14, "gamma")
g_fem3_14 <-fitdist(cd_fem3_14, "gamma")
g_fem4_14 <-fitdist(cd_fem4_14, "gamma")

g_male1_14 <- fitdist(cd_male1_14, "gamma")
g_male2_14 <- fitdist(cd_male2_14, "gamma")
g_male3_14 <-fitdist(cd_male3_14, "gamma")
g_male4_14 <-fitdist(cd_male4_14, "gamma")

cd_fem1_16 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, year == 2016) %>%
  pull(volajohdia)
cd_fem2_16 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, year == 2016) %>%
  pull(volajohdia)
cd_fem3_16 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, year == 2016) %>%
  pull(volajohdia)
cd_fem4_16 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, year == 2016) %>%
  pull(volajohdia)

cd_male1_16 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1, year == 2016) %>%
  pull(volajohdia)
cd_male2_16 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2, year == 2016) %>%
  pull(volajohdia)
cd_male3_16 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3, year == 2016) %>%
  pull(volajohdia)
cd_male4_16 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4, year == 2016) %>%
  pull(volajohdia)

g_fem1_16 <- fitdist(cd_fem1_16, "gamma")
g_fem2_16 <- fitdist(cd_fem2_16, "gamma")
g_fem3_16 <-fitdist(cd_fem3_16, "gamma")
g_fem4_16 <-fitdist(cd_fem4_16, "gamma")

g_male1_16 <- fitdist(cd_male1_16, "gamma")
g_male2_16 <- fitdist(cd_male2_16, "gamma")
g_male3_16 <-fitdist(cd_male3_16, "gamma")
g_male4_16 <-fitdist(cd_male4_16, "gamma")

cd_fem1_18 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, year == 2018) %>%
  pull(volajohdia)
cd_fem2_18 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, year == 2018) %>%
  pull(volajohdia)
cd_fem3_18 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, year == 2018) %>%
  pull(volajohdia)
cd_fem4_18 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, year == 2018) %>%
  pull(volajohdia)

cd_male1_18 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1, year == 2018) %>%
  pull(volajohdia)
cd_male2_18 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2, year == 2018) %>%
  pull(volajohdia)
cd_male3_18 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3, year == 2018) %>%
  pull(volajohdia)
cd_male4_18 <- data_input %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4, year == 2018) %>%
  pull(volajohdia)

g_fem1_18 <- fitdist(cd_fem1_18, "gamma")
g_fem2_18 <- fitdist(cd_fem2_18, "gamma")
g_fem3_18 <-fitdist(cd_fem3_18, "gamma")
g_fem4_18 <-fitdist(cd_fem4_18, "gamma")

g_male1_18 <- fitdist(cd_male1_18, "gamma")
g_male2_18 <- fitdist(cd_male2_18, "gamma")
g_male3_18 <-fitdist(cd_male3_18, "gamma")
g_male4_18 <-fitdist(cd_male4_18, "gamma")

# Listas de los parámetros de gamma y proporciones para cada año y grupo
gammas <- list(
  list(g_fem1_08, g_fem2_08, g_fem3_08, g_fem4_08),
  list(g_fem1_10, g_fem2_10, g_fem3_10, g_fem4_10),
  list(g_fem1_12, g_fem2_12, g_fem3_12, g_fem4_12),
  list(g_fem1_14, g_fem2_14, g_fem3_14, g_fem4_14),
  list(g_fem1_16, g_fem2_16, g_fem3_16, g_fem4_16),
  list(g_fem1_18, g_fem2_18, g_fem3_18, g_fem4_18)
)

p_abs_list <- list(
  list(input_female[[2,5]], input_female[[10,5]], input_female[[18,5]], input_female[[26,5]], input_female[[34,5]], input_female[[42,5]]),
  list(input_female[[4,5]], input_female[[12,5]], input_female[[20,5]], input_female[[28,5]], input_female[[36,5]], input_female[[44,5]]),
  list(input_female[[6,5]], input_female[[14,5]], input_female[[22,5]], input_female[[30,5]], input_female[[38,5]], input_female[[46,5]]),
  list(input_female[[8,5]], input_female[[16,5]], input_female[[24,5]], input_female[[32,5]], input_female[[40,5]], input_female[[48,5]])
)

p_form_list <- list(
  list(input_female[[1,5]], input_female[[9,5]], input_female[[17,5]], input_female[[25,5]], input_female[[33,5]], input_female[[41,5]]),
  list(input_female[[3,5]], input_female[[11,5]], input_female[[19,5]], input_female[[27,5]], input_female[[35,5]], input_female[[43,5]]),
  list(input_female[[5,5]], input_female[[13,5]], input_female[[21,5]], input_female[[29,5]], input_female[[37,5]], input_female[[45,5]]),
  list(input_female[[7,5]], input_female[[15,5]], input_female[[23,5]], input_female[[31,5]], input_female[[39,5]], input_female[[47,5]])
)

# Listas de los parámetros de gamma y proporciones para cada año y grupo masculino
gammas_male <- list(
  list(g_male1_08, g_male2_08, g_male3_08, g_male4_08),
  list(g_male1_10, g_male2_10, g_male3_10, g_male4_10),
  list(g_male1_12, g_male2_12, g_male3_12, g_male4_12),
  list(g_male1_14, g_male2_14, g_male3_14, g_male4_14),
  list(g_male1_16, g_male2_16, g_male3_16, g_male4_16),
  list(g_male1_18, g_male2_18, g_male3_18, g_male4_18)
)

p_abs_list_male <- list(
  list(input_male[[2,5]], input_male[[10,5]], input_male[[18,5]], input_male[[26,5]], input_male[[34,5]], input_male[[42,5]]),
  list(input_male[[4,5]], input_male[[12,5]], input_male[[20,5]], input_male[[28,5]], input_male[[36,5]], input_male[[44,5]]),
  list(input_male[[6,5]], input_male[[14,5]], input_male[[22,5]], input_male[[30,5]], input_male[[38,5]], input_male[[46,5]]),
  list(input_male[[8,5]], input_male[[16,5]], input_male[[24,5]], input_male[[32,5]], input_male[[40,5]], input_male[[48,5]])
)

p_form_list_male <- list(
  list(input_male[[1,5]], input_male[[9,5]], input_male[[17,5]], input_male[[25,5]], input_male[[33,5]], input_male[[41,5]]),
  list(input_male[[3,5]], input_male[[11,5]], input_male[[19,5]], input_male[[27,5]], input_male[[35,5]], input_male[[43,5]]),
  list(input_male[[5,5]], input_male[[13,5]], input_male[[21,5]], input_male[[29,5]], input_male[[37,5]], input_male[[45,5]]),
  list(input_male[[7,5]], input_male[[15,5]], input_male[[23,5]], input_male[[31,5]], input_male[[39,5]], input_male[[47,5]])
)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


####################################
# ESTIMATING AAF FOR BREAST CANCER #
####################################
b1_bcan <- 0.01018
var_bcan <- 0.0000007208
bcan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)
for (i in 1:length(bcan_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_bcan <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b1_bcan, 
                  var_beta = var_bcan, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_bcan) && !is.null(result_bcan$Point_Estimate)) {
      bcan_female[i, paste0("Fem", j, "_point")] <- result_bcan$Point_Estimate
      bcan_female[i, paste0("Fem", j, "_lower")] <- result_bcan$Lower_CI
      bcan_female[i, paste0("Fem", j, "_upper")] <- result_bcan$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", bcan_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}


# LIP AND ORAL CAVITY CANCER (BOTH SEXES)
b1_locan <- 0.02474
b2_locan <- -0.00004
rr_locan_fd <- 1.2
cov_matrix_locan <- matrix(c(0.000002953, -0.0000000127,
                             -0.0000000127, 0.000000000102), 
                           nrow = 2, ncol = 2, byrow = TRUE)

rr_locan_fun <- function(x, betas){
  # Example: assuming a quadratic model
  rr <- exp(betas[1] * x + betas[2] * x^2)
  return(rr)
}

# Define a data frame to store results
locan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Bucle sobre cada año y cada grupo femenino
for (i in 1:length(locan_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result <- confint_paf_vcov(gamma = gammas[[i]][[j]], 
                               betas = c(b1_locan, b2_locan), 
                               cov_matrix = cov_matrix_locan, 
                               p_abs = p_abs_list[[j]][[i]], 
                               p_form = p_form_list[[j]][[i]], 
                               rr_fd = rr_locan_fd, 
                               rr_function = rr_locan_fun)
    
    # Guarda los resultados en el data frame
    locan_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
    locan_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    locan_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}


# Display the results matrix
locan_female

# Define un data frame para almacenar los resultados para hombres
locan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(locan_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result_male <- confint_paf_vcov(gamma = gammas_male[[i]][[j]], 
                                    betas = c(b1_locan, b2_locan), 
                                    cov_matrix = cov_matrix_locan, 
                                    p_abs = p_abs_list_male[[j]][[i]], 
                                    p_form = p_form_list_male[[j]][[i]], 
                                    rr_fd = rr_locan_fd, 
                                    rr_function = rr_locan_fun)
    
    # Guarda los resultados en el data frame
    locan_male[i, paste0("Male", j, "_point")] <- result_male$point_estimate
    locan_male[i, paste0("Male", j, "_lower")] <- result_male$lower_ci
    locan_male[i, paste0("Male", j, "_upper")] <- result_male$upper_ci
  }
}


################################################
# ESTIMATING AAF FOR OTHER PHARINGEAL CANCER   #
################################################

# OTHER PHARINGEAL CANCER (BOTH SEXES)
b1_opcan <- 0.02474
b2_opcan <- -0.00004
rr_opcan_fd <- 1.2
cov_matrix_opcan <- matrix(c(0.000002953, -0.0000000127,
                             -0.0000000127, 0.000000000102), 
                           nrow = 2, ncol = 2, byrow = TRUE)


# Define a data frame to store results
opcan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Bucle sobre cada año y cada grupo femenino
for (i in 1:length(opcan_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result <- confint_paf_vcov(gamma = gammas[[i]][[j]], 
                               betas = c(b1_opcan, b2_opcan), 
                               cov_matrix = cov_matrix_opcan, 
                               p_abs = p_abs_list[[j]][[i]], 
                               p_form = p_form_list[[j]][[i]], 
                               rr_fd = rr_opcan_fd, 
                               rr_function = rr_locan_fun)
    
    # Guarda los resultados en el data frame
    opcan_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
    opcan_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    opcan_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}


# Display the results matrix
opcan_female

# Define un data frame para almacenar los resultados para hombres
opcan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(opcan_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result_male <- confint_paf_vcov(gamma = gammas_male[[i]][[j]], 
                                    betas = c(b1_opcan, b2_opcan), 
                                    cov_matrix = cov_matrix_opcan, 
                                    p_abs = p_abs_list_male[[j]][[i]], 
                                    p_form = p_form_list_male[[j]][[i]], 
                                    rr_fd = rr_opcan_fd, 
                                    rr_function = rr_locan_fun)
    
    # Guarda los resultados en el data frame
    opcan_male[i, paste0("Male", j, "_point")] <- result_male$point_estimate
    opcan_male[i, paste0("Male", j, "_lower")] <- result_male$lower_ci
    opcan_male[i, paste0("Male", j, "_upper")] <- result_male$upper_ci
  }
}


##########################################
# ESTIMATING AAF FOR OESOPHAFUS CANCER   #
##########################################
rr_oescan_cal <- function(b1,b2,x){
  exp(b1*x+b2*x**3)
}
b1_oescan <- 0.0132063596418668
b2_oescan <- -4.14801974664481*10^-08
rr_oescan_fd <- 1.16
cov_matrix_oescan <- matrix(c(1.525706e-07, -6.885205e-13,
                              -6.885205e-13, 2.953e-06), 
                            nrow = 2, ncol = 2, byrow = TRUE)

# Define a data frame to store results
oescan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Bucle sobre cada año y cada grupo femenino
for (i in 1:length(oescan_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result <- confint_paf_vcov(gamma = gammas[[i]][[j]], 
                               betas = c(b1_oescan, b2_oescan), 
                               cov_matrix = cov_matrix_opcan, 
                               p_abs = p_abs_list[[j]][[i]], 
                               p_form = p_form_list[[j]][[i]], 
                               rr_fd = rr_oescan_fd, 
                               rr_function = rr_locan_fun)
    
    # Guarda los resultados en el data frame
    oescan_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
    oescan_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    oescan_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}


# Display the results matrix
oescan_female

# Define un data frame para almacenar los resultados para hombres
oescan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(oescan_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result_male <- confint_paf_vcov(gamma = gammas_male[[i]][[j]], 
                                    betas = c(b1_oescan, b2_oescan), 
                                    cov_matrix = cov_matrix_opcan, 
                                    p_abs = p_abs_list_male[[j]][[i]], 
                                    p_form = p_form_list_male[[j]][[i]], 
                                    rr_fd = rr_oescan_fd, 
                                    rr_function = rr_locan_fun)
    
    # Guarda los resultados en el data frame
    oescan_male[i, paste0("Male", j, "_point")] <- result_male$point_estimate
    oescan_male[i, paste0("Male", j, "_lower")] <- result_male$lower_ci
    oescan_male[i, paste0("Male", j, "_upper")] <- result_male$upper_ci
  }
}

################################################
# ESTIMATING AAF FOR COLON AND RECTUM CANCER   #
################################################
b1_crcan_fem <- 0.006279
rr_crcan_fd_fem <- 1.05
var_rr_crcan_fd_fem <- 0.1459680025873172
b1_crcan_male <- 0.006279
var_crcan <- 0.000000907
rr_crcan_fd_male <- 2.19

crcan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)
for (i in 1:length(crcan_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_bcan <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b1_crcan_fem, 
                  var_beta = var_crcan, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = rr_crcan_fd_fem, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_bcan) && !is.null(result_bcan$Point_Estimate)) {
      crcan_female[i, paste0("Fem", j, "_point")] <- result_bcan$Point_Estimate
      crcan_female[i, paste0("Fem", j, "_lower")] <- result_bcan$Lower_CI
      crcan_female[i, paste0("Fem", j, "_upper")] <- result_bcan$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", bcan_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}

crcan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(crcan_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_bcan <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b1_crcan_male, 
                  var_beta = var_crcan, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = rr_crcan_fd_male, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_bcan) && !is.null(result_bcan$Point_Estimate)) {
      crcan_male[i, paste0("Male", j, "_point")] <- result_bcan$Point_Estimate
      crcan_male[i, paste0("Male", j, "_lower")] <- result_bcan$Lower_CI
      crcan_male[i, paste0("Male", j, "_upper")] <- result_bcan$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", bcan_female$Year[i], "y grupo Male", j, "\n")
    }
  }
}

######################################
# ESTIMATING AAF FOR LIVER CANCER   #
#####################################
b1_lican <- 0.005041
var_lican <- 0.000003097
rr_lican_fd_male <- 2.23
rr_lican_fd_fem <- 2.68

lican_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)
for (i in 1:length(lican_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_lican <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b1_lican, 
                  var_beta = var_lican, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = rr_lican_fd_fem, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_bcan) && !is.null(result_bcan$Point_Estimate)) {
      lican_female[i, paste0("Fem", j, "_point")] <- result_lican$Point_Estimate
      lican_female[i, paste0("Fem", j, "_lower")] <- result_lican$Lower_CI
      lican_female[i, paste0("Fem", j, "_upper")] <- result_lican$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", lican_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}

lican_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(lican_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_lican <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b1_lican, 
                  var_beta = var_lican, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = rr_lican_fd_male, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_bcan) && !is.null(result_bcan$Point_Estimate)) {
     lican_male[i, paste0("Male", j, "_point")] <- result_lican$Point_Estimate
      lican_male[i, paste0("Male", j, "_lower")] <- result_lican$Lower_CI
      lican_male[i, paste0("Male", j, "_upper")] <- result_lican$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", lican_female$Year[i], "y grupo Male", j, "\n")
    }
  }
}

######################################
# ESTIMATING AAF FOR LARYNX CANCER   #
#####################################
b1_lxcan <- 0.01462
b2_lxcan <- -0.00002
cov_matrix_lxcan <- matrix(c(3.585e-06, -1.62e-08,
                             -1.62e-08, 1.26e-10), 
                           nrow = 2, ncol = 2, byrow = TRUE)

rr_fd_lxcan <- 1.18
rr_lxcan_fun <- function(x, betas){
  exp(x * betas[1] + x^2 * betas[2])
}

# Define a data frame to store results
lxcan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Bucle sobre cada año y cada grupo femenino
for (i in 1:length(lxcan_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result <- confint_paf_vcov(gamma = gammas[[i]][[j]], 
                               betas = c(b1_lxcan, b2_lxcan), 
                               cov_matrix = cov_matrix_lxcan, 
                               p_abs = p_abs_list[[j]][[i]], 
                               p_form = p_form_list[[j]][[i]], 
                               rr_fd = rr_fd_lxcan, 
                               rr_function = rr_lxcan_fun)
    
    # Guarda los resultados en el data frame
    lxcan_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
    lxcan_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    lxcan_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}

# Define un data frame para almacenar los resultados para hombres
lxcan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(lxcan_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result_male <- confint_paf_vcov(gamma = gammas_male[[i]][[j]], 
                                    betas = c(b1_lxcan, b2_lxcan), 
                                    cov_matrix = cov_matrix_lxcan, 
                                    p_abs = p_abs_list_male[[j]][[i]], 
                                    p_form = p_form_list_male[[j]][[i]], 
                                    rr_fd = rr_fd_lxcan, 
                                    rr_function = rr_lxcan_fun)
    
    # Guarda los resultados en el data frame
    lxcan_male[i, paste0("Male", j, "_point")] <- result_male$point_estimate
    lxcan_male[i, paste0("Male", j, "_lower")] <- result_male$lower_ci
    lxcan_male[i, paste0("Male", j, "_upper")] <- result_male$upper_ci
  }
}


