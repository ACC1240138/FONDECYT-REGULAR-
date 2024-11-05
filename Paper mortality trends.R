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
  dplyr::select(year, sexo, exp, edad_tramo, volajohdia, cvolaj,dias_binge)


input <- data_input %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo,year, edad_tramo, cvolaj) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>% 
  filter(cvolaj == "fd" | cvolaj == "ltabs") %>% 
  dplyr::select(-weighted_count)

data_input %>% 
  mutate(hed = ifelse(dias_binge > 0, 1 , 0)) %>% 
  filter(!is.na(hed)) %>% 
  group_by(year, sexo, edad_tramo) %>% 
  count(hed) %>% 
  mutate(prop_hed = n/sum(n))

data_hed %>% 
  group_by(year, hed, sexo) %>% 
  count(aux = volajohdia > 0) %>% 
  filter(aux == TRUE, year == 2018)

input_male <- input %>% 
  filter(sexo == "Hombre")
input_female <- input %>% 
  filter(sexo == "Mujer")

# Function to fit gamma distribution
# Define the x_vals range (adjust as needed)
x_vals <- seq(0.1, 150, length.out = 1500)

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

# FOR ISCHAEMIC STROKE AND INJURIES WE NEED THE PROPORTION OF CURRENT DRINKERS WHO
# ENGAGE IN HEAVY EPISODIC DRINKING NHED and HED
paf_hed_function <- function(x_60,x_150, y_nhed, y_hed_60, 
                             y_hed_150, rr_nhed, rr_hed_60,
                             rr_hed_150, rr_form,p_abs, 
                             p_form, p_hed) {
  
  trap_int_hed <- function(x, y, rr, prop_abs, rr_form, prop_form, p_hed) {
    
    dx <- x[2] - x[1]
    ncgamma <- sum((y[-1] + y[-length(y)]) / 2) * dx
    
    normalized_y <- ((1 - (prop_abs + prop_form)) * p_hed) * y / ncgamma
    
    excess_rr <- rr - 1
    
    weighted_excess_rr <- normalized_y * excess_rr
    
    numerator <- (rr_form - 1) * prop_form + sum((weighted_excess_rr[-1] + weighted_excess_rr[-length(weighted_excess_rr)]) / 2) * dx
    
    denominator <- numerator + 1
    
    paf <- round(numerator / denominator, 3)
    
    return(paf)
  }
  
  
  int_ri_nhed <- trap_int_hed(
    x = x_60, 
    y = y_nhed, 
    rr = rr_nhed, 
    prop_abs = p_abs, 
    rr_form = rr_form, 
    prop_form = p_form, 
    p_hed = 1 - p_hed
  )
  
  int_ri_hed1 <- trap_int_hed(
    x = x_60, 
    y = y_hed_60, 
    rr = rr_hed_60, 
    prop_abs = p_abs, 
    rr_form = rr_form, 
    prop_form = p_form, 
    p_hed = p_hed
  )
  
  int_ri_hed2 <- trap_int_hed(
    x = x_150, 
    y = y_hed_150, 
    rr = rr_hed_150, 
    prop_abs = p_abs, 
    rr_form = rr_form, 
    prop_form = p_form, 
    p_hed = p_hed
  )
  
  num <- (int_ri_nhed + int_ri_hed1 + int_ri_hed2)
  den <- num + 1
  paf_ri_fem1 <- num / den
  
  return(paf_ri_fem1)
}

confint_paf_hed <- function(gammas, betas, cov_matrix, p_abs, p_form, rr_fd, rr_function_nhed,
                            rr_function_hed, p_hed) {
  set.seed(145)
  n_sim <- 10000
  simulated_pafs <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Simulate PCA mean and SD using the gamma distribution
    pca_sim_nhed <- rgamma(1000, shape = gammas[[1]]$estimate["shape"], 
                           rate = gammas[[1]]$estimate["rate"])
    pca_sim_hed <- rgamma(1000, shape = gammas[[2]]$estimate["shape"], 
                          rate = gammas[[2]]$estimate["rate"])
    
    # Calculate the shape and rate parameters from the simulated PCA values
    mean_sim_nhed <- mean(pca_sim_nhed)
    sd_sim_nhed <- sd(pca_sim_nhed)
    mean_sim_hed <- mean(pca_sim_hed)
    sd_sim_hed <- sd(pca_sim_hed)
    
    # Calculate shape and rate for the new gamma distribution
    shape_sim_nhed <- (mean_sim_nhed / sd_sim_nhed)^2
    rate_sim_nhed <- mean_sim_nhed / (sd_sim_nhed^2)
    shape_sim_hed <- (mean_sim_hed / sd_sim_hed)^2
    rate_sim_hed <- mean_sim_hed / (sd_sim_hed^2)
    
    # Simulate the gamma distribution
    y_gamma_sim_nhed <- dgamma(x_vals_nhed, shape = shape_sim_nhed, rate = rate_sim_nhed)
    y_gamma_sim_hed_60 <- dgamma(x_vals_nhed, shape = shape_sim_hed, rate = rate_sim_hed)
    y_gamma_sim_hed_150 <- dgamma(x_vals_hed, shape = shape_sim_hed, rate = rate_sim_hed)
    
    # Simulate the beta coefficients jointly from a multivariate normal distribution
    beta_sim <- MASS::mvrnorm(1, mu = betas, Sigma = cov_matrix)
    
    # Compute the relative risk based on the simulated betas
    rr_sim_nhed <- rr_function_nhed(x = x_vals_nhed, beta = beta_sim[1])
    rr_sim_hed_60 <- rr_function_hed(x = x_vals_nhed, betas = beta_sim)
    rr_sim_hed_150 <- rr_function_hed(x = x_vals_hed, betas = beta_sim)
    
    # Simulate proportions of lifetime abstainers and former drinkers
    prop_abs_sim <- max(rnorm(1000, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000)), 0.001)
    prop_form_sim <- max(rnorm(1000, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000)), 0.001)
    prop_hed_sim <- max(rnorm(1000, mean = p_hed, sd = sqrt(p_hed * (1 - p_hed) / 1000)), 0.001)
    
    # Calculate PAF using the customized function
    simulated_pafs[i] <- paf_hed_function(
      x_60 = x_vals_nhed, 
      x_150 = x_vals_hed, 
      y_nhed = y_gamma_sim_nhed, 
      y_hed_60 = y_gamma_sim_hed_60, 
      y_hed_150 = y_gamma_sim_hed_150, 
      rr_nhed = rr_sim_nhed, 
      rr_hed_60 = rr_sim_hed_60, 
      rr_hed_150 = rr_sim_hed_150,
      rr_form = rr_fd, 
      p_abs = prop_abs_sim, 
      p_form = prop_form_sim, 
      p_hed = prop_hed_sim
    )
  }
  
  # Remove NaN values from simulated PAFs
  simulated_pafs <- simulated_pafs[!is.nan(simulated_pafs)]
  
  if (length(simulated_pafs) == 0) {
    stop("All simulations resulted in NaN values. Please check your input parameters.")
  }
  
  # Calculate the 95% confidence interval
  paf_lower <- quantile(simulated_pafs, 0.025)
  paf_upper <- quantile(simulated_pafs, 0.975)
  paf_point_estimate <- mean(simulated_pafs)
  
  return(list(point_estimate = round(paf_point_estimate,3), lower_ci = round(paf_lower,3), upper_ci = round(paf_upper,3)))
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

cd_fem1_nhed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, hed == 0, year == 2008) %>%
  pull(volajohdia)
cd_fem1_hed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, hed == 1, year == 2008) %>%
  pull(volajohdia)

cd_fem2_nhed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, hed == 0, year == 2008) %>%
  pull(volajohdia)
cd_fem2_hed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, hed == 1, year == 2008) %>%
  pull(volajohdia)

cd_fem3_nhed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, hed == 0, year == 2008) %>%
  pull(volajohdia)
cd_fem3_hed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, hed == 1, year == 2008) %>%
  pull(volajohdia)

cd_fem4_nhed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, hed == 0, year == 2008) %>%
  pull(volajohdia)
cd_fem4_hed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, hed == 1, year == 2008) %>%
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

# TRAMO 1
cd_male1_nhed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1, hed == 0, year == 2008) %>%
  pull(volajohdia)
cd_male1_hed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 1, hed == 1, year == 2008) %>%
  pull(volajohdia)

# TRAMO 2
cd_male2_nhed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2, hed == 0, year == 2008) %>%
  pull(volajohdia)
cd_male2_hed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 2, hed == 1, year == 2008) %>%
  pull(volajohdia)

# TRAMO 3
cd_male3_nhed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3, hed == 0, year == 2008) %>%
  pull(volajohdia)
cd_male3_hed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 3, hed == 1, year == 2008) %>%
  pull(volajohdia)

# TRAMO 4
cd_male4_nhed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4, hed == 0, year == 2008) %>%
  pull(volajohdia)
cd_male4_hed_08 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Hombre", edad_tramo == 4, hed == 1, year == 2008) %>%
  pull(volajohdia)


g_fem1_08 <- fitdist(cd_fem1_08, "gamma")
g_fem2_08 <- fitdist(cd_fem2_08, "gamma")
g_fem3_08 <-fitdist(cd_fem3_08, "gamma")
g_fem4_08 <-fitdist(cd_fem4_08, "gamma")

g_fem1_nhed_08 <- fitdist(cd_fem1_nhed_08, "gamma")
g_fem2_nhed_08 <- fitdist(cd_fem2_nhed_08, "gamma")
g_fem3_nhed_08 <-fitdist(cd_fem3_nhed_08, "gamma")
g_fem4_nhed_08 <-fitdist(cd_fem4_nhed_08, "gamma")
g_fem1_hed_08 <- fitdist(cd_fem1_hed_08, "gamma")
g_fem2_hed_08 <- fitdist(cd_fem2_hed_08, "gamma")
g_fem3_hed_08 <-fitdist(cd_fem3_hed_08, "gamma")
g_fem4_hed_08 <-fitdist(cd_fem4_hed_08, "gamma")

g_male1_08 <- fitdist(cd_male1_08, "gamma")
g_male2_08 <- fitdist(cd_male2_08, "gamma")
g_male3_08 <-fitdist(cd_male3_08, "gamma")
g_male4_08 <-fitdist(cd_male4_08, "gamma")

g_male1_nhed_08 <- fitdist(cd_male1_nhed_08, "gamma")
g_male2_nhed_08 <- fitdist(cd_male2_nhed_08, "gamma")
g_male3_nhed_08 <-fitdist(cd_male3_nhed_08, "gamma")
g_male4_nhed_08 <-fitdist(cd_male4_nhed_08, "gamma")
g_male1_hed_08 <- fitdist(cd_male1_hed_08, "gamma")
g_male2_hed_08 <- fitdist(cd_male2_hed_08, "gamma")
g_male3_hed_08 <-fitdist(cd_male3_hed_08, "gamma")
g_male3_hed_08 <-fitdist(cd_male4_hed_08, "gamma")

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

cd_fem1_nhed_10 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, hed == 0, year == 2010) %>%
  pull(volajohdia)
cd_fem1_hed_10 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 1, hed == 1, year == 2010) %>%
  pull(volajohdia)

cd_fem2_nhed_10 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, hed == 0, year == 2010) %>%
  pull(volajohdia)
cd_fem2_hed_10 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 2, hed == 1, year == 2010) %>%
  pull(volajohdia)

cd_fem3_nhed_10 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, hed == 0, year == 2010) %>%
  pull(volajohdia)
cd_fem3_hed_10 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 3, hed == 1, year == 2010) %>%
  pull(volajohdia)

cd_fem4_nhed_10 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, hed == 0, year == 2010) %>%
  pull(volajohdia)
cd_fem4_hed_10 <- data_hed %>% 
  filter(volajohdia > 0, sexo == "Mujer", edad_tramo == 4, hed == 1, year == 2010) %>%
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

# Define years and age groups
years <- c(2008, 2010, 2012, 2014, 2016, 2018)
age_groups <- 1:4

# Initialize a list to store the results
cd_male_hed <- list()

# Loop through each year and age group to generate the data
for (year in years) {
  year_data <- list()
  for (age_group in age_groups) {
    # Filter data for non-heavy episodic drinkers (nhed) for each year and age group
    nhed_data <- data_hed %>% 
      filter(volajohdia > 0, sexo == "Hombre", edad_tramo == age_group, hed == 0, year == year) %>%
      pull(volajohdia)
    
    # Filter data for heavy episodic drinkers (hed) for each year and age group
    hed_data <- data_hed %>% 
      filter(volajohdia > 0, sexo == "Hombre", edad_tramo == age_group, hed == 1, year == year) %>%
      pull(volajohdia)
    
    # Store the results in the list for the given year and age group
    year_data[[paste0("tramo", age_group, "_nhed")]] <- nhed_data
    year_data[[paste0("tramo", age_group, "_hed")]] <- hed_data
  }
  cd_male_hed[[as.character(year)]] <- year_data
}

# Initialize a list to store gamma fits for each year, age group, and `hed` status
gamma_fits_male <- list()

# Loop through each year and age group to fit the gamma distribution
for (year in names(cd_male_hed)) {
  year_fits <- list()
  for (age_group in age_groups) {
    # Retrieve nhed data for the given year and age group
    nhed_data <- cd_male_hed[[year]][[paste0("tramo", age_group, "_nhed")]]
    if (length(nhed_data) > 1) { # Check if there is enough data to fit a distribution
      year_fits[[paste0("tramo", age_group, "_nhed")]] <- fitdist(nhed_data, "gamma")
    } else {
      year_fits[[paste0("tramo", age_group, "_nhed")]] <- NA  # Not enough data
    }
    
    # Retrieve hed data for the given year and age group
    hed_data <- cd_male_hed[[year]][[paste0("tramo", age_group, "_hed")]]
    if (length(hed_data) > 1) {
      year_fits[[paste0("tramo", age_group, "_hed")]] <- fitdist(hed_data, "gamma")
    } else {
      year_fits[[paste0("tramo", age_group, "_hed")]] <- NA
    }
  }
  gamma_fits_male[[year]] <- year_fits
}

g_fem1_10 <- fitdist(cd_fem1_10, "gamma")
g_fem2_10 <- fitdist(cd_fem2_10, "gamma")
g_fem3_10 <-fitdist(cd_fem3_10, "gamma")
g_fem4_10 <-fitdist(cd_fem4_10, "gamma")

g_fem1_nhed_10 <- fitdist(cd_fem1_nhed_10, "gamma")
g_fem2_nhed_10 <- fitdist(cd_fem2_nhed_10, "gamma")
g_fem3_nhed_10 <-fitdist(cd_fem3_nhed_10, "gamma")
g_fem4_nhed_10 <-fitdist(cd_fem4_nhed_10, "gamma")
g_fem1_hed_10 <- fitdist(cd_fem1_hed_10, "gamma")
g_fem2_hed_10 <- fitdist(cd_fem2_hed_10, "gamma")
g_fem3_hed_10 <-fitdist(cd_fem3_hed_10, "gamma")
g_fem4_hed_10 <-fitdist(cd_fem4_hed_10, "gamma")

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

cd_fem1_nhed_12 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 1, hed == 0, year == 2012) %>%
  pull(volajohdia)
cd_fem1_hed_12 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 1, hed == 1, year == 2012) %>%
  pull(volajohdia)

cd_fem2_nhed_12 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 2, hed == 0, year == 2012) %>%
  pull(volajohdia)
cd_fem2_hed_12 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 2, hed == 1, year == 2012) %>%
  pull(volajohdia)

cd_fem3_nhed_12 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 3, hed == 0, year == 2012) %>%
  pull(volajohdia)
cd_fem3_hed_12 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 3, hed == 1, year == 2012) %>%
  pull(volajohdia)

cd_fem4_nhed_12 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 4, hed == 0, year == 2012) %>%
  pull(volajohdia)
cd_fem4_hed_12 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 4, hed == 1, year == 2012) %>%
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

g_fem1_nhed_12 <- fitdist(cd_fem1_nhed_12, "gamma")
g_fem2_nhed_12 <- fitdist(cd_fem2_nhed_12, "gamma")
g_fem3_nhed_12 <-fitdist(cd_fem3_nhed_12, "gamma")
g_fem4_nhed_12 <-fitdist(cd_fem4_nhed_12, "gamma")
g_fem1_hed_12 <- fitdist(cd_fem1_hed_12, "gamma")
g_fem2_hed_12 <- fitdist(cd_fem2_hed_12, "gamma")
g_fem3_hed_12 <-fitdist(cd_fem3_hed_12, "gamma")
g_fem4_hed_12 <-fitdist(cd_fem4_hed_12, "gamma")

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

cd_fem1_nhed_14 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 1, hed == 0, year == 2014) %>%
  pull(volajohdia)
cd_fem1_hed_14 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 1, hed == 1, year == 2014) %>%
  pull(volajohdia)

cd_fem2_nhed_14 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 2, hed == 0, year == 2014) %>%
  pull(volajohdia)
cd_fem2_hed_14 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 2, hed == 1, year == 2014) %>%
  pull(volajohdia)

cd_fem3_nhed_14 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 3, hed == 0, year == 2014) %>%
  pull(volajohdia)
cd_fem3_hed_14 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 3, hed == 1, year == 2014) %>%
  pull(volajohdia)

cd_fem4_nhed_14 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 4, hed == 0, year == 2014) %>%
  pull(volajohdia)
cd_fem4_hed_14 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 4, hed == 1, year == 2014) %>%
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

g_fem1_nhed_14 <- fitdist(cd_fem1_nhed_14, "gamma")
g_fem2_nhed_14 <- fitdist(cd_fem2_nhed_14, "gamma")
g_fem3_nhed_14 <-fitdist(cd_fem3_nhed_14, "gamma")
g_fem4_nhed_14 <-fitdist(cd_fem4_nhed_14, "gamma")
g_fem1_hed_14 <- fitdist(cd_fem1_hed_14, "gamma")
g_fem2_hed_14 <- fitdist(cd_fem2_hed_14, "gamma")
g_fem3_hed_14 <-fitdist(cd_fem3_hed_14, "gamma")
g_fem4_hed_14 <-fitdist(cd_fem4_hed_14, "gamma")


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

cd_fem1_nhed_16 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 1, hed == 0, year == 2016) %>%
  pull(volajohdia)
cd_fem1_hed_16 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 1, hed == 1, year == 2016) %>%
  pull(volajohdia)

cd_fem2_nhed_16 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 2, hed == 0, year == 2016) %>%
  pull(volajohdia)
cd_fem2_hed_16 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 2, hed == 1, year == 2016) %>%
  pull(volajohdia)

cd_fem3_nhed_16 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 3, hed == 0, year == 2016) %>%
  pull(volajohdia)
cd_fem3_hed_16 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 3, hed == 1, year == 2016) %>%
  pull(volajohdia)

cd_fem4_nhed_16 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 4, hed == 0, year == 2016) %>%
  pull(volajohdia)
cd_fem4_hed_16 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 4, hed == 1, year == 2016) %>%
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

g_fem1_nhed_16 <- fitdist(cd_fem1_nhed_16, "gamma")
g_fem2_nhed_16 <- fitdist(cd_fem2_nhed_16, "gamma")
g_fem3_nhed_16 <-fitdist(cd_fem3_nhed_16, "gamma")
g_fem4_nhed_16 <-fitdist(cd_fem4_nhed_16, "gamma")
g_fem1_hed_16 <- fitdist(cd_fem1_hed_16, "gamma")
g_fem2_hed_16 <- fitdist(cd_fem2_hed_16, "gamma")
g_fem3_hed_16 <-fitdist(cd_fem3_hed_16, "gamma")
g_fem4_hed_16 <-fitdist(cd_fem4_hed_16, "gamma")

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

cd_fem1_nhed_18 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 1, hed == 0, year == 2018) %>%
  pull(volajohdia)
cd_fem1_hed_18 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 1, hed == 1, year == 2018) %>%
  pull(volajohdia)

cd_fem2_nhed_18 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 2, hed == 0, year == 2018) %>%
  pull(volajohdia)
cd_fem2_hed_18 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 2, hed == 1, year == 2018) %>%
  pull(volajohdia)

cd_fem3_nhed_18 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 3, hed == 0, year == 2018) %>%
  pull(volajohdia)
cd_fem3_hed_18 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 3, hed == 1, year == 2018) %>%
  pull(volajohdia)

cd_fem4_nhed_18 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 4, hed == 0, year == 2018) %>%
  pull(volajohdia)
cd_fem4_hed_18 <- data_hed %>% 
  filter(sexo == "Mujer", edad_tramo == 4, hed == 1, year == 2018) %>%
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

g_fem1_nhed_18 <- fitdist(cd_fem1_nhed_18, "gamma")
g_fem2_nhed_18 <- fitdist(cd_fem2_nhed_18, "gamma")
g_fem3_nhed_18 <-fitdist(cd_fem3_nhed_18, "gamma")
g_fem4_nhed_18 <-fitdist(cd_fem4_nhed_18, "gamma")
g_fem1_hed_18 <- fitdist(cd_fem1_hed_18, "gamma")
g_fem2_hed_18 <- fitdist(cd_fem2_hed_18, "gamma")
g_fem3_hed_18 <-fitdist(cd_fem3_hed_18, "gamma")
g_fem4_hed_18 <-fitdist(cd_fem4_hed_18, "gamma")

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

p_hed_list_male <- list(
  list(data_hed[[1,6]], data_hed[[9,6]], data_hed[[17,6]], data_hed[[25,6]],data_hed[[33,6]],data_hed[[41,6]]),
  list(data_hed[[2,6]], data_hed[[10,6]], data_hed[[18,6]], data_hed[[26,6]],data_hed[[34,6]],data_hed[[42,6]]),
  list(data_hed[[3,6]], data_hed[[11,6]], data_hed[[19,6]], data_hed[[27,6]],data_hed[[35,6]],data_hed[[43,6]]),
  list(data_hed[[4,6]], data_hed[[12,6]], data_hed[[20,6]], data_hed[[28,6]],data_hed[[36,6]],data_hed[[44,6]])
)
p_hed_list_fem <- list(
    list(data_hed[[5,6]], data_hed[[13,6]], data_hed[[21,6]], data_hed[[29,6]],data_hed[[37,6]],data_hed[[45,6]]),
    list(data_hed[[6,6]], data_hed[[14,6]], data_hed[[22,6]], data_hed[[30,6]],data_hed[[38,6]],data_hed[[46,6]]),
    list(data_hed[[7,6]], data_hed[[15,6]], data_hed[[23,6]], data_hed[[31,6]],data_hed[[39,6]],data_hed[[47,6]]),
    list(data_hed[[8,6]], data_hed[[16,6]], data_hed[[24,6]], data_hed[[32,6]],data_hed[[40,6]],data_hed[[48,6]])
  )

gammas_fem_hed <- list(
  # 2008 data
  list(
    list(g_fem1_nhed_08, g_fem1_hed_08),
    list(g_fem2_nhed_08, g_fem2_hed_08),
    list(g_fem3_nhed_08, g_fem3_hed_08),
    list(g_fem4_nhed_08, g_fem4_hed_08)
  ),
  # 2010 data
  list(
    list(g_fem1_nhed_10, g_fem1_hed_10),
    list(g_fem2_nhed_10, g_fem2_hed_10),
    list(g_fem3_nhed_10, g_fem3_hed_10),
    list(g_fem4_nhed_10, g_fem4_hed_10)
  ),
  # 2012 data
  list(
    list(g_fem1_nhed_12, g_fem1_hed_12),
    list(g_fem2_nhed_12, g_fem2_hed_12),
    list(g_fem3_nhed_12, g_fem3_hed_12),
    list(g_fem4_nhed_12, g_fem4_hed_12)
  ),
  # 2014 data
  list(
    list(g_fem1_nhed_14, g_fem1_hed_14),
    list(g_fem2_nhed_14, g_fem2_hed_14),
    list(g_fem3_nhed_14, g_fem3_hed_14),
    list(g_fem4_nhed_14, g_fem4_hed_14)
  ),
  # 2016 data
  list(
    list(g_fem1_nhed_16, g_fem1_hed_16),
    list(g_fem2_nhed_16, g_fem2_hed_16),
    list(g_fem3_nhed_16, g_fem3_hed_16),
    list(g_fem4_nhed_16, g_fem4_hed_16)
  ),
  # 2018 data
  list(
    list(g_fem1_nhed_18, g_fem1_hed_18),
    list(g_fem2_nhed_18, g_fem2_hed_18),
    list(g_fem3_nhed_18, g_fem3_hed_18),
    list(g_fem4_nhed_18, g_fem4_hed_18)
  )
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
                               cov_matrix = cov_matrix_oescan, 
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
                                    cov_matrix = cov_matrix_oescan, 
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
      confint_paf(gamma = gammas_male[[i]][[j]], 
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
    if (!is.null(lican_female) && !is.null(lican_female$Point_Estimate)) {
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
      confint_paf(gamma = gammas_male[[i]][[j]], 
                  beta = b1_lican, 
                  var_beta = var_lican, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = rr_lican_fd_male, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_lican) && !is.null(result_lican$Point_Estimate)) {
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

######################
# NEURO-PSYCHIATRIC #
#####################
# EPILEPSY (BOTH SEXES)
b1_epi <- 1.22861
var_epi <- 0.1391974
rr_epi_fun <- function(x, beta){
  exp(beta * (x+0.5)/100)
}

epi_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

for (i in 1:length(epi_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_epi <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b1_epi, 
                  var_beta = var_epi, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_epi_fun)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_bcan) && !is.null(result_bcan$Point_Estimate)) {
      epi_female[i, paste0("Fem", j, "_point")] <- result_epi$Point_Estimate
      epi_female[i, paste0("Fem", j, "_lower")] <- result_epi$Lower_CI
      epi_female[i, paste0("Fem", j, "_upper")] <- result_epi$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", epi_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}

epi_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(epi_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_epi <- tryCatch({
      confint_paf(gamma = gammas_male[[i]][[j]], 
                  beta = b1_epi, 
                  var_beta = var_epi, 
                  p_abs = p_abs_list_male[[j]][[i]], 
                  p_form = p_form_list_male[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_epi_fun)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_epi) && !is.null(result_epi$Point_Estimate)) {
      epi_male[i, paste0("Male", j, "_point")] <- result_epi$Point_Estimate
      epi_male[i, paste0("Male", j, "_lower")] <- result_epi$Lower_CI
      epi_male[i, paste0("Male", j, "_upper")] <- result_epi$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", epi_male$Year[i], "y grupo Male", j, "\n")
    }
  }
}


################
# OTHER CAUSES #
################

# Diabetes Mellitus (Male)
b1_dm_male <- 0.1763703
b2_dm_male <- -0.0728256
vcov_diabetes_male <- matrix(c(0.1681525, -0.2240129,
                               -0.2240129, 0.29964479),
                             nrow = 2, ncol = 2, byrow = TRUE)


rr_diabetes_male_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  
  # Compute the relative risk
  rr <- exp((x / 100)^2 * b1 + (x / 100)^3 * b2)
  
  return(rr)
}

dm_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(dm_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result_male <- confint_paf_vcov(gamma = gammas_male[[i]][[j]], 
                                    betas = c(b1_dm_male, b2_dm_male), 
                                    cov_matrix = vcov_diabetes_male, 
                                    p_abs = p_abs_list_male[[j]][[i]], 
                                    p_form = p_form_list_male[[j]][[i]], 
                                    rr_fd = 1.18, 
                                    rr_function = rr_diabetes_male_fun)
    
    # Guarda los resultados en el data frame
    dm_male[i, paste0("Male", j, "_point")] <- result_male$point_estimate
    dm_male[i, paste0("Male", j, "_lower")] <- result_male$lower_ci
    dm_male[i, paste0("Male", j, "_upper")] <- result_male$upper_ci
  }
}

# Diabetes Mellitus (Female)
b1_dm_fem <- -1.3133910
b2_dm_fem <- 1.0142390
vcov_diabetes_female <- matrix(c(0.1681525, -0.2240129,
                                 -0.2240129, 0.29964479),
                               nrow = 2, ncol = 2, byrow = TRUE)

rr_diabetes_female_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  
  # Compute the relative risk
  rr <- exp(sqrt(x / 100) * b1 + (x / 100) * b2)
  
  return(rr)
}

dm_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(dm_fem$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result <- confint_paf_vcov(gamma = gammas[[i]][[j]], 
                                    betas = c(b1_dm_fem, b2_dm_fem), 
                                    cov_matrix = vcov_diabetes_female, 
                                    p_abs = p_abs_list[[j]][[i]], 
                                    p_form = p_form_list[[j]][[i]], 
                                    rr_fd = 1.14, 
                                    rr_function = rr_diabetes_female_fun)
    
    # Guarda los resultados en el data frame
    dm_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
    dm_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    dm_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}

# TUBERCULOSIS (BOTH SEXES)
b_tb <- 0.0179695
var_tb <- 0.007215**2


tb_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

for (i in 1:length(tb_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b_tb, 
                  var_beta = var_tb, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      tb_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
      tb_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
      tb_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", tb_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}

tb_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(tb_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas_male[[i]][[j]], 
                  beta = b_tb, 
                  var_beta = var_tb, 
                  p_abs = p_abs_list_male[[j]][[i]], 
                  p_form = p_form_list_male[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      tb_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
      tb_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
      tb_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", tb_male$Year[i], "y grupo Male", j, "\n")
    }
  }
}

# HIV/AIDS (Male)
# ONLY FOR PPL OVER 61 GRS
b_hiv_male <- log(1.54)
# for consumption > 61 grams/day
var_hiv_male <- 0.078210772**2

rr_hiv_male <- function(x, beta) {
  b1 <- beta[1]
  
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
  # Case 1: x <= 61
  rr[x <= 61] <- 1
  
  # Case 2: x > 61
  rr[x > 61] <- exp(b1)
  
  return(rr)
}
# HIV/AIDS (Female)
# ONLY FOR PPL OVER 49
b_hiv_fem <- log(1.54)
# for consumption > 49 grams/day
var_hiv_fem <- 0.0782107722**2

rr_hiv_fem <- function(x, beta) {
  b1 <- beta[1]
  
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
  # Case 1: x <= 61
  rr[x <= 49] <- 1
  
  # Case 2: x > 61
  rr[x > 49] <- exp(b1)
  
  return(rr)
}

hiv_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

for (i in 1:length(hiv_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b_hiv_fem, 
                  var_beta = var_hiv_fem, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_hiv_fem)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      hiv_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
      hiv_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
      hiv_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", hiv_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}

hiv_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(hiv_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas_male[[i]][[j]], 
                  beta = b_hiv_male, 
                  var_beta = var_hiv_male, 
                  p_abs = p_abs_list_male[[j]][[i]], 
                  p_form = p_form_list_male[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_hiv_male)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      hiv_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
      hiv_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
      hiv_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", hiv_male$Year[i], "y grupo Male", j, "\n")
    }
  }
}

# LOWER RESPORATORY INFECTIONS (BOTH SEXES)
b_lri <- 0.4764038
var_lri <- 0.19220552^2
rr_lri_fun <- function(x, beta) {
  b1 <- beta[1]
  
  # Calculate y1
  y1 <- (x + 0.0399999618530273) / 100
  
  # Compute the relative risk
  rr <- exp(b1 * y1)
  
  return(rr)
}

lri_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

for (i in 1:length(lri_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b_lri, 
                  var_beta = var_lri, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_lri_fun)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      lri_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
      lri_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
      lri_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", lri_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}

lri_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(lri_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas_male[[i]][[j]], 
                  beta = b_lri, 
                  var_beta = var_lri, 
                  p_abs = p_abs_list_male[[j]][[i]], 
                  p_form = p_form_list_male[[j]][[i]], 
                  rr_form = 1, 
                  rr_function = rr_lri_fun)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      lri_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
      lri_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
      lri_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", lri_male$Year[i], "y grupo Male", j, "\n")
    }
  }
}

# LIVER CIRRHOSIS (MALES)
b1_lc_male <- 1.687111
b2_lc_male <- 1.106413
rr_fd_lc_male <- 3.26
vcov_lc_male <- matrix(c(0.0359478, -0.0359478,
                         -0.0359478, 0.07174495),
                       nrow = 2, ncol = 2, byrow = TRUE)

# LIVER CIRRHOSIS (FEMALES)
b1_lc_female <- 2.351821
b2_lc_female <- 0.9002139
vcov_lc_female <- matrix(c(0.05018842, -0.05018842,
                                        -0.05018842, 0.10270352),
                                      nrow = 2, ncol = 2, byrow = TRUE)


# For former drinkers (RRFD)
rr_fd_lc_female <- 3.26

rr_lc_male_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  
  # Compute the relative risk
  rr <- numeric(length(x))
  rr[x<=1] <- 1 + x * (exp((b1 + b2) * (1 + 0.1699981689453125) / 100))
  rr[x>1] <- exp((b1 + b2) * (x + 0.1699981689453125) / 100)
  
  return(rr)
}

rr_lc_female_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  
  # Compute the relative risk
  rr <- numeric(length(x))
  rr[x<=1] <- 1 + x * exp((b1 + b2) * sqrt((1 + 0.1699981689453125) / 100))
  rr[x>1] <- exp((b1 + b2) * sqrt((x + 0.1699981689453125) / 100))
  
  return(rr)
}

lc_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(lc_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result_male <- confint_paf_vcov(gamma = gammas_male[[i]][[j]], 
                                    betas = c(b1_lc_male, b2_lc_male), 
                                    cov_matrix = vcov_lc_male, 
                                    p_abs = p_abs_list_male[[j]][[i]], 
                                    p_form = p_form_list_male[[j]][[i]], 
                                    rr_fd = rr_fd_lc_male, 
                                    rr_function = rr_lc_male_fun)
    
    # Guarda los resultados en el data frame
    lc_male[i, paste0("Male", j, "_point")] <- result_male$point_estimate
    lc_male[i, paste0("Male", j, "_lower")] <- result_male$lower_ci
    lc_male[i, paste0("Male", j, "_upper")] <- result_male$upper_ci
  }
}

lc_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(lc_fem$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result <- confint_paf_vcov(gamma = gammas[[i]][[j]], 
                               betas = c(b1_lc_female, b2_lc_female), 
                               cov_matrix = vcov_lc_female, 
                               p_abs = p_abs_list[[j]][[i]], 
                               p_form = p_form_list[[j]][[i]], 
                               rr_fd = rr_fd_lc_female, 
                               rr_function = rr_lc_female_fun)
    
    # Guarda los resultados en el data frame
    lc_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
    lc_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    lc_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}

# PANCREATITIS 
b1_panc <- 0.0173451
var_panc <- 0.003803**2

# for former drinkers (both sexes)
rr_panc_fd <- 2.2

panc_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(panc_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas_male[[i]][[j]], 
                  beta = b1_panc, 
                  var_beta = var_panc, 
                  p_abs = p_abs_list_male[[j]][[i]], 
                  p_form = p_form_list_male[[j]][[i]], 
                  rr_form = rr_panc_fd, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      panc_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
      panc_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
      panc_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", lri_male$Year[i], "y grupo Male", j, "\n")
    }
  }
}


panc_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

for (i in 1:length(panc_fem$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b1_panc, 
                  var_beta = var_panc, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = rr_panc_fd, 
                  rr_function = rr_linear)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      panc_fem[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
      panc_fem[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
      panc_fem[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", lri_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}

############
# INJURIES #
############

# ROAD INJURIES
# Coefficients for Non-HED
b1_ri <- 0.00299550897979837
var_b1_ri <- 0.00050867822^2
# Coefficients for HED
b2_ri <- 0.959350221334602

# Variance-Covariance Matrix for HED
cov_matrix_ri <- matrix(c(0.00050867822^2, 0,
                          0, 0.227875857649849^2), 
                        nrow = 2, ncol = 2, byrow = TRUE)


ri_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)
# Loop over each year and each female age group
for (i in 1:length(ri_fem$Year)) {
  for (j in 1:4) {
    # Apply confint_paf_hed for each year and age group
    result <- confint_paf_hed(
      gammas = gammas_fem_hed[[i]][[j]], 
      betas = c(b1_ri,b2_ri), 
      cov_matrix = cov_matrix_ri, 
      p_abs = p_abs_list[[j]][[i]], 
      p_form = p_form_list[[j]][[i]], 
      rr_fd = 1, 
      rr_function_nhed = rr_linear,
      rr_function_hed = rr_ri_hed_fun,
      p_hed = p_hed_list_fem[[j]][[i]]
    )
    
    # Store results in the data frame
    ri_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
    ri_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    ri_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}

# Initialize a data frame to store results for males
ri_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# Loop over each year and each male age group
for (i in 1:nrow(ri_male)) {
  year <- as.character(ri_male$Year[i])  # Convert year to string for indexing
  
  for (j in 1:4) {
    # Get the gamma fits for nhed and hed for the current year and age group for males
    gamma_nhed <- gamma_fits_male[[year]][[paste0("tramo", j, "_nhed")]]
    gamma_hed <- gamma_fits_male[[year]][[paste0("tramo", j, "_hed")]]
    
    # Apply confint_paf_hed for each year and age group
    result <- confint_paf_hed(
      gammas = list(gamma_nhed, gamma_hed), 
      betas = c(b1_ri, b2_ri), 
      cov_matrix = cov_matrix_ri, 
      p_abs = p_abs_list_male[[j]][[i]], 
      p_form = p_form_list_male[[j]][[i]], 
      rr_fd = 1, 
      rr_function_nhed = rr_linear,
      rr_function_hed = rr_ri_hed_fun,
      p_hed = p_hed_list_male[[j]][[i]]
    )
    
    # Store results in the data frame
    ri_male[i, paste0("Male", j, "_point")] <- result$point_estimate
    ri_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
    ri_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
  }
}

# Display results
print(ri_male)

#POISONINGS, FALLS, FIRE, HEAT AND HOT SUBSTANCES, DROWING,
# EXPOSURE TO MECHANICAL FORCES, OTHER UNINTENTIONAL INJURIES

b1_injuries <- 0.0199800266267306
b2_injuries <- 0.647103242058538
vcov_injuries <- matrix(c(0.000509186^2,0,0,
                          0.155119431^2), nrow = 2, byrow = T)

# Initialize a data frame to store results for males
injuries_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# Loop over each year and each male age group
for (i in 1:nrow(injuries_male)) {
  year <- as.character(injuries_male$Year[i])  # Convert year to string for indexing
  
  for (j in 1:4) {
    # Get the gamma fits for nhed and hed for the current year and age group for males
    gamma_nhed <- gamma_fits_male[[year]][[paste0("tramo", j, "_nhed")]]
    gamma_hed <- gamma_fits_male[[year]][[paste0("tramo", j, "_hed")]]
    
    # Apply confint_paf_hed for each year and age group
    result <- confint_paf_hed(
      gammas = list(gamma_nhed, gamma_hed), 
      betas = c(b1_injuries, b2_injuries), 
      cov_matrix = vcov_injuries, 
      p_abs = p_abs_list_male[[j]][[i]], 
      p_form = p_form_list_male[[j]][[i]], 
      rr_fd = 1, 
      rr_function_nhed = rr_linear,
      rr_function_hed = rr_ri_hed_fun,
      p_hed = p_hed_list_male[[j]][[i]]
    )
    
    # Store results in the data frame
    injuries_male[i, paste0("Male", j, "_point")] <- result$point_estimate
    injuries_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
    injuries_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
  }
}

injuries_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)
# Loop over each year and each female age group
for (i in 1:length(injuries_fem$Year)) {
  for (j in 1:4) {
    # Apply confint_paf_hed for each year and age group
    result <- confint_paf_hed(
      gammas = gammas_fem_hed[[i]][[j]], 
      betas = c(b1_injuries,b2_injuries), 
      cov_matrix = vcov_injuries, 
      p_abs = p_abs_list[[j]][[i]], 
      p_form = p_form_list[[j]][[i]], 
      rr_fd = 1, 
      rr_function_nhed = rr_linear,
      rr_function_hed = rr_ri_hed_fun,
      p_hed = p_hed_list_fem[[j]][[i]]
    )
    
    # Store results in the data frame
    injuries_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
    injuries_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    injuries_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}

# self inflicted injuries and interpersonal violence

b1_violence <- 0.0199800266267306
b2_violence <- 0.647103242058538
vcov_violence <- matrix(c(0.000509186^2,0,0,
                          0.155119431^2), nrow = 2, byrow = T)

# Initialize a data frame to store results for males
violence_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# Loop over each year and each male age group
for (i in 1:nrow(violence_male)) {
  year <- as.character(violence_male$Year[i])  # Convert year to string for indexing
  
  for (j in 1:4) {
    # Get the gamma fits for nhed and hed for the current year and age group for males
    gamma_nhed <- gamma_fits_male[[year]][[paste0("tramo", j, "_nhed")]]
    gamma_hed <- gamma_fits_male[[year]][[paste0("tramo", j, "_hed")]]
    
    # Apply confint_paf_hed for each year and age group
    result <- confint_paf_hed(
      gammas = list(gamma_nhed, gamma_hed), 
      betas = c(b1_injuries, b2_injuries), 
      cov_matrix = vcov_injuries, 
      p_abs = p_abs_list_male[[j]][[i]], 
      p_form = p_form_list_male[[j]][[i]], 
      rr_fd = 1, 
      rr_function_nhed = rr_linear,
      rr_function_hed = rr_ri_hed_fun,
      p_hed = p_hed_list_male[[j]][[i]]
    )
    
    # Store results in the data frame
    violence_male[i, paste0("Male", j, "_point")] <- result$point_estimate
    violence_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
    violence_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
  }
}

violence_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)
# Loop over each year and each female age group
for (i in 1:length(violence_fem$Year)) {
  for (j in 1:4) {
    # Apply confint_paf_hed for each year and age group
    result <- confint_paf_hed(
      gammas = gammas_fem_hed[[i]][[j]], 
      betas = c(b1_injuries,b2_injuries), 
      cov_matrix = vcov_injuries, 
      p_abs = p_abs_list[[j]][[i]], 
      p_form = p_form_list[[j]][[i]], 
      rr_fd = 1, 
      rr_function_nhed = rr_linear,
      rr_function_hed = rr_ri_hed_fun,
      p_hed = p_hed_list_fem[[j]][[i]]
    )
    
    # Store results in the data frame
    violence_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
    violence_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    violence_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}



##################
# CARDIOVASCULAR #
##################

# HYPERTENSIVE HEART DISEASE (MALE)
b1_hhd_male <- 0.0150537
var_hhd_male <- 0.0024196**2
rr_fd_hhd_male <- 1.05

rr_hhd_male_fun <- function(x, beta) {
  rr = exp(beta * x)
  return(rr)
}


# Define un data frame para almacenar los resultados para hombres
hhd_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(hhd_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result <- tryCatch({
      confint_paf(gamma = gammas_male[[i]][[j]], 
                  beta = b1_hhd_male, 
                  var_beta = var_hhd_male, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = rr_fd_hhd_male, 
                  rr_function = rr_hhd_male_fun)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result) && !is.null(result$Point_Estimate)) {
      hhd_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
      hhd_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
      hhd_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", hhd_male$Year[i], "y grupo Male", j, "\n")
    }
  }
}


# HYPERTENSIVE HEART DISEASE (FEMALE)
b1_hhd_fem <- 0.0154196
var_hhd_fem <- 0.006459**2

rr_hhd_female_fun <- function(x, beta) {
  rr = exp(beta * x)
  return(rr)
}

# Define a data frame to store results
hhd_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Bucle sobre cada año y cada grupo femenino
for (i in 1:length(hhd_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result <- confint_paf_vcov(gamma = gammas[[i]][[j]], 
                               betas = b1_hhd_fem, 
                               cov_matrix = sqrt(var_hhd_fem), 
                               p_abs = p_abs_list[[j]][[i]], 
                               p_form = p_form_list[[j]][[i]], 
                               rr_fd = 1, 
                               rr_function = rr_hhd_female_fun)
    
    # Guarda los resultados en el data frame
    hhd_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
    hhd_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    hhd_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}

# ISCHAEMIC HEART DISEASE MALE
# ISCHAEMIC HEART DISEASE FEMALE

# INTRACEREBRAL HAEMORRHAGE (MALES)
b1_ich_male <- 0.6898937
var_b1_ich_male <- 0.1141980
# INTRACEREBRAL HAEMORRHAGE (FEMALES)
b1_ich_female <- 1.466406
var_b1_ich_female <- 0.3544172

rr_fd_ich_female <- 1.36

rr_ich_male <- function(x, beta) {
  b1 <- beta
  
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
  # Case 1: x <= 1
  rr[x <= 1] <- 1 - x[x <= 1] * (1 - exp(b1 * (1 + 0.0028572082519531) / 100))
  
  # Case 2: x > 1
  rr[x > 1] <- exp(b1 * (x[x > 1] + 0.0028572082519531) / 100)
  
  return(rr)
}
ich_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(ich_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_ich_male <- tryCatch({
      confint_paf(gamma = gammas_male[[i]][[j]], 
                  beta = b1_ich_male, 
                  var_beta = var_b1_ich_male, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = 1.25, 
                  rr_function = rr_ich_male)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_ich_male) && !is.null(result_ich_male$Point_Estimate)) {
      ich_male[i, paste0("Male", j, "_point")] <- result_ich_male$Point_Estimate
      ich_male[i, paste0("Male", j, "_lower")] <- result_ich_male$Lower_CI
      ich_male[i, paste0("Male", j, "_upper")] <- result_ich_male$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", ich_male$Year[i], "y grupo Male", j, "\n")
    }
  }
}

ich_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

for (i in 1:length(ich_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf para cada grupo
    result_ich <- tryCatch({
      confint_paf(gamma = gammas[[i]][[j]], 
                  beta = b1_ich_female, 
                  var_beta = var_b1_ich_female, 
                  p_abs = p_abs_list[[j]][[i]], 
                  p_form = p_form_list[[j]][[i]], 
                  rr_form = rr_fd_ich_female, 
                  rr_function = rr_ich_male)
    }, error = function(e) NULL)
    
    # Comprobación de que result_bcan tenga contenido antes de asignar
    if (!is.null(result_ich) && !is.null(result_ich$Point_Estimate)) {
      ich_female[i, paste0("Fem", j, "_point")] <- result_ich$Point_Estimate
      ich_female[i, paste0("Fem", j, "_lower")] <- result_ich$Lower_CI
      ich_female[i, paste0("Fem", j, "_upper")] <- result_ich$Upper_CI
    } else {
      # Opcional: mensaje para ver dónde ocurre el problema
      cat("Error al calcular para el año:", ich_female$Year[i], "y grupo Fem", j, "\n")
    }
  }
}

# ISCHAEMIC STROKE MALE
b1_stroke_male1 <- 1.111874
b2_stroke_male1 <- 0.4030081
b3_stroke_male1 <- 0.3877538

cov_matrix_stroke_male1 <- matrix(c(0.199819, 0, 0,
                                    0, 0.003720, 0.002018,
                                    0, 0.002018, 0.003539), 
                                  nrow = 3, ncol = 3, byrow = TRUE)

rr_stroke_male_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  b3 <- betas[3]
  
  # Calculate y1 and y2
  y1 <- (1 + 0.0028572082519531) / 100
  y2 <- (x + 0.0028572082519531) / 100
  
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
  # Case 1: x ≤ 1
  rr[x <= 1] <- 1 - x[x <= 1] * (1 - exp(b1 * (b2 * sqrt(y1) + b3 * sqrt(y1) * log(y1))))
  
  # Case 2: x > 1
  rr[x > 1] <- exp(b1 * (b2 * sqrt(y2[x > 1]) + b3 * sqrt(y2[x > 1]) * log(y2[x > 1])))
  
  return(rr)
}

stroke_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# Bucle sobre cada año y cada grupo masculino
for (i in 1:length(stroke_male$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result_male <- confint_paf_vcov(gamma = gammas_male[[i]][[j]], 
                                    betas = c(b1_stroke_male1, b2_stroke_male1,
                                              b3_stroke_male1), 
                                    cov_matrix = cov_matrix_stroke_male1, 
                                    p_abs = p_abs_list_male[[j]][[i]], 
                                    p_form = p_form_list_male[[j]][[i]], 
                                    rr_fd = 1, 
                                    rr_function = rr_stroke_male_fun)
    
    # Guarda los resultados en el data frame
    stroke_male[i, paste0("Male", j, "_point")] <- result_male$point_estimate
    stroke_male[i, paste0("Male", j, "_lower")] <- result_male$lower_ci
    stroke_male[i, paste0("Male", j, "_upper")] <- result_male$upper_ci
  }
}

# ISCHAEMIC STROKE FEMALE
b1_stroke_female1 <- 1.111874
b2_stroke_female1 <- -2.48768
b3_stroke_female1 <- 3.7087240

cov_matrix_stroke_female1 <- matrix(c(
  0.446999664429404,  0,                 0,              
  0,                 0.4875627,      -0.31633911,     
  0,                 -0.31633911,       0.6645197    
), nrow = 3, ncol = 3, byrow = TRUE)

rr_stroke_fem_fun1 <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  b3 <- betas[3]
  
  # Calculate y1 and y2
  y1 <- (1 + 0.0028572082519531) / 100
  y2 <- (x + 0.0028572082519531) / 100
  
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
  # Case 1: x ≤ 1
  rr[x <= 1] <- 1 - x[x <= 1] * (1 - exp(b1 * (b2 * sqrt(y1) + b3 * sqrt(y1))))
  
  # Case 2: x > 1
  rr[x > 1] <- exp(b1 * (b2 * sqrt(y2[x > 1]) + b3 * sqrt(y2[x > 1])))
  
  return(rr)
}

# Define a data frame to store results
stroke_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Bucle sobre cada año y cada grupo femenino
for (i in 1:length(stroke_female$Year)) {
  for (j in 1:4) {
    # Llama a confint_paf_vcov para cada grupo
    result <- confint_paf_hed(gamma = gammas[[i]][[j]], 
                              betas = c(b1_stroke_female1,
                                        b2_stroke_female1,
                                        b3_stroke_female1), 
                              cov_matrix = cov_matrix_stroke_female1, 
                              p_abs = p_abs_list[[j]][[i]], 
                              p_form = p_form_list[[j]][[i]], 
                              rr_fd = 1, 
                              rr_function = rr_stroke_fem_fun1)
    
    # Guarda los resultados en el data frame
    stroke_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
    stroke_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    stroke_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}
x_vals_nhed <- seq(0.1, 60, length.out = 1500)
x_vals_hed <- seq(60,150, length.out = 1500)

stroke_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# Stroke-specific parameters for females
betas_stroke_female <- c(b1_stroke_female1, b2_stroke_female1, b3_stroke_female1)
cov_matrix_stroke <- cov_matrix_stroke_female1

# Loop over each year and each female age group
for (i in 1:length(stroke_female$Year)) {
  for (j in 1:4) {
    # Apply confint_paf_hed for each year and age group
    result <- confint_paf_hed1(
      gammas = gammas_fem_hed[[i]][[j]], 
      betas = betas_stroke_female, 
      cov_matrix = cov_matrix_stroke, 
      p_abs = p_abs_list[[j]][[i]], 
      p_form = p_form_list[[j]][[i]], 
      rr_fd = 1, 
      rr_function_hed = rr_stroke_fem_fun1,
      p_hed = p_hed_list_fem[[j]][[i]]
    )
    
    # Store results in the data frame
    stroke_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
    stroke_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
    stroke_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
  }
}

confint_paf_hed1 <- function(gammas, betas, cov_matrix, p_abs, p_form, rr_fd,
                             rr_function_hed, p_hed) {
  set.seed(145)
  n_sim <- 10000
  simulated_pafs <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Simulate PCA mean and SD using the gamma distribution
    pca_sim_nhed <- rgamma(1000, shape = gammas[[1]]$estimate["shape"], 
                           rate = gammas[[1]]$estimate["rate"])
    pca_sim_hed <- rgamma(1000, shape = gammas[[2]]$estimate["shape"], 
                          rate = gammas[[2]]$estimate["rate"])
    
    # Calculate the shape and rate parameters from the simulated PCA values
    mean_sim_nhed <- mean(pca_sim_nhed)
    sd_sim_nhed <- sd(pca_sim_nhed)
    mean_sim_hed <- mean(pca_sim_hed)
    sd_sim_hed <- sd(pca_sim_hed)
    
    # Calculate shape and rate for the new gamma distribution
    shape_sim_nhed <- (mean_sim_nhed / sd_sim_nhed)^2
    rate_sim_nhed <- mean_sim_nhed / (sd_sim_nhed^2)
    shape_sim_hed <- (mean_sim_hed / sd_sim_hed)^2
    rate_sim_hed <- mean_sim_hed / (sd_sim_hed^2)
    
    # Simulate the gamma distribution
    y_gamma_sim_nhed <- dgamma(x_vals_nhed, shape = shape_sim_nhed, rate = rate_sim_nhed)
    y_gamma_sim_hed_60 <- dgamma(x_vals_nhed, shape = shape_sim_hed, rate = rate_sim_hed)
    y_gamma_sim_hed_150 <- dgamma(x_vals_hed, shape = shape_sim_hed, rate = rate_sim_hed)
    
    # Simulate the beta coefficients jointly from a multivariate normal distribution
    beta_sim <- MASS::mvrnorm(1, mu = betas, Sigma = cov_matrix)
    
    # Compute the relative risk based on the simulated betas
    rr_sim_nhed <- 1
    rr_sim_hed_60 <- rr_function_hed(x = x_vals_nhed, betas = beta_sim)
    rr_sim_hed_150 <- rr_function_hed(x = x_vals_hed, betas = beta_sim)
    
    # Simulate proportions of lifetime abstainers and former drinkers
    prop_abs_sim <- max(rnorm(1000, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000)), 0.001)
    prop_form_sim <- max(rnorm(1000, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000)), 0.001)
    prop_hed_sim <- max(rnorm(1000, mean = p_hed, sd = sqrt(p_hed * (1 - p_hed) / 1000)), 0.001)
    
    # Calculate PAF using the customized function
    simulated_pafs[i] <- paf_hed_function(
      x_60 = x_vals_nhed, 
      x_150 = x_vals_hed, 
      y_nhed = y_gamma_sim_nhed, 
      y_hed_60 = y_gamma_sim_hed_60, 
      y_hed_150 = y_gamma_sim_hed_150, 
      rr_nhed = rr_sim_nhed, 
      rr_hed_60 = rr_sim_hed_60, 
      rr_hed_150 = rr_sim_hed_150,
      rr_form = rr_fd, 
      p_abs = prop_abs_sim, 
      p_form = prop_form_sim, 
      p_hed = prop_hed_sim
    )
  }
  
  # Remove NaN values from simulated PAFs
  simulated_pafs <- simulated_pafs[!is.nan(simulated_pafs)]
  
  if (length(simulated_pafs) == 0) {
    stop("All simulations resulted in NaN values. Please check your input parameters.")
  }
  
  # Calculate the 95% confidence interval
  paf_lower <- quantile(simulated_pafs, 0.025)
  paf_upper <- quantile(simulated_pafs, 0.975)
  paf_point_estimate <- mean(simulated_pafs)
  
  return(list(point_estimate = round(paf_point_estimate,3), lower_ci = round(paf_lower,3), upper_ci = round(paf_upper,3)))
}
