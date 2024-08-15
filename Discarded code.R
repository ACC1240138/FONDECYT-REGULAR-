oh_level_mid <- oh_level_median %>% 
  mutate(oh_level = case_when(sexo == "Hombre" & cvolaj == "cat1" ~ 10,
                              sexo == "Hombre" & cvolaj == "cat2" ~ 30,
                              sexo == "Hombre" & cvolaj == "cat3" ~ 105,
                              sexo == "Mujer" & cvolaj == "cat1" ~ 10,
                              sexo == "Mujer" & cvolaj == "cat2" ~ 30,
                              sexo == "Mujer" & cvolaj == "cat3" ~ 95,
                              TRUE ~ 0))

filtered_data <- data %>% 
  filter(volajohdia <= 150) %>% 
  mutate(edad_tramo = case_when(between(edad, 15, 29)~1,
                                between(edad, 30,44)~2,
                                between(edad,45,59)~3,
                                between(edad,60,65)~4),
         cvolaj = case_when(oh1 == "No" ~ "ltabs",
                            oh2 == ">30" | oh2 == ">1 año" ~ "fd",
                            sexo == "Mujer" & volajohdia > 0 & volajohdia <= 19.99 ~ "cat1",
                            sexo == "Mujer" & volajohdia >= 20 & volajohdia <= 39.99 ~ "cat2",
                            sexo == "Mujer" & volajohdia > 40  ~ "cat3",
                            sexo == "Hombre" & volajohdia > 0 & volajohdia <= 39.99 ~ "cat1",
                            sexo == "Hombre" & volajohdia >= 40 & volajohdia <= 59.99 ~ "cat2",
                            sexo == "Hombre" & volajohdia >= 60  ~ "cat3",
                            TRUE ~ NA))
paf_base_mid <- oh_level_mid %>% 
  left_join(prop_level, by = c("sexo", "edad_tramo","cvolaj"))
paf_bc_mid <- paf_base_mid %>% 
  filter(sexo=="Mujer") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b1_bcan*oh_level),1),NA))
paf_tb_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs", cvolaj != "fd") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b_tb*oh_level),1),NA))
paf_lri_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs", cvolaj != "fd") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b_lri*((oh_level+0.0399999618530273)/100)),1),NA))
paf_locan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_locan*oh_level+b2_locan*oh_level**2),1)))
paf_opcan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_opcan*oh_level+b2_opcan*oh_level**2),1)))
paf_oescan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.16,round(exp(b1_oescan*oh_level+b2_oescan*oh_level**3),1)))
paf_crcan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.19,
                        cvolaj == "fd" & sexo == "Mujer"~1.05,
                        sexo == "Hombre"~exp(b1_crcan_male*oh_level),
                        sexo == "Mujer"~exp(b1_crcan_fem*oh_level)))
paf_lican_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.23,
                        cvolaj == "fd" & sexo == "Mujer"~2.68,
                        sexo == "Hombre"~exp(b1_lican_male*oh_level),
                        sexo == "Mujer"~exp(b1_lican_fem*oh_level)))
paf_lxcan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.18,exp(b1_lxcan*oh_level)))
paf_dm_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~1.18,
                        cvolaj == "fd" & sexo == "Mujer"~1.14,
                        sexo == "Hombre"~exp((b1_dm_male*(oh_level/100)**2)+(b2_dm_male*(oh_level/100)**3)),
                        sexo == "Mujer"~exp(b1_dm_fem*sqrt(oh_level/100)+b2_dm_fem*(oh_level/100))))
paf_epi_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1,exp(b1_epi*((oh_level+0.5)/100))))

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


# EMPIRICAL DISTRIBUTION
# Define the relative risk coefficient
beta_1 <- 0.0179695
data08 <- data08 %>%
  filter(!is.na(volajohdia) & is.finite(volajohdia))
# Define the relative risk function
rr_function <- function(x) {
  exp(b1_panc * x)
}

# Define the prevalence function based on the empirical data
prevalence_function <- function(x, data) {
  # Compute the density estimate of the empirical data
  dens <- density(data$volajohdia, from = 0, to = 150, n = 512)
  # Interpolate the density estimate to get the value at x
  approx(dens$x, dens$y, xout = x)$y
}

# Define the integrand function that combines both the prevalence and relative risk functions
integrand <- function(x, data) {
  prevalence_function(x, data) * (rr_function(x) - 1)
}

# Perform numerical integration from 1 to 150
result <- integrate(integrand, lower = 1, upper = 150, data = data08)





data %>% 
  group_by(year) %>% 
  summarise(mean(volajohdia, na.rm = T),
            sd(volajohdia, na.rm = T))


data08 <- data %>% 
  filter(year == 2008)
table(!is.na(data08$voltotdia))

data08 %>% 
  summarise(sum(voltotdia*exp, na.rm = T)/sum(exp))

gamma_data <- data08 %>%
  filter(volajohdia > 0) %>% 
  mutate(edad_tramo = case_when(between(edad, 15, 29)~1,
                                between(edad, 30,44)~2,
                                between(edad,45,59)~3,
                                between(edad,60,65)~4)) %>% 
  dplyr::select(edad_tramo, volajohdia, sexo)
fit_gamma_model <- function(subset_data) {
  if (nrow(subset_data) < 2) return(NULL)  # Avoid fitting if too few data points
  
  # Remove NAs from the subset data
  subset_data <- na.omit(subset_data)
  
  # Fit the gamma regression model
  gamma_model <- gamlss(volajohdia ~ 1, family = GA, data = subset_data)
  
  # Predict the proportion of consumption for each observed level of volajohdia
  predicted_proportions <- predict(gamma_model, type = "response")
  
  # Create a data frame with the original rows and predicted proportions
  result <- data.frame(volajohdia = subset_data$volajohdia, predicted_proportions = predicted_proportions)
  
  return(result)
}


# Apply the function to each stratum of sex and age group
predicted_proportions_list <- gamma_data %>%
  group_by(sexo, edad_tramo) %>%
  mutate(predicted_proportions = ifelse(volajohdia > 0, {
    # Filter out zero values for fitting the gamma model
    fit_data <- cur_data() %>% filter(volajohdia > 0)
    
    if (nrow(fit_data) < 2) {
      rep(NA, n())
    } else {
      # Fit the gamma regression model
      gamma_model <- gamlss(volajohdia ~ 1, family = GA, data = fit_data)
      
      # Predict the proportion of consumption for each observed level of volajohdia
      predict(gamma_model, type = "response")
    }
  }, NA_real_))

gamma_model <- gamlss(volajohdia ~ 1, family = GA, data = gamma_data[gamma_data$edad_tramo == 1,])
fitted_values <- fitted(gamma_model, what = "mu")
sigma_values <- fitted(gamma_model, what = "sigma")
gamma_data <- gamma_data %>%
  mutate(predicted_proportions = dGA(volajohdia, mu = fitted_values, sigma = sigma_values))

hist(gamma_data$volajohdia)
sum(gamma_data$predicted_proportions)

predict_tramo1 <- predict(gamma_model, type = "response")
gamma_data <- gamma_data %>%
  mutate(predicted_values = predict(gamma_model, type = "response"))
hist(gamma_data$predicted_values)
# Calculate the total of the predicted values for normalizing
total_predicted <- sum(gamma_data$predicted_values, na.rm = TRUE)

# Calculate the predicted proportions as the predicted values divided by the total predicted values
gamma_data <- gamma_data %>%
  mutate(predicted_proportions = predicted_values / total_predicted)







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



























# PROPORTION OF LIFETIME ABSTAINERS
prop_abs_male <- capped_data %>% 
  filter(sexo == "Hombre") %>% 
  count(oh1) %>% 
  mutate(perc = round(n/sum(n),3)) %>% 
  filter(oh1 == "No") %>% 
  pull(perc)

# 0.173
prop_abs_fem <- capped_data %>% 
  filter(sexo == "Mujer") %>% 
  count(oh1) %>% 
  mutate(perc = round(n/sum(n),3)) %>% 
  filter(oh1 == "No") %>% 
  pull(perc)
# 0.271

# PROPORTION OF FORMER DRINKERS
prop_fd_male <- capped_data %>% 
  filter(sexo == "Hombre", !is.na(oh2)) %>% 
  count(oh2) %>% 
  mutate(prop = round(n/sum(n),3)) %>% 
  filter(oh2 == ">30" | oh2 == ">1 año") %>%
  mutate(prop_sum = sum(prop)) %>% 
  pull(prop_sum) %>% 
  unique()
# 0.377
prop_fd_fem <- capped_data %>% 
  filter(sexo == "Mujer", !is.na(oh2)) %>% 
  count(oh2) %>% 
  mutate(prop = round(n/sum(n),3)) %>% 
  filter(oh2 == ">30" | oh2 == ">1 año") %>%
  mutate(prop_sum = sum(prop)) %>% 
  pull(prop_sum) %>% 
  unique()
# 0.535

# PROPORTION OF CURRENT DRINKERS MALE
prop_cd_male <- capped_data %>% 
  filter(!is.na(cvolaj), cvolaj != "Abstinentes") %>% 
  group_by(sexo,edad_tramo) %>% 
  count(cvolaj) %>% 
  mutate(prop = round(n/sum(n),3)) %>% 
  filter(sexo == "Hombre") %>% 
  left_join(oh_level, by = c("sexo","edad_tramo","cvolaj"))

prop_cd_fem <- capped_data %>% 
  filter(!is.na(cvolaj), cvolaj != "Abstinentes") %>% 
  group_by(sexo,edad_tramo) %>% 
  count(cvolaj) %>% 
  mutate(prop = round(n/sum(n),3)) %>% 
  filter(sexo == "Mujer") %>% 
  left_join(oh_level, by = c("sexo","edad_tramo","cvolaj"))
# tuberculosis hombres

tb_male <- prop_cd_male %>% 
  mutate(rr = exp(0.0179695*prop)) %>% 
  dplyr::select(-n)
tb_fem <- prop_cd_fem %>% 
  mutate(rr = exp(0.0179695*prop)) %>% 
  dplyr::select(-n)

numerador <- 0.718*0.01+0.153*0+0.098*0+0.031*0
denominador <- 0.718*1.01+0.153*1.00+0.098*1+0.031*1
(prop_abs_fem+numerador)/(prop_abs_fem+denominador)
