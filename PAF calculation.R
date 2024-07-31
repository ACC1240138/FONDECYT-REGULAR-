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



####### CATEGORICAL VERSION
table(data$cvolaj)
capped_data <- data %>% 
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
                            sexo == "Hombre" & volajohdia >= 60 ~ "cat3",
                            TRUE ~ NA))
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

# STRATIFIED CONSUMPTION LEVEL
oh_level_median <- capped_data %>% 
  filter(!is.na(volajohdia),!is.na(cvolaj)) %>% 
  group_by(sexo, edad_tramo,cvolaj) %>%
  summarise(oh_level = round(median(volajohdia),1))

oh_level_mid <- oh_level_median %>% 
  mutate(oh_level = case_when(sexo == "Hombre" & cvolaj == "cat1" ~ 10,
                              sexo == "Hombre" & cvolaj == "cat2" ~ 30,
                              sexo == "Hombre" & cvolaj == "cat3" ~ 105,
                              sexo == "Mujer" & cvolaj == "cat1" ~ 10,
                              sexo == "Mujer" & cvolaj == "cat2" ~ 30,
                              sexo == "Mujer" & cvolaj == "cat3" ~ 95,
                              TRUE ~ 0))
# STRATIFIED PROPORTION OF CONSUMPTION LEVEL
prop_level <- capped_data %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo, edad_tramo, cvolaj) %>%
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>%
  group_by(sexo, edad_tramo) %>%
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>%
  dplyr::select(-weighted_count)

# JOIN BOTH TABLES
paf_base_median <- oh_level_median %>% 
  left_join(prop_level, by = c("sexo", "edad_tramo","cvolaj"))
paf_base_mid <- oh_level_mid %>% 
  left_join(prop_level, by = c("sexo", "edad_tramo","cvolaj"))

####################################
# ESTIMATING AAF FOR BREAST CANCER #
####################################

paf_bc_median <- paf_base_median %>% 
  filter(sexo=="Mujer") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b1_bcan*oh_level),1),NA))

paf_bc_mid <- paf_base_mid %>% 
  filter(sexo=="Mujer") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b1_bcan*oh_level),1),NA))

cat_paf_calculator <- function(data){
  num1 = data[1,5]*(data[1,6]-1)+data[2,5]*(data[2,6]-1)+data[3,5]*(data[3,6]-1)
  den1 = num1+1
  num2 = data[6,5]*(data[6,6]-1)+data[7,5]*(data[7,6]-1)+data[8,5]*(data[8,6]-1)
  den2 = num2+1
  num3 = data[11,5]*(data[11,6]-1)+data[12,5]*(data[12,6]-1)+data[13,5]*(data[13,6]-1)
  den3 = num3+1
  num4 = data[16,5]*(data[16,6]-1)+data[17,5]*(data[17,6]-1)+data[18,5]*(data[18,6]-1)
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
cat_paf_calculator(paf_bc_mid)


####################################
# ESTIMATING AAF FOR TUBERCULOSIS  #
####################################
paf_tb_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs", cvolaj != "fd") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b_tb*oh_level),1),NA))
paf_tb_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs", cvolaj != "fd") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b_tb*oh_level),1),NA))

cat_paf <- function(data){
  num1 = data[1,5]*(data[1,6]-1)+data[2,5]*(data[2,6]-1)+data[3,5]*(data[3,6]-1)
  den1 = num1+1
  num2 = data[4,5]*(data[4,6]-1)+data[5,5]*(data[5,6]-1)+data[6,5]*(data[6,6]-1)
  den2 = num2+1
  num3 = data[7,5]*(data[7,6]-1)+data[8,5]*(data[8,6]-1)+data[9,5]*(data[9,6]-1)
  den3 = num3+1
  num4 = data[10,5]*(data[10,6]-1)+data[11,5]*(data[11,6]-1)+data[12,5]*(data[12,6]-1)
  den4 = num4+1
  
  num5 = data[13,5]*(data[13,6]-1)+data[14,5]*(data[14,6]-1)+data[15,5]*(data[15,6]-1)
  den5 = num5+1
  num6 = data[16,5]*(data[16,6]-1)+data[17,5]*(data[17,6]-1)+data[18,5]*(data[18,6]-1)
  den6 = num6+1
  num7 = data[19,5]*(data[19,6]-1)+data[20,5]*(data[20,6]-1)+data[21,5]*(data[21,6]-1)
  den7 = num7+1
  num8 = data[22,5]*(data[22,6]-1)+data[23,5]*(data[23,6]-1)+data[24,5]*(data[24,6]-1)
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
cat_paf(paf_tb_mid)

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
paf_lri_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs", cvolaj != "fd") %>% 
  mutate(rr = ifelse(cvolaj != "ltabs" & cvolaj != "fd",round(exp(b_lri*((oh_level+0.0399999618530273)/100)),1),NA))

cat_paf(paf_lri_median)
cat_paf(paf_lri_mid)

####################################################
# ESTIMATING AAF FOR LIP AND ORAL CAVITY CANCER   #
####################################################
paf_locan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_locan*oh_level+b2_locan*oh_level**2),1)))
paf_locan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_locan*oh_level+b2_locan*oh_level**2),1)))

cat_paf_fd <- function(data){
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
cat_paf_fd(paf_locan_mid)
cat_paf_fd(paf_locan_median)
cat_paf(paf_locan_median %>% filter(cvolaj != "fd"))
cat_paf(paf_locan_mid %>% filter(cvolaj != "fd"))

################################################
# ESTIMATING AAF FOR OTHER PHARINGEAL CANCER   #
################################################
paf_opcan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_locan*oh_level+b2_locan*oh_level**2),1)))
paf_opcan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.2,round(exp(b1_opcan*oh_level+b2_opcan*oh_level**2),1)))

cat_paf_fd(paf_opcan_median)
cat_paf_fd(paf_opcan_mid)
cat_paf(paf_opcan_median %>% filter(cvolaj != "fd"))
cat_paf(paf_opcan_mid %>%  filter(cvolaj != "fd"))

##########################################
# ESTIMATING AAF FOR OESOPHAFUS CANCER   #
##########################################
paf_oescan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.16,round(exp(b1_oescan*oh_level+b2_oescan*oh_level**3),1)))
paf_oescan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.16,round(exp(b1_oescan*oh_level+b2_oescan*oh_level**3),1)))

cat_paf_fd(paf_oescan_median)
cat_paf_fd(paf_oescan_mid)

cat_paf(paf_oescan_median %>% filter(cvolaj != "fd"))
cat_paf(paf_oescan_mid %>% filter(cvolaj != "fd"))
################################################
# ESTIMATING AAF FOR COLON AND RECTUM CANCER   #
################################################
paf_crcan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.19,
                        cvolaj == "fd" & sexo == "Mujer"~1.05,
                        sexo == "Hombre"~exp(b1_crcan_male*oh_level),
                        sexo == "Mujer"~exp(b1_crcan_fem*oh_level)))
paf_crcan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.19,
                        cvolaj == "fd" & sexo == "Mujer"~1.05,
                        sexo == "Hombre"~exp(b1_crcan_male*oh_level),
                        sexo == "Mujer"~exp(b1_crcan_fem*oh_level)))

cat_paf_fd(paf_crcan_median)
cat_paf_fd(paf_crcan_mid)
cat_paf(paf_crcan_median %>% filter(cvolaj != "fd"))
cat_paf(paf_crcan_mid %>% filter(cvolaj != "fd"))

######################################
# ESTIMATING AAF FOR LIVER CANCER   #
#####################################
paf_lican_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.23,
                        cvolaj == "fd" & sexo == "Mujer"~2.68,
                        sexo == "Hombre"~exp(b1_lican_male*oh_level),
                        sexo == "Mujer"~exp(b1_lican_fem*oh_level)))
paf_lican_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~2.23,
                        cvolaj == "fd" & sexo == "Mujer"~2.68,
                        sexo == "Hombre"~exp(b1_lican_male*oh_level),
                        sexo == "Mujer"~exp(b1_lican_fem*oh_level)))

cat_paf_fd(paf_lican_median)
cat_paf_fd(paf_lican_mid)
cat_paf(paf_lican_median %>% filter(cvolaj != "fd"))
cat_paf(paf_lican_mid %>% filter(cvolaj != "fd"))


######################################
# ESTIMATING AAF FOR LARYNX CANCER   #
#####################################
paf_lxcan_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.18,exp(b1_lxcan*oh_level)))
                       
paf_lxcan_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1.18,exp(b1_lxcan*oh_level)))


cat_paf_fd(paf_lxcan_median)
cat_paf_fd(paf_lxcan_mid)
cat_paf(paf_lxcan_median %>% filter(cvolaj != "fd"))
cat_paf(paf_lxcan_mid %>% filter(cvolaj != "fd"))

########################################
# ESTIMATING AAF FOR DIABETES MELLITUS #
########################################
paf_dm_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~1.18,
                        cvolaj == "fd" & sexo == "Mujer"~1.14,
                        sexo == "Hombre"~exp((b1_dm_male*(oh_level/100)**2)+(b2_dm_male*(oh_level/100)**3)),
                        sexo == "Mujer"~exp(b1_dm_fem*sqrt(oh_level/100)+b2_dm_fem*(oh_level/100))))
paf_dm_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = case_when(cvolaj == "fd" & sexo == "Hombre"~1.18,
                        cvolaj == "fd" & sexo == "Mujer"~1.14,
                        sexo == "Hombre"~exp((b1_dm_male*(oh_level/100)**2)+(b2_dm_male*(oh_level/100)**3)),
                        sexo == "Mujer"~exp(b1_dm_fem*sqrt(oh_level/100)+b2_dm_fem*(oh_level/100))))


cat_paf_fd(paf_dm_median)
cat_paf_fd(paf_dm_mid)
cat_paf(paf_dm_median %>% filter(cvolaj != "fd"))
cat_paf(paf_dm_mid %>% filter(cvolaj != "fd"))

################################
# ESTIMATING AAF FOR EPILEPSY #
###############################
paf_epi_median <- paf_base_median %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1,exp(b1_epi*((oh_level+0.5)/100))))

paf_epi_mid <- paf_base_mid %>% 
  filter(cvolaj != "ltabs") %>% 
  mutate(rr = ifelse(cvolaj == "fd",1,exp(b1_epi*((oh_level+0.5)/100))))

cat_paf_fd(paf_epi_median)
cat_paf_fd(paf_epi_mid)
cat_paf(paf_epi_median %>% filter(cvolaj != "fd"))
cat_paf(paf_epi_mid %>% filter(cvolaj != "fd"))

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
cat_paf(paf_hhd)








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
