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
  dplyr::select(year, sexo, exp, edad_tramo, volajohdia, cvolaj)


input <- data_input %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo,year, edad_tramo, cvolaj) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>% 
  filter(cvolaj == "fd" | cvolaj == "ltabs")

p_abs_male <- input %>% 
  filter(sexo == "Hombre", cvolaj == "ltabs")
p_form_male <- input %>% 
  filter(sexo == "Hombre", cvolaj == "fd")
p_abs_fem <- input %>% 
  filter(sexo == "Mujer", cvolaj == "ltabs")
p_form_fem <- input %>% 
  filter(sexo == "Mujer", cvolaj == "fd")

# PROPORTION OF LIFE TIME ABSTAINERS IN MALES
p_abs_male1_08 <- p_abs_male[1,6] %>% pull()
p_abs_male2_08 <- p_abs_male[2,6] %>% pull()
p_abs_male3_08 <- p_abs_male[3,6] %>% pull()
p_abs_male4_08 <- p_abs_male[4,6] %>% pull()

p_abs_male1_10 <- p_abs_male[5,6] %>% pull()
p_abs_male2_10 <- p_abs_male[6,6] %>% pull()
p_abs_male3_10 <- p_abs_male[7,6] %>% pull()
p_abs_male4_10 <- p_abs_male[8,6] %>% pull()

p_abs_male1_12 <- p_abs_male[9,6] %>% pull()
p_abs_male2_12 <- p_abs_male[10,6] %>% pull()
p_abs_male3_12 <- p_abs_male[11,6] %>% pull()
p_abs_male4_12 <- p_abs_male[12,6] %>% pull()

p_abs_male1_14 <- p_abs_male[13,6] %>% pull()
p_abs_male2_14 <- p_abs_male[14,6] %>% pull()
p_abs_male3_14 <- p_abs_male[15,6] %>% pull()
p_abs_male4_14 <- p_abs_male[16,6] %>% pull()

p_abs_male1_16 <- p_abs_male[17,6] %>% pull()
p_abs_male2_16 <- p_abs_male[18,6] %>% pull()
p_abs_male3_16 <- p_abs_male[19,6] %>% pull()
p_abs_male4_16 <- p_abs_male[20,6] %>% pull()

p_abs_male1_18 <- p_abs_male[21,6] %>% pull()
p_abs_male2_18 <- p_abs_male[22,6] %>% pull()
p_abs_male3_18 <- p_abs_male[23,6] %>% pull()
p_abs_male4_18 <- p_abs_male[24,6] %>% pull()

# PROPORTION OF FORMER DRINKERS IN MALES
p_form_male1_08 <- p_form_male[1,6] %>% pull()
p_form_male2_08 <- p_form_male[2,6] %>% pull()
p_form_male3_08 <- p_form_male[3,6] %>% pull()
p_form_male4_08 <- p_form_male[4,6] %>% pull()

p_form_male1_10 <- p_form_male[5,6] %>% pull()
p_form_male2_10 <- p_form_male[6,6] %>% pull()
p_form_male3_10 <- p_form_male[7,6] %>% pull()
p_form_male4_10 <- p_form_male[8,6] %>% pull()

p_form_male1_12 <- p_form_male[9,6] %>% pull()
p_form_male2_12 <- p_form_male[10,6] %>% pull()
p_form_male3_12 <- p_form_male[11,6] %>% pull()
p_form_male4_12 <- p_form_male[12,6] %>% pull()

p_form_male1_14 <- p_form_male[13,6] %>% pull()
p_form_male2_14 <- p_form_male[14,6] %>% pull()
p_form_male3_14 <- p_form_male[15,6] %>% pull()
p_form_male4_14 <- p_form_male[16,6] %>% pull()

p_form_male1_16 <- p_form_male[17,6] %>% pull()
p_form_male2_16 <- p_form_male[18,6] %>% pull()
p_form_male3_16 <- p_form_male[19,6] %>% pull()
p_form_male4_16 <- p_form_male[20,6] %>% pull()

p_form_male1_18 <- p_form_male[21,6] %>% pull()
p_form_male2_18 <- p_form_male[22,6] %>% pull()
p_form_male3_18 <- p_form_male[23,6] %>% pull()
p_form_male4_18 <- p_form_male[24,6] %>% pull()

# PROPORTION OF LIFE TIME ABSTAINERS IN FEMALES
p_abs_fem1_08 <- p_abs_fem[1,6] %>% pull()
p_abs_fem2_08 <- p_abs_fem[2,6] %>% pull()
p_abs_fem3_08 <- p_abs_fem[3,6] %>% pull()
p_abs_fem4_08 <- p_abs_fem[4,6] %>% pull()

p_abs_fem1_10 <- p_abs_fem[5,6] %>% pull()
p_abs_fem2_10 <- p_abs_fem[6,6] %>% pull()
p_abs_fem3_10 <- p_abs_fem[7,6] %>% pull()
p_abs_fem4_10 <- p_abs_fem[8,6] %>% pull()

p_abs_fem1_12 <- p_abs_fem[9,6] %>% pull()
p_abs_fem2_12 <- p_abs_fem[10,6] %>% pull()
p_abs_fem3_12 <- p_abs_fem[11,6] %>% pull()
p_abs_fem4_12 <- p_abs_fem[12,6] %>% pull()

p_abs_fem1_14 <- p_abs_fem[13,6] %>% pull()
p_abs_fem2_14 <- p_abs_fem[14,6] %>% pull()
p_abs_fem3_14 <- p_abs_fem[15,6] %>% pull()
p_abs_fem4_14 <- p_abs_fem[16,6] %>% pull()

p_abs_fem1_16 <- p_abs_fem[17,6] %>% pull()
p_abs_fem2_16 <- p_abs_fem[18,6] %>% pull()
p_abs_fem3_16 <- p_abs_fem[19,6] %>% pull()
p_abs_fem4_16 <- p_abs_fem[20,6] %>% pull()

p_abs_fem1_18 <- p_abs_fem[21,6] %>% pull()
p_abs_fem2_18 <- p_abs_fem[22,6] %>% pull()
p_abs_fem3_18 <- p_abs_fem[23,6] %>% pull()
p_abs_fem4_18 <- p_abs_fem[24,6] %>% pull()

# PROPORTION OF FORMER DRINKERS IN FEMALES
p_form_fem1_08 <- p_form_fem[1,6] %>% pull()
p_form_fem2_08 <- p_form_fem[2,6] %>% pull()
p_form_fem3_08 <- p_form_fem[3,6] %>% pull()
p_form_fem4_08 <- p_form_fem[4,6] %>% pull()

p_form_fem1_10 <- p_form_fem[5,6] %>% pull()
p_form_fem2_10 <- p_form_fem[6,6] %>% pull()
p_form_fem3_10 <- p_form_fem[7,6] %>% pull()
p_form_fem4_10 <- p_form_fem[8,6] %>% pull()

p_form_fem1_12 <- p_form_fem[9,6] %>% pull()
p_form_fem2_12 <- p_form_fem[10,6] %>% pull()
p_form_fem3_12 <- p_form_fem[11,6] %>% pull()
p_form_fem4_12 <- p_form_fem[12,6] %>% pull()

p_form_fem1_14 <- p_form_fem[13,6] %>% pull()
p_form_fem2_14 <- p_form_fem[14,6] %>% pull()
p_form_fem3_14 <- p_form_fem[15,6] %>% pull()
p_form_fem4_14 <- p_form_fem[16,6] %>% pull()

p_form_fem1_16 <- p_form_fem[17,6] %>% pull()
p_form_fem2_16 <- p_form_fem[18,6] %>% pull()
p_form_fem3_16 <- p_form_fem[19,6] %>% pull()
p_form_fem4_16 <- p_form_fem[20,6] %>% pull()

p_form_fem1_18 <- p_form_fem[21,6] %>% pull()
p_form_fem2_18 <- p_form_fem[22,6] %>% pull()
p_form_fem3_18 <- p_form_fem[23,6] %>% pull()
p_form_fem4_18 <- p_form_fem[24,6] %>% pull()

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
  summarise(gamma_fit = list(fit_gamma_with_density(volajohdia, x_vals))) %>%
  mutate(shape = map_dbl(gamma_fit, "shape"),
         rate = map_dbl(gamma_fit, "rate")) %>%
  ungroup() %>%
  dplyr::select(-gamma_fit)

# TRAPEZOIDAL INTEGRATION FUNCTION
trap_int <- function(x, shape, rate, rr, prop_abs, rr_form, prop_form) {
  
  # Interval width
  dx <- x[2] - x[1]
  y = dgamma(x_vals, shape = shape, rate = rate)
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

####################################
# ESTIMATING AAF FOR BREAST CANCER #
####################################
# Calculate the relative risk for each x
rr_bcan <- rr_linear(x_vals, b1_bcan)

paf_bc1_08 <- trap_int(x = x_vals, shape = input_gamma[[25,4]],
                    rate = input_gamma[[25,5]], rr = rr_bcan, 
                    prop_abs = female_props_df[1,4],
                    rr_form = 1, prop_form = female_props_df[2,4])

paf_bc2_08 <- trap_int(x = x_vals, 
                       shape = input_gamma[[26,4]],
                       rate = input_gamma[[26,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[3,4],
                       rr_form = 1,
                       prop_form = female_props_df[4,4])

paf_bc3_08 <- trap_int(x = x_vals, 
                       shape = input_gamma[[27,4]],
                       rate = input_gamma[[27,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[5,4],
                       rr_form = 1,
                       prop_form = female_props_df[6,4])

paf_bc4_08 <- trap_int(x = x_vals, 
                       shape = input_gamma[[28,4]],
                       rate = input_gamma[[28,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[7,4],
                       rr_form = 1,
                       prop_form = female_props_df[8,4])

paf_bc1_10 <- trap_int(x = x_vals, shape = input_gamma[[29,4]],
                       rate = input_gamma[[29,5]], rr = rr_bcan, 
                       prop_abs = female_props_df[9,4],
                       rr_form = 1, prop_form = female_props_df[10,4])

paf_bc2_10 <- trap_int(x = x_vals, 
                       shape = input_gamma[[30,4]],
                       rate = input_gamma[[30,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[11,4],
                       rr_form = 1,
                       prop_form = female_props_df[12,4])

paf_bc3_10 <- trap_int(x = x_vals, 
                       shape = input_gamma[[31,4]],
                       rate = input_gamma[[31,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[13,4],
                       rr_form = 1,
                       prop_form = female_props_df[14,4])

paf_bc4_10 <- trap_int(x = x_vals, 
                       shape = input_gamma[[32,4]],
                       rate = input_gamma[[32,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[15,4],
                       rr_form = 1,
                       prop_form = female_props_df[16,4])

paf_bc1_12 <- trap_int(x = x_vals, shape = input_gamma[[33,4]],
                       rate = input_gamma[[33,5]], rr = rr_bcan, 
                       prop_abs = female_props_df[17,4],
                       rr_form = 1, prop_form = female_props_df[18,4])

paf_bc2_12 <- trap_int(x = x_vals, 
                       shape = input_gamma[[34,4]],
                       rate = input_gamma[[34,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[19,4],
                       rr_form = 1,
                       prop_form = female_props_df[20,4])

paf_bc3_12 <- trap_int(x = x_vals, 
                       shape = input_gamma[[35,4]],
                       rate = input_gamma[[35,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[21,4],
                       rr_form = 1,
                       prop_form = female_props_df[22,4])

paf_bc4_12 <- trap_int(x = x_vals, 
                       shape = input_gamma[[36,4]],
                       rate = input_gamma[[36,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[23,4],
                       rr_form = 1,
                       prop_form = female_props_df[24,4])

paf_bc1_14 <- trap_int(x = x_vals, shape = input_gamma[[37,4]],
                       rate = input_gamma[[37,5]], rr = rr_bcan, 
                       prop_abs = female_props_df[25,4],
                       rr_form = 1, prop_form = female_props_df[26,4])

paf_bc2_14 <- trap_int(x = x_vals, 
                       shape = input_gamma[[38,4]],
                       rate = input_gamma[[38,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[27,4],
                       rr_form = 1,
                       prop_form = female_props_df[28,4])

paf_bc3_14 <- trap_int(x = x_vals, 
                       shape = input_gamma[[39,4]],
                       rate = input_gamma[[39,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[29,4],
                       rr_form = 1,
                       prop_form = female_props_df[30,4])

paf_bc4_14 <- trap_int(x = x_vals, 
                       shape = input_gamma[[40,4]],
                       rate = input_gamma[[40,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[31,4],
                       rr_form = 1,
                       prop_form = female_props_df[32,4])

paf_bc1_16 <- trap_int(x = x_vals, shape = input_gamma[[41,4]],
                       rate = input_gamma[[41,5]], rr = rr_bcan, 
                       prop_abs = female_props_df[33,4],
                       rr_form = 1, prop_form = female_props_df[34,4])

paf_bc2_16 <- trap_int(x = x_vals, 
                       shape = input_gamma[[42,4]],
                       rate = input_gamma[[42,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[35,4],
                       rr_form = 1,
                       prop_form = female_props_df[36,4])

paf_bc3_16 <- trap_int(x = x_vals, 
                       shape = input_gamma[[43,4]],
                       rate = input_gamma[[43,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[37,4],
                       rr_form = 1,
                       prop_form = female_props_df[38,4])

paf_bc4_16 <- trap_int(x = x_vals, 
                       shape = input_gamma[[44,4]],
                       rate = input_gamma[[44,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[39,4],
                       rr_form = 1,
                       prop_form = female_props_df[40,4])

paf_bc1_18 <- trap_int(x = x_vals, shape = input_gamma[[45,4]],
                       rate = input_gamma[[45,5]], rr = rr_bcan, 
                       prop_abs = female_props_df[41,4],
                       rr_form = 1, prop_form = female_props_df[42,4])

paf_bc2_18 <- trap_int(x = x_vals, 
                       shape = input_gamma[[46,4]],
                       rate = input_gamma[[46,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[43,4],
                       rr_form = 1,
                       prop_form = female_props_df[44,4])

paf_bc3_18 <- trap_int(x = x_vals, 
                       shape = input_gamma[[47,4]],
                       rate = input_gamma[[47,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[45,4],
                       rr_form = 1,
                       prop_form = female_props_df[46,4])

paf_bc4_18 <- trap_int(x = x_vals, 
                       shape = input_gamma[[48,4]],
                       rate = input_gamma[[48,5]],
                       rr = rr_bcan, 
                       prop_abs = female_props_df[47,4],
                       rr_form = 1,
                       prop_form = female_props_df[48,4])

