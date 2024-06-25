#########################################################
# CLEANING AND MERGING NATIONAL DRUG AND ALCOHOL SURVEY #
#########################################################

# LIBRARIES
library(dplyr)
library(haven)
library(readr)

################
# FIRST FILTER #
################

# 2008
enpg08 <- read_dta("Data/Encuesta Nacional de Drogas/BDENPG2008/BD ENPG 2008/Base de datos ENPG 2008 (PG).DTA.dta")
data08 <- enpg08 %>% mutate(year = 2008)
data08 <- data08 %>% 
  select(id, year, region, comuna,exp, sexo, edad,  
         trabajo = q264, religion = q277, nedu = q279, ecivil = q281, 
         p_originario = q282, ingreso = q284,oh1 = q12, audit1 = q18, audit2 = q19, audit3 = q20)
write_rds(data08, "Data/Encuesta Nacional de Drogas/BDENPG2008/BD ENPG 2008/BDENPG2008.RDS", compress = "gz")

# 2010
enpg10 <- read_dta("Data/Encuesta Nacional de Drogas/BDENPG2010/BD ENPG 2010 (Stata y SPSS)/Base de datos ENPG 2010 (PG).DTA.dta")
data10 <- enpg10 %>% mutate(year = 2010) %>% 
  select(id = folio, year,region = pregion, nombre_comuna = pnomcomun, comuna = pcodcom, 
         exp = factor_ajustado_com, sexo = psexoent,edad = pedadent, nedu = p288a, 
         trabajo = p268, religion = p285,
         ecivil = p282, p_originario = p283, ingreso = p291,
         oh1 = p013, audit1 = p019, audit2 = p020, audit3 = p021)
write_rds(data10, "Data/Encuesta Nacional de Drogas/BDENPG2010/BD ENPG 2010 (Stata y SPSS)/BDENPG2010.RDS", compress = "gz")

# BASE 2012 
enpg12 <- read_dta("Data/Encuesta Nacional de Drogas/BDENPG2012/Base de datos ENPG 2012 (PG).dta")
data12 <- enpg12 %>% mutate(year = 2012) %>% 
  select(id = idencuesta, year, region, nombre_comuna = comuna, comuna = codigo_comuna, 
         exp = PONDERADOR, sexo, edad, nedu1 = p184_1, nedu2 = p185_1, trabajo = p174, 
         religion = p180, ecivil = p177, p_originario = p179, ingreso = p188,
         oh1 = p10, audit1 = p19, audit2 = p20, audit3 = p21)
write_rds(data12, "Data/Encuesta Nacional de Drogas/BDENPG2012/BDENPG2012.RDS", compress = "gz")

# BASE 2014
enpg14 <- read_sav("Data/Encuesta Nacional de Drogas/BDENPG2014/Base de datos ENPG 2014 (PG).sav")
data14 <- enpg14 %>% mutate(year = 2014) %>% 
  select(id = idencuesta, year, region = Region, comuna = Nombre_Comuna, exp = RND_F2_MAY_AJUS_com, sexo, edad,
         nedu1 = dp9_1, nedu2 = dp10_1, trabajo = co7, religion = dp5, ecivil = dp2, 
         p_originario = dp4, ingreso = dp13, oh1 = oh1, audit1 = oh10, audit2 = oh11,
         audit3 = oh12)
write_rds(data14, "Data/Encuesta Nacional de Drogas/BDENPG2014/BDENPG2014.RDS", compress = "gz")

# BASE 2016
enpg16 <- read_sav("Data/Encuesta Nacional de Drogas/BDENPG2016/BD ENPG 2016 (Stata y SPSS)/base ENPG 2016 publico general.sav")
exp16 <- read_sav("Data/Encuesta Nacional de Drogas/BDENPG2016/BD ENPG 2016 (Stata y SPSS)/factoresdeexpansion.sav")
data16 <- enpg16 %>% 
  left_join(exp16, by = "idencuesta") %>% 
  mutate(year = 2016) %>% 
  select(id = idencuesta, year, region = región, comuna, exp = Fexp, sexo, edad,
         nedu1 = dp_9_a, nedu2 = dp_10_a, trabajo = co_7, religion = dp_5, ecivil = dp_2, 
         p_originario = dp_4, ingreso = dp_13, oh1 = oh_1, audit1 = oh_14, audit2 = oh_15,
         audit3 = oh_16)
write_rds(data16, "Data/Encuesta Nacional de Drogas/BDENPG2016/BDENPG2016.RDS", compress = "gz")

# BASE 2018
enpg18 <- read_sav("Data/Encuesta Nacional de Drogas/BDENPG2018/BD ENPG 2018 (Stata y SPSS)/Base de datos ENPG 2018 (PG).sav")
data18 <- enpg18 %>% 
  mutate(year = 2018) %>% 
  select(id = SbjNum, year, region = Region, comuna, exp = Fexp, sexo = S01, edad = S02,
         nedu1 = T_DP_12_1, nedu2 = T_DP_13_1, trabajo = CO_7, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, oh1 = OH_1, audit1 = OH_12, audit2 = OH_13,
         audit3 = OH_14)
write_rds(data18, "Data/Encuesta Nacional de Drogas/BDENPG2018/BDENPG2018.RDS", compress = "gz")

# BASE 2020
enpg20 <- read_dta("Data/Encuesta Nacional de Drogas/BDENPG2020/Base publico general enpg 2020 (stata 12 o inferior).dta")
data20 <- enpg20 %>% 
  mutate(year = 2020) %>% 
  select(id = SbjNum, year, region = REGION, comuna = Nom_comuna, exp = FACT_PERS_COMUNA, sexo = S01, edad = S02,
         nedu1 = DP_12, nedu2 = DP_13, trabajo = CO_6, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, oh1 = OH_1, audit1 = OH_8, audit2 = OH_9,
         audit3 = OH_10)
write_rds(data20, "Data/Encuesta Nacional de Drogas/BDENPG2020/BDENPG2020.RDS", compress = "gz")

# BASE 2022
enpg22 <- read_dta("Data/Encuesta Nacional de Drogas/BDENPG2022/Base de datos publica ENPG 2022 (Stata 13).dta")
data22 <- enpg22 %>% 
  mutate(year = 2022) %>% 
  select(id = FOLIO, year, region = REGION, comuna = COD_COMUNA, exp = FACTOR_EXPANSION, sexo = SEXO, edad = EDAD,
         nedu1 = DP_12, nedu2 = DP_13, trabajo = CO_6, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, oh1 = OH_1, audit1 = OH_8, audit2 = OH_9,
         audit3 = OH_10)

write_rds(data22, "Data/Encuesta Nacional de Drogas/BDENPG2022/BDENPG2022.RDS", compress = "gz")

#########
# MERGE #
#########
data08 <- readRDS("Data/Encuesta Nacional de Drogas/BDENPG2008/BD ENPG 2008/BDENPG2008.RDS") 
data10 <- readRDS("C:/Users/Jose/Desktop/RESEARCH/Álvaro/FONDECYT/Data/Encuesta Nacional de Drogas/BDENPG2010/BD ENPG 2010 (Stata y SPSS)/BDENPG2010.RDS")
data12 <- readRDS("C:/Users/Jose/Desktop/RESEARCH/Álvaro/FONDECYT/Data/Encuesta Nacional de Drogas/BDENPG2012/BDENPG2012.RDS")
data14 <- readRDS("C:/Users/Jose/Desktop/RESEARCH/Álvaro/FONDECYT/Data/Encuesta Nacional de Drogas/BDENPG2014/BDENPG2014.RDS")
data16 <- readRDS("C:/Users/Jose/Desktop/RESEARCH/Álvaro/FONDECYT/Data/Encuesta Nacional de Drogas/BDENPG2016/BDENPG2016.RDS")
data18 <- readRDS("C:/Users/Jose/Desktop/RESEARCH/Álvaro/FONDECYT/Data/Encuesta Nacional de Drogas/BDENPG2018/BDENPG2018.RDS")
data20 <- readRDS("C:/Users/Jose/Desktop/RESEARCH/Álvaro/FONDECYT/Data/Encuesta Nacional de Drogas/BDENPG2020/BDENPG2020.RDS")
data22 <- readRDS("C:/Users/Jose/Desktop/RESEARCH/Álvaro/FONDECYT/Data/Encuesta Nacional de Drogas/BDENPG2022/BDENPG2022.RDS")

# BASE 2008
data08 <- readRDS("Data/Encuesta Nacional de Drogas/BDENPG2008/BD ENPG 2008/BDENPG2008.RDS")
new_names <- paste("oh", 1:21, sep = "_")  
names_to_change <- paste("q", 12:32, sep = "")  
all_names <- names(data08)
all_names[which(all_names %in% names_to_change)] <- new_names
names(data08) <- all_names

data08 <- data08 %>%
  mutate(across(where(is.numeric), ~na_if(., 99)), 
         across(where(is.character), ~na_if(., "99")),
         across(.cols = starts_with("oh_17"):starts_with("oh_21"), 
                       .fns = ~na_if(., 9)),
         across(.cols = c("oh_1", "oh_2", "oh_4", "oh_6"), 
                .fns = ~na_if(., 9)),
         oh_5 = ifelse(oh_5 >30, NA, oh_5)) 
data08 <- data08[,-c(27:37)]
write_rds(data08, "Data/Encuesta Nacional de Drogas/BDENPG2008/BD ENPG 2008/BDENPG2008_clean.RDS", compress = "gz")

# BASE 2010
data10 <- readRDS("C:/Users/Jose/Desktop/RESEARCH/Álvaro/FONDECYT/Data/Encuesta Nacional de Drogas/BDENPG2010/BD ENPG 2010 (Stata y SPSS)/BDENPG2010.RDS")
new_names <- paste("oh", 1:20, sep = "_")  
names_to_change <- paste("p",0, 13:32, sep = "")  
all_names <- names(data10)
all_names[which(all_names %in% names_to_change)] <- new_names
names(data10) <- all_names
data10 <- data10 %>% 
  mutate(across(.cols = c("oh_1", "oh_2", "oh_3", "oh_4","oh_5","oh_6",
                          "oh_17","oh_18","oh_19","oh_20"), 
                .fns = ~na_if(., 0))) %>% 
  rename(sexo = psexoent, edad = pedadent,region = pregion, exp = factor_ajustado_com)
write_rds(data10, "Data/Encuesta Nacional de Drogas/BDENPG2010/BD ENPG 2010 (Stata y SPSS)/BDENPG2010_clean.RDS", compress = "gz")

# BASE 2012

