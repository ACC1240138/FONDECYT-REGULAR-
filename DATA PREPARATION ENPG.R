#########################################################
# CLEANING AND MERGING NATIONAL DRUG AND ALCOHOL SURVEY #
#########################################################

# LIBRARIES
library(dplyr)
library(haven)
library(readr)

#############
# ENPG 2008 #
#############

enpg08 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2008.RDS")
enpg08 <- enpg08 %>% mutate(year = 2008)
data08 <- enpg08 %>% 
  select(id, year, area, region,comuna, n_per, exp, sexo, edad,  
         religion = q277, nedu = q279, ecivil = q281, 
         p_originario = q282, ingreso = q284,oh1 = q12, oh2 = q15, oh3 = q16,
         tab1 = q5,tab2 = q8, tab3 = q9, tab4 = q10, mar1 = q33, mar2 = q36,
         coc1 = q79, coc2 = q82, q133a, q133b, q133c, q133d, q133e, q133f, q133g,
         q133h, q133i, q133j, q134, q135, q136, q137, q138, q139, q140, q141, q142, q143,
         audit1 = q18, audit2 = q19, audit3 = q20) %>% 
  mutate(id = factor(id),
         mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
         mar2 = factor(mar2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
         tab2 = factor(tab2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
    sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(ecivil <= 2 ~ "soltero",
                                   ecivil == 3 ~ "casado",
                                   ecivil == 4 | ecivil == 5 ~ "separado",
                                   ecivil == 6 | ecivil == 7 ~ "viudo")),
         religion = factor(religion, levels = c(1:4), labels = c("Catolico",
                                                                 "evangelico/protestante",
                                                                 "otra",
                                                                 "ninguna")),
         p_originario = factor(ifelse(p_originario == 9,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu == 1 ~ "basica incompleta",
                                 nedu == 2 ~ "basica completa",
                                 nedu == 3 ~ "media incompleta",
                                 nedu == 4 ~ "media completa",
                                 nedu == 5 | nedu == 7 ~ "superior incompleta",
                                 nedu == 6 | nedu == 8 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 9),
         ingreso = na_if(ingreso, 0))
         

#############
# ENPG 2010 #
#############

enpg10 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2010.RDS")
data10 <- enpg10 %>% mutate(year = 2010) %>% 
  select(id = folio, year,region = pregion,
         exp = factor_ajustado_com, sexo = psexoent,edad = pedadent, nedu = p288a, 
         religion = p285,ecivil = p282, p_originario = p283, ingreso = p291,
         oh1 = p013, oh2 = p016, audit1 = p019, audit2 = p020, audit3 = p021) %>% 
  mutate(id = factor(id),
    sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(ecivil == 0 ~ NA,
                                   ecivil <= 2 ~ "soltero",
                                   ecivil == 3 ~ "casado",
                                   ecivil == 4 | ecivil == 5 ~ "separado",
                                   ecivil == 6 | ecivil == 7 ~ "viudo")),
         religion = factor(religion, levels = c(1:4), labels = c("Catolico",
                                                                 "evangelico/protestante",
                                                                 "otra",
                                                                 "ninguna")),
         p_originario = factor(ifelse(p_originario == 9,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu < 2 | nedu == 16 ~ "basica incompleta",
                                 nedu > 2 & nedu < 5 ~ "basica completa",
                                 nedu >= 5 & nedu <= 8 ~ "media incompleta",
                                 nedu == 4 ~ "media completa",
                                 nedu == 9 | nedu == 11 | nedu == 13 ~ "superior incompleta",
                                 nedu == 10 | nedu == 12 | nedu == 14 | nedu == 15 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 9))


#############
# ENPG 2012 #
#############
enpg12 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2012.RDS")
data12 <- enpg12 %>% mutate(year = 2012) %>% 
  select(id = idencuesta, year, region = "región",
         exp = PONDERADOR, sexo, edad, nedu1 = p184_1, nedu2 = p185_1, 
         religion = p180, ecivil = p177, p_originario = p179, ingreso = p188,
         oh1 = p10, oh2 = p13,audit1 = p19, audit2 = p20, audit3 = p21) %>% 
mutate(id = factor(id),
  sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
       oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
       oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                            "2 a 3 veces a la semana","4 o mas veces a la semana")),
       audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                            "7-8", "9 o mas")),
       audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                            "mensualmente",
                                                            "semanalmente",
                                                            "todos o casi todos los dias")),
       ecivil = factor(case_when(
                                 ecivil == 1 ~ "soltero",
                                 ecivil == 2 ~ "casado",
                                 ecivil == 4  ~ "viudo",
                                 ecivil == 3 | ecivil == 5 ~ "separado")),
       religion = factor(case_when(religion == 1 ~ "Catolico",
                                   religion == 2 ~ "evangelico/protestante",
                                   religion >= 3 & religion <= 11 ~ "otra",
                                   religion == 12 ~ "ninguna")),
       p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
       nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                               nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                               nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                               TRUE~NA)),
       ingreso = na_if(ingreso, 88),
       ingreso = na_if(ingreso, 99)) %>% 
  select(-nedu1, -nedu2)

#############
# ENPG 2014 #
#############
enpg14 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2014.RDS")
data14 <- enpg14 %>% mutate(year = 2014) %>% 
  select(id = idencuesta, year, region = Region, exp = RND_F2_MAY_AJUS_com, sexo, edad,
         nedu1 = dp9_1, nedu2 = dp10_1, religion = dp5, ecivil = dp2, 
         p_originario = dp4, ingreso = dp13, oh1,oh2 =oh4,  audit1 = oh10, audit2 = oh11,
         audit3 = oh12) %>% 
mutate(id = factor(id),
  sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
       oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
       oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                            "2 a 3 veces a la semana","4 o mas veces a la semana")),
       audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                            "7-8", "9 o mas")),
       audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                            "mensualmente",
                                                            "semanalmente",
                                                            "todos o casi todos los dias")),
       ecivil = factor(case_when(
         ecivil == 1 ~ "soltero",
         ecivil == 2 ~ "casado",
         ecivil == 4  ~ "viudo",
         ecivil == 3 | ecivil == 5 ~ "separado")),
       religion = factor(case_when(religion == 1 ~ "Catolico",
                                   religion == 2 ~ "evangelico/protestante",
                                   religion >= 3 & religion <= 11 ~ "otra",
                                   religion == 12 ~ "ninguna")),
       p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
       nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                               nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                               nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                               TRUE~NA)),
       ingreso = na_if(ingreso, 88),
       ingreso = na_if(ingreso, 99)) %>% 
  select(-nedu1, -nedu2)

#############
# ENPG 2016 #
#############
enpg16 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2016.RDS")
exp16 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/Expansion16.RDS")

data16 <- enpg16 %>% 
  left_join(exp16, by = "idencuesta") %>% 
  mutate(year = 2016) %>% 
  select(id = idencuesta, year, region = región, exp = Fexp.x, sexo, edad,
         nedu1 = dp_9_a, nedu2 = dp_10_a, religion = dp_5, ecivil = dp_2, 
         p_originario = dp_4, ingreso = dp_13, oh1 = oh_1, oh2 = oh_4,audit1 = oh_14, audit2 = oh_15,
         audit3 = oh_16) %>% 
mutate(id = factor(id),
  sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
       oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
       oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                            "2 a 3 veces a la semana","4 o mas veces a la semana")),
       audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                            "7-8", "9 o mas")),
       audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                            "mensualmente",
                                                            "semanalmente",
                                                            "todos o casi todos los dias")),
       ecivil = factor(case_when(
         ecivil == 1 ~ "soltero",
         ecivil == 2 & ecivil == 6 ~ "casado",
         ecivil == 4  ~ "viudo",
         ecivil == 3 | ecivil == 5 ~ "separado")),
       religion = factor(case_when(religion == 1 ~ "Catolico",
                                   religion == 2 ~ "evangelico/protestante",
                                   religion >= 3 & religion <= 11 ~ "otra",
                                   religion == 12 ~ "ninguna")),
       p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
       nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                               nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                               nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                               TRUE~NA)),
       ingreso = na_if(ingreso, 88),
       ingreso = na_if(ingreso, 99)) %>% 
  select(-nedu1, -nedu2)

#############
# ENPG 2018 #
#############
enpg18 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2018.RDS")
data18 <- enpg18 %>% 
  mutate(year = 2018) %>% 
  select(id = SbjNum, year, region = Region, exp = Fexp, sexo = S01, edad = S02,
         nedu1 = T_DP_12_1, nedu2 = T_DP_13_1, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, oh1 = OH_1, oh2 = OH_4, audit1 = OH_12, audit2 = OH_13,
         audit3 = OH_14) %>% 
  mutate(id = factor(id),
    sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(
           ecivil == 1 ~ "soltero",
           ecivil == 2 & ecivil == 6 ~ "casado",
           ecivil == 4  ~ "viudo",
           ecivil == 3 | ecivil == 5 ~ "separado")),
         religion = factor(case_when(religion == 1 ~ "Catolico",
                                     religion == 2 ~ "evangelico/protestante",
                                     religion >= 3 & religion <= 11 ~ "otra",
                                     religion == 12 ~ "ninguna")),
         p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                                 nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                                 nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 88),
         ingreso = na_if(ingreso, 99)) %>% 
  select(-nedu1, -nedu2)

#############
# ENPG 2020 #
#############
enpg20 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2020.RDS")
data20 <- enpg20 %>% 
  mutate(year = 2020) %>% 
  select(id = SbjNum, year, region = REGION, exp = FACT_PERS_COMUNA, sexo = S01, edad = S02,
         nedu1 = DP_12, nedu2 = DP_13, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, oh1 = OH_1, oh2 = OH_4, audit1 = OH_8, audit2 = OH_9,
         audit3 = OH_10) %>% 
  mutate(id = factor(id),
    sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(
           ecivil == 1 ~ "soltero",
           ecivil == 2 & ecivil == 6 ~ "casado",
           ecivil == 4  ~ "viudo",
           ecivil == 3 | ecivil == 5 ~ "separado")),
         religion = factor(case_when(religion == 1 ~ "Catolico",
                                     religion == 2 ~ "evangelico/protestante",
                                     religion >= 3 & religion <= 11 ~ "otra",
                                     religion == 12 ~ "ninguna")),
         p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                                 nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                                 nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 88),
         ingreso = na_if(ingreso, 99)) %>% 
  select(-nedu1, -nedu2)

#############
# ENPG 2022 #
#############
enpg22 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2022.RDS")
data22 <- enpg22 %>% 
  mutate(year = 2022) %>% 
  select(id = FOLIO, year, region = REGION, exp = FACTOR_EXPANSION, sexo = SEXO, edad = EDAD,
         nedu1 = DP_12, nedu2 = DP_13, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, oh1 = OH_1, oh2 = OH_4, audit1 = OH_8, audit2 = OH_9,
         audit3 = OH_10) %>% 
  mutate(id = factor(id),
    sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(
           ecivil == 1 ~ "soltero",
           ecivil == 2 & ecivil == 6 ~ "casado",
           ecivil == 4  ~ "viudo",
           ecivil == 3 | ecivil == 5 ~ "separado")),
         religion = factor(case_when(religion == 1 ~ "Catolico",
                                     religion == 2 ~ "evangelico/protestante",
                                     religion >= 3 & religion <= 11 ~ "otra",
                                     religion == 12 ~ "ninguna")),
         p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                                 nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                                 nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 88),
         ingreso = na_if(ingreso, 99)) %>% 
  select(-nedu1, -nedu2)

#########
# MERGE #
#########

enpg_full <- bind_rows(data08, data10, data12,data14,data16,data18,data20,data22)

write_rds(enpg_full, "ENPG_FULL.RDS", compress = "gz")


########################
# ESTIMATE ALCOHOL USE #
########################

enpg_full <- 