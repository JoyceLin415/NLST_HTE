# 072720
rm(list = ls())

library(here)
library(survival)

dat0<-read.csv("NSLT_participants_dat0.csv")
dat<-dat0

#---- Data management & summary statistics ----#
names(dat)
table(dat$rndgroup, useNA = "ifany") # 26722 CT, 26730 X-ray

# Num of ineligible obs
table(dat$ineligible, dat$elig, useNA = "ifany") # 53248 elig, 204 inelig
table(dat$elig, dat$rndgroup, useNA = "ifany")
chisq.test(dat$elig, dat$rndgroup) # p = 0.3629
dat<-dat[dat$elig == "Eligible Participant",]
table(dat$rndgroup, useNA = "ifany") # 26627 CT, 26621 X-ray2


# Check num of obs with fup_days = 0 & remove them
str(dat$fup_days)
dat$fup_is_0<-ifelse(dat$fup_days==0, 1, 0)
table(dat$fup_is_0) # 58 fup = 0 (they all dead)
table(dat$deathstat, dat$fup_is_0)
chisq.test(dat$fup_is_0, dat$rndgroup) # p = 0.8623
dat<-dat[dat$fup_days!=0,] 
table(dat$rndgroup, useNA = "ifany") # 26593 CT, 26578 x-ray

# Num who received T0 screenig
table(dat$scr_res0)
table(dat$scr_res0, dat$rndgroup, useNA = "ifany")
dat$screen_t0_yn<-ifelse(substr(dat$scr_res0, 0, 3) %in% c("Neg", "Pos"), 1,
                         ifelse(substr(dat$scr_res0, 0, 3) %in% c("Not", "Ina"), 0, NA))
table(dat$screen_t0_yn, useNA = "ifany") # 52174 left, 1074 excluded
table(dat$screen_t0_yn, dat$rndgroup, useNA = "ifany")
chisq.test(dat$screen_t0_yn, dat$rndgroup) # p < 0.0001
dat<-dat[dat$screen_t0_yn == 1,]
table(dat$rndgroup, useNA = "ifany") # 26203 CT, 25913 x-ray


# time to T0 scr
str(dat$scr_days0)
dat$scr_days0_num<-as.numeric(as.character(dat$scr_days0))
tapply(dat$scr_days0_num, dat$rndgroup, summary)
wilcox.test(dat$scr_days0_num~dat$rndgroup) # p < 0.0001
t.test(log(dat$scr_days0_num+1)~dat$rndgroup, na.rm = T) # p < 0.0001
table(is.na(dat$scr_days0_num), dat$screen_t0_yn) # All patients receive screening at T0 have scr_days0_num

tmp<-dat[dat$scr_days0_num>=0 & is.na(dat$scr_days0_num) == 0,]
tapply(tmp$scr_days0_num, tmp$rndgroup, summary)
wilcox.test(tmp$scr_days0_num~tmp$rndgroup) # p < 0.0001

# Outcome
## All cause mortality
### Adjust follow up time (see NLST user guide)
dat$fup_days_adj_allcausemort<-ifelse(dat$deathstat == "No report of death", 
                                      ifelse(dat$rndgroup == "Spiral CT", dat$fup_days - 58.1772, # CT arm
                                             ifelse(dat$rndgroup == "X-ray", dat$fup_days - 58.7590, NA)), # x-ray arm
                                      dat$fup_days) # death (no change)
summary(dat$fup_days_adj_allcausemort)
tapply(dat$fup_days_adj_allcausemort, dat$deathstat, summary)
dat$fup_days_adj_allcausemort<-ifelse(dat$fup_days_adj_allcausemort < 0 & is.na(dat$fup_days_adj_allcausemort) == 0, 0, dat$fup_days_adj_allcausemort)
tapply(dat$fup_days_adj_allcausemort, dat$deathstat, summary)

### All cause death
dat$allcause_death<-ifelse(dat$deathstat %in% c("EVP certified", "Death Certificate coded", "Death Certificate received but not coded"), 1, 0)
table(dat$deathstat, dat$allcause_death)

### 8-yr mortality
maxtm<-360*8
dat$allcause_death_8yr<-ifelse(dat$fup_days_adj_allcausemort <= maxtm, dat$allcause_death,
                               ifelse(dat$fup_days_adj_allcausemort > maxtm, 0, NA))
table(dat$allcause_death_8yr, useNA = "ifany") # 3791 deaths

dat$allcause_death_8yr_time<-ifelse(dat$fup_days_adj_allcausemort <= maxtm, dat$fup_days_adj_allcausemort,
                               ifelse(dat$fup_days_adj_allcausemort > maxtm, maxtm, NA))
summary(dat$allcause_death_8yr_time, useNA = "ifany")

## Lung cancer mortality
### Adjust follow up time (see NLST user guide)
dat$fup_days_adj_lungmort<-ifelse(dat$deathcutoff %in% c("No death or no date of death", "Death Not Included"), 
                                      ifelse(dat$rndgroup == "Spiral CT", dat$fup_days - 394.6020, # CT arm
                                             ifelse(dat$rndgroup == "X-ray", dat$fup_days - 392.1746, NA)), # x-ray arm
                                      dat$fup_days) # death (no change)
summary(dat$fup_days_adj_lungmort)
tapply(dat$fup_days_adj_lungmort, dat$deathcutoff, summary)
dat$fup_days_adj_lungmort<-ifelse(dat$fup_days_adj_lungmort < 0 & is.na(dat$fup_days_adj_lungmort) == 0, 0, dat$fup_days_adj_lungmort)
tapply(dat$fup_days_adj_lungmort, dat$deathcutoff, summary)

### Lung cancer deaths
dat$lung_death<-ifelse(dat$finaldeathlc %in% c("Death due to lung cancer or work-up of suspected lung cancer") & 
                         dat$deathcutoff == "Death Included", 1, 0)
table(dat$finaldeathlc, dat$lung_death, dat$deathcutoff)

### 8-yr mortality
maxtm<-360*8
dat$lung_death_8yr<-ifelse(dat$fup_days_adj_lungmort <= maxtm, dat$lung_death,
                           ifelse(dat$fup_days_adj_lungmort > maxtm, 0, NA))
table(dat$lung_death_8yr, useNA = "ifany") # 769 deaths

dat$lung_death_8yr_time<-ifelse(dat$fup_days_adj_lungmort <= maxtm, dat$fup_days_adj_lungmort,
                                ifelse(dat$fup_days_adj_lungmort > maxtm, maxtm, NA))
summary(dat$lung_death_8yr_time, useNA = "ifany")

# Check num of lost FU
table(dat$lung_death, dat$allcause_death, useNA = "ifany")

# Cov
## Trt
dat$CT<-ifelse(dat$rndgroup=="Spiral CT", 1,
                         ifelse(dat$rndgroup=="X-ray",0,NA))
table(dat$rndgroup, dat$CT, useNA = "ifany")

## Age
summary(dat$age)

## Sex
table(dat$gender, useNA = "ifany") # 21795 female, 31376 male
dat$female<-ifelse(dat$gender == "Female", 1,
                   ifelse(dat$gender == "Male", 0, NA))
table(dat$female, dat$gender, useNA = "ifany")

## Race/Ethnic
table(dat$race, useNA = "ifany");table(dat$ethnic, useNA = "ifany")
fisher.test(dat$race, dat$rndgroup, simulate.p.value=TRUE) # p = 0.08
fisher.test(dat$ethnic, dat$rndgroup, simulate.p.value=TRUE) # p = 0.0009

dat$race_condensed<-ifelse(substr(dat$race,1,7) == "Missing" | 
                             dat$race %in% c("Participant refused to answer","Unknown/ decline to answer"), 
                           NA, as.character(dat$race))
table(dat$race_condensed, dat$race, useNA = "ifany")
table(dat$race_condensed, dat$rndgroup, useNA = "ifany")
fisher.test(dat$race_condensed, dat$rndgroup, simulate.p.value=TRUE) # p = 0.4053

dat$race_condensed2<-factor(ifelse(dat$race_condensed %in% c("American Indian or Alaskan Native",
                                                                       "Native Hawaiian or Other Pacific Islander",
                                                                       "More than one race"), "Other", 
                                        ifelse(dat$race_condensed == "Black or African-American", "Black", as.character(dat$race_condensed))))
dat$race_condensed2_num<-ifelse(dat$race_condensed2 == "White", 1,
                                ifelse(dat$race_condensed2 == "Black", 2,
                                       ifelse(dat$race_condensed2 == "Asian", 3,
                                              ifelse(dat$race_condensed2 == "Other", 4, NA))))
table(dat$race_condensed2, dat$race_condensed2_num, useNA = "ifany")

dat$ethnic_condensed<-ifelse(substr(dat$ethnic,1,7) == "Missing" | 
                               dat$ethnic %in% c("Participant refused to answer","Unknown/ decline to answer"), 
                             NA, as.character(dat$ethnic))
table(dat$ethnic_condensed)
table(dat$ethnic, dat$ethnic_condensed, useNA = "ifany")
fisher.test(dat$ethnic_condensed, dat$rndgroup, simulate.p.value=TRUE) # p = 0.0020

table(dat$race_condensed, dat$ethnic_condensed, useNA = "ifany")

## Smoking
### Current vs former
table(dat$cigsmok, dat$rndgroup, useNA = "ifany")

### Age quit
str(dat$age_quit)
dat$age_quit_num<-as.numeric(as.character(dat$age_quit))
tapply(dat$age_quit_num, dat$rndgroup, summary)
tapply(dat$age_quit_num, dat$cigsmok, summary)

### Yrs from cessation
dat$yrs_from_cessation0<-ifelse(dat$cigsmok=="Current",0,
                               ifelse(dat$cigsmok=="Former",dat$age-dat$age_quit_num, NA))
summary(dat$yrs_from_cessation0)
tapply(dat$yrs_from_cessation0, dat$cigsmok, summary)

dat$yrs_from_cessation<-ifelse(dat$yrs_from_cessation0<0 & is.na(dat$yrs_from_cessation0)==0, 0,
                               ifelse(dat$yrs_from_cessation0>15 & is.na(dat$yrs_from_cessation0)==0, NA, as.numeric(dat$yrs_from_cessation0)))
summary(dat$yrs_from_cessation)
tapply(dat$yrs_from_cessation, dat$cigsmok, summary)

tmp<-dat[dat$yrs_from_cessation<0 & is.na(dat$yrs_from_cessation)==0,] # 12 pt

tmp<-dat[dat$yrs_from_cessation>15 & is.na(dat$yrs_from_cessation)==0,] # 288 pt


### Age start
str(dat$smokeage)
dat$smokeage_num<-as.numeric(as.character(dat$smokeage))
tapply(dat$smokeage_num, dat$rndgroup, summary)



### Pipe
table(dat$pipe, dat$rndgroup, useNA = "ifany")
dat$pipe_YNNA<-ifelse(dat$pipe == "Missing", NA, as.character(dat$pipe))
table(dat$pipe_YNNA, dat$pipe, useNA = "ifany")

chisq.test(dat$pipe_YNNA, dat$rndgroup) # p = 0.0653

### Cigar
table(dat$cigar, dat$rndgroup, useNA = "ifany")
dat$cigar_YNNA<-ifelse(dat$cigar == "Missing", NA, as.character(dat$cigar))
table(dat$cigar_YNNA, dat$cigar, useNA = "ifany") # 153 NA

### pack per yr
str(dat$pkyr)
tapply(dat$pkyr, dat$rndgroup, summary)
tapply(dat$pkyr, dat$rndgroup, hist)
wilcox.test(dat$pkyr~dat$rndgroup) # p = 0.9653

### Smoke live
table(dat$smokelive, dat$rndgroup, useNA = "ifany")
dat$smokelive_YNNA<-ifelse(dat$smokelive == "Missing", NA, as.character(dat$smokelive))
table(dat$smokelive_YNNA, dat$smokelive, useNA = "ifany") # 158 NA

### Smoke work
table(dat$smokework, dat$rndgroup, useNA = "ifany")
dat$smokework_YNNA<-ifelse(dat$smokework == "Missing", NA, as.character(dat$smokework))
table(dat$smokework_YNNA, dat$smokework, useNA = "ifany") # 259 NA

### Family hist
dat$fam_hist_lung<-ifelse(dat$famchild=="Yes" |dat$famfather=="Yes" | dat$fammother=="Yes"| 
                            dat$famsister=="Yes" | dat$fambrother=="Yes", "Yes",
                          ifelse(dat$famchild=="No" & dat$famfather=="No" & dat$fammother=="No" & 
                                   dat$famsister=="No" & dat$fambrother=="No", "No", NA))
table(dat$fam_hist_lung, useNA = "ifany") # 1544 missing
dat$fam_hist_lung_num<-ifelse(dat$fam_hist_lung=="Yes",1,
                              ifelse(dat$fam_hist_lung=="No",0,NA))
table(dat$fam_hist_lung_num, dat$fam_hist_lung, useNA = "ifany")


## Working hist and masks
dat$wrk_hist<-ifelse(dat$wrkasbe=="Yes"|dat$wrkbaki=="Yes"|dat$wrkbutc=="Yes"|
                       dat$wrkchem=="Yes"|dat$wrkcoal=="Yes"|dat$wrkcott=="Yes"|dat$wrkfarm=="Yes"|
                       dat$wrkfire=="Yes"|dat$wrkflou=="Yes"|dat$wrkfoun=="Yes"|dat$wrkhard=="Yes"|
                       dat$wrkpain=="Yes"|dat$wrksand=="Yes"|dat$wrkweld=="Yes", "Yes", 
                     ifelse(dat$wrkasbe=="No"&dat$wrkbaki=="No"&dat$wrkbutc=="No"&
                              dat$wrkchem=="No"&dat$wrkcoal=="No"&dat$wrkcott=="No"&dat$wrkfarm=="No"&
                              dat$wrkfire=="No"&dat$wrkflou=="No"&dat$wrkfoun=="No"&dat$wrkhard=="No"&
                              dat$wrkpain=="No"&dat$wrksand=="No"&dat$wrkweld=="No", "No", NA))
table(dat$wrk_hist, useNA = "ifany") # 88 NA

# BMI
dat$weight<-as.numeric(as.character(dat$weight))
dat$height<-as.numeric(as.character(dat$height))
dat$BMI<-as.numeric(as.character(dat$weight))/(as.numeric(as.character(dat$height))^2)*703

## Diag hist
dat[,c("diagadas_YNNA", "diagasbe_YNNA", "diagbron_YNNA", "diagchas_YNNA", "diagchro_YNNA", "diagcopd_YNNA", "diagdiab_YNNA", "diagemph_YNNA", "diagfibr_YNNA", 
       "diaghear_YNNA", "diaghype_YNNA", "diagpneu_YNNA", "diagsarc_YNNA", "diagsili_YNNA", "diagstro_YNNA", "diagtube_YNNA", # medical hist
       "cancblad_YNNA", "cancbrea_YNNA", "canccerv_YNNA", "canccolo_YNNA", "cancesop_YNNA", "canckidn_YNNA", "canclary_YNNA", "canclung_YNNA", 
       "cancnasa_YNNA", "cancoral_YNNA", "cancpanc_YNNA", "cancphar_YNNA", "cancstom_YNNA", "cancthyr_YNNA", "canctran_YNNA")]<-
  apply(dat[,c("diagadas", "diagasbe", "diagbron", "diagchas", "diagchro", "diagcopd", "diagdiab", "diagemph", "diagfibr", 
                                                        "diaghear", "diaghype", "diagpneu", "diagsarc", "diagsili", "diagstro", "diagtube", # medical hist
                                                        "cancblad", "cancbrea", "canccerv", "canccolo", "cancesop", "canckidn", "canclary", "canclung", 
                                                        "cancnasa", "cancoral", "cancpanc", "cancphar", "cancstom", "cancthyr", "canctran")], 2, 
                                                 function(x){ifelse(x == "Missing", NA, as.character(x))})
dat$diagcopd_num<-ifelse(dat$diagcopd_YNNA=="Yes",1,
                         ifelse(dat$diagcopd_YNNA=="No",0,NA))
dat$diagemph_num<-ifelse(dat$diagemph_YNNA=="Yes",1,
                         ifelse(dat$diagemph_YNNA=="No",0,NA))

table(dat$diagcopd_num, dat$diagcopd_YNNA, useNA = "ifany")
table(dat$diagemph_num, dat$diagemph_YNNA, useNA = "ifany")


# Number of confirmed cancer
dat$num_confirmed_condensed<-ifelse(dat$num_confirmed>=2, 2, dat$num_confirmed)
table(dat$num_confirmed_condensed, dat$num_confirmed)

# Alcohol per day
dat$acrin_drinknum_curr_num<-as.numeric(as.character(dat$acrin_drinknum_curr))
dat$acrin_drinknum_form_num<-as.numeric(as.character(dat$acrin_drinknum_form))

dat$acfin_num_alc_perday<-ifelse(is.na(dat$acrin_drinknum_curr_num)==0 & is.na(dat$acrin_drinknum_form_num)==0, NA,
                           ifelse(is.na(dat$acrin_drinknum_curr_num)==1 | is.na(dat$acrin_drinknum_form_num)==1,
                                  ifelse(is.na(dat$acrin_drinknum_curr_num)==0, dat$acrin_drinknum_curr_num/7,
                                         ifelse(is.na(dat$acrin_drinknum_form_num)==0, dat$acrin_drinknum_form_num/7, NA)),NA))


# Married
table(dat$marital, useNA = "ifany")
dat$marital_condensed<-ifelse(dat$marital %in% c("Missing", "Not Ascertained", "Participant refused to answer"), NA, 
                              as.character(dat$marital))
table(dat$marital_condensed, dat$marital, useNA = "ifany") # 154 NA

# Education
table(dat$educat, dat$rndgroup, useNA = "ifany")
dat$educat_condensed<-ifelse(substr(dat$educat, 1, 7) %in% c("Missing", "Unknown"), NA, 
                             as.character(dat$educat))
table(dat$educat_condensed, dat$educat, useNA = "ifany")

dat1<-dat[,c("cen", "pid", "rndgroup", "study", # Study
             "age", "educat_condensed",  "ethnic", "gender", "height", "marital_condensed" , "race", "weight", "race_condensed", "ethnic_condensed","BMI", # Demo
             "age_quit_num", "cigar_YNNA", "cigsmok", "pipe_YNNA", "pkyr", "smokeage_num", "smokeday", 
             "smokelive_YNNA", "smokework_YNNA", "smokeyr", "yrs_from_cessation",# Smk
             "scr_days0", 
             "fup_days_adj_allcausemort", "allcause_death", 
             "allcause_death_8yr", "allcause_death_8yr_time", 
             "fup_days_adj_lungmort", "lung_death", 
             "lung_death_8yr", "lung_death_8yr_time", # Outcome
             "diagadas_YNNA", "diagasbe_YNNA", "diagbron_YNNA", "diagchas_YNNA", "diagchro_YNNA", "diagcopd_YNNA", "diagdiab_YNNA", "diagemph_YNNA", "diagfibr_YNNA", 
             "diaghear_YNNA", "diaghype_YNNA", "diagpneu_YNNA", "diagsarc_YNNA", "diagsili_YNNA", "diagstro_YNNA", "diagtube_YNNA", # medical hist
             "cancblad_YNNA", "cancbrea_YNNA", "canccerv_YNNA", "canccolo_YNNA", "cancesop_YNNA", "canckidn_YNNA", "canclary_YNNA", "canclung_YNNA", 
             "cancnasa_YNNA", "cancoral_YNNA", "cancpanc_YNNA", "cancphar_YNNA", "cancstom_YNNA", "cancthyr_YNNA", "canctran_YNNA", # cancer hist
             "fam_hist_lung", # family history
              "wrk_hist", #work hist,
             "race_condensed2_num","female","diagcopd_num","diagemph_num","fam_hist_lung_num","CT"
             )] # 67 var
write.csv(dat1, "NSLT_baseline_072320.csv")

