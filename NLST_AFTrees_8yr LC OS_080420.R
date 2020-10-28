# 9 cov; calculate That

rm(list = ls())
library(here)
library(BART)
library(AFTrees)

dat0<-read.csv("/Users/linjungyi/Documents/mount sinai/LH PCORI/NLST/NSLT_baseline_072320.csv")
dat0$current_smk<-ifelse(dat0$cigsmok=="Current",1,
                         ifelse(dat0$cigsmok=="Former",0,NA))
table(dat0$current_smk, useNA = "ifany")

dat<-dat0[,c("pid", "CT",
             "fup_days_adj_allcausemort", "allcause_death",
             "fup_days_adj_lungmort", "lung_death",  
             "age",  "female", "race_condensed2_num", "pkyr", "BMI",
             "diagcopd_num", "diagemph_num", "fam_hist_lung_num", "current_smk")]
dat_noNA<-na.omit(dat) # N = 50062
table(dat_noNA$CT, useNA = "ifany") # 25204 CT, 24858 X-ray

dat_noNA$fup_days_adj_allcausemort_new<-ifelse(dat_noNA$fup_days_adj_allcausemort==0,0.1,dat_noNA$fup_days_adj_allcausemort)
dat_noNA$fup_days_adj_lungmort_new<-ifelse(dat_noNA$fup_days_adj_lungmort==0,0.1,dat_noNA$fup_days_adj_lungmort)

#---- Surv BART 7 cov ----#
tmp<-dat_noNA

# Categorical variables were labeled as 0,1,... so need to change them to factor
tmp$female<-factor(tmp$female)
tmp$race_condensed2_num<-factor(tmp$race_condensed2_num)
tmp$fam_hist_lung_num<-factor(tmp$fam_hist_lung_num)
tmp$diagcopd_num<-factor(tmp$diagcopd_num)
tmp$diagemph_num<-factor(tmp$diagemph_num)
tmp$current_smk<-factor(tmp$current_smk)

# To calculate ITE, assume each person receive both CT and Xray
## Train (trt as original received one)
train<-tmp[,c("CT","age","female","race_condensed2_num","pkyr","diagcopd_num","diagemph_num",
              "fam_hist_lung_num", "BMI", "current_smk")]

## Test data (reverse the trt in original data)
test_CT2X<-tmp[,c("CT","age","female","race_condensed2_num","pkyr","diagcopd_num","diagemph_num",
                  "fam_hist_lung_num", "BMI", "current_smk")] 
test_CT2X$CT<- 1-test_CT2X$CT # Change CT to X-ray and X-ray to CT for ITE

train$CT<-factor(train$CT)
test_CT2X$CT<-factor(test_CT2X$CT)

# LC AFTrees
## Modeling
Sys.time()
set.seed(1234)
bart_LC_8yr<-AFTrees(x.train = train, 
                     y.train = tmp$fup_days_adj_lungmort_new, status = tmp$lung_death,  
                     x.test = test_CT2X, ndpost=1000, nskip = 300)
Sys.time() 
saveRDS(bart_LC_8yr, "/Users/linjungyi/Documents/mount sinai/LH PCORI/NLST/bart_LC_8yr_7covBMICursmk_300_1000.rds")

#bart_LC_8yr<-readRDS("/Users/linjungyi/Documents/mount sinai/LH PCORI/NLST/bart_LC_8yr_7covBMICursmk_300_1000.rds")
#str(bart_LC_8yr)
#head(bart_LC_8yr$m.train[,c(1:100)])

## Calculate That for CT and Xray (ITE is calculated in XXX_plots.r)
AFT_model<-bart_LC_8yr
num_post<-nrow(AFT_model$locations) # number of posterior

tmp[,c(paste("yhat_train_MCMC_",seq(1,num_post),sep=""))]<-t(AFT_model$m.train)
tmp[,c(paste("yhat_test_MCMC_",seq(1,num_post),sep=""))]<-t(AFT_model$m.test)

### If the obs received CT, That for CT will be in train data; if received Xray, That for CT will be in test data
for(i in 1:num_post){
  tmp[,c(paste("yhat_CT_MCMC_",i,sep=""))]<-ifelse(tmp$CT==1, 
                                                   tmp[,c(paste("yhat_train_MCMC_",i,sep=""))],
                                                   tmp[,c(paste("yhat_test_MCMC_",i,sep=""))])
  tmp[,c(paste("yhat_Xray_MCMC_",i,sep=""))]<-ifelse(tmp$CT==0, 
                                                     tmp[,c(paste("yhat_train_MCMC_",i,sep=""))],
                                                     tmp[,c(paste("yhat_test_MCMC_",i,sep=""))])
}

tmp<-tmp[,-which(names(tmp)%in%c(paste("yhat_train_MCMC_",seq(1,num_post),sep=""),
                                 paste("yhat_test_MCMC_",seq(1,num_post),sep="")))]
write.csv(tmp,"/Users/linjungyi/Documents/mount sinai/LH PCORI/NLST/080420/survdat_yhat_MCMC_8yrLC9cov_080420.csv",
          row.names = F)

# OS
## Modeling
Sys.time()
set.seed(1234)
bart_OS_8yr<-AFTrees(x.train = train, 
                     y.train = tmp$fup_days_adj_allcausemort_new, status = tmp$allcause_death,  
                     x.test = test_CT2X, ndpost=1000, nskip = 300)
Sys.time() 
saveRDS(bart_OS_8yr, "/Users/linjungyi/Documents/mount sinai/LH PCORI/NLST/bart_OS_8yr_7covBMICursmk_300_1000.rds")

#bart_OS_8yr<-readRDS("/Users/linjungyi/Documents/mount sinai/LH PCORI/NLST/bart_OS_8yr_7covBMICursmk_300_1000.rds")
#str(bart_OS_8yr)
#head(bart_OS_8yr$m.train[,c(1:100)])

## Calculate That for CT and Xray (ITE is calculated in XXX_plots.r)
AFT_model<-bart_OS_8yr
num_post<-nrow(AFT_model$locations)

tmp[,c(paste("yhat_train_MCMC_",seq(1,num_post),sep=""))]<-t(AFT_model$m.train)
tmp[,c(paste("yhat_test_MCMC_",seq(1,num_post),sep=""))]<-t(AFT_model$m.test)
for(i in 1:num_post){
  tmp[,c(paste("yhat_CT_MCMC_",i,sep=""))]<-ifelse(tmp$CT==1, 
                                                   tmp[,c(paste("yhat_train_MCMC_",i,sep=""))],
                                                   tmp[,c(paste("yhat_test_MCMC_",i,sep=""))])
  tmp[,c(paste("yhat_Xray_MCMC_",i,sep=""))]<-ifelse(tmp$CT==0, 
                                                     tmp[,c(paste("yhat_train_MCMC_",i,sep=""))],
                                                     tmp[,c(paste("yhat_test_MCMC_",i,sep=""))])
}

tmp<-tmp[,-which(names(tmp)%in%c(paste("yhat_train_MCMC_",seq(1,num_post),sep=""),
                                 paste("yhat_test_MCMC_",seq(1,num_post),sep="")))]
write.csv(tmp,"/Users/linjungyi/Documents/mount sinai/LH PCORI/NLST/080420/survdat_yhat_MCMC_8yrOS9cov_080320.csv",
          row.names = F)







