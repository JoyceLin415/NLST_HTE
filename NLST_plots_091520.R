rm(list=ls())

library(here)
library(survival)
library(tableone)

options(scipen = 10000)

dat0<-read.csv(here("NSLT_baseline_072320.csv"))
dat0$current_smk<-ifelse(dat0$cigsmok=="Current",1,
                         ifelse(dat0$cigsmok=="Former",0,NA))
table(dat0$current_smk, useNA = "ifany")

dat<-dat0[,c("pid", "CT",
             "fup_days_adj_allcausemort", "allcause_death",
             "fup_days_adj_lungmort", "lung_death",  
             "age",  "female", "race_condensed2_num", "pkyr", "BMI",
             "diagcopd_num", "diagemph_num", "fam_hist_lung_num", "current_smk")]
dat$trt<-ifelse(dat$CT==1, "CT", "X-ray")

dat$fup_days_adj_allcausemort_new<-ifelse(dat$fup_days_adj_allcausemort==0,0.1,dat$fup_days_adj_allcausemort)
dat$fup_days_adj_lungmort_new<-ifelse(dat$fup_days_adj_lungmort==0,0.1,dat$fup_days_adj_lungmort)

dat_noNA<-na.omit(dat) # N = 50062
table(dat_noNA$CT, useNA = "ifany") # 25204 CT, 24858 X-ray

#--- LC Mort and OS KM (not included in manuscript) ----
osplot<-function(osdays, os, group, title, datax, ylab, xlab, surv_legend, maxtm, line_col, show_p, 
                 legend_x, legend_y, show_numatrisk){
  if(sum(is.na(osdays), is.na(os), is.na(group)) > 0){
    stop("Use data without NA!!!!!!!")
  }

  groupname<-substr(names(survfit(Surv(osdays, os)~group)$strata), 7,
                    nchar(names(survfit(Surv(osdays, os)~group)$strata)))
  n_grp<-length(groupname)
  
  os0 <- summary(survfit(Surv(osdays, os)~group))
  data0 <- cbind(os0$time, os0$surv, rep(1:n_grp, table(os0$strata)))
  data0 <- rbind(data0, cbind(rep(0, n_grp), rep(1, n_grp), 1:n_grp))  
  data0 <- rbind(data0, cbind(rep(max(maxtm, max(osdays, na.rm = T)), n_grp), tapply(data0[,2], data0[,3], min), 1:n_grp))
  data0 <- data0[order(data0[,3], data0[,1]),]
  tt<-survdiff((Surv(osdays, os)~group), rho=0)
  pval<-1-pchisq(tt$chisq, length(tt$n)-1)
  pval<-ifelse(pval<0.001, paste("<0.001"), paste(round(pval,3)))
  
  par(mar = c(7,4.5,2,0.5))
  plot(c(0, max(osdays, na.rm = T)), c(0,1), pch=" ", ylab=ylab, xlab=xlab,font.lab=2,
       cex.lab=1.4, mgp=c(2.5,1,0), xaxt="n", yaxt="n", bty="n", main=title, xlim=c(0,maxtm), ylim=c(0.7,1))
  axis(1, at=seq(0, maxtm, 12*30), labels=seq(0, (maxtm)/30, 12), pos=0.701, cex.axis=1.4, font=2, lwd.ticks=3, lwd=3)
  axis(2, at=0:20*0.05, labels=0:20*5, pos=0.001, cex.axis=1.4, font=2, lwd.ticks=3, lwd=3, las=2)
  
  for(i in 1:n_grp){
    #lines(data0[data0[,3]==i,1], data0[data0[,3]==i,2], lwd=4, type="s", col=line_col[i], lty = line_lty[i])
    lines(data0[data0[,3]==i,1], data0[data0[,3]==i,2], lwd=2, type="s", col=line_col[i])
    text(legend_x*maxtm, legend_y-0.02*i, groupname[i], col = line_col[i], cex = 1.4)
    text((legend_x+0.25)*maxtm, legend_y-0.02*i, paste(round(min(data0[data0[,3]==i & data0[,1]<=maxtm,][,2]), 3)), cex = 1.4)
  }
  
  text((legend_x+0.25)*maxtm, legend_y, paste(surv_legend), cex = 1.4)

  
  if(show_p == T){
    text((legend_x+0.68)*maxtm, legend_y, paste("log-rank test p-value"), cex = 1.4)
    text((legend_x+0.68)*maxtm, legend_y-0.02, paste(pval), cex = 1.4)
  }
  
  if(show_numatrisk == T){
    mtext("Patients at risk", side = 1, at = 1, line = 3, cex = 1.3)
    
    for(i in 1:n_grp){
      mtext(groupname[i], side = 1, at = -340, line = 3+i*1.2, col = line_col[i], cex = 1.3)
    }
    for(i in seq(0,maxtm,12*30)){
      nar<-tapply(osdays, group, function(x){sum(x>=i, na.rm = T)})
      for(j in 1:length(nar)){
        mtext(nar[j], side = 1, at = i, line = (3+j*1.2), col = line_col[j], cex = 1.3)
      }
    }
  }
}

# Max FU
summary(dat$fup_days_adj_allcausemort_new)/360
summary(dat$fup_days_adj_lungmort_new)/360

table(dat$lung_death, dat$trt, useNA = "ifany")

jpeg(here("091520","Fig 1_1 LCS OS by trt.jpg"), height = 7, width = 14, units = "in", res = 500)
par(mfrow = c(1,2))
# LC Mort by trt
osplot(osdays = dat_noNA$fup_days_adj_lungmort_new, os = dat_noNA$lung_death, group = dat_noNA$trt, 
       title = "",
       datax = dat_noNA, ylab = "Lung Cancer Survival (%)", xlab = "Months from Randomization", 
       surv_legend = "8-year LC survival", maxtm = 360*8, 
       line_col = c("deepskyblue3","deeppink3"), show_p = T, legend_x = 0.1, legend_y = 0.78, show_numatrisk = T)

# OS by trt
os1<-survfit(Surv(time=dat_noNA$fup_days_adj_allcausemort, event = dat_noNA$allcause_death)~dat_noNA$trt)
os1_dat_noNA<-cbind(summary(os1)$strata, summary(os1)$surv, summary(os1)$time)
tapply(os1_dat_noNA[,2],os1_dat_noNA[,1],min)
tapply(os1_dat_noNA[,3],os1_dat_noNA[,1],max)/365
#plot(os1, ylim=c(0.8,1))

osplot(osdays = dat_noNA$fup_days_adj_allcausemort, os = dat_noNA$allcause_death, group = dat_noNA$trt, 
       title = "",
       datax = dat_noNA, ylab = "Overall Survival (%)", xlab = "Months from Randomization", 
       surv_legend = "8-year OS", maxtm = 360*8, 
       line_col = c("deepskyblue3","deeppink3"), show_p = T, legend_x = 0.1, legend_y = 0.78, show_numatrisk = T)

mtext(text = "A. Lung Cancer Survival", side = 3, line = -1.3, at = 0.11, outer = T, font=2, cex=1.4)
mtext(text = "B. Overall Survival", side = 3, line = -1.3, at = 0.59, outer = T, font=2, cex=1.4)
dev.off()

# Median FU
survfit(Surv(time=dat_noNA$fup_days_adj_allcausemort/360, event=(1-dat_noNA$allcause_death))~1)

# Freq of death
table(dat_noNA$allcause_death, dat_noNA$trt, useNA = "ifany")
23429+1775; 1775/(23429+1775)
23004+1854; 1854/(23004+1854)

table(dat_noNA$lung_death, dat_noNA$trt, useNA = "ifany")
24872+332; 332/(24872+332)
24453+405; 405/(24453+405)

#--- Table 1 ----
table1<-CreateTableOne(vars=c("age",  "female", "race_condensed2_num", "pkyr", "BMI",
                              "diagcopd_num", "diagemph_num", "fam_hist_lung_num", "current_smk"),
                       strata = "trt", data=dat_noNA, 
                       factorVars = c("female", "race_condensed2_num", "diagcopd_num", "diagemph_num", 
                                      "fam_hist_lung_num", "current_smk"), test=T, addOverall = T)
write.csv(print(table1, nonnormal = c("age", "pkyr", "BMI"), minMax = T, showAllLevels = T),
          here("091520", "Table 1.csv"))

hist(dat_noNA$pkyr)
hist(dat_noNA$BMI)
hist(dat_noNA$age)

#--- Fig 2 ----
# LC ITE (diff of log time)
tmp_lc<-read.csv(here("080420","survdat_yhat_MCMC_8yrLC9cov_080420.csv")) # Skip 300, keep 1000
tmp_lc[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))]<-
  (tmp_lc[,c(paste("yhat_CT_MCMC_",seq(1,1000),sep=""))]-tmp_lc[,c(paste("yhat_Xray_MCMC_",seq(1,1000),sep=""))])
tmp_lc$difflogT_mean_overMCMC<-apply(tmp_lc[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],1,mean)
summary(tmp_lc$difflogT_mean_overMCMC)
tmp_lc$difflogT_025_overMCMC<-apply(tmp_lc[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],1,function(x){quantile(x,0.025)})
tmp_lc$difflogT_975_overMCMC<-apply(tmp_lc[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],1,function(x){quantile(x,0.975)})

## Order data by ITE
tmp_lc_ITE<-tmp_lc[order(tmp_lc$difflogT_mean_overMCMC),]
tmp_lc_ITE$id_new<-c(1:nrow(tmp_lc_ITE))

# OS ITE (diff of log time)
tmp_os<-read.csv(here("080420","survdat_yhat_MCMC_8yrOS9cov_080320.csv")) # Skip 300, keep 1000
tmp_os[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))]<-
  (tmp_os[,c(paste("yhat_CT_MCMC_",seq(1,1000),sep=""))]-tmp_os[,c(paste("yhat_Xray_MCMC_",seq(1,1000),sep=""))])
tmp_os$difflogT_mean_overMCMC<-apply(tmp_os[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],1,mean)
summary(tmp_os$difflogT_mean_overMCMC)
tmp_os$difflogT_025_overMCMC<-apply(tmp_os[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],1,function(x){quantile(x,0.025)})
tmp_os$difflogT_975_overMCMC<-apply(tmp_os[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],1,function(x){quantile(x,0.975)})

## Order data by ITE
tmp_os_ITE<-tmp_os[order(tmp_os$difflogT_mean_overMCMC),]
tmp_os_ITE$id_new<-c(1:nrow(tmp_os_ITE))

# Plot ITE
jpeg(here("091520","Fig 2_1 Ratio of T (individual) from AFTrees 8yr LC OS 9 cov.jpg"), height = 7, width = 14, 
     units = "in", res = 500)
par(mar = c(5,5,4.1,0.8), mgp = c(3.7,1,0), mfrow=c(1,2))
## LC
plot(exp(tmp_lc_ITE$difflogT_mean_overMCMC), tmp_lc_ITE$id_new, yaxt="n", type = "n", xlim = c(0.6,1.5),
     main = "", ylab = "Patient index", 
     xlab = paste("Ratio of failure time (CT/CXR)", sep = ""), 
     cex.axis = 1.4, cex.main = 1.4, cex.lab = 1.4, font.lab=2, font=2)
axis(side = 2, at = seq(0,50000,5000), lab = c(paste(seq(0,50,5),"K",sep = "")), 
     las = 1, cex.axis = 1.4, cex.lab = 1.4, font=2)
for(i in 1:nrow(tmp_lc_ITE)){
  lines(x = c(exp(tmp_lc_ITE$difflogT_025_overMCMC[i]), exp(tmp_lc_ITE$difflogT_975_overMCMC[i])), 
        y = rep(tmp_lc_ITE$id_new[i],2), lty = 1, col = "grey", lwd = 0.5)
  #cat(paste("i=",i,sep=""))
}
abline(v = 1, col = "red", lty = 2, lwd = 1.5)
points(exp(tmp_lc_ITE$difflogT_mean_overMCMC), tmp_lc_ITE$id_new, pch = 20, cex = 0.6)

## OS
plot(exp(tmp_os_ITE$difflogT_mean_overMCMC), tmp_os_ITE$id_new, yaxt="n", type = "n", xlim = c(0.6,1.5),
     main = "", ylab = "Patient index", 
     xlab = paste("Ratio of failure time (CT/CXR)", sep = ""), 
     cex.axis = 1.4, cex.main = 1.4, cex.lab = 1.4, font.lab=2, font=2)
axis(side = 2, at = seq(0,50000,5000), lab = c(paste(seq(0,50,5),"K",sep = "")), 
     las = 1, cex.axis = 1.4, cex.lab = 1.4, font=2)
for(i in 1:nrow(tmp_os_ITE)){
  lines(x = c(exp(tmp_os_ITE$difflogT_025_overMCMC[i]), exp(tmp_os_ITE$difflogT_975_overMCMC[i])), 
        y = rep(tmp_os_ITE$id_new[i],2), lty = 1, col = "grey", lwd = 0.5)
  #cat(paste("i=",i,sep=""))
}
abline(v = 1, col = "red", lty = 2, lwd = 1.5)
points(exp(tmp_os_ITE$difflogT_mean_overMCMC), tmp_os_ITE$id_new, pch = 20, cex = 0.6)

mtext(text = "A. Lung Cancer Survival", side = 3, line = -2.0, at = 0.11, outer = T, font=2, cex=1.4)
mtext(text = "B. Overall Survival", side = 3, line = -2.0, at = 0.58, outer = T, font=2, cex=1.4)

dev.off()


#--- difflogT HTE (subgroup) LC (not in manuscipt) ----
# Create 300 subgroups based on quantile of difflogT
n_subgrp<-300 
for(i in 1:n_subgrp){
  cut0<-quantile(tmp_lc$difflogT_mean_overMCMC,1/n_subgrp*(i-1))
  cut1<-quantile(tmp_lc$difflogT_mean_overMCMC,1/n_subgrp*(i))
  
  if(i == 1){
    tmp_lc$difflogT_grp<-ifelse(tmp_lc$difflogT_mean_overMCMC<=cut1,i,0)
  }else if(i >1 & i <= n_subgrp){
    tmp_lc$difflogT_grp<-ifelse(tmp_lc$difflogT_mean_overMCMC<=cut1 & tmp_lc$difflogT_mean_overMCMC>cut0,i,tmp_lc$difflogT_grp)
  }else {
    tmp_lc$difflogT_grp<-NA
  }
}
table(tmp_lc$difflogT_grp, tmp_lc$CT)
length(table(tmp_lc$difflogT_grp)) # 300 out of 300
tmp_lc$difflogT_grp0<-tmp_lc$difflogT_grp

## If all observations in subgroup i all received CT/Xray, change their subgroup to 999
difflogT_grp_0_ID<-rownames(table(tmp_lc$difflogT_grp0, tmp_lc$CT))[
  table(tmp_lc$difflogT_grp0, tmp_lc$CT)[,1]==0 | table(tmp_lc$difflogT_grp0, tmp_lc$CT)[,2]==0]

tmp_lc$difflogT_grp<-ifelse(tmp_lc$difflogT_grp0%in%as.numeric(difflogT_grp_0_ID), 999, tmp_lc$difflogT_grp0)
length(table(tmp_lc$difflogT_grp)) # 280 out of 300

#write.csv(table(tmp_lc$difflogT_grp0, tmp_lc$CT), here("080420/plot","num of obs in each difflogT subgroup_9cov_8yrLC.csv"))

# difflogT Plot
## Remove group 999
tmp_lc2<-tmp_lc[tmp_lc$difflogT_grp!=999,]
difflogT_dat<-data.frame("grp_id" = c(1:length(table(tmp_lc2$difflogT_grp))) ,
                         "difflogT_mean_in_subgrp" = tapply(tmp_lc2$difflogT_mean_overMCMC,tmp_lc2$difflogT_grp,mean))
difflogT_dat[,paste("difflogT_mean_in_subgrp_MCMC_",seq(1,1000),sep="")] <- apply(tmp_lc2[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,
                                                                                  function(x){tapply(x,tmp_lc2$difflogT_grp,mean)})
difflogT_dat[,c("difflogT_025_in_subgrp_MCMC_")] <- apply(difflogT_dat[,c(paste("difflogT_mean_in_subgrp_MCMC_",seq(1,1000),sep=""))],1,quantile,probs = 0.025)
difflogT_dat[,c("difflogT_975_in_subgrp_MCMC_")] <- apply(difflogT_dat[,c(paste("difflogT_mean_in_subgrp_MCMC_",seq(1,1000),sep=""))],1,quantile,probs = 0.975)

min(difflogT_dat[,c("difflogT_025_in_subgrp_MCMC_")]) # -0.07
max(difflogT_dat[,c("difflogT_975_in_subgrp_MCMC_")]) # 0.09

jpeg(here("080420/plot","Ratio of T subgrp results from BART 8yr LC 9 cov.jpg"), 
     height = 8, width = 6, units = "in", res = 2500)
difflogT_dat<-difflogT_dat[order(difflogT_dat$difflogT_mean_in_subgrp),]
difflogT_dat$grp_id2<-seq(1,nrow(difflogT_dat),1)
par(mar = c(5,5,4.1,2.1), mgp = c(3.7,1,0))
plot(exp(difflogT_dat$difflogT_mean_in_subgrp), difflogT_dat$grp_id2, yaxt="n", type = "n", cex = 0.6, xlim = c(0.9,1.1),
     main = "", ylab = "Subgroup index", 
     xlab = paste("Ratio of T (CT/X-ray)\n",nrow(difflogT_dat), " groups", sep = ""), 
     cex.axis = 1.3, cex.main = 1.4, cex.lab = 1.3)
axis(side = 2, at = seq(0,300,50), las = 1, cex.axis = 1.3, cex.lab = 1.3)
for(i in 1:nrow(difflogT_dat)){
  lines(x = c(exp(difflogT_dat[i,c("difflogT_025_in_subgrp_MCMC_")]), 
              exp(difflogT_dat[i,c("difflogT_975_in_subgrp_MCMC_")]) ), 
        y = rep(i,2), lty = 1, col = "grey", lwd = 0.5)
}
abline(v = 1, col = "red", lty = 2, lwd = 1.5)
points(exp(difflogT_dat$difflogT_mean_in_subgrp), difflogT_dat$grp_id2, pch = 20, cex = 0.6)
dev.off()

#--- difflogT HTE (subgroup) OS (not in manuscipt) ----
# Create subgroups
n_subgrp<-300 
for(i in 1:n_subgrp){
  cut0<-quantile(tmp_os$difflogT_mean_overMCMC,1/n_subgrp*(i-1))
  cut1<-quantile(tmp_os$difflogT_mean_overMCMC,1/n_subgrp*(i))
  
  if(i == 1){
    tmp_os$difflogT_grp<-ifelse(tmp_os$difflogT_mean_overMCMC<=cut1,i,0)
  }else if(i >1 & i <= n_subgrp){
    tmp_os$difflogT_grp<-ifelse(tmp_os$difflogT_mean_overMCMC<=cut1 & tmp_os$difflogT_mean_overMCMC>cut0,i,tmp_os$difflogT_grp)
  }else {
    tmp_os$difflogT_grp<-NA
  }
}
table(tmp_os$difflogT_grp, tmp_os$CT)
length(table(tmp_os$difflogT_grp)) # 300 out of 300
tmp_os$difflogT_grp0<-tmp_os$difflogT_grp
difflogT_grp_0_ID<-rownames(table(tmp_os$difflogT_grp0, tmp_os$CT))[
  table(tmp_os$difflogT_grp0, tmp_os$CT)[,1]==0 | table(tmp_os$difflogT_grp0, tmp_os$CT)[,2]==0]

tmp_os$difflogT_grp<-ifelse(tmp_os$difflogT_grp0%in%as.numeric(difflogT_grp_0_ID), 999, tmp_os$difflogT_grp0)
length(table(tmp_os$difflogT_grp)) # 297 out of 300
tapply(tmp_os$difflogT_mean_overMCMC, tmp_os$difflogT_grp, summary)
#write.csv(table(tmp_os$difflogT_grp0, tmp_os$CT), here("080420/plot","num of obs in each difflogT subgroup_9cov_8yrOS.csv"))

# difflogT Plot
tmp_os2<-tmp_os[tmp_os$difflogT_grp!=999,]
difflogT_dat<-data.frame("grp_id" = c(1:length(table(tmp_os2$difflogT_grp))) ,
                         "difflogT_mean_in_subgrp" = tapply(tmp_os2$difflogT_mean_overMCMC,tmp_os2$difflogT_grp,mean))
difflogT_dat[,paste("difflogT_mean_in_subgrp_MCMC_",seq(1,1000),sep="")] <- apply(tmp_os2[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,
                                                                                  function(x){tapply(x,tmp_os2$difflogT_grp,mean)})
difflogT_dat[,c("difflogT_025_in_subgrp_MCMC_")] <- apply(difflogT_dat[,c(paste("difflogT_mean_in_subgrp_MCMC_",seq(1,1000),sep=""))],1,quantile,probs = 0.025)
difflogT_dat[,c("difflogT_975_in_subgrp_MCMC_")] <- apply(difflogT_dat[,c(paste("difflogT_mean_in_subgrp_MCMC_",seq(1,1000),sep=""))],1,quantile,probs = 0.975)

min(difflogT_dat[,c("difflogT_025_in_subgrp_MCMC_")]) # -0.25
max(difflogT_dat[,c("difflogT_975_in_subgrp_MCMC_")]) # 0.32

jpeg(here("080420/plot","Ratio of T subgrp results from BART 8yr OS 9 cov.jpg"), 
     height = 8, width = 6, units = "in", res = 2500)
difflogT_dat<-difflogT_dat[order(difflogT_dat$difflogT_mean_in_subgrp),]
difflogT_dat$grp_id2<-seq(1,nrow(difflogT_dat),1)
par(mar = c(5,5,4.1,2.1), mgp = c(3.7,1,0))
plot(exp(difflogT_dat$difflogT_mean_in_subgrp), difflogT_dat$grp_id2, yaxt="n", type = "n", cex = 0.6, xlim = c(0.7,1.4),
     main = "", ylab = "Subgroup index", 
     xlab = paste("Ratio of T (CT/X-ray)\n",nrow(difflogT_dat), " groups", sep = ""), 
     cex.axis = 1.3, cex.main = 1.4, cex.lab = 1.3)
axis(side = 2, at = seq(0,300,50), las = 1, cex.axis = 1.3, cex.lab = 1.3)
for(i in 1:nrow(difflogT_dat)){
  lines(x = c(exp(difflogT_dat[i,c("difflogT_025_in_subgrp_MCMC_")]), 
              exp(difflogT_dat[i,c("difflogT_975_in_subgrp_MCMC_")]) ), 
        y = rep(i,2), lty = 1, col = "grey", lwd = 0.5)
}
abline(v = 1, col = "red", lty = 2, lwd = 1.5)
points(exp(difflogT_dat$difflogT_mean_in_subgrp), difflogT_dat$grp_id2, pch = 20, cex = 0.6)
dev.off()

#--- Compare difflogT between groups LC (not in manuscipt) ----
mean_CI<-function(id){
  n_Y<-length(id); n_N<-50062-length(id)
  p_Y<-round(n_Y/50062, digits = 2); p_N<-round(n_N/50062, digits = 2)
  
  node_dat_Y<-tmp_lc[tmp_lc$pid%in%id==1,]
  (u_Y<-round(exp(mean(node_dat_Y$difflogT_mean_overMCMC)), digits = 4))
  (L95_Y<-round(exp(quantile(apply(node_dat_Y[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,mean), 0.025)), digits = 4))
  (U95_Y<-round(exp(quantile(apply(node_dat_Y[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,mean), 0.975)), digits = 4))
  
  node_dat_N<-tmp_lc[tmp_lc$pid%in%id==0,]
  (u_N<-round(exp(mean(node_dat_N$difflogT_mean_overMCMC)), digits = 4))
  (L95_N<-round(exp(quantile(apply(node_dat_N[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,mean), 0.025)), digits = 4))
  (U95_N<-round(exp(quantile(apply(node_dat_N[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,mean), 0.975)), digits = 4))
  
  mat<-matrix(c(n_Y,p_Y,u_Y,L95_Y,U95_Y,
                n_N,p_N,u_N,L95_N,U95_N), nrow = 1)
  colnames(mat)<-c("n_Y","p_Y","u_Y","L95_Y","U95_Y",
                   "n_N","p_N","u_N","L95_N","U95_N")
  return(mat)
  
}

# Race
tapply(tmp_lc$difflogT_mean_overMCMC, tmp_lc$race_condensed2_num, summary)

## White+Other 
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(1,4)]
(white<-mean_CI(node_id))

## Black+Asian
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(2,3)]
(black_asian<-mean_CI(node_id))

# Age
summary(tmp_lc$age)

## Find cut from 1 pct to 99 pct by 5
quantile(tmp_lc$age,0.25) # 57
quantile(tmp_lc$age,0.75) # 65

a<-(65-56)
u_diff<-matrix(NA,nrow=1,ncol=a)
for(i in 1:a){
  cutpt<-56+1*i
  id2<-tmp_lc$pid[tmp_lc$age>cutpt]
  u_diff[,i]<-abs(mean(tmp_lc$difflogT_mean_overMCMC[tmp_lc$pid%in%id2==1])-mean(tmp_lc$difflogT_mean_overMCMC[tmp_lc$pid%in%id2==0]))
}
colnames(u_diff)<-seq(57,65,1)
max(u_diff) # cut = 65
u_diff

##
node_id<-tmp_lc$pid[tmp_lc$age>65]
(age65<-mean_CI(node_id))

# Female
node_id<-tmp_lc$pid[tmp_lc$female==1]
(female<-mean_CI(node_id))

# Pkyr
summary(tmp_lc$pkyr)

## Find cut from 1 pct to 99 pct by 5
quantile(tmp_lc$pkyr,0.01) # 30
quantile(tmp_lc$pkyr,0.99) # 138

a<-(140-25)/5
u_diff<-matrix(NA,nrow=1,ncol=a)
for(i in 1:a){
  cutpt<-25+5*i
  id2<-tmp_lc$pid[tmp_lc$pkyr>cutpt]
  u_diff[,i]<-abs(mean(tmp_lc$difflogT_mean_overMCMC[tmp_lc$pid%in%id2==1])-mean(tmp_lc$difflogT_mean_overMCMC[tmp_lc$pid%in%id2==0]))
}
colnames(u_diff)<-seq(30,140,5)
max(u_diff) # cut = 30
u_diff

## 
node_id<-tmp_lc$pid[tmp_lc$pkyr>30]
(pkyr30<-mean_CI(node_id))

# COPD
node_id<-tmp_lc$pid[tmp_lc$diagcopd_num==1]
(copd<-mean_CI(node_id))

# Emphy
node_id<-tmp_lc$pid[tmp_lc$diagemph_num==1]
(emph<-mean_CI(node_id))

# Fam Hist
node_id<-tmp_lc$pid[tmp_lc$fam_hist_lung_num==1]
(famhist<-mean_CI(node_id))

# BMI
summary(tmp_lc$BMI)

## Find cut from 1 pct to 99 pct by 5
quantile(tmp_lc$BMI,0.25) # 24
quantile(tmp_lc$BMI,0.75) # 31

a<-(31-23)/1
u_diff<-matrix(NA,nrow=1,ncol=a)
for(i in 1:a){
  cutpt<-23+1*i
  id2<-tmp_lc$pid[tmp_lc$BMI>cutpt]
  u_diff[,i]<-abs(mean(tmp_lc$difflogT_mean_overMCMC[tmp_lc$pid%in%id2==1])-mean(tmp_lc$difflogT_mean_overMCMC[tmp_lc$pid%in%id2==0]))
}
colnames(u_diff)<-seq(24,31,1)
max(u_diff) # cut = 31
u_diff

## 
node_id<-tmp_lc$pid[tmp_lc$BMI>31]
(BMI31<-mean_CI(node_id))

# Current/former smoker
node_id<-tmp_lc$pid[tmp_lc$current_smk==1]
(cur_smk<-mean_CI(node_id))

# Combine
node_diff<-rbind(white, black_asian, age65, female,
                 pkyr30, copd, emph, famhist, BMI31, cur_smk)
rownames(node_diff)<-c("White or Other", "Black or Asian", "Age > 65", "Female",
                       "Pack-year > 30", "COPD", "Emphysema", "Family history", "BMI > 31", 
                       "Current smoker")
node_diff
write.csv(node_diff, here("080420/plot", "node_diff_LC_AFT_9cov.csv"))

#--- Compare difflogT between groups OS (not in manuscipt) ----
mean_CI<-function(id){
  n_Y<-length(id); n_N<-50062-length(id)
  p_Y<-round(n_Y/50062, digits = 2); p_N<-round(n_N/50062, digits = 2)
  
  node_dat_Y<-tmp_os[tmp_os$pid%in%id==1,]
  (u_Y<-round(exp(mean(node_dat_Y$difflogT_mean_overMCMC)), digits = 4))
  (L95_Y<-round(exp(quantile(apply(node_dat_Y[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,mean), 0.025)), digits = 4))
  (U95_Y<-round(exp(quantile(apply(node_dat_Y[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,mean), 0.975)), digits = 4))
  
  node_dat_N<-tmp_os[tmp_os$pid%in%id==0,]
  (u_N<-round(exp(mean(node_dat_N$difflogT_mean_overMCMC)), digits = 4))
  (L95_N<-round(exp(quantile(apply(node_dat_N[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,mean), 0.025)), digits = 4))
  (U95_N<-round(exp(quantile(apply(node_dat_N[,c(paste("difflogT_MCMC_",seq(1,1000),sep=""))],2,mean), 0.975)), digits = 4))
  
  mat<-matrix(c(n_Y,p_Y,u_Y,L95_Y,U95_Y,
                n_N,p_N,u_N,L95_N,U95_N), nrow = 1)
  colnames(mat)<-c("n_Y","p_Y","u_Y","L95_Y","U95_Y",
                   "n_N","p_N","u_N","L95_N","U95_N")
  return(mat)
}

# Race
tapply(tmp_os$difflogT_mean_overMCMC, tmp_os$race_condensed2_num, summary)

## White+Other 
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1,4)]
(white<-mean_CI(node_id))

## Black+Asian
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(2,3)]
(black_asian<-mean_CI(node_id))

# Age
summary(tmp_os$age)

## Find cut from 1 pct to 99 pct by 5
quantile(tmp_os$age,0.25) # 57
quantile(tmp_os$age,0.75) # 65

a<-(65-56)
u_diff<-matrix(NA,nrow=1,ncol=a)
for(i in 1:a){
  cutpt<-56+1*i
  id2<-tmp_os$pid[tmp_os$age>cutpt]
  u_diff[,i]<-abs(mean(tmp_os$difflogT_mean_overMCMC[tmp_os$pid%in%id2==1])-mean(tmp_os$difflogT_mean_overMCMC[tmp_os$pid%in%id2==0]))
}
colnames(u_diff)<-seq(57,65,1)
max(u_diff) # cut = 57
u_diff

##
node_id<-tmp_os$pid[tmp_os$age>57]
(age57<-mean_CI(node_id))

# Female
node_id<-tmp_os$pid[tmp_os$female==1]
(female<-mean_CI(node_id))

# Pkyr
summary(tmp_os$pkyr)

## Find cut from 1 pct to 99 pct by 5
quantile(tmp_os$pkyr,0.01) # 30
quantile(tmp_os$pkyr,0.99) # 138

a<-(140-25)/5
u_diff<-matrix(NA,nrow=1,ncol=a)
for(i in 1:a){
  cutpt<-25+5*i
  id2<-tmp_os$pid[tmp_os$pkyr>cutpt]
  u_diff[,i]<-abs(mean(tmp_os$difflogT_mean_overMCMC[tmp_os$pid%in%id2==1])-mean(tmp_os$difflogT_mean_overMCMC[tmp_os$pid%in%id2==0]))
}
colnames(u_diff)<-seq(30,140,5)
max(u_diff) # cut = 35
u_diff

## 
node_id<-tmp_os$pid[tmp_os$pkyr>35]
(pkyr35<-mean_CI(node_id))

# COPD
node_id<-tmp_os$pid[tmp_os$diagcopd_num==1]
(copd<-mean_CI(node_id))

# Emphy
node_id<-tmp_os$pid[tmp_os$diagemph_num==1]
(emph<-mean_CI(node_id))

# Fam Hist
node_id<-tmp_os$pid[tmp_os$fam_hist_lung_num==1]
(famhist<-mean_CI(node_id))

# BMI
summary(tmp_os$BMI)

## Find cut from 1 pct to 99 pct by 5
quantile(tmp_os$BMI,0.25) # 24
quantile(tmp_os$BMI,0.75) # 31

a<-(31-23)/1
u_diff<-matrix(NA,nrow=1,ncol=a)
for(i in 1:a){
  cutpt<-23+1*i
  id2<-tmp_os$pid[tmp_os$BMI>cutpt]
  u_diff[,i]<-abs(mean(tmp_os$difflogT_mean_overMCMC[tmp_os$pid%in%id2==1])-mean(tmp_os$difflogT_mean_overMCMC[tmp_os$pid%in%id2==0]))
}
colnames(u_diff)<-seq(24,31,1)
max(u_diff) # cut = 30
u_diff

## 
node_id<-tmp_os$pid[tmp_os$BMI>30]
(BMI30<-mean_CI(node_id))

# Current/former smoker
node_id<-tmp_os$pid[tmp_os$current_smk==1]
(cur_smk<-mean_CI(node_id))

# Combine
node_diff<-rbind(white, black_asian, age57, female,
                 pkyr35, copd, emph, famhist, BMI30,cur_smk )
rownames(node_diff)<-c("White or Other", "Black or Asian", "Age > 57", "Female",
                       "Pack-year > 35", "COPD", "Emphysema", "Family history", "BMI > 30", 
                       "Current smoker")
node_diff
write.csv(node_diff, here("080420/plot", "node_diff_OS_9cov_AFT_080420.csv"))



#--- Fit the fit LC ----
tmp_lc_tree<-tmp_lc[,c("difflogT_mean_overMCMC", "age", "race_condensed2_num", "pkyr", "BMI",
                       "female", "diagcopd_num","diagemph_num","fam_hist_lung_num", "current_smk")]

tmp_lc_tree$race_condensed2_num<-factor(tmp_lc_tree$race_condensed2_num)
tmp_lc_tree$female<-factor(tmp_lc_tree$female)
tmp_lc_tree$diagcopd_num<-factor(tmp_lc_tree$diagcopd_num)
tmp_lc_tree$diagemph_num<-factor(tmp_lc_tree$diagemph_num)
tmp_lc_tree$fam_hist_lung_num<-factor(tmp_lc_tree$fam_hist_lung_num)
tmp_lc_tree$current_smk<-factor(tmp_lc_tree$current_smk)

# Fit the Fit
## Step 1: Uni CART
var_list0<-c("age", "race_condensed2_num", "pkyr", "BMI",
             "female", "diagcopd_num","diagemph_num","fam_hist_lung_num", "current_smk")
R2<-0
names(R2)<-"Empty"

for(i in 1:length(var_list0)){
  var_i<-var_list0[i]
  cart_uni<-rpart(tmp_lc_tree$difflogT_mean_overMCMC~tmp_lc_tree[,var_i], data = tmp_lc_tree, method = "anova")
  
  R2<-c(R2, 1-min(cart_uni$cptable[,c("rel error")]))
  
  names(R2)[i+1]<-var_i
  cat(i)
}
(var_uni<-names(R2)[R2==max(R2)]) # race
R2_uni<-R2
(maxR2_0<-max(R2_uni)) # 0.82

# Step 2: Add a var to cart with max R2
d<-1
var_list<-var_list0[-which(var_list0==var_uni)]
var_old<-var_uni
R2<-0
maxR2_1<-maxR2_0

while(d>=0.01){
  for(i in 1:length(var_list)){
    var_i<-var_list[i]
    var_dat<-tmp_lc_tree[,c("difflogT_mean_overMCMC", var_old, var_i)]
    cart_2<-rpart(difflogT_mean_overMCMC~., data = var_dat, method = "anova")
    R2<-c(R2, 1-min(cart_2$cptable[,c("rel error")]))
    names(R2)[i+1]<-var_i
  }
  
  (var_max<-names(R2)[R2==max(R2)] )
  (maxR2_2<-max(R2))
  (d<-maxR2_2-maxR2_1)
  
  if(d>=0.01){
    maxR2_1<-maxR2_2
    var_old<-c(var_old, var_max)
    var_list<-var_list[-which(var_list==var_max)]
    R2<-0
  }
  
}
var_old
cart1<-rpart(difflogT_mean_overMCMC~race_condensed2_num, data = tmp_lc_tree, method = "anova")
cart2<-rpart(difflogT_mean_overMCMC~race_condensed2_num+diagcopd_num, data = tmp_lc_tree, method = "anova")
cart3<-rpart(difflogT_mean_overMCMC~race_condensed2_num+diagcopd_num+age, data = tmp_lc_tree, method = "anova")
cart4<-rpart(difflogT_mean_overMCMC~race_condensed2_num+diagcopd_num+age+diagemph_num, data = tmp_lc_tree, method = "anova")

max(1-cart1$cptable[,c("rel error")]) # 0.82
max(1-cart2$cptable[,c("rel error")]) # 0.94
max(1-cart3$cptable[,c("rel error")]) # 0.97
max(1-cart4$cptable[,c("rel error")]) # 0.98

pdf(here("080420/plot", "LC_BART_9cov_fit the fit_080420.pdf"))
prp(cart4)
dev.off()

# White
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(1)]
(values<-mean_CI(node_id))

# White & No COPD
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(1) & tmp_lc$diagcopd_num==0]
(values<-mean_CI(node_id))

# White & No COPD & Age>=72
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(1) & tmp_lc$diagcopd_num==0 & tmp_lc$age>=72]
(values<-mean_CI(node_id))

# White & No COPD & Age<72
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(1) & tmp_lc$diagcopd_num==0 & tmp_lc$age<72]
(values<-mean_CI(node_id))

# White & No COPD & Age<72 & No emph
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(1) & tmp_lc$diagcopd_num==0 & tmp_lc$age<72 & tmp_lc$diagemph_num==0]
(values<-mean_CI(node_id))

# White & No COPD & Age<72 & Emph
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(1) & tmp_lc$diagcopd_num==0 & tmp_lc$age<72 & tmp_lc$diagemph_num==1]
(values<-mean_CI(node_id))

# White & COPD
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(1) & tmp_lc$diagcopd_num==1]
(values<-mean_CI(node_id))

# Black+Asian+Other
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(2,3,4)]
(values<-mean_CI(node_id))

# Asian
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(3)]
(values<-mean_CI(node_id))

# Black+Other
node_id<-tmp_lc$pid[tmp_lc$race_condensed2_num%in%c(2,4)]
(values<-mean_CI(node_id))

#--- Fit the fit OS ----
tmp_os_tree<-tmp_os[,c("difflogT_mean_overMCMC", "age", "race_condensed2_num", "pkyr", "BMI",
                       "female", "diagcopd_num","diagemph_num","fam_hist_lung_num", "current_smk")]

tmp_os_tree$race_condensed2_num<-factor(tmp_os_tree$race_condensed2_num)
tmp_os_tree$female<-factor(tmp_os_tree$female)
tmp_os_tree$diagcopd_num<-factor(tmp_os_tree$diagcopd_num)
tmp_os_tree$diagemph_num<-factor(tmp_os_tree$diagemph_num)
tmp_os_tree$fam_hist_lung_num<-factor(tmp_os_tree$fam_hist_lung_num)
tmp_os_tree$current_smk<-factor(tmp_os_tree$current_smk)

# Fit the Fit
## Step 1: Uni CART
var_list0<-c("age", "race_condensed2_num", "pkyr", "BMI",
             "female", "diagcopd_num","diagemph_num","fam_hist_lung_num","current_smk")
R2<-0
names(R2)<-"Empty"

for(i in 1:length(var_list0)){
  var_i<-var_list0[i]
  cart_uni<-rpart(tmp_os_tree$difflogT_mean_overMCMC~tmp_os_tree[,var_i], data = tmp_os_tree, method = "anova")
  
  R2<-c(R2, 1-min(cart_uni$cptable[,c("rel error")]))
  
  names(R2)[i+1]<-var_i
  cat(i)
}
(var_uni<-names(R2)[R2==max(R2)]) # race
R2_uni<-R2
(maxR2_0<-max(R2_uni)) # 0.74

# Step 2: Add a var to cart with max R2
d<-1
var_list<-var_list0[-which(var_list0==var_uni)]
var_old<-var_uni
R2<-0
maxR2_1<-maxR2_0

while(d>=0.01){
  for(i in 1:length(var_list)){
    var_i<-var_list[i]
    var_dat<-tmp_os_tree[,c("difflogT_mean_overMCMC", var_old, var_i)]
    cart_2<-rpart(difflogT_mean_overMCMC~., data = var_dat, method = "anova")
    R2<-c(R2, 1-min(cart_2$cptable[,c("rel error")]))
    names(R2)[i+1]<-var_i
  }
  
  (var_max<-names(R2)[R2==max(R2)] )
  (maxR2_2<-max(R2))
  (d<-maxR2_2-maxR2_1)
  
  if(d>=0.01){
    maxR2_1<-maxR2_2
    var_old<-c(var_old, var_max)
    var_list<-var_list[-which(var_list==var_max)]
    R2<-0
  }
  
}
var_old
cart1<-rpart(difflogT_mean_overMCMC~race_condensed2_num, data = tmp_os_tree, method = "anova")
cart2<-rpart(difflogT_mean_overMCMC~race_condensed2_num+diagemph_num, data = tmp_os_tree, method = "anova")
cart3<-rpart(difflogT_mean_overMCMC~race_condensed2_num+diagemph_num+pkyr, data = tmp_os_tree, method = "anova")

max(1-cart1$cptable[,c("rel error")]) # 0.74
max(1-cart2$cptable[,c("rel error")]) # 0.85
max(1-cart3$cptable[,c("rel error")]) # 0.96


pdf(here("080420/plot", "OS_BART_9cov_fit the fit_080420.pdf"))
prp(cart3)
dev.off()


# White+Black+Other
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1,2,4)]
(values<-mean_CI(node_id))

# Other
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(4)]
(values<-mean_CI(node_id))

# White+Black
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1,2)]
(values<-mean_CI(node_id))

# White+Black & Emph
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1,2) & tmp_os$diagemph_num==1]
(values<-mean_CI(node_id))

# White+Black & No Emph
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1,2) & tmp_os$diagemph_num==0]
(values<-mean_CI(node_id))

# White+Black & No Emph & pkyr<37
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1,2) & tmp_os$diagemph_num==0 & tmp_os$pkyr<37]
(values<-mean_CI(node_id))

# White & No Emph & pkyr<37
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1) & tmp_os$diagemph_num==0 & tmp_os$pkyr<37]
(values<-mean_CI(node_id))

# Black & No Emph & pkyr<37
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(2) & tmp_os$diagemph_num==0 & tmp_os$pkyr<37]
(values<-mean_CI(node_id))

# White+Black & No Emph & pkyr>=37
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1,2) & tmp_os$diagemph_num==0 & tmp_os$pkyr>=37]
(values<-mean_CI(node_id))

# White & No Emph & pkyr>=37
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(1) & tmp_os$diagemph_num==0 & tmp_os$pkyr>=37]
(values<-mean_CI(node_id))

# Black & No Emph & pkyr>=37
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(2) & tmp_os$diagemph_num==0 & tmp_os$pkyr>=37]
(values<-mean_CI(node_id))

# Asian
node_id<-tmp_os$pid[tmp_os$race_condensed2_num%in%c(3)]
(values<-mean_CI(node_id))


