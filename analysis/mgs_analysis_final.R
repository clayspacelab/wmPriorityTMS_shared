rm(list=ls())

library(dplyr)
library(tidyr)
library(permuco)

SCRIPTS_DIR = '/datd/TMS_Priority/analysis_Nathan'
DATA_FILE = '../data/TMS_priority_MGS_data_11-Oct-2024_allconds.csv'

setwd(SCRIPTS_DIR)

datadf_allconds <- read.table(DATA_FILE,header=TRUE,sep=',')

subj_valid = c(1,2,3,4,5,6,7,11,12,16,17,21,22,24)
subj_excl <- c(8,9,20)

#double check we have the right subjects
stopifnot("Wrong subjects detected!"=sort(unique(datadf_allconds$subject))==sort(c(subj_valid,subj_excl)))

#compute averages for ANOVA
datadf_ave_allconds <- summarise(group_by(datadf_allconds,subject,hemi,condition,priority),
                         mean_i_sacc_err=mean(i_sacc_err),
                        mean_i_sacc_rt=mean(i_sacc_rt))

#contrasts constant across analyses
datadf_ave_allconds$priority.f <- as.factor(datadf_ave_allconds$priority)
contrasts(datadf_ave_allconds$priority.f) <- contr.sum
contrasts(datadf_ave_allconds$priority.f)

datadf_ave_allconds$hemi.f <- as.factor(datadf_ave_allconds$hemi)
contrasts(datadf_ave_allconds$hemi.f) <- contr.sum
contrasts(datadf_ave_allconds$hemi.f)

head(datadf_ave_allconds)


#### ANALYSIS ####

### NoTMS vs. sPCS

#primary analyses focus on PCS vs. noTMS, so remove IPS
datadf_ave_all <- subset(datadf_ave_allconds,condition != 'l_ips2')
stopifnot("Wrong conditions detected"=unique(datadf_ave_all$condition)==c("l_spcs","noTMS"))

#set factors
datadf_ave_all$condition.f <- as.factor(datadf_ave_all$condition)
contrasts(datadf_ave_all$condition.f) <- contr.sum
contrasts(datadf_ave_all$condition.f)

head(datadf_ave_all)

#subset based on exclusion criteria
datadf_ave <- subset(datadf_ave_all,!(subject %in% subj_excl))
#unique(datadf_ave$subject)
stopifnot("Wrong subjects detected"=sort(unique(datadf_ave$subject))==subj_valid)


#let's separate by hemi for follow-up analyses
datadfR_ave <- subset(datadf_ave,hemi=='RVF')
datadfL_ave <- subset(datadf_ave,hemi=='LVF')

head(datadfR_ave)
head(datadfL_ave)


## Sac Error


#Right Hemi
aov.mgs.perm.R <- aovperm(mean_i_sacc_err ~ condition.f * priority.f + 
                          Error(subject/(condition.f*priority.f)),
                        data = datadfR_ave, np = 10000)

summary(aov.mgs.perm.R)

#Left Hemi
aov.mgs.perm.L <- aovperm(mean_i_sacc_err ~ condition.f * priority.f + 
                            Error(subject/(condition.f*priority.f)),
                          data = datadfL_ave, np = 10000)

summary(aov.mgs.perm.L)

#3-way ANOVA
aov.mgs.perm <- aovperm(mean_i_sacc_err ~ condition.f * priority.f * hemi.f + 
                          Error(subject/(condition.f*priority.f*hemi.f)),
                        data = datadf_ave, np = 10000)

summary(aov.mgs.perm)


## Sac RT

#Right Hemi
aov.mgs.perm.R.rt <- aovperm(mean_i_sacc_rt ~ condition.f * priority.f + 
                               Error(subject/(condition.f*priority.f)),
                             data = datadfR_ave, np = 10000)

summary(aov.mgs.perm.R.rt)

#Left Hemi
aov.mgs.perm.L.rt <- aovperm(mean_i_sacc_rt ~ condition.f * priority.f + 
                               Error(subject/(condition.f*priority.f)),
                             data = datadfL_ave, np = 10000)

summary(aov.mgs.perm.L.rt)


### IPS2 analyses (RVF only)

#subset to RVF and based on exclusion criteria
datadfR_ips_ave <- subset(datadf_ave_allconds,hemi == 'RVF' & !(subject %in% subj_excl))
stopifnot("Wrong hemi detected"=unique(datadfR_ips_ave$hemi)=="RVF")
stopifnot("Wrong subjects detected"=sort(unique(datadfR_ips_ave$subject))==subj_valid)

#we just want to check for difference in priority effect for ips vs. notms 
#set factors
datadfR_ips_noTMS_ave <- subset(datadfR_ips_ave,condition!='l_spcs')
stopifnot("Wrong conds detected"=unique(datadfR_ips_noTMS_ave$condition)==c('l_ips2','noTMS'))

#set condition factors
datadfR_ips_noTMS_ave$condition.f <- as.factor(datadfR_ips_noTMS_ave$condition)
contrasts(datadfR_ips_noTMS_ave$condition.f) <- contr.sum
contrasts(datadfR_ips_noTMS_ave$condition.f)
head(datadfR_ips_noTMS_ave)


## Sac error

#ips vs notms
aov.mgs.perm.ips_notms <- aovperm(mean_i_sacc_err ~ condition.f * priority.f + 
                            Error(subject/(condition.f*priority.f)),
                          data = datadfR_ips_noTMS_ave, np = 10000)

summary(aov.mgs.perm.ips_notms)


## Sac RT

#ips vs notms
aov.mgs.perm.ips_notms.rt <- aovperm(mean_i_sacc_rt ~ condition.f * priority.f + 
                                    Error(subject/(condition.f*priority.f)),
                                  data = datadfR_ips_noTMS_ave, np = 10000)

summary(aov.mgs.perm.ips_notms.rt)

