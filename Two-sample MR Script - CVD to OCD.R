#Analysis by Robyn Wootton
#Begun April 2024
#R version 4.4.0 (2024-04-24) -- "Puppy Cup"
#TwoSampleMR v0.6.1
################################################################################
## Contents
# 1. Load required packages 
# 2. Read in Exposures 
# 3. Read in outcome data 
# 4. Harmonising Exposure and Outcome Data
# 5. The MR Analysis
# 6. MR Egger Sensitivity Checks
# 7. Test for reverse causation
# 8. Visual Inspection
# 9. Plot the Results

################################################################################
#########               1. Load required packages                   ############
################################################################################
rm(list=ls())

#install.packages("devtools")
library(devtools)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)      
#install.packages("ggplot2")
library(ggplot2)
#install.packages("knitr")
library(knitr)
#install_github("MRCIEU/MRInstruments")
library(MRInstruments)
#install.packages("ieugwasr")
library(ieugwasr)
#install.packages("dplyr")
library(dplyr)

################################################################################
#########                  2. Read in Exposures                     ############
################################################################################
#### First extract the traits from the IEU repository
#The data structure for the two-sample MR package is explained here:
# https://mrcieu.github.io/TwoSampleMR/articles/index.html

#We now need a token to be able to extract from IEU GWAS repository
#https://api.opengwas.io 

#Binary traits
#Coronary Artery Disease ebi-a-GCST003116
#MI ieu-a-798 - MR BASE says mixed, not available in european only?
#stroke - ebi-a-GCST006906
#T2D - ebi-a-GCST005047

#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

#Extract binary instruments in the IEU repository
exposure_dat_bin <- TwoSampleMR::extract_instruments(
  c("ebi-a-GCST003116", "ieu-a-798", "ebi-a-GCST006906", "ebi-a-GCST005047", "finn-b-I9_MI"),
  p1 = 5e-08,
  clump = TRUE,
  p2 = 5e-08,
  r2 = 0.001,
  kb = 10000)

#get N snp
n <- subset(exposure_dat_bin, exposure_dat_bin$id.exposure=="ebi-a-GCST003116")
#CAD = 42
n <- subset(exposure_dat_bin, exposure_dat_bin$id.exposure=="ebi-a-GCST006906")
#stroke = 8
n <- subset(exposure_dat_bin, exposure_dat_bin$id.exposure=="ebi-a-GCST005047")
#T2D = 33
n <- subset(exposure_dat_bin, exposure_dat_bin$id.exposure=="ieu-a-798")
#MI mixed ancestry = 26
n <- subset(exposure_dat_bin, exposure_dat_bin$id.exposure=="finn-b-I9_MI")
#MI finn genn = 13

#Continuous traits in the IEU repository
#Total Cholesterol - ebi-a-GCST002221
#HDL Cholesterol - ebi-a-GCST002223
#LDL Cholesterol - ebi-a-GCST002222
#Triglycerides - ebi-a-GCST002216
#BMI - ebi-a-GCST002783

#Extract continuous instruments
exposure_dat_cont <- TwoSampleMR::extract_instruments(
  c("ebi-a-GCST002223", "ebi-a-GCST002222", "ebi-a-GCST002221", "ebi-a-GCST002783", "ebi-a-GCST002216"),
  p1 = 5e-08,
  clump = TRUE,
  p2 = 5e-08,
  r2 = 0.001,
  kb = 10000)

#get N snp
n <- subset(exposure_dat_cont, exposure_dat_cont$id.exposure=="ebi-a-GCST002221")
#total chol = 87
n <- subset(exposure_dat_cont, exposure_dat_cont$id.exposure=="ebi-a-GCST002223")
#HDL = 87
n <- subset(exposure_dat_cont, exposure_dat_cont$id.exposure=="ebi-a-GCST002222")
#LDL = 75
n <- subset(exposure_dat_cont, exposure_dat_cont$id.exposure=="ebi-a-GCST002216")
#trig = 54
n <- subset(exposure_dat_cont, exposure_dat_cont$id.exposure=="ebi-a-GCST002783")
#bmi = 79

#### Second, extract the traits not available in the IEU repository

##Heart failure
#Perform clumping
hf_exp <- read.table("/pathname/heartfailure_sumstats.txt", header=T)
head(hf_exp)
#name the columns as per ld_clump
hf_exp$rsid <- hf_exp$SNP
hf_exp$pval <- hf_exp$p
#Subset to instrument
hf_exp <- subset(hf_exp, hf_exp$pval<5e-8)
hf_exp <- ld_clump(hf_exp, clump_kb = 1000, clump_r2 = 0.01, opengwas_jwt = XXXX)
#Save the instrument
write.csv(hf_exp, "/pathname/HeartFailure_instrument.csv", quote=F, row.names=F )

#Read in using TwoSampleMR Package
hf_exp_dat <- TwoSampleMR::read_exposure_data(
  filename = "/pathname/HeartFailure_instrument.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p")
hf_exp_dat$exposure <- "Heart Failure"


##heart rate variability
#Perform clumping
exp <- read.table("/pathname/Nolte_28613276_SDNN.txt", header=T)
head(exp)
#name the columns as per ld_clump
exp$rsid <- exp$rs_number
exp$pval <- exp$p.value
#Subset to instrument
exp <- subset(exp, exp$pval<5e-8)
exp <- ld_clump(exp, clump_kb = 1000, clump_r2 = 0.01, opengwas_jwt = XXXX)
#Save the instrument
write.csv(exp, "/pathname/HRV_instrument.csv", quote=F, row.names=F )

#Read in using TwoSampleMR Package
hrv_exp_dat <- TwoSampleMR::read_exposure_data(
  filename = "/pathname/HRV_instrument.csv",
  sep = ",",
  snp_col = "rs_number",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "reference_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p.value")
hrv_exp_dat$exposure <- "Heart Rate Variability"


##Systolic Blood pressure 
#Perform clumping
exp <- read.csv("/pathname/SBP_summary.csv", header=T)
head(exp)
#name the columns as per ld_clump
exp$rsid <- exp$V3
exp$pval <- exp$meta.pval
#Subset to instrument
exp <- subset(exp, exp$pval<5e-8)
exp <- ld_clump(exp, clump_kb = 1000, clump_r2 = 0.01, opengwas_jwt = XXXX)
#Save the instrument
write.csv(exp, "/pathname/SBP_instrument.csv", quote=F, row.names=F )

#Read in using TwoSampleMR Package
sbp_exp_dat <- TwoSampleMR::read_exposure_data(
  sep=",",
  filename = "/pathname/SBP_instrument.csv",
  snp_col = "V3",
  beta_col = "meta.beta",
  se_col = "meta.se",
  effect_allele_col = "new_REF",
  other_allele_col = "new_ALT",
  eaf_col = "meta.w_eaf",
  pval_col = "meta.pval")
sbp_exp_dat$exposure <- "Systolic BP"

#Format effect sizes to be in SD units (to match the others)
#I have calculated the meta-analysed total sample SD to be 19.10
sbp_exp_dat$beta.exposure <- sbp_exp_dat$beta.exposure/19.1
sbp_exp_dat$se.exposure <- sbp_exp_dat$se.exposure/19.1

##Diastolic Blood pressure 
#Perform clumping
exp <- read.csv("/pathname/DBP_summary.csv", header=T)
head(exp)
#name the columns as per ld_clump
exp$rsid <- exp$V3
exp$pval <- exp$meta.pval
#Subset to instrument
exp <- subset(exp, exp$pval<5e-8)
exp <- ld_clump(exp, clump_kb = 1000, clump_r2 = 0.01, opengwas_jwt = XXXX)
#Save the instrument
write.csv(exp, "/pathname/DBP_instrument.csv", quote=F, row.names=F )

#Read in using TwoSampleMR Package
dbp_exp_dat <- TwoSampleMR::read_exposure_data(
  filename = "/pathname/DBP_instrument.csv",
  sep = ",",
  snp_col = "V3",
  beta_col = "meta.beta",
  se_col = "meta.se",
  effect_allele_col = "new_REF",
  other_allele_col = "new_ALT",
  eaf_col = "meta.w_eaf",
  pval_col = "meta.pval")
dbp_exp_dat$exposure <- "Diastolic BP"

#Format effect sizes to be in SD units
#I have calculated the meta-analysed total sample SD to be 10.86
dbp_exp_dat$beta.exposure <- dbp_exp_dat$beta.exposure/10.86
dbp_exp_dat$se.exposure <- dbp_exp_dat$se.exposure/10.86

##Thrombosis
#Perform clumping
exp <- read.table("/pathname/thrombosis.txt", header=T)
head(exp)
#name the columns as per ld_clump
exp$rsid <- exp$rsid
exp$pval <- exp$P.value
#Subset to instrument
exp <- subset(exp, exp$pval<5e-8)
exp <- ld_clump(exp, clump_kb = 1000, clump_r2 = 0.01, opengwas_jwt = XXXX)
#Save the instrument
write.csv(exp, "/pathname/thrombosis_instrument.csv", quote=F, row.names=F )

#Read in using TwoSampleMR Package
thromb_exp_dat <- TwoSampleMR::read_exposure_data(
  filename = "/pathname/thrombosis_instrument.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "Freq1",
  pval_col = "pval")
thromb_exp_dat$exposure <- "Thrombosis"

##################################################################################
###########                  3. Read in Outcome Data                   ###########  
##################################################################################

#Read in the OCD outcome data (excluding 23andMe)
ocd <- read.table("/pathname/daner_OCD_full_wo23andMe_190522", header=T)
head(ocd)
#make into log odds
ocd$beta <- log(ocd$OR)
#save
write.csv(ocd, "/pathname/ocd_outcomedata.csv", quote=F, row.names=F)

##Read in the outcome data for each exposure
#Binary exposures
out_dat_bin <- TwoSampleMR::read_outcome_data(
  snps=exposure_dat_bin$SNP,
  filename = "/pathname/ocd_outcomedata.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ_U_1114613",
  pval_col = "P")
out_dat_bin$outcome <- "OCD"

#add the outcome columns later for Steiger
out_dat_bin$samplesize.outcome<- 1138106
out_dat_bin$ncase.outcome<-23493
out_dat_bin$ncontrol.outcome<-1114613
out_dat_bin$units.outcome<-"log odds"
out_dat_bin$prevalence.outcome<-0.01

#continuous exposures
out_dat_cont <- TwoSampleMR::read_outcome_data(
  snps=exposure_dat_cont$SNP,
  filename = "/pathname/ocd_outcomedata.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ_U_1114613",
  pval_col = "P")
out_dat_cont$outcome <- "OCD"

#add the outcome columns later for steiger
out_dat_cont$samplesize.outcome<- 1138106
out_dat_cont$ncase.outcome<-23493
out_dat_cont$ncontrol.outcome<-1114613
out_dat_cont$units.outcome<-"log odds"
out_dat_cont$prevalence.outcome<-0.01

#heart failure
out_dat_hf <- TwoSampleMR::read_outcome_data(
  snps=hf_exp_dat$SNP,
  filename = "/pathname/ocd_outcomedata.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ_U_1114613",
  pval_col = "P")
out_dat_hf$outcome <- "OCD"

#add the columns later for steiger
out_dat_hf$samplesize.outcome<- 1138106
out_dat_hf$ncase.outcome<-23493
out_dat_hf$ncontrol.outcome<-1114613
out_dat_hf$units.outcome<-"log odds"
out_dat_hf$prevalence.outcome<-0.01

#heart rate variability
out_dat_hrv <- TwoSampleMR::read_outcome_data(
  snps=hrv_exp_dat$SNP,
  filename = "/pathname/ocd_outcomedata.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ_U_1114613",
  pval_col = "P")
out_dat_hrv$outcome <- "OCD"

#add the columns later for steiger
out_dat_hrv$samplesize.outcome<- 1138106
out_dat_hrv$ncase.outcome<-23493
out_dat_hrv$ncontrol.outcome<-1114613
out_dat_hrv$units.outcome<-"log odds"
out_dat_hrv$prevalence.outcome<-0.01


#SBP
out_dat_sbp <- TwoSampleMR::read_outcome_data(
  snps=sbp_exp_dat$SNP,
  filename = "/pathname/ocd_outcomedata.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ_U_1114613",
  pval_col = "P")
out_dat_sbp$outcome <- "OCD"

#add the columns later for steiger
out_dat_sbp$samplesize.outcome<- 1138106
out_dat_sbp$ncase.outcome<-23493
out_dat_sbp$ncontrol.outcome<-1114613
out_dat_sbp$units.outcome<-"log odds"
out_dat_sbp$prevalence.outcome<-0.01


#DBP
out_dat_dbp <- TwoSampleMR::read_outcome_data(
  snps=dbp_exp_dat$SNP,
  filename = "/pathname/ocd_outcomedata.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ_U_1114613",
  pval_col = "P")
out_dat_dbp$outcome <- "OCD"

#add the columns later for steiger
out_dat_dbp$samplesize.outcome<- 1138106
out_dat_dbp$ncase.outcome<-23493
out_dat_dbp$ncontrol.outcome<-1114613
out_dat_dbp$units.outcome<-"log odds"
out_dat_dbp$prevalence.outcome<-0.01

##Thrombosis
out_dat_thromb <- TwoSampleMR::read_outcome_data(
  snps=thromb_exp_dat$SNP,
  filename = "/pathname/ocd_outcomedata.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ_U_1114613",
  pval_col = "P")
out_dat_thromb$outcome <- "OCD"

#add the columns later for steiger
out_dat_thromb$samplesize.outcome<- 1138106
out_dat_thromb$ncase.outcome<-23493
out_dat_thromb$ncontrol.outcome<-1114613
out_dat_thromb$units.outcome<-"log odds"
out_dat_thromb$prevalence.outcome<-0.01

################################################################################
#########      4. Harmonising Exposure and Outcome Data             ############
################################################################################

#harmonise the data, IEU binary traits
data_bin <- TwoSampleMR::harmonise_data( 
  exposure_dat = exposure_dat_bin ,
  outcome_dat = out_dat_bin,
  action = 2)

#harmonise the data, IEU continuous traits
data_cont <- TwoSampleMR::harmonise_data( 
  exposure_dat = exposure_dat_cont ,
  outcome_dat = out_dat_cont,
  action = 2)

#harmonise the data, HF
data_hf <- TwoSampleMR::harmonise_data( 
  exposure_dat = hf_exp_dat ,
  outcome_dat = out_dat_hf,
  action = 2)

#harmonise the data, HRV
data_hrv <- TwoSampleMR::harmonise_data( 
  exposure_dat = hrv_exp_dat ,
  outcome_dat = out_dat_hrv,
  action = 2)

#harmonise the data, SBP
data_sbp <- TwoSampleMR::harmonise_data( 
  exposure_dat = sbp_exp_dat ,
  outcome_dat = out_dat_sbp,
  action = 2)

#harmonise the data, DBP
data_dbp <- TwoSampleMR::harmonise_data( 
  exposure_dat = dbp_exp_dat ,
  outcome_dat = out_dat_dbp,
  action = 2)

#harmonise the data, thrombosis
data_thromb <- TwoSampleMR::harmonise_data( 
  exposure_dat = thromb_exp_dat ,
  outcome_dat = out_dat_thromb,
  action = 2)

################################################################################
#########                   5. The MR Analysis                      ############
################################################################################

#Cholesterol contains comma so rename
data_cont$exposure[data_cont$exposure=="Cholesterol, total || id:ebi-a-GCST002221"] <- "Cholesterol total || id:ebi-a-GCST002221"

#Separate into each trait for the analyses
data_cad <- subset(data_bin, data_bin$exposure=="Coronary artery disease || id:ebi-a-GCST003116")
data_t2d <- subset(data_bin, data_bin$exposure=="Type 2 diabetes || id:ebi-a-GCST005047")
data_mi <- subset(data_bin, data_bin$exposure=="Myocardial infarction || id:ieu-a-798")
data_mi_finn <- subset(data_bin, data_bin$exposure=="Myocardial infarction || id:finn-b-I9_MI")
data_stroke <- subset(data_bin, data_bin$exposure=="Stroke || id:ebi-a-GCST006906")
data_bmi <- subset(data_cont, data_cont$exposure=="Body mass index || id:ebi-a-GCST002783")
data_ldl <- subset(data_cont, data_cont$exposure=="LDL cholesterol || id:ebi-a-GCST002222")
data_hdl <- subset(data_cont, data_cont$exposure=="HDL cholesterol || id:ebi-a-GCST002223")
data_chol <- subset(data_cont, data_cont$exposure=="Cholesterol total || id:ebi-a-GCST002221")
data_trig <- subset(data_cont, data_cont$exposure=="Triglycerides || id:ebi-a-GCST002216")

#Create lists of the dataframes (separately for binary and continuous)
all_bin <- list(data_cad, data_t2d, data_mi, data_stroke, data_hf, data_thromb, data_mi_finn)
all_cont <- list(data_bmi, data_ldl, data_hdl, data_chol, data_trig, data_hrv, data_sbp, data_dbp)

#loop over Continuous outcomes
output_het <- c()
output_res <- c()
output_egger <- c()

for(i in 1: length(all_cont)){
  name_list <- c("data_bmi", "data_ldl", "data_hdl", "data_chol", "data_trig", "data_hrv", "data_sbp", "data_dbp")
  name <- name_list[i]
  df <- data.frame(all_cont[i])
  mr_het <- TwoSampleMR::mr_heterogeneity(df)
  res <- TwoSampleMR::mr(df, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  or <- TwoSampleMR::generate_odds_ratios(res)
  mr_egger_int <- TwoSampleMR::mr_pleiotropy_test(df)
  mr_egger_int$lci <- mr_egger_int$egger_intercept-(1.96*mr_egger_int$se)
  mr_egger_int$uci <- mr_egger_int$egger_intercept+(1.96*mr_egger_int$se)
  output_het <-rbind(output_het, mr_het)
  output_res <-rbind(output_res, or)
  output_egger <-rbind(output_egger, mr_egger_int)
  write.csv(output_het, "/pathname/heterogeneity_stats_cont_rev.csv", quote=FALSE, row.names=FALSE)
  write.csv(output_res, "/pathname/MRresults_cont_rev.csv", quote=FALSE, row.names=FALSE)
  write.csv(output_egger, "/pathname/egger_cont_rev.csv", quote=FALSE, row.names=FALSE)
  
  ## Scatter Plot
  # Shows the distribution of the SNP-exposure effects against the SNP-outcome effects
  # Includes the main MR methods plotted
  p1 <- TwoSampleMR::mr_scatter_plot(res, df)
  filename_scatter <- paste0("/pathname/", name, "_scatter.png")
  ggsave(p1[[1]], file=filename_scatter, width=7, height=7)
}


#loop over Binary outcomes
output_het <- c()
output_res <- c()
output_egger <- c()

for(i in 1: length(all_bin)){
  name_list <- c("data_cad", "data_t2d", "data_mi", "data_stroke", "data_hf", "data_thromb", "data_mi_finn")
  name <- name_list[i]
  df <- data.frame(all_bin[i])
  mr_het <- mr_heterogeneity(df)
  res <- mr(df, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  res$b <- res$b*0.693 #following Stephen Burgess suggestion so estimates are per doubling in odds for exposure
  or <- generate_odds_ratios(res)
  mr_egger_int <- mr_pleiotropy_test(df)
  mr_egger_int$lci <- mr_egger_int$egger_intercept-(1.96*mr_egger_int$se)
  mr_egger_int$uci <- mr_egger_int$egger_intercept+(1.96*mr_egger_int$se)
  output_het <-rbind(output_het, mr_het)
  output_res <-rbind(output_res, or)
  output_egger <-rbind(output_egger, mr_egger_int)
  write.csv(output_het, "/pathname/heterogeneity_stats_bin_rev.csv", quote=FALSE, row.names=FALSE)
  write.csv(output_res, "/pathname/MRresults_bin_rev.csv", quote=FALSE, row.names=FALSE)
  write.csv(output_egger, "/pathname/egger_bin_rev.csv", quote=FALSE, row.names=FALSE)
  
  ## Scatter Plot
  # Shows the distribution of the SNP-exposure effects against the SNP-outcome effects
  # Includes the main MR methods plotted
  p1 <- mr_scatter_plot(res, df)
  filename_scatter <- paste0("/pathname/", name, "_scatter.png")
  ggsave(p1[[1]], file=filename_scatter, width=7, height=7)
}


################################################################################
#########               6. MR Egger Sensitivity Checks              ############
################################################################################

# I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

ll <- list( data_cad, data_t2d, data_mi, data_stroke, data_bmi,  data_ldl,  data_hdl, 
            data_chol,  data_trig, data_hf, data_cont, data_hrv, data_sbp, data_dbp, data_thromb, data_mi_finn)

output <- c()
for(i in 1: length(ll)){
  name_list <- c("data_cad", "data_t2d", "data_mi", "data_stroke", "data_bmi",  "data_ldl",  "data_hdl", 
                 "data_chol",  "data_trig", 'data_hf', "data_cont", "data_hrv", "data_sbp", "data_dbp", "data_thromb", "data_mi_finn")
  name <- name_list[i]
  df <- data.frame(ll[i])
  BetaXG   = df$beta.exposure
  seBetaXG = df$se.exposure 
  seBetaYG <- df$se.outcome
  BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
  
  # Calculate F statistics
  F   = BXG^2/seBetaXG^2
  count_10 = sum(F > 10)
  perc_F = (count_10/length(F))*100
  mF  = mean(F)
  # and I-squared statistics
  Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
  Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
  
  results <- cbind(name, mF, perc_F, Isq_unweighted, Isq_weighted)
  colnames(results) <- c("data", "mF", "perc_F", "Isq_unweighted", "Isq_weighted") #following the MRBase naming convention
  output <-rbind(output, results)
  colnames(output) <- c("data", "mF", "perc_F", "Isq_unweighted", "Isq_weighted")
  write.csv(output, "/pathname/regressiondilution_rev_revision.csv", quote=FALSE, row.names=FALSE)
} 

################################################################################
#########              7. Test for reverse causation                ############
################################################################################

## Steiger filtering
#continuous exposures/outcomes must contain: samplesize and units="SD"
#binary exposures/outcomes must contain: samplesize, ncase, ncontrol, prevalence and units="log odds"


#CAD
data_cad $samplesize.exposure<- 86995
data_cad $ncase.exposure<-22233
data_cad $ncontrol.exposure<- 64762
data_cad $units.exposure<-"log odds"
data_cad $prevalence.exposure<-0.05

#T2D
data_t2d$samplesize.exposure<- 69033
data_t2d$ncase.exposure <- 12171
data_t2d$ncontrol.exposure<- 	56862
data_t2d$units.exposure<-"log odds"
data_t2d$prevalence.exposure<-0.06
data_t2d $eaf.exposure<-data_t2d $eaf.outcome

#MI
data_mi $samplesize.exposure<- 171875
data_mi $ncase.exposure<-48371
data_mi $ncontrol.exposure<- 123504
data_mi $units.exposure<-"log odds"
data_mi $prevalence.exposure<-0.017

#Stroke
data_stroke $samplesize.exposure<- 521612
data_stroke $ncase.exposure<-67162
data_stroke $ncontrol.exposure<- 454450
data_stroke $units.exposure<-"log odds"
data_stroke $prevalence.exposure<-0.25

#Heart Failure
data_hf $samplesize.exposure<- 977323
data_hf $ncase.exposure<-47309
data_hf $ncontrol.exposure<- 930014
data_hf $units.exposure<-"log odds"
data_hf $prevalence.exposure<-0.023

#BMI
data_bmi $samplesize.exposure<- 339224
data_bmi $units.exposure<-"SD"

#LDL
data_ldl $samplesize.exposure<- 188578
data_ldl $units.exposure<-"SD"

#HDL
data_hdl $samplesize.exposure<- 188578
data_hdl $units.exposure<-"SD"

#Cholesterol
data_chol $samplesize.exposure<- 188578
data_chol $units.exposure<-"SD"

#Heart Rate Variability
data_hrv $samplesize.exposure<- 53174
data_hrv $units.exposure<-"SD"

#Triglicerides
data_trig $samplesize.exposure<- 188578
data_trig $units.exposure<-"SD"

#Systolic Blood Pressure
data_sbp $samplesize.exposure<- 150134
data_sbp $units.exposure<-"SD"

#Diastolic Blood Pressure
data_dbp $samplesize.exposure<- 150134
data_dbp $units.exposure<-"SD"

#Thrombosis
data_thromb $samplesize.exposure<- 995864
data_thromb $ncase.exposure<-62879
data_thromb $ncontrol.exposure<- 932985
data_thromb $units.exposure<-"log odds"
data_thromb $prevalence.exposure<-0.002

#Loop to run Steiger filtering for each exposure-outcome 
ll <- list(data_cad, data_t2d,  data_mi, data_stroke,  data_hf,   data_bmi,  data_ldl, 
           data_hdl, data_chol,  data_hrv,  data_trig,  data_sbp, data_dbp, data_thromb) 

output <- c()
for(i in 1: length(ll)){
  name_list <- c("data_cad", "data_t2d",  "data_mi", "data_stroke",  "data_hf",   "data_bmi",  "data_ldl", 
                 "data_hdl", "data_chol",  "data_hrv",  "data_trig",  "data_sbp", "data_dbp", "data_thromb")
  name <- name_list[i]
  df <- data.frame(ll[i])
  analysis_df <- subset(df, df$mr_keep==T)
  steiger <- TwoSampleMR::steiger_filtering(analysis_df)
  false <- length(steiger$steiger_dir[steiger$steiger_dir == FALSE])
  true <- length(steiger$steiger_dir[steiger$steiger_dir == TRUE])
  percent <- (true/(true+false))*100
  filter <- cbind(name, false, true, percent)
  colnames(filter) <- c("data", "falseN", "trueN", "percent_true") #following the MRBase naming convention
  output <-rbind(output, filter)
  colnames(output) <- c("data", "falseN", "trueN", "percent_true")
  write.csv(output, "/pathname/steiger_rev.csv", quote=FALSE, row.names=FALSE)
}

################################################################################
#########                   8. Visual Inspection                    ############
################################################################################

ll <- list( data_cad, data_t2d, data_mi, data_stroke, data_bmi,  data_ldl,  data_hdl, 
            data_chol,  data_trig, data_hf, data_cont, data_hrv, data_sbp, data_dbp, data_thromb)

for(i in 1: length(ll)){
  name_list <- c("data_cad", "data_t2d", "data_mi", "data_stroke", "data_bmi",  "data_ldl",  "data_hdl", 
                 "data_chol",  "data_trig", 'data_hf', "data_cont", "data_hrv", "data_sbp", "data_dbp", "data_thromb")
  name <- name_list[i]
  df <- data.frame(ll[i])
  
  ## Leave-one-out Analysis
  # Re-runs the IVW analysis, leaving one SNP out at a time
  # Check effects are not driven by a single SNP, estimate should stay relatively consistent 
  res_loo <- mr_leaveoneout(df)
  p2 <- mr_leaveoneout_plot(res_loo)
  loo_file <- paste0("/pathname/Visual Inspection/rev_", name, "loo.png")
  ggsave(p2[[1]], file=loo_file, width=7, height=10)
  
  ## Single SNP Analysis
  # Calculates a wald ratio for each individual SNP
  # Check for outlier SNPs, estimates should stay relatively consistent across SNPs
  res_single <- mr_singlesnp(df, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  p3 <- mr_forest_plot(res_single)
  forest_file <- paste0("/pathname/", name, "_single_forest.png")
  ggsave(p3[[1]], file=forest_file, width=7, height=10)
} 


################################################################################
#########                   9. Plot the Results                     ############
################################################################################
#install.packages("forestplot")
library(forestplot)

##Continuous outcomes
cont_res <- read.csv("/pathname/MRresults_cont_rev.csv", header=T)
head(cont_res)
cont_res$exposure <- gsub("\\|.*","",cont_res$exposure)
bin_res <- read.csv("/pathnames/MRresults_bin_rev.csv", header=T)
head(bin_res)
tail(bin_res)
#remove finn genn
table(bin_res$exposure)
bin_res <- subset(bin_res, bin_res$exposure != "Myocardial infarction || id:finn-b-I9_MI")
bin_res <- bin_res[-which(rownames(bin_res) == "14"), ] #check that this removes the MR Egger result for stroke because of low I2
bin_res$exposure <- gsub("\\|.*","",bin_res$exposure)
res <- rbind(cont_res, bin_res)

##Plot of just IVW results(combined binary and cont)
cont_res_ivw <- subset(cont_res, cont_res$method=="Inverse variance weighted")
bin_res_ivw <- subset(bin_res, bin_res$method=="Inverse variance weighted")
res_ivw <- rbind(cont_res_ivw, bin_res_ivw)
table(res_ivw$exposure)
res_ivw$exposure <- gsub("\\|.*","",res_ivw$exposure)

#Create a dataframe of the text
table_text <- cbind(
  c("Exposure", res_ivw$exposure),
  c("OR", sprintf("%.2f", res_ivw$or)),
  c("95% CI", paste(sprintf("%.2f", res_ivw$or_lci95), "-", sprintf("%.2f", res_ivw$or_uci95)))
)

#work out limits to set the axes
summary(res_ivw$or_lci95)
summary(res_ivw$or_uci95)

#Create the plot
pdf.options(reset = TRUE, onefile = FALSE)
pdf("cvd-ocd_IVW.pdf", width=6,height=7)
# Create the forest plot
forestplot(
  labeltext = table_text,
  mean = c(NA, res_ivw$or), # Add NA for the header row
  lower = c(NA, res_ivw$or_lci95),
  upper = c(NA, res_ivw$or_uci95),
  new_page = TRUE,
  zero=1,
  
  # Custom appearance
  is.summary = c(TRUE, rep(FALSE, 14)), # Summary at first and last rows
  xlab = "Odds Ratio",
  clip =c(0.5, 1.3), 
  xticks = seq(0.5, 1.4, by = 0.2),
  
  # Colors and styles
  col = fpColors(box = "darkblue", line = "darkblue", summary = "royalblue"),
  boxsize = 0.25, # Box size for individual studies
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 12),
    ticks = gpar(fontsize = 14),
    xlab = gpar(fontsize = 14, fontface = "italic"),
    title = gpar(fontsize = 14, fontface = "bold")
  ),
  
  # Grid lines
  grid = structure(c(1), gp = gpar(lty = 2, col = "gray"))
)
dev.off()

##Plot to Compare MI before and after
bin_res <- read.csv("/pathname/MRresults_bin_rev.csv", header=T)
head(bin_res)
#Subset to MI
table(bin_res$exposure)
bin_res <- subset(bin_res, bin_res$exposure== "Myocardial infarction || id:finn-b-I9_MI" | bin_res$exposure=="Myocardial infarction || id:ieu-a-798")
bin_res$GWAS[bin_res$exposure== "Myocardial infarction || id:finn-b-I9_MI"] <- "FinnGen"
bin_res$GWAS[bin_res$exposure== "Myocardial infarction || id:ieu-a-798"] <- "Nikpay et al. (2015)"
bin_res$GWAS <- factor(bin_res$GWAS, levels = c("FinnGen", "Nikpay et al. (2015)"))
bin_res$method <- factor(bin_res$method, levels = c("Weighted mode", "Weighted median", "MR Egger", "Inverse variance weighted"))

library(ggplot2)
ggplot(bin_res, aes(x = or, y = method, color = GWAS)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(xmin = or_lci95, xmax = or_uci95, y=method, color = GWAS), position = position_dodge(width = 0.5), width = 0.4) +
  labs(
    x = "Odds Ratio (95% Confidence Interval)",
    y=""
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("FinnGen" = "mediumpurple4", "Nikpay et al. (2015)" = "darkorange2"))

