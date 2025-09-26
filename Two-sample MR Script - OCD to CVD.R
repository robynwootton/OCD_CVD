#Analysis by Robyn Wootton
#Analysis begun March 2024
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
#########             1. Load Required Packages                     ############
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
#########                 2. Read in Exposures                      ############
################################################################################
#The data structure for the two-sample MR package is explained here: https://mrcieu.github.io/TwoSampleMR/articles/index.html

#Read in the OCD instrument (including 23andMe)
ocd <- read.csv("/pathname/OCD_GWS_inc23andMe.csv", header=T)
#Check they are all GWS
summary(ocd$P)
#convert the OR to log odds
ocd$beta <- log(ocd$OR)
#Perform clumping
#name the columns as per ld_clump
ocd$rsid <- ocd$SNP
ocd$pval <- ocd$P
#Subset to instrument
ocd <- subset(ocd, ocd$pval<5e-8)
ocd <- ld_clump(ocd, clump_kb = 1000, clump_r2 = 0.01, opengwas_jwt = XXXX)
#Save the instrument
write.csv(ocd, "/pathname/OCD_GWS_inc23andMe_clumped.csv", quote=F, row.names=F )

#Read in the Instrument using TwoSampleMR
exp_dat <- TwoSampleMR::read_exposure_data(
  filename = "/pathname/OCD_GWS_inc23andMe_clumped.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ_U_2044417",
  pval_col = "P")
exp_dat$exposure <- "OCD"

#add the columns later for steiger
exp_dat$samplesize.exposure<- 2098077
exp_dat$ncase.exposure<-53660
exp_dat$ncontrol.exposure<-2044417
exp_dat$units.exposure<-"log odds"
exp_dat$prevalence.exposure<-0.01

#Calculate variance explained by the instrument in discovery sample
N <- 2098077
exp_dat$r2 <- (2*(exp_dat$beta^2)* exp_dat$eaf*(1-exp_dat$eaf))/
  (2*(exp_dat$beta^2)* exp_dat$eaf*(1-exp_dat$eaf)+(exp_dat$se^2)*2*N* exp_dat$eaf*(1-exp_dat$eaf))
#total variance explained by the instrument:
sum(exp_dat$r2)
sum(exp_dat$r2)*100

##################################################################################
#######                #3. Read in outcome data                            ####### 
##################################################################################
#First of all, I will do the outcomes which are available on MR Base
#ao <- available_outcomes()
#head(ao)

#Tell it where to read from
#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

#Binary
#Coronary Artery Disease ebi-a-GCST003116
#MI - ieu-a-798 
#stroke - ebi-a-GCST006906
#T2D - ebi-a-GCST005047

#Continuous
#Total Cholesterol - ebi-a-GCST002221
#HDL Cholesterol - ebi-a-GCST002223
#LDL Cholesterol - ebi-a-GCST002222
#Triglycerides - ebi-a-GCST002216
#BMI - ebi-a-GCST002783

#extract from MR Base
outcome_dat_bin <- TwoSampleMR::extract_outcome_data(exp_dat$SNP, c("ebi-a-GCST003116", "ieu-a-798", "ebi-a-GCST006906", "ebi-a-GCST005047", "finn-b-I9_MI"))
outcome_dat_cont <- TwoSampleMR::extract_outcome_data(exp_dat$SNP, c("ebi-a-GCST002223", "ebi-a-GCST002222", "ebi-a-GCST002221", "ebi-a-GCST002783", "ebi-a-GCST002216"))

##And extract the outcomes that aren't on MR Base
#Heart failure, shah
hf_out_dat <- TwoSampleMR::read_outcome_data(
  snps = exp_dat$SNP,
  filename = "/pathname/HeartFailure.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p")
hf_out_dat$outcome <- "Heart Failure"

#heart rate variability

hrv_out_dat <- TwoSampleMR::read_outcome_data(
  snps = exp_dat$SNP,
  filename = "/pathname/Nolte_28613276_SDNN.txt",
  sep = " ",
  snp_col = "rs_number",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "reference_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p.value")
hrv_out_dat$outcome <- "Heart Rate Variability"


##Systolic Blood pressure 
sbp_out_dat <- TwoSampleMR::read_outcome_data(
  sep=",",
  snps = exp_dat$SNP,
  filename = "/pathname/SBP_summary.csv",
  snp_col = "V3",
  beta_col = "meta.beta",
  se_col = "meta.se",
  effect_allele_col = "new_REF",
  other_allele_col = "new_ALT",
  eaf_col = "meta.w_eaf",
  pval_col = "meta.pval")
sbp_out_dat$outcome <- "Systolic BP"

#Format effect sizes to be in SD units
#I have calculated the meta-analysed total sample SD to be 19.10
sbp_out_dat$beta.outcome <- sbp_out_dat$beta.outcome/19.1
sbp_out_dat$se.outcome <- sbp_out_dat$se.outcome/19.1

#Diastolic Blood Pressure
dbp_out_dat <- TwoSampleMR::read_outcome_data(
  snps = exp_dat$SNP,
  filename = "/pathname/DBP_summary.csv",
  sep = ",",
  snp_col = "V3",
  beta_col = "meta.beta",
  se_col = "meta.se",
  effect_allele_col = "new_REF",
  other_allele_col = "new_ALT",
  eaf_col = "meta.w_eaf",
  pval_col = "meta.pval")
dbp_out_dat$outcome <- "Diastolic BP"

#Format effect sizes to be in SD units
#I have calculated the meta-analysed total sample SD to be 10.86
dbp_out_dat$beta.outcome <- dbp_out_dat$beta.outcome/10.86
dbp_out_dat$se.outcome <- dbp_out_dat$se.outcome/10.86

##Suicide as a positive control
suicide_out_dat <- TwoSampleMR::read_outcome_data(
  snps = exp_dat$SNP,
  filename = "/pathname/suicide_summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "OR",
  se_col = "SE",
  effect_allele_col = "A1.y",
  other_allele_col = "A2.y",
  eaf_col = "FRQ_U_492022",
  pval_col = "P")
suicide_out_dat$outcome <- "Suicide"
#make into log odds
suicide_out_dat$beta.outcome <- log(suicide_out_dat$beta.outcome)

###Thrombosis
thromb_out_dat <- TwoSampleMR::read_outcome_data(
  snps = exp_dat$SNP,
  filename = "/pathname/thrombosis_summary.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "Freq1",
  pval_col = "P.value")
thromb_out_dat$outcome <- "Thrombosis"

################################################################################
#########      4. Harmonising Exposure and Outcome Data             ############
################################################################################

#harmonise the data, IEU binary
data_bin <- TwoSampleMR::harmonise_data( 
  exposure_dat = exp_dat ,
  outcome_dat = outcome_dat_bin,
  action = 2)

#harmonise the data, IEU cont
data_cont <- TwoSampleMR::harmonise_data( 
  exposure_dat = exp_dat ,
  outcome_dat = outcome_dat_cont,
  action = 2)

#harmonise the data, HF
data_hf <- TwoSampleMR::harmonise_data( 
  exposure_dat = exp_dat ,
  outcome_dat = hf_out_dat,
  action = 2)

#harmonise the data, HRV
data_hrv <- TwoSampleMR::harmonise_data( 
  exposure_dat = exp_dat ,
  outcome_dat = hrv_out_dat,
  action = 2)

#harmonise the data, SBP
data_sbp <- TwoSampleMR::harmonise_data( 
  exposure_dat = exp_dat ,
  outcome_dat = sbp_out_dat,
  action = 2)

#harmonise the data, DBP
data_dbp <- TwoSampleMR::harmonise_data( 
  exposure_dat = exp_dat ,
  outcome_dat = dbp_out_dat,
  action = 2)

#harmonise the data, suicide
data_sui <- TwoSampleMR::harmonise_data( 
  exposure_dat = exp_dat ,
  outcome_dat = suicide_out_dat,
  action = 2)

#harmonise the data, thrombosis
data_thromb <- TwoSampleMR::harmonise_data( 
  exposure_dat = exp_dat ,
  outcome_dat = thromb_out_dat,
  action = 2)

################################################################################
#########                   5. The MR Analysis                      ############
################################################################################

#Cholesterol contains comma so rename
data_cont$outcome[data_cont$outcome=="Cholesterol, total || id:ebi-a-GCST002221"] <- "Cholesterol total || id:ebi-a-GCST002221"

#Separate into each of the analyses
data_cad <- subset(data_bin, data_bin$outcome=="Coronary artery disease || id:ebi-a-GCST003116")
data_t2d <- subset(data_bin, data_bin$outcome=="Type 2 diabetes || id:ebi-a-GCST005047")
data_mi <- subset(data_bin, data_bin$outcome=="Myocardial infarction || id:ieu-a-798")
data_mi_finn <- subset(data_bin, data_bin$outcome=="Myocardial infarction || id:finn-b-I9_MI")
data_stroke <- subset(data_bin, data_bin$outcome=="Stroke || id:ebi-a-GCST006906")
data_bmi <- subset(data_cont, data_cont$outcome=="Body mass index || id:ebi-a-GCST002783")
data_ldl <- subset(data_cont, data_cont$outcome=="LDL cholesterol || id:ebi-a-GCST002222")
data_hdl <- subset(data_cont, data_cont$outcome=="HDL cholesterol || id:ebi-a-GCST002223")
data_chol <- subset(data_cont, data_cont$outcome=="Cholesterol total || id:ebi-a-GCST002221")
data_trig <- subset(data_cont, data_cont$outcome=="Triglycerides || id:ebi-a-GCST002216")

#Create lists of the results
all_bin <- list(data_cad, data_t2d, data_mi, data_stroke, data_hf, data_thromb, data_sui, data_mi_finn)
all_cont <- list(data_bmi, data_ldl, data_hdl, data_chol, data_trig, data_hrv, data_sbp, data_dbp)

library(ggplot2)

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
  res$b <- res$b*0.693 #following Stephen Burgess suggestion so estimates are per doubling in odds for exposure
  or <- TwoSampleMR::generate_odds_ratios(res)
  mr_egger_int <- TwoSampleMR::mr_pleiotropy_test(df)
  mr_egger_int$lci <- mr_egger_int$egger_intercept-(1.96*mr_egger_int$se)
  mr_egger_int$uci <- mr_egger_int$egger_intercept+(1.96*mr_egger_int$se)
  output_het <-rbind(output_het, mr_het)
  output_res <-rbind(output_res, or)
  output_egger <-rbind(output_egger, mr_egger_int)
  write.csv(output_het, "/pathname/heterogeneity_stats_cont.csv", quote=FALSE, row.names=FALSE)
  write.csv(output_res, "/pathname/MRresults_cont.csv", quote=FALSE, row.names=FALSE)
  write.csv(output_egger, "/pathname/egger_cont.csv", quote=FALSE, row.names=FALSE)
  
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
  name_list <- c("data_cad", "data_t2d", "data_mi", "data_stroke", "data_hf", "data_thromb", "data_sui", "data_mi_finn")
  name <- name_list[i]
  df <- data.frame(all_bin[i])
  mr_het <- TwoSampleMR::mr_heterogeneity(df)
  res <- TwoSampleMR::mr(df, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  res$b <- res$b*0.693 #following Stephen Burgess suggestion so estimates are per doubling in odds for exposure
  or <- TwoSampleMR::generate_odds_ratios(res)
  mr_egger_int <- TwoSampleMR::mr_pleiotropy_test(df)
  mr_egger_int$lci <- mr_egger_int$egger_intercept-(1.96*mr_egger_int$se)
  mr_egger_int$uci <- mr_egger_int$egger_intercept+(1.96*mr_egger_int$se)
  output_het <-rbind(output_het, mr_het)
  output_res <-rbind(output_res, or)
  output_egger <-rbind(output_egger, mr_egger_int)
  write.csv(output_het, "/pathname/heterogeneity_stats_bin.csv", quote=FALSE, row.names=FALSE)
  write.csv(output_res, "/pathname/MRresults_bin.csv", quote=FALSE, row.names=FALSE)
  write.csv(output_egger, "/pathname/egger_bin.csv", quote=FALSE, row.names=FALSE)
  
  ## Scatter Plot
  # Shows the distribution of the SNP-exposure effects against the SNP-outcome effects
  # Includes the main MR methods plotted
  p1 <- TwoSampleMR::mr_scatter_plot(res, df)
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

ll <- list( data_cad, data_t2d, data_mi, data_stroke, data_sui, data_bmi,  data_ldl,  data_hdl, 
            data_chol,  data_trig, data_hf, data_cont, data_hrv, data_sbp, data_dbp, data_thromb, data_mi_finn)

output <- c()
for(i in 1: length(ll)){
  name_list <- c("data_cad", "data_t2d", "data_mi", "data_stroke", "data_sui", "data_bmi",  "data_ldl",  "data_hdl", 
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
  write.csv(output, "/pathname/regressiondilution_revision.csv", quote=FALSE, row.names=FALSE)
} 


################################################################################
#########              7. Test for reverse causation                ############
################################################################################

## Steiger filtering
#continuous exposures/outcomes must contain: samplesize and units="SD"
#binary exposures/outcomes must contain: samplesize, ncase, ncontrol, prevalence and units="log odds"

#CAD
data_cad $samplesize.outcome<- 86995
data_cad $ncase.outcome<-22233
data_cad $ncontrol.outcome<- 64762
data_cad $units.outcome<-"log odds"
data_cad $prevalence.outcome<-0.05

#T2D
data_t2d $samplesize.outcome<- 69033
data_t2d $ncase.outcome <- 12171
data_t2d $ncontrol.outcome<- 	56862
data_t2d $units.outcome<-"log odds"
data_t2d $prevalence.outcome<-0.06
data_t2d $eaf.outcome <- data_t2d $eaf.exposure

#MI
data_mi $samplesize.outcome<- 171875
data_mi $ncase.outcome<-48371
data_mi $ncontrol.outcome<- 123504
data_mi $units.outcome<-"log odds"
data_mi $prevalence.outcome<-0.017

#Stroke
data_stroke $samplesize.outcome<- 521612
data_stroke $ncase.outcome<-67162
data_stroke $ncontrol.outcome<- 454450
data_stroke $units.outcome<-"log odds"
data_stroke $prevalence.outcome<-0.25

#Heart Failure
data_hf $samplesize.outcome<- 977323
data_hf $ncase.outcome<-47309
data_hf $ncontrol.outcome<- 930014
data_hf $units.outcome<-"log odds"
data_hf $prevalence.outcome<-0.023

#BMI
data_bmi $samplesize.outcome<- 339224
data_bmi $units.outcome<-"SD"

#LDL
data_ldl $samplesize.outcome<- 188578
data_ldl $units.outcome<-"SD"

#HDL
data_hdl $samplesize.outcome<- 188578
data_hdl $units.outcome<-"SD"

#Cholesterol
data_chol $samplesize.outcome<- 188578
data_chol $units.outcome<-"SD"

#Heart Rate Variability
data_hrv $samplesize.outcome<- 53174
data_hrv $units.outcome<-"SD"

#Triglicerides
data_trig $samplesize.outcome<- 188578
data_trig $units.outcome<-"SD"

#Systolic Blood Pressure
data_sbp $samplesize.outcome<- 150134
data_sbp $units.outcome<-"SD"

#Diastolic Blood Pressure
data_dbp $samplesize.outcome<- 150134
data_dbp $units.outcome<-"SD"

#Thrombosis
data_thromb $samplesize.outcome<- 995864
data_thromb $ncase.outcome<-62879
data_thromb $ncontrol.outcome<- 932985
data_thromb $units.outcome<-"log odds"
data_thromb $prevalence.outcome<-0.002

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
  steiger <- steiger_filtering(analysis_df)
  false <- length(steiger$steiger_dir[steiger$steiger_dir == FALSE])
  true <- length(steiger$steiger_dir[steiger$steiger_dir == TRUE])
  percent <- (true/(true+false))*100
  filter <- cbind(name, false, true, percent)
  colnames(filter) <- c("data", "falseN", "trueN", "percent_true") #following the MRBase naming convention
  output <-rbind(output, filter)
  colnames(output) <- c("data", "falseN", "trueN", "percent_true")
  write.csv(output, "/pathname/steiger.csv", quote=FALSE, row.names=FALSE)
}

################################################################################
#########                   8. Visual Inspection                    ############
################################################################################

ll <- list( data_cad, data_t2d, data_mi, data_stroke, data_bmi,  data_ldl,  data_hdl, 
            data_chol,  data_trig, data_hf, data_cont, data_hrv, data_sbp, data_dbp, data_thromb, data_sui)

for(i in 1: length(ll)){
  name_list <- c("data_cad", "data_t2d", "data_mi", "data_stroke", "data_bmi",  "data_ldl",  "data_hdl", 
                 "data_chol",  "data_trig", 'data_hf', "data_cont", "data_hrv", "data_sbp", "data_dbp", "data_thromb", "data_sui")
  name <- name_list[i]
  df <- data.frame(ll[i])
  
  ## Leave-one-out Analysis
  # Re-runs the IVW analysis, leaving one SNP out at a time
  # Check effects are not driven by a single SNP, estimate should stay relatively consistent 
  res_loo <- mr_leaveoneout(df)
  p2 <- mr_leaveoneout_plot(res_loo)
  loo_file <- paste0("/pathname/", name, "_loo.png")
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

##Continuous outcomes
cont_res <- read.csv("/pathname/MRresults_cont.csv", header=T)
head(cont_res)
#remove Egger because the I2 is very low
cont_res <- subset(cont_res, cont_res$method!="MR Egger")

#Binary Outcomes
bin_res <- read.csv("/pathname/MRresults_bin.csv", header=T)
head(bin_res)
#remove the MR Egger Results due to low I2
bin_res <- subset(bin_res, bin_res$method!="MR Egger")
table(bin_res$outcome)
bin_res <- subset(bin_res, bin_res$outcome!="Myocardial infarction || id:finn-b-I9_MI")

##Plot of just IVW results
##Continuous
table(cont_res$method)
cont_res_ivw <- subset(cont_res, cont_res$method=="Inverse variance weighted")
bin_res_ivw <- subset(bin_res, bin_res$method=="Inverse variance weighted")
cont_res_ivw$outcome <- gsub("\\|.*","",cont_res_ivw$outcome)
bin_res_ivw$outcome <- gsub("\\|.*","",bin_res_ivw$outcome)

table_text <- cbind(
  c("Outcome", cont_res_ivw$outcome),
  c("Beta", sprintf("%.3f", cont_res_ivw$b)),
  c("95% CI", paste(sprintf("%.3f", cont_res_ivw$lo_ci), "-", sprintf("%.3f", cont_res_ivw$up_ci)))
)

summary(cont_res_ivw$lo_ci)
summary(cont_res_ivw$up_ci)

pdf.options(reset = TRUE, onefile = FALSE)
pdf("ocd-cvd_cont_IVW.pdf", width=6,height=7)
# Create the forest plot
forestplot(
  labeltext = table_text,
  mean = c(NA, cont_res_ivw$b), # Add NA for the header row
  lower = c(NA, cont_res_ivw$lo_ci),
  upper = c(NA, cont_res_ivw$up_ci),
  new_page = TRUE,
  
  # Custom appearance
  is.summary = c(TRUE, rep(FALSE, 8)), # Summary at first and last rows
  xlab = "Beta (95% CI)",
  clip =c(-0.2, 0.2), 
  
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

##Binary Outcomes
table_text <- cbind(
  c("Outcomes", bin_res_ivw$outcome),
  c("OR", sprintf("%.2f", bin_res_ivw$or)),
  c("95% CI", paste(sprintf("%.2f", bin_res_ivw$or_lci95), "-", sprintf("%.2f", bin_res_ivw$or_uci95)))
)

summary(bin_res_ivw$or_lci95)
summary(bin_res_ivw$or_uci95)

pdf.options(reset = TRUE, onefile = FALSE)
pdf("ocd-cvd_bin_IVW.pdf", width=6,height=6)
# Create the forest plot
forestplot(
  labeltext = table_text,
  mean = c(NA, bin_res_ivw$or), # Add NA for the header row
  lower = c(NA, bin_res_ivw$or_lci95),
  upper = c(NA, bin_res_ivw$or_uci95),
  new_page = TRUE,
  zero=1,
  
  # Custom appearance
  is.summary = c(TRUE, rep(FALSE, 7)), # Summary at first and last rows
  xlab = "Odds Ratio",
  clip =c(0.8, 1.5), 
  
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

#Compare MI before and after
bin_res <- read.csv("/pathname/MRresults_bin.csv", header=T)
head(bin_res)
#Subset to MI
table(bin_res$outcome)
bin_res <- subset(bin_res, bin_res$outcome== "Myocardial infarction || id:finn-b-I9_MI" | bin_res$outcome=="Myocardial infarction || id:ieu-a-798")
bin_res <- subset(bin_res, bin_res$method!="MR Egger")
bin_res$GWAS[bin_res$outcome== "Myocardial infarction || id:finn-b-I9_MI"] <- "FinnGen"
bin_res$GWAS[bin_res$outcome== "Myocardial infarction || id:ieu-a-798"] <- "Nikpay et al. (2015)"
bin_res$GWAS <- factor(bin_res$GWAS, levels = c("FinnGen", "Nikpay et al. (2015)"))
bin_res$method <- factor(bin_res$method, levels = c("Weighted mode", "Weighted median", "Inverse variance weighted"))

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

