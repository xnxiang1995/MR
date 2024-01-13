#install.packages("shiny")
#install.packages("devtools")
#devtools::install_github("explodecomputer/plinkbinr")
#get_plink_exe()
#devtools::install_github("MRCIEU/TwoSampleMR")
#install.packages("data.table")
#install.packages("dplyr")
#install.packages("tidyverse")
#install.packages("ieugwasr")
#install.packages("readxl")
#install.packages("writexl")

rm(list = ls())
setwd(dir="G:\\Others\\Program\\MR\\Exposure\\immune\\731\\ALL")  #Set work path

library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(ieugwasr)
library(plinkbinr)
library(stringr)

exp <- fread("Finngen_R9_HCC_exp.txt",header = T)
exp$trait <- 'Hepatocellular Carcinoma'  #Exposure name
exp$n <- as.numeric(exp$n)          
head(exp)

#p threshold
exp<-subset(exp,p<5e-8)  
write.csv(exp, file="exp_clumped.csv")

##Transfer to TSMR format
exp_data<-read_exposure_data(
  filename = "exp_clumped.csv",
  sep = ",",
  phenotype_col = "trait",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  samplesize_col = "n",
  eaf_col = "freq",
  pval_col = "p")

#Online LD
#exp_data<-clump_data(exp,clump_kb =10000,clump_r2= 0.001,clump_p1= 1,clump_p2= 1,pop ="EUR")

#local clump
d4<- ld_clump(
  #dat = X1,
  clump_kb = 500,
  clump_r2 = 0.01,
  pop = "EUR",
  dplyr::tibble(rsid=exp_data$SNP, pval=exp_data$pval.exposure, id=exp_data$id.exposure),
  #get_plink_exe()
  plink_bin = "D:/R/R-4.3.1/library/plinkbinr/bin/plink_Windows.exe",
  #Ref Genome
  bfile = "G:\\Others\\Program\\MR\\Exposure\\immune\\1kg.v3\\EUR"
)

exp_data<-subset(exp_data,SNP %in% d4$rsid) 


#Set outcome
FileNames <-list.files(paste0(getwd()),pattern=".rds")
outcomeids <- str_sub(FileNames,start = 1,end = 21)
out_comes <- outcomeids


#############loop code######################
qaq <- 1
for (qaq in 1:length(outcomeids)) { 
outcomeid <- outcomeids[qaq]
out_come <- out_comes[qaq]

if(length(exp_data[,1])>2){
  
  out<- try(read_rds(paste0(getwd(),"/",FileNames[qaq])),silent = T)
  out$PHENO <- FileNames[qaq]
  outcomeid <- out
  outcome_dat<-merge(exp_data,outcomeid,by.x = "SNP",by.y = "SNP")
  head(out)
  out_data<-format_data(outcome_dat,
                        type="outcome", snps = NULL,
                        phenotype_col = "PHENO",
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "standard_error",
                        pval_col = "P-value",
                        samplesize_col = "n",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",
                        eaf_col = "effect_allele_frequency")
  
  if(length(out_data[,1])>0){  
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exp_data,
      outcome_dat = out_data)
    
    ####Remove palindromic sequences
    dat <-subset(dat,mr_keep==TRUE)
    
    #Calculate F value and R2, delete weak IV
    get_f<-function(dat,F_value=10){
      log<-is.na(dat$eaf.exposure)
      log<-unique(log)
      if(length(log)==1)
      {if(log==TRUE){
        print("data don't have eaf，can not calculate F")
        return(dat)}
      }
      if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("data don't have beta，can not calculate F")
        return(dat)}
      if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("data don't have se，can not calculate F")
        return(dat)}
      if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("data don't have sample size，can not calculate F")
        return(dat)}
      if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
        R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
        F<- (dat$samplesize.exposure-2)*R2/(1-R2)
        dat$R2<-R2
        dat$F<-F
        dat<-subset(dat,F>F_value)
        return(dat)
      }
    }
    dat <- get_f(dat, F_value = 10)
    
    #Search confounders via Phenoscanner
    colnames(dat)
    grp_size<-100
    grps<-split(dat$SNP,ceiling(seq_along(dat$SNP)/grp_size))
    result<-map_dfr(grps,~phenoscanner(snpquery=.x,
                                       catalogue = "GWAS",
                                       pvalue=1e-05,
                                       proxies = "EUR",
                                       r2=0.8,
                                       build = 38)$results)
    #Delete confounders by key words
    delete_data_confounder <- function(confounder,pheno_data,data,exp,out){
      confounder_SNP <- list()
      for (i in 1:length(confounder)) {
        confounder_SNP[[i]]<- subset(pheno_data,pheno_data$trait==confounder[i])
      }
      confounder_SNP <-dplyr::bind_rows(confounder_SNP)%>%
        dplyr::select(snp)%>%
        dplyr::distinct()
      print(paste0("number of confounder SNP:",nrow(confounder_SNP)))
      write.csv(confounder_SNP,file = paste0('~/',deparse(substitute(exp)),"_",deparse(substitute(out)),'_confounder_SNP.csv'),row.names=FALSE,
                quote = FALSE)
      dat <- data%>%
        dplyr::filter(!SNP %in% confounder_SNP$snp)
    }
    dat <- delete_data_confounder(confounder=c("Body mass index",
                                               "Body mass index in females less than or equal to 50 years of age",
                                               "Body mass index in females greater than 50 years of age",
                                               "Body mass index in males greater than 50 years of age",
                                               "Body mass index in males",
                                               "Body mass index in females",
                                               "Body mass index in non-smokers",
                                               "Body mass index adjusted for smoking",
                                               "Body mass index adjusted for physical activity",
                                               "Past tobacco smoking",
                                               "Current tobacco smoking",
                                               "Pack years adult smoking as proportion of life span exposed to smoking",
                                               "Ever smoked",
                                               "Tobacco smoking: smokes on most or all days",
                                               "Smoking status: previous",
                                               "Tobacco smoking: ex-smoker",
                                               "Number of cigarettes previously smoked daily",
                                               "Nicotine dependence smoking cigarettes per day",
                                               "Pack years of smoking preview only",
                                               "Type II diabetes",
                                               "Type 2 diabetes",
                                               "Type 2 diabetes combined control dataset",
                                               "Insulin-dependent diabetes mellitus",
                                               "Diabetes diagnosed by doctor",
                                               "Self-reported diabetes",
                                               "viral infection"),
                                  pheno_data=result,
                                  data = dat,
                                  exp = exp_data,
                                  out = out_data)

    #MR analysis
    res=TwoSampleMR::mr(dat,method_list= c("mr_ivw" ,
                                           "mr_weighted_median" ,
                                           "mr_egger_regression",
                                           "mr_simple_mode",
                                           "mr_weighted_mode",
                                           "mr_wald_ratio"))
    print(paste0(out_come,"_SNP数_",res$nsnp[1]))
    results <- TwoSampleMR::generate_odds_ratios(res)
    results$estimate <- paste0(
      format(round(results$or, 2), nsmall = 2), " (", 
      format(round(results$or_lci95, 2), nsmall = 2), "-",
      format(round(results$or_uci95, 2), nsmall = 2), ")")
    resdata <- dat
    openxlsx::write.xlsx(dat, file = paste0(out_come,"-dat",".xlsx"), rowNames = FALSE)
    
    names(resdata)
    Assumption13 <- subset(resdata,mr_keep==TRUE,
                           select = c("SNP","pval.exposure",
                                      "pval.outcome", # "F_statistic",
                                      "mr_keep"))
    openxlsx::write.xlsx(x = list(
      "main"=results,
      "Assumption13"=Assumption13),
      overwrite = TRUE,
      paste0(out_come,"-res.xlsx"))
    
  }}


if(length(dat[,1])>2){
  res_hete <- TwoSampleMR::mr_heterogeneity(dat)
  res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
  res_leaveone <- mr_leaveoneout(dat)   
  

  #MR—presso 
            #res_presso <- TwoSampleMR::run_mr_presso(dat,NbDistribution = 1000)
                         #  [["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
            #sink(paste0(out_come,"_PRESSO.txt"),append=FALSE,split = FALSE) 
            #print(res_presso)
            #sink()
            #print(res_presso)
  p1 <- mr_scatter_plot(res, dat)
  p1[[1]]
  pdf(paste0(out_come,"_scatter.pdf"))
  print(p1[[1]])
  dev.off()
  
  res_single <- mr_singlesnp(dat)
  p2 <- mr_forest_plot(res_single)
  pdf(paste0(out_come,"_forest.pdf"))
  print(p2[[1]])
  dev.off()
  
  p3 <- mr_funnel_plot(res_single)
  pdf(paste0(out_come,"_funnel.pdf"))
  print(p3[[1]])
  dev.off()
  
  res_loo <- mr_leaveoneout(dat)
  pdf(paste0(out_come,"_leave_one_out.pdf"))
  print(mr_leaveoneout_plot(res_loo))
  dev.off()
  
  library(magrittr)
  res3 <- res[1:3,]
  
  #Data format transform
  library(magrittr)
  # Main result 
  res4 <- tidyr::pivot_wider(
    res3,names_from ="method",names_vary = "slowest",
    values_from = c("b","se","pval") )
  # Heterogeneity statistics
  res_hete2 <- tidyr::pivot_wider(
    res_hete,names_from ="method",names_vary = "slowest",
    values_from = c("Q","Q_df","Q_pval") ) %>% 
    dplyr::select( -id.exposure,-id.outcome,-outcome,-exposure)
  # Horizontal pleiotropy
  res_plei2 <- dplyr::select(res_plei,
                             egger_intercept,se,pval)
  
  # Merge
  res_ALL <- cbind(res4, res_hete2, res_plei2)
  
  write.csv(res_ALL,file = paste0(out_come,".csv"), row.names = FALSE)
}}


#Result output
file.remove("exp_clumped.csv")
fs=list.files("./", pattern = ".csv",full.names = TRUE)   
df = map_dfr(fs, read.csv)
write.csv(df,"FinnGen_HCC_to_immune_5E-8_500_0.01.csv")

