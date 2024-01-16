rm(list = ls())

#install.packages("shiny")
#install.packages("devtools")
#devtools::install_github("explodecomputer/plinkbinr",force = TRUE)
#devtools::install_github("MRCIEU/TwoSampleMR")
#remotes::install_github("phenoscanner/phenoscanner")
#install.packages("data.table")
#install.packages("dplyr")
#install.packages("tidyverse")
#install.packages("ieugwasr")
#install.packages("readxl")
#install.packages("writexl")
#install.packages("MendelianRandomization")
#install.packages("purrr")
#install.packages("readr")

setwd(dir="G:\\Others\\Program\\MR\\Exposure\\immune\\data\\Expresso") #set work path

library(phenoscanner)
library(friendly2MR)
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(ieugwasr)
library(plinkbinr)
library(pacman)
library(MendelianRandomization)
library(purrr)
library(readr)
get_plink_exe()

#set exposure
FileNames <-list.files(paste0(getwd()),pattern=".txt.gz")
exp_dat_ids <- str_sub(FileNames,start = 1,end = 12)
exps <- str_sub(FileNames,start = 1,end = 12)

#set outcome
out<-fread("Finngen_R9_HCC.txt",header = T) #local outcome
out$trait <- 'Hepatocellular Carcinoma'     #outcome name

outcomeid <- out
rm(out)
head(outcomeid)

###############loop code######################
qaq <- 1
for (qaq in 1:length(exp_dat_ids)) {  
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  
  d3<- try(fread(paste0(getwd(),"/",FileNames[qaq]),fill=TRUE),silent = T)
  
  d3$PHENO <- exps[qaq]
  d3<-format_data(d3,
                 type="exposure",
                 phenotype_col = "PHENO",
                 snp_col = "SNP",
                 beta_col = "beta",
                 se_col = "standard_error",
                 pval_col = "p_value",
                 samplesize_col = "n",
                 eaf_col = "effect_allele_frequency",
                 effect_allele_col = "effect_allele",
                 other_allele_col = "other_allele")
 
# exp_data <- clump_data(d3,clump_kb = 500,clump_r2 = 0.01)

#local clump     
 d4<- ld_clump(
    #dat = X1,
    clump_kb = 500,
    clump_r2 = 0.01,
    pop = "EUR",
    dplyr::tibble(rsid=d3$SNP, pval=d3$pval.exposure, id=d3$id.exposure),
    #get_plink_exe()
    plink_bin = "D:/R/R-4.3.1/library/plinkbinr/bin/plink_Windows.exe",
    #Genome Ref
    bfile = "G:/Others/Program/MR/Exposure/immune/1kg.v3/EUR"
  )
  
  exp_data<-subset(d3,SNP %in% d4$rsid) 
    
  if(length(exp_data[,1])>0){
    outcome_dat<-merge(exp_data,outcomeid,by.x = "SNP",by.y = "SNP")
    write.csv(outcome_dat,file = "d.csv")
    out_data <- read_outcome_data(
      snps = exp_data$SNP, 
      filename = "d.csv",
      sep = ",",
      phenotype_col = "trait",
      snp_col = "SNP",
      beta_col = "b",
      se_col = "se",
      eaf_col="freq",
      samplesize_col = "n",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      pval_col = "p")

    if(length(out_data[,1])>0){  
      dat <- TwoSampleMR::harmonise_data(
        exposure_dat = exp_data,
        outcome_dat = out_data)
      
      ####Delete Palindromic sequence
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
                                                 "Diabetes mellitus type 2",
                                                 "Type II diabetes",
                                                 "Type 2 diabetes",
                                                 "Type 2 diabetes combined control dataset",
                                                 "Insulin-dependent diabetes mellitus",
                                                 "Diabetes diagnosed by doctor",
                                                 "Self-reported diabetes",
                                                 "Alcohol dependence symptom count",
                                                 "Alcohol consumption drinkers vs non drinkers",
                                                 "Alcohol consumption heavy vs lightnon drinkers",
                                                 "Alcohol consumption",
                                                 "Alcohol consumption in current drinkers",
                                                 "Alcohol consumption drinks per week",
                                                 "Alcohol consumption over the past year",
                                                 "Cause of death: fatty liver",
                                                 "Nonalcoholic fatty liver disease",
                                                 "Fibrosis and cirrhosis of liver",
                                                 "Liver cirrhosis",
                                                 "Hepatitis C induced liver fibrosis",
                                                 "Hepatitis C virus HCV induced liver fibrosis",
                                                 "Self-reported hepatitis c",
                                                 "Hepatitis c",
                                                 "Hepatitis c chronic",
                                                 "Chronic hepatitis C infection",
                                                 "Liver fibrosis severity in HIVhepatitis C co infection",
                                                 "Self-reported hepatitis b",
                                                 "Chronic hepatitis B progressio",
                                                 "Chronic hepatitis B progression in individuals age 35 or older",
                                                 "Chronic hepatitis B infection",
                                                 "Hepatitis b chronic"),
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
      
      print(paste0(exp_dat_ids[qaq],"_SNP_No_",res$nsnp[1]))
      
      results <- TwoSampleMR::generate_odds_ratios(res)
      
      results$estimate <- paste0(
        format(round(results$or, 2), nsmall = 2), " (", 
        format(round(results$or_lci95, 2), nsmall = 2), "-",
        format(round(results$or_uci95, 2), nsmall = 2), ")")
      resdata <- dat
      openxlsx::write.xlsx(dat,file = paste0(exp_dat_ids[qaq],"-dat.xlsx"), rowNames = FALSE)
      
      names(resdata)
      Assumption13 <- subset(resdata,mr_keep==TRUE,
                             select = c("SNP","pval.exposure",
                                        "pval.outcome", # "F_statistic",
                                        "mr_keep"))
      
      openxlsx::write.xlsx(x = list(
        "main"=results,
        "Assumption13"=Assumption13),
        overwrite = TRUE,
        paste0(exp_dat_ids[qaq],"-res.xlsx"))
      
    }}
  
  if(length(dat[,1])>2){
    res_hete <- TwoSampleMR::mr_heterogeneity(dat)
    res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    res_leaveone <- mr_leaveoneout(dat)   
 
        #MR—presso 
            res_presso <- TwoSampleMR::run_mr_presso(dat,NbDistribution = 1000)#[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
            sink(paste0(exp_dat_ids[qaq],"_PRESSO.txt"),append=FALSE,split = FALSE) 
            print(res_presso)
            sink()
            print(res_presso)
    
    p1 <- mr_scatter_plot(res, dat)
    p1[[1]]
    pdf(paste0(exp_dat_ids[qaq],"_scatter.pdf"))
    print(p1[[1]])
    dev.off()
    
    res_single <- mr_singlesnp(dat)
    p2 <- mr_forest_plot(res_single)
    pdf(paste0(exp_dat_ids[qaq],"_forest.pdf"))
    print(p2[[1]])
    dev.off()
    
    p3 <- mr_funnel_plot(res_single)
    pdf(paste0(exp_dat_ids[qaq],"_funnel.pdf"))
    print(p3[[1]])
    dev.off()
    
    res_loo <- mr_leaveoneout(dat)
    pdf(paste0(exp_dat_ids[qaq],"_leave_one_out.pdf"))
    print(mr_leaveoneout_plot(res_loo))
    dev.off()
    
    library(magrittr)
    res3 <- res[1:3,]
    
    #Format transform
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
    
    write.csv(res_ALL,file = paste0(exp_dat_ids[qaq],".csv"), row.names = FALSE)
  }}

#Results output 
file.remove("d.csv")
fs=list.files("./", pattern = "csv",full.names = TRUE) 
df = map_dfr(fs, read.csv)
write.csv(df,"FinnGen_Immune_to_HCC_500_0.01_1E-5_Expresso.csv")
