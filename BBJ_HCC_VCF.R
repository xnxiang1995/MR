#install.packages("vcfR")
#install.packages("gwasvcf")

#Load VCF file
rm(list = ls())
workpath="G:\\Others\\Program\\MR\\Outcome\\PLC-JENGER";setwd(workpath)  #set work path
data=VariantAnnotation::readVcf("bbj-a-158.vcf.gz")
df=gwasvcf::vcf_to_tibble(data, id="bbj-a-158")
df1<-dat[1:10,]

# get colnames
colnames(df)

#preprocessing, insert phenotye, case number, transfer pval etc.
df$Phenotype <- "Hepatocellular Carcinoma"
df$P <- 10^(-df$LP)
df$SS <- 197611
df$Cases <- 1866
df$Control <- 195745

# As exposure
BBJ.HCC_exp <- subset(df, df$P<0.0005)
BBJ.HCC_exp<- TwoSampleMR::format_data(
  BBJ.HCC_exp,
  type = "exposure",
  phenotype_col = "Phenotype",
  snp_col = "rsid",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P",
  gene_col = "nearest_genes",
  chr_col = "seqnames",
  pos_col = "start",
  samplesize_col = "SS",
  ncase_col = "Cases",
  ncontrol_col = "Control"
)

# As outcome
BBJ.HCC_out <-df
BBJ.HCC_out <- TwoSampleMR::format_data(
  BBJ.HCC_out,
  type = "outcome",
  phenotype_col = "Phenotype",
  snp_col = "rsid",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P",
  gene_col = "nearest_genes",
  chr_col = "seqnames",
  pos_col = "start",
  samplesize_col = "SS",
  ncase_col = "Cases",
  ncontrol_col = "Control"
)

# Change name
library(dplyr)
colnames(BBJ.HCC_exp)
colnames(BBJ.HCC_exp)=c("chr.exposure","pos.exposure","A2","A1","b","se","freq","n","id.exposure","SNP","p","ncase.exposure","ncontrol.exposure","exposure","mr_keep.exposure","pval_origin.exposure")

colnames(BBJ.HCC_out)
colnames(BBJ.HCC_out)=c("chr.outcome","pos.outcome","A2","A1","b","se","freq","n","id.outcome","SNP","p","ncase.outcome","ncontrol.outcome","outcome","mr_keep.outcome","pval_origin.outcome")

# Result output
write.table(BBJ.HCC_exp,file = "BBJ.HCC_exp.txt", quote = F, sep = "\t",row.names = F)
write.csv(BBJ.HCC_exp,file = "BBJ.HCC_exp.csv")
write.table(BBJ.HCC_out,file = "BBJ.HCC_out.txt", quote = F, sep = "\t",row.names = F)
write.csv(BBJ.HCC_out,file = "BBJ.HCC_out.csv")
