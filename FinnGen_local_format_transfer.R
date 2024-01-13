rm(list = ls())
workpath="G:\\Others\\Program\\MR\\Outcome\\PLC-FinnGen";setwd(workpath)  #set work path
df <- vroom::vroom("finngen_R9_C3_HEPATOCELLU_CARC_EXALLC.gz")
df2 <- df[1:10,]
names(df2)
View(df2)

######As Exposure
exp_finngen_R9 <- subset(df, df$pval<0.0005)
exp_finngen_R9 <- TwoSampleMR::format_data(
  exp_finngen_R9,
  type = "exposure",
  phenotype_col = "Phenotype",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  eaf_col = "af_alt_cases",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  gene_col = "nearest_genes",
  chr_col = "#chrom",
  pos_col = "pos"
)

#####As Outcome
out_finngen_R9 <- TwoSampleMR::format_data(
  df,
  type = "outcome",
  phenotype_col = "Phenotype",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  eaf_col = "af_alt_cases",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  gene_col = "nearest_genes",
  chr_col = "#chrom",
  pos_col = "pos"
)

# Add manifest information
# load manifest
manifest <- read.delim("R9_manifest.tsv")
# Extract file name information
manifest2 <- manifest[
  grep(pattern = "finngen_R9_C3_HEPATOCELLU_CARC_EXALLC.gz", #your file
       manifest$path_https),]
manifest2
out_finngen_R9$ncase.outcome <- manifest2$num_cases
out_finngen_R9$ncontrol.outcome <- manifest2$num_controls
out_finngen_R9$samplesize.outcome <- manifest2$num_cases + manifest2$num_controls
out_finngen_R9$outcome <- manifest2$name
df3 <- out_finngen_R9[1:10,]

#exp_finngen_R9$ncase.outcome <- manifest2$num_cases
#exp_finngen_R9$ncontrol.outcome <- manifest2$num_controls
#exp_finngen_R9$samplesize.outcome <- manifest2$num_cases + manifest2$num_controls
#exp_finngen_R9$outcome <- manifest2$name
#df3 <- exp_finngen_R9[1:10,]

library(dplyr)
colnames(out_finngen_R9)
colnames(out_finngen_R9)=c("chr.outcome","pos.outcome","A2","A1","SNP","gene.outcome","p","b","se","freq","outcome","mr_keep.outcome","pval_origin.outcome","id.outcome","ncase.outcome","ncontrol.outcome","n")

colnames(exp_finngen_R9)
colnames(exp_finngen_R9)=c("chr.exposure","pos.exposure","A2","A1","SNP","gene.outcome","p","b","se","freq","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure","ncase.exposure","ncontrol.exposure","n")

#Result output
write.table(exp_finngen_R9,file = "Finngen_R9_HCC_exp.txt", quote = F, sep = "\t",row.names = F)
write.csv(exp_finngen_R9,file = "Finngen_R9_HCC_exp.csv")
