library(dplyr)

extract_numbers <- function(x) {
  x <- gsub("[^0-9,]", "", x)
  # Use regular expression to find digits and capture them in a group
  strsplit(x, ",", fixed = TRUE)[[1]] %>%
    # Coerce to numeric, handling NA for non-numeric values
    as.numeric(., na.rm = T) %>%
    na.omit(.) %>%
    # Sum the extracted numbers
    sum(.)
}

df_GWAS<-read.table("GWAS_AD_catalog_extracted.tsv",sep="\t",header=T)
df_study<-read.table("GWAS_AD_catalog_studies.tsv",sep="\t",header=T)

df_tmp<-df_GWAS %>% filter(grepl("^rs",SNPS))
df_GWAS_filt<-df_tmp %>% filter(!grepl(" ",SNPS))


df_GWAS_filt<-df_GWAS_filt[order(df_GWAS_filt[,1],df_GWAS_filt[,2]),]
df_GWAS_filt<-df_GWAS_filt %>% distinct(SNPS, .keep_all=T)
write.table(df_GWAS_filt[,1],"GWAS_AD_catalog_rsID.txt",quote=F,col.names=F,row.names=F)
colnames(df_GWAS_filt)<-c("SNP","P","STUDY.ACCESSION")

df_study<-df_study[,c("accessionId","discoverySampleAncestry")]
colnames(df_study)<-c("STUDY.ACCESSION","N_raw")

df_study$N <- sapply(as.character(df_study[, 2]), extract_numbers)
df_study$N[df_study$N==0]<-NA

df_m<-merge(df_GWAS_filt,df_study,all.x=T)
write.table(df_m[,-4],"GWAS_AD_catalog.summary",sep="\t",quote=F,row.names=F)
