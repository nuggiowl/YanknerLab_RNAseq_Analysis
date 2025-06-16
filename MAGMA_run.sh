input_DEG="DEG_gene_set_input.tsv"
input_DEG_FC="DEG_gene_set_input_FCapplied.tsv"
input_KAI_DEG="DEG_gene_set_input_Kai.tsv"
input_KAI_DEG_FC="DEG_gene_set_input_Kai_FCapplied.tsv"
input_KAI_DEG_FC_up="DEG_gene_set_input_Kai_FCapplied_up.tsv"
input_KAI_DEG_FC_down="DEG_gene_set_input_Kai_FCapplied_down.tsv"
input_KAI_DEG_inc_common="DEG_gene_set_input_Kai_inc_common.tsv"
input_KAI_DEG_FC_inc_common="DEG_gene_set_input_Kai_FCapplied_inc_common.tsv"
Macosko_DEG_Abeta="DEG_gene_set_Macosko_Abeta_upreg.tsv"
Macosko_DEG_AbetaTau="DEG_gene_set_Macosko_AbetaTau_upreg.tsv"

#gunzip -c All_20180423.vcf.gz | vcf2bed --sort-tmpdir=tmp --max-mem=30G - > hg19.dbSNP151.bed

awk -F'\t' '{print $22"\t"$28"\t"$37}' GWAS_AD_catalog.tsv > GWAS_AD_catalog_extracted.tsv
#manual curation
#Run R

#Requires long time & memory
#awk 'NR==FNR {h[$1] = 1; next} {if(h[$4]==1) print$0}' GWAS_AD_catalog_rsID.txt hg19.dbSNP151.bed > GWAS_AD_catalog_hg19.bed
awk -F'[\t;=]' 'BEGIN {OFS="\t"}; {print $4,$1,$12}' GWAS_AD_catalog_hg19.bed > GWAS_AD_catalog_LOC_SNPs

~/lib/magma_v1.10/magma --annotate window=35,10 --snp-loc GWAS_AD_catalog_LOC_SNPs --gene-loc NCBI37.3.gene.loc.converted --out GWAS_AD_catalog
~/lib/magma_v1.10/magma --bfile g1000_eur --gene-annot GWAS_AD_catalog.genes.annot --pval GWAS_AD_catalog.summary ncol=N --out GWAS_AD_catalog
~/lib/magma_v1.10/magma --gene-results GWAS_AD_catalog.genes.raw --set-annot $input_DEG col=1,2 --out GWAS_AD_catalog
~/lib/magma_v1.10/magma --gene-results GWAS_AD_catalog.genes.raw --set-annot $input_DEG_FC col=1,2 --out GWAS_AD_catalog_FCapplied
~/lib/magma_v1.10/magma --gene-results GWAS_AD_catalog.genes.raw --set-annot $input_KAI_DEG col=1,2 --out GWAS_AD_catalog_Kai
~/lib/magma_v1.10/magma --gene-results GWAS_AD_catalog.genes.raw --set-annot $input_KAI_DEG_FC col=1,2 --out GWAS_AD_catalog_Kai_FCapplied
~/lib/magma_v1.10/magma --gene-results GWAS_AD_catalog.genes.raw --set-annot $input_KAI_DEG_FC_up col=1,2 --out GWAS_AD_catalog_Kai_FCapplied_up
~/lib/magma_v1.10/magma --gene-results GWAS_AD_catalog.genes.raw --set-annot $input_KAI_DEG_FC_down col=1,2 --out GWAS_AD_catalog_Kai_FCapplied_down
~/lib/magma_v1.10/magma --gene-results GWAS_AD_catalog.genes.raw --set-annot $input_KAI_DEG_inc_common col=1,2 --out GWAS_AD_catalog_Kai_inc_common
~/lib/magma_v1.10/magma --gene-results GWAS_AD_catalog.genes.raw --set-annot $input_KAI_DEG_FC_inc_common col=1,2 --out GWAS_AD_catalog_Kai_FCapplied_inc_common

awk -F'[\t:]' 'NR!=1{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$9+$10}' ../GWAS_summary/nallsEtAl2019_excluding23andMe_allVariants.tab > GWAS_PD_2019.avinput
annotate_variation.pl GWAS_PD_2019.avinput humandb/ -filter -build hg19 -dbtype avsnp150 -out GWAS_PD_2019
awk -F'\t' '{print $2"\t"$3"\t"$4}' GWAS_PD_2019.hg19_avsnp150_dropped | sed 's/chr//g' > GWAS_PD_2019_LOC_SNPs
head -n 1 ../GWAS_summary/nallsEtAl2019_excluding23andMe_allVariants.tab | awk '{print "Build\t"$1"\tCHR\tBEG\tEND\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tP\t"$8"\t"$9"\tN"}' > GWAS_PD_2019.summary
cat GWAS_PD_2019.hg19_avsnp150_dropped >> GWAS_PD_2019.summary
~/lib/magma_v1.10/magma --annotate window=35,10 --snp-loc GWAS_PD_2019_LOC_SNPs --gene-loc NCBI37.3.gene.loc.converted --out GWAS_PD_2019
~/lib/magma_v1.10/magma --bfile g1000_eur --gene-annot GWAS_PD_2019.genes.annot --pval GWAS_PD_2019.summary ncol=N --out GWAS_PD_2019
~/lib/magma_v1.10/magma --gene-results GWAS_PD_2019.genes.raw --set-annot $input_DEG col=1,2 --out GWAS_PD_2019
~/lib/magma_v1.10/magma --gene-results GWAS_PD_2019.genes.raw --set-annot $input_DEG_FC col=1,2 --out GWAS_PD_2019_FCapplied
~/lib/magma_v1.10/magma --gene-results GWAS_PD_2019.genes.raw --set-annot $Macosko_DEG_Abeta col=1,2 --out GWAS_PD_2019_Macosko_Abeta
~/lib/magma_v1.10/magma --gene-results GWAS_PD_2019.genes.raw --set-annot $Macosko_DEG_AbetaTau col=1,2 --out GWAS_PD_2019_Macosko_AbetaTau

awk -F'\t' 'NR!=1{print $2"\t"$1"\t"$3}' ../GWAS_summary/iPSYCH-PGC_ASD_Nov2017 > GWAS_ASD_2017_LOC_SNPs
awk -F'\t' '{print $0"\t46350"}' ../GWAS_summary/iPSYCH-PGC_ASD_Nov2017 > GWAS_ASD_2017.summary
#change column names to P and N
~/lib/magma_v1.10/magma --annotate window=35,10 --snp-loc GWAS_ASD_2017_LOC_SNPs --gene-loc NCBI37.3.gene.loc.converted --out GWAS_ASD_2017
~/lib/magma_v1.10/magma --bfile g1000_eur --gene-annot GWAS_ASD_2017.genes.annot --pval GWAS_ASD_2017.summary ncol=N --out GWAS_ASD_2017
~/lib/magma_v1.10/magma --gene-results GWAS_ASD_2017.genes.raw --set-annot $input_DEG col=1,2 --out GWAS_ASD_2017
~/lib/magma_v1.10/magma --gene-results GWAS_ASD_2017.genes.raw --set-annot $input_DEG_FC col=1,2 --out GWAS_ASD_2017_FCapplied
~/lib/magma_v1.10/magma --gene-results GWAS_ASD_2017.genes.raw --set-annot $Macosko_DEG_Abeta col=1,2 --out GWAS_ASD_2017_Macosko_Abeta
~/lib/magma_v1.10/magma --gene-results GWAS_ASD_2017.genes.raw --set-annot $Macosko_DEG_AbetaTau col=1,2 --out GWAS_ASD_2017_Macosko_AbetaTau

zcat ../GWAS_summary/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz | awk '/^[^#]/ {print}' | awk -F'\t' 'NR!=1{print $2"\t"$1"\t"$3}' > GWAS_SCZ_2022_LOC_SNPs
zcat ../GWAS_summary/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz | awk '/^[^#]/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11"\t"$14+$15}' > GWAS_SCZ_2022.summary
#change column names to SNP, P and N
~/lib/magma_v1.10/magma --annotate window=35,10 --snp-loc GWAS_SCZ_2022_LOC_SNPs --gene-loc NCBI37.3.gene.loc.converted --out GWAS_SCZ_2022
~/lib/magma_v1.10/magma --bfile g1000_eur --gene-annot GWAS_SCZ_2022.genes.annot --pval GWAS_SCZ_2022.summary ncol=N --out GWAS_SCZ_2022
~/lib/magma_v1.10/magma --gene-results GWAS_SCZ_2022.genes.raw --set-annot $input_DEG col=1,2 --out GWAS_SCZ_2022
~/lib/magma_v1.10/magma --gene-results GWAS_SCZ_2022.genes.raw --set-annot $input_DEG_FC col=1,2 --out GWAS_SCZ_2022_FCapplied

zcat ../GWAS_summary/pgc-bip2021-all.vcf.tsv.gz | awk '/^[^#]/ {print}' | awk -F'\t' '{print $3"\t"$1"\t"$2}' > GWAS_BIP_2021_LOC_SNPs
zcat ../GWAS_summary/pgc-bip2021-all.vcf.tsv.gz | awk '/^[^#]/ || /^#[^#]/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$14+$15}' > GWAS_BIP_2021.summary
#change column names to SNP, P and N
~/lib/magma_v1.10/magma --annotate window=35,10 --snp-loc GWAS_BIP_2021_LOC_SNPs --gene-loc NCBI37.3.gene.loc.converted --out GWAS_BIP_2021
~/lib/magma_v1.10/magma --bfile g1000_eur --gene-annot GWAS_BIP_2021.genes.annot --pval GWAS_BIP_2021.summary ncol=N --out GWAS_BIP_2021
~/lib/magma_v1.10/magma --gene-results GWAS_BIP_2021.genes.raw --set-annot $input_DEG col=1,2 --out GWAS_BIP_2021
~/lib/magma_v1.10/magma --gene-results GWAS_BIP_2021.genes.raw --set-annot $input_DEG_FC col=1,2 --out GWAS_BIP_2021_FCapplied
