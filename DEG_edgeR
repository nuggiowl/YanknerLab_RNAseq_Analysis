library("edgeR")
library("ggplot2")
library("ggrepel")

directory <- "."
sampleFiles <- list.files(directory,pattern="\\.htseq$")
sampleFilesTable <- cbind(sampleName=sapply(strsplit(sampleFiles,'_'),function(x) paste(x[1:4], collapse = "_")), fileName=sampleFiles)
colnames(sampleFilesTable)[1]<-"ID"
sampleTable <- read.table("./sampledata.csv",sep="\t",header=T)
colnames(sampleTable)[1]<-"ID"
mergedTable <- merge(sampleFilesTable, sampleTable[,c("ID","group","sequencingBatch")])
colnames(mergedTable)[3]<-"Group"
colnames(mergedTable)[4]<-"Batch"
d<-readDGE(sampleFiles,path=directory,columns=c(1,2),header=F,group=mergedTable$Group)
x=DGEList(counts=d$counts,group=d$samples$group,samples=mergedTable)
y<-x[,which(!x$samples$Group=="wt.water")]
y$samples$group<-factor(y$samples$group,levels=c("ad3xtg.water","ad3xtg.lithiumorotate0.7"))

# filter for expressed genes
minCPM=1
# determine the min number of samples in a group.
GroupSize=3 #set this value
keep=(rowSums(cpm(y)>=minCPM)>=GroupSize)
y=y[keep, , keep.lib.sizes=TRUE]

#normalize data
y=calcNormFactors(object=y, method="TMM")

logcpm=cpm(y=y, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1)
write.table(logcpm,"res_normalizedCPM_EdgeR_1_3xTG_LiO_water.txt",sep="\t")

#RUVR method `doRUVr` with specified k `RUVrK`
library(RUVSeq)
RUVrK=1
design.tmp=model.matrix(as.formula("~Batch+group"), data=y$samples)
colnames(design.tmp)=gsub(":", ".", colnames(design.tmp))
colnames(design.tmp)=gsub("\\(|\\)", "", colnames(design.tmp))
colnames(design.tmp)=make.names(colnames(design.tmp))
z=estimateDisp(y=y, design=design.tmp, robust=FALSE)
#fit model
fit.tmp=glmFit(y=z, design=design.tmp, robust=FALSE)

contrasts=makeContrasts(
    ad3xtg.lithiumorotate0.7vad3xtg.water=groupad3xtg.lithiumorotate0.7,
    levels=design.tmp)

res=glmLRT(glmfit=fit.tmp, contrast=contrasts)
tab_raw=as.data.frame(topTags(object=res, n=Inf, adjust.method="fdr", sort.by="PValue")$table)
write.table(tab_raw,"res_EdgeR_3xTG_LiO_water.txt",sep="\t")

r=residuals(object=fit.tmp, type="deviance")
result=RUVr(x=logcpm, cIdx=rownames(logcpm), k=RUVrK, residuals=r, isLog=TRUE) # covariates in result[["W"]]
rownames(result[["W"]])=colnames(result[["normalizedCounts"]]) #rows are not named
#add covariates to y
tmp=result[["W"]][rownames(y$samples), , drop=FALSE]
colnames(tmp)<-"W_RUVr"
y$samples=data.frame(y$samples, tmp)

contrasts=makeContrasts(
    ad3xtg.lithiumorotate0.7vad3xtg.water=groupad3xtg.lithiumorotate0.7,
    levels=design.tmp)

res=glmLRT(glmfit=fit.tmp, contrast=contrasts)
tab=as.data.frame(topTags(object=res, n=Inf, adjust.method="fdr", sort.by="PValue")$table)
empirical <- rownames(y$counts)[which(!(rownames(y$counts) %in% rownames(tab)[1:5000]))]
RUVg_res <- RUVg(y$counts, empirical, k=RUVrK)
rownames(RUVg_res[["W"]])=colnames(RUVg_res[["normalizedCounts"]]) #rows are not named
#add covariates to y
tmp=RUVg_res[["W"]][rownames(y$samples), , drop=FALSE]
colnames(tmp)<-"W_RUVg"
y$samples=data.frame(y$samples, tmp)

design=model.matrix(as.formula("~Batch+W_RUVr+group"), data=y$samples)
colnames(design)=gsub(":", ".", colnames(design))
colnames(design)=gsub("\\(|\\)", "", colnames(design))
colnames(design)=make.names(colnames(design))

y=estimateDisp(y=y, design=design, robust=FALSE)
fit=glmFit(y=y, design=design, robust=FALSE)

contrasts=makeContrasts(
    ad3xtg.lithiumorotate0.7vad3xtg.water=groupad3xtg.lithiumorotate0.7,
    levels=design)

res_RUVr=glmLRT(glmfit=fit, contrast=contrasts)
tab_RUVr=as.data.frame(topTags(object=res_RUVr, n=Inf, adjust.method="fdr", sort.by="PValue")$table)
write.table(tab_RUVr,"res_EdgeR_RUVr_3xTG_LiO_water.txt",sep="\t")

fdr_threshold <- 0.05 # Significance threshold for Adjusted P-value (FDR)
logfc_threshold <- 0 # Significance threshold for log2 Fold Change (e.g., 1.0 means a 2-fold change)
n_top_genes_to_label <- 20 # Number of top significant genes to label on the plot

tab_RUVr$gene_symbol <- rownames(tab_RUVr) # Use row names (often Ensembl IDs)
tab_RUVr$gene_symbol[is.na(tab_RUVr$gene_symbol) | tab_RUVr$gene_symbol == ""] <- rownames(tab_RUVr)[is.na(tab_RUVr$gene_symbol) | tab_RUVr$gene_symbol == ""]
tab_RUVr$Significance <- "Not Significant"
# Upregulated: FDR < threshold and logFC > threshold
tab_RUVr$Significance[tab_RUVr$FDR < fdr_threshold & tab_RUVr$logFC > logfc_threshold] <- "Upregulated"
# Downregulated: FDR < threshold and logFC < -threshold
tab_RUVr$Significance[tab_RUVr$FDR < fdr_threshold & tab_RUVr$logFC < -logfc_threshold] <- "Downregulated"
tab_RUVr$Significance[tab_RUVr$gene_symbol %in% c("Apoe")] <- "Important"
tab_RUVr$Significance <- factor(tab_RUVr$Significance, levels = c("Important","Upregulated", "Downregulated", "Not Significant"))
tab_RUVr_ordered <- tab_RUVr[order(tab_RUVr$FDR, -abs(tab_RUVr$logFC)), ]
top_genes_df<-rbind(head(tab_RUVr_ordered[tab_RUVr_ordered$Significance!="Not Significant",],n_top_genes_to_label),
        tab_RUVr_ordered[tab_RUVr_ordered$Significance != "Not Significant" & abs(tab_RUVr_ordered$logFC)>1,],
        tab_RUVr_ordered[tab_RUVr_ordered$gene_symbol %in% c("Apoe"),])

volcano_plot <- ggplot(tab_RUVr[order(tab_RUVr$Significance,decreasing=T),], aes(x = logFC, y = -log10(FDR))) +
    geom_point(aes(color = Significance)) +
    scale_color_manual(values = c("Important" = "black","Upregulated" = "red","Downregulated" = "blue",
                                  "Not Significant" = "grey")) +
    # Add labels ONLY for the top genes stored in top_genes_df
    geom_text_repel(data = top_genes_df, # Use the filtered data frame for labels
                    aes(label = gene_symbol),
                    size = 3,
                    #box.padding = 0.5,
                    #point.padding = 0.3,
                    max.overlaps = Inf,
                    segment.color = 'grey50',
                    show.legend = FALSE) +
    # Labels and Theme
    labs(title = "Volcano Plot of Differential Gene Expression",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value (FDR)",
         color = "Significance") + # Legend title
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size=9, face="italic"),
      legend.position = "bottom"
    )
print(volcano_plot)
ggsave("volcano_plot_RUVr_Apoe_add_point.pdf", plot = volcano_plot, width = 8, height = 6, dpi = 300)
