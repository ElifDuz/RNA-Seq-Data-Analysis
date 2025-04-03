#                      
#             ****RNA-Seq data analysis pipeline****
#                ****from raw counts to DEGs****

# Elif Duz
# 28.03.25

# libraries
library(affy)
library(GEOquery)
library(tidyverse)
library(readxl)
library(limma)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
library(biomaRt)
library(tidyr)
library(DESeq2)
library(dplyr)


# option 1: get the normalized data 
# get the data from pre normalized counts
# these type of data were uploaded by analyzer to GEO website or given as supplementary files.
# generally given in EnsembleIDs.
# go to conversion part at line 88 and unique the data

data <- read.delim("GSE222515_output_tpm_all_samples.txt")
#------------------------------------------------------

# option 2: get the count data
# some count files given as zip files
# first unzip the files into counts folder than import them.

gseid= "GSE205457"
folder_name <- paste0(gseid, "data")# ---> create a folder same name with data
untar(paste0(folder_name, "/", gseid, "_RAW.tar"), exdir = folder_name)

counts <- paste0(folder_name, "/counts")
dir.create(counts, showWarnings = FALSE)

list.files(folder_name)
gz_files <- list.files(folder_name, pattern = "\\.out.tab.gz$", full.names = TRUE)

for (file in gz_files) {
  new_file <- gsub(".gz$", "", basename(file))  
  gunzip(file, destname = file.path(counts, new_file), remove = FALSE)
}

#upload the files from folder
foldername <- "GSE205457data/counts/"
files <- list.files(foldername, pattern = "*.tab", full.names = TRUE)
data <- lapply(files, function(samples) {
  df <- read_delim(samples, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 4)
  colnames(df) <- c("Gene", basename(samples))  
  return(df)
})

data <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), data)
#rearrange the colnames
colnames(data) <- sub("_.*", "", colnames(data))

# Deseq2 normalization

metadata= read_excel("merve/GEO Verileri.xlsx", sheet = "GSE205457")
normdata= data[, metadata$Accession]
rownames(normdata)= data$Gene
matr <- data.frame(status=factor(c(rep("Tumor",6),rep("Control",3))),row.names=colnames(normdata))
mat_deseq <- as.matrix(round(normdata))

deseq <- DESeqDataSetFromMatrix(countData = mat_deseq, colData= matr, design= ~status)
deseqsize <- estimateSizeFactors(deseq)
normalized_deseq <- counts(deseqsize, normalized=T)
normdata= data.frame(geneID=rownames(normdata),normalized_deseq )

# if the normalization does not seem proper:
# apply variance stabilization
sizeFactors(deseqsize)
head(cbind(mat_deseq[,1:3], normalized_deseq[,1:3]))

vsd <- vst(deseqsize, blind=FALSE)
normdata <- assay(vsd)

#sometimes EnsembleID's have versions, remove them for conversion

rownames(normdata) <- sub("\\..*", "", rownames(normdata))

#RNA-Seq datasets generally have Ensemble ID's. Convert them into gene symbols.
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr= listAttributes(ensembl)
convIDs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"), 
                 filters = "ensembl_gene_id", 
                 values = rownames(normdata),
                 mart = ensembl)

data= data.frame(names=rownames(normdata),normdata )
data= merge(data, convIDs, by.x= "names", by.y= "ensembl_gene_id", all.x= TRUE)
data2= data.frame(genesym=data$hgnc_symbol, data[,2:(length(data)-2)]) # last two columns are gene biotype and gene symbol

data2$genesym=trimws(toupper(data2$genesym))

#we can sum the counts belongs to same gene symbols
df_unique <- data2 %>%
  group_by(genesym) %>% 
  summarise(across(everything(), sum)) #sum the expression comes from same gene

df_unique <- df_unique[!(is.na(df_unique$genesym) | df_unique$genesym == ""), ] #remove the nameless row
df_unique= unique(df_unique)

# if column names include unnecessary characters ----> check the column names for consistency
# colnames(df_unique)= gsub("^X", "", colnames(df_unique))

# type column represents the disease status, arrange it based on status----> check the disease status column for comparison
df_unique= data.frame(genesym= toupper(df_unique$genesym),df_unique[, metadata[metadata$type2=="Tumor",]$Accession],df_unique[, metadata[metadata$type2=="Control",]$Accession])
df_unique <- df_unique[rowSums(df_unique != 0, na.rm = TRUE) > 0, ] #at least one row (>0) that total expression is different than 0

df_unique <- df_unique %>%
  rowwise() %>% 
  filter(sum(c_across(-1)) >= 9) %>%  #--------> you can also filter the genes based on total counts, we have 9 samples so we used 9 in here
  ungroup()

#colnames(df_unique)= gsub("^X", "", colnames(df_unique)) #---> arrange the column names if necessary

rnames= df_unique$genesym
df_unique$genesym <- NULL
rownames(df_unique) <- rnames

# PCA analysis
sigmatrix= data.frame(status=factor(metadata$type2),rownames=colnames(df_unique))

pcaPlot <- function(data, Key, title){
  require(ggplot2)
  df <- as.matrix(data)
  # remove the NA values
  df[is.infinite(df)] <- NA
  df <- na.omit(df)  
  # remove the rows with 0 in all samples
  df <- df[rowSums(df != 0, na.rm = TRUE) > 0, ]
  # remove the 0 varience
  df <- df[apply(df, 1, var, na.rm = TRUE) > 0, ]
  
  pcaanal <- prcomp(t(df))
  
  pcaplot <- data.frame(PC1 = pcaanal$x[,1], 
                        PC2 = pcaanal$x[,2], 
                        Sample = rownames(pcaanal$x),
                        Key = Key)  
  
  ggplot(data = pcaplot, aes(x = PC1, y = PC2, color = Key)) +
    geom_point(size = 3) +
    geom_text(aes(label = Sample), hjust = 0.6, vjust = 0, size = 4) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_bw() +
    theme(legend.title = element_blank())  
}
pcaPlot(df_unique,sigmatrix$status,paste0(gseid,"_PCA"))
# ----> to check the cumulative proportions:
df=as.matrix(df_unique)
pcaanal=prcomp(t(df))
summary(pcaanal)

# save the PCA plot
png(paste0(gseid, "_PCA.png"), width = 3000, height = 1500, res = 300)
pcaPlot(df_unique,sigmatrix$status,paste0(gseid,"_PCA"))
dev.off()

# save the normalized data
save(metadata, df_unique, file= paste0(gseid, "_normalized_unique.RData"))

# option 
# if any of the samples removed from the dataset, rearrange both data and metadata:
metadata= metadata[metadata$type!= "TMZ2",]
rnames= rownames(df_unique)
df_unique= df_unique[, metadata$Accession]
rownames(df_unique)= rnames
#------------------------------------------------------

# DEG analysis with limma 
# option 1:unpaired samples
sigmatrix= data.frame(status=factor(metadata$type2),rownames=colnames(df_unique))
disease_control_list <- sigmatrix

group <- factor(metadata$type2, levels = c("Tumor", "Control")) #-----> check the levels and rearrange them (disease vs control)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(df_unique, design)
contrast_matrix <- makeContrasts(Tumor_vs_Control = Tumor - Control, levels = design) #----> rearrange this part, disease vs control
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, adjust = "BH", number = Inf)
#------------------------------------------------------

# option 2: paired samples

samples <- data.frame(
  patient <- factor(c(1:3 ,1:3)),  # match the samples first three samples are recurrent, last three samples from primary 
  group <- factor(metadat$type2))

design= model.matrix(~ patient + group, data=samples)

fit <- lmFit(mesenchymaldata, design)
fit <- eBayes(fit)
results <- topTable(fit, coef="grouprecurrent", adjust="BH", number=Inf) #----> rearrange this part (grouprecurrent) recurrent is diseased group 
#------------------------------------------------------

DEGs= result[result$P.Value <0.05,]#----> define p-value for significance, you can choose adjusted or non-adjusted p-values based on number of DEGs
result$FC= 2^(result$logFC)

write.table(result, file= paste0(gseid, "_DEGs.xls"))

results2 <- data.frame(
  gene = rownames(result),      # Gene names
  logFC = result$logFC,         # Log2 fold change
  pvalue = result$P.Value       # P-values
)

results2$negLogP <- -log10(results2$pvalue)

results2$color <- ifelse(results2$logFC > 0.5849625 & results2$pvalue< 0.05, "upregulated", 
                         ifelse(results2$logFC < -0.5849625 & results2$pvalue < 0.05, "downregulated", "not significant"))
top_genes <- results2[order(results2$pvalue), ][1:1000, ] #---> you can select the number of genes shown in graph

# Volcano plot
png(paste0(gseid, "_volcano.png"), width = 2000, height = 1300, res = 300) #-------> check the graph and arrange the size

ggplot(results2, aes(x = logFC, y = negLogP, color = color)) +
  geom_point(alpha = 0.6, size = 3) +  
  scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "not significant" = "grey")) + 
  labs(title = paste0(gseid," Volcano Plot"),
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = c(-0.5849625, 0.5849625), linetype = "dashed", color = "blue") +  
  geom_text(data = top_genes, aes(label = gene), vjust = 1.5, size = 3, check_overlap = TRUE)  

dev.off()

# enrichment analysis
upgenes= na.omit(rownames(result[result$P.Val<0.05 & result$FC>1.5,]))
downgenes= na.omit(rownames(result[result$P.Val<0.05 & result$FC<0.6666,]))

# convert gene symbols to ENSID for enrichment analysis
ENSID= gconvert(as.vector(upgenes))$target
module= gost(as.vector(ENSID), evcodes=TRUE)$result
module2up= module[(module$source != "TF") & (module$source != "HPA") & (module$source != "MIRNA") ,]
GOtermsup= module2up[grep("GO:", module2up$source),]
othersup= module2up[grep("GO:", module2up$source, invert=TRUE),]
moduleup= module2up[,c(9,10,11,3,6,16)]

ENSID= gconvert(as.vector(downgenes))$target
module= gost(as.vector(ENSID), evcodes=TRUE)$result
module2down= module[(module$source != "TF") & (module$source != "HPA") & (module$source != "MIRNA") ,]
GOtermsdown= module2down[grep("GO:", module2down$source),]
othersdown= module2down[grep("GO:", module2down$source, invert=TRUE),]
moduledown= module2down[,c(9,10,11,3,6,16)]

moduleup$sign= "upregulated"
moduledown$sign= "downregulated"
resultenrichment= rbind(moduleup, moduledown)
write.table(resultenrichment, paste0( gseid,"_enrichment", ".xls"))

# if you want to check the expression value of selected genes in disease and control grupes:
# for two conditions:

genes = data.frame(genes=c("SYTL2", "SPINK6", "MMP11", "ESRRA","TPRG1", "NSMF")) #---> get the genes from excel table or just write their name
genes= unique(genes$genes)
df_unique= data.frame(genesym= rownames(df_unique), df_unique)
selected= df_unique[df_unique$genesym %in% genes, ]
metadata_selected= metadata[, c("Accession","type2")] #------> Accession and type

data_long <- selected %>%
  pivot_longer(cols = -genesym, names_to = "Sample", values_to = "Expression")

data_long <- data_long %>%
  left_join(metadata_selected, by = c("Sample" = "Accession"))

colnames(data_long)[4]= "type2" #--------> one of the column represents the disease status

# if you want to rearrange the disease status names
data_long$type2[which(data_long$type2=="TMZ1")]="TMZ2"

ggplot(data_long, aes(x = type2, y = Expression, fill = type2)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  stat_compare_means(comparisons = list(c("Control", "Tumor")), #-------> change the status names
                     method = "t.test", 
                     label = "p.signif", 
                     tip.length = 0.05) +
  facet_wrap(~ genesym, scales = "free_y") +
  theme_minimal() +
  theme(strip.text = element_text(size = 15, color = "black")) + # the size and the color of the gene symbols
  labs(title = gseid,
       x = "Group",
       y = "Expression Level")



