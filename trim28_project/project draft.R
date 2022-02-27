library(dplyr)
library(ggplot2)
library(MetBrewer)
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(ggpubr)
library(limma)
library(reshape2)
library(RColorBrewer)
library(corrplot)


#downloaded files from expression atlas

#open files

counts_file <- "trim28_project/E-MTAB-9809-raw-counts.tsv"
coldata_file <- "trim28_project/E-MTAB-9809-experiment-design.tsv"

counts <- read.table(counts_file, sep = '\t',header = T, row.names = 1, stringsAsFactors = TRUE)
metadata <-  read.table(coldata_file, row.names = 1, header = T, sep = '\t', 
                        stringsAsFactors = TRUE)

#keep a 'dictionary of gene ids and names, then remove gene names from counts file

geneNames <- as.data.frame(counts[,"Gene.Name"])
geneNames$Gene.ID <- rownames(counts)
geneNames <- geneNames %>% rename("Gene.Name" = 'counts[, "Gene.Name"]')
counts <- subset( counts, select = -Gene.Name)

# get only informative columns of metadata, sex and group

metadata <- metadata %>% select('Sample.Characteristic.genotype.', 'Sample.Characteristic.sex.') %>%
        rename("Sex" = 'Sample.Characteristic.sex.' ) %>%
        mutate(Group = ifelse(grepl("Trim28", Sample.Characteristic.genotype.), "Trim28KO", "WT"))



metadata$lib_size = colSums(counts)


ggplot(metadata, aes(fill = Group)) +
        geom_col(aes(x = rownames(metadata), y = lib_size)) +
        labs(title = "Barplot of library sizes") +
        #scale_x_discrete(limits = rownames(metadata)) +
        theme_pubr() +
        xlab(NULL) +
        ylab("Library size") +
        scale_fill_manual(values=met.brewer("Isfahan1", 2)) +
        facet_wrap(~Sex, scales = "free") +
        theme(axis.text.x = element_text(angle=90, vjust=0.4))


cpm <- apply(counts, 2, 
             function(x) as.numeric(x)/sum(as.numeric(x)) * 10^6)


colSums(cpm)


keep <- rowMeans(cpm) > 30
table(keep)


counts.filtered = counts[keep,]
dim(counts.filtered)



metadata_filtered_samples <- subset(metadata, row.names(metadata) != ("ERR4873315"))
subset_Samples <- as.vector(rownames(metadata_filtered_samples))
subset_reads <- counts.filtered[, subset_Samples]

# create deseq object comapring Group and Sex

countdata.deseq <- DESeq2::DESeqDataSetFromMatrix(countData = subset_reads, 
                                                  colData = metadata_filtered_samples, 
                                                  design = ~ Group + Sex)



# Determine the size factors to use for normalization
dds <- estimateSizeFactors(countdata.deseq)


# Extract the normalized counts
countdata_normalized_counts_all <- counts(dds, normalized = TRUE)




V <- apply(countdata_normalized_counts_all , 1, var)
selectedGenes <- names(V[order(V, decreasing = T)][1:1000])

countdata.deseq1 <- countdata.deseq[selectedGenes, ]



# Determine the size factors to use for normalization
dds <- estimateSizeFactors(countdata.deseq1)

# Extract the normalized counts
countdata_normalized_counts <- counts(dds, normalized = TRUE)


# Transform the normalized counts 
vsd <- vst(dds, blind=TRUE, nsub = 1000)

# Extract the matrix of transformed counts
vsd_mat <- assay(vsd)



mat <- limma::removeBatchEffect(vsd_mat, vsd$Sex)

assay(vsd) <- mat
counts_batch_corrected <- assay(vsd)


# Compute the correlation values between samples
vsd_cor <- cor(vsd_mat) 

# Plot the heatmap

breaksList <- seq(0.8, 1, by = 0.001)


heatmap.corr <- pheatmap(vsd_cor, annotation = select(metadata, c(Group, Sex)), 
                         #annotation_colors = list(Group = c('none' = "lightblue3", immunized = "red3", 
                          #                                  "immunized_cross" = "salmon")),
                         breaks = breaksList,
                         # cluster_rows = FALSE,
                         #cluster_cols = FALSE,
                         
                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                         main = "Correlations Heatmap")



colnames(vsd_mat) <- paste(metadata$Group, metadata$Sex, sep = "_")


correlationMatrix <- cor(vsd_mat)
correlationMatrix

#https://melyssaminto.com/blog-posts/how%20to%20change%20corrplot()%20color%20limits%20and%20color%20sc%20c5f87712fa7448948f8cf81ad1ecff56

corrplot(correlationMatrix,
         type = "upper",
         title = "\n\n\n",
         is.corr = TRUE,
         order = "hclust",
         col.lim = c(0.3,1),
         addrect = 2,
         #col = colorRampPalette(c("white", "deepskyblue", "blue4"))(100),
         col = colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(50),
         #addCoef.col = TRUE,
         method = "shade",
         tl.col = "black",
         diag = TRUE)

# Plot PCA 
pcaData <- plotPCA(vsd, intgroup = c("Group", "Sex"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"),digits = 1)

ggplot(pcaData, aes(x = PC1, y = PC2, color = Group)) +
        coord_fixed()+
        geom_point(size = 2, aes(shape = Sex)) +
        xlab(paste0("PC1: ", percentVar[1], "%")) + 
        ylab(paste0("PC2: ", percentVar[2], "%")) +
        geom_text(label = rownames(pcaData), size = 4, nudge_y = 0.2, show.legend = FALSE)+
        #scale_shape_manual(values = c(17,16,15))+
        #stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=exp)) +
        labs(color = "Group") +
        guides(color = guide_legend(order = 1),
               shape = guide_legend(order = 2))+
        theme_pubr(legend = "right")


# Run analysis

dds_deseq <- DESeq(countdata.deseq1)

plotDispEsts(dds_deseq)

# Extract the results of the differential expression analysis

dds_res <- results(dds_deseq, 
                   contrast = c("Group", "Trim28KO", 'WT'), 
                   alpha = 0.05, 
                   lfcThreshold = 0.0)


resultsNames(dds_deseq)

summary(dds_res)


# Shrink the log2 fold changes

shrink_res <- lfcShrink(dds=dds_deseq, coef=2, type="apeglm")



# Save results as a data frame
dds_res_all <- data.frame(shrink_res)

# Subset the results to only return the significant genes with p-adjusted values less than 0.05
dds_res_sig <- subset(dds_res_all, padj < 0.05 & abs(log2FoldChange) >= 0.5)




# Generate logical column 
dds_res_all <- data.frame(dds_res) %>% mutate(threshold = (padj < 0.05 & abs(log2FoldChange) >= 0.58)) %>%
        filter(!is.na(padj))



dds_res_all_named <- left_join(rownames_to_column(dds_res_all), geneNames, by=c("rowname" = "Gene.ID"))

# Create the volcano plot
ggplot(dds_res_all_named) + 
        geom_point(aes(x = log2FoldChange, y = -log10(pvalue), color = threshold)) + 
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") + 
        theme(legend.position = "none", 
              plot.title = element_text(size = rel(1.5), hjust = 0.5), 
              axis.title = element_text(size = rel(1.25)))+
        theme_pubr(legend = "none")


ggplot() +
        geom_point(data = dds_res_all_named, mapping = aes(x = log2FoldChange, y = -log10(pvalue)), 
                   col = "grey80", pch = 21, size = 4) + 
        geom_point(data = subset(dds_res_all_named, log2FoldChange < -0.58 & padj < 0.05),
                   aes(log2FoldChange, -log10(pvalue)), color = "red", size = 4)  +
        geom_point(data = subset(dds_res_all_named, log2FoldChange > 0.58 & padj < 0.05),
                   aes(log2FoldChange, -log10(pvalue)), color = "green3", size = 4) +
        #ggtitle("in-vitro -CD3 (red) vs none (green) \nsignificant genes are colored") +
        geom_text(data = subset(dds_res_all_named, abs(log2FoldChange) > 1.5 & padj < 0.05),
                  aes(log2FoldChange, -log10(pvalue), label = Gene.Name), 
                  vjust =1, size = 4)+
        theme_pubr()

# Subset normalized counts to significant genes
sig_norm_counts_dds <- countdata_normalized_counts[rownames(dds_res_sig), ]



# Choose heatmap color palette
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

# Plot heatmap
heatmap.numClusters <- pheatmap(sig_norm_counts_dds, 
                                #color = heat_colors, 
                                cluster_rows = T, 
                                show_rownames = F,
                                #labels_row = T,
                                annotation = select(metadata, c(Group, Sex)), 
                                #annotation_colors = list(Condition = c('none' = "lightblue3", immunized = "red3")),
                                scale = "row")



alpha = 0.05
sigtab = dds_res[which(dds_res$padj < alpha), ]
#order results by pvalues
sigtab = sigtab[order(sigtab$pvalue),]
head(sigtab)


dds_res_all_full <- left_join(rownames_to_column(dds_res_all), geneNames, by=c("rowname" = "Gene.ID")) %>% arrange(padj)

dds_res_all_full[1,]



#normalize with deseq

#PCA

#volcanos


#random forest



