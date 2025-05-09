# 1. Load Required Libraries
#---------------------------
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(AnnotationDbi)
library(biomaRt)
library(ReactomePA)
library(EnhancedVolcano)
library(pheatmap)
library(GSVA)               # For immune signature enrichment
library(GSEABase)           # To handle gene sets
library(msigdbr)            # Immune signatures

#---------------------------
# 2. Load Data
#---------------------------
setwd("D:/rheumatoid/")

rh_count <- read.table("D:/rheumatoid/GSE89408_GEO_count_matrix_rename.txt.gz", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata <- read.csv("D:/rheumatoid/metadata_1.csv", row.names = 1)
colnames(rh_count) <- rownames(metadata)
metadata$disease <- gsub("[^a-zA-Z0-9._]", "_", metadata$disease)
metadata$disease <- factor(metadata$disease)

#---------------------------
# 3. Preprocessing
#---------------------------
keep <- rowSums(rh_count >= 10) >= 2
rh_count <- rh_count[keep, ]
rh_count <- round(rh_count)

#---------------------------
# 4. DESeq2 Analysis
#---------------------------
dds <- DESeqDataSetFromMatrix(countData = rh_count, colData = metadata, design = ~ disease)
vsd <- vst(dds, blind = TRUE)
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef = 2, type = "apeglm")
res_sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1), ]
res_sig <- na.omit(res_sig)
write.csv(as.data.frame(res), file = "DEGs_deseq2.csv")

#---------------------------
# 5. PCA Plot
#---------------------------
pcaData <- plotPCA(vsd, intgroup = "disease", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = disease)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA of RNA-Seq Samples")

#---------------------------
# 6. Volcano Plot
#---------------------------
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'Differential Expression')

#---------------------------
# 7. Heatmap of Top 50 DEGs
#---------------------------
top_genes <- head(rownames(res_sig[order(res_sig$padj), ]), 50)
write.csv(top_genes, file = "top50_genes.csv", row.names = TRUE)
mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))
annotation_col <- as.data.frame(metadata["disease"])
pheatmap(mat,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         main = "Top 50 DEGs Heatmap")

#---------------------------
# 8. Functional Enrichment (GO & KEGG)
#---------------------------
gene_symbols <- rownames(res_sig)
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pvalueCutoff = 0.05,
                readable = TRUE)

ekegg <- enrichKEGG(gene = entrez_ids,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

write.csv(as.data.frame(ego), "GO_Enrichment.csv")
write.csv(as.data.frame(ekegg), "KEGG_Enrichment.csv")

barplot(ego, showCategory = 10, title = "GO Enrichment")
dotplot(ekegg, showCategory = 10, title = "KEGG Enrichment")

#---------------------------
# 9. Immune Checkpoint Gene Expression
#---------------------------
immune_checkpoints <- c("PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1")
checkpoint_expr <- assay(vsd)[immune_checkpoints, ]
pheatmap(checkpoint_expr,
         annotation_col = annotation_col,
         main = "Immune Checkpoint Genes")

#---------------------------
# 10. Immune Signature Enrichment (GSVA)
#---------------------------
# Get hallmark or immune signature gene sets
msig_immune <- msigdbr(species = "Homo sapiens", category = "C7")  # C7 = immunologic signatures
immune_list <- split(msig_immune$gene_symbol, msig_immune$gs_name)

gsva_result <- gsva(as.matrix(assay(vsd)), immune_list, method = "gsva", kcdf = "Poisson")
write.csv(gsva_result, "Immune_Signature_Enrichment_ssGSEA.csv")

# Heatmap of enriched immune signatures (top 30)
top_sets <- head(rownames(gsva_result[order(rowMeans(gsva_result), decreasing = TRUE), ]), 30)
pheatmap(gsva_result[top_sets, ],
         annotation_col = annotation_col,
         main = "Immune Signature Enrichment (ssGSEA)")
