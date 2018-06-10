# |-----------------------------------------------------------------------------------|
# | Project:  Study of Diabetic Nephropathy in HK2 Kidney Cells                       |
# | Script:   RNA-seq data analysis and visualization using DEGSeq                    |
# | Author:   Davit Sargsyan                                                          |
# | Created:  05/31/2018                                                              |
# | Modified: 06/10/2018(DS) Completed analysis                                       |
# |-----------------------------------------------------------------------------------|
# sink(file = "tmp/log_hk2_rnaseq_degseq_v1.txt")
date()

# https://bioconductor.org/packages/release/bioc/html/DEGseq.html
# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

# Header----
require(data.table)
require(ggplot2)
require(DEGseq)
require(knitr)
require(ggdendro)

# Treatment legend----
trt.names <- c("LG",
               "HG",
               "MITC")

# Load data----
dt1 <- fread("data/featurecounts.results.human.csv",
             skip = 1)
dt1 <- droplevels(dt1[, c("Geneid",
                          "HG.dedup.bam", 
                          "LG.dedup.bam",
                          "MIC1.dedup.bam")])
colnames(dt1) <- c("gene",
                trt.names)

dt1

# Remove genes with low counts----
summary(dt1[, -1])
tmp <- rowSums(dt1[, -1])
# Remove if total across 3 samples is no more than 10
dt1 <- droplevels(subset(dt1,
                         tmp > 10))
dt1
# 17,147 genes left, down from 26,364 genes

# DEGseq----
# a. (HG - LG)----
DEGexp(geneExpMatrix1 = dt1,
       geneCol1 = 1, 
       expCol1 = 3, 
       groupLabel1 = colnames(dt1)[3],
       
       geneExpMatrix2 = dt1,
       geneCol2 = 1, 
       expCol2 = 2,
       groupLabel2 = colnames(dt1)[2],
       
       foldChange = 2,
       qValue = 0.1,
       thresholdKind = 5, 
       rawCount = TRUE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp")

hg_lg <- fread("tmp/output_score.txt")
hg_lg
hg_lg[hg_lg$`Signature(q-value(Storey et al. 2003) < 0.1)`,]

# Write as CSV----
write.csv(hg_lg,
          file = "tmp/HK2_RNAseq_DEGseq_HG-LG.csv",
          row.names = FALSE)

# MA Plot----
hg_lg[, mu := (log2(value1) + log2(value2))/2]
hg_lg[, diff := log2(value1) - log2(value2)]
hg_lg$color <- NA
hg_lg$color[hg_lg$`q-value(Storey et al. 2003)` < 0.05 & hg_lg$diff > 0] <- "green"
hg_lg$color[hg_lg$`q-value(Storey et al. 2003)` < 0.05 & hg_lg$diff < 0] <- "red"
tmp <- hg_lg[!is.na(color), ]
tmp <- tmp[sample(x = 1:nrow(tmp), 
                  size = nrow(tmp),
                  replace = FALSE), ]

p1 <- ggplot(hg_lg) +
  geom_point(aes(x = mu,
                 y = diff),
             shape = ".") +
  geom_point(data = tmp,
             aes(x = mu,
                 y = diff,
                 fill = color),
             shape = 21,
             size = 2,
             alpha = 0.5) +
  geom_hline(yintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual("Gene Expression Change",
                    values = c("green", 
                               "red"),
                    labels = c("Upregulated",
                               "Downregulated")) +
  scale_x_continuous("Mean") +
  scale_y_continuous("Difference",
                     breaks = seq(from = -5,
                                  to = 5,
                                  by = 1)) +
  ggtitle("HK2 Gene Expression, HG-LG, FDR < 0.05") +
  theme(plot.title = element_text(hjust = 0.5))
p1
  
tiff(filename = "tmp/hk2_rnaseq_DEGseq_HG-LG_maplot.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# b. (MITC - HG)----
DEGexp(geneExpMatrix1 = dt1,
       geneCol1 = 1, 
       expCol1 = 4, 
       groupLabel1 = colnames(dt1)[4],
       
       geneExpMatrix2 = dt1,
       geneCol2 = 1, 
       expCol2 = 3,
       groupLabel2 = colnames(dt1)[3],
       
       foldChange = 2,
       qValue = 0.1,
       thresholdKind = 5, 
       rawCount = TRUE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp")

mitc_hg <- fread("tmp/output_score.txt")
mitc_hg
mitc_hg[mitc_hg$`Signature(q-value(Storey et al. 2003) < 0.1)`, ]

# Write as CSV----
write.csv(mitc_hg,
          file = "tmp/HK2_RNAseq_DEGseq__MITC-HG.csv",
          row.names = FALSE)

# MA Plot----
mitc_hg[, mu := (log2(value1) + log2(value2))/2]
mitc_hg[, diff := log2(value1) - log2(value2)]
mitc_hg$color <- NA
mitc_hg$color[mitc_hg$`q-value(Storey et al. 2003)` < 0.05 & mitc_hg$diff > 0] <- "green"
mitc_hg$color[mitc_hg$`q-value(Storey et al. 2003)` < 0.05 & mitc_hg$diff < 0] <- "red"
tmp <- mitc_hg[!is.na(color), ]
tmp <- tmp[sample(x = 1:nrow(tmp), 
                  size = nrow(tmp),
                  replace = FALSE), ]

p2 <- ggplot(mitc_hg) +
  geom_point(aes(x = mu,
                 y = diff),
             shape = ".") +
  geom_point(data = tmp,
             aes(x = mu,
                 y = diff,
                 fill = color),
             shape = 21,
             size = 2,
             alpha = 0.5) +
  geom_hline(yintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual("Gene Expression Change",
                    values = c("green", 
                               "red"),
                    labels = c("Upregulated",
                               "Downregulated")) +
  scale_x_continuous("Mean") +
  scale_y_continuous("Difference",
                     breaks = seq(from = -5,
                                  to = 5,
                                  by = 1)) +
  ggtitle("HK2 Gene Expression, MITC-HG, FDR < 0.05") +
  theme(plot.title = element_text(hjust = 0.5))
p2

tiff(filename = "tmp/hk2_rnaseq_DEGseq_MITC-HG_maplot.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Venn diagram----
g1 <- hg_lg[`q-value(Storey et al. 2003)` < 0.05 & 
              `log2(Fold_change) normalized` > 1,]$GeneNames
# 197 genes
g2 <- hg_lg[`q-value(Storey et al. 2003)` < 0.05 & 
              `log2(Fold_change) normalized` < -1,]$GeneNames
# 223 genes

g3 <- mitc_hg[`q-value(Storey et al. 2003)` < 0.05 & 
                `log2(Fold_change) normalized` > 1,]$GeneNames
# 284 genes
g4 <- mitc_hg[`q-value(Storey et al. 2003)` < 0.05 & 
                `log2(Fold_change) normalized` < -1,]$GeneNames
# 291 genes

up.dn <- g1[g1 %in% g4]
# 57 genes
dn.up <- g2[g2 %in% g3]
# 86 genes

# Combine and save the lists----
all.genes <- Reduce(f = function(a, b){
  merge(a, b, all = TRUE)
},
x = list(data.table(gene = g1,
                    `log2(HG/LG) > 1` = g1),
         data.table(gene = g2,
                    `log2(HG/LG) < -1` = g2),
         data.table(gene = g3,
                    `log2(MITC/HG) > 1` = g3),
         data.table(gene = g4,
                    `log2(MITC/HG) < -1` = g4)))
all.genes
write.csv(all.genes,
          file = "tmp/HK2_MITC_all_sign_genes.csv",
          row.names = FALSE)
  
# Heatmap----
ll <- unique(c(up.dn,
               dn.up))

t1 <- merge(hg_lg[hg_lg$GeneNames %in% ll,
                  c("GeneNames",
                    "log2(Fold_change) normalized")],
            mitc_hg[mitc_hg$GeneNames %in% ll,
                    c("GeneNames",
                      "log2(Fold_change) normalized")],
            by = "GeneNames")
colnames(t1) <- c("Gene",
                  "HG-LG",
                  "MITC-HG")
t1 <- t1[order(t1$`HG-LG`,
               decreasing = TRUE), ]
t1
write.csv(t1,
          file = "tmp/hk2_mitc_genes_q-0.5_log2-1.csv",
          row.names = FALSE)

ll <- melt.data.table(data = t1,
                      id.vars = 1,
                      measure.vars = 2:3,
                      variable.name = "Comparison",
                      value.name = "Gene Expression Diff")
ll$Comparison <- factor(ll$Comparison,
                        levels = c("MITC-HG",
                                   "HG-LG"))
lvls <- ll[ll$Comparison == "HG-LG", ]
ll$Gene <- factor(ll$Gene,
                  levels = lvls$Gene[order(lvls$`Gene Expression Diff`)])
# Keep all 143 genes for the plot----
gene.keep <- unique(ll$Gene[order(abs(ll$`Gene Expression Diff`))])
ll <- droplevels(subset(ll,
                        Gene %in% gene.keep))
ll

p3 <- ggplot(data = ll) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(Comparison),
                y = Gene,
                fill = `Gene Expression Diff`),
            color = "white") +
  geom_text(data = ll[Comparison == "HG-LG", ],
            aes(x = rep(1.75,
                        nlevels(Gene)),
                y = Gene,
                label = unique(Gene),
                angle = 90 + seq(from = 0,
                            to = 360,
                            length.out = nlevels(Gene))[as.numeric(Gene)]),
            hjust = 0) +
  scale_fill_gradient2(low = "red",
                       high = "green",
                       mid = "grey",
                       midpoint = 0,
                       name = "Gene Expr Diff") +
  scale_x_continuous(limits = c(0,
                                max(as.numeric(ll$Comparison)) + 0.5),
                     expand = c(0, 0)) +
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("Changes in Gene Expression
          Fold-Change > 0.3 and q-Value < 0.5") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p4

tiff(filename = "tmp/hk2_rnaseq_DEGseq_heatmap.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Heatmap with clustering----
dt2 <- data.frame(dt1[gene %in% t1$Gene, ])
rownames(dt2) <- dt2$gene
dt2 <- dt2[, -1]
dt2 <- log2(dt2 + 0.5)
head(dt2)

# Compute distances between colNamess----
sampleDists <- dist(dt2)
sampleDists

# Make dendrogram data----
dhc <- as.dendrogram(hclust(d = sampleDists),
                     horiz = TRUE)
ddata <- dendro_data(dhc, 
                     type = "rectangle")

# Segment data----
dtp1 <- segment(ddata)
head(dtp1)

# Hitmap data----
dtp2 <- melt.data.table(data.table(gene = rownames(dt2),
                                   dt2),
                        id.vars = 1,
                        measure.vars = 2:4,
                        variable.name = "trt")
dtp2$gene <- factor(dtp2$gene,
                        levels = ddata$labels$label)

p3 <- ggplot(data = dtp2) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(trt),
                y = gene, 
                fill = value),
            color = "white") +
  geom_text(data = dtp2[trt == "LG", ],
            aes(x = rep(2.75,
                        nlevels(gene)),
                y = gene,
                angle = 90 + seq(from = 0,
                                 to = 360,
                                 length.out = nlevels(gene))[as.numeric(gene)],
                label = unique(gene)),
            hjust = 0) +
  scale_fill_gradient2(low = "white", 
                       high = "red", 
                       name = "Log2 Expr") +
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("HK2: Selected Genes Clustered by Log2 RNA Expressions") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_segment(data = dtp1,
               aes(x = -sqrt(y) + 0.5,
                   y = x, 
                   xend = -sqrt(yend) + 0.5,
                   yend = xend),
               size = 1) 
print(p3)

tiff(filename = "tmp/hk2_rnaseq_DEGseq_clustered_heatmap.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()

# sessionInfo()
# sink()