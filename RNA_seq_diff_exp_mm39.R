# Make sure these packages are locally installed on your user account to use them in ARC
library(DESeq2)
library(Rsubread)

setwd("/work/abc123/RNA_seq")

# Step 1: prepare your own reference genome index.
buildindex(basename = "jindex", reference = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.fna")


# Step 2: read mapping
align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR1783798_1.fastq",
      readfile2 = "/work/abc123/RNA_seq/SRR1783798_2.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR1783799_1.fastq",
      readfile2 = "/work/abc123/RNA_seq/SRR1783799_2.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR1783800_1.fastq",
      readfile2 = "/work/abc123/RNA_seq/SRR1783800_2.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR3391140.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR3391141.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR3391143.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR5008775.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR5008776.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR5008777.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR5008783.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR5008784.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR5008786.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR6847855_1.fastq",
      readfile2 = "/work/abc123/RNA_seq/SRR6847855_2.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR6847856_1.fastq",
      readfile2 = "/work/abc123/RNA_seq/SRR6847856_2.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR6847857_1.fastq",
      readfile2 = "/work/abc123/RNA_seq/SRR6847857_2.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR7852754_1.fastq",
      readfile2 = "/work/abc123/RNA_seq/SRR7852754_2.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)

align(index = "jindex",
      readfile1 = "/work/abc123/RNA_seq/SRR7852755_1.fastq",
      readfile2 = "/work/abc123/RNA_seq/SRR7852755_2.fastq",
      type = "rna",
      input_format = "FASTQ",
      output_format = "BAM",
      sortReadsByCoordinates = TRUE,
      useAnnotation = TRUE,
      annot.ext = "/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
      isGTF = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      chrAliases = NULL)


filesbam_2021 <- list.files(pattern = "\\.subread.BAM$")
filesbam_2021

fc_2021_mm39 <- featureCounts(filesbam_2021,
                              
                              # annotation
                              annot.ext="/work/abc123/RNA_seq/GCF_000001635.27_GRCm39_genomic.gtf",
                              isGTFAnnotationFile=TRUE,
                              GTF.featureType="exon",
                              GTF.attrType="gene_id",
                              chrAliases=NULL,
                              
                              # level of summarization
                              useMetaFeatures=TRUE,
                              
                              # overlap between reads and features
                              allowMultiOverlap=FALSE,
                              minOverlap=1,
                              largestOverlap=FALSE,
                              readExtension5=0,
                              readExtension3=0,
                              read2pos=NULL,
                              
                              # multi-mapping reads
                              countMultiMappingReads=FALSE,
                              fraction=FALSE,
                              
                              # read filtering
                              minMQS=0,
                              splitOnly=FALSE,
                              nonSplitOnly=FALSE,
                              primaryOnly=FALSE,
                              ignoreDup=FALSE,
                              
                              # strandness
                              strandSpecific=0,
                              
                              # exon-exon junctions
                              juncCounts=FALSE,
                              genome=NULL,
                              
                              # parameters specific to paired end reads
                              isPairedEnd=TRUE,
                              requireBothEndsMapped=TRUE,
                              checkFragLength=FALSE,
                              minFragLength=50,
                              maxFragLength=600,
                              countChimericFragments=TRUE,
                              autosort=TRUE)


save.image("RNA_seq_diff_exp5-24-21_mm39.RData")

q()


#///////////////////////////////////////////////////////////////////////////////////////////////////

#Start the RNA work flow 

setwd("C:/UTSA/Lab/Coding/RNA-seq_Pipeline")
load("RNA_seq_diff_exp5-24-21_mm39.RData")
save.image("RNA_seq_diff_exp5-24-21_mm39.RData")

# You can check that your fc object was set correctly using head() and specificng to return the information in the counts 
head(fc_2021_mm39$counts)
head(fc_2021_mm39$targets)
ncol(fc_2021_mm39$counts)

BiocManager::install("edgeR")
library(edgeR)

#edgeR analysis
group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,6,7,8,8,8,9,9,9,9,9,10,10,10))
group_1 <- c("Sertoli", "Sertoli", "Sertoli", "Granulosa", "Granulosa", "Granulosa", "PGC Male", "PGC Male", "PGC Male", "PGC Female", "PGC Female", "PGC Female", "Endoderm", "Mesoderm", "Epiblast", "ESC", "ESC", "ESC", "iPSC", "iPSC", "iPSC", "iPSC", "iPSC", "SSC", "SSC", "SSC")
y <- DGEList(counts=fc_2021_mm39$counts,group=group)

#Filter out genes that have no expression
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

#Normalize the library size using TMM normalization
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

plotMD(cpm(y, log=TRUE), column=26)
abline(h=0, col="red", lty=2, lwd=2)

#To perform quasi-likelihood F-tests. The coef value is for each of the groups (so in this workflow it could be 1-10). The FDR value indicates false discovery rate so this should be very low same for the p-value. 
fit <- glmQLFit(y,design, robust = TRUE)
qlf <- glmQLFTest(fit,coef=10)
topTags(qlf)

#Estimating the dispersion
install.packages("statmod")
library(statmod)
plotBCV(y)

####### This wont effect donwstream pipeline it's just to show the outliers after this estimated common dispersion
DISPLAYOutlier <- glmQLFit(y, design, robust=TRUE)
head(DISPLAYOutlier$coefficients)
plotQLDisp(DISPLAYOutlier)
DGEListFitTest <- glmQLFTest(DISPLAYOutlier, coef=10)
topTags(DGEListFitTest)

# Heatmap analysis
install.packages("ellipsis")
library(ellipsis)
install.packages("vctrs")
library(vctrs)
install.packages("tidyverse")
library(tidyverse)
install.packages("magrittr")
library(magrittr)
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
BiocManager::install("TxDb.Mmusculus.UCSC.mm39.refGene")
library(TxDb.Mmusculus.UCSC.mm39.refGene)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
install.packages("pheatmap")
library(pheatmap)
install.packages("RColorBrewer")
library(RColorBrewer)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
require(org.Mm.eg.db)
GeneSymbol <- mapIds(org.Mm.eg.db, keys = rownames(y), keytype="SYMBOL", column="ENTREZID")
class(GeneSymbol)

#Lets take this object and add it to our DGEListgroup_1 object as a data frame and then add on the ENTREZID as a ne column

y$gene <- data.frame(ENTREZID = GeneSymbol)
head(y$gene)
SYMBOL <- rownames(y$gene)
ENTREZID <- y$gene[,1] 
y$gene <- data.frame(Symbol = SYMBOL, ENTREZID = ENTREZID)

gene.length <- transcriptsBy(TxDb.Mmusculus.UCSC.mm39.refGene, by = "gene")
gene.length_2
gene.length_2 <- unlist(gene.length)
gene.length_2$ENTREZID <- names(gene.length_2)

names(gene.length_2) <- gene.length_2$tx_name
gene.length <- relist(gene.length_2, gene.length)
gene.length.df <- as.data.frame(gene.length)
gene.length.df
gene.length.df <- gene.length.df[ -c(1:2) ]
gene.length.df.2 = gene.length.df %>% group_by(ENTREZID) %>% top_n(n = 1, wt = width) %>% distinct( ENTREZID, .keep_all = TRUE)
gene.length.df.2
gene.length.df.2$length = gene.length.df.2$width

y$gene = y$gene %>% left_join(dplyr::select(gene.length.df.2, c("length", "ENTREZID")), by = "ENTREZID")

head(y$gene)
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#*************************************************************************************
# to plot heatmap you have to scale the data using the edgeR cpm() function so the difference in reads will be indicated by log(cpm).
DGEList.cpm.2021 <- cpm(y, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
gene.length.df.2021 = data.frame(ENTREZID = y$gene[,2]) %>%  left_join(dplyr::select(gene.length.df.2, c("length", "ENTREZID")), by = "ENTREZID")

# genes to plot : 
genes.x = c("Hprt", "Gusb", "Actb", "Pgk1", "Gapdh", "Cryaa", "Pgk2", "Sox2", "Pou5f1", "Nanog", "Fut4", "Myc", "Pou2f3", "Klf4", "Tfcp2l1", "Dppa5a", "Lck", "Smad1", "Lin28a", "Utf1", "Eras", "Zfp42", "Hsp90b1", "Nanos3", "Tbx4", "Prdm1", "Cdx2", "Ddx4", "Dazl", "Ddx5", "Sohlh1", "Sohlh2", "Dppa3", "Hey1", "Itgb3", "Dnd1", "Dppa4", "Ifitm3", "Piwil2", "Tbx2", "Eomes", "Sox3", "Ngn3", "Rarg", "Id4", "Nano2", "Sall4", "Gfra1", "Wt1", "Sox9", "Gdnf", "Inha", "Fshr", "Foxo1", "Dnmt1", "Dnmt3a", "Dnmt3b", "Tet1", "Tet2", "Tet3", "Uhrf1", "Tdgf1", "En1", "Gfap", "Isl1", "Lhx1", "Nes", "T", "Gsc", "Lef1", "Meox1", "Tie1", "Gata4", "Foxa2", "Pdx1", "Nodal", "Sox7", "Sox17", "Cdx2", "Fgfr2", "Cga", "Hand1", "5t4")
df.cpm = as.data.frame(DGEList.cpm.2021)
df.cpm$ENTREZID <- y$gene[,2]
df.cpm = df.cpm %>% left_join(dplyr::select(y$gene, c("Symbol",  "ENTREZID")), by = "ENTREZID")
df.cpm.subset = filter(df.cpm, Symbol %in% genes.x )
mex.cpm.subset = as.matrix(df.cpm.subset[, 1:26])
rownames(mex.cpm.subset) = df.cpm.subset$Symbol
colnames(mex.cpm.subset) = group_1
# 
pdf("Heatmap_2021_mm10.pdf", width = 20, height = 18)
pheatmap(mex.cpm.subset, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         border_color = NA,
         cellwidth = NA, 
         cellheight = NA, 
         scale = "none", 
         display_numbers = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         show_rownames = T, 
         show_colnames = T, 
         main = NA)
dev.off()


points <- c(0,1,2,3,4,5,6,7,8,9)
colors = brewer.pal(10, "Paired")

pdf("PCA_2021_mm10.pdf")
plotMDS(y, col=colors[group], pch=points[group])
legend("bottomleft", legend=c("Sertoli", "Granulosa", "PGC Male", "PGC Female", "Endoderm", "Mesoderm", "Epiblast", "ESC", "iPSC", "SSC"), pch=points, col=colors, ncol=2)
dev.off()
#////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#This concludes what I needed to do for RNA seq while at UTSA. There is obviously much more you can do with pathway analysis, etc.. However, this should indicate the ease at which you can select RNA seq experiments that have been published previously and are available on GEO and use that data to quickly determine DEG in different cell types.
#Thank you -jake
q()
