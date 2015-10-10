setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
br_ar <- read.table("Brassica_rapa_final.cds_Arabi_blastn", sep = "\t")

head(br_ar, 50)
dim(br_ar)

ar_met <- read.table("Athal_metabolism.csv", sep = ",")
head(ar_met)
ar_met
?sub

#remove .1 on AT names
br_ar$V2 <- sub("(\\.)(\\d)", "", br_ar$V2)


br_sub <- as.data.frame(br_ar[which(br_ar$V2 %in% ar_met$value),])
head(br_sub)
dim(br_sub)
unique(br_sub$V1)
?unique

br_ar_mr <- merge(ar_met, br_sub, by.x = "value", by.y = "V2")
head(br_ar_mr)
dim(br_ar_mr)
str(br_ar_mr)
br_ar_mr$value <- as.character(br_ar_mr$value)
br_ar_mr$V1 <- as.character(br_ar_mr$V1)
unique(br_ar_mr$V1)

# go to data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
write.table(br_sub, "br_met_genes.csv", sep = ",", col.names = TRUE, row.names = FALSE)
head(br_sub)


# shade genes
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
read.table(strings)
?read.table
shade_genes <- read.table("athal_shade_genes.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(shade_genes)
dim(shade_genes)

# find duplicates
dups <- duplicated(shade_genes$AGI)
dups
shade_genes <- shade_genes[!duplicated(shade_genes$AGI),]
dim(shade_genes)
# 94 arabidopsis genes

br_sub <- as.data.frame(br_ar[which(br_ar$V2 %in% shade_genes$AGI),])
dim(br_sub) #258
head(br_sub)

# ideally based on fit scores, but just want to take a quick look
# need to figure out fit scores
br_sub <- br_sub[!(duplicated(br_sub[c("V1","V2")]) | duplicated(br_sub[c("V1","V2")], fromLast = TRUE)), ]
dim(br_sub) # 189
head(br_sub)
br_sub <- br_sub[1:2]
br_sub

# go to data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
write.table(br_sub, "br_shade_genes.csv", sep = ",", col.names = TRUE, row.names = FALSE)
head(br_sub)




# go to data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

# large file (112 MB) takes about 1 minute to get into memory
mapped_counts2 <- read.delim("RIL_v1.5_mapping.tsv", header = TRUE, sep = "\t")
dim(mapped_counts2)
dim(mapped_counts)
# [1] 43151   843

colnames(mapped_counts2)
rownames(mapped_counts2)[1:10]

#replace all NA values with 0 
mapped_counts[is.na(mapped_counts)] <- 0
head(mapped_counts)
tail(mapped_counts)

#remove first row
mapped_counts <- mapped_counts[-1,]
mapped_counts2 <- mapped_counts2[-1,]
head(mapped_counts)[,1:10]

#remove these columns from full dataset because either CR or UN are missing
# Br_trtCR:Br_RIL234 columns 340:343
# Br_trtCR:Br_RIL311 columns 556:558
# causing failures to converge below in lmFit()
# need to remove now and redue the design matrix
colnames(mapped_counts)
dim(mapped_counts)
mapped_counts <- mapped_counts[,-c(556:558)]
dim(mapped_counts)
mapped_counts <- mapped_counts[,-c(340:343)]
dim(mapped_counts)
colnames(mapped_counts)

library(DESeq2)
rownames(mapped_counts) <- mapped_counts$mapped_counts
mapped_counts <- mapped_counts[-1]
dim(mapped_counts)
head(mapped_counts[1:10])
head(mapped_counts2[1:10])

samples <- as.data.frame(names(mapped_counts))
colnames(samples) <- paste("samples")
head(samples)
dim(samples)

samples$condition <- sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\4", samples$samples)
head(samples)
samples$type <- paste("single-read")
rownames(samples) <- samples$samples
dim(samples)
# samples <- samples[-1]

# countData <- counts(mapped_counts)
dds <- DESeqDataSetFromMatrix(countData = mapped_counts, colData = samples, design = ~ condition)

vsd <- varianceStabilizingTransformation(dds)
rld <- rlog(dds, fast = TRUE)
rld2 <- rlog(dds, blind = FALSE, fast = TRUE)
?rlog
head(rld2)
str(rld2)
assays <- as.data.frame(assay(rld2))
str(assays)
assays[1]
assays$gene <- mapped_counts2$gene
head(assays)[1:10]

br_genes <- br_sub$V1
vst_shade <- assays[assays$gene %in% br_genes,]
dim(vst_shade)
head(vst_shade)[1:10]
vst_shade <- vst_shade[c(836,1:835)]
# go to data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
?write.table
write.table(vst_shade, "rlog2_shade_brassica_shade.csv", sep = ",", row.names = FALSE)


setwd("/Users/Cody_2/git.repos/brassica_genetic_map/Input")
gene_contrasts <- read.table("gene_marker_contrast_matrix_long.csv", sep = ",", header = TRUE)
head(gene_contrasts)
dim(gene_contrasts)

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
vst_shade <- read.table("rlog2_shade_brassica_shade.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(vst_shade)[1:10]
dim(vst_shade)


setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
gene_genotype <- gene_contrasts[gene_contrasts$tx_name %in% vst_shade$gene,]
dim(gene_genotype)
head(gene_genotype)
write.table(gene_genotype, "shade_gene_cis_genotypes_v1.csv", sep = ",", row.names = FALSE)

# infile shade genes
# join athal blast N table
# remove duplicates
# subset counts file
# VST with DEseq
# dists <- dist(t(assay(rld2)))
# plot(hclust(dists))




