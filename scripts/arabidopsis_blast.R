setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
br_ar <- read.table("Brassica_rapa_final.cds_Arabi_blastn", sep = "\t")

head(br_ar, 50)
dim(br_ar)

ar_met <- read.table("Athal_metabolism.csv", sep = ",")
head(ar_met)
ar_met
?sub
br_ar$V2 <- sub("(\\.)(\\d)", "", br_ar$V2)
br_ar_t

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

# gene expression of each gene
# index back into list data structure
# reformat object
# run model



# shade genes
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
read.table(strings)
?read.table
shade_genes <- read.table("athal_shade_genes.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(shade_genes)
dim(shade_genes)
duplicated
dups <- duplicated(shade_genes$AGI)
dups
shade_genes <- shade_genes[!duplicated(shade_genes$AGI),]
dim(shade_genes)


br_sub <- as.data.frame(br_ar[which(br_ar$V2 %in% shade_genes$AGI),])
dim(br_sub)
head(br_sub)

shade_genes <- shade_genes[!duplicated(shade_genes$AGI),]
dim(shade_genes)

# ideally based on fit scores, but just want to take a quick look
br_sub <- br_sub[!(duplicated(br_sub[c("V1","V2")]) | duplicated(br_sub[c("V1","V2")], fromLast = TRUE)), ]
dim(br_sub)
head(br_sub)

# go to data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

# large file (112 MB) takes about 1 minute to get into memory
mapped_counts <- read.delim("RIL_v1.5_mapping.tsv", header = TRUE, sep = "\t")
dim(mapped_counts)
# [1] 43151   843

colnames(mapped_counts)

#replace all NA values with 0 
mapped_counts[is.na(mapped_counts)] <- 0
head(mapped_counts)
tail(mapped_counts)

#remove first row
mapped_counts <- mapped_counts[-1,]
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


# infile shade genes
# join athal blast N table
# remove duplicates
# subset counts file
# VST with DEseq

