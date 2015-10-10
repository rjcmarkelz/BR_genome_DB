
library(GenomicFeatures)
library(ChIPpeakAnno)

#used the BEDtoGFF.pl script to convert the BED file to GTF
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
?system.file
?makeTranscriptDbFromGFF
?read.table
# test <- read.table("brassica_rapa_final.gtf", header = FALSE)
# head(test)
# str(test)

# update with new genome mapping
br_bed <- read.table("Brassica_rapa_v1.5_final.bed", header = FALSE)
head(br_bed)
str(br_bed)

#format for entry in database
#subset columns required by makeTranscriptDb

#transcripts data frame
transcripts <- br_bed[,c("V4","V4","V1","V6","V2","V3")]
head(transcripts, 10)
tail(transcripts, 20)
dup_transcripts <- as.numeric(duplicated(transcripts$tx_name))
head(dup_transcripts)
sum(dup_transcripts) # 0 


scaf_num <- transcripts[ grep("Scaffold", transcripts$tx_chrom) , ]  
dim(scaf_num)
# 1760 scaffold genes
head(transcripts)

# new genes from upendra
upendra_genes <- transcripts[ grep("Bra1", transcripts$tx_name) , ]  
dim(upendra_genes)
#2732    6
head(upendra_genes)

upendra_genes2 <- transcripts[ grep("Bra0", transcripts$tx_name) , ]  
dim(upendra_genes2)

# # Explore gaps on chromosome 5
# A05     2926071 2914390 2937752
# A05     8139901 8117310 8162493
# A05     9725452 9523975 9926930 R500    R500    I
# A05     18811891 18811483        18812300 
dim(transcripts[transcripts$tx_chrom == "A05" & transcripts$tx_end > 2937752 & transcripts$tx_start < 8117310, ])
head(transcripts[transcripts$tx_chrom == "A05" & transcripts$tx_end > 2937752 & transcripts$tx_start < 8117310, ], 20)
# [1] 875   6

# Explore gaps on chromosome 2
# A02     2434565 2381099 2488032 
# A02     9344913 9334845 9354981 
dim(transcripts[transcripts$tx_chrom == "A02" & transcripts$tx_end > 2488032 & transcripts$tx_start < 9334845, ])
head(transcripts[transcripts$tx_chrom == "A02" & transcripts$tx_end > 2488032 & transcripts$tx_start < 9334845, ], 20)
# [1] 1201    6

# 2015_02_19
# It appears that the gene number reduction from the known marker bins below and the transcript database are caused by 
# genes falling into: telemeric regions, marker gaps, and scaffolds where we have little information
# from ~43000 we drop down to ~25000 that we can conservatively call allele specific expression

#make unique id for each gene based on row number for now 1-44239
transcripts$V4 <- rownames(transcripts)
head(transcripts)

#replace names of new 
colnames(transcripts) <- c("tx_id", "tx_name", "tx_chrom", "tx_strand",
	                         "tx_start", "tx_end")
head(transcripts)
str(transcripts)
transcripts$tx_id <- as.integer(transcripts$tx_id)
transcripts$tx_id

# write table for eQTL cis/trans
write.table(transcripts, "transcripts_eqtl_start_stop_eqtl.csv", sep = ",", col.names = TRUE, row.names = FALSE)

# splicings dataframe
# will not include splice sites to simplify this for demonstration purposes
# splice sites will all be length of gene until I figure out how to parse BED

splicings <- transcripts
head(splicings)

# splicings number is 1 for each gene
splicings$exon_rank <- rep(1, nrow(splicings))
head(splicings)

# rearrange splicings with only mandatory columns
# rename columns
splicings <- splicings[, c(1,7,5,6)]
head(splicings)
colnames(splicings) <- c("tx_id", "exon_rank", "exon_start", "exon_end")
head(splicings)


brassica_db <- makeTranscriptDb(transcripts, splicings)
brassica_db
# success
class(brassica_db)
?RangedData
?GRanges

save(brassica_db, file = "brassica_gene_db.RData")
#example QTL output
# > bayesint(out_em_2, "A10", prob = 0.98, expandtomarkers = TRUE)
#              chr      pos      lod
# A10_16230621 A10 82.90340 13.10773
# A10_16069633 A10 83.30983 14.08732
# A10_16032467 A10 84.54993 12.45285

#simple test for QTL range
gr <- GRanges(seqnames = "A10", ranges = IRanges(start=16032467, end=16230621))
gr

#look for overlaps between QTL support interval and transcripts
transcripts_qtl <- transcriptsByOverlaps(brassica_db, gr)
# GRanges object with 48 ranges and 2 metadata columns:
transcripts_qtl$tx_name

##############
# 2015_22_01 brassica_UV QTL
##############
QTL.7.3 
#              chr      pos      lod
# A07x19348283 A07 64.30808 6.919690
# A07x19521092 A07 65.12285 9.415914
# A07x20152281 A07 67.16263 7.023449
gr7_3 <- GRanges(seqnames = "A07", ranges = IRanges(start=19348283, end=20152281))
gr7_3 
gr7_3_candidates <- transcriptsByOverlaps(brassica_db, gr7_3)
gr7_3_candidates$tx_name
gr7_3_candidates
write.table(gr7_3_candidates$tx_name, "gr7_3_candidates.txt" )

QTL.8.2
#              chr      pos      lod
# A08x20128470 A08 73.08586 4.377882
# cA08.loc77   A08 77.00000 6.608733
# A08x20766708 A08 80.61716 4.991970
gr8_2 <- GRanges(seqnames = "A08", ranges = IRanges(start=20128470, end=20766708))
gr8_2 
gr8_2_candidates <- transcriptsByOverlaps(brassica_db, gr8_2)
gr8_2_candidates$tx_name
gr8_2_candidates
write.table(gr8_2_candidates$tx_name, "gr8_2_candidates.txt" )

QTL.7.3.wider.LOD
#              chr      pos      lod
# A07x18905857 A07 63.49163 4.694500
# A07x19521092 A07 65.12285 7.008322
# A07x20306343 A07 68.38861 4.493748
gr7_3_wide <- GRanges(seqnames = "A07", ranges = IRanges(start=18905857, end=20306343))
gr7_3_wide 
gr7_3_wide_candidates <- transcriptsByOverlaps(brassica_db, gr7_3_wide)
gr7_3_wide_candidates$tx_name
gr7_3_wide_candidates
?write.table
write.table(gr7_3_wide_candidates$tx_name, "gr7_3_wide_candidates.txt" )


?read.table
br_annotation <- read.table("Brassica_rapa_v1.5_final_annotation.txt", header = FALSE, sep = "\t")
head(br_annotation)
tail(br_annotation)
str(br_annotation)

##############
# 2015_02_19 marker ranges
##############
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/input")
snpmap <- read.delim("bin-genotypes_ref1.5_v0.1.1_tab.txt", header = TRUE, sep = "\t")
head(snpmap)[,1:20]
dim(snpmap)

marker_ranges <- GRanges(seqnames = snpmap$chr, IRanges(start=snpmap$bin.start, end=snpmap$bin.end))
head(marker_ranges)
marker_ranges
(marker_ranges)

brassica_db
#subset genes to remove scaffolds in db
binned_genes <- transcriptsByOverlaps(brassica_db, marker_ranges)
head(binned_genes)
tail(binned_genes)
binned_genes
# noticed that binned genes is only 25276 genes long, not 43456 like the brassica_db
head(binned_genes$tx_name)

# get overlap indices
binned_genes2 <- findOverlaps(binned_genes, marker_ranges)
str(binned_genes2)
binned_genes2 

#replace transcript ranges with bin range they fall into
ranges(binned_genes)[queryHits(binned_genes2)] = ranges(marker_ranges)[subjectHits(binned_genes2)]

# double check
brassica_genes
marker_ranges
binned_genes

#print these rows
binned_genes$tx_id
binned_genes$ranges


df <- data.frame(chr         = seqnames(binned_genes),
                 bin.start   = start(binned_genes),
                 bin.end     = end(binned_genes),
                 tx_name     = binned_genes$tx_name)
head(df)
dim(df)
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/input")
write.table(df, file="gene_marker_ranges.csv", sep=",", row.names = FALSE, col.names = TRUE)

# for tiffany app

df <- read.table("gene_marker_ranges.csv", sep=",", header = TRUE)
head(df)

dim(df)
dim(snpmap)
?merge
df_snpmap <- merge(df, snpmap, all.x = TRUE)
head(df_snpmap)
dim(df_snpmap)
write.table(df_snpmap, file="gene_marker_ranges_merge.csv", sep=",", row.names = FALSE, col.names = TRUE)

df_snpmap_red <- df_snpmap[,-c(1:3,5)]
head(df_snpmap_red)
#replace R500 with 0 and IMB211 with 1
str(df_snpmap_red)

df_snpmap_red[] <- lapply(df_snpmap_red, as.character)

df_snpmap_red[df_snpmap_red == "R500"] <- 0
df_snpmap_red[df_snpmap_red == "IMB211"] <- 1
df_snpmap_red[df_snpmap_red == "HET"] <- NA
head(df_snpmap_red)[,1:20]
tail(df_snpmap_red)
# R runs better with a long table compared to a wide one
write.table(df_snpmap_red, file="gene_marker_contrast_matrix_long.csv", sep=",", row.names = TRUE, col.names = TRUE)

#save col names as vector
df_snpmap_names <- df_snpmap_red$tx_name
df_snpmap_red <- df_snpmap_red[,-c(1)]

#transpose for contrast matrix
df_snpmap_red <- as.data.frame(t(df_snpmap_red))
head(df_snpmap_red)[,1:20]
colnames(df_snpmap_red) <- df_snpmap_names
head(df_snpmap_red)[,1:20]


write.table(df_snpmap_red, file="gene_marker_contrast_matrix.csv", sep=",", row.names = TRUE, col.names = TRUE)

##############
# 2015_04_3 branching 
##############
QTL1_range <- GRanges(seqnames = "A06", ranges = IRanges(start=6674642, end=18162101))
QTL1_range 
QTL1_range_candidates <- transcriptsByOverlaps(brassica_db, QTL1_range)
QTL1_range_candidates$tx_name
QTL1_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL1_range_candidates$tx_name, "Branching_QTL1.txt" )

QTL2_range <- GRanges(seqnames = "A09", ranges = IRanges(start=10979146, end=20118570))
QTL2_range 
QTL2_range_candidates <- transcriptsByOverlaps(brassica_db, QTL2_range)
QTL2_range_candidates$tx_name
QTL2_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL2_range_candidates$tx_name, "Branching_QTL2.txt" )

QTL3_range <- GRanges(seqnames = "A10", ranges = IRanges(start=1631569, end=13422306))
QTL3_range 
QTL3_range_candidates <- transcriptsByOverlaps(brassica_db, QTL3_range)
QTL3_range_candidates$tx_name
QTL3_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL3_range_candidates$tx_name, "Branching_QTL3.txt" )

QTL4_range <- GRanges(seqnames = "A07", ranges = IRanges(start=23803155, end=24855082))
QTL4_range 
QTL4_range_candidates <- transcriptsByOverlaps(brassica_db, QTL4_range)
QTL4_range_candidates$tx_name
QTL4_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL4_range_candidates$tx_name, "Branching_QTL4.txt" )

QTL5_range <- GRanges(seqnames = "A07", ranges = IRanges(start=12178600, end=13111909))
QTL5_range 
QTL5_range_candidates <- transcriptsByOverlaps(brassica_db, QTL5_range)
QTL5_range_candidates$tx_name
QTL5_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL5_range_candidates$tx_name, "Branching_QTL5.txt" )

QTL6_range <- GRanges(seqnames = "A08", ranges = IRanges(start=19588369, end=20679206))
QTL6_range 
QTL6_range_candidates <- transcriptsByOverlaps(brassica_db, QTL6_range)
QTL6_range_candidates$tx_name
QTL6_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL6_range_candidates$tx_name, "Branching_QTL6.txt" )

QTL7_range <- GRanges(seqnames = "A04", ranges = IRanges(start=3248299, end=12327483))
QTL7_range 
QTL7_range_candidates <- transcriptsByOverlaps(brassica_db, QTL7_range)
QTL7_range_candidates$tx_name
QTL7_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL7_range_candidates$tx_name, "Branching_QTL7.txt" )

QTL8_range <- GRanges(seqnames = "A06", ranges = IRanges(start=5582169, end=8655497))
QTL8_range 
QTL8_range_candidates <- transcriptsByOverlaps(brassica_db, QTL8_range)
QTL8_range_candidates$tx_name
QTL8_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL8_range_candidates$tx_name, "Branching_QTL8.txt" )

QTL9_range <- GRanges(seqnames = "A03", ranges = IRanges(start=11791866, end=22151769))
QTL9_range 
QTL9_range_candidates <- transcriptsByOverlaps(brassica_db, QTL9_range)
QTL9_range_candidates$tx_name
QTL9_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL9_range_candidates$tx_name, "Branching_QTL9.txt" )

QTL10_range <- GRanges(seqnames = "A06", ranges = IRanges(start=11791866, end=22151769))
QTL10_range 
QTL10_range_candidates <- transcriptsByOverlaps(brassica_db, QTL10_range)
QTL10_range_candidates$tx_name
QTL10_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL10_range_candidates$tx_name, "Branching_QTL10.txt" )

QTL11_range <- GRanges(seqnames = "A03", ranges = IRanges(start=24029011, end=24556705))
QTL11_range 
QTL11_range_candidates <- transcriptsByOverlaps(brassica_db, QTL11_range)
QTL11_range_candidates$tx_name
QTL11_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL11_range_candidates$tx_name, "Branching_QTL11.txt" )

QTL12_range <- GRanges(seqnames = "A03", ranges = IRanges(start=33896102, end=38599148))
QTL12_range 
QTL12_range_candidates <- transcriptsByOverlaps(brassica_db, QTL12_range)
QTL12_range_candidates$tx_name
QTL12_range_candidates
setwd("/Users/Cody_2/git.repos/brassica_branching/")
write.table(QTL12_range_candidates$tx_name, "Branching_QTL12.txt" )





