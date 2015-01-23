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
head(transcripts)
tail(transcripts)

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
# succes
class(brassica_db)
?RangedData
?GRanges

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




