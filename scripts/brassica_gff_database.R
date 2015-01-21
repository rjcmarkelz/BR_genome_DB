library(GenomicFeatures)
library(ChIPpeakAnno)

#used the BEDtoGFF.pl script to convert the BED file to GTF
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
?system.file
?makeTranscriptDbFromGFF
?read.table
test <- read.table("brassica_rapa_final.gtf", header = FALSE)
head(test)
str(test)

test2 <- read.table("brassica_rapa_final.bed", header = FALSE)
head(test2)
str(test2)

#format for entry in database
#subset columns required by makeTranscriptDb

#transcripts data frame
transcripts <- test2[,c("V4","V4","V1","V6","V2","V3")]
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
splicings$exon_rank <- rep(1,nrow(splicings))
head(splicings)

# rearrange splicings with only mandatory columns
# rename columns
splicings <- splicings[,c(1,7,5,6)]
head(splicings)
colnames(splicings) <- c("tx_id", "exon_rank", "exon_start", "exon_end")



brassica_db <- makeTranscriptDb(transcripts, splicings)
brassica_db
# succes
head(brassica_db)

?extractTranscriptSeqs
gene <- extractTranscriptSeqs()


class(brassica_db)
?RangedData

?GRanges
gr <- GRanges(seqnames = "A10",
              ranges = IRanges(start=16032467, end=16230621))

gr
transcriptsByOverlaps(brassica_db,gr)


# > bayesint(out_em_2, "A10", prob = 0.98, expandtomarkers = TRUE)
#              chr      pos      lod
# A10_16230621 A10 82.90340 13.10773
# A10_16069633 A10 83.30983 14.08732
# A10_16032467 A10 84.54993 12.45285



