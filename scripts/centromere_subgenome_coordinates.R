library(ggplot2)
#infile genomic coordinates of genes
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")


transcripts <- read.table("transcripts_eqtl_start_stop_eqtl.csv", sep = ",", header = TRUE)
head(transcripts)
str(transcripts)
dim(transcripts)

centromeres <- read.table("centromere_subgenome_coordinates.csv", sep = ",", header = TRUE)
head(centromeres)
str(centromeres)
dim(centromeres)

trans_cent <- merge(centromeres, transcripts, by.x = "geneID", by.y = "tx_name", all.x = TRUE)
head(trans_cent)
tail(trans_cent)
dim(trans_cent)

plot(trans_cent$Start, trans_cent$tx_start)
# one to one (nearly)

trans_cent$Mbp <- trans_cent$tx_start/
trans_cent$y_hold <- as.numeric(paste(1))
head(trans_cent)
str(trans_cent)

# not finished yet
trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + geom_tile(aes(x = y_hold, group = Block, color = Block, y = Mbp)) +
                        facet_grid(~ tx_chrom) 

trans_cent_plot







