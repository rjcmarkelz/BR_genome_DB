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



