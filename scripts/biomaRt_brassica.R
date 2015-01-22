library("biomaRt")
#look at available marts
listMarts()

plant_genomes <- useMart("plants_mart_24")
listDatasets(plant_genomes)

brapa <- useDataset("brapa_eg_gene", mart = plant_genomes)

# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")
library("GenomicFeatures")

brapa_db <- makeTranscriptDbFromBiomart(biomart = "plants_mart_24",
	           dataset = "brapa_eg_gene")
brapa_db
