install.packages("devtools")
devtools::install_github("nicolewhite/RNeo4j")
library(RNeo4j)

?startGraph

graph <-  startGraph("http://localhost:7474/db/data/")
startGraph

?createNode
addConstraint(graph, "gene", "name")
gene1 <- createNode(graph, "gene", name = "Bra_10001")
gene2 <- createNode(graph, "gene", name = "Bra_10002") 

getConstraint(graph)

addConstraint(graph, "genotype", "name")
genotype3 <- createNode(graph, c("genotype", "environment"), name = "genotype3")
genotype2 <- createNode(graph, "genotype", name = "genotype2")
genotype1 <- createNode(graph, "genotype", name = "genotype1")

qtl3 <- createNode(graph, c("qtl", "environment"), name = "qtl3")
qtl2 <- createNode(graph, "qtl", name = "qtl2")
qtl1 <- createNode(graph, "qtl", name = "qtl1")

clear(graph)
Y
browse(graph)

###########
###########
# TAGS
###########
###########

#species: name:
#genotype name: 
?addConstraint
addConstraint(graph, "genotype", "genotype")
RIL_2 <- createNode(graph, c("genotype","Brassica"), genotype = "RIL_2")
RIL_1 <- createNode(graph, c("genotype","Brassica"), genotype = "RIL_1")
RIL_3 <- createNode(graph, c("genotype","Brassica"), genotype = "RIL_3")

#marker name:
addConstraint(graph, "marker", "name")
mk_1 <- createNode(graph, c("marker","Brassica"), name = "A01x3323536", chr = "A01", cM = 3.323536)
mk_2 <- createNode(graph, c("marker","Brassica"), name = "A02x222925", chr = "A02", cM = 56.48621943)
mk_3 <- createNode(graph, c("marker","Brassica"), name = "A03x9953796", chr = "A03", cM = 54.25056711)


#individual
##investigation
##treatment
#individual
addConstraint(graph, "Observation", "observation")
obs_1 <- createNode(graph, c("Observation"), observation = "1", value = 6, units = "cm")
obs_2 <- createNode(graph, c("Observation"), observation = "2 ", value = 5, units = "cm")
obs_3 <- createNode(graph, c("Observation"), observation = "3", value = 10, units = "cm")

obs_4 <- createNode(graph, c("Observation"), observation = "4", value = 20, units = "cm")
obs_5 <- createNode(graph, c("Observation"), observation = "5", value = 25, units = "cm")
obs_6 <- createNode(graph, c("Observation"), observation = "6", value = 25, units = "cm")

obs_7 <- createNode(graph, c("Observation"), observation = "7", value = 50, units = "ret time")
obs_8 <- createNode(graph, c("Observation"), observation = "8", value = 100, units = "ret time")
obs_9 <- createNode(graph, c("Observation"), observation = "9", value = 34, units = "ret time")

obs_10 <- createNode(graph, c("Observation"), observation = "10", value = 50, units = "ret time")
obs_11 <- createNode(graph, c("Observation"), observation = "11", value = 100, units = "ret time")
obs_12 <- createNode(graph, c("Observation"), observation = "12", value = 34, units = "ret time")

#time varying observations
obs_13 <- createNode(graph, c("Observation"), observation = "13", value = 50, units = "cm")
obs_14 <- createNode(graph, c("Observation"), observation = "14", value = 100, units = "cm")
obs_15 <- createNode(graph, c("Observation"), observation = "15", value = 1500, units = "cm")
obs_13
#investigations
addConstraint(graph, "Investigation", "Investigation")
exp_1 <- createNode(graph, c("Investigation"), Investigation = "edwards")
exp_2 <- createNode(graph, c("Investigation"), Investigation = "brock")
exp_3 <- createNode(graph, c("Investigation"), Investigation = "duchaine")

#treatments
addConstraint(graph, "Treatment", "Treatment")
trt_1 <- createNode(graph, c("Treatment"), Treatment = "Uncrowded")
trt_2 <- createNode(graph, c("Treatment"), Treatment = "Crowded")
trt_3 <- createNode(graph, c("Treatment"), Treatment = "Drought")

#structures
addConstraint(graph, "Structure", "Structure")
leaf <- createNode(graph, c("Structure"), Structure = "leaf")
petiole <- createNode(graph, c("Structure"), Structure = "petiole")
root <- createNode(graph, c("Structure"), Structure = "root")


#traits
leaflength <- createNode(graph, c("leaflength","trait"), trait = "leaflength")
petiolelength <- createNode(graph, c("petiolelength","trait"),
                   trait = "petiolelength")
metab <- createNode(graph, c("metabolite","trait"),
                     name = "3-Hydroxypropyl", URL = "http://tinyurl.com/q8zrhsq")
createRel(leaflength, "IN_STURCTURE", leaf)
createRel(petiolelength, "IN_STURCTURE", petiole)


######
#RELATIONSHIPS
######
#obs 1 to 6
createRel(obs_1, "IN_EXPERIMENT", exp_1)
createRel(obs_2, "IN_EXPERIMENT", exp_1)
createRel(obs_3, "IN_EXPERIMENT", exp_1)
createRel(obs_4, "IN_EXPERIMENT", exp_1)
createRel(obs_5, "IN_EXPERIMENT", exp_1)
createRel(obs_6, "IN_EXPERIMENT", exp_1)

createRel(obs_1, "IN_TREATMENT", trt_1)
createRel(obs_2, "IN_TREATMENT", trt_1)
createRel(obs_3, "IN_TREATMENT", trt_1)
createRel(obs_4, "IN_TREATMENT", trt_2)
createRel(obs_5, "IN_TREATMENT", trt_2)
createRel(obs_6, "IN_TREATMENT", trt_2)

createRel(obs_1, "IN_STRUCTURE", petiole)
createRel(obs_2, "IN_STRUCTURE", petiole)
createRel(obs_3, "IN_STRUCTURE", petiole)
createRel(obs_4, "IN_STRUCTURE", petiole)
createRel(obs_5, "IN_STRUCTURE", petiole)
createRel(obs_6, "IN_STRUCTURE", petiole)

createRel(obs_1, "TRAIT", leaflength)
createRel(obs_2, "TRAIT", leaflength)
createRel(obs_3, "TRAIT", leaflength)
createRel(obs_4, "TRAIT", leaflength)
createRel(obs_5, "TRAIT", leaflength)
createRel(obs_6, "TRAIT", leaflength)

createRel(obs_1, "IS_GENOTYPE", RIL_1)
createRel(obs_2, "IS_GENOTYPE", RIL_2)
createRel(obs_3, "IS_GENOTYPE", RIL_3)
createRel(obs_4, "IS_GENOTYPE", RIL_1)
createRel(obs_5, "IS_GENOTYPE", RIL_2)
createRel(obs_6, "IS_GENOTYPE", RIL_3)

######
#RELATIONSHIPS
######
#obs 7 to 12
createRel(obs_7, "IN_EXPERIMENT", exp_2)
createRel(obs_8, "IN_EXPERIMENT", exp_2)
createRel(obs_9, "IN_EXPERIMENT", exp_2)
createRel(obs_10, "IN_EXPERIMENT", exp_2)
createRel(obs_11, "IN_EXPERIMENT", exp_2)
createRel(obs_12, "IN_EXPERIMENT", exp_2)

createRel(obs_7, "IN_TREATMENT", trt_1)
createRel(obs_8, "IN_TREATMENT", trt_1)
createRel(obs_9, "IN_TREATMENT", trt_1)
createRel(obs_10, "IN_TREATMENT", trt_2)
createRel(obs_11, "IN_TREATMENT", trt_2)
createRel(obs_12, "IN_TREATMENT", trt_2)

createRel(obs_7, "IN_STRUCTURE", leaf)
createRel(obs_8, "IN_STRUCTURE", leaf)
createRel(obs_9, "IN_STRUCTURE", leaf)
createRel(obs_10, "IN_STRUCTURE", leaf)
createRel(obs_11, "IN_STRUCTURE", leaf)
createRel(obs_12, "IN_STRUCTURE", leaf)

createRel(obs_7, "TRAIT", metab)
createRel(obs_8, "TRAIT", metab)
createRel(obs_9, "TRAIT", metab)
createRel(obs_10, "TRAIT", metab)
createRel(obs_11, "TRAIT", metab)
createRel(obs_12, "TRAIT", metab)

createRel(obs_7, "IS_GENOTYPE", RIL_1)
createRel(obs_8, "IS_GENOTYPE", RIL_2)
createRel(obs_9, "IS_GENOTYPE", RIL_3)
createRel(obs_10, "IS_GENOTYPE", RIL_1)
createRel(obs_11, "IS_GENOTYPE", RIL_2)
createRel(obs_12, "IS_GENOTYPE", RIL_3)


#####
#RELATIONSHIPS
######
#obs 13 15
createRel(obs_13, "IN_EXPERIMENT", exp_3)
createRel(obs_14, "IN_EXPERIMENT", exp_3)
createRel(obs_15, "IN_EXPERIMENT", exp_3)

createRel(obs_13, "IN_TREATMENT", trt_3)
createRel(obs_14, "IN_TREATMENT", trt_3)
createRel(obs_15, "IN_TREATMENT", trt_3)

createRel(obs_13, "IN_STRUCTURE", leaf)
createRel(obs_14, "IN_STRUCTURE", leaf)
createRel(obs_15, "IN_STRUCTURE", leaf)

createRel(obs_13, "TRAIT", metab)
createRel(obs_14, "TRAIT", metab)
createRel(obs_15, "TRAIT", metab)

createRel(obs_13, "TIME", obs_14)
createRel(obs_14, "TIME", obs_15)

createRel(obs_13, "IS_GENOTYPE", RIL_1)
createRel(obs_14, "IS_GENOTYPE", RIL_2)
createRel(obs_15, "IS_GENOTYPE", RIL_3)



#genotype combo
rel_geno1 <- createRel(RIL_1, "HAS_GENOTYPE", mk_1, genotype = "AA")
rel_geno2 <- createRel(RIL_1, "HAS_GENOTYPE", mk_2, genotype = "BB")
rel_geno3 <- createRel(RIL_1, "HAS_GENOTYPE", mk_3, genotype = "AA")

rel_geno1_2 <- createRel(RIL_2, "HAS_GENOTYPE", mk_1, genotype = "AA")
rel_geno2_2 <- createRel(RIL_2, "HAS_GENOTYPE", mk_2, genotype = "AA")
rel_geno3_2 <- createRel(RIL_2, "HAS_GENOTYPE", mk_3, genotype = "BB")

rel_geno1_2 <- createRel(RIL_3, "HAS_GENOTYPE", mk_1, genotype = "BB")
rel_geno2_2 <- createRel(RIL_3, "HAS_GENOTYPE", mk_2, genotype = "BB")
rel_geno3_2 <- createRel(RIL_3, "HAS_GENOTYPE", mk_3, genotype = "BB")


QTL_met <- createNode(graph, c("QTL"), name = "QTL_met")
QTL_leaf <- createNode(graph, c("QTL"), name = "QTL_leaf")

createRel(metab, "HAS_QTL", QTL_met)
createRel(QTL_met, "QTL_LOCATION", mk_1, lod = 4.6, R500_allele = 5, IMB_allele = 0.1)
createRel(QTL_met, "QTL_LOCATION", mk_3, lod = 10, R500_allele = 4, IMB_allele = 2)

createRel(leaflength, "HAS_QTL", QTL_leaf)
createRel(QTL_leaf, "QTL_LOCATION", mk_1, lod = 20, R500_allele = 5, IMB_allele = 0.1)

browse(graph)

#time
#leaf
#petiole
#root
#shoot
#vegetative
#reproductive
##name
##value

#qtl

#environment
#time
##datatype:light
##value:
##units:












# okay now for a whole new db from a dataframe
?data.frame
data <- data.frame(Origin = c("SFO", "AUS", "MCI"),
                  FlightNum = c(1, 2, 3),
                  Destination = c("PDX", "MCI", "LGA"))
head(data)
str(data)
data$Origin <- as.character(data$Origin)
data$Destination <- as.character(data$Destination)


query <- "
MERGE (origin:Airport {name:{origin_name}})
MERGE (destination:Airport {name:{dest_name}})
CREATE (origin)<-[:ORIGIN]-(:Flight {number:{flight_num}})-[:DESTINATION]->(destination)
"
query
?newTransaction
newTransaction

trans <- newTransaction(graph)
trans

for (i in 1:nrow(data)) {
  origin_name <- data[i, ]$Origin
  dest_name <- data[i, ]$Destination
  flight_num <- data[i, ]$FlightNum

  appendCypher(trans, 
               query, 
               origin_name = origin_name, 
               dest_name = dest_name, 
               flight_num = flight_num)
}

head(data)

appendCypher(trans, 
            query, 
  origin_name <- data$Origin[1])
origin_name

origin_name = data[1, ]$Origin
str(origin_name)              



?appendCypher
?commit
commit(t)

cypher(graph, "MATCH (o:Airport)<-[:ORIGIN]-(f:Flight)-[:DESTINATION]->(d:Airport)
               RETURN o.name, f.number, d.name")


