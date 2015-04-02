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
addConstraint(graph, "genotype", "name")
RIL_2 <- createNode(graph, c("genotype","Brassica"), name = "RIL_2", species = "Brassica rapa")
RIL_1 <- createNode(graph, c("genotype","Brassica"), name = "RIL_1", species = "Brassica rapa")
RIL_3 <- createNode(graph, c("genotype","Brassica"), name = "RIL_3", species = "Brassica rapa")

#marker name:
addConstraint(graph, "marker", "name")
mk_1 <- createNode(graph, c("marker","Brassica"), name = "A01x3323536", chr = "A01", cM = 3.323536)
mk_2 <- createNode(graph, c("marker","Brassica"), name = "A02x222925", chr = "A02", cM = 56.48621943)
mk_3 <- createNode(graph, c("marker","Brassica"), name = "A03x9953796", chr = "A03", cM = 54.25056711)


#individual
##investigation
##treatment
#individual
RIL_1_cr <- createNode(graph, c("investigation","treatment"), name = "RIL1_cr",
              genotype = "RIL_1", trt = "cr", exp = "edwards")
RIL_2_cr <- createNode(graph, c("investigation","treatment"), name = "RIL2_cr",
               genotype = "RIL_2", trt = "cr", exp = "edwards")
RIL_3_cr <- createNode(graph, c("investigation","treatment"), name = "RIL3_cr",
                      genotype = "RIL_3", trt = "cr", exp = "edwards")

RIL_1_un <- createNode(graph, c("investigation","treatment"), name = "RIL1_un",
                 genotype = "RIL_1", trt = "un", exp = "edwards")
RIL_2_un <- createNode(graph, c("investigation","treatment"), name = "RIL2_un",
                 genotype = "RIL_2", trt = "un", exp = "edwards")
RIL_3_un <- createNode(graph, c("investigation","treatment"), name = "RIL3_un",
                 genotype = "RIL_3", trt = "un", exp = "edwards")

#trait
met_1_1 <- createNode(graph, c("trait","metabolite","leaf","shoot"),
                   name = "3-Hydroxypropyl")
met_1_2 <- createNode(graph, c("trait","metabolite","root"),
                   name = "3-Hydroxypropyl")
met_1_3 <- createNode(graph, c("metabolite"), name = "3-Hydroxypropyl",
                       URL = "http://tinyurl.com/q8zrhsq")
#delete(met_1_1)

#create some genotype trait relationships
rel_1 <- createRel(RIL_1_cr, "HAS_TRAIT", met_1_1, value = 1.2)
rel_2 <- createRel(RIL_2_cr, "HAS_TRAIT", met_1_1, value = 2)
rel_3 <- createRel(RIL_3_cr, "HAS_TRAIT", met_1_1, value = 10)
rel_4 <- createRel(RIL_1_un, "HAS_TRAIT", met_1_1, value = 1)
rel_5 <- createRel(RIL_2_un, "HAS_TRAIT", met_1_1, value = 5)
rel_6 <- createRel(RIL_3_un, "HAS_TRAIT", met_1_1, value = 6)

rel_2_1 <- createRel(RIL_1_cr, "HAS_TRAIT", met_1_2, value = 1.5)
rel_2_2 <- createRel(RIL_2_cr, "HAS_TRAIT", met_1_2, value = 6)
rel_2_3 <- createRel(RIL_3_cr, "HAS_TRAIT", met_1_2, value = 5)
rel_2_4 <- createRel(RIL_1_un, "HAS_TRAIT", met_1_2, value = 4)
rel_2_5 <- createRel(RIL_2_un, "HAS_TRAIT", met_1_2, value = 1)
rel_2_6 <- createRel(RIL_3_un, "HAS_TRAIT", met_1_2, value = 8)

rel_met_1 <- createRel(met_1_1, "IS_METABOLITE", met_1_3)
rel_met_1 <- createRel(met_1_2, "IS_METABOLITE", met_1_3)

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

createRel(RIL_1_cr, "IS_GENOTYPE", RIL_1)
createRel(RIL_1_un, "IS_GENOTYPE", RIL_1)
createRel(RIL_2_cr, "IS_GENOTYPE", RIL_2)
createRel(RIL_2_cr, "IS_GENOTYPE", RIL_2)
createRel(RIL_3_cr, "IS_GENOTYPE", RIL_3)
createRel(RIL_3_cr, "IS_GENOTYPE", RIL_3)



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


