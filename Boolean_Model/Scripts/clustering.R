######################################
#####Identification of attractors#####
######################################


#####Boolnet#####

install.packages("BoolNet")
library(BoolNet)
library(igraph)

#To download the network move to the respective directory
getwd()
setwd("C:/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th")
getwd()

#To convert the network to an object of class BooleanNetwork
network_mcts <- loadNetwork("network_mcts_c1.txt")

#Attractors are stable cycles of states in a Boolean network. BoolNet is able to identify attractors 
#in synchronous and asynchronous Boolean networks.


########################
###Synchronous search###
########################

#This means that the software starts from all possible states of the network and performs
#synchronous state transitions until a simple or steady-state attractor is reached.
attractors_syn <- getAttractors(network_mcts)
attractors_syn

#Typing attractors_syn calls a special print method that presents the attractor in a human readable way. 
#Here, a state in an attractor is represented by a binary vector, where each entry of the vector codes for one gene. 
#An alternative is to print only the names of the active genes (i.e., the genes that are set to 1) instead of the
#full vector by calling the print() method explicitly with a changed parameter.
active <- print(attractors_syn, activeOnly=TRUE)

#The function getAttractorSequence() can be used to obtain the sequence of states
#that constitute a specific attractor as a table, in this case the attractor 3.
getAttractorSequence(attractors_syn, 3)


#########################
###Asynchronous search###
#########################

#In this search at each point of time t, only one of the transition functions is chosen at random, and the corresponding Boolean
#variable is updated.
attractors_asy <- getAttractors(network_mcts,
                                type="asynchronous",
                                method="random",
                                startStates=10000)
attractors_asy
#conducts an asynchronous search with 10000 random start states on the mcts network.

#With BoolNet is not possible to obtain attractor binary vectors, for that reason it was used BoolNetPerturb.


#####BoolNetPerturb##### 
install.packages("devtools")
install_github("mar-esther23/boolnet-perturb")
library(devtools)
library(BoolNetPerturb)


#################
###Synchronous###
#################

attr <- getAttractors(network_mcts)
#Convert a BoolNet attractor to dataframe.
attr.df <- attractorToDataframe(attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/atractores_syn.csv")


##################
###Asynchronous###
##################

attr <- getAttractors(network_mcts, 
                      type="asynchronous",
                      method="random",
                      startStates=10000)
attr.df <- attractorToDataframe(attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/atractores_asy.csv")

verifySyncronousVsAsyncronous(network_mcts, attr)



########################
#####Edit dataframe#####
########################

#The dataframe obtained before was slightly edited. The column of attractors
#was numeric so it was transformed into a no numeric.  
df <- read.csv("/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/atractores_asy.csv")
df <- select(df, -state, -X)
df$ID <- seq.int(nrow(df)) 
metadata <- select(df, attractor) 
metadata$muestra <- 'Atractor'
metadata$nombre <- paste(metadata$muestra, metadata$attractor) 
metadata <- select(metadata, nombre, attractor) 
df$ID <- paste(metadata$nombre, df$ID) 
metadata$ID <- df$ID 
write.csv(df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/atractores_asy_a.csv", row.names = FALSE)
write.csv(metadata, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/met_asy_a.csv", row.names = FALSE)
#These files were used to do the projections.
