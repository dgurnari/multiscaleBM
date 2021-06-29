setwd("~/GitHub//multiscaleBM//")

library(data.table) # to read csv faster
library(Rcpp)

sourceCpp('R/BallMapper.cpp')

#######################################
MultiScaleBallMapperGroupActionCpp <- function( points , values , epsilon , orbit )
{
  output <- SimplifiedMultiScaleBallMapperCppInterfaceGroupActionAndSparseRepresentation( points , values , epsilon , orbit )
  return_list <- output
}#BallMapperCpp


storeMultiScaleBallMapperGraphInFile <- function( outputFromBallMapper , filename = "BM_graph" )
{
  #Writing vertices
  #fwrite( as.matrix(outputFromBallMapper$vertices) , file=paste(filename,"_vertices",sep=""), col.names = F, row.names = F)
  
  #Writing landmarks
  #fwrite( as.matrix(outputFromBallMapper$landmarks), file=paste(filename,"_landmarks",sep=""), col.names = F, row.names = F)
  
  #Writing edges
  fwrite(matrix(unlist(outputFromBallMapper$edges), nrow=length(outputFromBallMapper$edges), byrow = T), 
         file=paste(filename,"_edges",sep=""), col.names = F, row.names = F)
  
  #Writing points_in_order_from_landmarks
  #fwrite(matrix(unlist(outputFromBallMapper$points_in_order_from_landmarks), nrow=length(outputFromBallMapper$points_in_order_from_landmarks), byrow = T) ,
  #       file=paste(filename,"_points_in_order_from_landmarks",sep=""), col.names = F, row.names = F)  
  
  #Writing distance_of_points_in_order_from_landmarks
  fwrite(matrix(unlist(outputFromBallMapper$distance_of_points_in_order_from_landmarks), nrow=length(outputFromBallMapper$distance_of_points_in_order_from_landmarks), byrow = T) ,
         file=paste(filename,"_distance_of_points_in_order_from_landmarks",sep=""), col.names = F, row.names = F)
  
  #Writing coloration_in_order_from_landmarks
  fwrite(matrix(unlist(outputFromBallMapper$coloration_in_order_from_landmarks), nrow=length(outputFromBallMapper$coloration_in_order_from_landmarks), byrow = T) ,
         file=paste(filename,"_coloration_in_order_from_landmarks",sep=""), col.names = F, row.names = F)
  
}#storeMultiScaleBallMapperGraphInFile



#######################################
# SMALL EXAMPLE
pts <- rbind( c(0,0), c(1,0), c(10,0), c(12,0), c(30,0), c(33,0) )
values = as.data.frame( c(1,2,3,4,5,6) )
# identity orbit, each element is a singleton. 
orbit = as.data.frame( c(1,2,3,4,5,6) )   

epsilon = 1
bm <- MultiScaleBallMapperGroupActionCpp( pts , values , epsilon , orbit )

storeMultiScaleBallMapperGraphInFile(bm, filename = paste0("output/test_", epsilon))



#######################################
# KNOTS EXAMPLE
jones_13n_MIRRORS <- fread('data/Jones_fromK_upto_13n_MIRRORS.csv', sep=',', header = TRUE)
sapply(jones_13n_MIRRORS[, 1:10], class)

# the table has 57 columns, the same structure as before, but twice the rows
colors <- jones_13n_MIRRORS[, 1:5]
coeff <- jones_13n_MIRRORS[, 6:40]
# add the norm of the coefficients
colors$norm <- wordspace::rowNorms(as.matrix(coeff))

# we will again use jones_13_MIRRORS but we need also this file
# this files links each point (rows in jones_13_MIRRORS) with its mirror 
orbit <- read.csv('data/Jones_fromK_upto_13n_MIRRORS_orbs.csv',header=FALSE)

# we can then compute a Symmetric BM
# create a bm of radius epsilon, and color by the signature 
epsilon <- 20


jones_bm <- MultiScaleBallMapperGroupActionCpp( coeff , colors$signature , epsilon , orbit )

storeMultiScaleBallMapperGraphInFile(jones_bm, filename = paste0("output/jones13n_", epsilon))

