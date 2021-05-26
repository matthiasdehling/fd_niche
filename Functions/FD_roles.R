#####################################################################################################################
# Function to calculate FD from traits and resource use (network data)                                              #
# D.M. Dehling                                                                                                      #
#                                                                                                                   #
# Reference:                                                                                                        #
# Dehling DM, Stouffer DB (2018) Bringing the Eltonian niche into functional diversity. Oikos 127: 1711-1723        #
#                                                                                                                   #
#####################################################################################################################

fd.niche <- function(networks, traits, 
                     splits=30, 
                     dimensions = dim(pre.coord)[2], 
                     do.pcoa = TRUE,
                     distance = "euclidean",
                     calculate.centroid = TRUE
                     ){

#########################################################################
# Data

# Load list of networks in table format, two trophic levels, interaction frequency (can be all 1 if not available)
# rename columns in networks
for (i in 1:length(networks))  { colnames(networks[[i]]) = c("CONSUMER", "RESOURCE", "INTSTR")}

# rename species so that they can be identified unambiguously:
for (i in 1:length(networks)) { 
  networks[[i]]$CONSUMER = paste("nw", i, "_",networks[[i]]$CONSUMER, sep="")
  networks[[i]]$RESOURCE = paste("nw", i, "_",networks[[i]]$RESOURCE, sep="")
}


# Load traits in table format (as a list)
# rename species so that they can be identified unambiguously:
for (i in 1:length(traits)) { 
  rownames(traits[[i]]) = paste("nw", i, "_", rownames(traits[[i]]), sep="")
}
  

#########################################################################
# Libraries and functions

# Load libraries
library(FD)
library(vegan)
library(StatMatch)
library(MASS)
library(bipartite)
library(geometry)      
library(Matrix)

# Load functions (make sure inhull.R is in the same directory as FD_niche.R)
source(paste(source_directory,"/inhull.R", sep=""))
  
  
#########################################################################
# [1] List of niche coordinates per consumer species per site

niche.coord      <- list()  
orig.niche.coord <- list()

# 1.1 Combine trait data from all sites
pre.coord <- do.call("rbind", traits)

# get resource coordinates from species (either from PCoA coordinates or directly from traits)
if(do.pcoa == TRUE){    
  print("Niche coordinates from PCoA coordinates")  

  
  if(distance == "mahalanobis"){  # 1.2 Calculate distance matrix (Mahalanobis distance)
    m <- dim(pre.coord)[1]
    e <- matrix(0,m,m)
    rownames(e) <- rownames(pre.coord)
    colnames(e) <- rownames(pre.coord)

    for (i in 1:m) { e[i,]<-sqrt(mahalanobis(pre.coord, colMeans(pre.coord[i,]),var(pre.coord)))
    }

    # 1.3 get PCoA coordinates
    plant.coord <- dbFD(as.dist(e), print.pco = TRUE)$x.axes
    plant.coord <- plant.coord[,1:dimensions]
    plant.coord <- as.data.frame(plant.coord)

  } else{
  
    pcoa.plant.coord <- pcoa(dist(pre.coord, "euclidean"))$vectors
    pcoa.plant.coord <- pcoa.plant.coord[,1:dimensions]
    plant.coord      <- as.data.frame(pcoa.plant.coord)

    }  
  
# 1.4 create a data.frame with new header, get rid of row names
    plant.coord$RESOURCE = row.names(plant.coord)
    row.names(plant.coord) <- NULL

# 1.5a Create list of coordinates of resource species per consumer per site from PCoA

    for (i in seq(1, length(traits))){

      unord1 <- merge(plant.coord, networks[[i]])
      orig.niche.coord[[i]] <- unord1[order(unord1$CONSUMER),]
    
      unord2 <- merge(plant.coord, networks[[i]])[,-1]
      niche.coord[[i]] <- unord2[order(unord2$CONSUMER),]
    }

} else {
  print("Niche coordinates from traits")  
  
  plant.coord <- pre.coord
  for (i in 1:dimensions){names(plant.coord)[i] <- paste("Axis.",i,sep = "")}
  
  plant.coord$RESOURCE = row.names(plant.coord)
  row.names(plant.coord) <- NULL

# 1.5b Create list of coordinates of resources per consumer directly from traits    
  for (i in seq(1, length(traits))){
    
    unord1 <- merge(plant.coord, networks[[i]])
    orig.niche.coord[[i]] <- unord1[order(unord1$CONSUMER),]
    
    unord2 <- merge(plant.coord, networks[[i]])[,-1]
    niche.coord[[i]] <- unord2[order(unord2$CONSUMER),]
    
  }                        
}
  
print("Niche coordinates calculated")


#########################################################################
# [2]

# Calculate minimum niche among the species with enough interaction partners.
# Then add an artificial niche to the species for which there are not enough interaction partners.
#
# Determine number of dimensions to be used (default = number of traits)
#
# Identify the species that consume more plant species than number of dimensions.
# Based on those, the minimum niche is calculated (the minimum range in any dimension).
# The minimum range of each dimension is added around each species' niche midpoint.

###################
# Identify smallest range per dimension

Axes <- names(niche.coord[[i]])[!names(niche.coord[[i]]) %in% c("CONSUMER", "INTSTR")]

niche.range.sel = NULL
for (i in seq(1, length(niche.coord))){

  # identify the species with more interaction partners than dimensions
  sp.w.enough <- unique(niche.coord[[i]]$CONSUMER)[table(niche.coord[[i]]$CONSUMER) > dimensions]
  # determine the range of each niche dimension for each species of the reduced species pool
  if(length(sp.w.enough) > 0) niche.range.sel = rbind(niche.range.sel, 
                                                      (aggregate(.~CONSUMER, niche.coord[[i]][niche.coord[[i]]$CONSUMER%in%sp.w.enough,], FUN = max)[,Axes] -
                                                       aggregate(.~CONSUMER, niche.coord[[i]][niche.coord[[i]]$CONSUMER%in%sp.w.enough,], FUN = min)[,Axes])
                                                )  

}  

# get the minimum niche range per axis
min.range <- apply(niche.range.sel, 2, min, na.rm = TRUE)



######################  
# Calculate niche centroid

if(calculate.centroid==TRUE){  
  ORD <- vector("list", length = length(networks))
  
  for (i in 1:length(networks)){ 
    ORD[[i]] <- matrix(data = NA, nrow = length(levels(as.factor(networks[[i]]$CONSUMER))), ncol = length(Axes), dimnames = list(levels(as.factor(networks[[i]]$CONSUMER)), Axes))
    
    # Centroid of resources consumed by each consumer species (weighted mean, if available)
    
    for (j in 1:dim(ORD[[i]])[1]){
      
      INTSPP   <- as.matrix(plant.coord[plant.coord$RESOURCE %in% as.vector(networks[[i]][networks[[i]]$CONSUMER==levels(as.factor(networks[[i]]$CONSUMER))[j],2]),Axes])
      
      INTSTR   <- networks[[i]][networks[[i]]$CONSUMER==levels(as.factor(networks[[i]]$CONSUMER))[j],3]
      
      CENTRSPP <- INTSPP * INTSTR
      
      #print(levels(as.factor(networks[[i]]$CONSUMER))[j])
      
      #print(INTSPP)
      #print(INTSTR)
      #print(CENTRSPP)
      #print(CENTRSPP/INTSPP)
      
      if (length(CENTRSPP)>dimensions) CENTRSPP <- colSums(CENTRSPP)
      ORD[[i]][j,] <- CENTRSPP
    }
    
  }
  niche.centroids <- ORD
  
}else{
  niche.centroids <- rep(list(NA), length(networks))
}  


############################
# add a grid of points around the niche centre of each species with length = the minimum range per axis

for (i in seq(1, length(niche.coord))){
    
  for (k in 1:length(unique(niche.coord[[i]]$CONSUMER))){
    niche.coord[[i]] <- rbind(niche.coord[[i]][c(Axes,"CONSUMER")],         # orig. niche coordinates
                              cbind(expand.grid(as.data.frame(rbind(ORD[[i]][k,]+(min.range/2), ORD[[i]][k,]-(min.range/2)))),
                                    CONSUMER = unique(niche.coord[[i]]$CONSUMER)[k]  # add the species names again to the new data
                              )
                        )
  }
}


#########################################################################
# [3] Build the grid
  
# Get max range of plant coordinates in any of the axes
max.plant.coord   <- apply(plant.coord[,-which(names(plant.coord)=="RESOURCE")],2, max)
min.plant.coord   <- apply(plant.coord[,-which(names(plant.coord)=="RESOURCE")],2, min)
range.plant.coord <- max.plant.coord - min.plant.coord

# Get max range of niche coordinates (can be outside of range of plant coordinates!)  
max.niche.coord   <- apply(do.call(rbind, niche.coord)[,1:dimensions], 2, max)
min.niche.coord   <- apply(do.call(rbind, niche.coord)[,1:dimensions], 2, min)
range.niche.coord <- max.niche.coord - min.niche.coord
  
grid.coord <- vector("list", dimensions)
  
if(length(splits)==dimensions){new.splits=splits} else{new.splits = rep(splits, dimensions)} 

new.splits <- round( range.plant.coord / (max(range.plant.coord)/splits) )
  
for (i in 1:dimensions){
  grid.coord[[i]] <- seq(min.plant.coord[i], max.plant.coord[i], length.out=new.splits[i])
}  
  
print(paste(c("Building grid, splits = ", new.splits), collapse=" "))    

# build the grid from the coordinats of each axis
point.grid <- expand.grid(grid.coord)
dim(point.grid)
  
# As a test:
# Only the coordinates of the grid that fall into the convex hull of all resource species
# This should be faster.
  
test       <- matrix(inhull(point.grid, as.matrix(do.call(rbind, niche.coord)[,1:dimensions])), ncol=1)
point.grid <- point.grid[which(test>-1),] 

print("Grid created")
  
  
  
#########################################################################
# [4] Calculate the indices

# Test whether the points of the point.grid fall into the convex hull of species' niches
# Add columns to the matrix that show whether points are outside (-1), inside (1) or on (0) the hull

FD_1 = list()
FD_2 = list( fd.sum          = list(),
             fd.base         = list(), 
             contrib.sum     = list(),
             contrib.base    = list(),
             niche.points    = list(),
             niche.coord     = list(),
             niche.centroids = list(),
             orig.niche.coord= list(),
             row.sum         = list())

for (i in seq(1, length(niche.coord))){
  j=1
  niche.points <- matrix(inhull(point.grid[,1:dimensions], 
                              as.matrix(niche.coord[[i]][niche.coord[[i]]$CONSUMER==unique(niche.coord[[i]]$CONSUMER)[j],1:dimensions])), ncol=1)

  if(length(unique(niche.coord[[i]]$CONSUMER)) > 1){  
  
    for (j in 2:length(unique(niche.coord[[i]]$CONSUMER))){
  
      niche.points <- cbind(niche.points, 
                            inhull(point.grid[,1:dimensions], 
                            as.matrix(niche.coord[[i]][niche.coord[[i]]$CONSUMER==unique(niche.coord[[i]]$CONSUMER)[j],1:dimensions]))
                    )
      print(paste("Network", i, ":  Species", j, "/", length(unique(niche.coord[[i]]$CONSUMER))))
    }

  }  

  # head(niche.points)
  species = unique(niche.coord[[i]]$CONSUMER)
  colnames(niche.points) <- species
  
  # Inhull gives +1 (inside), -1 (outside), and 0 (on the hull). For convenience, it is better to change it to +1 (inside and on), and 0 (outside)
  # There are hardly any (or no) points ON the hull. Convert all zeros into ones:
  niche.points[which(niche.points==0)]<-1

  # Then convert all -1 into into 0.
  niche.points[which(niche.points==-1)]<-0

  # calculate the niche sizes = number of points that fall into the convex hull
  niche.sizes <- matrix(rep(9999,3*length(species)), nrow = 3, 
                        dimnames = list(c("InOn", "Out", "Sum"), paste("bird", as.character(seq(1, length(species), 1)), sep = "")))

  for (m in 1:length(species)){
    niche.sizes[1,m] <- length(which(niche.points[,m]==1)) 
    niche.sizes[2,m] <- length(which(niche.points[,m]==0))
    niche.sizes[3,m] <- niche.sizes[1,m]+niche.sizes[2,m]
  } 

  colnames(niche.sizes) <- species



  # Calculate functional role diversity

  FD.sum  <- length(which(niche.points != 0))         # FD sum
  FD.base <- length(which(rowSums(niche.points) > 0)) # FD.base
  
  # Relative contribution to FDsum
  ctrb.sum <- niche.sizes["InOn",]/sum(niche.sizes["InOn",])
  
  # Contribution to FD.base
  # 1. calculate the relative niche sizes = number of points that fall into the convex hull divided by number of species for which point also falls in the niche
  wghtd.niche.sizes <- rep(9999,length(species))

  for (lm in 1:length(species)){
    wghtd.niche.sizes[lm] <- sum( niche.points[which(niche.points[,lm]=="1"),lm]/
                                      rowSums(matrix(niche.points[which(niche.points[,lm]=="1"),], ncol = length(species))) )
  }

#  sum(wghtd.niche.sizes) # should be identical to FD.base


  # 2. Relative individual contribution to FD.base [= weighted niche size / FD.base]
  ctrb.base <- wghtd.niche.sizes/FD.base
  names(ctrb.base) <- as.vector(species)
#  sum(ctrb.base)                           # =! 1
  

  
  # Compile the results in two ways: per assemblage (FD_1) and per index (FD_2)
  FD_1[[i]] = list(
    fd.sum          = FD.sum,
    fd.base         = FD.base,
    contrib.sum     = ctrb.sum,
    contrib.base    = ctrb.base,
    niche.points    = niche.points,
    niche.coord     = niche.coord[[i]],
    niche.centroids = niche.centroids[[i]],
    orig.niche.coord= orig.niche.coord[[i]],
    row.sum         = rowSums(niche.points)
  )

  FD_2$fd.sum[[i]]          = FD.sum 
  FD_2$fd.base[[i]]         = FD.base
  FD_2$contrib.sum[[i]]     = ctrb.sum 
  FD_2$contrib.base[[i]]    = ctrb.base
  FD_2$niche.points[[i]]    = niche.points
  FD_2$niche.coord[[i]]     = niche.coord[[i]]
  FD_2$niche.centroids[[i]] = niche.centroids[[i]]
  FD_2$orig.niche.coord[[i]]= orig.niche.coord[[i]]
  FD_2$row.sum[[i]]         = rowSums(niche.points)
  
}

return(list(FD_1, FD_2))

}
