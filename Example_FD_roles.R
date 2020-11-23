# Data

# Load networks as list, each network in table format: consumer, resource, interaction strength (if available)
networks <- list(read.table("Data/Network1.txt", header=TRUE),
                 read.table("Data/Network2.txt", header=TRUE)
)

# Load traits as list, each data frame with species names as row names
traits <-  list(read.table("Data/PlantTraitsNW1.txt", header=TRUE, row.names=1),
                read.table("Data/PlantTraitsNW2.txt", header=TRUE, row.names=1)
)


# Functions

source_directory <- paste(getwd(),"/Functions", sep="")

source("Functions/FD_roles.R")            # Main function  (FDbase, FDsum, species' contributions, niche centroids)
source("Functions/beta_FD_networks.R")    # Community level beta diversity (beta FDsum, beta FDbase)
source("Functions/beta_FD_species.R")     # Species level beta diversity (overlap between the PRNs of individual species from all communities)


# Tests

# Functional role diversity
test <-fd.niche(networks, traits, splits=40)
test[[2]]$fd.sum
test[[2]]$fd.base
test[[2]]$contrib.sum
test[[2]]$contrib.base
test[[2]]$niche.centroids

# Beta functional roles (beta FDsum, beta FDbase, each with turnover and nestedness component)
# Jaccard dissimilarity based on overlap in PRNs (i.e. dissimilarity in functional roles) on the network level, with turnover and nestedness
beta.fd.network(test[[2]])

# Jaccard dissimilarity based on overlap in PRNs (i.e. dissimilarity in functional roles) on the species level, with turnover and nestedness
beta.fd.sp(test[[2]])
