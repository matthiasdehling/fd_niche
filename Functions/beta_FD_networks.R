beta.fd.network <- function(FD.output){

  #Network level
  
  # beta FDsum
  fd.beta.sum      <- matrix(0, nrow = length(FD.output$niche.points), ncol = length(FD.output$niche.points))
  fd.beta.sum.turn <- matrix(0, nrow = length(FD.output$niche.points), ncol = length(FD.output$niche.points))
  fd.beta.sum.nest <- matrix(0, nrow = length(FD.output$niche.points), ncol = length(FD.output$niche.points))
  
  for ( i in 1:length(FD.output$niche.points) ) {
    for (j in setdiff(1:length(FD.output$niche.points), i)){
          min.ab           <- pmin(rowSums(FD.output$niche.points[[i]]), rowSums(FD.output$niche.points[[j]]))
          shared.ab        <- sum(min.ab)
          only.a           <- sum(rowSums(FD.output$niche.points[[i]]) - min.ab)
          only.b           <- sum(rowSums(FD.output$niche.points[[j]]) - min.ab) 
          
          fd.beta.sum[i,j]      <- (only.a + only.b) / (only.a + only.b + shared.ab)
          fd.beta.sum.turn[i,j] <- (2 * min(only.a, only.b)) / (shared.ab + 2 * min(only.a, only.b))
          fd.beta.sum.nest[i,j] <- (abs(only.a - only.b) / (shared.ab + only.a + only.b)) * (shared.ab / (shared.ab + 2 * min(only.a, only.b)))
    }  
  }


  # beta FDbase
  fd.beta.base      <- matrix(0, nrow = length(FD.output$niche.points), ncol = length(FD.output$niche.points))
  fd.beta.base.turn <- matrix(0, nrow = length(FD.output$niche.points), ncol = length(FD.output$niche.points))
  fd.beta.base.nest <- matrix(0, nrow = length(FD.output$niche.points), ncol = length(FD.output$niche.points))
  
  for ( i in 1:length(FD.output$niche.points) ) {
    for (j in setdiff(1:length(FD.output$niche.points), i)){
      rowsum_i <- rowSums(FD.output$niche.points[[i]])
      rowsum_i[rowsum_i > 1] <- 1
      rowsum_j <- rowSums(FD.output$niche.points[[j]])
      rowsum_j[rowsum_j > 1] <- 1
      
      shared.ab        <- length(which(rowsum_i + rowsum_j == 2)) 
      only.a           <- length(which(rowsum_i - rowsum_j == 1))
      only.b           <- length(which(rowsum_i - rowsum_j == -1))
      
      fd.beta.base[i,j]      <- (only.a + only.b) / (only.a + only.b + shared.ab)
      fd.beta.base.turn[i,j] <- (2 * min(only.a, only.b)) / (shared.ab + 2 * min(only.a, only.b))
      fd.beta.base.nest[i,j] <- (abs(only.a - only.b) / (shared.ab + only.a + only.b)) * (shared.ab / (shared.ab + 2 * min(only.a, only.b)))
    }  
  }
  
  
  FD.beta.networks <- list(fd.beta.sum       = fd.beta.sum, 
                           fd.beta.sum.turn  = fd.beta.sum.turn, 
                           fd.beta.sum.nest  = fd.beta.sum.nest,
                           fd.beta.base      = fd.beta.base, 
                           fd.beta.base.turn = fd.beta.base.turn, 
                           fd.beta.base.nest = fd.beta.base.nest)
  
  return(FD.beta.networks) 
}
