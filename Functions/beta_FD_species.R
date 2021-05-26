beta.fd.sp <- function(FD.output){
  
  #Species level
  f.beta.sp      = vector("list", length(FD.output$niche.points))
  f.beta.turn.sp = vector("list", length(FD.output$niche.points))
  f.beta.nest.sp = vector("list", length(FD.output$niche.points))
  
  for ( i in 1:length(FD.output$niche.points) ) {
 
   Indexzahl = 1
  
       for (j in 1:length(FD.output$niche.points)){
        sp.a        <- dim(FD.output$niche.points[[i]])[2]
        sp.b        <- dim(FD.output$niche.points[[j]])[2]
        f.beta      <- matrix(nrow = sp.a, ncol = sp.b)
        f.beta.turn <- matrix(nrow = sp.a, ncol = sp.b)
        f.beta.nest <- matrix(nrow = sp.a, ncol = sp.b)
        
        for (k in seq(1, sp.a, 1)){
          for (l in seq(1, sp.b, 1)){
            only.a           <- length(which(FD.output$niche.points[[i]][,k] - FD.output$niche.points[[j]][,l] ==  1 ))
            only.b           <- length(which(FD.output$niche.points[[i]][,k] - FD.output$niche.points[[j]][,l] == -1 ))
            shared.ab        <- length(which(FD.output$niche.points[[i]][,k] + FD.output$niche.points[[j]][,l] ==  2 ))
            f.beta[k,l]      <- (only.a + only.b) / (only.a + only.b + shared.ab)
            f.beta.turn[k,l] <- (2 * min(only.a, only.b)) / (shared.ab + 2 * min(only.a, only.b))
            f.beta.nest[k,l] <- (abs(only.a - only.b) / (shared.ab + only.a + only.b)) * (shared.ab / (shared.ab + 2 * min(only.a, only.b)))
          }  
        }
        
        rownames(f.beta)      <- colnames(FD.output$niche.points[[i]])
        rownames(f.beta.turn) <- colnames(FD.output$niche.points[[i]])
        rownames(f.beta.nest) <- colnames(FD.output$niche.points[[i]])
        colnames(f.beta)      <- colnames(FD.output$niche.points[[j]])
        colnames(f.beta.turn) <- colnames(FD.output$niche.points[[j]])
        colnames(f.beta.nest) <- colnames(FD.output$niche.points[[j]])
        
        f.beta.sp[[i]][[Indexzahl]]      <- f.beta
        f.beta.turn.sp[[i]][[Indexzahl]] <- f.beta.turn
        f.beta.nest.sp[[i]][[Indexzahl]] <- f.beta.nest
        
        Indexzahl = Indexzahl + 1
      }
  }
  
  FD.beta.sp <- list(f.beta.sp      = f.beta.sp, 
                     f.beta.turn.sp = f.beta.turn.sp, 
                     f.beta.nest.sp = f.beta.nest.sp)
  
  return(FD.beta.sp) 
}
