library(network)
library(ggplot2)



CreateClusterList <- function(d,Conn1,Conn2,p0){
  ClusterList <- NULL


  for (i in 1:d){
    TempCluster <- i
    for(j in 1:d){

      if(i!=j){
        if(Conn2[i,j]>p0){

          TempCluster <- c(TempCluster,j)
        }
      }
    }
    ClusterList <- c(ClusterList, list(TempCluster))
  }

  return(ClusterList)
}


LengthClusterList <- function(d, ClusterList){
  
  
  Length <- numeric(d)
  
  for(i in 1:d){
    
    Length[i] <-length(ClusterList[[i]])
  }
  
  return(Length)
}


FindOneCluster <- function(d, ClusterList){
  ClusterList <- ClusterList
  DifCluster <- NULL
  
  LengthCluster <- LengthClusterList(d,ClusterList)
  
  MaxI <- which.max(LengthCluster)
  
  if(max(LengthCluster)== 1){
    cat("No more cluster")
    return(0)
    
  }
  
  
  
  LengthMaxI <- length(ClusterList[[MaxI]])
  
  for(j in ClusterList[[MaxI]]){
    
    DifCluster[[1]] <- union(DifCluster[[1]],ClusterList[[j]])
    ClusterList[[j]] <- 0
    
    
  }
  
  n <- 1
  
  while( (length(DifCluster[[1]])>LengthMaxI & length(DifCluster[[1]])<d) | n >= 1){
    n <- 0
    
    LengthMaxI <- length(DifCluster[[1]])
    for(j in DifCluster[[1]]){
      if(ClusterList[[j]][1] == 0){
        ### LengthMaxI <- LengthMaxI +1
        
        next
      }
      
      if(ClusterList[[j]][1]  != 0){
        
        cat("j=",j,"\n")
        cat(ClusterList[[j]])
        
        
        
        
        for(k in ClusterList[[j]]){
          
          if(ClusterList[[k]][1]!=0){
            n <- n +1
            DifCluster[[1]] <- union(DifCluster[[1]],ClusterList[[k]] )
            
            ClusterList[[k]] <- 0
          }
        }
      }
    }
  }
  
  if (length(DifCluster[[1]]) == d){
    for (l in 1:d){
      ClusterList[[l]] <- 0
    }
  }
  
  return(list(DifCluster[[1]],ClusterList))
  
}




CreateDifCluster <- function( d, ClusterList){
  
  NoCluster <- 1
  DifCluster <- NULL
  List_FindCluster <- FindOneCluster(d, ClusterList)
  
  while(List_FindCluster[[1]][1]!=0){
    ClusterList <- List_FindCluster[[2]]
    
    DifCluster[[NoCluster]] <- List_FindCluster[[1]]
    
    NoCluster <- NoCluster+1
    
    List_FindCluster <- FindOneCluster(d, ClusterList)
  }
  
  for( j in 1:d){
    
    if(ClusterList[[j]][1] != 0){
      DifCluster[[NoCluster]] <- ClusterList[[j]]
      
      NoCluster <- NoCluster+1
    }
  }
  
  return(DifCluster)
}



CalculateFDA <- function(typeI,conn){
  prob <- as.vector(conn)
  
  return( sum(1-prob[prob>typeI])/length(which(prob>typeI)))
  
  
}


label2mat <- function(clusterlabel){
  d <- length(clusterlabel)
  clustermat <- matrix(0, nrow = d, ncol = d)
  for (i in 1:max(clusterlabel)){
    iInd <- clusterlabel == i
    clustermat[iInd,iInd] <- 1
  }
  return(clustermat)
}



