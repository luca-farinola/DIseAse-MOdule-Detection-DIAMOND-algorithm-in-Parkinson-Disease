createAdjMatrix <- function(interactome){

  source <- as.character(interactome$gene_ID_1)
  target <- as.character(interactome$gene_ID_2)
  
  nodi <- unique(c(source,target))
  
  N <- length(nodi)
  
  adj <- matrix(0,N,N)
  colnames(adj) <- nodi
  row.names(adj) <- nodi
  
  for(i in 1:N){
    
    ind_s <- which(source %in% nodi[i])
    if (length(ind_s) > 0 ){
      j <- which(nodi %in% target[ind_s])
        adj[i,j] = 1
    }
    
    ind_t <- which(target %in% nodi[i])
    if (length(ind_t) > 0 ){
      j <- which(nodi %in% source[ind_t])
        adj[i,j] = 1
    }
    
  }

  network <- list(adj=adj,nodi=nodi)
  
return(network)

}
