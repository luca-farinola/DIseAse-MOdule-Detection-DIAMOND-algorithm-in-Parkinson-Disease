findCandidate <- function(adj_nn,adj,N,s0,conversion_table){
  
  output <- sapply(rownames(adj_nn), function(y){
    
    k <- sum(adj[y,])
    ks <- sum(adj_nn[y,])
    
    pval <- sum(dhyper(ks:k,s0,(N-s0),k))
    
    return(c(pval,k,ks))
    
    })
                   
  output = t(output)
  colnames(output) <- c("p-value", "k", "ks")
  
  list_candidates <- row.names(output)
  pval <- output[,"p-value"]
  k <- output[,"k"]
  ks <- output[,"ks"]
  
  ind <- which( pval == min(pval))
  
  candidate <- list_candidates[ind]
  pval_candidate <- pval[ind]
  k_candidate <- k[ind]
  ks_candidate <- ks[ind]
  
  symbol <- makeConversion(candidate,conversion_table)
  
  node2add <- data.frame(
    ID = candidate,
    symbol = symbol,
    pval = pval_candidate,
    k = k_candidate,
    ks = ks_candidate)
  
  return(node2add)
  
}
