createCandidateMatrix <- function(sp_list,adj){
  
  tmp <- adj[,(colnames(adj) %in% sp_list)]

  ind_zero_row <- which(rowMeans(tmp) == 0)
  if (length(ind_zero_row) > 0 ) tmp <- tmp[-ind_zero_row,]
  
  ind <- which( rownames(tmp) %in% sp_list )
  if (length(ind) > 0 ){
    adj_nn <- tmp[-ind,] 
  } else {
    adj_nn <- tmp
  }
  
  return(adj_nn)
  
}

