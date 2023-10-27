findDiseaseGenes <- function(disease,adj,nodi,seed_proteins,conversion_table,N_iter){
  
  
  ##########################################################
  sp_list <- seed_proteins[disease,]
  sp_list <- sp_list[sp_list != ""]
  sp_list <- sp_list[(sp_list %in% nodi)]
  s0 <- length(sp_list)
  
  adj_nn <- createCandidateMatrix(sp_list,adj)
  ##########################################################
  
  N <- dim(adj)[1]
  
  new_disease_genes <- NULL
  iter = 0
  while (iter < N_iter){
    
    node2add <- findCandidate(adj_nn,adj,N,s0,conversion_table)
    
    new_disease_genes <- rbind(new_disease_genes,node2add)
    
    ##########################################################
    sp_list <- c(sp_list, node2add$ID)
    s0 <- length(sp_list)
    adj_nn <- createCandidateMatrix(sp_list,adj)
    ##########################################################
    
    count = length(node2add$ID)
    iter = iter + count
    print(paste("#:",iter, disease))
    
  }
  
  return(new_disease_genes)
  
}


