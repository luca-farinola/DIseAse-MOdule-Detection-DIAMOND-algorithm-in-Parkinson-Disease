getInputData <- function(filename_interactome,filename_seed_proteins){
  
  interactome <- read.table(filename_interactome, comment.char = "#", header = F, sep = "\t", check.names = F)
  colnames(interactome) <- c("gene_ID_1", "gene_ID_2", "gene_symbol_1", "gene_symbol_2", "sources")
  
  # eliminate self-loop (i.e. diag(adj)=0 and adj has not col/row with all zero values)
  ind <- which(interactome$gene_ID_1 == interactome$gene_ID_2)
  interactome <- interactome[-ind,]
  
  # create conversion table
  tmp1 <- cbind(interactome$gene_ID_1,interactome$gene_symbol_1)
  tmp2 <- cbind(interactome$gene_ID_2,interactome$gene_symbol_2)
  conversion_table = data.frame(unique(rbind(tmp1,tmp2)))
  colnames(conversion_table) <- c("ID", "Symbol")
  
  seed_proteins <- read.table(filename_seed_proteins, colClasses = "character", header = F, sep = "\t", check.names = F)
  
  disease <- seed_proteins[,1]
  seed_proteins <- seed_proteins[,-1]
  rownames(seed_proteins) <- disease
  
  input <- list(
    conversion_table = conversion_table,
    interactome = interactome,
    seed_proteins = seed_proteins,
    disease = disease
  )
  
  return(input)
  
}

