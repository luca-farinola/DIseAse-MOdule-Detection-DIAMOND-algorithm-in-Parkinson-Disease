makeConversion <- function(ID,conversion_table){
  
  ind <- which(conversion_table$ID %in% ID)
  
  symbol <- conversion_table$Symbol[ind]
  
  return(symbol)
  
}

