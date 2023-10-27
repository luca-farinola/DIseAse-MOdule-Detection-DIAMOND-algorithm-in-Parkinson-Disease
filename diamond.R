options(stringsAsFactors = F)
setwd("C:/Users/lukfa/OneDrive/Desktop/Netbior")

##########################
source("script/getInputData.R")
source("script/createAdjMatrix.R")
source("script/makeConversion.R")
source("script/createCandidateMatrix.R")
source("script/findCandidate.R")
source("script/findDiseaseGenes.R")
##########################
filename_interactome<- "files/interactome.txt"
filename_seed_proteins <- "files/seedsid.txt"

dirRes <- 'Results/'
if( !file.exists(dirRes) ){
  dir.create(dirRes)
}

N_iter <- 500

##########################
input <- getInputData(filename_interactome,filename_seed_proteins)

conversion_table <- input$conversion_table
interactome <- input$interactome
interactome <- interactome[!(interactome$sources == 'literature') &
                             !(interactome$sources == 'literature,signaling') & 
                             !(interactome$sources == 'literature,metabolic') &
                             !(interactome$sources == 'literature,kinase')
                           ,]

seed_proteins <- input$seed_proteins
disease <- 'PD'

##########################
network <- createAdjMatrix(interactome)
adj <- network$adj
nodi <- network$nodi

lapply(disease, function(x){
  
  new_disease_genes <- findDiseaseGenes(x,adj,nodi,seed_proteins,conversion_table,N_iter)
  
  num <- which(disease %in% x)
  
  write.table(new_disease_genes, paste(dirRes,paste(num,"_",x,".txt",sep=""), sep = ""), sep = "\t", 
              row.names = F, col.names = c("gene ID","gene name","p-value","K","Ks"), quote = F)
  
})



