############################################################################

#importing libraries 

library(RCy3)
library(GENIE3)
library(igraph)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)


#importing files 

setwd("C:/Users/lukfa/OneDrive/Desktop/Netbior")

interactome <- read.table("files/interactome.txt", header = T, sep = "\t", check.names = F)
interactome <- interactome[!(interactome$sources == 'literature') &
                             !(interactome$sources == 'literature,signaling') & 
                             !(interactome$sources == 'literature,metabolic') &
                             !(interactome$sources == 'literature,kinase')
                             ,]

seeds <- read.table("files/seedsid.txt", header = F, sep = "\t", check.names = F)
seeds <- t(seeds)
seeds <- as.integer(seeds[-1])
seeds <- as.data.frame(seeds)
colnames(seeds) <- "gene ID"
seeds$type <- "seed"
table_diamond <- read.table("Results/500_associated_pruned_diamond.txt", header = T, sep = "\t", check.names = F)
table_diamond <- table_diamond[order(table_diamond$`p-value`), ]
#table_diamond <- head(table_diamond, 10)
table_diamond$type <- "diamond"
diamond_genes <- table_diamond[,c(1,6)]

bestgenesNE <- rbind(seeds, diamond_genes)
colnames(bestgenesNE) <- c('id','type')

#subset interactome 

subsetinteractome <- interactome[interactome$gene_ID_1 %in% bestgenesNE$id,]
subsetinteractome <- subsetinteractome[subsetinteractome$gene_ID_2 %in% bestgenesNE$id,]

#removing self loops 

duplicate_rows <- subsetinteractome$gene_ID_1 == subsetinteractome$gene_ID_2
subsetinteractome <- subsetinteractome[!duplicate_rows, ]

#getting only diamond or seed 

Joint1 <- left_join(subsetinteractome ,bestgenesNE , by = join_by("gene_ID_1" == "id"))
Joint2 <- left_join(subsetinteractome ,bestgenesNE , by = join_by("gene_ID_2" == "id"))
colnames(Joint2) <- c("--", "--","--","gene_symbol_1","--","type") # i want column names to match to perform rbind
Joint <- rbind(Joint1[,c(3,6)],Joint2[,c(4,6)])
duplicated_rows <- duplicated(Joint$gene_symbol_1)
Joint <- Joint[!duplicated_rows, ]

diff_expr <-readr::read_tsv("files/E-GEOD-20168-A-AFFY-33-query-results.tsv")
colnames(diff_expr) <- c("ID", "Name","Element","foldc","pval","stat")

# subsetnetworks 

regulatory <- subsetinteractome[subsetinteractome$sources == 'regulatory',c(3,4,5)]

signaling <- subsetinteractome[subsetinteractome$sources == 'signaling',c(3,4,5)]

complexes <- subsetinteractome[subsetinteractome$sources == 'complexes',c(3,4,5)]

kinase <- subsetinteractome[subsetinteractome$sources == 'kinase',c(3,4,5)]

#network visualizzation in cytoskape

unique(subsetinteractome$sources)

cytoscapePing ()
cytoscapeVersionInfo ()

graph <- igraph::graph_from_data_frame(as.matrix(subsetinteractome[,c(3,4,5)]))

createNetworkFromIgraph(graph,title="50 Iteration")

copyVisualStyle("default","GENIE3")

setVisualStyle("GENIE3")

loadTableData(Joint, data.key.column = "gene_symbol_1", table.key.column = "id", table = "node")

loadTableData(diff_expr, data.key.column = "Name", table.key.column = "id", table = "node")

setNodeColorMapping("foldc", c(-1,0,1),style.name = "GENIE3", default.color = '#ffff99', 
                    colors=paletteColorBrewerRdBu, mapping.type = "c")

setNodeShapeMapping(style.name = "GENIE3",table.column = "type", table.column.values = c("seed", "diamond"), 
                    shapes = c('ellipse','diamond'))

# gene ontology seecting up and down regulated 

union_result <- union(subsetinteractome$gene_symbol_1,subsetinteractome$gene_symbol_2)

up <- diff_expr[diff_expr$foldc > 0,c(2)]

UNION_UP <- intersect(union_result, up$Name)

down <- diff_expr[diff_expr$foldc > 0,c(2)]

UNION_DOWN <- intersect(union_result, down$Name)

enrich <- enrichGO(UNION_UP, keyType = 'SYMBOL' ,ont = 'BP', OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1)

dotplot(enrich , showCategory =10)


