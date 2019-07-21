
setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/data')

file_path <- list.files(pattern='MCT.*.txt')

file_path <- file_path[grep("MCTs05", file_path, invert = T)]
file_path <- file_path[grep("MCTs06", file_path, invert = T)]
file_path <- file_path[grep("MCTs28", file_path, invert = T)]
file_path <- file_path[grep("MCTs29", file_path, invert = T)]

# total <- NULL
# 
# for (i in file_path) {
#   temp <- read.delim(i,sep='\t',header = T)
#   total <- rbind(total,temp)
# }

cazyme_data <- read.table('./Cazyme_total_byperson_qiime.txt',sep='\t',header = T)
rownames(cazyme_data) <- cazyme_data$Cazyme
cazyme_data <- cazyme_data[-1]
rownames(cazyme_data) <- gsub('.*;L4_','',rownames(cazyme_data))


node.target <- c()
node.source <- c()
for (i in 1:ncol(cazyme_data)) {
  for(j in 1:nrow(cazyme_data)){
    if(cazyme_data[j,i]>0){
      node.source <- c(node.source,colnames(cazyme_data)[i])
      node.target <- c(node.target,rownames(cazyme_data)[j])
    }
  }
}

output.edge <- data.frame(source=node.source, target=node.target, source_node_type=rep('MCT',length(node.source)),
                          target_node_type=ifelse(grepl('AA',node.target),'AA',
                                                  ifelse(grepl('CE',node.target),'CE',
                                                         ifelse(grepl('CBM',node.target),'CBM',
                                                                ifelse(grepl('GH',node.target),'GH',
                                                                       ifelse(grepl('GT',node.target),'GT',
                                                                              ifelse(grepl('PL',node.target),'PL','')))))))
write.table(output.edge,"./network/cazy_output_edge.txt", quote=FALSE, sep="\t", row.names=FALSE)

