
setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/Cazyme_output/')

file_path <- list.files(pattern='MCT.*.txt')

file_path <- file_path[grep("MCTs05", file_path, invert = T)]
file_path <- file_path[grep("MCTs06", file_path, invert = T)]
file_path <- file_path[grep("MCTs28", file_path, invert = T)]
file_path <- file_path[grep("MCTs29", file_path, invert = T)]

total <- NULL

for (i in file_path) {
  temp <- read.delim(i,sep='\t',header = T)
}
