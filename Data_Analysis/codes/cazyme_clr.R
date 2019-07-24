require(robCompositions)
require(tibble)
require(dplyr)

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

#clr of cazyme data

cazyme <- read.table('./data/Cazyme_total.txt',sep='\t',header = T,row.names = 1)
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")


cazyme_clr <- cenLR(cazyme+0.001)$x

cazyme_filter <- as.data.frame(t(cazyme_clr))
cazyme_filter <- rownames_to_column(cazyme_filter,var = 'X.SampleID')
cazyme_filter <- merge(cazyme_filter,map[c('X.SampleID','UserName','StudyDayNo')],by = 'X.SampleID')

co_vars <- aggregate(cazyme_filter[,!names(cazyme_filter)%in%c('X.SampleID','UserName','StudyDayNo')],by=list(cazyme_filter$UserName),FUN=function(x){sd(x)/mean(x)})

co_vars_cazyme <- colMeans(co_vars[-1])

hist(co_vars_cazyme)

