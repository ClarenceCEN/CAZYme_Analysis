require(robCompositions)
require(tibble)
require(dplyr)

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

#clr of cazyme data

cazyme <- read.table('./data/Cazyme_total.txt',sep='\t',header = T,row.names = 1)
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

##normalization
#cazyme <- sweep(cazyme,1,rowSums(cazyme),'/')
#cazyme_t <- as.data.frame(t(cazyme*1.0))

#myimpR = impRZilr(cazyme_t, maxit = 3, method = "lm", dl = rep(0.01,ncol(cazyme_t)), verbose = T) # pick a good detection limit...


cazyme_clr <- cenLR(t(cazyme)+0.001)$x

cazyme_filter <- as.data.frame((cazyme_clr))
cazyme_filter <- rownames_to_column(cazyme_filter,var = 'X.SampleID')
cazyme_filter <- merge(cazyme_filter,map[c('X.SampleID','UserName','StudyDayNo')],by = 'X.SampleID')

co_vars <- aggregate(cazyme_filter[,!names(cazyme_filter)%in%c('X.SampleID','UserName','StudyDayNo')],by=list(cazyme_filter$UserName),FUN=function(x){sd(x)/mean(x)})

co_vars_cazyme <- colMeans(co_vars[-1])

cazyme_clr_trans <- t(cazyme_clr)
save(cazyme_clr_trans, file = "data/cazyme_clr.RData")

hist(co_vars_cazyme)

cazyme_mean <- aggregate(t(cazyme),by=list(cazyme_filter$UserName),FUN=mean)
rownames(cazyme_mean) <- cazyme_mean$Group.1;cazyme_mean <- cazyme_mean[-1]

myimpR = impRZilr(cazyme_mean, maxit = 3, method = "lm", dl = rep(0.000001,ncol(cazyme_mean)), verbose = T) # pick a good detection limit...
cazyme_mean_clr <- cenLR(myimpR$x)$x.clr

rownames(cazyme_mean_clr) <- unique(cazyme_filter$UserName)
cazyme_mean_clr_trans <- as.data.frame(t(cazyme_mean_clr))
save(cazyme_mean_clr_trans, file = "data/cazyme_mean_clr.RData")
