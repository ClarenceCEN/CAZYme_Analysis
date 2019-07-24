require(robCompositions)
require(tibble)
require(dplyr)

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

#clr of cazyme data

cazyme <- read.table('./data/Cazyme_total.txt',sep='\t',header = T,row.names = 1)
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

cazyme_filter <- as.data.frame(t(cazyme))
cazyme_filter <- rownames_to_column(cazyme_filter,var = 'X.SampleID')
cazyme_filter <- merge(cazyme_filter,map[c('X.SampleID','UserName','StudyDayNo')],by = 'X.SampleID')
co_vars <- aggregate(cazyme_filter[,!names(cazyme_filter)%in%c('X.SampleID','UserName','StudyDayNo')],by=list(cazyme_filter$UserName),FUN=function(x){sd(x)/mean(x)})
cazyme_sum <- aggregate(cazyme_filter[,!names(cazyme_filter)%in%c('X.SampleID','UserName','StudyDayNo')],by=list(cazyme_filter$UserName),FUN=sum)
cazyme_sum[cazyme_sum>0] <- 1
cazyme_sum_to_keep <- names(cazyme_sum[-1])[colSums(cazyme_sum[-1]) > 34*0.5]
hist(co_vars_cazyme)

##normalization
#cazyme <- sweep(cazyme+0.000001,1,rowSums(cazyme+0.000001),'/')

#cazyme_t <- as.data.frame(t(cazyme*1.0))

#myimpR = impRZilr(cazyme_t, maxit = 3, method = "lm", dl = rep(0.00000000001,ncol(cazyme_t)), verbose = T) # pick a good detection limit...


#cazyme_clr <- cenLR(t(cazyme)+0.000001)$x

#cazyme_filter <- as.data.frame((cazyme_clr))
co_vars[is.na(co_vars)] <- 0
co_vars_cazyme <- colMeans(co_vars[-1])
co_vars_cazyme <- co_vars_cazyme[order(co_vars_cazyme,decreasing = T)]
cazyme_cor_to_keep <- names(co_vars_cazyme)[1:100]
rownames(cazyme_filter) <- cazyme_filter$X.SampleID
cazyme_to_keep <- intersect(cazyme_cor_to_keep,cazyme_sum_to_keep)
cazyme_cleaned <- select(cazyme_filter,cazyme_to_keep)
hist(co_vars_cazyme[cazyme_to_keep])

#cazyme_filter <- sweep(cazyme_filter,2,colSums(cazyme_filter),'/')
myimpR = impRZilr(cazyme_filter*1.0, maxit = 3, method = "lm", dl = rep(1,ncol(cazyme_filter*1.0)), verbose = T) # pick a good detection limit...





cazyme_clr_trans <- t(cazyme_clr)
save(cazyme_clr_trans, file = "data/cazyme_clr.RData")

hist(co_vars_cazyme[cazyme_to_keep])


cazyme_mean <- aggregate(cazyme_cleaned,by=list(cazyme_filter$UserName),FUN=mean)
rownames(cazyme_mean) <- cazyme_mean$Group.1;cazyme_mean <- cazyme_mean[-1]

myimpR = impRZilr(cazyme_mean, maxit = 3, method = "lm", dl = rep(0.000001,ncol(cazyme_mean)), verbose = T) # pick a good detection limit...
cazyme_mean_clr <- cenLR(myimpR$x)$x.clr

rownames(cazyme_mean_clr) <- unique(cazyme_filter$UserName)
cazyme_mean_clr_trans <- as.data.frame(t(cazyme_mean_clr))
save(cazyme_mean_clr_trans, file = "data/cazyme_mean_clr.RData")
