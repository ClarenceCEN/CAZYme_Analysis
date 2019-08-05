require(robCompositions)
require(tibble)
require(dplyr)
require(vegan)

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

#clr of cazyme data

cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
#cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

cazyme <- as.data.frame(t(sweep(cazyme,2,colSums(cazyme),'/')))

cazyme_filter <- as.data.frame(cazyme)

cazyme_filter <- rownames_to_column(cazyme_filter,var = 'X.SampleID')
cazyme_filter <- merge(cazyme_filter,map[c('X.SampleID','UserName','StudyDayNo')],by = 'X.SampleID')

cazyme_mean <- aggregate(cazyme,by=list(cazyme_filter$UserName),FUN=mean)
rownames(cazyme_mean) <- cazyme_mean$Group.1;cazyme_mean <- cazyme_mean[-1]
#co_vars <- sapply(cazyme_mean,FUN = function(x){sd(x)/mean(x)})
#co_vars <- co_vars[order(co_vars,decreasing = T)]
#hist(co_vars)

cazyme_mean_for_fil <- cazyme_mean
cazyme_mean_for_fil[cazyme_mean_for_fil>0] <- 1
cazyme_mean <- cazyme_mean[colSums(cazyme_mean_for_fil)>34*0.3]
myimpR = impRZilr(cazyme_mean, maxit = 3, method = "lm", dl = rep(0.00000001,ncol(cazyme_mean)), verbose = T) # pick a good detection limit...

cazyme_mean_clr <- cenLR(myimpR$x)$x
rownames(cazyme_mean_clr) <- rownames(cazyme_mean)
save(cazyme_mean_clr, file = "data/cazyme_mean_clr_try.RData")
save(cazyme_mean_clr, file = "data/cazyme_median_clr.RData")

cazymes_to_keep <- colnames(cazyme_mean)
save(cazymes_to_keep,file='./data/cazymes_to_keep.RData')
#co_vars <- aggregate(cazyme_filter[,!names(cazyme_filter)%in%c('X.SampleID','UserName','StudyDayNo')],by=list(cazyme_filter$UserName),FUN=function(x){sd(x)/mean(x)})
#cazyme_sum <- aggregate(cazyme_filter[,!names(cazyme_filter)%in%c('X.SampleID','UserName','StudyDayNo')],by=list(cazyme_filter$UserName),FUN=sum)
#cazyme_sum[cazyme_sum>0] <- 1
#cazyme_sum_to_keep <- names(cazyme_sum[-1])[colSums(cazyme_sum[-1]) > 34*0.5]


##normalization
#cazyme <- sweep(cazyme+0.000001,1,rowSums(cazyme+0.000001),'/')

#cazyme_t <- as.data.frame(t(cazyme*1.0))

#myimpR = impRZilr(cazyme_t, maxit = 3, method = "lm", dl = rep(0.00000000001,ncol(cazyme_t)), verbose = T) # pick a good detection limit...


#cazyme_clr <- cenLR(t(cazyme)+0.000001)$x

#cazyme_filter <- as.data.frame((cazyme_clr))
co_vars[is.na(co_vars)] <- 0
co_vars_cazyme <- colMeans(co_vars[-1])
co_vars_cazyme <- co_vars_cazyme[order(co_vars_cazyme,decreasing = T)]
cazyme_cor_to_keep <- names(co_vars_cazyme)[1:150]
rownames(cazyme_filter) <- cazyme_filter$X.SampleID
cazyme_to_keep <- intersect(cazyme_cor_to_keep,cazyme_sum_to_keep)
cazyme_cleaned <- select(cazyme_filter,cazyme_to_keep)
cazyme_cleaned <- select(cazyme_filter,rownames(cazyme_filter))
hist(co_vars_cazyme[cazyme_to_keep])


cazyme_cleaned <- sweep(cazyme_cleaned,1,rowSums(cazyme_cleaned),'/')
cazyme_cleaned <- as.data.frame(t(cazyme_filter))
#myimpR = impRZilr(cazyme_cleaned, maxit = 3, method = "lm", dl = rep(0.00000001,ncol(cazyme_cleaned*1.0)), verbose = T) # pick a good detection limit...
cazyme_clr <- cenLR(cazyme_cleaned +0.00000001)$x




cazyme_clr_trans <- as.data.frame(t(cazyme_clr))
save(cazyme_clr_trans, file = "data/cazyme_clr.RData")

hist(co_vars_cazyme[cazyme_to_keep])


cazyme_mean <- aggregate(cazyme_cleaned,by=list(cazyme_filter$UserName),FUN=mean)
rownames(cazyme_mean) <- cazyme_mean$Group.1;cazyme_mean <- cazyme_mean[-1]

cazyme_mean <- sweep(cazyme_mean,2,colSums(cazyme_mean),'/')


myimpR = impRZilr(cazyme_mean, maxit = 3, method = "lm", dl = rep(0.000001,ncol(cazyme_mean)), verbose = T) # pick a good detection limit...
cazyme_mean_clr <- cenLR(myimpR$x)$x.clr
rownames(cazyme_mean_clr) <- unique(cazyme_filter$UserName)


cazyme_mean_clr_trans <- as.data.frame(cazyme_mean_clr)
save(cazyme_mean_clr_trans, file = "data/cazyme_mean_clr.RData")
