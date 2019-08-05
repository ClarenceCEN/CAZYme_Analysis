require(robCompositions)
require(tibble)
require(dplyr)
require(vegan)

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

#clr of cazyme data
cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
#cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")


cazyme_filter <- as.data.frame(t(cazyme))
cazyme_filter <- rownames_to_column(cazyme_filter,var = 'X.SampleID')
cazyme_filter <- merge(cazyme_filter,map[c('X.SampleID','UserName','StudyDayNo')],by = 'X.SampleID')

cazyme_sum <- aggregate(t(cazyme),by=list(cazyme_filter$UserName),FUN=sum)
rownames(cazyme_sum) <- cazyme_sum$Group.1;cazyme_sum <- cazyme_sum[-1]

cazyme_sum_for_fil <- cazyme_sum
cazyme_sum_for_fil[cazyme_sum_for_fil>0] <- 1
cazyme_sum <- cazyme_sum[colSums(cazyme_sum_for_fil)>34*0.25]

cazyme_to_keep <- colnames(cazyme_sum)

cazyme_co_var <- aggregate(t(cazyme[cazyme_to_keep,]),by=list(cazyme_filter$UserName),FUN=function(x){sd(x)/(mean(x))})
rownames(cazyme_co_var) <- cazyme_co_var$Group.1;cazyme_co_var <- cazyme_co_var[-1]

cazyme_co_var[is.na(cazyme_co_var)] <- 0
cazyme_co_var <- cazyme_co_var[,order(colMeans(cazyme_co_var),decreasing = T)]

hist(colMeans(cazyme_co_var))
cazyme_co_var <- cazyme_co_var[,1:100]

cazyme_to_keep <- colnames(cazyme_co_var)

save(cazyme_to_keep,file='./data/cazymes_to_keep_0801.RData')



cazyme <- cazyme[cazyme_to_keep,]*1.0
cazyme_trans <- t(cazyme)
cazyme_trans <- cazyme_trans[rowSums(cazyme_trans)!=0,]
cazyme_trans <- as.data.frame(sweep(cazyme_trans,1,rowSums(cazyme_trans),'/'))


cazyme_filter <- rownames_to_column(cazyme_trans,var = 'X.SampleID')
cazyme_filter <- merge(cazyme_filter,map[c('X.SampleID','UserName','StudyDayNo')],by = 'X.SampleID')


myimpR = impRZilr(cazyme_trans, maxit = 3, method = "lm", dl = rep(0.0000001,ncol(cazyme_trans)), verbose = T) # pick a good detection limit...
cazyme_trans_imp <- myimpR$x
cazyme_trans_clr <- cenLR(cazyme_trans_imp)$x

cazyme_mean <- aggregate(cazyme_trans_clr,by=list(cazyme_filter$UserName),FUN=mean)
cazyme_mean <- aggregate(cazyme_trans,by=list(cazyme_filter$UserName),FUN=mean)
rownames(cazyme_mean) <- cazyme_mean$Group.1;cazyme_mean <- cazyme_mean[-1]
