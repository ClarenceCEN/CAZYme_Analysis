require(robCompositions)
require(tibble)
require(dplyr)
require(vegan)
require(ggplot2)
require(ggsci)
require(mice)

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

#clr of cazyme data
cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
#cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

cazyme_filter <- as.data.frame(t(cazyme))
cazyme_filter <- rownames_to_column(cazyme_filter,var = 'X.SampleID')
cazyme_filter <- merge(cazyme_filter,map[c('X.SampleID','UserName','StudyDayNo')],by = 'X.SampleID')

cazyme_filter <- cazyme_filter[,-grep('Other',colnames(cazyme_filter))]
colnames(cazyme_filter) <- gsub(".*;L4_",'',colnames(cazyme_filter))
#cazyme_filter <- cazyme_filter[,colnames(cazyme_filter)!='X.SampleID']

cazyme_list <- list()
for(i in unique(cazyme_filter$UserName)){
  print(i)
  temp_cazyme <- cazyme_filter[cazyme_filter$UserName==i,!colnames(cazyme_filter)%in%c('X.SampleID','UserName','StudyDayNo')]
  #temp_cazyme <- sweep(temp_cazyme,1,rowSums(temp_cazyme),'/')
  temp_cazyme <- sweep(temp_cazyme,1,rowSums(temp_cazyme),'/')
  temp <- as.data.frame(apply(temp_cazyme, 2, function(x){
    iqr <- IQR(x)
    q_25 <-quantile(x)['25%']
    q_75 <- quantile(x)['75%']
    x[x==0] <- NA
    x[x>=q_75+1.5*iqr||x<=q_25-1.5*iqr] <- NA
    x
  }))
  
  temp2 <- sapply(temp,function(x){
    sum(is.na(x))
  })
  
  temp <- temp[,temp2<nrow(temp_cazyme)*0.25]
  temp$StudyDayNo <- cazyme_filter[cazyme_filter$UserName==i,'StudyDayNo']
  
  imputed_Data <- mice(temp, m=3, maxit = 30, method = 'pmm', seed = 500)
  temp_imp <- complete(imputed_Data,2)
  
  temp3 <- sapply(temp_imp,function(x){
    sum(is.na(x))
  })
  
  temp_imp <- temp_imp[,temp3==0]
  
  temp_imp$X.SampleID <- cazyme_filter[cazyme_filter$UserName==i,'X.SampleID']
  
  cazyme_list[[i]] <- temp_imp
}

save(cazyme_list,file = './data/cazyme_list.RData')

load("./data/cazyme_list.RData")
cazyme_list_clr <- lapply(cazyme_list, function(x){
  b <- x[,!colnames(x)%in%c("StudyDayNo",'X.SampleID')]
  a <- as.data.frame(cenLR(b)$x)
  a$X.SampleID <- x$X.SampleID
  a
})
save(cazyme_list_clr,file = './data/cazy_list_clr.RData')


example.data <- data.frame(y=cazyme_list$MCTs37[['PL8']],x=1:nrow(cazyme_list$MCTs37))
example.data$color1 <- ifelse(example.data$y==0,'zero','non-zero')

#example.data$y <- example.data$y/sum(example.data$y)

g <- ggplot(example.data,aes(x=x,y=y)) + geom_point(aes(color=color1),size=5) + theme_classic() +
  scale_colour_npg() + ylab("MCTs10_CBM20") +xlab("Days")
g

quantile(example.data$y)
IQR(example.data$y) * 1.5
boxplot(example.data$y)

example.data[example.data==0] <- NA

#impute missing value with time
imputed_Data <- mice(example.data[1:2], m=5, maxit = 50, method = 'pmm', seed = 500)
summary(imputed_Data)

example.data_imp <- complete(imputed_Data,2)
example.data_imp$color1 <- example.data$color1

g <- ggplot(example.data_imp,aes(x=x,y=y)) + geom_point(aes(color=color1),size=5) + theme_classic() +
  scale_colour_npg() + ylab("MCTs10_CBM20") +xlab("Days")
g

#impute missing value without time
imputed_Data <- mice(example.data[1], m=5, maxit = 50, method = 'pmm', seed = 500)
summary(imputed_Data)

example.data_imp <- complete(imputed_Data,2)
example.data_imp$color1 <- example.data$color1

g <- ggplot(example.data_imp,aes(x=x,y=y)) + geom_point(aes(color=color1),size=5) + theme_classic() +
  scale_colour_npg() + ylab("MCTs10_CBM20") +xlab("Days")
g


temp <- as.data.frame(apply(temp_cazyme, 2, function(x){
  print(x)
  iqr <- IQR(x)
  q_25 <-quantile(x)['25%']
  q_75 <- quantile(x)['75%']
  x[x==0] <- NA
  x[x>=q_75+1.5*iqr||x<=q_25-1.5*iqr] <- NA
  x
}))

temp <- temp[,sapply(temp, function)]

temp2 <- sapply(temp,function(x){
  sum(is.na(x))
})
