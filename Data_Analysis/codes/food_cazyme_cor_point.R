require(robCompositions)
require(tibble)
require(psych)
require(reshape2)
require(dplyr)

# load colors
source(file = "G:/Dan_Lab/dietstudy_analyses-master/lib/colors/UserNameColors.R")
UserNameColors['MCTs05'] <- '#99ff00'

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

food <- read.table('./data/masterdecaydiet.txt',header = T)
food$taxonomy <- as.character(food$taxonomy)

split <- strsplit(food$taxonomy,";") 

foodStrings <- sapply(split,function(x) paste(x[1:2],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",foodStrings,perl=T) # clean tips
food_L2 <- rowsum(food[-ncol(food)],foodStrings) 
rownames(food_L2) <- gsub(".*;L2_",'',rownames(food_L2))
#rownames(food_L2) = sapply(strsplit(food$taxonomy,";"),function(x) paste(x[1:2],collapse=";"));

load("./data/cazymes_to_keep.RData")
cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
#cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
cazyme <- cazyme[cazymes_to_keep,]

cazyme <- cazyme[-grep('Other',rownames(cazyme)),]
rownames(cazyme) <- gsub(".*;L4_",'',rownames(cazyme))

map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# identify soylent samples
soylent <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
soylent <- as.character(soylent$X.SampleID)

cazyme <- cazyme[,!colnames(cazyme)%in%soylent]
food_L2 <- food_L2[,!colnames(food_L2)%in%soylent]

samples <- intersect(colnames(cazyme),colnames(food_L2))
cazyme <- cazyme[samples]
food_L2 <- food_L2[samples]

cazyme_c <- as.data.frame(t(cazyme))
cazyme_c <- rownames_to_column(cazyme_c,var = 'X.SampleID')
cazyme_c <- merge(cazyme_c,map[c('X.SampleID','UserName')],by='X.SampleID')

food_c <- as.data.frame(t(food_L2))
food_c <- rownames_to_column(food_c,var='X.SampleID')
food_c <- merge(food_c,map[c('X.SampleID','UserName')],by='X.SampleID')


load('./data/cazyme_food_cor.RData')
sigs <- lapply(cazyme_food_list, function(x) subset(x, fdr_p <= 0.2))
allsigs <- do.call("rbind", sigs)

allsigs$cazy_cat <- ifelse(grepl('AA',allsigs$cazyme),'AA',
                           ifelse(grepl('CBM',allsigs$cazyme),'CBM',
                                  ifelse(grepl('GT',allsigs$cazyme),'GT',
                                         ifelse(grepl('GH',allsigs$cazyme),'GH',
                                                ifelse(grepl('PL',allsigs$cazyme),'PL','CE')))))

allsigs$bin <- ifelse(allsigs$coef < 0, "Negative (+/-)", 
                      ifelse(allsigs$coef > 0, "Positive (+/+ or -/-)", "NA"))
allsigs$bin <- factor(allsigs$bin, levels = c("Positive (+/+ or -/-)", "Negative (+/-)"))

allsigs$pairs <- paste0(allsigs$food,'_',allsigs$cazyme)
cazyme_set <- as.character(unique(allsigs$pairs))


select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))>1){
    print(i)
    select_cazymes <- c(select_cazymes,i)
  }
}

select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))==1 && length(unique(temp_data$id))>1){
    print(i)
    select_cazymes <- c(select_cazymes,i)
  }
}

selected_data <- allsigs[allsigs$pairs%in%select_cazymes,]

selected_food_df <- list()
for(i in unique(selected_data$pairs)){
  print(i)

  selected_food_df[[i]] <- subset(selected_data,selected_data$pairs==i)
}

food_c$UserName <- as.character(food_c$UserName)
cazyme_c$UserName <- as.character(cazyme_c$UserName)

selected_food <- list()
selected_food <- lapply(selected_food_df,function(x){
    food_c[food_c$UserName%in%as.character(x$id),colnames(food_c)%in%c("UserName",as.character(x$food))] %>%
    melt(id.vars = 'UserName',variable.name = 'food',value.name = 'weight')
})

selected_cazyme <- list()
selected_cazyme <- lapply(selected_food_df,function(x){
  cazyme_c[cazyme_c$UserName%in%as.character(x$id),colnames(cazyme_c)%in%c("UserName",as.character(x$cazyme))] %>%
    melt(id.vars = 'UserName',variable.name = 'cazyme',value.name = 'count')
})


selected_cazyme_plot <- do.call("rbind", selected_cazyme)
selected_food_plot <- do.call("rbind", selected_food)

#selected_plot <- merge(selected_food,selected_cazyme,by='UserName',all=F)
selected_plot <- bind_cols(selected_food_plot,selected_cazyme_plot)
selected_plot <- selected_plot[,names(selected_plot)!='UserName1']

selected_plot <- mutate_if(selected_plot,is.factor,as.character)
unique(selected_plot$food)

#Organ_meats_sausages_and_lunchmeats_CBM4
plot1 <- subset(selected_plot,selected_plot$food=='Organ_meats_sausages_and_lunchmeats')
plot1 <- filter(plot1,plot1$cazyme=='CBM4')

g <- ggplot(plot1,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(~UserName,scales = 'free') +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  guides(size=NULL) +
  ylab('CBM4') + xlab("Organ_meats_sausages_and_lunchmeats")
g
ggsave("./result/linear_plots/L2_food/Organ_meats_sausages_and_lunchmeats_CBM4.pdf",
       width = 5,height = 5)

#Organ_meats_sausages_and_lunchmeats_GH26
plot2 <- subset(selected_plot,selected_plot$food=='Organ_meats_sausages_and_lunchmeats')
plot2 <- filter(plot2,plot2$cazyme=='GH26')

g <- ggplot(plot2,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(~UserName,scales = 'free') +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  guides(size=NULL) +
  ylab('GH26') + xlab("Organ_meats_sausages_and_lunchmeats")
g
ggsave("./result/linear_plots/L2_food/Organ_meats_sausages_and_lunchmeats_GH26.pdf",
       width = 5,height = 5)


#Deepyellow_vegetables_GH137
plot3 <- subset(selected_plot,selected_plot$food=='Deepyellow_vegetables')
plot3 <- filter(plot3,plot3$cazyme=='GH137')

g <- ggplot(plot3,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(.~UserName,scales = "free") +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GH137') + xlab("Deepyellow_vegetables_GH137")
g
ggsave("./result/linear_plots/L2_food/Deepyellow_vegetables_GH137.pdf",
       width = 5,height = 5)

#Citrus_fruits_juices_CBM48
plot4 <- subset(selected_plot,selected_plot$food=='Citrus_fruits_juices')
plot4 <- filter(plot4,plot4$cazyme=='CBM48')

g <- ggplot(plot4,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(.~UserName,scales = "free") +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('CBM48') + xlab("Citrus_fruits_juices")
g
ggsave("./result/linear_plots/L2_food/Citrus_fruits_juices_CBM48.pdf",
       width = 5,height = 5)

#Citrus_fruits_juices_GT5
plot5 <- subset(selected_plot,selected_plot$food=='Citrus_fruits_juices')
plot5 <- filter(plot5,plot5$cazyme=='GT5')

g <- ggplot(plot5,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(.~UserName,scales = "free") +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GT5') + xlab("Citrus_fruits_juices")
g
ggsave("./result/linear_plots/L2_food/Citrus_fruits_juices_GT5.pdf",
       width = 5,height = 5)

#Cereals_not_cooked_or_NS_as_to_cooked_GH51
plot6 <- subset(selected_plot,selected_plot$food=='Cereals_not_cooked_or_NS_as_to_cooked')
plot6 <- filter(plot6,plot6$cazyme=='GH51')

g <- ggplot(plot6,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(.~UserName,scales = "free") +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GH51') + xlab("Cereals_not_cooked_or_NS_as_to_cooked")
g
ggsave("./result/linear_plots/L2_food/Cereals_not_cooked_or_NS_as_to_cooked_GH51.pdf",
       width = 5,height = 5)


#Quick_breads_GH137
plot7 <- subset(selected_plot,selected_plot$food=='Quick_breads')
plot7 <- filter(plot7,plot7$cazyme=='GH137')

g <- ggplot(plot7,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(.~UserName,scales = "free") +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GH137') + xlab("Quick_breads")
g
ggsave("./result/linear_plots/L2_food/Quick_breads_GH137.pdf",
       width = 5,height = 5)

#Nonalcoholic_beverages_GT51
plot8 <- subset(selected_plot,selected_plot$food=='Nonalcoholic_beverages')
plot8 <- filter(plot8,plot8$cazyme=='GT51')

g <- ggplot(plot8,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(.~UserName,scales = "free") +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GT51') + xlab("Nonalcoholic_beverages")
g
ggsave("./result/linear_plots/L2_food/Nonalcoholic_beverages_GT51.pdf",
       width = 5,height = 5)

#L3
setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

food <- read.table('./data/masterdecaydiet.txt',header = T)
food$taxonomy <- as.character(food$taxonomy)

split <- strsplit(food$taxonomy,";") 

foodStrings <- sapply(split,function(x) paste(x[1:3],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",foodStrings,perl=T) # clean tips
food_L3 <- rowsum(food[-ncol(food)],foodStrings) 
rownames(food_L3) <- gsub(".*;L3_",'',rownames(food_L3))
#rownames(food_L2) = sapply(strsplit(food$taxonomy,";"),function(x) paste(x[1:2],collapse=";"));

load("./data/cazymes_to_keep.RData")
cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
#cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
cazyme <- cazyme[cazymes_to_keep,]

cazyme <- cazyme[-grep('Other',rownames(cazyme)),]
rownames(cazyme) <- gsub(".*;L4_",'',rownames(cazyme))

map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# identify soylent samples
soylent <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
soylent <- as.character(soylent$X.SampleID)

cazyme <- cazyme[,!colnames(cazyme)%in%soylent]
food_L3 <- food_L3[,!colnames(food_L3)%in%soylent]

samples <- intersect(colnames(cazyme),colnames(food_L3))
cazyme <- cazyme[samples]
food_L3 <- food_L3[samples]

cazyme_c <- as.data.frame(t(cazyme))
cazyme_c <- rownames_to_column(cazyme_c,var = 'X.SampleID')
cazyme_c <- merge(cazyme_c,map[c('X.SampleID','UserName')],by='X.SampleID')

food_c <- as.data.frame(t(food_L3))
food_c <- rownames_to_column(food_c,var='X.SampleID')
food_c <- merge(food_c,map[c('X.SampleID','UserName')],by='X.SampleID')

cazyme_food_list <- list()

for(i in unique(food_c$UserName)){
  print(i)
  coef_v <- c()
  fdr_v <- c()
  food_v <- c()
  cazyme_v <- c()
  p_v <- c()
  id <- c()
  temp_food <- food_c[food_c$UserName==i,!colnames(food_c)%in%c('X.SampleID','UserName')]
  temp_cazyme <- cazyme_c[cazyme_c$UserName==i,!colnames(cazyme_c)%in%c('X.SampleID','UserName')]
  temp_cor_r <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$r
  temp_cor_p <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$p
  for(a in 1:nrow(temp_cor_p)){
    #print(a)
    for (b in 1:ncol(temp_cor_p)){
      if(!is.na(temp_cor_p[a,b])){
        coef_v <- c(coef_v,temp_cor_r[a,b])
        fdr_v <- c(fdr_v,temp_cor_p[a,b])
        food_v <- c(food_v,colnames(temp_food)[a])
        cazyme_v <- c(cazyme_v,colnames(temp_cazyme)[b])
        id <- c(id,i)
      }
    }
  }
  #fdr_v <- p.adjust(p_v,method = 'fdr')
  temp_cor_df <- data.frame(food=food_v,cazyme=cazyme_v,coef=coef_v,id=id,fdr_p=fdr_v)
  cazyme_food_list[[i]] <- temp_cor_df
}
save(cazyme_food_list,file = './data/cazyme_food_cor_L3.RData')
#save(cazyme_food_list,file = './data/cazyme_food_cor_try_L3.RData')

load('./data/cazyme_food_cor_L3.RData')
sigs <- lapply(cazyme_food_list, function(x) subset(x, fdr_p <= 0.2))
allsigs <- do.call("rbind", sigs)

allsigs$cazy_cat <- ifelse(grepl('AA',allsigs$cazyme),'AA',
                           ifelse(grepl('CBM',allsigs$cazyme),'CBM',
                                  ifelse(grepl('GT',allsigs$cazyme),'GT',
                                         ifelse(grepl('GH',allsigs$cazyme),'GH',
                                                ifelse(grepl('PL',allsigs$cazyme),'PL','CE')))))

allsigs$bin <- ifelse(allsigs$coef < 0, "Negative (+/-)", 
                      ifelse(allsigs$coef > 0, "Positive (+/+ or -/-)", "NA"))
allsigs$bin <- factor(allsigs$bin, levels = c("Positive (+/+ or -/-)", "Negative (+/-)"))

allsigs$pairs <- paste0(allsigs$food,'_',allsigs$cazyme)
cazyme_set <- as.character(unique(allsigs$pairs))


select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))>1){
    print(i)
    select_cazymes <- c(select_cazymes,i)
  }
}

selected_data <- allsigs[allsigs$pairs%in%select_cazymes,]

selected_food_df <- list()
for(i in unique(selected_data$pairs)){
  print(i)
  
  selected_food_df[[i]] <- subset(selected_data,selected_data$pairs==i)
}

food_c$UserName <- as.character(food_c$UserName)
cazyme_c$UserName <- as.character(cazyme_c$UserName)

selected_food <- list()
selected_food <- lapply(selected_food_df,function(x){
  food_c[food_c$UserName%in%as.character(x$id),colnames(food_c)%in%c("UserName",as.character(x$food))] %>%
    melt(id.vars = 'UserName',variable.name = 'food',value.name = 'weight')
})

selected_cazyme <- list()
selected_cazyme <- lapply(selected_food_df,function(x){
  cazyme_c[cazyme_c$UserName%in%as.character(x$id),colnames(cazyme_c)%in%c("UserName",as.character(x$cazyme))] %>%
    melt(id.vars = 'UserName',variable.name = 'cazyme',value.name = 'count')
})


selected_cazyme_plot <- do.call("rbind", selected_cazyme)
selected_food_plot <- do.call("rbind", selected_food)

#selected_plot <- merge(selected_food,selected_cazyme,by='UserName',all=F)
selected_plot <- bind_cols(selected_food_plot,selected_cazyme_plot)
selected_plot <- selected_plot[,names(selected_plot)!='UserName1']

selected_plot <- mutate_if(selected_plot,is.factor,as.character)
unique(selected_plot$food)

#Carrots_GH137
plot1 <- subset(selected_plot,selected_plot$food=='Carrots')
plot1 <- filter(plot1,plot1$cazyme=='GH137')

g <- ggplot(plot1,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(~UserName,scales = 'free') +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GH137') + xlab("Carrots")
g
ggsave("./result/linear_plots/L3_food/Carrots_GH137.pdf",
       width = 5,height = 5)

#Citrus_fruits_CBM48
plot2 <- subset(selected_plot,selected_plot$food=='Citrus_fruits')
plot2 <- filter(plot2,plot2$cazyme=='CBM48')

g <- ggplot(plot2,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(~UserName,scales = 'free') +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('CBM48') + xlab("Citrus_fruits")
g
ggsave("./result/linear_plots/L3_food/Citrus_fruits_CBM48.pdf",
       width = 5,height = 5)

#Citrus_fruits_GT5
plot3 <- subset(selected_plot,selected_plot$food=='Citrus_fruits')
plot3 <- filter(plot3,plot3$cazyme=='GT5')

g <- ggplot(plot3,aes(x=weight,y=count))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(~UserName,scales = 'free') +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GT5') + xlab("Citrus_fruits")
g
ggsave("./result/linear_plots/L3_food/Citrus_fruits_GT5.pdf",
       width = 5,height = 5)


#L1
setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

food <- read.table('./data/masterdecaydiet.txt',header = T)
food$taxonomy <- as.character(food$taxonomy)

split <- strsplit(food$taxonomy,";") 

foodStrings <- sapply(split,function(x) paste(x[1:1],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",foodStrings,perl=T) # clean tips
food_L1 <- rowsum(food[-ncol(food)],foodStrings) 
rownames(food_L1) <- gsub(".*;L1_",'',rownames(food_L1))
#rownames(food_L2) = sapply(strsplit(food$taxonomy,";"),function(x) paste(x[1:2],collapse=";"));

load("./data/cazymes_to_keep.RData")
cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
#cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
cazyme <- cazyme[cazymes_to_keep,]

cazyme <- cazyme[-grep('Other',rownames(cazyme)),]
rownames(cazyme) <- gsub(".*;L4_",'',rownames(cazyme))

map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# identify soylent samples
soylent <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
soylent <- as.character(soylent$X.SampleID)

cazyme <- cazyme[,!colnames(cazyme)%in%soylent]
food_L1 <- food_L1[,!colnames(food_L1)%in%soylent]

samples <- intersect(colnames(cazyme),colnames(food_L1))
cazyme <- cazyme[samples]
food_L1 <- food_L1[samples]

cazyme_c <- as.data.frame(t(cazyme))
cazyme_c <- rownames_to_column(cazyme_c,var = 'X.SampleID')
cazyme_c <- merge(cazyme_c,map[c('X.SampleID','UserName')],by='X.SampleID')

food_c <- as.data.frame(t(food_L1))
food_c <- rownames_to_column(food_c,var='X.SampleID')
food_c <- merge(food_c,map[c('X.SampleID','UserName')],by='X.SampleID')

cazyme_food_list <- list()

for(i in unique(food_c$UserName)){
  print(i)
  coef_v <- c()
  fdr_v <- c()
  food_v <- c()
  cazyme_v <- c()
  p_v <- c()
  id <- c()
  temp_food <- food_c[food_c$UserName==i,!colnames(food_c)%in%c('X.SampleID','UserName')]
  temp_cazyme <- cazyme_c[cazyme_c$UserName==i,!colnames(cazyme_c)%in%c('X.SampleID','UserName')]
  temp_cor_r <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$r
  temp_cor_p <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$p
  for(a in 1:nrow(temp_cor_p)){
    #print(a)
    for (b in 1:ncol(temp_cor_p)){
      if(!is.na(temp_cor_p[a,b])){
        coef_v <- c(coef_v,temp_cor_r[a,b])
        fdr_v <- c(fdr_v,temp_cor_p[a,b])
        food_v <- c(food_v,colnames(temp_food)[a])
        cazyme_v <- c(cazyme_v,colnames(temp_cazyme)[b])
        id <- c(id,i)
      }
    }
  }
  #fdr_v <- p.adjust(p_v,method = 'fdr')
  temp_cor_df <- data.frame(food=food_v,cazyme=cazyme_v,coef=coef_v,id=id,fdr_p=fdr_v)
  cazyme_food_list[[i]] <- temp_cor_df
}
save(cazyme_food_list,file = './data/cazyme_food_cor_L1.RData')
#save(cazyme_food_list,file = './data/cazyme_food_cor_try_L3.RData')

load('./data/cazyme_food_cor_L1.RData')
sigs <- lapply(cazyme_food_list, function(x) subset(x, fdr_p <= 0.2))
allsigs <- do.call("rbind", sigs)

allsigs$cazy_cat <- ifelse(grepl('AA',allsigs$cazyme),'AA',
                           ifelse(grepl('CBM',allsigs$cazyme),'CBM',
                                  ifelse(grepl('GT',allsigs$cazyme),'GT',
                                         ifelse(grepl('GH',allsigs$cazyme),'GH',
                                                ifelse(grepl('PL',allsigs$cazyme),'PL','CE')))))

allsigs$bin <- ifelse(allsigs$coef < 0, "Negative (+/-)", 
                      ifelse(allsigs$coef > 0, "Positive (+/+ or -/-)", "NA"))
allsigs$bin <- factor(allsigs$bin, levels = c("Positive (+/+ or -/-)", "Negative (+/-)"))

allsigs$pairs <- paste0(allsigs$food,'_',allsigs$cazyme)
cazyme_set <- as.character(unique(allsigs$pairs))


select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))>1){
    print(i)
    select_cazymes <- c(select_cazymes,i)
  }
}

