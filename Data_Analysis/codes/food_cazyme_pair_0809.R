require(robCompositions)
require(tibble)
require(psych)
require(reshape2)
require(dplyr)

load("./data/cazyme_food_cor_ind_0808.RData")
load("./data/cazyme_food_cor_ind_0815.RData")
load('./data/cazy_list_clr_2.RData')

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

map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# identify soylent samples
soylent <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
soylent <- as.character(soylent$X.SampleID)

food_L2 <- food_L2[,!colnames(food_L2)%in%soylent]


food_daily <- as.data.frame(t(food_L2))
food_daily <- rownames_to_column(food_daily,var='X.SampleID')
food_daily <- merge(food_daily,map[c('X.SampleID','UserName')],by='X.SampleID')

food_daily_list <- list()
for(i in unique(food_daily$UserName)){
  food_daily_list[[i]] <- select(food_daily[food_daily$UserName==i,],-UserName)
}

food_daily_list <- lapply(food_daily_list, function(x){
  temp <- select(x,-X.SampleID)
  temp <- temp[,colSums(temp)!=0]
  temp$X.SampleID <- x$X.SampleID
  temp
})

save(food_daily_list,file = './data/food_daily_L2.RData')

food_long <- lapply(food_daily_list, function(x){
  #print(x)
  x %>% melt(id.vars='X.SampleID',variable.name = 'food',value.name = 'weight') -> x
  x
})

food_long <- do.call('rbind',food_long) %>% mutate_if(is.factor,as.character)

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

#different
select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))>1&& length(unique(temp_data$id))>1){
    print(i)
    print(temp_data$id)
    select_cazymes <- c(select_cazymes,i)
  }
}

#same
select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))==1 && length(unique(temp_data$id))>1){
    print(i)
    print(temp_data$id)
    select_cazymes <- c(select_cazymes,i)
  }
}

selected_data <- allsigs[allsigs$pairs%in%select_cazymes,]

selected_food_df <- list()
for(i in unique(selected_data$pairs)){
  print(i)

  selected_food_df[[i]] <- subset(selected_data,selected_data$pairs==i)
}

#food_daily$UserName <- as.character(food_daily$UserName)
#food_daily <- melt(food_daily,id.vars = c('X.SampleID','UserName'),variable.name = 'food',
#                   value.name = 'weight')

cazyme_long <- lapply(cazyme_list_clr, function(x){
  #print(x)
  x %>% melt(id.vars='X.SampleID',variable.name = 'cazyme',value.name = 'relative_abundance_clr') -> x
  x
})

cazyme_long <- do.call('rbind',cazyme_long) %>% mutate_if(is.factor,as.character)

food_cazyme_long <- merge(food_long,cazyme_long,by='X.SampleID')
food_cazyme_long <- merge(food_cazyme_long,map[,c('X.SampleID','UserName')],by='X.SampleID')


# selected_food <- list()
# selected_food <- lapply(selected_food_df,function(x){
#     food_daily[food_daily$UserName%in%as.character(x$id),colnames(food_daily)%in%c("UserName",as.character(x$food))] %>%
#     melt(id.vars = 'UserName',variable.name = 'food',value.name = 'weight')
# })
# 
# 
# selected_cazyme_plot <- do.call("rbind", selected_cazyme)
# selected_food_plot <- do.call("rbind", selected_food)
# 
# #selected_plot <- merge(selected_food,selected_cazyme,by='UserName',all=F)
# selected_plot <- bind_cols(selected_food_plot,selected_cazyme_plot)
# selected_plot <- selected_plot[,names(selected_plot)!='UserName1']
# 
# selected_plot <- mutate_if(selected_plot,is.factor,as.character)
# unique(selected_plot$food)

#diff
#Deepyellow_vegetables_GH92
plot1 <- subset(food_cazyme_long,food_cazyme_long$food=='Deepyellow_vegetables')
plot1 <- filter(plot1,plot1$cazyme=='GH92',plot1$UserName%in%c("MCTs07",'MCTs16'))

g <- ggplot(plot1,aes(x=weight,y=relative_abundance_clr))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(~UserName,scales = 'free') +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  guides(size=NULL) +
  ylab('GH92 (clr-adjusted relative abundance)') + xlab("Deepyellow_vegetables")
g
ggsave("./result/linear_plots/L2_ind/Deepyellow_vegetables_GH92.pdf",
       width = 5,height = 5)

#Cakes_cookies_pies_pastries_bars_GH106
plot2 <- subset(food_cazyme_long,food_cazyme_long$food=='Cakes_cookies_pies_pastries_bars')
plot2 <- filter(plot2,plot2$cazyme=='GH106',plot2$UserName%in%c("MCTs32",'MCTs33'))

g <- ggplot(plot2,aes(x=weight,y=relative_abundance_clr))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(~UserName,scales = 'free') +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  guides(size=NULL) +
  ylab('GH106 (clr-adjusted relative abundance)') + xlab("Cakes_cookies_pies_pastries_bars")
g
ggsave("./result/linear_plots/L2_ind/Cakes_cookies_pies_pastries_bars_GH106.pdf",
       width = 5,height = 5)


#same
#White_potatoes_and_Puerto_Rican_starchy_vegetables_GH43
plot4 <- subset(food_cazyme_long,food_cazyme_long$food=='White_potatoes_and_Puerto_Rican_starchy_vegetables')
plot4 <- filter(plot4,plot4$cazyme=='GH43',plot4$UserName%in%c("MCTs18",'MCTs32'))

g <- ggplot(plot4,aes(x=weight,y=relative_abundance_clr))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(.~UserName,scales = "free") +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GH43 (clr-adjusted relative abundance)') + xlab("White_potatoes_and_Puerto_Rican_starchy_vegetables")
g
ggsave("./result/linear_plots/L2_ind/White_potatoes_and_Puerto_Rican_starchy_vegetables_GH43.pdf",
       width = 5,height = 5)

#Organ_meats_sausages_and_lunchmeats_GH28
plot5 <- subset(food_cazyme_long,food_cazyme_long$food=='Organ_meats_sausages_and_lunchmeats')
plot5 <- filter(plot5,plot5$cazyme=='GH28',plot5$UserName%in%c("MCTs18",'MCTs32'))

g <- ggplot(plot5,aes(x=weight,y=relative_abundance_clr))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(.~UserName,scales = "free") +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GH28 (clr-adjusted relative abundance)') + xlab("Organ_meats_sausages_and_lunchmeats")
g
ggsave("./result/linear_plots/L2_ind/Organ_meats_sausages_and_lunchmeats_GH28.pdf",
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


map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# identify soylent samples
soylent <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
soylent <- as.character(soylent$X.SampleID)

food_L3 <- food_L3[,!colnames(food_L3)%in%soylent]


food_daily <- as.data.frame(t(food_L3))
food_daily <- rownames_to_column(food_daily,var='X.SampleID')
food_daily <- merge(food_daily,map[c('X.SampleID','UserName')],by='X.SampleID')

food_daily_list <- list()
for(i in unique(food_daily$UserName)){
  food_daily_list[[i]] <- select(food_daily[food_daily$UserName==i,],-UserName)
}

food_daily_list <- lapply(food_daily_list, function(x){
  temp <- select(x,-X.SampleID)
  temp <- temp[,colSums(temp)!=0]
  temp$X.SampleID <- x$X.SampleID
  temp
})

food_long <- lapply(food_daily_list, function(x){
  #print(x)
  x %>% melt(id.vars='X.SampleID',variable.name = 'food',value.name = 'weight') -> x
  x
})

food_long <- do.call('rbind',food_long) %>% mutate_if(is.factor,as.character)

cazyme_food_list <- list()


for(i in names(food_daily_list)){
  print(i)
  coef_v <- c()
  fdr_v <- c()
  food_v <- c()
  cazyme_v <- c()
  p_v <- c()
  id <- c()
  temp_food <- food_daily_list[[i]] #%>% select(-X.SampleID)
  temp_cazyme <- cazyme_list_clr[[i]] #%>% select(-X.SampleID)
  temp_food <- temp_food[temp_food$X.SampleID%in%intersect(temp_food$X.SampleID,temp_cazyme$X.SampleID),] 
  temp_cazyme <- temp_cazyme[temp_cazyme$X.SampleID%in%intersect(temp_food$X.SampleID,temp_cazyme$X.SampleID),]%>% select(-X.SampleID)
  temp_food %>%select(-X.SampleID) ->temp_food
  temp_cor_r <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$r
  temp_cor_p <- corr.test(temp_food,temp_cazyme,adjust = 'none',method = 'spearman')$p
  temp_cor_fdr <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$p
  for(a in 1:nrow(temp_cor_p)){
    #print(a)
    for (b in 1:ncol(temp_cor_p)){
      if(!is.na(temp_cor_p[a,b])){
        coef_v <- c(coef_v,temp_cor_r[a,b])
        fdr_v <- c(fdr_v,temp_cor_fdr[a,b])
        food_v <- c(food_v,colnames(temp_food)[a])
        cazyme_v <- c(cazyme_v,colnames(temp_cazyme)[b])
        p_v <- c(p_v,temp_cor_p[a,b])
        id <- c(id,i)
      }
    }
  }
  temp_cor_df <- data.frame(food=food_v,cazyme=cazyme_v,coef=coef_v,id=id,fdr_p=fdr_v,p=p_v)
  cazyme_food_list[[i]] <- temp_cor_df
}
save(cazyme_food_list,file = './data/cazyme_food_cor_ind_L3.RData')
#save(cazyme_food_list,file = './data/cazyme_food_cor_try_L3.RData')

load('./data/cazyme_food_cor_indL3.RData')
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


#different
select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))>1){
    print(i)
    print(temp_data$id)
    select_cazymes <- c(select_cazymes,i)
  }
}

#same
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


cazyme_long <- lapply(cazyme_list_clr, function(x){
  #print(x)
  x %>% melt(id.vars='X.SampleID',variable.name = 'cazyme',value.name = 'relative_abundance_clr') -> x
  x
})

cazyme_long <- do.call('rbind',cazyme_long) %>% mutate_if(is.factor,as.character)

food_cazyme_long <- merge(food_long,cazyme_long,by='X.SampleID')
food_cazyme_long <- merge(food_cazyme_long,map[,c('X.SampleID','UserName')],by='X.SampleID')


#different
#Sour_cream_GH139
plot1 <- subset(food_cazyme_long,food_cazyme_long$food=='Sour_cream')
plot1 <- filter(plot1,plot1$cazyme=='GH139',plot1$UserName%in%c("MCTs28",'MCTs33'))

g <- ggplot(plot1,aes(x=weight,y=relative_abundance_clr))+geom_point(aes(fill=UserName),alpha=0.5,size=4,shape=21)+
  facet_grid(~UserName,scales = 'free') +theme_classic() +
  stat_smooth(method = lm,color="black",size=1,se=FALSE)+
  scale_fill_manual(values = UserNameColors) +
  ylab('GH139 (clr-adjusted relative abundance)') + xlab("Sour_cream") 
g
ggsave("./result/linear_plots/L3_ind/Sour_cream_GH139.pdf",
       width = 5,height = 5)


#L1
setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

food <- read.table('./data/masterdecaydiet.txt',header = T)
food$taxonomy <- as.character(food$taxonomy)

split <- strsplit(food$taxonomy,";") 

foodStrings <- sapply(split,function(x) paste(x[1:1],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",foodStrings,perl=T) # clean tips
food_L1 <- rowsum(food[-ncol(food)],foodStrings) 
rownames(food_L1) <- gsub("L1_",'',rownames(food_L1))
#rownames(food_L2) = sapply(strsplit(food$taxonomy,";"),function(x) paste(x[1:2],collapse=";"));


map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# identify soylent samples
soylent <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
soylent <- as.character(soylent$X.SampleID)

food_L1 <- food_L1[,!colnames(food_L1)%in%soylent]


food_daily <- as.data.frame(t(food_L1))
food_daily <- rownames_to_column(food_daily,var='X.SampleID')
food_daily <- merge(food_daily,map[c('X.SampleID','UserName')],by='X.SampleID')

food_daily_list <- list()
for(i in unique(food_daily$UserName)){
  food_daily_list[[i]] <- select(food_daily[food_daily$UserName==i,],-UserName)
}

food_daily_list <- lapply(food_daily_list, function(x){
  temp <- select(x,-X.SampleID)
  temp <- temp[,colSums(temp)!=0]
  temp$X.SampleID <- x$X.SampleID
  temp
})

food_long <- lapply(food_daily_list, function(x){
  #print(x)
  x %>% melt(id.vars='X.SampleID',variable.name = 'food',value.name = 'weight') -> x
  x
})

food_long <- do.call('rbind',food_long) %>% mutate_if(is.factor,as.character)

cazyme_food_list <- list()


for(i in names(food_daily_list)){
  print(i)
  coef_v <- c()
  fdr_v <- c()
  food_v <- c()
  cazyme_v <- c()
  p_v <- c()
  id <- c()
  temp_food <- food_daily_list[[i]] #%>% select(-X.SampleID)
  temp_cazyme <- cazyme_list_clr[[i]] #%>% select(-X.SampleID)
  temp_food <- temp_food[temp_food$X.SampleID%in%intersect(temp_food$X.SampleID,temp_cazyme$X.SampleID),] 
  temp_cazyme <- temp_cazyme[temp_cazyme$X.SampleID%in%intersect(temp_food$X.SampleID,temp_cazyme$X.SampleID),]%>% select(-X.SampleID)
  temp_food %>%select(-X.SampleID) ->temp_food
  temp_cor_r <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$r
  temp_cor_p <- corr.test(temp_food,temp_cazyme,adjust = 'none',method = 'spearman')$p
  temp_cor_fdr <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$p
  for(a in 1:nrow(temp_cor_p)){
    #print(a)
    for (b in 1:ncol(temp_cor_p)){
      if(!is.na(temp_cor_p[a,b])){
        coef_v <- c(coef_v,temp_cor_r[a,b])
        fdr_v <- c(fdr_v,temp_cor_fdr[a,b])
        food_v <- c(food_v,colnames(temp_food)[a])
        cazyme_v <- c(cazyme_v,colnames(temp_cazyme)[b])
        p_v <- c(p_v,temp_cor_p[a,b])
        id <- c(id,i)
      }
    }
  }
  temp_cor_df <- data.frame(food=food_v,cazyme=cazyme_v,coef=coef_v,id=id,fdr_p=fdr_v,p=p_v)
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

#different   #None
select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))>1){
    print(i)
    print(as.character(temp_data$id))
    select_cazymes <- c(select_cazymes,i)
  }
}

#same   #None
select_cazymes <- c()
for(i in cazyme_set){
  temp_data = allsigs[allsigs$pairs==i,]
  if(length(unique(temp_data$bin))==1 && length(unique(temp_data$id))>1){
    print(i)
    print(as.character(temp_data$id))
    select_cazymes <- c(select_cazymes,i)
  }
}


