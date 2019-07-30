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
cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
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


load('./data/cazyme_food_cor_try.RData')
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

#Citrus_fruits_juices
plot1 <- subset(selected_plot,selected_plot$food=='Citrus_fruits_juices')
plot1 <- filter(plot1,plot1$UserName!='MCTs05')

g <- ggplot(plot1,aes(x=weight,y=count))+geom_point(aes(color=UserName),alpha=0.8)+
  facet_grid(~UserName,scales = 'free_x')
g

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
#cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
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
save(cazyme_food_list,file = './data/cazyme_food_cor_try_L3.RData')

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


selected_food <- food_c[food_c$UserName%in%as.character(selected_data$id),colnames(food_c)%in%c("UserName",as.character(selected_data$food))]
selected_cazyme <- cazyme_c[cazyme_c$UserName%in%as.character(selected_data$id),colnames(cazyme_c)%in%c("UserName",as.character(selected_data$cazyme))]

selected_food <- melt(selected_food,id.vars = 'UserName',variable.name = 'food',value.name = 'weight')
selected_cazyme <- melt(selected_cazyme,id.vars = 'UserName',variable.name = 'cazyme',value.name = 'count')

#selected_plot <- merge(selected_food,selected_cazyme,by='UserName',all=F)
selected_plot <- bind_cols(selected_food,selected_cazyme)
selected_plot <- selected_plot[,names(selected_plot)!='UserName1']

g <- ggplot(selected_plot,aes(x=weight,y=count))+geom_point(aes(color=UserName))
g


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

selected_food <- food_c[food_c$UserName%in%as.character(selected_data$id),colnames(food_c)%in%c("UserName",as.character(selected_data$food))]
selected_cazyme <- cazyme_c[cazyme_c$UserName%in%as.character(selected_data$id),colnames(cazyme_c)%in%c("UserName",as.character(selected_data$cazyme))]

selected_food$UserName <- as.character(selected_food$UserName)
selected_cazyme$UserName <- as.character(selected_cazyme$UserName)

selected_plot <- merge(selected_food,selected_cazyme,by='UserName',all=F)
selected_plot <- bind_cols(selected_food,selected_cazyme)
selected_plot <- selected_plot[-4]

g <- ggplot(selected_plot,aes(x=Citrus_fruits_juices,y=CBM48))+geom_point(aes(color=UserName))
g
