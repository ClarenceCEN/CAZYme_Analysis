require(reshape2)
require(randomForest)
require(pheatmap)
require(dplyr)
require(ape)
require(e1071)

setwd("G:/Dan_Lab/dietstudy_analyses-master/")
load(file = "data/test_personal_diet_dat.Rdata")

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
load(file = './data/cazy_list_clr_2.RData')

food <- read.table('./data/masterdecaydiet.txt',header = T)
food$taxonomy <- as.character(food$taxonomy)

split <- strsplit(food$taxonomy,";") 

foodStrings <- sapply(split,function(x) paste(x[1:2],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",foodStrings,perl=T) # clean tips
food_L2 <- rowsum(food[-ncol(food)],foodStrings) 
rownames(food_L2) <- gsub(".*;L2_",'',rownames(food_L2))
rownames(food_L2) <- paste0('food_',rownames(food_L2))
food_L2 <- as.data.frame(t(food_L2))



soylent <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
soylentIDs <- droplevels(soylent$X.SampleID)

cazyme_to_keep <- colnames(cazyme_list_clr[['MCTs01']])
for(i in names(cazyme_list_clr)){
  print(i)
  cazyme_to_keep <- base::intersect(cazyme_to_keep,colnames(cazyme_list_clr[[i]]))
}
cazyme_to_keep <- setdiff(cazyme_to_keep,'X.SampleID')
save(cazyme_to_keep,file='./data/cazyme_to_keep_for_pred.RData')

all_cazyme <- lapply(cazyme_list_clr, function(x){
  y <- melt(x,id.vars='X.SampleID',variable.name = 'cazyme',value.name='relative_abundance')
  y
})

all_cazyme <- do.call('rbind',all_cazyme)
all_cazyme <- dcast(all_cazyme,X.SampleID~cazyme,value.var = 'relative_abundance')
all_cazyme <- all_cazyme[,c(cazyme_to_keep,'X.SampleID')]
all_cazyme <- merge(all_cazyme,map[c('X.SampleID','UserName','StudyDayNo')],by='X.SampleID')
all_cazyme <- all_cazyme[!all_cazyme$X.SampleID%in%soylentIDs,]

cazyme_dist <- dist(select(all_cazyme,-UserName,-StudyDayNo,-X.SampleID))
cazymepc <- as.data.frame(pcoa(cazyme_dist)$vectors)
rownames(cazymepc) <- all_cazyme$X.SampleID

cazyme_predict_dat <- list()
for(i in unique(all_cazyme$UserName)){

  temp_data <- NULL
  #temp_sampleid
  temp_data$X.SampleID.x <- dat[[i]]$X.SampleID.x
  temp_data$X.SampleID.y <- dat[[i]]$X.SampleID.y
  interested_X.SampleID <- Reduce(base::intersect,list(v1=temp_data$X.SampleID.y,
                                                       v2=temp_data$X.SampleID.x,
                                                       v3=all_cazyme$X.SampleID,
                                                       v4=rownames(food_L2)))
  temp_data$X.SampleID.x <- temp_data$X.SampleID.x[temp_data$X.SampleID.x%in%interested_X.SampleID]
  temp_data$X.SampleID.y <- temp_data$X.SampleID.y[temp_data$X.SampleID.y%in%interested_X.SampleID]
  temp_cazyme.x <- all_cazyme[all_cazyme$X.SampleID%in%interested_X.SampleID,] %>% select(-UserName,-StudyDayNo,-X.SampleID)
  temp_cazyme.y <- all_cazyme[all_cazyme$X.SampleID%in%interested_X.SampleID,] %>% select(-UserName,-StudyDayNo,-X.SampleID)
  colnames(temp_cazyme.y) <- paste0(colnames(temp_cazyme.x),'.y')
  colnames(temp_cazyme.x) <- paste0(colnames(temp_cazyme.x),'.x')
  temp_cazymepc <- cazymepc[rownames(cazymepc)%in%interested_X.SampleID,][1:5]
  colnames(temp_cazymepc) <- paste0('cazyme.',colnames(temp_cazymepc))
  temp_food <- food_L2[rownames(food_L2)%in%interested_X.SampleID,]
  temp_data <- cbind(temp_data,temp_cazyme.x,temp_cazyme.y,temp_cazymepc,temp_food)
  if(nrow(temp_data)<=5){
    next()
  }
  print(i)
  cazyme_predict_dat[[i]] <- temp_data
}


save(cazyme_predict_dat,file = './data/cazyme_predict_list.RData')

load(file = './data/cazyme_predict_list.RData')

# cazyme_rf_importance <- list()
# for(i in names(cazyme_predict_dat)){
#   print(i)
#   today.cazymes <- paste0(cazyme_to_keep,'.x')
#   tomorrow.cazymes <- paste0(cazyme_to_keep,'.y')
#   
#   rf.importance.list <- list()
#   for(j in seq_along(tomorrow.cazymes)){
#     today.cazyme <- today.cazymes[j]
#     tomorrow.cazyme <- tomorrow.cazymes[j]
#     cazyme.features <- c("cazyme.Axis.1","cazyme.Axis.2", "cazyme.Axis.3", "cazyme.Axis.4", "cazyme.Axis.5")
#     food.features <- colnames(food_L2)
#     
#     predictor.feats <- cazyme_predict_dat[[i]][,c(cazyme.features,food.features)]
#     result.cazyme.tomorrow <- cazyme_predict_dat[[i]][tomorrow.cazyme]
#     
#     rf.formula <- as.formula(paste0(tomorrow.cazyme,'~',paste(cazyme.features,collapse = '+'),'+',
#                                     paste(food.features,collapse = '+')))
#     
#     model <- randomForest(rf.formula,data = cazyme_predict_dat[[i]],importance=TRUE)
# 
#     rf.importance <- as.data.frame(model$importanceSD)
#     colnames(rf.importance) <- paste0(gsub('.y','',tomorrow.cazyme))
#     rf.importance.list[[tomorrow.cazyme]] <- rf.importance
#   }
#   rf.importance.sub <- do.call('cbind',rf.importance.list)
#   rf.importance.sub[is.na(rf.importance.sub)] <- 0
#   cazyme_rf_importance[[i]] <- rf.importance.sub
# }
# 
# cazyme_rf_importance.all <- do.call('cbind',cazyme_rf_importance)
# cazyme_rf_importance.all <- cazyme_rf_importance.all[apply(cazyme_rf_importance.all,1,sd)!=0,]
# 
# apply(cazyme_rf_importance.all,1,sd)!=0

pdf(file = './result/import_try.pdf',width = 50,height = 5)
pheatmap(cazyme_rf_importance.all,scale = 'row')
dev.off()


#average for each subject
set.seed(666)
cazyme_rf_importance <- list()
cazyme_rf_importance_person <- list()
for(i in names(cazyme_predict_dat)){
  print(i)
  today.cazymes <- paste0(cazyme_to_keep,'.x')
  tomorrow.cazymes <- paste0(cazyme_to_keep,'.y')
  
  rf.importance.list <- list()
  #rf.importance.person.list <- list()
  for(j in seq_along(tomorrow.cazymes)){
    today.cazyme <- today.cazymes[j]
    tomorrow.cazyme <- tomorrow.cazymes[j]
    cazyme.features <- c("cazyme.Axis.1","cazyme.Axis.2", "cazyme.Axis.3", "cazyme.Axis.4", "cazyme.Axis.5")
    food.features <- colnames(food_L2)
    
    predictor.feats <- cazyme_predict_dat[[i]][,c(cazyme.features,food.features)]
    result.cazyme.tomorrow <- cazyme_predict_dat[[i]][tomorrow.cazyme]
    
    rf.formula <- as.formula(paste0(tomorrow.cazyme,'~',paste(cazyme.features,collapse = '+'),'+',
                                    paste(food.features,collapse = '+')))
    a <- cazyme_predict_dat[[i]]
    model <- randomForest(rf.formula,data = cazyme_predict_dat[[i]],importance=F,ntree=800)
    rf.importance <- as.data.frame(model$importance)
    colnames(rf.importance) <- paste0(gsub('.y','',tomorrow.cazyme))
    rf.importance.list[[tomorrow.cazyme]] <- rf.importance
  }
  rf.importance.sub <- do.call('cbind',rf.importance.list)
  rf.importance.sub[is.na(rf.importance.sub)] <- 0
  cazyme_rf_importance_person[[i]] <- rf.importance.sub
  rf.importance.sub <- as.data.frame(rowMeans(rf.importance.sub))
  
  #colnames(rf.importance.sub) <- paste0(i,'.avg.importance')
  colnames(rf.importance.sub) <- i
  cazyme_rf_importance[[i]] <- rf.importance.sub
}

cazyme_rf_importance.all <- do.call('cbind',cazyme_rf_importance)
cazyme_rf_importance.all <- cazyme_rf_importance.all[apply(cazyme_rf_importance.all,1,sd)!=0,]
rownames(cazyme_rf_importance.all) <- gsub('food_','',rownames(cazyme_rf_importance.all))

anno_food <- as.data.frame(rownames(cazyme_rf_importance.all))
rownames(anno_food) <- rownames(cazyme_rf_importance.all)
colnames(anno_food) <- 'category'
anno_food$category <- gsub('cazyme.Axis.*','cazyme',anno_food$category)

food_name <- food$taxonomy

  for(j in seq_along(anno_food$category)){
    if(sum(grepl(anno_food$category[j],food_name))>0){
      #anno_food_name <- c(anno_food_name,sample(food_name[grepl(anno_food$category[j],food_name)],1))
      anno_food$category[j] <- sample(food_name[grepl(anno_food$category[j],food_name)],1)
    }
  }

anno_food$category <- gsub(';L2.*','',anno_food$category)
anno_food$category <- gsub('L1_','',anno_food$category)




pdf(file = './result/import_try_add_group.pdf',width = 15,height = 8)
pheatmap(cazyme_rf_importance.all,scale = 'row', color=colorRampPalette(c("green","black","red"))(100),
         border_color = NA, annotation_row = anno_food, cluster_rows = F)
dev.off()

pdf(file = './result/import_try_add_group_col.pdf',width = 15,height = 8)
pheatmap(cazyme_rf_importance.all,scale = 'column', color=colorRampPalette(c("green","black","red"))(100),
         border_color = NA, annotation_row = anno_food, cluster_rows = T)
dev.off()


