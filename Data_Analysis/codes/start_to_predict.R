require(reshape2)


setwd("G:/Dan_Lab/dietstudy_analyses-master/")
load(file = "data/test_personal_diet_dat.Rdata")

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
load(file = './data/cazy_list_clr_2.RData')

cazyme_to_keep <- colnames(cazyme_list_clr[['MCTs01']])
for(i in names(cazyme_list_clr)){
  print(i)
  cazyme_to_keep <- base::intersect(cazyme_to_keep,colnames(cazyme_list_clr[[i]]))
}

all_cazyme <- lapply(cazyme_list_clr, function(x){
  y <- melt(x,id.vars='X.SampleID',variable.name = 'cazyme',value.name='relative_abundance')
  y
})

all_cazyme <- do.call('rbind',all_cazyme)
all_cazyme <- dcast(all_cazyme,X.SampleID~cazyme,value.var = 'relative_abundance')
all_cazyme <- all_cazyme[,cazyme_to_keep]
all_cazyme <- merge(all_cazyme,map[c('X.SampleID','UserName','StudyDayNo')],by='X.SampleID')


cazyme_predict_dat <- list()
for(i in unique(all_cazyme$UserName)){
  print(i)
  temp_data <- NULL
  #temp_sampleid
  temp_data$X.SampleID.x <- dat[[i]]$X.SampleID.x
  temp_data$X.SampleID.y <- dat[[i]]$X.SampleID.y
  temp_data$X.SampleID.x <- temp_data$X.SampleID.x[temp_data$X.SampleID.x%in%all_cazyme$X.SampleID]
  temp_data$X.SampleID.y <- temp_data$X.SampleID.y[temp_data$X.SampleID.y%in%all_cazyme$X.SampleID]
  temp_cazyme.x <- all_cazyme[all_cazyme$X.SampleID%in%dat[[i]]$X.SampleID.x,] %>% select(-UserName,-StudyDayNo,-X.SampleID)
  temp_cazyme.y <- all_cazyme[all_cazyme$X.SampleID%in%dat[[i]]$X.SampleID.y,] %>% select(-UserName,-StudyDayNo,-X.SampleID)
  colnames(temp_cazyme.y) <- paste0(colnames(temp_cazyme.x),'.y')
  colnames(temp_cazyme.x) <- paste0(colnames(temp_cazyme.x),'.x')
  temp_data <- cbind(temp_data,temp_cazyme.x,temp_cazyme.y)
  if(nrow(temp_data)==0){
    next()
  }
  cazyme_predict_dat[[i]] <- temp_data
}
