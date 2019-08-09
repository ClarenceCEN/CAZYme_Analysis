load("./data/cazyme_food_cor_ind_0808.RData")
load('./data/cazy_list_clr.RData')
require(robCompositions)
require(tibble)
require(psych)
require(reshape2)
require(dplyr)

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

food_c <- as.data.frame(t(food_L2))
food_c <- rownames_to_column(food_c,var='X.SampleID')
food_c <- merge(food_c,map[c('X.SampleID','UserName')],by='X.SampleID')


sigs <- lapply(cazyme_food_list, function(x) subset(x, fdr_p <= 0.1))
allsigs <- do.call("rbind", sigs)

allsigs$cazy_cat <- ifelse(grepl('AA',allsigs$cazyme),'AA',
                           ifelse(grepl('CBM',allsigs$cazyme),'CBM',
                                  ifelse(grepl('GT',allsigs$cazyme),'GT',
                                         ifelse(grepl('GH',allsigs$cazyme),'GH',
                                                ifelse(grepl('PL',allsigs$cazyme),'PL','CE')))))

allsigs$bin <- ifelse(allsigs$coef < 0, "Negative (+/-)", 
                      ifelse(allsigs$coef > 0, "Positive (+/+ or -/-)", "NA"))
allsigs$bin <- factor(allsigs$bin, levels = c("Positive (+/+ or -/-)", "Negative (+/-)"))


#numbers of cazyme
length(unique(allsigs$cazyme))  #43

for(i in unique(allsigs$cazyme)){
  temp <- allsigs[allsigs$cazyme==i,]
  #print(paste(i,length(unique(temp$id))))
  if(length(unique(temp$id))>1){
    print(paste(i,length(unique(temp$id)),temp$id,temp$food))
  }

}

allsigs <- mutate_if(allsigs,is.factor,as.character)

#GH88
sub_cazyme_id <- allsigs[allsigs$cazyme=="GH88",]$id
sub_cazyme_food <- allsigs[allsigs$cazyme=="GH88",]$food

sub_cazyme <- cazyme_list_clr[sub_cazyme_id]
sub_food <- food_c[food_c$UserName%in%names(sub_cazyme),c("X.SampleID","UserName",sub_cazyme_food)]
sub_cazyme <- lapply(sub_cazyme, function(x){select(x,c("GH88",'X.SampleID'))})
sub_cazyme <- do.call('rbind',sub_cazyme)

sub_plot <- merge(sub_cazyme,sub_food,by='X.SampleID')

sub_plot <- melt(sub_plot,id.vars = c("X.SampleID","UserName",'GH88'),variable.name = 'food',
                 value.name = 'weight')


g <- ggplot(sub_plot,aes(x=weight,y=GH88)) + geom_point(aes(color=UserName),size=4) +
  facet_grid(~UserName+food,scales = 'free') +theme_classic() + stat_smooth(method = lm,se=FALSE) +
  scale_color_manual(values = UserNameColors)
g
