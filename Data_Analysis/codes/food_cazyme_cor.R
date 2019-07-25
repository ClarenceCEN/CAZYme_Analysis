
require(robCompositions)
require(tibble)
require(psych)

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

cazyme_food_list <- list()
food_list <- list()
cazyme_list <- list()

for(i in unique(food_c$UserName)){
  coef_v <- c()
  fdr_v <- c()
  food_v <- c()
  cazyme_v <- c()
  temp_food <- food_c[food_c$UserName==i,!colnames(food_c)%in%c('X.SampleID','UserName')]
  temp_cazyme <- cazyme_c[cazyme_c$UserName==i,!colnames(cazyme_c)%in%c('X.SampleID','UserName')]
  temp_cor_r <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$r
  temp_cor_p <- corr.test(temp_food,temp_cazyme,adjust = 'fdr',method = 'spearman')$p
  for(a in 1:nrow(temp_cor_p)){
    for (b in 1:ncol(temp_cor_p)){
      if(!is.na(temp_cor_p[a,b])){
        coef_v <- c(coef_v,temp_cor_r)
        fdr_v <- c(fdr_v,temp_cor_p)
        food_v <- c(food_v,colnames(temp_food)[a])
        cazyme_v <- c(cazyme_v,colnames(temp_cazyme)[b])
      }
    }
  }
  
}
