require(robCompositions)
require(tibble)
require(psych)
require(reshape2)
require(dplyr)

load("./data/cazyme_food_cor_ind_0808.RData")
load('./data/cazy_list_clr.RData')
load("./data/food_daily_L2.RData")

cazyme_food_cor <- do.call('rbind',cazyme_food_list) %>% filter(fdr_p<0.2)
cazyme_food_cor <- mutate_if(cazyme_food_cor,is.factor,as.character)


cazyme_food_cor$cazy_cat <- ifelse(grepl('AA',cazyme_food_cor$cazyme),'AA',
                           ifelse(grepl('CBM',cazyme_food_cor$cazyme),'CBM',
                                  ifelse(grepl('GT',cazyme_food_cor$cazyme),'GT',
                                         ifelse(grepl('GH',cazyme_food_cor$cazyme),'GH',
                                                ifelse(grepl('PL',cazyme_food_cor$cazyme),'PL','CE')))))

food <- read.table('./data/masterdecaydiet.txt',header = T)
food_names <- as.character(food$taxonomy)

#cazyme_food_cor$L1 <- gsub(pastecazyme_food_cor$food)
#grep(cazyme_food_cor$food,food_names)

full_name <- sapply(cazyme_food_cor$food,function(x){
  full_name <- unique(food_names[grep(x,paste0(";L2_",food_names))])[1]
  full_name
},simplify = T)

cazyme_selected <- c()
cazyme_selected_list <- list()
for(i in unique(cazyme_food_cor$cazyme)){
  temp <- cazyme_food_cor[cazyme_food_cor$cazyme==i,]
  if(nrow(temp)>1){
    cazyme_selected <- c(cazyme_selected,i)
    cazyme_selected_list[[i]] <- temp %>% mutate_if(is.factor,as.character)
  }
}
