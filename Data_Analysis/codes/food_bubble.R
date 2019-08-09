
require(robCompositions)
require(tibble)
require(psych)
require(reshape2)

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

food <- read.table('./data/masterdecaydiet.txt',header = T)
food$taxonomy <- as.character(food$taxonomy)

split <- strsplit(food$taxonomy,";") 

foodStrings <- sapply(split,function(x) paste(x[1:2],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",foodStrings,perl=T) # clean tips
food_L2 <- rowsum(food[-ncol(food)],foodStrings) 
rownames(food_L2) <- gsub(".*;L2_",'',rownames(food_L2))
#rownames(food_L2) = sapply(strsplit(food$taxonomy,";"),function(x) paste(x[1:2],collapse=";"));

cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)


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

#cazyme_c <- as.data.frame(t(cazyme))
#cazyme_c <- rownames_to_column(cazyme_c,var = 'X.SampleID')
#cazyme_c <- merge(cazyme_c,map[c('X.SampleID','UserName')],by='X.SampleID')

food_c <- as.data.frame(t(food_L2))
food_c <- rownames_to_column(food_c,var='X.SampleID')
food_c <- merge(food_c,map[c('X.SampleID','UserName')],by='X.SampleID')

cazyme_food_list <- list()
load('./data/cazy_list_clr.RData')

for(i in unique(food_c$UserName)){
  print(i)
  coef_v <- c()
  fdr_v <- c()
  food_v <- c()
  cazyme_v <- c()
  p_v <- c()
  id <- c()
  temp_food <- food_c[food_c$UserName==i,!colnames(food_c)%in%c('X.SampleID','UserName')]
  #temp_cazyme <- cazyme_c[cazyme_c$UserName==i,!colnames(cazyme_c)%in%c('X.SampleID','UserName')]
  temp_cazyme <- cazyme_list_clr[[i]]
  temp_cazyme <- temp_cazyme[temp_cazyme$X.SampleID%in%food_c[food_c$UserName==i,'X.SampleID'],]
  temp_cazyme <- temp_cazyme[,!colnames(temp_cazyme)=='X.SampleID']
  temp_food <- temp_food[,colSums(temp_food)!=0]

  #temp_cazyme_pre <- temp_cazyme
  #temp_cazyme <- sweep(temp_cazyme,1,rowSums(temp_cazyme),'/')
  #temp_cazyme_pre[temp_cazyme_pre>=1] <- 1
  #temp_cazyme_pre <- temp_cazyme_pre[,colSums(temp_cazyme_pre)>0.75*nrow(temp_cazyme_pre)]
  #temp_cazyme <- temp_cazyme[,colnames(temp_cazyme)%in%colnames(temp_cazyme_pre)]
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
  #fdr_v <- p.adjust(p_v,method = 'fdr')
  temp_cor_df <- data.frame(food=food_v,cazyme=cazyme_v,coef=coef_v,id=id,fdr_p=fdr_v,p=p_v)
  cazyme_food_list[[i]] <- temp_cor_df
}

save(cazyme_food_list,file = "./data/cazyme_food_cor_ind_0808.RData")


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

allsigs$food <- gsub("Dry_Beans_Peas_Other_Legumes_Nuts_and_Seeds", "Legumes", allsigs$food)
allsigs$food <- gsub("Fats_Oils_and_Salad_Dressings", "Fats", allsigs$food)
allsigs$food <- gsub("Grain_Product", "Grains", allsigs$food)
allsigs$food <- gsub("Milk_and_Milk_Products", "Milks", allsigs$food)
allsigs$food <- gsub("Meat_Poultry_Fish_and_Mixtures", "Meats", allsigs$food)
allsigs$food <- gsub("Sugars_Sweets_and_Beverages", "Sweets and Beverages", allsigs$food)
allsigs$food <- gsub("Cereals_not_cooked_or_NS_as_to_cooked", "Grains", allsigs$food)
allsigs$food <- gsub("Grain_mixtures_frozen_plate_meals_soups", "Grains", allsigs$food)
allsigs$food <- gsub("Nuts_nut_butters_and_nut_mixtures", "Nuts", allsigs$food)
allsigs$food <- gsub("Darkgreen_vegetables", "Vegetables", allsigs$food)
allsigs$food <- gsub("Alcoholic_beverages", "Sweets and Beverages", allsigs$food)
allsigs$food <- gsub("Citrus_fruits_juices", "Fruits", allsigs$food)
allsigs$food <- gsub("Deepyellow_vegetables", "Vegetables", allsigs$food)
allsigs$food <- gsub("Cakes_cookies_pies_pastries_bars", "Grains", allsigs$food)
allsigs$food <- gsub("Crackers_and_salty_snacks_from_grain", "Grains", allsigs$food)
allsigs$food <- gsub("Creams_and_cream_substitutes", "Milks", allsigs$food)
allsigs$food <- gsub("Organ_meats_sausages_and_lunchmeats", "Meats", allsigs$food)
allsigs$food <- gsub("Meatpoultry_fish_with_nonmeat", "Meats", allsigs$food)
allsigs$food <- gsub("Fish_and_shellfish", "Meats", allsigs$food)
allsigs$food <- gsub("Other_vegetables", "Vegetables", allsigs$food)
allsigs$food <- gsub("Other_fruits", "Fruits", allsigs$food)
allsigs$food <- gsub("Formulated_nutrition_beverages_energy_drinks_sports_drinks_function", "Sweets and Beverages", allsigs$food)
allsigs$food <- gsub("Nonalcoholic_beverages", "Sweets and Beverages", allsigs$food)
allsigs$food <- gsub("Pancakes_waffles_French_toast_other", "Grains", allsigs$food)
allsigs$food <- gsub("Fruit_juices_and_nectars_excluding_citrus", "Fruits", allsigs$food)
allsigs$food <- gsub("Frozen_and_shelfstable_plate_meals_soups_and_gravies", "Meats", allsigs$food)
allsigs$food <- gsub("Milk_desserts_sauces_gravies", "Milks", allsigs$food)
allsigs$food <- gsub("Milks_and_milk_drinks", "Milks", allsigs$food)
allsigs$food <- gsub("Lamb_veal_game_other", "Meats", allsigs$food)
allsigs$food <- gsub("Pastas_cooked_cereals_rice", "Grains", allsigs$food)
allsigs$food <- gsub("Quick_breads", "Grains", allsigs$food)
allsigs$food <- gsub("Salad_dressings", "Fats", allsigs$food)
allsigs$food <- gsub("Sugars_and_sweets", "Sweets and Beverages", allsigs$food)
allsigs$food <- gsub("Seeds_and_seed_mixtures", "Legumes", allsigs$food)
allsigs$food <- gsub("Tomatoes_and_tomato_mixtures", "Vegetables", allsigs$food)
allsigs$food <- gsub("Egg_mixtures", "Eggs", allsigs$food)
allsigs$food <- gsub("Dried_fruits", "Fruits", allsigs$food)
allsigs$food <- gsub("Yeast", "Grains", allsigs$food)
allsigs$food <- gsub("White_potatoes_and_Puerto_Rican_starchy_vegetables", "Vegetables", allsigs$food)
allsigs$food <- gsub("Vegetables_with_meat_poultry_fish", "Vegetables", allsigs$food)
allsigs$food <- gsub("Water_noncarbonated", "Sweets and Beverages", allsigs$food)
allsigs$food <- gsub("Poultry", "Meats", allsigs$food)
allsigs$food <- gsub("Pork", "Meats", allsigs$food)
allsigs$food <- gsub("Cheeses", "Milks", allsigs$food)
allsigs$food <- gsub("Beef", "Meats", allsigs$food)
allsigs$food <- gsub("Fruits_and_juices_baby_food", "Fruits", allsigs$food)
allsigs$food <- gsub("Meat_substitutes_mainly_cereal_protein", "Meats", allsigs$food)

#write.csv(allsigs,"./data/network/food_cazyme_cor_try.csv",row.names = F,quote = F)

# load colors
source(file = "G:/Dan_Lab/dietstudy_analyses-master/lib/colors/UserNameColors.R")
UserNameColors['MCTs05'] <- '#99ff00'

## get information to label key values repeated in more than one person for labeling
allsigs_names <- allsigs[colnames(allsigs) %in% c("food", "cazyme")]
allsigs_names_pairs <- paste(allsigs_names$food, allsigs_names$cazyme)

names_repeated<-as.data.frame(table(allsigs_names_pairs))
names_repeated<- subset(names_repeated, names_repeated$Freq >2)
table(names_repeated$Freq)
names_repeated <- colsplit(names_repeated$allsigs_names_pairs, " ", c("food", "cazyme"))
names_repeated$label <- "*"




allsigs <- merge(allsigs, names_repeated, all.x = T)

allsigs$id <- gsub("MCTs", "", allsigs$id)

table(allsigs$id)

length(table(allsigs$id))


names(UserNameColors) <- gsub("MCTs", "", names(UserNameColors))
pal <- rev(c("#ff40f2", "#ff0000", "#008c4b", "#00138c", "#8c235b", "#ffbfbf", "#8c7723", "#468c75", "#8091ff", "#ff80c4", "#8c3123", "#fff2bf", "#40fff2", "#69698c", "#ff0044", "#ff9180", "#e5ff80", "#bffbff", "#5940ff", "#8c696e", "#8c7369", "#858c69", "#40d9ff", "#c480ff", "#ff8c40", "#4b8c00", "#23698c", "#69238c", "#8c4b00", "#bfffbf", "#004b8c", "#eabfff", "#ffc480", "#40ff59", "#80c4ff", "#ffd940" ))

require(ggplot2)


allsigs <- allsigs[allsigs$fdr_p>0.00001,]

myplot <- ggplot(data = allsigs, aes(x = coef, y = cazyme, size = -log(fdr_p), color = id)) +
  geom_point(alpha = 0.6) +
  #geom_point(alpha = 0.8, color = "darkgrey", pch = 21) +
  facet_grid(food~bin, scales = "free", space = "free_y")+
  scale_color_manual(values = UserNameColors) +
  theme_classic() +
  guides(color = guide_legend(nrow = 10, title = "Subject", title.position = "top"),
         size = guide_legend(title.position = "top", title = "-log(FDR p-value)")) +
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9, color = "black"),
        panel.grid.major = element_line(colour = "lightgrey"),
        strip.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  ylab("Cazyme") +
  xlab("Spearman correlation") +
  scale_x_continuous(trans = "reverse")


myplot

ggsave('./result/food_cor_imp_ind_fdr0.3.pdf',height = 20,width = 8,limitsize = F)
