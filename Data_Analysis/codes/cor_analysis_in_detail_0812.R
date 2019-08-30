require(robCompositions)
require(tibble)
require(psych)
require(reshape2)
require(dplyr)
require(ggplot2)

load("./data/cazyme_food_cor_ind_0815.RData")
load('./data/cazy_list_clr.RData')
load("./data/food_daily_L2.RData")

source(file = "G:/Dan_Lab/dietstudy_analyses-master/lib/colors/UserNameColors.R")
UserNameColors['MCTs05'] <- '#99ff00'

cazyme_food_cor <- do.call('rbind',cazyme_food_list) %>% filter(fdr_p<0.2)
cazyme_food_cor <- mutate_if(cazyme_food_cor,is.factor,as.character)


cazyme_food_cor$cazy_cat <- ifelse(grepl('AA',cazyme_food_cor$cazyme),'AA',
                                   ifelse(grepl('CBM',cazyme_food_cor$cazyme),'CBM',
                                          ifelse(grepl('GT',cazyme_food_cor$cazyme),'GT',
                                                 ifelse(grepl('GH',cazyme_food_cor$cazyme),'GH',
                                                        ifelse(grepl('PL',cazyme_food_cor$cazyme),'PL','CE')))))


id_count <- as.data.frame(table(cazyme_food_cor$id))
colnames(id_count) <- c('id','freq')
id_count <- arrange(id_count,desc(freq));id_count$id <- as.character(id_count$id)

g <- ggplot(id_count,aes(x=reorder(id,desc(freq)),y=freq,fill=id)) + geom_bar(stat = 'identity') + theme_classic()+
  scale_fill_manual(values = UserNameColors) +
  geom_text(aes(label=freq),vjust=-0.2) +
  theme(axis.title.x = element_blank())
g
ggsave('./result/cor_id_bar.jpg',width = 12,height = 6)

cazyme_count <- as.data.frame(table(cazyme_food_cor$cazyme))
colnames(cazyme_count) <- c('cazyme','freq')
cazyme_count$cat <- ifelse(grepl('AA',cazyme_count$cazyme),'AA',
                           ifelse(grepl('CBM',cazyme_count$cazyme),'CBM',
                                  ifelse(grepl('GT',cazyme_count$cazyme),'GT',
                                         ifelse(grepl('GH',cazyme_count$cazyme),'GH',
                                                ifelse(grepl('PL',cazyme_count$cazyme),'PL','CE')))))
cazyme_count <- as.data.frame(table(cazyme_count$cat))
colnames(cazyme_count) <- c('cazyme','freq')

g <- ggplot(cazyme_count,aes(x=reorder(cazyme,desc(freq)),y=freq,fill=cazyme)) + geom_bar(stat = 'identity') + theme_classic()+
  #scale_fill_manual(values = UserNameColors) +
  geom_text(aes(label=freq),vjust=-0.2) +
  theme(axis.title.x = element_blank()) 
  #facet_grid(cat~.,scales='free') +
  #geom_segment(aes(yend=cazyme),xend=0,color='grey50')
  
g
ggsave('./result/cor_cazyme_bar.jpg',width = 6,height = 6)

#for network
cazyme_food_cor$bin <- ifelse(cazyme_food_cor$coef < 0, "Negative (+/-)", 
                      ifelse(cazyme_food_cor$coef > 0, "Positive (+/+ or -/-)", "NA"))
cazyme_food_cor <- filter(cazyme_food_cor,fdr_p<0.1)
cazyme_food_cor$food <- gsub("Dry_Beans_Peas_Other_Legumes_Nuts_and_Seeds", "Legumes", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Fats_Oils_and_Salad_Dressings", "Fats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Grain_Product", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Milk_and_Milk_Products", "Milks", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Meat_Poultry_Fish_and_Mixtures", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Sugars_Sweets_and_Beverages", "Sweets and Beverages", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Cereals_not_cooked_or_NS_as_to_cooked", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Grain_mixtures_frozen_plate_meals_soups", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Nuts_nut_butters_and_nut_mixtures", "Nuts", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Darkgreen_vegetables", "Vegetables", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Alcoholic_beverages", "Sweets and Beverages", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Citrus_fruits_juices", "Fruits", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Deepyellow_vegetables", "Vegetables", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Cakes_cookies_pies_pastries_bars", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Crackers_and_salty_snacks_from_grain", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Creams_and_cream_substitutes", "Milks", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Organ_meats_sausages_and_lunchmeats", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Meatpoultry_fish_with_nonmeat", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Fish_and_shellfish", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Other_vegetables", "Vegetables", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Other_fruits", "Fruits", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Formulated_nutrition_beverages_energy_drinks_sports_drinks_function", "Sweets and Beverages", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Nonalcoholic_beverages", "Sweets and Beverages", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Pancakes_waffles_French_toast_other", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Fruit_juices_and_nectars_excluding_citrus", "Fruits", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Frozen_and_shelfstable_plate_meals_soups_and_gravies", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Milk_desserts_sauces_gravies", "Milks", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Milks_and_milk_drinks", "Milks", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Lamb_veal_game_other", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Pastas_cooked_cereals_rice", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Quick_breads", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Salad_dressings", "Fats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Sugars_and_sweets", "Sweets and Beverages", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Seeds_and_seed_mixtures", "Legumes", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Tomatoes_and_tomato_mixtures", "Vegetables", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Egg_mixtures", "Eggs", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Dried_fruits", "Fruits", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Yeast", "Grains", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("White_potatoes_and_Puerto_Rican_starchy_vegetables", "Vegetables", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Vegetables_with_meat_poultry_fish", "Vegetables", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Water_noncarbonated", "Sweets and Beverages", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Poultry", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Pork", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Cheeses", "Milks", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Beef", "Meats", cazyme_food_cor$food)
cazyme_food_cor$food <- gsub("Fruits_and_juices_baby_food", "Fruits", cazyme_food_cor$food)

write.table(cazyme_food_cor,'./data/cor_network_L2.txt',sep = '\t',quote = F,row.names = F)
