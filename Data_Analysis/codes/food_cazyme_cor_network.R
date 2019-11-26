load('./data/cazyme_food_cor.RData')
load('./data/cazyme_food_cor_ind_0815.RData')

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


#for network
allsigs$pairs <- paste0(allsigs$food,'_',allsigs$cazyme)
allsigs <- arrange(allsigs,allsigs$fdr_p)
allsigs$detect_dup <- paste0(allsigs$food,'_',allsigs$cazyme,"_",allsigs$id)
allsigs <- allsigs[!duplicated(allsigs$detect_dup),]

write.csv(allsigs,"./data/network/food_cazyme_cor_0.2.csv",row.names = F,quote = F)

cazyme_cate <- c()
food_cazyme_pair <- c()
allsigs <- mutate_if(allsigs,is.factor,as.character)
for(i in unique(allsigs$pairs)){
  print(i)
  temp_cazyme <- allsigs[allsigs$pairs==i,'cazy_cat'][1]
  food_cazyme_pair <- c(food_cazyme_pair,i)
  cazyme_cate <- c(cazyme_cate,temp_cazyme)
  print(temp_cazyme)
  if(length(food_cazyme_pair)!=length(cazyme_cate)){
    break()
  }
}

subjects <- c()
subjects_type <- c()
for(i in unique(allsigs$id)){
  print(i)
  subjects <- c(subjects,i)
  subjects_type <- c(subjects_type,"subjects")

}

nodes <- c(food_cazyme_pair,subjects)
nodes_type <- c(cazyme_cate,subjects_type)

output.node <- data.frame(nodes=nodes, type=nodes_type)
write.table(output.node, "./data/network/output_node.txt", quote=FALSE, sep="\t", row.names=FALSE)



load('./data/cazyme_food_cor_L3.RData')
sigs <- lapply(cazyme_food_list, function(x) subset(x, fdr_p <= 0.05))
allsigs <- do.call("rbind", sigs)

allsigs$cazy_cat <- ifelse(grepl('AA',allsigs$cazyme),'AA',
                           ifelse(grepl('CBM',allsigs$cazyme),'CBM',
                                  ifelse(grepl('GT',allsigs$cazyme),'GT',
                                         ifelse(grepl('GH',allsigs$cazyme),'GH',
                                                ifelse(grepl('PL',allsigs$cazyme),'PL','CE')))))

allsigs$bin <- ifelse(allsigs$coef < 0, "Negative (+/-)", 
                      ifelse(allsigs$coef > 0, "Positive (+/+ or -/-)", "NA"))
allsigs$bin <- factor(allsigs$bin, levels = c("Positive (+/+ or -/-)", "Negative (+/-)"))

#for network
allsigs$pairs <- paste0(allsigs$food,'_',allsigs$cazyme)
write.csv(allsigs,"./data/network/food_cazyme_cor_L3.csv",row.names = F,quote = F)
