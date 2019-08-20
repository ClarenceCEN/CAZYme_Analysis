require(rmarkdown)
require(knitr)
require(tidyverse)
require(RColorBrewer) 
require(cowplot)
require(reshape2)
require(ggdendro)
require(vegan)
require(ape)

setwd("G:/Dan_Lab/dietstudy_analyses-master/")

map_sample <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
map_username <- read.table("data/maps/UserName_map.txt", sep = "\t", header = T, comment = "")
map_sample$StudyDate <- as.Date.factor(map_sample$StudyDate, format = "%m/%d/%y")

# identify soylent samples
soylent <- map_sample[map_sample$UserName %in% c("MCTs11", "MCTs12"),]
soylent <- as.character(soylent$X.SampleID)

##########food##############
# load the food distance matrix, unweighted unifrac
food_un <- read.delim("data/diet/processed_food/dhydrt_smry_no_soy_beta/unweighted_unifrac_dhydrt.smry.no.soy.txt", row = 1) # weighted in not significant
food_dist <- as.dist(food_un)

# make non-tree distance matrix for food (for supplemental)
food <- read.delim("data/diet/processed_food/dhydrt.smry.no.soy.txt", row = 1)
food <- food[,colnames(food) %in% colnames(food_un)]
no_tree_dist <- dist(t(food))

###############nutrients############
# load nutrition data
nutr <- read.delim("data/diet/processed_nutr/nutr_65_smry_no_soy.txt", row = 1)

# normalize nutrition data across features (rows)
nutr_n <- sweep(nutr, 1, rowSums(nutr), "/")

# make nutrition distance matrix (euclidean)
nutr_dist <- dist(t(nutr_n))


setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
cazyme <- cazyme[,!colnames(cazyme)%in%soylent]


map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
cazyme_filter <- as.data.frame(t(cazyme))
cazyme_filter <- rownames_to_column(cazyme_filter,var = 'X.SampleID')
cazyme_filter <- merge(cazyme_filter,map[c('X.SampleID','UserName','StudyDayNo')],by = 'X.SampleID')

cazyme_filter <- cazyme_filter[,-grep('Other',colnames(cazyme_filter))]
colnames(cazyme_filter) <- gsub(".*;L4_",'',colnames(cazyme_filter))


#jaccard
cazyme_count <- cazyme_filter[,!colnames(cazyme_filter)%in%c('X.SampleID','StudyDayNo','UserName')]
cazyme_count[cazyme_count>0] <- 1
cazyme_present <- aggregate(cazyme_count,by = list(cazyme_filter$UserName),FUN = sum)
rownames(cazyme_present) <- cazyme_present$Group.1;cazyme_present <- cazyme_present[-1]
cazyme_present[cazyme_present<3] <- 0
cazyme_present[cazyme_present>=2] <- 1


cazyme_dist <- vegdist(cazyme_present,method = 'jaccard')
cazyme_dist <- dist(cazyme_present)

#Aitchison
load('./data/cazy_list_clr.RData')
cazyme_mean_list <- list()
for(i in names(cazyme_list_clr)){
  temp <- cazyme_list_clr[[i]]
  temp <- select(temp,-X.SampleID)
  temp$UserName <- i
  temp_mean <- aggregate(temp[,colnames(temp)!='UserName'],by=list(temp$UserName),FUN=mean)
  temp_mean_long <- melt(temp_mean, id.vars='Group.1',variable.name = 'cazyme',value.name = 'relative_abundance_clr')
  cazyme_mean_list[[i]] <- temp_mean_long
}
cazyme_mean <- do.call('rbind',cazyme_mean_list) %>% mutate_if(is.factor,as.character)
cazyme_mean <- dcast(cazyme_mean,Group.1~cazyme)
cazyme_mean[is.na(cazyme_mean)] <- 0
cazyme_mean <- cazyme_mean[!cazyme_mean$Group.1%in%c("MCTs11", "MCTs12"),]
cazyme_dist <- dist(cazyme_mean[-1])




# make pcoa
cazyme_all_pcoa <- data.frame(pcoa(cazyme_all_dist)$vectors)

eigen <- pcoa(cazyme_all_dist)$values$Eigenvalues

percent_var <- signif(eigen/sum(eigen), 4)*100

# move rownames to a column
cazyme_all_pcoa <- rownames_to_column(cazyme_all_pcoa, var = "X.SampleID")

#merge map and PCOA by SampleID - select both lines and run together - won't work otherwise
cazyme_all_pcoa <- inner_join(cazyme_all_pcoa, map_sample, by = 'X.SampleID')

## PERMANOVA for within subject grouping

mytest <- adonis(cazyme_all_dist~cazyme_all_pcoa$UserName, permutations = 999)
plot_p <- mytest$aov.tab$`Pr(>F)`[1]

mytest


cazyme_all_plot_one_color <- ggplot(cazyme_all_pcoa, aes(x = Axis.1, y = Axis.2, fill = UserName)) +
  geom_point(color =  "#5d47ff", alpha=0.75, size = 2, alpha = 0.6) +
  stat_ellipse(color = "dark grey", type = "norm", linetype = 2, level = 0.95, alpha = 0.5) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) +
  theme_classic() +                                                                      
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 4),
        axis.title = element_text(size=9)) +
  guides(color = guide_legend(ncol = 4)) +
annotate("text", x = 40, y = -50, label = paste0("p-value = ",plot_p))

cazyme_all_plot_one_color + theme(legend.position = 'none')
ggsave('./result/ado_test.pdf')

#### MAKE 2E #####

# make pcoas 
pcoa_f <- as.data.frame(pcoa(food_dist)$vectors)
pcoa_c <- as.data.frame(pcoa(cazyme_dist)$vectors)

# procrustes
pro <- procrustes(pcoa_f, pcoa_c)
pro_test <- protest(pcoa_f, pcoa_c, perm = 999)
mantel_test <- mantel(food_dist,cazyme_dist,method = 'spearman')

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Food (Unweighted Unifrac)"
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Cazyme (Aitchison's)"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif, 1)

plot <- rbind(beta_pro, trans_pro)

food_cazyme <- ggplot(plot) +
  geom_point(size = 3, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#ff0f17", "#ffb405")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.3, y = -0.27, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 


food_cazyme_leg <- get_legend(food_cazyme)


food_cazyme + theme(legend.position = "none")

ggsave('./result/pro_test_selected.pdf')


# make pcoas 
pcoa_n <- as.data.frame(pcoa(nutr_dist)$vectors)
pcoa_c <- as.data.frame(pcoa(cazyme_dist)$vectors)

# procrustes
pro <- procrustes(pcoa_n, pcoa_c)
pro_test <- protest(pcoa_n, pcoa_c, perm = 999)
mantel_test <- mantel(nutr_dist,cazyme_dist,method = 'spearman')

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Nutrition (Euclidean)"
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Cazyme (Aitchison's)"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif, 1)

plot <- rbind(beta_pro, trans_pro)

nutr_cazyme <- ggplot(plot) +
  geom_point(size = 3, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#1a661a", "#00becc")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.2, y = -0.2, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 


nutr_cazyme_leg <- get_legend(nutr_cazyme)


nutr_cazyme + theme(legend.position = "none")

ggsave('./result/pro_test_nutr_try.pdf')

# Import the beta diversity table from grains_beta (weighted or unweighted)
setwd("G:/Dan_Lab/dietstudy_analyses-master/")
food_beta <- read.table(file="data/diet/fiber/grains_data/grains_beta/unweighted_unifrac_grains_fiber.txt")
cazyme_pro <- cazyme_pro[,colnames(cazyme_pro)%in%colnames(food_beta)]
food_dist <- as.dist(food_beta)
cazyme_dist <- dist(t(cazyme_pro))


pcoa_f <- as.data.frame(pcoa(food_dist)$vectors)
pcoa_c <- as.data.frame(pcoa(cazyme_dist)$vectors)

# procrustes
pro <- procrustes(pcoa_f, pcoa_c)
pro_test <- protest(pcoa_f, pcoa_c, perm = 999)
mantel_test <- mantel(food_dist,cazyme_dist,method = 'spearman')

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Food (Unweighted Unifrac)"
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Cazyme (Aitchison's)"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif, 1)

plot <- rbind(beta_pro, trans_pro)

food_cazyme <- ggplot(plot) +
  geom_point(size = 3, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#ff0f17", "#ffb405")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.3, y = -0.27, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 


food_cazyme_leg <- get_legend(food_cazyme)


food_cazyme + theme(legend.position = "none")

ggsave('./result/pro_test_try.pdf')