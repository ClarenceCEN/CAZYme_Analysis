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


setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')
load('./data/cazyme_mean_clr.RData')
cazyme <- cazyme_mean_clr_trans
cazyme <- cazyme[,colnames(cazyme)%in%colnames(food_un)]
cazyme_dist <- dist(t(cazyme))


load('./data/cazyme_clr.RData')
cazyme_all <- as.data.frame(cazyme_clr_trans)
cazyme_all <- cazyme_all[,!colnames(cazyme_all)%in%soylent]
cazyme_all_dist <- dist(t(cazyme_all))

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

ggsave('./result/pro_test.pdf')


