require(vegan)

indiv_prodecay <- as.list(NULL)
pval.list.decay <- as.list(NULL)

setwd("G:/Dan_Lab/dietstudy_analyses-master")

map_sample <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

set.seed(7)
source("lib/results_scripts/4a_Microbiome and food pair longitudinally/1_Figure3/load_data_lists.R")

setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')
load('./data/cazy_list_clr_2.RData')

cazyme.pcoa <- list()
for(i in names(diet.pcoa.decay)){
  temp_cazy <- cazyme_list_clr[[i]]
  rownames(temp_cazy) <- temp_cazy$X.SampleID
  temp_X.Sample <- intersect(temp_cazy$X.SampleID,rownames(diet.pcoa.decay[[i]]))
  temp_cazy <- temp_cazy[temp_cazy$X.SampleID%in%temp_X.Sample,] %>% select(-X.SampleID)
  diet.pcoa.decay[[i]] <- diet.pcoa.decay[[i]][rownames(diet.pcoa.decay[[i]])%in%temp_X.Sample,]
  temp_pcoa <- as.data.frame(pcoa(dist(temp_cazy))$vectors)
  cazyme.pcoa[[i]] <- temp_pcoa
}


source(file = "G:/Dan_Lab/dietstudy_analyses-master/lib/colors/UserNameColors.R")
UserNameColors['MCTs05'] <- '#99ff00'

for (n in names(cazyme.pcoa)) {
  
  
  pcoa_f <- diet.pcoa.decay[[n]]
  pcoa_c <- cazyme.pcoa[[n]]
  
  # procrustes
  pro <- procrustes(pcoa_c, pcoa_f)
  pro_test <- protest(pcoa_c, pcoa_f, perm = 999) # may want to skip this for plotting?
  
  eigen <- sqrt(pro$svd$d)
  percent_var <- signif(eigen/sum(eigen), 4)*100
  
  beta_pro <- data.frame(pro$X)
  trans_pro <- data.frame(pro$Yrot)
  beta_pro$UserName <- rownames(beta_pro)
  beta_pro$type <- "Food Distance (Unweighted Unifrac)"
  trans_pro$UserName <- rownames(trans_pro)
  trans_pro$type <- "Cazyme Distance (Aitchison's)"
  
  colnames(trans_pro) <- colnames(beta_pro)
  
  pval <- pro_test$signif
  
  plot <- rbind(beta_pro, trans_pro)
  
  plot$X.SampleID <- plot$UserName
  plot$UserName <- unique(map_sample[map_sample$X.SampleID%in%rownames(plot),'UserName'])
  
  
  indiv_prodecay[[n]] <- ggplot(plot) +
    geom_point(shape=21,size = 3, alpha=0.75, aes(x = X.SampleID, y = Axis.1, fill = type, color=UserName)) +
    scale_fill_manual(values = c("#ff0f17", "#ffb405")) +
    scale_color_manual(values = UserNameColors) +
    theme_bw() +
    geom_line(aes(x = X.SampleID, y = Axis.1, group=X.SampleID), col = "darkgrey", alpha = 0.6) +
    labs(title=gsub("MCTs", "", n),caption = paste0("pval=",pval)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #plot.margin = margin(0.02, .01, .01, .01, "in"),
          plot.title = element_text(size = 9, hjust = 0.5),
          legend.position = 'none',
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1,
          plot.caption = element_text(size = (rel(0.8)))) +
    
    #annotate("text", x = max(rownames(plot)), y = min(plot$Axis.1), label = paste0("pval=",pval), size = 2) +
    NULL
  
  pval.list.decay[[n]] <- pval
  
}

# sort for plotting A #
new_order <- names(sort(unlist(pval.list.decay[pval.list.decay < 0.05])))

order_indiv_prodecay <- indiv_prodecay[new_order]

new_order

names(order_indiv_prodecay)



### PLOT 3A  fig.height=2.7, fig.width=7 ####

grid.arrange(grobs = order_indiv_prodecay, nrow=3)

fig <- arrangeGrob(grobs = order_indiv_prodecay, nrow=3)

ggsave("./result/personal_pro.pdf", fig, width = 12, height = 6)
