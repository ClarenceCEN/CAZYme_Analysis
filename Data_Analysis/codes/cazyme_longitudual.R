require(rmarkdown)
require(knitr)
require(tidyverse)
require(RColorBrewer) 
require(cowplot)
require(reshape2)
require(ggdendro)
require(vegan)
require(gridExtra)
require(pheatmap)




setwd('G:/Dan_Lab/codes/CAZyme/CAZYme_Analysis/Data_Analysis/')

load("./codes/order.RData")

cazyme <- read.table('./data/Cazyme_total.txt',sep='\t',header = T,row.names = 1)
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
pal <- rev(c("#ff40f2", "#ff0000", "#008c4b", "#00138c", "#8c235b", "#ffbfbf", "#8c7723", "#468c75", "#8091ff", "#ff80c4", "#8c3123", "#fff2bf", "#40fff2", "#69698c", "#ff0044", "#ff9180", "#e5ff80", "#bffbff", "#5940ff", "#8c696e", "#8c7369", "#858c69", "#40d9ff", "#c480ff", "#ff8c40", "#4b8c00", "#23698c", "#69238c", "#8c4b00", "#bfffbf", "#004b8c", "#eabfff", "#ffc480", "#40ff59", "#80c4ff", "#ffd940" ))


rownames(cazyme) <- gsub(".*;L4_",'',rownames(cazyme))


#draw the heatmap
heat1 <- as.data.frame(t(cazyme))
heat1 <- rownames_to_column(heat1,var='X.SampleID')
heat1 <- merge(heat1,map,by = 'X.SampleID')

cazyme_family <- rownames(cazyme)

for(i in unique(heat1$UserName)){
  print(i)
  heat_temp <- heat1[heat1$UserName==i,] %>% select('StudyDayNo',cazyme_family)
  rownames(heat_temp) <- heat_temp$StudyDayNo
  heat_temp <- heat_temp[-1]
  heat_temp[heat_temp>0] <- 1
  heat_temp <- as.data.frame(t(heat_temp))
  pdf(paste0('./result/',i,'.pdf'),width = 5,height = 12)
  pheatmap(heat_temp,cluster_cols = F)
  dev.off()
}














#nomalization
cazyme <- sweep(cazyme,2,colSums(cazyme),'/')

#order by mean value
cazyme <- cazyme[order(rowMeans(cazyme),decreasing = F),]

#generate plotting table
plot1 <- as.data.frame(t(cazyme))
plot1 <- rownames_to_column(plot1,var='X.SampleID')
plot1 <- melt(plot1,id = 'X.SampleID',variable.name = 'cazyme')
plot1 <- merge(plot1,map,by = 'X.SampleID')

# get right number of colors for plotting
no_cols <- length(unique(plot1$cazyme))
#colors_func <- sample(cols_func(no_cols))
colors_func <- sample(pal, no_cols,replace = T)
mycols <- sample(colors(),no_cols)   

#change the order
plot1$UserName <- factor(plot1$UserName, levels = ord_factor)


cazyme_plot <- ggplot(data = plot1, aes(x=StudyDayNo, y = value, fill=cazyme)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = mycols) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_blank(), 
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        panel.spacing.x=unit(0.05, "lines")) +
  guides(fill = guide_legend(reverse = TRUE, 
                             keywidth = 0.5, 
                             keyheight = 0.5, 
                             nrow = 5)) +
  #nrow = 5)) +
  ylab("Relative Abundance") 

cazyme_plot

ggsave("./result/cazy.pdf",width=20,height=10,dpi = 600)



