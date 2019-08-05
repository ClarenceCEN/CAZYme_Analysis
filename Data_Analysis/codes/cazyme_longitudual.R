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

load("./data/order.RData")

cazyme <- read.table('./data/Cazyme_total2.txt',sep='\t',header = T,row.names = 1)
#cazyme <- read.table('./data/Cazyme_total_try.txt',sep='\t',header = T,row.names = 1)
map <- read.table("./maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
pal <- rev(c("#ff40f2", "#ff0000", "#008c4b", "#00138c", "#8c235b", "#ffbfbf", "#8c7723", "#468c75", "#8091ff", "#ff80c4", "#8c3123", "#fff2bf", "#40fff2", "#69698c", "#ff0044", "#ff9180", "#e5ff80", "#bffbff", "#5940ff", "#8c696e", "#8c7369", "#858c69", "#40d9ff", "#c480ff", "#ff8c40", "#4b8c00", "#23698c", "#69238c", "#8c4b00", "#bfffbf", "#004b8c", "#eabfff", "#ffc480", "#40ff59", "#80c4ff", "#ffd940" ))
pal <- c("#e54545", "#735050", "#731406", "#ff4400", "#cc8970", "#f27724", "#331c03", "#ffe7cc", "#664605", "#e5b85c", "#8c8169", "#cad900", "#3c400a", "#668000", "#bce6a1", "#40f224", "#1d7334", "#40ffa6", "#1d402f", "#7ee6d1", "#7ca6a3", "#00add9", "#1b4d59", "#00294d", "#3d9df2", "#cce7ff", "#0c63e6", "#69778c", "#15358c", "#333640", "#00008c", "#070033", "#6026ff", "#695980", "#dcc2f2", "#b56cd9", "#3c0640", "#ff26d4", "#b33e7c", "#ff1a75", "#590929", "#e6a1b3", "#331418")





rownames(cazyme) <- gsub(".*;L4_",'',rownames(cazyme))
rownames(cazyme) <- gsub("L1_Others;Others;Others;",'',rownames(cazyme))
#order by mean value
cazyme <- cazyme[order(rowMeans(cazyme),decreasing = F),]


plot2 <- as.data.frame(t(cazyme))
plot2 <- rownames_to_column(plot2,var='X.SampleID')
plot2 <- melt(plot2,id = 'X.SampleID',variable.name = 'cazyme')
plot2 <- merge(plot2,map,by = 'X.SampleID')

plot2$present <- ifelse(plot2$value>0,1,0)

#change the order
plot2$UserName <- factor(plot2$UserName, levels = ord_factor)

# get right number of colors for plotting
no_cols <- length(unique(plot2$cazyme))
#colors_func <- sample(cols_func(no_cols))
colors_func <- sample(pal, no_cols,replace = T)
colors_func <- sample(colors(),no_cols)   


cazyme_plot <- ggplot(data = plot2, aes(x=StudyDayNo, y = cazyme)) +
  geom_point(aes(size=present,color=cazyme)) +
  facet_grid(.~UserName, scales = "free_y") +
  scale_color_manual(values = colors_func) +
  #scale_x_discrete(drop = FALSE) +
  theme_classic() + scale_x_continuous(breaks = seq(1,17,1))+
  theme(strip.text.x = element_text(angle = 0, size = 15, face = "italic"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_rect(color = "grey"), 
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        panel.spacing.x=unit(1.5, "lines")) +
  guides(color = guide_legend(reverse = TRUE, 
                             keywidth = 0.5, 
                             keyheight = 0.5, 
                             nrow = 5),size=FALSE)

cazyme_plot

ggsave("./result/cazy_pre.pdf",width=150,height=80,dpi = 600,limitsize = F)




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
colors_func <- sample(colors(),no_cols)  

#change the order
plot1$UserName <- factor(plot1$UserName, levels = ord_factor)


cazyme_plot <- ggplot(data = plot1, aes(x=StudyDayNo, y = value, fill=cazyme)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_func) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(strip.text.x = element_text(angle = 0, size = 9, face = "italic"),
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





