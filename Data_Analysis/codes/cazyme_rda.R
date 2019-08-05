
load('./data/cazyme_median_clr.RData')
load("./data/cazymes_to_keep_0801.RData")

cazyme_mean_clr <- cazyme_mean_clr[,colnames(cazyme_mean_clr)%in%cazyme_to_keep]
map <- read.table("./maps/food_map.txt", sep = "\t", header = T, comment = "")

cazyme_mean <- cazyme_mean[rownames(cazyme_mean)%in%rownames(food_L3),]

map_t <- aggregate(map[c('Gender','Age','BMI')],list(map$UserName),unique)
rownames(map_t) <- map_t$Group.1;map_t <- map_t[-1]

food <- read.delim("./data/dhydrt.smry.no.soy.txt",row.names = 1)
food$taxonomy <- as.character(food$taxonomy)
split <- strsplit(food$taxonomy,";") 
foodStrings <- sapply(split,function(x) paste(x[1:1],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",foodStrings,perl=T) # clean tips
food_L1 <- rowsum(food[-ncol(food)],foodStrings) 
rownames(food_L1) <- gsub("L1_",'',rownames(food_L1))
food_L1 <- as.data.frame(t(food_L1))

#food <- food[,colnames(food)%in%rownames(map_t)]




source(file = "G:/Dan_Lab/dietstudy_analyses-master/lib/colors/UserNameColors.R")
UserNameColors['MCTs05'] <- '#99ff00'

map_t <- map_t[rownames(map_t)%in%rownames(food),]
cazyme_mean_clr <- cazyme_mean_clr[rownames(cazyme_mean_clr)%in%rownames(food),]


#pca
library(vegan);library(ape)
cazyme_dist <- dist(cazyme_mean_clr)

cazyme_pcoa <- pcoa(cazyme_dist)
cazyme_plots <- as.data.frame(cazyme_pcoa$vectors[,1:2])
pc1 <- cazyme_pcoa$values[1,2];pc2 <- cazyme_pcoa$values[2,2]

#rda
cazyme_rda <- rda(cazyme_mean_clr~.,map_t)
cazyme_rda_s <- summary(cazyme_rda)
anova.cca(cazyme_rda)
rda1 <- cazyme_rda_s[["cont"]][["importance"]][2,1]
rda2 <- cazyme_rda_s[["cont"]][["importance"]][2,2]
basic_explanation <- (rda1+rda2)/(pc1+pc2)   #34% 

cazyme_food_rda <- rda(cazyme_mean_clr~.,food_L1)
#a <- step(cazyme_food_rda)
cazyme_food_rda_s <- summary(cazyme_food_rda)
anova.cca(cazyme_food_rda)
plot(cazyme_food_rda)
rda1 <- cazyme_food_rda_s[["cont"]][["importance"]][2,1]
rda2 <- cazyme_food_rda_s[["cont"]][["importance"]][2,2]

food_explanation <- (rda1+rda2)/(pc1+pc2)   #55%, higher than food-microbiome

cazyme.food <- as.data.frame(cazyme_food_rda_s$biplot)
cazyme.subject <- as.data.frame(cazyme_food_rda_s$sites)



library(ggplot2)
g <- ggplot(cazyme_plots,aes(x=Axis.1,y=Axis.2))+geom_point(aes(color=rownames(cazyme_plots)),size=3)+theme_classic()+
  xlab(paste("PC1",round(pc1*100,2),"%")) +ylab(paste("PC2",round(pc2*100,2),"%")) +
  scale_color_manual(values = UserNameColors) + labs(color='Subject')
g
ggsave("./result/pcoa.pdf",width = 7,height = 5)

library(ggplot2)
g <- ggplot(cazyme.subject,aes(x=RDA1,y=RDA2))+geom_point(aes(color=rownames(cazyme.subject)),size=3)+theme_classic()+
  xlab(paste("RDA1",round(rda1*100,2),"%")) +ylab(paste("RDA2",round(rda2*100,2),"%")) +
  scale_color_manual(values = UserNameColors) + labs(color='Subject') 
for (i in 1:nrow(cazyme.food)){
  g <- g  + annotate("text",x = cazyme.food[i,1]*7+0.4, y = cazyme.food[i,2]*7+0.4, label = rownames(cazyme.food)[i], color = "darkorange", size = 4)
  g <- g  + annotate("segment", x = 0, y = 0, xend = cazyme.food[i,1]*7, yend = cazyme.food[i,2]*7, color = "darkorange", arrow = arrow(),size=1)
}
g
ggsave("./result/rda.pdf",width = 7,height = 5)
