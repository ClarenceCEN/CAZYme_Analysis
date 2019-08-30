require(dplyr)
require(ggplot2)
require(ggsci)

load('./data/cazy_list_clr_2.RData')

cazy_dist <- lapply(cazyme_list_clr, function(x){
  x <- select(x,-X.SampleID)
  
  d <- dist(x)
  d <- as.matrix(d)
  d[upper.tri(d)] = 0
  d <- as.vector(d[d!=0])
  d
})

#cazy_dist <- do.call('cbind',cazy_dist)

soylent.dist <- c(cazy_dist[['MCTs11']],cazy_dist[['MCTs12']])

soylent.dist <- data.frame(Group='soylent',distance=soylent.dist)

regular.dist <- c()
for(i in names(cazy_dist)){
  if(!i%in%c("MCTs11", "MCTs12")){
    regular.dist <- c(regular.dist,cazy_dist[[i]])
  }
}

regular.dist <- data.frame(Group='regular',distance=regular.dist)

dist.data <- rbind(regular.dist,soylent.dist)

g <- ggplot(dist.data,aes(x=Group,y=distance,fill=Group)) + geom_boxplot() +
  theme_classic() + scale_fill_ucscgb() + ylab("Aichison's distance") + 
  theme(axis.title.x = element_blank()) + guides(fill=FALSE)
g
