#############################
#####UMAP and clustering#####
#############################

library(umap)
library(factoextra)
library(dplyr)
library(ggrepel)

met <- read.csv("/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/met_asy_a.csv", header = TRUE, row.names = 'ID')
attrac <- read.csv("/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/atractores_asy_a.csv", header = TRUE, row.names = 'ID')
attrac <- select(attrac, -attractor)


##########
###UMAP###
##########

attrac <- as.matrix(attrac)
umap_vis <- umap(attrac, controlscale = TRUE, scale = 3, dotsize = 1, random_state = 123, negative_sample_rate = 25L, n_epochs = 300) 
plot(umap_vis$layout)


############
#Clustering#
############

library(NbClust)

#This package provides 30 indices for determining the number of clusters.
Res <- NbClust(data = umap_vis$layout, distance = "euclidean",diss = NULL, 
               min.nc = 2, max.nc = 5, 
               method = "kmeans", index = "all")
Res$Best.nc

#############################
#To visualize the clustering#
#############################

umap_plot <- as.data.frame(umap_vis$layout)
umap_plot <- cbind(umap_plot,Res$Best.partition)
colnames(umap_plot) <- c("umap1","umap2","Cluster")
ggplot(umap_plot, aes(x = umap1,y = umap2))+
  geom_point( aes(color = factor(Cluster)))+
  scale_colour_manual("Clusters", values = c("red", "blue", "green"))+
  theme_bw()

#To get attractor tags 
write.csv(umap_plot, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/atractores_asy_t.csv", row.names = TRUE)





###################################
#####Corroborating projections#####
###################################

plot_cluster <- function(layout, var_cluster, palette, etiquetas)
{
  ggplot(layout, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=2) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) +
    geom_text_repel(aes(label=etiquetas$attractor),hjust=0, vjust=0)
}


###Hierarchical
umap_plot <- as.data.frame(umap_vis$layout)
#plot(umap_plot)
umap_original <- umap_plot
fit_cluster_hierarchical <- hclust(dist(scale(umap_vis$layout)))
plot(fit_cluster_hierarchical) 
k = max(Res$Best.partition)
umap_original$clusters <- factor(cutree(fit_cluster_hierarchical, k = k)) 
plot_th <- plot_cluster(umap_original, "clusters", "Paired", met)  
plot_th


###Kmeans
set.seed(123)
umap_plot <- as.data.frame(umap_vis$layout)
umap_original <- umap_plot
k = max(Res$Best.partition)
fit_cluster_kmeans <- kmeans(umap_plot, k)
umap_original$clusters_kmeans <- factor(fit_cluster_kmeans$cluster)
plot <- plot_cluster(umap_original, "clusters_kmeans", "Paired", met)
plot


sessionInfo()