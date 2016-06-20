library(dplyr)
library(wgcna)
library(ggplot)
library(cowplot)

# test different minModuleSize - choose one with maximum modularity
sizes=seq(5,100, by=1)
eig=list()
eig2011=list()
eig2012=list()

for (s in sizes){
# Run with module size that give maximum modularity
minModuleSize = s;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );

# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedLabels)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);

# Cluster consensus modules - CONSIDER DIFFERENT CLUSTERING METHOD
consMETree = hclust(as.dist(consMEDiss), method = "average");

# merge close modules
merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)

# Eigengenes of the new merged modules:
consMEs = merge$newMEs;

# Examine variance explained by eigenvector (PC1)

eig2011[[s]]=moduleEigengenes(multiExpr[[1]]$data, 
                         merge$colors, 
                         impute = TRUE, 
                         nPC = 1, 
                         align = "along average", 
                         excludeGrey = TRUE, 
                         grey = if (is.numeric(merge$colors)) 0 else "grey",
                         subHubs = TRUE,
                         trapErrors = FALSE, 
                         softPower = 10,
                         scale = TRUE,
                         verbose = 0, indent = 0)

eig2012[[s]]=moduleEigengenes(multiExpr[[2]]$data, 
                              merge$colors, 
                              impute = TRUE, 
                              nPC = 1, 
                              align = "along average", 
                              excludeGrey = TRUE, 
                              grey = if (is.numeric(merge$colors)) 0 else "grey",
                              subHubs = TRUE,
                              trapErrors = FALSE, 
                              softPower = 10,
                              scale = TRUE,
                              verbose = 0, indent = 0)

eig[[s]]=as.data.frame(t(rbind(round(eig2011[[s]]$varExplained, digits=3), round(eig2012[[s]]$varExplained, digits=3), rep(s, times=length(eig2011[[s]]$varExplained)))))
}

eig=eig[!sapply(eig, is.null)]

save(eig, file="optimize.min.mod.size.eigen.var.RData")

eig.df=ldply(eig)
colnames(eig.df)<-c("y2011", "y2012", "min.module.size")
eig.mdf=melt(eig.df, id.var="min.module.size")
p1=ggplot(eig.mdf, aes(x=factor(min.module.size), y=value, fill=factor(variable))) + 
  geom_boxplot() +
  xlab("Minimum module size") +
  ylab("% variance explained (PC1)") +
  guides(fill=guide_legend(title=NULL)) +
  scale_x_discrete(breaks = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100))
  
eig.mean=ddply(eig.mdf, "min.module.size", summarize, mean(value))
colnames(eig.mean)<-c("min.module.size", "mean.var")
p2=ggplot(eig.mean, aes(x=min.module.size, y=mean.var)) + geom_point(aes(x=min.module.size, y=mean.var)) +
  geom_vline(xintercept = c(18), linetype="solid", color = "grey", size=0.5) +
  xlab("Minimum module size") +
  ylab("Mean % variance explained (PC1)") +
  scale_x_continuous(breaks = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100))

p=plot_grid(p1, p2, labels=c("A", "B"), nrow=2, ncol=1)
save_plot("opt.module.size.eigen.var.expl.pdf", p, ncol = 1, nrow = 2, base_height = 4,
          base_width = 12)


