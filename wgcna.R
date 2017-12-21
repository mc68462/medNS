library(igraph)
library(plyr)
library(vegan)
library(gridExtra)
library(gplots)
library(compositions)
source("taxplot_functions.R")

# set working directory
workingDir = "."
setwd(workingDir)
dir.create("results")
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# load filtered phyloseq table from data.distribution.R saved in results
filt2011 <- t(as.matrix(as.data.frame(otu_table(readRDS("results/filtered.phyloseq.2011.rds")))))
filt2012 <- t(as.matrix(as.data.frame(otu_table(readRDS("results/filtered.phyloseq.2012.rds")))))

# pseudocount + clr transformation for downstream module eigenegene calculation
filt2011.p <- filt2011 + 1e-6
clr2011 <- compositions::clr(filt2011.p)

filt2012.p <- filt2012 + 1e-6
clr2012 <- compositions::clr(filt2012.p)

commonOTUs <- intersect(colnames(clr2011), colnames(clr2012))
clr2011.com <- clr2011[,commonOTUs]
clr2012.com <- clr2012[,commonOTUs]

# subset only interesected OTUs for untransformed, raw abundance matrices
OTU.2011.common <- filt2011[,commonOTUs]
OTU.2012.common <- filt2011[,commonOTUs]

# abundance measures of common OTUs
OTU.2011.common.ra <- vegan::decostand(OTU.2011.common[,commonOTUs], "total", MARGIN=1)
OTU.2012.common.ra <- vegan::decostand(OTU.2011.common[,commonOTUs], "total", MARGIN=1)
rowSums(OTU.2011.common.ra)
rowSums(OTU.2012.common.ra)

# calculate mean OTU abundance
OTU.2011.common.mean=colMeans(OTU.2011.common.ra)
OTU.2012.common.mean=colMeans(OTU.2012.common.ra)
commonOTUs.means=cbind(OTU.2011.common.mean, OTU.2012.common.mean)
mean=rowMeans(commonOTUs.means)

# calcuate median OTU obundance
OTU.2011.common.med=apply(OTU.2011.common.ra, 2, median)
OTU.2012.common.med=apply(OTU.2012.common.ra, 2, median)
commonOTUs.meds=cbind(OTU.2011.common.med, OTU.2012.common.med)
median=rowMeans(commonOTUs.meds)

# calculate geometric mean OTU abundance
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
OTU.2011.common.gmean=apply(OTU.2011.common.ra, 2, gm_mean)
OTU.2012.common.gmean=apply(OTU.2012.common.ra, 2, gm_mean)
commonOTUs.gmeans=cbind(OTU.2011.common.gmean, OTU.2012.common.gmean)
gmean=rowMeans(commonOTUs.gmeans)

# calculate max OTU abundance
OTU.2011.common.max=apply(OTU.2011.common.ra, 2, max)
OTU.2012.common.max=apply(OTU.2012.common.ra, 2, max)
commonOTUs.max=cbind(OTU.2011.common.max, OTU.2012.common.max)
max=rowMeans(commonOTUs.max)

# merge OTU abundance metrics into a single dataframe
commonOTUs.abund=data.frame(mean, gmean, max)
colnames(commonOTUs.abund) <- c("mean", "gmean", "max")

# Place data in to list of 2 datasets
nSets = 2;
setLabels = c("2011", "2012")
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(clr2011.com));
names(multiExpr[[1]]$data) = colnames(clr2011.com);
rownames(multiExpr[[1]]$data) = rownames(clr2011.com);
multiExpr[[2]] = list(data = as.data.frame(clr2012.com));
names(multiExpr[[2]]$data) = colnames(clr2012.com);
rownames(multiExpr[[2]]$data) = rownames(clr2012.com);


# Check that the data has the correct format for many functions operating on multiple sets
exprSize = checkSets(multiExpr)

# Create list of environmental trait data
md = read.table("data/TableS4_metadata.txt", header=T, sep='\t', as.is=T);

traitData = subset(md, md$year!=2010)
dim(traitData)
names(traitData)

# remove columns with unneeded information
# NOX and DIN show multicolinearity with vifcor() - exclude NOX
allTraits = traitData[, -c(2:7,15,16,17,20:23)];

# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)

# Form a multi-set structure that will hold the sample data
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$sample_name);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();

# Define data set dimensions
Dims = vector(mode ="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data); 
  Dims[[set]] =  list(data = length(rownames(multiExpr[[set]])))
}
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

# Import sparcc correlation matrix and pseudo pvalue matrix - use each value as a node attribute - for Antonio's script and for comparison of cor and pval from sparCC
cor2011 <- as.matrix(read.table("data/cor_sparcc.2011.out", sep='\t', header=T, check.names = F, row.names = 1))
cor2011.S <-  0.5*(cor2011 + 1)

cor2012 <- as.matrix(read.table("data/cor_sparcc.2012.out", sep='\t', header=T, check.names = F, row.names = 1))
cor2012.S <- 0.5*(cor2012 + 1)

# select power that achieves a scale-free topology R2 of > 0.8
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set - TRY SPEARMAN?
adjacencies[1,,] <- cor2011.S
adjacencies[2,,] <- cor2011.S

dimnames(adjacencies)[[1]] <- c("2011", "2012")
dimnames(adjacencies)[[2]] <- as.vector(colnames(cor2011))
dimnames(adjacencies)[[3]] <- as.vector(colnames(cor2011))

# Run TOM (Topological Overlap Matrix) analysis - initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ],TOMType="unsigned");

# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale 2012 TOM
  if (set>1)
  {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
}

# find consensus TOM using parallel minimum
consensusTOM = pmin(TOM[1, , ], TOM[2, , ]);
# Clustering - consider trying ward.D or ward.D2 for the method
consTree = hclust(as.dist(1-consensusTOM), method = "average");

# Run with module size that give maximum modularity
minModuleSize = 10;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );

# plot dendrogram of unmerged modules
c=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#e41a1c','#377eb8','#4daf4a','#984ea3','#ffff33','#a65628','#f781bf','#ffffb3','#fb8072')
c=c[c(1:length(unique(unmergedLabels)))]
unmergedColors = labels2colors(unmergedLabels, zeroIsGrey=T, colorSeq=c)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedLabels)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);

# Cluster consensus modules - CONSIDER DIFFERENT CLUSTERING METHOD
consMETree = hclust(as.dist(consMEDiss), method = "average")

# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")

# merge close modules
merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)

# Numeric module labels
moduleLabels = merge$colors;
# make vector of nicer colors
c=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#e41a1c','#377eb8','#4daf4a','#984ea3','#ffff33','#a65628','#f781bf','#ffffb3','#fb8072')
c=c[c(1:length(unique(unmergedLabels)))]
moduleColors=labels2colors(labels=moduleLabels, zeroIsGrey=FALSE, colorSeq=c)

# Eigengenes of the new merged modules:
consMEs = merge$newMEs;

# Perform redundancy analysis (RDA) with module eigengenes (MEs) and trait data (explanatory variables)
X.2011.o <- as.matrix(Traits[[1]]$data)
X.2011 <- X.2011.o[,c(1,2,3,4,5,6,7,8)]
X.2011.s <- vegan::decostand(X.2011, method = "standardize", na.rm = TRUE)

X.2012.o <- as.matrix(Traits[[2]]$data)
X.2012 <- X.2012.o[,c(1,2,3,4,5,6,7,8)]
X.2012.s <- vegan::decostand(X.2012, method = "standardize", na.rm = TRUE)

# extract consensus modeule eigenegenes from consMEs object
Y.2011 <- as.matrix(consMEs[[1]]$data)
rownames(Y.2011) <- gsub("FL_11", "", rownames(X.2011))
Y.2012 <- as.matrix(consMEs[[2]]$data)
rownames(Y.2012) <- gsub("FL_12", "",rownames(X.2012))

# prepare data for RDA - remove NAs
which(is.na(X.2011.s), arr.ind=TRUE)
XO.2011 <- X.2011.s[-c(45,35,65),]
YO.2011 <- Y.2011[-c(45,35,65),]

which(is.na(X.2012), arr.ind=TRUE)
XO.2012 <- X.2012[-c(17,27,32,37,41,44),]
YO.2012 <- Y.2012[-c(17,27,32,37,41,44),]

# RDA
par(mfrow=c(1,2))
rda.2011 <- vegan::rda(YO.2011, XO.2011)
plot(rda.2011)
rda.2012 <- vegan::rda(YO.2012, XO.2012)
plot(rda.2012)

par(mfrow=c(1,2))
hist(residuals(rda.2011))
hist(residuals(rda.2012))

# Assign seasonal traits to modules and set up metadata for plotting relative abundance barplots
# merged modules
df.trait.merged <- data.frame(module=as.character(c(1,2,3,4,5,7,8)), 
                       trait_cor=c("spring", "autumn", "summer", "spring.summer", "winter.spring", "winter", "summer.autumn"))

# unmerged modules
df.trait.unmerged <- data.frame(module=as.character(c(1,2,3,4,5,6,7,8,9,10,11)), 
                       trait_cor=c("spring.bloom", "autumn.winter", "summer.bloom", "summer.bloom", "spring.bloom", "autumn.winter", "autumn.winter", "summer.bloom", "spring.bloom", "autumn.winter", "autumn.winter"))

# Node metadata 
tax=as.data.frame(tax_table(cogito.physeq))
cols=c("domain", "phylum", "class", "order", "family", "genus")
tax$ltax=apply(tax[,cols], 1, paste , collapse = ";")
nodes=colnames(multiExpr[[1]]$data)

nd <- data.frame(node = colnames(multiExpr[[1]]$data),
                 module = as.character(moduleLabels),
                 domain = tax$domain[match(nodes, rownames(tax))],
                 phylum = tax$phylum[match(nodes, rownames(tax))],
                 class = tax$class[match(nodes, rownames(tax))],
                 order = tax$order[match(nodes, rownames(tax))],
                 family = tax$family[match(nodes, rownames(tax))],
                 genus = tax$genus[match(nodes, rownames(tax))],
                 ltax = tax$ltax[match(nodes, rownames(tax))],
                 mean.2011 = OTU.2011.common.mean,
                 mean.2012 = OTU.2012.common.mean,
                 median.2011 = OTU.2011.common.med,
                 median.2012 = OTU.2012.common.med,
                 gmean.2011 = OTU.2011.common.gmean,
                 gmean.2012 = OTU.2012.common.gmean,
                 max.2011 = OTU.2011.common.max,
                 max.2012 = OTU.2012.common.max)
megamodule <- df.trait.merged$trait_cor[match(nd$module, df.trait.merged$module)]
nodeData <- cbind(nd, megamodule)

nodeData <- nodeData[order(as.numeric(nodeData$module)),]

# module data
module <- sort(as.numeric(unique(nodeData$module)))
OF <- ddply(nodeData, .(module), summarize, freq=length(node))
OTU_freq <- OF[order(as.numeric(sort(OF$module))),]
PO <- ddply(nodeData, .(module), summarize, round((length(node)/nGenes)*100, digits = 1))
pct_OTUs <- PO[order(as.numeric(sort(PO$module))),]
mean.2011 <- ddply(nodeData, .(module), summarize, mean=sum(mean.2011))
mean_OTU_abund.2011 <- mean.2011[order(as.numeric(sort(mean.2011$module))),]
mean.2012 <- ddply(nodeData, .(module), summarize, mean=sum(mean.2012))
mean_OTU_abund.2012 <- mean.2012[order(as.numeric(sort(mean.2012$module))),]
mean_OTU_abund <- rowMeans(cbind(mean_OTU_abund.2011$mean, mean_OTU_abund.2012$mean))

gmean.2011 <- ddply(nodeData, .(module), summarize, gmean=sum(gmean.2011))
gmean_OTU_abund.2011 <- gmean.2011[order(as.numeric(sort(gmean.2011$module))),]
gmean.2012 <- ddply(nodeData, .(module), summarize, gmean=sum(gmean.2012))
gmean_OTU_abund.2012 <- gmean.2012[order(as.numeric(sort(gmean.2012$module))),]
gmean_OTU_abund <- rowMeans(cbind(gmean_OTU_abund.2011$gmean, gmean_OTU_abund.2012$gmean))

max.2011 <- ddply(nodeData, .(module), summarize, max=sum(max.2011))
max_OTU_abund.2011 <- max.2011[order(as.numeric(sort(max.2011$module))),]
max.2012 <- ddply(nodeData, .(module), summarize, max=sum(max.2012))
max_OTU_abund.2012 <- max.2012[order(as.numeric(sort(max.2012$module))),]
max_OTU_abund <- rowMeans(cbind(max_OTU_abund.2011$max, max_OTU_abund.2012$max))
TC <- as.data.frame(df.trait.merged)

tdf=data.frame(module=module, OTU_freq=OTU_freq$freq, 
               pct_OTUs=pct_OTUs[,2], 
               mean_OTU_abund=mean_OTU_abund, 
               max_OTU_abund=max_OTU_abund,
               gmean_OTU_abund=gmean_OTU_abund, 
               trait_cor=TC$trait_cor) #[-1,]
tdf$module <- factor(tdf$module, levels = as.numeric(tdf$module))
tdf$trait_cor <- factor(tdf$trait_cor, levels = c("spring", "autumn", "summer", "spring.summer", "winter.spring", "winter", "summer.autumn"))


# Plot OTUs within each consensus module
lmod=dlply(nodeData, c("module"), function(x) x=x)
lmod.df=data.frame(module=ldply(names(lmod), function(x) paste("ME",unique(lmod[[x]]$module), sep='')))
# make list of grepapable characters for each module
lgrep=lapply(lmod, function(x) gsub("(\\|$)", "", capture.output(cat(gsub("(.*)", "^\\1_|", (x$node)), sep=''))))
names(lgrep)<-lmod.df$V1

# merged modules 0.25
c1=c("ME5", "ME7", "ME2", "ME1", "ME8", "ME4", "ME3")

# SELECT one megamodule at a time for plotting - change 'm' and 'mod'- run until dev.off()
m="c1"
mod=c1

pdf(paste("results/", m, ".otu.abund.2011-12.n10.pdf", sep=''), height=8, width=10*length(mod))

p2011=list()
p2012=list()

plotyear="2011"
p2011=lapply(
  mod,
  function(i) taxplot_subnets(otu_matrix="data/TableS4.txt", 
                              sample_data="data/TableS4_metadata.txt", 
                              taxa=lgrep[[i]], n=10, abund=0.0, main=names(lgrep)[which(names(lgrep)==i)], Y=as.numeric(plotyear), legend="right")
)

plotyear="2012"
p2012=lapply(
  mod,
  function(i) taxplot_subnets(otu_matrix="data/TableS4.txt", 
                              sample_data="data/TableS4_metadata.txt", 
                              taxa=lgrep[[i]], n=10, abund=0.0, main=names(lgrep)[which(names(lgrep)==i)], Y=as.numeric(plotyear), legend="right")
)
do.call(grid.arrange, c(p2011, p2012, nrow=2))
dev.off()


# Convert Sparcc graphs into igraph objects with annotations
cor2011 <- cor2011
cor2011m <- melt(cor2011)
pval2011 <- as.matrix(read.table("data/pvals.2011.two-sided.txt", sep='\t', header=T, check.names = F, row.names = 1))
pval2011m <- melt(pval2011)

g1.df=cbind(cor2011m, pval2011m)[,c(1:3,6)]
colnames(g1.df) <- c("X1", "X2", "cor", "pval")
g1 <- igraph::simplify(graph_from_data_frame(g1.df, directed = F), remove.multiple = T, remove.loops = T, edge.attr.comb = list(cor=function(x) sum(x)/2, pval=function(x) sum(x)/2))

cor2012 <- cor2012
cor2012m <- melt(cor2012)
pval2012 <- as.matrix(read.table("data/pvals.2011.two-sided.txt", sep='\t', header=T, check.names = F, row.names = 1))
pval2012m <- melt(pval2012)

g2.df=cbind(cor2012m, pval2012m)[,c(1:3,6)]
colnames(g2.df) <- c("X1", "X2", "cor", "pval")
g2 <- igraph::simplify(graph_from_data_frame(g2.df, directed = F), remove.multiple = T, remove.loops = T, edge.attr.comb = list(cor=function(x) sum(x)/2, pval=function(x) sum(x)/2))

# 2011
V(g1)$module <- nodeData$module[match(V(g1)$name, nodeData$node)]
V(g1)$megamodule <- nodeData$megamodule[match(V(g1)$name, nodeData$node)]
V(g1)$domain <- as.character(nodeData$domain[match(V(g1)$name, nodeData$node)])
V(g1)$phylum <- as.character(nodeData$phylum[match(V(g1)$name, nodeData$node)])
V(g1)$class <- as.character(nodeData$class[match(V(g1)$name, nodeData$node)])
V(g1)$order <- as.character(nodeData$order[match(V(g1)$name, nodeData$node)])
V(g1)$family <- as.character(nodeData$family[match(V(g1)$name, nodeData$node)])
V(g1)$genus <- as.character(nodeData$genus[match(V(g1)$name, nodeData$node)])
V(g1)$mean <- as.numeric(nodeData$mean.2011[match(V(g1)$name, nodeData$node)])
V(g1)$max <- as.numeric(nodeData$max.2011[match(V(g1)$name, nodeData$node)])
V(g1)$gmean <- as.numeric(nodeData$gmean.2011[match(V(g1)$name, nodeData$node)])
V(g1)$color <- "#878787"
V(g1)$color[which(V(g1)$class == "Flavobacteriia")] <- "#d53e4f"
V(g1)$color[which(V(g1)$class == "Gammaproteobacteria")] <- "#5e4fa2"
V(g1)$color[which(V(g1)$class == "Alphaproteobacteria")] <- "#3288bd"
V(g1)$color[which(V(g1)$class == "Betaproteobacteria")] <- "#e6f598"
V(g1)$color[which(V(g1)$class == "Acidimicrobiia")] <- "#fdae61"
V(g1)$color[which(V(g1)$class == "Thermoplasmata")] <- "#66c2a5"
V(g1)$color[which(V(g1)$class == "Cytophagia")] <- "#f46d43"
V(g1)$color[which(V(g1)$class == "Planctomycetacia")] <- "#abdda4"
V(g1)$color[which(V(g1)$class == "Epsilonproteobacteria")] <- "#fee08b"
E(g1)$weight <- E(g1)$cor
E(g1)$color[which(E(g1)$weight > 0)] <- "#bdbdbd"
E(g1)$color[which(E(g1)$weight < 0)] <- "#d53e4f"

# 2012
V(g2)$module <- nodeData$module[match(V(g2)$name, nodeData$node)]
V(g2)$megamodule <- nodeData$megamodule[match(V(g2)$name, nodeData$node)]
V(g2)$domain <- as.character(nodeData$domain[match(V(g2)$name, nodeData$node)])
V(g2)$phylum <- as.character(nodeData$phylum[match(V(g2)$name, nodeData$node)])
V(g2)$class <- as.character(nodeData$class[match(V(g2)$name, nodeData$node)])
V(g2)$order <- as.character(nodeData$order[match(V(g2)$name, nodeData$node)])
V(g2)$family <- as.character(nodeData$family[match(V(g2)$name, nodeData$node)])
V(g2)$genus <- as.character(nodeData$genus[match(V(g2)$name, nodeData$node)])
V(g2)$mean <- as.numeric(nodeData$mean.2012[match(V(g2)$name, nodeData$node)])
V(g2)$max <- as.numeric(nodeData$max.2012[match(V(g2)$name, nodeData$node)])
V(g2)$gmean <- as.numeric(nodeData$gmean.2012[match(V(g2)$name, nodeData$node)])
V(g2)$color <- "#878787"
V(g2)$color[which(V(g2)$class == "Flavobacteriia")] <- "#d53e4f"
V(g2)$color[which(V(g2)$class == "Gammaproteobacteria")] <- "#5e4fa2"
V(g2)$color[which(V(g2)$class == "Alphaproteobacteria")] <- "#3288bd"
V(g2)$color[which(V(g2)$class == "Betaproteobacteria")] <- "#e6f598"
V(g2)$color[which(V(g2)$class == "Acidimicrobiia")] <- "#fdae61"
V(g2)$color[which(V(g2)$class == "Thermoplasmata")] <- "#66c2a5"
V(g2)$color[which(V(g2)$class == "Cytophagia")] <- "#f46d43"
V(g2)$color[which(V(g2)$class == "Planctomycetacia")] <- "#abdda4"
V(g2)$color[which(V(g2)$class == "Epsilonproteobacteria")] <- "#fee08b"
E(g2)$weight <- E(g2)$cor
E(g2)$color[which(E(g2)$weight > 0)] <- "#bdbdbd"
E(g2)$color[which(E(g2)$weight < 0)] <- "#d53e4f"

module_name_list_2011=unique(V(g1)$module)
module_name_list_2012=unique(V(g2)$module)

# make list of individual module igraphs 
module_graphs_2011=list()
module_dfs_2011=list()
for (m in module_name_list_2011){
  gv=delete_vertices(g1, which(V(g1)$module != m))
  gv$name <- paste("mod.", m, sep='')
  gvf <- delete.vertices(gv,  which(V(gv)$max < 0))
  gvef <- delete.edges(gvf, which(abs(E(gvf)$cor) < 0))
  gvefp <- delete.edges(gvef, which(E(gvef)$pval > 1))
  V(gvefp)$size <- (V(gvefp)$max)*200
  V(gvefp)$label.cex <- 0.8
  E(gvefp)$weight <- E(gvefp)$cor
  module_graphs_2011[[m]] = gvefp
  module_dfs_2011[[m]] = get.data.frame(gvefp, what="edges")
}

module_graphs_2012=list()
module_dfs_2012=list()
for (m in module_name_list_2012){
  gv=delete_vertices(g2, which(V(g2)$module != m))
  gv$name <- paste("mod.", m, sep='')
  gvf <- delete.vertices(gv,  which(V(gv)$max < 0.0))
  gvef <- delete.edges(gvf, which(abs(E(gvf)$cor) < 0))
  gvefp <- delete.edges(gvef, which(E(gvef)$pval > 1))
  V(gvefp)$size <- (V(gvefp)$max)*200
  V(gvefp)$label.cex <- 0.8
  E(gvefp)$weight <- E(gvefp)$cor
  module_graphs_2012[[m]] = gvefp
  module_dfs_2012[[m]] = get.data.frame(gvefp, what="edges")
}

# List of individual unfiltered module network igraph objects
saveRDS(module_graphs_2011, file = "results/unfiltered_2011_module_igraph_list.rds")
saveRDS(module_graphs_2012, file = "results/unfiltered_2012_module_igraph_list.rds")
# List of individual unfiltered module network edge data frames 
saveRDS(module_dfs_2012, file = "results/unfiltered_2011_module_edgelist.rds")
saveRDS(module_dfs_2011, file = "results/unfiltered_2012_module_edgelist.rds")

# generate filtered module networks based on analysis of network graph deconstruction (Figure S12)
p2011g <- list()
p2011f <- list()
for (m1 in 1:length(module_name_list_2011)){
  test_graph <- module_graphs_2011[[m1]]
  m.name1 <- test_graph$name
  filt_list <- c(0.51, 0.53, 0.51, 0.61, 0.62, 0.52, 0.5)  # based on analysis of network graph deconstruction
  filt_graph_e <- igraph::delete_edges(test_graph, which(abs(E(test_graph)$weight) < filt_list[m1]))
  filt_graph_v <- igraph::delete_vertices(filt_graph_e, V(filt_graph_e)[degree(filt_graph_e) == 0])
  V(filt_graph_v)$degree <- degree(filt_graph_v, mode='all', normalized = F)
  V(filt_graph_v)$strength <- graph.strength(filt_graph_v, mode='all', weights=E(filt_graph_v)$cor)
  V(filt_graph_v)$betweenness <- betweenness(filt_graph_v, directed = F, weights = abs(E(filt_graph_v)$cor), normalized = T)
  E(filt_graph_v)$color[which(E(filt_graph_v)$weight > 0)] <- "#bdbdbd"
  E(filt_graph_v)$color[which(E(filt_graph_v)$weight < 0)] <- "#d53e4f"
  E(filt_graph_v)$dir[which(E(filt_graph_v)$weight > 0)] <- "positive"
  E(filt_graph_v)$dir[which(E(filt_graph_v)$weight < 0)] <- "negative"
  p2011f[[m1]] <- filt_graph_v
  write.graph(filt_graph_v, format = "gml", file = paste("results/", m.name1, ".2011.filt.gml", sep = ''))
}

p2012g <- list()
p2012f <- list()
for (m2 in 1:length(module_name_list_2012)){
  test_graph <- module_graphs_2012[[m2]]
  m.name2 <- test_graph$name
  filt_list <- c(0.58, 0.55, 0.56, 0.61, 0.57, 0.59, 0.55) # based on analysis of network graph deconstruction
  filt_graph_e <- igraph::delete_edges(test_graph, which(abs(E(test_graph)$weight) < filt_list[m2]))
  filt_graph_v <- igraph::delete_vertices(filt_graph_e, V(filt_graph_e)[degree(filt_graph_e) == 0])
  V(filt_graph_v)$degree <- degree(filt_graph_v, mode='all', normalized = F)
  V(filt_graph_v)$strength <- graph.strength(filt_graph_v, mode='all', weights=E(filt_graph_v)$cor)
  V(filt_graph_v)$betweenness <- betweenness(filt_graph_v, directed = F, weights = abs(E(filt_graph_v)$cor), normalized = T)
  E(filt_graph_v)$color[which(E(filt_graph_v)$weight > 0)] <- "#bdbdbd"
  E(filt_graph_v)$color[which(E(filt_graph_v)$weight < 0)] <- "#d53e4f"
  E(filt_graph_v)$dir[which(E(filt_graph_v)$weight > 0)] <- "positive"
  E(filt_graph_v)$dir[which(E(filt_graph_v)$weight < 0)] <- "negative"
  p2012f[[m2]] <- filt_graph_v
  write.graph(filt_graph_v, format = "gml", file = paste("results/", m.name2, ".2012.filt.gml", sep = ''))
}

saveRDS(p2011f, file = "results/filtered_2011_module_igraph_objects.rds")
saveRDS(p2012f, file = "results/filtered_2012_module_igraph_objects.rds")

# intersect module networks to evaluate network conservation across 2011 and 2012
g.inters <- list()

for (i in 1:length(module_name_list_2011)){
  net1 <- p2011f[[i]]
  net2 <- p2012f[[i]]
  
  g <- graph.intersection(net1, net2, keep.all.vertices = FALSE)
  E(g)$weight <- (E(g)$cor_1 + E(g)$cor_2)/2
  E(g)$pval <- (E(g)$pval_1 + E(g)$pval_2)/2
  V(g)$max <- (V(g)$max_1 + V(g)$max_2)/2
  V(g)$color <- V(g)$color_1
  E(g)$color <- E(g)$color_1
  V(g)$degree <- degree(g, mode='all', normalized = F)
  V(g)$strength <- strength(g, mode='all')
  V(g)$betweenness <- betweenness(g, directed = F, weights = abs(E(g)$weight), normalized = T)
  g$name <- g$name_1
  gd <- igraph::delete_vertices(g, V(g)[degree(g) == 0])
  g.inters[[i]] <- gd
}

saveRDS(g.inters, file = "results/intersected_module_igraph_objects.rds")
