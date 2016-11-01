library(igraph)
library(ggplot2)
library(plyr)
library(reshape)
library(splitstackshape)
library(phyloseq)
library(vegan)
library(cowplot)
library(gridExtra)
library(gplots)
library(WGCNA)

# set working directory
workingDir <- ""
setwd(workingDir); 

# set data directories to variables
meta <- "data/TableS4.sample.metadata.txt"

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# take filtered phyloseq objects cogito_physeq_2011_otu_sample_physeq_filt_prev and cogito_physeq_2012_otu_sample_physeq_filt_prev from OTU.pre_processing.R
cogito.physeq.2011.filt <- cogito_physeq_2011_otu_sample_physeq_filt_prev
cogito.physeq.2012.filt <- cogito_physeq_2012_otu_sample_physeq_filt_prev
save(cogito.physeq.2011.filt, cogito.physeq.2012.filt, file="filtered.cogito.data.Rata")

# extract OTU tables from filtered phyloseq objects
filt2011 <- t(as.matrix(as.data.frame(otu_table(cogito.physeq.2011.filt))))
filt2012 <- t(as.matrix(as.data.frame(otu_table(cogito.physeq.2012.filt))))

# Perform transformation: Hellinger transformation of log scaled abundances + pseudocount
hellog2011 <- decostand(log(filt2011 + 1 ), method="hellinger", MARGIN=1)
hellog2012 <- decostand(log(filt2012 + 1 ), method="hellinger", MARGIN=1)

# set OTU table to untransformed relative abundance matrix for calculations of OTU abundances for downstream plotting
OTU.2011 <- t(as.data.frame(as.matrix(otu_table(cogito.physeq.2011.filt))))
OTU.2012 <- t(as.data.frame(as.matrix(otu_table(cogito.physeq.2012.filt))))

# intersect 2011 and 2012 OTUs to include only those OTUs common in both datasets and then take the proportion
commonOTUs <- intersect (colnames(OTU.2011), colnames(OTU.2012))
commonOTUs2011.ra <- decostand(OTU.2011[,commonOTUs], "total", MARGIN=1)
commonOTUs2012.ra <- decostand(OTU.2012[,commonOTUs], "total", MARGIN=1)
# rowSums should == 1
rowSums(commonOTUs2011.ra)
rowSums(commonOTUs2012.ra)

# abundance measures of common OTUs for downstream plotting
# mean
commonOTUs2011.mean <- colMeans(commonOTUs2011.ra)
commonOTUs2012.mean <- colMeans(commonOTUs2012.ra)
commonOTUs.means <- cbind(commonOTUs2011.mean, commonOTUs2012.mean)
mean <- rowMeans(commonOTUs.means)

# median
commonOTUs2011.med <- apply(commonOTUs2011.ra, 2, median)
commonOTUs2012.med <- apply(commonOTUs2012.ra, 2, median)
commonOTUs.meds <- cbind(commonOTUs2011.med, commonOTUs2012.med)
median <- rowMeans(commonOTUs.meds)

# geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
commonOTUs2011.gmean <- apply(commonOTUs2011.ra, 2, gm_mean)
commonOTUs2012.gmean <- apply(commonOTUs2012.ra, 2, gm_mean)
commonOTUs.gmeans <- cbind(commonOTUs2011.gmean, commonOTUs2012.gmean)
gmean <- rowMeans(commonOTUs.gmeans)

# max
commonOTUs2011.max <- apply(commonOTUs2011.ra, 2, max)
commonOTUs2012.max <- apply(commonOTUs2012.ra, 2, max)
commonOTUs.max <- cbind(commonOTUs2011.max, commonOTUs2012.max)
max <- rowMeans(commonOTUs.max)

# combine abundance measures into a single dataframe
commonOTUs.abund <- data.frame(mean, gmean, max)
colnames(commonOTUs.abund) <- c("mean", "gmean", "max")

# set OTU transformed matrices to new variable
OTU.2011.trans <- hellog2011
OTU.2012.trans <- hellog2012

# intersect 2011 and 2012 transformed OTUs to include only those OTUs common in both datasets
commonOTUs <- intersect(colnames(OTU.2011.trans), colnames(OTU.2012.trans))
commonOTUs2011.trans <- OTU.2011.trans[,commonOTUs]
commonOTUs2012.trans <- OTU.2012.trans[,commonOTUs]

# place data in to list of 2 datasets for WGCNA
nSets = 2;
setLabels <- c("2011", "2012")
multiExpr <- vector(mode = "list", length = nSets)
multiExpr[[1]] <- list(data = as.data.frame(commonOTUs2011.trans));
names(multiExpr[[1]]$data) <- colnames(commonOTUs2011.trans);
rownames(multiExpr[[1]]$data) <- rownames(commonOTUs2011.trans);
multiExpr[[2]] <- list(data = as.data.frame(commonOTUs2012.trans));
names(multiExpr[[2]]$data) <- colnames(commonOTUs2012.trans);
rownames(multiExpr[[2]]$data) <- rownames(commonOTUs2012.trans);

# check that the data has the correct format for many functions operating on multiple sets
exprSize <- checkSets(multiExpr)

# load sample metadata
md <- read.table(meta, header=T, sep='\t', as.is=T);
# subset sample metadata to include only 2011 and 2012 to be used in WGCNA
traitData <- subset(md, md$year!=2010)
dim(traitData)
names(traitData)

# remove columns that hold information we do not need (i.e., non-continuous data)
# NOX and DIN show multicolinearity with vifcor() - exclude NOX 
allTraits = traitData[, c("sample_name", "Secchi.meter", "Temperature", "Salinity", "SiO2", "PO4", "NO2", "NO3", "BBE.Chla", "DIN")];
dim(allTraits)
names(allTraits)

# create a list containings subsampled sample metadata
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$sample_name);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();

# define data set dimensions
Dims = vector(mode ="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data); 
  Dims[[set]] =  list(data = length(rownames(multiExpr[[set]])))
}
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

# Save input objects
save(multiExpr, Traits, nGenes, nSamples, nSets, setLabels, exprSize, commonOTUs.abund, commonOTUs, commonOTUs2011.trans, commonOTUs2012.trans,
     file = "Consensus-dataInput.2011.2012.hellog.RData")

# Choose a set of soft-thresholding powers
powers = c(seq(4,20,by=1));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, 
                                                     RsquaredCut = 0.8,
                                                     networkType = 'signed',
                                                     corFnc = cor, 
                                                     corOptions = list(use = 'p'), 
                                                     powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();

# Plot mean connectivity and power
power.df=data.frame(powerTables[[1]]$data[,1],powerTables[[1]]$data[,2], powerTables[[2]]$data[,2])
colnames(power.df)<- c("power", "y2011", "y2012")
mean.df=data.frame(powerTables[[1]]$data[,1],powerTables[[1]]$data[,5], powerTables[[2]]$data[,5])
colnames(mean.df)<- c("power", "y2011", "y2012")

mpower=melt(power.df, id.vars="power")
sft=ggplot(mpower, aes(x = power, y = value, group=variable, colour=variable)) + 
  geom_line(size=0.5) + 
  geom_point(size=1) +
  ylab("Scale Free Topology R2") +
  xlab("Power") +
  geom_hline(yintercept=0.8, linetype='solid', color='#969696', size=0.35) +
  geom_vline(xintercept=15, linetype='solid', color='#969696', size=0.35) +
  labs(color = "Year")
mmean=melt(mean.df, id.vars="power")
mean=ggplot(mmean, aes(x = power, y = value, group=variable, colour=variable)) + 
  geom_line(size=0.5) + 
  geom_point(size=1) +
  ylab("Mean connectivity") +
  xlab("Power") +
  labs(color = "Year")

p=plot_grid(sft, mean, labels=c("A", "B"), ncol = 2, nrow = 1)
save_plot("scale.free.topology.power.table.pdf", p, ncol = 1, nrow = 1, base_height = 4,
          base_aspect_ratio = 2.4)

# run auto pipeline for detection of consensus modules 
# use a softPower that achieves a scale free-topology from above plot (i.e., R2 > 0.8) and minimum module size 10
softPower=15
minModule=10

net = blockwiseConsensusModules(
  multiExpr, 
  maxBlockSize=5000,
  power = softPower, 
  minModuleSize = minModule,
  randomSeed = 12345,
  TOMType = "signed",
  TOMDenom="min",
  networkType = "signed", 
  corType = "pearson",
  saveIndividualTOMs = TRUE,
  individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
  saveConsensusTOMs = TRUE,
  consensusTOMFileNames = "consensusTOM-block.%b.RData",
  consensusQuantile = 0,
  deepSplit = 2, 
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)

# numeric module labels
moduleLabels = net$colors;
moduleColors=labels2colors(labels=moduleLabels, zeroIsGrey=TRUE)

# set eigengenes of the consensus modules to a variable
consMEs = net$multiMEs

save(consMEs, moduleColors, merge, moduleLabels, file=paste("soft_power.",softPower,"min_module_size.", minModule,"EigenCons.RData", sep=''))

# examine variance explained by eigenvector (PC1)
eig2011=moduleEigengenes(multiExpr[[1]]$data, 
                         net$colors, 
                         impute = TRUE, 
                         nPC = 1, 
                         align = "along average", 
                         excludeGrey = FALSE, 
                         grey = if (is.numeric(net$colors)) 0 else "grey",
                         subHubs = TRUE,
                         trapErrors = FALSE, 
                         softPower = softPower,
                         scale = TRUE,
                         verbose = 0, indent = 0)

eig2012=moduleEigengenes(multiExpr[[2]]$data, 
                         net$colors,
                         impute = TRUE, 
                         nPC = 1, 
                         align = "along average", 
                         excludeGrey = FALSE, 
                         grey = if (is.numeric(net$colors)) 0 else "grey",
                         subHubs = TRUE,
                         trapErrors = FALSE, 
                         softPower = softPower,
                         scale = TRUE,
                         verbose = 0, indent = 0)

eig=as.data.frame(t(rbind(eig2011$varExplained, eig2012$varExplained)))
module<-as.vector(substring(names(eig2011$eigengenes), 3))
eig.df=data.frame(module,eig)
colnames(eig.df)<-c("module", "y2011", "y2012")
eig.df$module <- factor(eig.df$module, levels = eig.df$module)

meig=melt(eig.df, id.vars="module")
pvar=ggplot(data=meig, aes(x=module, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +
  scale_fill_manual(values=c("deeppink3", "#1f78b4")) +
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.key = element_blank()) +
  xlab("Module") +
  scale_x_discrete(breaks=module, expand=c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0.005, 0.005)) +
  ylab("PC1 % variance explained")

pdf("var.expl.pc1.pdf", height=5, width=10)
plot(pvar)
dev.off()

#### find correlations between modules and environmental traits
# set up a list to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# calculate Pearson correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}

# convert numerical lables to colors for labeling of modules in the plot and color palette
MEColorNames = names(consMEs[[1]]$data);
my_palette <- colorRampPalette(c("#08519c", "white", "#cb181d"))(n = 199)

# initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
# find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative]);
# find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive]);
rownames(consensusCor)=rownames(moduleTraitCor[[1]])
colnames(consensusCor)=colnames(moduleTraitCor[[1]])
consensusCor=consensusCor[substring(rownames(moduleTraitCor[[1]]), 3) != 0,]

# plot consensus correlations from 2011 and 2012 (consensus positive in red, consensus negative in blue, no consensus between years in black)
pdf("Consensus.trait.env.plot.auto.wardD2.pdf", height=15, width=20)
heatmap.2(consensusCor,
          hclustfun = function(x) hclust(x, method = 'ward.D2'), 
          distfun = function(x) dist(x,method = 'euclidean'),
          Rowv = TRUE,
          na.color="black",
          main=paste("Consensus correlation"),
          cexCol=1 + 1/log10(dim(consensusCor)[1]),
          cexRow=1 + 1/log10(dim(consensusCor)[1]),
          margins = c(30, 20),
          symm=F,
          scale="column",
          trace="none",
          dendrogram="row",
          col = my_palette,
          key = TRUE,
          keysize = 0.6,
          density.info="none")
dev.off()

# separate module-trait correlation heatmaps per year                                                                                                                                                                                                                                                                                                                                        # plot correlation separately per year
moduleTraitCor.sub=moduleTraitCor[[1]][rownames(moduleTraitCor) != "ME0",]
names(moduleTraitCor)<-c("2011", "2012")
sets=c("2011", "2012")

for (set in sets){
  pdf(paste("Module.trait.plot.",set,".pdf", sep=''), height=15, width=20)
  heatmap.2(moduleTraitCor[[set]],
            na.color="white",
            main=paste("Module-trait relationship",set,sep=' '),
            cexCol=1 + 1/log10(dim(consensusCor)[1]),
            cexRow=1 + 1/log10(dim(consensusCor)[1]),
            margins = c(30, 20),
            symm=F,
            scale="column", 
            trace="none", 
            dendrogram="row",   
            col = my_palette,   
            key = TRUE,
            keysize = 0.6,
            density.info="none")
  dev.off()
}

#### generate barplots of OTUs within each consensus modules
# determine classification of modules according to consensus module-trait correlations
autumn-winter=c("ME11", "ME2","ME3", "ME20", "ME12", "ME21", "ME14", "ME1", "ME15")
early-bloom=c("ME10", "ME16")
spring-bloom=c("ME7", "ME19", "ME13", "ME25", "ME9", "ME17", "ME8", "ME22")
summer-bloom=c("ME5", "ME24", "ME18", "ME6")
late-bloom=c("ME23", "ME4")
no.module=c("ME0")
# place module classifications into a dataframe
df.trait <- data.frame(module=as.character(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)), 
                       trait_cor=c("no.module","autumn-winter","autumn-winter","autumn-winter","late-bloom", "summer-bloom", "summer-bloom",
                                   "spring-bloom", "spring-bloom", "spring-bloom", "early-bloom", "autumn-winter", "autumn-winter",
                                   "spring-bloom", "autumn-winter", "autumn-winter", "early-bloom", "spring-bloom", "summer-bloom",
                                   "spring-bloom", "autumn-winter", "autumn-winter", "spring-bloom", "late-bloom", "summer-bloom", "spring-bloom"))
# first generate a table of OTU attributes (i.e., taxonomy, abundances)
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
                 mean.2011 = commonOTUs2011.mean,
                 mean.2012 = commonOTUs2012.mean,
                 median.2011 = commonOTUs2011.med,
                 median.2012 = commonOTUs2012.med,
                 gmean.2011 = commonOTUs2011.gmean,
                 gmean.2012 = commonOTUs2012.gmean,
                 max.2011 = commonOTUs2011.max,
                 max.2012 = commonOTUs2012.max)
megamodule <- df.trait$trait_cor[match(nd$module, df.trait$module)]
nodeData <- cbind(nd, megamodule)
nodeData <- nodeData[order(as.numeric(nodeData$module)),]

# agglomerate OTU data (i.e., abundances) across each module for plotting
module <- sort(as.numeric(unique(nodeData$module)))
OF <- ddply(nodeData, .(module), summarize, length(node))
OTU_freq <- OF[order(as.numeric(sort(OF$module))),]
PO <- ddply(nodeData, .(module), summarize, round((length(node)/nGenes)*100, digits = 1))
pct_OTUs <- PO[order(as.numeric(sort(PO$module))),]
mean.2011 <- ddply(nodeData, .(module), summarize, sum(mean.2011))
mean_OTU_abund.2011 <- mean.2011[order(as.numeric(sort(mean.2011$module))),]
mean.2012 <- ddply(nodeData, .(module), summarize, sum(mean.2012))
mean_OTU_abund.2012 <- mean.2012[order(as.numeric(sort(mean.2012$module))),]
mean_OTU_abund <- rowMeans(cbind(mean_OTU_abund.2011$`sum(mean.2011)`, mean_OTU_abund.2012$`sum(mean.2012)`))
gmean.2011 <- ddply(nodeData, .(module), summarize, sum(gmean.2011))
gmean_OTU_abund.2011 <- gmean.2011[order(as.numeric(sort(gmean.2011$module))),]
gmean.2012 <- ddply(nodeData, .(module), summarize, sum(gmean.2012))
gmean_OTU_abund.2012 <- gmean.2012[order(as.numeric(sort(gmean.2012$module))),]
gmean_OTU_abund <- rowMeans(cbind(gmean_OTU_abund.2011$`sum(gmean.2011)`, gmean_OTU_abund.2012$`sum(gmean.2012)`))
max.2011 <- ddply(nodeData, .(module), summarize, sum(max.2011))
max_OTU_abund.2011 <- max.2011[order(as.numeric(sort(max.2011$module))),]
max.2012 <- ddply(nodeData, .(module), summarize, sum(max.2012))
max_OTU_abund.2012 <- max.2012[order(as.numeric(sort(max.2012$module))),]
max_OTU_abund <- rowMeans(cbind(max_OTU_abund.2011$`sum(max.2011)`, max_OTU_abund.2012$`sum(max.2012)`))
TC <- as.data.frame(df.trait)

# place module agglomerations into a dataframe
tdf=data.frame(module=module, OTU_freq=OTU_freq$`length(node)`, 
               pct_OTUs=pct_OTUs$`round((length(node)/nGenes) * 100, di...`, 
               mean_OTU_abund=mean_OTU_abund, 
               max_OTU_abund=max_OTU_abund,
               gmean_OTU_abund=gmean_OTU_abund, 
               trait_cor=TC$trait_cor) #[-1,]
tdf$module <- factor(tdf$module, levels = as.numeric(tdf$module))
tdf$trait_cor <- factor(tdf$trait_cor, levels = c("autumn-winter","early-bloom", "spring-bloom","summer-bloom", "late-bloom", "no.module"))

# extract module membership for each OTU
lmod=dlply(nodeData, c("module"), function(x) x=x)
lmod.df=data.frame(module=ldply(names(lmod), function(x) paste("ME",unique(lmod[[x]]$module), sep='')))

# make list of OTU ids for each module that are grep-compatible for plotting
lgrep=lapply(lmod, function(x) gsub("(\\|$)", "", capture.output(cat(gsub("(.*)", "^\\1_|", (x$node)), sep=''))))
names(lgrep)<-lmod.df$V1

# SELECT one megamodule (i.e., autumn-winter, early-bloom, spring-bloom, summer-bloom, late-bloom, no.module) at a time for plotting of 2011 and 2012 - change 'm' and 'mod' with desired module- run until dev.off()
# requires loading taxplot.functions.R for taxplot_subsets script
# e.g.
m="autumn-winter"
mod=autumn-winter

pdf(paste(m, ".otu.abund.2011-12.n10.pdf", sep=''), height=8, width=10*length(mod))
p2011=list()
p2012=list()

plotyear="2011"
p2011=lapply(
  mod, 
  function(i) taxplot_subnets(otu_matrix ="/Users/mchafee/Documents/TextFiles/FINAL.2010-12.MED/cogito.fl.2010_2012.med.m100.d1.silva1.3.txt", 
                              sample_data="/Users/mchafee/Documents/TextFiles/FINAL.2010-12.MED/FL.metadata.2010-12.FINAL.txt", 
                              taxa=lgrep[[i]], n=12, abund=0.0, main=names(lgrep)[which(names(lgrep)==i)], Y=as.numeric(plotyear), name="long", legend="right")
)
plotyear="2012"
p2012=lapply(
  mod,
  function(i) taxplot_subnets(otu_matrix ="/Users/mchafee/Documents/TextFiles/FINAL.2010-12.MED/cogito.fl.2010_2012.med.m100.d1.silva1.3.txt", 
                              sample_data="/Users/mchafee/Documents/TextFiles/FINAL.2010-12.MED/FL.metadata.2010-12.FINAL.txt", 
                              taxa=lgrep[[i]], n=12, abund=0.0, main=names(lgrep)[which(names(lgrep)==i)], Y=as.numeric(plotyear), name="long", legend="right")
)
do.call(grid.arrange, c(p2011, p2012, nrow=2))
dev.off()
