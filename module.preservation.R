### Module Preservation
library(scales)

# first run wgcna.R to obtain data for 2011 and 2012 formatted into a list structure (multiExpr)
# use same softPower as in wgcna.R to calculate adjacency matrices, which are generated during blockwiseConsensusModules but not saved
softPower = 15;
# initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# calculate adjacencies in each individual data set - TRY SPEARMAN?
for (set in 1:nSets)
  adjacencies[set, , ] = adjacency(multiExpr[[set]]$data, 
                                   selectCols = NULL, 
                                   type = "signed", 
                                   power = softPower,
                                   corFnc = "cor", corOptions = "use = 'p'")

dimnames(adjacencies)[[1]] <- c("2011", "2012")
dimnames(adjacencies)[[2]] <- as.vector(colnames(multiExpr[[1]]$data))
dimnames(adjacencies)[[3]] <- as.vector(colnames(multiExpr[[1]]$data))

adj2011=adjacencies["2011",,]
adj2012=adjacencies["2012",,]

# Place correlation adjacency data in to list of 2 datasets
nSets = 2;
setLabels = c("2011", "2012")
multiAdj = vector(mode = "list", length = nSets)
multiAdj[[1]] = list(data=as.matrix(adj2011))
colnames(multiAdj[[1]]$data) = rownames(adj2011)
rownames(multiAdj[[1]]$data) = rownames(adj2011)
multiAdj[[2]] = list(data=as.matrix(adj2012))
colnames(multiAdj[[2]]$data) = rownames(adj2012)
rownames(multiAdj[[2]]$data) = rownames(adj2012)
names(multiAdj)<-c("2011", "2012")
multiCol=list(net$colors, net$colors)

# run modulePreservation and set nPermutations to << 1000 to reduce runtime
mp=modulePreservation(multiAdj, multiCol,
                      greyName=0,
                      maxModuleSize=nGenes,
                      maxGoldModuleSize=mean(table(net$colors)),
                      referenceNetworks=1,
                      dataIsExpr=F, 
                      networkType='signed', 
                      nPermutations=1000, 
                      randomSeed=1,
                      quickCor = 0,
                      verbose=3)

save(mp,  file="modulePreservation.RData")

# isolate the observed statistics and their Z scores
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

# place preservation statistics into a dataframe (plotMods to remove unassinged grey module '0' and random gold module '0.1')
module = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]
plotMods = !(module %in% c("0", "0.1"));
mp.df=data.frame(modules=module[plotMods], 
                 moduleSizes=moduleSizes[plotMods], 
                 statsObs[, c("medianRank.pres", "medianRank.qual")][plotMods,],
                 signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)[plotMods,],
                 trait.cor=tdf$trait_cor[match(module, tdf$module)][plotMods],
                 gmean=tdf$gmean_OTU_abund[match(module, tdf$module)][plotMods],
                 max=tdf$max_OTU_abund[match(module, tdf$module)][plotMods])
megamodule_factor = c("autumn.winter","early.bloom", "spring.bloom","summer.bloom", "late.bloom")

c=c("autumn.winter"="#7570b3",
    "early.bloom"="#d95f02", 
    "spring.bloom"="#1b9e77", 
    "summer.bloom"="#386cb0", 
    "late.bloom"="#e41a1c", 
    "no.module"="#999999")

# plot Zsummary and medianRank preservation (Figure S9)
p1=ggplot(mp.df, aes(mp.df$moduleSizes, mp.df$Zsummary.pres)) +
  geom_point(aes(x=mp.df$moduleSizes, y=mp.df$Zsummary.pres, color=mp.df$trait.cor, size=mp.df$gmean, name="Geometric mean")) +
  geom_text(size=6, aes(x=mp.df$moduleSizes, y=mp.df$Zsummary.pres), label=mp.df$modules) +
  ylab("Zsummary preservation") +
  xlab("log2(Module size)") +
  scale_color_manual(values=alpha(c, 0.8), name="Megamodule") +
  scale_x_continuous(breaks=c(25,50,100,200,400,800,1600), trans=log2_trans()) +
  scale_size_continuous(name="Geometric mean", range=c(1,20)) +
  geom_hline(yintercept = c(2), linetype='longdash', colour=c("grey"), size=1.5) +
  geom_hline(yintercept = c(10), linetype='longdash', colour=c("black"), size=1.5)
p2=ggplot(mp.df, aes(mp.df$moduleSizes, mp.df$medianRank.pres)) +
  geom_point(aes(x=mp.df$moduleSizes, y=mp.df$medianRank.pres, color=mp.df$trait.cor, size=mp.df$gmean, name="Geometric mean")) +
  geom_text(size=6, aes(x=mp.df$moduleSizes, y=mp.df$medianRank.pres), label=mp.df$modules) +
  ylab("medianRank preservation") +
  xlab("log2(Module size)") +
  scale_color_manual(values=alpha(c, 0.8), name="Megamodule") +
  scale_x_continuous(breaks=c(25,50,100,200,400,800), trans=log2_trans()) +
  scale_y_reverse() +
  scale_size_continuous(name="Geometric mean", range=c(5,20)) +
  geom_hline(yintercept = c(8,16),linetype="longdash", colour=c("#386cb0","#1b9e77"), size=1.5)
plot_grid(p1,p2, nrow=1, ncol=2)

# Statistical tests between trait groups
sp=mp.df[c(which(mp.df$trait.cor=="spring.bloom")),]
su=mp.df[c(which(mp.df$trait.cor=="summer.bloom")),]
late=mp.df[c(which(mp.df$trait.cor=="late.bloom")),]
aw=mp.df[c(which(mp.df$trait.cor=="autumn.winter")),]
eb=mp.df[c(which(mp.df$trait.cor=="early.bloom")),]

# bloom vs nonbloom
bloom=rbind(eb, sp, su)
bloom$trait.cor <- "bloom"
nonbloom=rbind(aw, late)
nonbloom$trait.cor <- "nonbloom"
mp.df.agg=rbind(bloom, nonbloom)

# medianRank ttest bloom vs nonbloom
t.test(x=bloom[,3], y=nonbloom[,3])

# medianRank ttest spring vs summer
t.test(x=sp[,3], y=su[,3])






