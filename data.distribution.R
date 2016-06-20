library(DESeq2)
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
library(dplyr)
library(vsn)
library(metagenomeSeq)
# set working directory
workingDir = "/Users/mchafee/Documents/Local_data_for_Rstudio/WGCNA/final/"
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#COGITO data
# from 2010-2012 (142 samples)
cogito.otus<-read.table(file="/Users/mchafee/Documents/TextFiles/FINAL.2010-12.MED/cogito.fl.2010_2012.med.m100.d1.silva1.3.txt", sep="\t", header=T, row.names = 1)

#massage cogito data, split node and taxonomy
cogito.otus.taxo<-as.data.frame(stringr::str_split_fixed(rownames(cogito.otus), "_", 2))
rownames(cogito.otus)<-cogito.otus.taxo[,1]
cogito.otus.taxo<-as.matrix((stringr::str_split_fixed(cogito.otus.taxo[,2], ";", 6)))
rownames(cogito.otus.taxo)<-rownames(cogito.otus)
colnames(cogito.otus.taxo)<-c("domain",
                              "phylum",
                              "class",
                              "order",
                              "family",
                              "genus")
cogito.otus.taxo[is.na(cogito.otus.taxo)]<-NA

#read metadata file 
cogito.sample.metadata<-read.table(file="/Users/mchafee/Documents/TextFiles/FINAL.2010-12.MED/FL.metadata.2010-12.FINAL.txt" , sep="\t", header=T, row.names=1)
sample_name=rownames(cogito.sample.metadata)
cogito.sample.metadata=cbind(sample_name, cogito.sample.metadata)
# merge all into physeq object
cogito.otus.taxo<-tax_table(cogito.otus.taxo)
cogito.otus.counts<-otu_table(cogito.otus,taxa_are_rows = TRUE)
cogito.sample.metadata<-sample_data(cogito.sample.metadata)
cogito.physeq<-phyloseq(cogito.otus.counts,cogito.otus.taxo)
cogito.physeq<-merge_phyloseq(cogito.physeq,cogito.sample.metadata)

# filter organelles
cogito.physeq.no.chloro.mito.norel<-subset_taxa(cogito.physeq, class!="Chloroplast" & family!="Mitochondria" & domain!="NoRelative")
cogito.physeq<-cogito.physeq.no.chloro.mito.norel

#Split data set by year
cogito.physeq.2011<-subset_samples(cogito.physeq, year=="2011")
cogito.physeq.2012<-subset_samples(cogito.physeq, year=="2012")

# remove OTU rows that sum to 0
#cogito.physeq.2011.filt=filter_taxa(cogito.physeq.2011, function(x) sum(x)>0, prune=T)
#cogito.physeq.2012.filt=filter_taxa(cogito.physeq.2012, function(x) sum(x)>0, prune=T)

# Look sparsity
pdf("count.distribution.sum.counts.vs.max.in.sample.pdf", height=5, width=10)
par(mfrow=c(1,2))

physeq_object <- cogito.physeq.2011
deseq_object <- phyloseq_to_deseq2(physeq_object, ~1)
rs <- rowSums(counts(deseq_object))
rmx <- apply(counts(deseq_object), 1, max)
plot(rs+1, rmx/rs, log="x", main="2011", xlab="Sum of counts per OTU", ylab="Max count/sum count")

physeq_object <- cogito.physeq.2012
deseq_object <- phyloseq_to_deseq2(physeq_object, ~1)
rs <- rowSums(counts(deseq_object))
rmx <- apply(counts(deseq_object), 1, max)
df <- data.frame(rs, rmx)
plot(rs+1, rmx/rs, log="x", main="2012", xlab="Sum of counts per OTU", ylab="Max count/sum count")

dev.off()

# OTU summary
cogito_physeq_2011_summary <- melt(otu_table(cogito.physeq.2011)) %>% 
  tbl_df %>% magrittr::set_colnames(c("otu_name","label","count")) %>% 
  group_by(label) %>% mutate(rel_abun = count/sum(count)) %>% 
  arrange(otu_name, label) %>% 
  filter(count > 0)

cogito_physeq_2012_summary <- melt(otu_table(cogito.physeq.2012)) %>% 
  tbl_df %>% magrittr::set_colnames(c("otu_name","label","count")) %>% 
  group_by(label) %>% mutate(rel_abun = count/sum(count)) %>% 
  arrange(otu_name, label) %>% 
  filter(count > 0)

# Prevalence, singletons and doubletons
cogito_physeq_2011_prevalence <- cogito_physeq_2011_summary %>% dplyr::select(otu_name, count) %>% group_by(otu_name) %>% summarise(prev = sum(count >0), total_counts = sum(count)) 
cogito_physeq_2011_abs_singletons <- cogito_physeq_2011_prevalence %>% dplyr::filter(total_counts <= 1, prev <= 1)
cogito_physeq_2011_abun_singletons <- cogito_physeq_2011_prevalence %>% dplyr::filter(total_counts > 1, prev <= 1)
cogito_physeq_2011_abs_singletons_names <- cogito_physeq_2011_abs_singletons$otu_name
cogito_physeq_2011_abun_doubletons <- cogito_physeq_2011_prevalence %>% filter(!(otu_name %in% cogito_physeq_2011_abs_singletons_names)) %>% dplyr::filter( prev <= 2)

cogito_physeq_2011_prevalence <- cogito_physeq_2011_prevalence %>% filter(!(otu_name %in% cogito_physeq_2011_abs_singletons_names))
cogito_physeq_2011_counts <- data.frame(class = c("Total OTUs","Absolute singletons", "Abundant singletons"), counts = c(length(unique(cogito_physeq_2011_summary$otu_name)) ,dim(cogito_physeq_2011_abs_singletons)[1], dim(cogito_physeq_2011_abun_singletons)[1]))
print(cogito_physeq_2011_counts)

cogito_physeq_2011_counts$class <- factor(cogito_physeq_2011_counts$class, levels=c("Total OTUs","Absolute singletons", "Abundant singletons"))
p1=ggplot(cogito_physeq_2011_counts, aes(class, counts)) +
  geom_bar(stat = "identity") + theme_bw() +
  xlab("") + ylab("# OTUs") + ggtitle("2011")

cogito_physeq_2012_prevalence <- cogito_physeq_2012_summary %>% dplyr::select(otu_name, count) %>% group_by(otu_name) %>% summarise(prev = sum(count >0), total_counts = sum(count)) 
cogito_physeq_2012_abs_singletons <- cogito_physeq_2012_prevalence %>% dplyr::filter(total_counts <= 1, prev <= 1)
cogito_physeq_2012_abun_singletons <- cogito_physeq_2012_prevalence %>% dplyr::filter(total_counts > 1, prev <= 1)
cogito_physeq_2012_abs_singletons_names <- cogito_physeq_2012_abs_singletons$otu_name
cogito_physeq_2012_abun_doubletons <- cogito_physeq_2012_prevalence %>% filter(!(otu_name %in% cogito_physeq_2012_abs_singletons_names)) %>% dplyr::filter( prev <= 2)

cogito_physeq_2012_prevalence <- cogito_physeq_2012_prevalence %>% filter(!(otu_name %in% cogito_physeq_2012_abs_singletons_names))
cogito_physeq_2012_counts <- data.frame(class = c("Total OTUs","Absolute singletons", "Abundant singletons"), counts = c(length(unique(cogito_physeq_2012_summary$otu_name)) ,dim(cogito_physeq_2012_abs_singletons)[1], dim(cogito_physeq_2012_abun_singletons)[1]))
print(cogito_physeq_2012_counts)

cogito_physeq_2012_counts$class <- factor(cogito_physeq_2012_counts$class, levels=c("Total OTUs","Absolute singletons", "Abundant singletons"))
p2=ggplot(cogito_physeq_2012_counts, aes(class, counts)) +
  geom_bar(stat = "identity") + theme_bw() +
  xlab("") + ylab("# OTUs") + ggtitle("2012")

pdf("singles.doubles.pdf", height=5, width=10)
plot_grid(p1, p2)
dev.off()

# Evaluate prevalence cut-off
cogito_physeq_2011_n_samples <- dim(sample_data(cogito.physeq.2011))[1]
cogito_physeq_2011_prevalence_dist <- plyr::ldply(seq(0.00,1,0.1), function (X) {
  cogito_physeq_2011_prevalence %>% 
    filter(prev >= cogito_physeq_2011_n_samples * X) %>%
    summarise(N=n()) %>% mutate(N = N, prev=paste(100*X, "%", sep = ""), nsamples = round(cogito_physeq_2011_n_samples * X))
})
cogito_physeq_2011_prevalence_dist$prev <- factor(cogito_physeq_2011_prevalence_dist$prev , levels=paste(100*seq(0.00,1,0.1), "%", sep = ""))
p1 <- ggplot(cogito_physeq_2011_prevalence_dist, aes(prev, N, group=1)) +
  geom_line() + geom_point()+ theme_bw() +
  xlab("Prevalence") + ylab("# OTUs") + ggtitle("2011")

cogito_physeq_2012_n_samples <- dim(sample_data(cogito.physeq.2012))[1]
cogito_physeq_2012_prevalence_dist <- plyr::ldply(seq(0.00,1,0.1), function (X) {
  cogito_physeq_2012_prevalence %>% 
    filter(prev >= cogito_physeq_2012_n_samples * X) %>%
    summarise(N=n()) %>% mutate(N = N, prev=paste(100*X, "%", sep = ""), nsamples = round(cogito_physeq_2012_n_samples * X))
})
cogito_physeq_2012_prevalence_dist$prev <- factor(cogito_physeq_2012_prevalence_dist$prev , levels=paste(100*seq(0.00,1,0.1), "%", sep = ""))
p2 <- ggplot(cogito_physeq_2012_prevalence_dist, aes(prev, N, group=1)) +
  geom_line() + geom_point()+ theme_bw() +
  xlab("Prevalence") + ylab("# OTUs") + ggtitle("2012")

pdf("prevalence.pdf", height=5, width=10)
plot_grid(p1, p2)
dev.off()

# Filter OTUs based on prevalence - 10% cutoff
cogito_physeq_2011_prevalence_filt <- cogito_physeq_2011_prevalence %>% 
  filter(prev >= cogito_physeq_2011_n_samples * 0.1)
summary(cogito_physeq_2011_prevalence_filt$prev)
cogito_physeq_2011_otu_sample_physeq_filt_prev <- prune_taxa(cogito_physeq_2011_prevalence_filt %>% .$otu_name %>% as.vector %>% as.character, cogito.physeq.2011)

cogito_physeq_2012_prevalence_filt <- cogito_physeq_2012_prevalence %>% 
  filter(prev >= cogito_physeq_2012_n_samples * 0.1)
summary(cogito_physeq_2012_prevalence_filt$prev)
cogito_physeq_2012_otu_sample_physeq_filt_prev <- prune_taxa(cogito_physeq_2012_prevalence_filt %>% .$otu_name %>% as.vector %>% as.character, cogito.physeq.2012)

# Plot heteroscedasticity
# larger average expression have on average larger observed variances across samples

# No transformation - rows are samples
filt2011=t(as.matrix(as.data.frame(otu_table(cogito_physeq_2011_otu_sample_physeq_filt_prev))))
filt2012=t(as.matrix(as.data.frame(otu_table(cogito_physeq_2012_otu_sample_physeq_filt_prev))))

# Log Hellinger
# decostand tranformations : 
# All methods have a default margin. MARGIN=1 means rows (sites in a normal data set) and MARGIN=2 means columns (species in a normal data set).
hellog2011 <- decostand(log(filt2011 + 1), method="hellinger", MARGIN=1)
hellog2012 <- decostand(log(filt2012 + 1), method="hellinger", MARGIN=1)

# Variance stabilization
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

filt2011.deseq <- phyloseq_to_deseq2(cogito_physeq_2011_otu_sample_physeq_filt_prev, ~1)

geoMeans <- apply(counts(filt2011.deseq), 1, gm_mean)
diagdds <- DESeq2::estimateSizeFactors(filt2011.deseq, geoMeans=geoMeans)
diagdds <- DESeq2::estimateDispersions(diagdds)
vst <- DESeq2::varianceStabilizingTransformation(diagdds, blind = FALSE)
vstMat2011 <- GenomicRanges::assay(vst)
vstMat2011[vstMat2011 < 0] <- 0

filt2012.deseq <- phyloseq_to_deseq2(cogito_physeq_2012_otu_sample_physeq_filt_prev, ~1)

geoMeans <- apply(counts(filt2012.deseq), 1, gm_mean)
diagdds <- DESeq2::estimateSizeFactors(filt2012.deseq, geoMeans=geoMeans)
diagdds <- DESeq2::estimateDispersions(diagdds)
vst <- DESeq2::varianceStabilizingTransformation(diagdds, blind = FALSE)
vstMat2012 <- GenomicRanges::assay(vst)
vstMat2012[vstMat2012 < 0] <- 0

# meanSDplot: Standard deviation and mean are calculated row-wise from the expression matrix (in) x. 
pdf("transformation.comparison.pdf", height=10, width=15)
par(mfrow=c(2,3))
meanSdPlot(t(filt2011),main="Raw counts - 2011", ylab = "SD", xlab="Rank(mean)")
meanSdPlot(t(hellog2011),main="Hellinger of log(counts) + 1 - 2011", ylab = "SD", xlab="Rank(mean)")
meanSdPlot(vstMat2011, main="VST - 2011", ylab = "SD", xlab="Rank(mean)")

meanSdPlot(t(filt2012),main="Raw counts - 2011", ylab = "SD", xlab="Rank(mean)")
meanSdPlot(t(hellog2012),main="Hellinger of log(counts) + 1 - 2011", ylab = "SD", xlab="Rank(mean)")
meanSdPlot(vstMat2012, main="VST - 2011", ylab = "SD", xlab="Rank(mean)")
dev.off()
