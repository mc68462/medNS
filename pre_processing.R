library(ggplot2)
library(reshape)
library(dplyr)
library(phyloseq)
library(cowplot)

#COGITO data
# from 2010-2012 (142 samples)
cogito.otus<-read.table(file="data/TableS4.txt", sep="\t", header=T, row.names = 1)

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
cogito.sample.metadata<-read.table(file="data/TableS4_metadata.txt" , sep="\t", header=T, row.names=1)
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

plot_grid(p1, p2)

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

plot_grid(p1, p2)

# Remove OTUs with < 10% Prevalence
cogito_physeq_2011_prevalence_filt <- cogito_physeq_2011_prevalence %>%
  filter(prev >= cogito_physeq_2011_n_samples * 0.1)
summary(cogito_physeq_2011_prevalence_filt$prev)
cogito_physeq_2011_otu_sample_physeq_filt_prev <- prune_taxa(cogito_physeq_2011_prevalence_filt %>% .$otu_name %>% as.vector %>% as.character, cogito.physeq.2011)

cogito_physeq_2012_prevalence_filt <- cogito_physeq_2012_prevalence %>%
  filter(prev >= cogito_physeq_2012_n_samples * 0.1)
summary(cogito_physeq_2012_prevalence_filt$prev)
cogito_physeq_2012_otu_sample_physeq_filt_prev <- prune_taxa(cogito_physeq_2012_prevalence_filt %>% .$otu_name %>% as.vector %>% as.character, cogito.physeq.2012)

saveRDS(cogito_physeq_2012_otu_sample_physeq_filt_prev, file = "results/filtered.phyloseq.2012.rds")
saveRDS(cogito_physeq_2011_otu_sample_physeq_filt_prev, file = "results/filtered.phyloseq.2011.rds")
