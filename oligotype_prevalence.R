library(ggplot2)

#import OTU matrix
df=read.table("data/Table_S3.OTU_matrix.txt", header=TRUE, row.names='nodes', sep='\t')
OTU_table=df[rownames(df)[grep("Chloroplast|Mitochondria|NoRelative", rownames(df), ignore.case=TRUE, invert=TRUE)],]
dim(OTU_table)

OTU.table.all <- OTU_table

# Endemism vs cosmopolitanism
library(ggplot2)

# set min abundance for count calculations (use all for plotting)
minabu <- 500
otus <- t(OTU.table.all)

# We count for each OTU in how many sites are present
site.otu <- apply(otus > minabu, 2, sum)
# Take percent total sample prevalence
site.otu.prev <- apply(otus > minabu, 2, mean)

# mean
if (minabu == 0) {
  # If minabu = 0 we can calculate the mean straightforward
  mean.abu <- apply(otus, 2, sum)/site.otu
} else {
  mean.abu <- rep(0, ncol(otus))
  # If not we have to find where there are the OTU with a value larger than minabu
  for (i in 1:ncol(otus)) {
    mask <- otus[, i] > minabu
    # And calculate the mean abundance
    mean.abu[i] <- sum(otus[mask, i])/max(1, site.otu[i])
  }
}
mean.abu[is.na(mean.abu)] <- 0


b<-as.data.frame(cbind(site.otu.prev[mean.abu > minabu], mean.abu[mean.abu > minabu]))
colnames(b)<-c("counts", "mean")
b$com[b$counts >= 0.75 & b$mean >= minabu]<-"Generalist"
b$com[b$counts <= 0.1  & b$mean >= minabu]<-"Specialist"
b$com[is.na(b$com)]<-"NA"

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

pdf("oligotype_prevalence.pdf", height=6, width=7)
ggplot(data=b) +
  geom_point(aes(counts, mean, colour=com)) +
  xlab("% Sample prevalence") +
  ylab("Mean abundance per sample") +
  scale_y_continuous(trans = log_trans(), breaks = base_breaks()) +
  scale_color_manual(values=c("red", "grey70", "blue")) +
  guides(colour=FALSE) +
  theme(legend.position = "bottom") +
  theme_bw()
dev.off()

#Return generalist & specialsits oligotypes
broad <- rownames(b)[b$counts >= 0.75 & b$mean >= minabu]

narrow <- rownames(b)[b$counts <= 0.1  & b$mean >= minabu]
