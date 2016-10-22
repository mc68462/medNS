# plotting functions based on grep searches of taxnoomy and lists of OTU ids

taxplot_grep <- function(otu_matrix, sample_data, taxa, main, abund, position=c("stack", "fill"), name=c("long", "short"), legend=c("right", "bottom", "none")){
  # import matrix
  df <- read.table(otu_matrix, header=TRUE, row.names='nodes', sep='\t')
  # import sample metadata
  md <- read.table(sample_data, header=TRUE, sep='\t', as.is=TRUE)

  # remove organelle sequences
  df <- df[rownames(df)[grep("Chloroplast|Mitochondria|NoRelative", rownames(df),  ignore.case=TRUE, invert=TRUE)],]
  dfRA <- scale(df, center=FALSE, scale=colSums(df))
  dfRA <- as.data.frame(dfRA)
  
  # subset based on minimum observed value (i.e., minimum abundance in a single sample required for plotting)
  dfRA$max <- apply(dfRA, 1, max)
  dfRAmax <- subset(dfRA, dfRA$max>=abund)
  dfRAmax <- t(dfRAmax[,1:length(colnames(dfRA))-1])
  
  # format names for use in for loop with grep to remove (), spaces and -
  colnames(dfRAmax) <- gsub("\\(|\\)","", colnames(dfRAmax))
  colnames(dfRAmax) <- gsub("-","_", colnames(dfRAmax))
  colnames(dfRAmax) <- gsub(" ","_", colnames(dfRAmax))
  colnames(dfRAmax) <- gsub(";uncultured","", colnames(dfRAmax))
  
  # set long or short format for OTU legend names
  if (name == "short"){
    target_oligotypes <- gsub("^(.*?)(_.*;)(.*$)", "\\1",grep(taxa, colnames(dfRAmax), value=TRUE))  #oligotype id only
} else{ 
    if (name == "long"){
      target_oligotypes <- gsub("^(.*?)(_.*;)(.*$)", "\\1 \\3",grep(taxa, colnames(dfRAmax), value=TRUE)) #long name
      }else{
        stop("name argument must == 'short' or 'long'")
      }
}

  # shade grey plot areas with spring and summer Julian days
  rect.spring <- data.frame(xmin=79, xmax=172, ymin=-Inf, ymax=Inf)                      
  rect.summer <- data.frame(xmin=172, xmax=265, ymin=-Inf, ymax=Inf)  
  
  # make plot using grep search results
  dfm <- melt(cbind(md,dfRAmax)[,c("julian_day", "year", colnames(dfRAmax)[grep(taxa, colnames(dfRAmax))])], id.vars=c("julian_day", "year"))  
  colourCount <- length(target_oligotypes)
  ggplot(dfm, aes(x=julian_day, y=value, fill=variable))  +    
    geom_bar(stat="identity", position=position, width=2, colour="black", size=0.1) +   
    theme_bw() %+replace% theme(panel.background=element_rect(fill=NA)) +               
    theme(legend.position=legend) +
    ggtitle(main) +
    guides(fill=guide_legend(keywidth = 1, keyheight = 1, override.aes=list(colour=NULL), title="Oligotype")) +  
    theme(plot.title = element_text(size=12,family="Helvetica", face="italic")) +
    theme(legend.text=element_text(size=10,family="Helvetica", face="plain")) +               
    theme(legend.title = element_text(size=12,family="Helvetica", face="plain")) +
    theme(axis.title = element_text(size=12,family="Helvetica", face="plain")) +
    theme(axis.text = element_text(size=8, family="Helvetica", face="plain")) +
    theme(axis.title.y=element_text(vjust=1)) +
    theme(axis.title.x=element_text(vjust=-0.2)) +
    theme(strip.text.y = element_text(size=12)) +
    theme(strip.text.x = element_text(size=12)) +
    scale_fill_manual(labels=target_oligotypes, values=rainbow(colourCount))  +   
    xlab("Julian Day") +    
    ylab("Frequency") +
    theme(panel.grid.minor=element_blank()) +     
    theme(panel.grid.major=element_blank()) +     
    theme(panel.border=element_rect(color="black", size=0.2))  + 
    geom_rect(data=rect.spring, aes (xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.15, inherit.aes = FALSE)  +
    geom_rect(data=rect.summer, aes (xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1,  inherit.aes = FALSE) +
    facet_grid(year ~ ., scale="free_y") +
    theme(plot.background = element_rect(fill = "transparent",colour = NA))
}


taxplot_subnets <- function(otu_matrix, sample_data, taxa, Year, main, n, abund, name=c("long", "short"), legend=c("right", "bottom", "none")){
                   
  # import matrix
  df <- read.table(otu_matrix, header=TRUE, row.names='nodes', sep='\t')
  # import sample metadata
  md <- read.table(sample_data, header=TRUE, sep='\t', as.is=TRUE)
  
  # remove organelle sequences
  df <- df[rownames(df)[grep("Chloroplast|Mitochondria|NoRelative", rownames(df),  ignore.case=TRUE, invert=TRUE)],]
  dfRA <- scale(df, center=FALSE, scale=colSums(df))
  dfRA <- as.data.frame(dfRA)
  
  # subset based on minimum observed value (i.e., minimum abundance in a single sample required for plotting)
  dfRA$max <- apply(dfRA, 1, max)
  dfRAmax <- subset(dfRA, dfRA$max>=abund)
  dfRAmax <- t(dfRAmax[,1:length(colnames(dfRA))-1])
  
  # format names for use in for loop with grep to remove (), spaces and -
  colnames(dfRAmax) <- gsub("\\(|\\)","", colnames(dfRAmax))
  colnames(dfRAmax) <- gsub("-","_", colnames(dfRAmax))
  colnames(dfRAmax) <- gsub(" ","_", colnames(dfRAmax))
  colnames(dfRAmax) <- gsub(";uncultured","", colnames(dfRAmax))
  
  # set long or short format for OTU legend names
  if (name == "short"){
    target_oligotypes <- gsub("^(.*?)(_.*;)(.*$)", "\\1",grep(taxa, colnames(dfRAmax), value=TRUE))  #oligotype id only
  } else{ 
    if (name == "long"){
      target_oligotypes <- gsub("^(.*?)(_.*;)(.*$)", "\\1 \\3",grep(taxa, colnames(dfRAmax), value=TRUE)) #long name
    }else{
      stop("name argument must == 'short' or 'long'")
    }
  }
  
  # shade grey plot areas with spring and summer Julian days
  rect.spring <- data.frame(xmin=79, xmax=172, ymin=-Inf, ymax=Inf)                      
  rect.summer <- data.frame(xmin=172, xmax=265, ymin=-Inf, ymax=Inf)  
  
  dfc=cbind(md,dfRAmax)
  dfs=subset(dfc, dfc$year==Year)[,c("julian_day", colnames(dfRAmax)[grep(taxa, colnames(dfRAmax))[1:n]])]
  dfm=melt(dfs, id.vars=c("julian_day"))  
  
  colourCount=length(colnames(dfRAmax)[grep(taxa, colnames(dfRAmax))[1:n]])
  ggplot(dfm, aes(x=julian_day, y=value, fill=variable))  +    
    geom_bar(stat="identity", position="stack", width=2, colour="black", size=0.1) +   
    theme_bw() %+replace% theme(panel.background=element_rect(fill=NA)) +               
    theme(legend.position=legend) +
    ggtitle(paste(main, Year)) +
    guides(fill=guide_legend(keywidth = 1, keyheight = 1, override.aes=list(colour=NULL), title="Oligotype")) +  
    theme(plot.title = element_text(size=16,family="Helvetica", face="plain")) +
    theme(legend.text=element_text(size=11,family="Helvetica", face="plain")) +               
    theme(legend.title = element_text(size=11,family="Helvetica", face="plain")) +
    theme(axis.title = element_text(size=12,family="Helvetica", face="plain")) +
    theme(axis.text = element_text(size=9, family="Helvetica", face="plain")) +
    theme(axis.title.y=element_text(vjust=1)) +
    theme(axis.title.x=element_text(vjust=-0.2)) +
    theme(strip.text.y = element_text(size=10)) +
    theme(strip.text.x = element_text(size=10)) +
    scale_fill_manual(labels=target_oligotypes, values=rainbow(colourCount))  +   
    xlab("Julian Day") +
    ylab("Frequency") +
    theme(panel.grid.minor=element_blank()) +     
    theme(panel.grid.major=element_blank()) +     
    theme(panel.border=element_rect(color="black", size=0.2))  + 
    geom_rect(data=rect.spring, aes (xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.15, inherit.aes = FALSE)  +
    geom_rect(data=rect.summer, aes (xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1,  inherit.aes = FALSE) +
    theme(plot.background = element_rect(fill = "transparent",colour = NA)) 
  
}
