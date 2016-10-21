### Identification of consensus modules between two OTU datsets using R package WGCNA (Langfelder & Horvath *BMC Bioinformatics* 2008)

#### Input files: an OTU count matrix and sample metadata file

otu matrix: TableS3.OTU_matrix.txt

sample metadata: TableS4.sample.metadata.txt

#### Workflow

##### OTU.pre_processing.R 
```
- Construct phyloseq object
- Calculate OTU singletons, doubletons and sample prevalence
- Filter OTUs < 10% sample prevalence
- Evaluate transformation of OTU matrix - Hellinger transformation of log scaled abundance (+ pseudocount 1)
```

##### wgcna.R
```
- Import filtered data from phyloseq object generated in OTU.pre_processing.R
- Collect basic abundance measures of OTUs (mean, geometric mean, maximum observed abundance)
- Perform OTU matrix transformation - Hellinger transformation of log scaled abundance (+ pseudocount 1)

- Run WGCNA 
- identify optimal softPower to achieve a scale free topology (scale free toplogy model fit R squared > 0.8)
   - run blockwiseConsensusModules
   - examine variance explained by eigenvector (PC1)
   - find correlations and consensus correlations of module eigengenes and sample environmental metadata (Figure 5)
   - identify 'megamodules' accoriding to hierarchical clustering patterns
   - assign module, megamodule, taxonomic and abundance data to each OTU in a dataframe
   - plot OTU barblots within each megamodule using taxplot_subnets (Figures 6, 7, S3-S7)
```

##### module.preservation.R
```
- calculate module preservation using Zsummary statistic and medianRank (Figure S9)
```
### OTU abundance barplots using taxplot_grep in taxplot.functions.R
```
- the taxplot_grep function uses a basic grep search to plot the abundance of OTUs and taxonomic paths of 
   interest across 2010-2012

- source(taxplot.function.R)
  usage: 
   otu_matrix = path to OTU table [TableS3.OTU_matrix.txt]
   sample_data = path to metadata table [TableS4.sample.metadata.txt]
   main = title 
   abund = minimum observed abundance in at least 1 sample
   taxa = taxonomic grep search - search feature must match a portion of the taxonomic path as in tax.txt 
   legend position = c("left", "right", "none")
   name = c("short", "long") [short = otuID only, long = otuID + lowest taxnomic level]
   position = c("stack", "fill") [stack = % total abundance, fill = scale abundance to 100%] 
 ```
#### Examples - taxa should be formatted according to data/taxonomy.txt for successful grep search

##### Class search: 
```
taxplot_grep(otu_matrix="TableS3.OTU_matrix.txt", sample_data="TableS4.sample.metadata.txt", 
   main = "Alphaproteobacteria", abund = 0.05, taxa = "Alphaproteobacteria", legend = "right", name="long", 
   position="stack")
  ```
  ![alt text](https://github.com/genomewalker/medNS/blob/master/plot_examples/Alphaproteobacteria.png)
  
##### Genera search: 
```
taxplot_grep(otu_matrix="TableS3.OTU_matrix.txt", sample_data="TableS4.sample.metadata.txt", 
   main = "Polaribacter", abund = 0.05, taxa = "Polaribacter", legend = "bottom", name="short", 
   position="stack")
```
  ![alt text](https://github.com/genomewalker/medNS/blob/master/plot_examples/Polaribacter.png)

##### OTU search: 
```
taxplot_grep(otu_matrix="TableS3.OTU_matrix.txt", sample_data="TableS4.sample.metadata.txt", main = "BD1_7 clade", 
   abund = 0.0, taxa = "^11259_", legend = "bottom", name="short", position="stack")
```
 ![alt text](https://github.com/genomewalker/medNS/blob/master/plot_examples/BD1_7_clade.png)
