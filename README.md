### Identification of consensus modules between two OTU datsets using R package WGCNA

#### Input files: an OTU count matrix and sample metadata file

otu matrix: TableS4.txt

sample metadata: TableS4.metadata.txt

#### Workflow

##### pre_processing.R 
```
- Construct phyloseq object
- Calculate OTU singletons, doubletons and sample prevalence
- Save filtered abundance matrices where oligotypes have > 10% sample prevalence
```

##### wgcna.R (Langfelder & Horvath *BMC Bioinformatics* 2008)
```
- Import filtered abundance data from phyloseq object generated in pre_processing.R
- Import SparCC correlation matrixes
- Run WGCNA to idenfity consensus modules in SparCC networks
- Perform redundancy analysis (RDA)
- Plot oligotype abundance for each consensus module using taxplot_functions.R
- Filter oligotypes based on patterns graph deconstruction
- Intersect network graphs
```

##### taxplot_grep
```
- the taxplot_grep function uses a basic grep search to plot the abundance of oligotypes and taxonomic paths of 
   interest across 2010-2012 datasets

- source(taxplot_functions.R)
  usage: 
   otu_matrix = path to OTU table [TableS4.txt]
   sample_data = path to metadata table [TableS4.metadata.txt]
   taxa = taxonomic grep search - search feature must match a portion of the taxonomic path as in tax.txt 
   main = title 
   abund = minimum observed abundance in at least 1 sample
   position = c("stack", "fill") [stack = % total abundance, fill = scale abundance to 100%] 
   name = c("short", "long") [short = otuID only, long = otuID + lowest taxnomic level]
   legend = c("left", "right", "none")
   
 ```
#### Examples - taxa should be formatted according to data/taxonomy.txt for successful grep search

##### Class search: 
```
taxplot_grep(otu_matrix="data/TableS4.txt", sample_data="data/TableS4.metadata.txt", 
   main = "Alphaproteobacteria", abund = 0.05, taxa = "Alphaproteobacteria", legend = "right", name="long", 
   position="stack")
  ```
  ![alt text](https://github.com/genomewalker/medNS/blob/master/plot_examples/Alphaproteobacteria.png)
  
##### Genera search: 
```
taxplot_grep(otu_matrix="data/TableS4.txt", sample_data="data/TableS4.metadata.txt", 
   main = "Polaribacter", abund = 0.05, taxa = "Polaribacter", legend = "bottom", name="short", 
   position="stack")
```
  ![alt text](https://github.com/genomewalker/medNS/blob/master/plot_examples/Polaribacter.png)

##### OTU search: 
```
taxplot_grep(otu_matrix="data/TableS4.txt", sample_data="data/TableS4.metadata.txt", main = "BD1_7 clade", 
   abund = 0.0, taxa = "^11259_", legend = "bottom", name="short", position="stack")
```
 ![alt text](https://github.com/genomewalker/medNS/blob/master/plot_examples/BD1_7_clade.png)
