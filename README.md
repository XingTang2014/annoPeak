# annoPeakR
This application was developed to help biologists visualize their ChIP-seq/ChIP-exo results and generate or validate their hypothesis eventually. To access this application, users can download this package and deploy it within R. Users can also access our online service from http://ccc-annopeak.osumc.edu/annoPeakR
## Features
* This application allows users to generate multiple types plots to compare multiple ChIP-seq/ChIP-exo experiments. 
* Peak sets from ChIP-seq/ChIP-exo experiments can be uploaded in a bed format. 
* annoPeakR is constituted by 5 analysis modules: Peak associated gene structures, Peak size distribution, Peak to nearest peak distances, Overlapping peak identifications, Overlapping peak-associated genes.

## Required R packages
* Shiny
* GenomicRanges 
* VennDiagram
* ggplot2
* RColorBrewer
* reshape2
* ChIPpeakAnno 
* xtable
* biomaRt 
* GSEABase 
* GO.db 
* GOstats
* TxDb.Mmusculus.UCSC.mm9.knownGene
* TxDb.Mmusculus.UCSC.mm10.knownGene
* TxDb.Hsapiens.UCSC.hg18.knownGene
* TxDb.Hsapiens.UCSC.hg19.knownGene
* TxDb.Hsapiens.UCSC.hg38.knownGene
* org.Hs.eg.db 
* org.Mm.eg.db 

## Installation
* Install R and RStudio
* Install R packages within RStudio as below.
```r
install.packages(c("shiny", "ggplot2", "VennDiagram", "RColorBrewer", "reshape2", "xtable", "gplots"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges", "ChIPpeakAnno", "biomaRt", "GSEABase", "GO.db", "GOstats", "TxDb.Mmusculus.UCSC.mm9.knownGene", "TxDb.Mmusculus.UCSC.mm10.knownGene", "TxDb.Hsapiens.UCSC.hg18.knownGene", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "org.Mm.eg.db"))
```
