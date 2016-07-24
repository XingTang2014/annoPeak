# annoPeakR
This application was developed to help biologists visualize their ChIP-seq/ChIP-exo results and generate or validate their hypothesis eventually. 
## Features
* This application allows users to generate multiple types plots to compare multiple ChIP-seq/ChIP-exo experiments. 
* Peak sets from ChIP-seq/ChIP-exo experiments can be uploaded in a bed format. 
* annoPeakR is constituted by 5 analysis modules: Peak associated gene structures, Peak size distribution, Peak to nearest peak distances, Overlapping peak identifications, Overlapping peak-associated genes.
* Tested R version to run this application is  

## Required R packages
* Shiny version 0.5.0
* GenomicRanges version 1.20.5
* VennDiagram version 1.6.9
* ggplot2 version 1.0.1
* RColorBrewer version 1.1-2
* reshape2 version 1.4.1 
* ChIPpeakAnno version 3.2.2 
* xtable version 1.7-4 
* biomaRt version 2.24.0 
* GSEABase version 1.30.2
* GO.db version 3.1.2 
* GOstats version 2.34.0
* TxDb.Mmusculus.UCSC.mm9.knownGene version 3.1.2
* TxDb.Mmusculus.UCSC.mm10.knownGene version 3.1.2 
* TxDb.Hsapiens.UCSC.hg18.knownGene version 3.1.2 
* TxDb.Hsapiens.UCSC.hg19.knownGene version 3.1.2 
* TxDb.Hsapiens.UCSC.hg38.knownGene version 3.1.2 
* org.Hs.eg.db version 3.1.2 
* org.Mm.eg.db version 3.1.2 

## Required R data 
* GeneSigDB [reactive](http://www.genesigdb.org/genesigdb/) 

## Installation
```r
install.packages(c("shiny", "ggplot2", "VennDiagram", "RColorBrewer", "reshape2", "xtable"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges", "ChIPpeakAnno", "biomaRt", "GSEABase", "GO.db", "GOstats", "TxDb.Mmusculus.UCSC.mm9.knownGene", "TxDb.Mmusculus.UCSC.mm10.knownGene", "TxDb.Hsapiens.UCSC.hg18.knownGene", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "org.Mm.eg.db"))
```
