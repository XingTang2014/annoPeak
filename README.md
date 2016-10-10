# annoPeak
This application was developed to help biologists visualize their ChIP-seq/ChIP-exo results and generate or validate their hypothesis eventually. To access this application, users can download this package and deploy it within R. Users can also access our online service from http://ccc-annopeak.osumc.edu/annoPeakR
## Features
* This application allows users to generate multiple types plots to compare multiple ChIP-seq/ChIP-exo experiments. 
* Peak sets from ChIP-seq/ChIP-exo experiments can be uploaded in a bed format. 
* annoPeak is constituted by 5 analysis modules: Peak associated gene structures, Peak size distribution, Peak to nearest peak distances, Overlapping peak identifications, Overlapping peak-associated genes.

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

Please download the right version of R for your own system from https://cloud.r-project.org/ and install it. 
RStudio is an integrated development environment (IDE) for R. Please download RStudio from https://www.rstudio.com/products/rstudio/#Desktop and install it.
* Install R packages within RStudio as below.
```r
install.packages(c("shiny", "ggplot2", "VennDiagram", "RColorBrewer", "reshape2", "xtable", "gplots"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges", "ChIPpeakAnno", "biomaRt", "GSEABase", "GO.db", "GOstats", "TxDb.Mmusculus.UCSC.mm9.knownGene", "TxDb.Mmusculus.UCSC.mm10.knownGene", "TxDb.Hsapiens.UCSC.hg18.knownGene", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "org.Mm.eg.db"))
```
* Download the source code for annoPeak from GitHub https://github.com/XingTang2014/annoPeak/.

Open ui.R from the download annoPeak folder. Click the button named "Run App" on the upper right corner of code editing window, an web page will be automatically invoked and the annoPeak application is ready to use. It may take a few minutes to load the required packages. 
