


run_GSEA <- function(genes, genome){
  if( genome %in% c("mm8", "mm9", "mm10")){ 
    specie <- "Mouse"
    ensembl = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
  }
  
  if( genome %in% c("hg18", "hg19", "hg38")){ 
    specie <- "Human"
    ensembl = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl")
  }
  
  adj.cutoff <- 0.05
  
  # Prepare gene list: map gene symbol to ensembl id
  symbol2ensembl <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id'), 
        values = genes, filters='external_gene_name', 
        mart = ensembl, uniqueRows=TRUE)
  all.genes <- unique(getBM(attributes=c('ensembl_gene_id'), 
                            mart = ensembl)[, 'ensembl_gene_id'])
  gene_list <- unique(symbol2ensembl[, 'ensembl_gene_id'])
  
  # prepare gene signatures: select gene signature with 10 more genes in the selected genome
  selectedSig <- unlist(lapply(genesigdb@.Data, function(x){x@organism})) == specie
  selectedSig <- genesigdb[selectedSig]
  ok <- sapply(geneIds(selectedSig), length)>10
  selectedSig <- selectedSig[ok]
  
  # start test
  res <- t(as.data.frame(lapply(geneIds(selectedSig), function(SigGenes){
    gene_list.in <- length( intersect( gene_list, SigGenes) )
    gene_list.not <- length( setdiff( gene_list, SigGenes) )
    univ.in <- length( intersect( all.genes, SigGenes) )
    univ.not <- length( setdiff( all.genes, SigGenes) )
    test.res <- fisher.test(matrix(c(gene_list.in, gene_list.not, univ.in, univ.not), nrow=2, byrow=TRUE), alternative="greater")
    c(test.res$p.value, test.res$estimate) 
  })))
  colnames(res) <- c("Pvalue", "OddsRatio")
  res <- data.frame(SigID=rownames(res), res, pubMedIds=unlist(lapply(selectedSig@.Data, function(x){x@pubMedIds})))
   
  Size <- unlist(lapply(geneIds(selectedSig), length))
  Count <- unlist(lapply(geneIds(selectedSig), function(x){ 
    g <- intersect(gene_list, x)
    s <- unique(symbol2ensembl[ symbol2ensembl[,2] %in% g, 1])
    length( s ) } ) )
  ExpCount <- format(length(gene_list) / length(all.genes) * Size, digits=1,nsmall=2)
  enrichedGenesInSig <- unlist(lapply( geneIds(selectedSig), function(x){ 
    g <- intersect(x, gene_list)
    s <- unique(symbol2ensembl[ symbol2ensembl[,2] %in% g, 1])
    paste(s, collapse=";") } ))
  BH.adj.Pvalue <- p.adjust(res[,"Pvalue"], method="BH")
  # columns of res include: pvalue, oddsratio, pmid, expcount, count, size, padjust, genes
  res <- data.frame(res, ExpCount=ExpCount, Count=Count, Size=Size, BH.adj.Pvalue=BH.adj.Pvalue, Gene=enrichedGenesInSig)
  res <- res[order(res$Pvalue), ]
  res<- res[ res$BH.adj.Pvalue < adj.cutoff, ]
  subset(res, Size < 1000)
}


# perform GO enrichment analysis with 
run_GO <- function(genes, genome){
  if( genome %in% c("mm8", "mm9", "mm10")){ 
    goDB <- "org.Mm.eg"
    specie <- "Mouse"
  }
  
  if( genome %in% c("hg18", "hg19", "hg38")){ 
    goDB <- "org.Hs.eg"
    specie <- "Human"
  }
  
  adj.cutoff <- 0.05
  
  # Prepare gene list: map gene symbol to ensembl id
  annoDB <- paste(goDB, ".db", sep="")
  SYMBOL2EG <- eval(parse(text = paste(goDB, "SYMBOL2EG", sep="")))
  GO2ALLEGS <- eval(parse(text = paste(goDB, "GO2ALLEGS", sep="")))
  SYMBOL <- eval(parse(text = paste(goDB, "SYMBOL", sep="")))
  all.genes <- mappedkeys( SYMBOL2EG )
  univ <- unique(unlist(mget(all.genes, SYMBOL2EG  )))
  gene_list <- intersect(genes, all.genes)
  entrezList <- unique(unlist(mget(gene_list, SYMBOL2EG)))
  
  # start test
  ParamObjs <- new("GOHyperGParams", geneIds = entrezList, universeGeneIds = univ, annotation = annoDB, ontology = "BP", testDirection = "over")
  hyper_res <- summary( hyperGTest(ParamObjs), pvalue=1)
  adjpValueList <- p.adjust( hyper_res[,2], method="BH" )
  output <- cbind( hyper_res, adjpValueList )[ adjpValueList < adj.cutoff, ]
  if( nrow(output) > 0 ) {
    colnames( output ) <- c( colnames( hyper_res), "BH-adj.Pvalue" )
    goMaps <- lapply(output[["GOBPID"]], function(x) unlist(mget(x, GO2ALLEGS)))
    goSelected <- lapply(goMaps, function(x) { 
      temp <- entrezList[entrezList %in% x]; unlist(mget(temp, SYMBOL))} )
    res <- data.frame( output, Genes=unlist(lapply(goSelected, function(x) paste(x, collapse=";"))) )
  }
  
  res <- res[order(res$Pvalue), ]
  subset(res, Size < 1000)
}


