
# output: peak_id, tss distances, nearest gene, gene region
target_gene_prediction <- function(inputBed, genome, promoter.left, promoter.right, enhancer.left){
  if(genome=="mm9"){
    gregions <- genes(TxDb.Mmusculus.UCSC.mm9.knownGene)
    orgAnn <- org.Mm.egSYMBOL
  }
  if(genome=="mm10"){
    gregions <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
    orgAnn <- org.Mm.egSYMBOL
  }
  if(genome=="hg19"){
    gregions <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    orgAnn <- org.Hs.egSYMBOL
  }
  if(genome=="hg18"){
    gregions <- genes(TxDb.Hsapiens.UCSC.hg18.knownGene)
    orgAnn <- org.Hs.egSYMBOL
  }
  if(genome=="hg38"){
    gregions <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
    orgAnn <- org.Hs.egSYMBOL
  }
  bed.anno <- annotatePeakInBatch(toGRanges(inputBed, format="BED"), 
                                  AnnotationData=gregions, 
                                  output="nearestLocation", 
                                  PeakLocForDistance="middle",
                                  FeatureLocForDistance="TSS",
                                  select="first")
  yy <-data.frame(insideFeature=bed.anno$insideFeature, tss_distances=as.numeric(bed.anno$distancetoFeature)) 
  gene_region <- apply(yy, 1, function(x){
                         if( any(is.na(x))) {x[1]} 
                         else if (as.numeric(x[2]) >= -enhancer.left & as.numeric(x[2]) < -promoter.left) {"enhancer"} 
                         else if( as.numeric(x[2]) >= -promoter.left & as.numeric(x[2]) < promoter.right) { "promoter"} 
                         else {x[1]} 
                         })
  d <- data.frame(peak_id=bed.anno$peak, 
                  peak_chr=bed.anno@seqnames,
                  peak_start=bed.anno@ranges@start,
                  peak_end=bed.anno@ranges@start + bed.anno@ranges@width,
             tss_distances=bed.anno$distancetoFeature, 
             nearest_gene=bed.anno$feature,
             nearest_gene_symbol=entrez2symbol(bed.anno$feature, orgAnn),
             gene_region=gene_region)
  d[!is.na(d$gene_region),]
} 


peak_distribution <- function(geneTargetList){
  regions <- c("upstream", "enhancer", "promoter", "inside", "downstream")
  d <- as.data.frame(lapply(geneTargetList, function(x){
    as.vector(table(x[,"gene_region"])[regions])}))
  colnames(d) <- names(geneTargetList)
  total <- apply(d, 2, sum)
  d[["region"]] <- factor(regions, levels=regions, ordered=TRUE)
  toPlot <- melt(d)
  toPlot[["rate"]] <- paste(round(100 * apply(toPlot, 1, function(x){ as.numeric(x["value"]) / as.numeric( total[x["variable"]])}) , 1), "%", sep="")
  colnames(toPlot) <- c("Region", "Peakset", "PeakCounts", "Percentage")
  toPlot
}

