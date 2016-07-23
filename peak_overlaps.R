
generate_peak_overlaps <- function(dataList, minimumDist4overlap){
  # format the original bed file, keep chr, start, end
  # ignor lenght of peak, use the center
  # name the columns as required by GRange
  dataList <- lapply(dataList, function(x){
    x <- x[,1:3]
    colnames(x) <- c("space", "start", "end")
    mid_pos <- round((as.numeric(x[,2]) + as.numeric(x[,3])) /2)
    x[,2] <- mid_pos
    x[,3] <- mid_pos
    x
  })
  minimumDist4overlap <- as.numeric(minimumDist4overlap)
  
  overlaps_table <- matrix(0,ncol=length(dataList), nrow=length(dataList))
  peak_sets <- names(dataList)
  colnames(overlaps_table) <- peak_sets
  rownames(overlaps_table) <- peak_sets
  for(i in 1:length(peak_sets)){
    for(j in i:length(peak_sets)){
      if(i == j){
        overlaps_table[i,j] <- nrow(dataList[[ peak_sets[i] ]])
      } else {
        gr_i <- toGRanges(dataList[[ peak_sets[i] ]], format="BED")
        gr_j <- toGRanges(dataList[[ peak_sets[j] ]], format="BED")
        overlaps <- findOverlaps(gr_i, gr_j, maxgap=minimumDist4overlap)
        overlaps_table[i,j] <- length(unique(overlaps@queryHits))
        overlaps_table[j,i] <- length(unique(overlaps@subjectHits))
      }
    }
  }
  rates_table <- overlaps_table
  for(i in peak_sets){
    rates_table[i,] <- overlaps_table[i,] / overlaps_table[i,i]
  }
  rates_table <- apply(rates_table, c(1,2), function(x){
    x <- format(x * 100, nsmall=2, digits=2)
    paste(x, "%", sep="")
  })
  overlaps_table <- apply(overlaps_table, c(1,2), as.character)
  res <- list(counts=overlaps_table, rates=rates_table)
  res
}

generate_gene_overlaps <- function(targets){
  overlaps_table <- matrix(0,ncol=length(targets), nrow=length(targets))
  peak_sets <- names(targets)
  colnames(overlaps_table) <- peak_sets
  rownames(overlaps_table) <- peak_sets
  for(i in 1:length(peak_sets)){
    for(j in i:length(peak_sets)){
      if(i == j){
        overlaps_table[i,j] <- length(unique(targets[[i]]))
      } else {
        overlaps_table[i,j] <- length(unique(intersect(targets[[i]], targets[[j]])))
        overlaps_table[j,i] <- length(unique(intersect(targets[[i]], targets[[j]])))
      }
    }
  }
  rates_table <- overlaps_table
  for(i in peak_sets){
    rates_table[i,] <- overlaps_table[i,] / overlaps_table[i,i]
  }
  rates_table <- apply(rates_table, c(1,2), function(x){
    x <- format(x * 100, nsmall=2, digits=2)
    paste(x, "%", sep="")
  })
  overlaps_table <- apply(overlaps_table, c(1,2), as.character)

  all_genes <- unique(unlist(targets))  
  merged_table <- matrix(FALSE, nrow=length(all_genes), ncol=length(targets))
  rownames(merged_table) <- all_genes
  colnames(merged_table) <- names(targets)
  for(i in names(targets)){
    merged_table[ targets[[i]],i] <- TRUE 
  }
  
  res <- list(counts=overlaps_table, rates=rates_table, merged_table=merged_table)
  res
}

