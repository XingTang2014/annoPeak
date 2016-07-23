
# concatenate elements of a list into a data frame 
listToDataframe <- function(Lis){
  n <- length(Lis)
  lables <- vector()
  if( length(unique(lapply(Lis, class))) != 1 ){
    stop("not the same type of data")
  }
  # when elements of the list are matrix or data.frame
  if( all(unlist(lapply(Lis, function(x){ is.data.frame(x) | is.matrix(x) }))) ){
    if(length(unique(lapply(Lis, ncol))) != 1 ){
      stop("unequal number of columns")
    }
    d <- matrix(nrow=0, ncol=ncol(Lis[[1]]))
    for(i in names(Lis)){
      lables <- c(lables, rep(i, nrow(Lis[[i]])))
      d <- rbind(d, Lis[[i]])
    }
    d <- data.frame(d, lables=lables)
  }
  # when elements of the list are vector
  if( all(unlist(lapply(Lis, is.vector))) ){
    d <- vector()
    for(i in names(Lis)){
      lables <- c(lables, rep(i, length(Lis[[i]])))
      d <- c(d, Lis[[i]])
    }
    d <- data.frame(values=d, lables=lables)
  }
  d
}


# convert entrez gene id to gene symbol
entrez2symbol <- function(entrez_ids, orgAnn){
  mapped_genes <- mappedkeys(orgAnn)
  orgAnn <- as.list(orgAnn)
  symbols <- rep(NA, length(entrez_ids))
  has_symbol <- entrez_ids %in% mapped_genes
  symbols[ has_symbol ] <- unlist(orgAnn[ entrez_ids[ has_symbol ] ])
  symbols
}