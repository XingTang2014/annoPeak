
pair_peak_distance <- function(peaks1, peaks2, name1, name2){
  if(name1 == name2){ 
    return("identical peakset")
    } # identical peakset
  
  merged_peaks <- data.frame(
    chromosome=c(as.character(peaks1[,1,drop=TRUE]), as.character(peaks2[,1,drop=TRUE])),
    pos=c(round( (as.numeric(peaks1[,2,drop=TRUE]) + as.numeric(peaks1[,3, drop=TRUE])) / 2 ),
          round( (as.numeric(peaks2[,2,drop=TRUE]) + as.numeric(peaks2[,3, drop=TRUE])) / 2 )),
    group=c(rep(name1, nrow(peaks1)), rep(name2, nrow(peaks2))))

  distances_list <- list()
  distances_list[[name1]] <- vector()
  distances_list[[name2]] <- vector()
  
  for(chr in unique(merged_peaks[,"chromosome"])){
    merged_peaks_in_chr <- subset(merged_peaks, chromosome==chr)
    merged_peaks_in_chr <- merged_peaks_in_chr[order( as.numeric(merged_peaks_in_chr[,"pos"]) ),]
    
    for(i in 1:nrow(merged_peaks_in_chr)){
      if(i ==1){ 
        before_group <- merged_peaks_in_chr[i,"group"]
        nClust <- 1
        clusters <- nClust} 
      else {
        if( before_group != merged_peaks_in_chr[i,"group"] ){
          nClust <- nClust + 1
          before_group <- merged_peaks_in_chr[i,"group"]
        }
          clusters <- c(clusters, nClust)
      }
    }
    
    if( nClust == 1 ){
      distances_list[[ merged_peaks_in_chr[1,"group"] ]] <- c(distances_list[[ merged_peaks_in_chr[1,"group"] ]],
                                                              rep(Inf, nrow(merged_peaks_in_chr)) )
      next;
    } # if all the peaks come from the same peak set
    if(nClust == 2){
      h <- merged_peaks_in_chr[1,"group"] # the first cluster
      t <- merged_peaks_in_chr[ nrow(merged_peaks_in_chr),"group"] # the last cluster
      distances_list[[ h ]]<- c(distances_list[[ h ]], 
                                abs( merged_peaks_in_chr[ merged_peaks_in_chr[,"group"]==h,"pos"] - 
                                       merged_peaks_in_chr[ which(merged_peaks_in_chr[,"group"]==t)[1],"pos"] ) )
      distances_list[[ t ]]<- c(distances_list[[ t ]], 
                                abs( merged_peaks_in_chr[ merged_peaks_in_chr[,"group"]==t,"pos"] - 
                                       merged_peaks_in_chr[ rev(which(merged_peaks_in_chr[,"group"]==h))[1],"pos"] ) )
      next;
    }
    
    
    up <- unlist(tapply(seq(1, nrow(merged_peaks_in_chr)), clusters, 
                 function(x){ rep(x[1] - 1, length(x)) }) ) # up stream cloest peak
    down <- unlist( tapply(seq(1, nrow(merged_peaks_in_chr)), clusters, 
                 function(x){ rep(x[length(x)] + 1, length(x)) }) ) # down stream cloest peak
    
    for(i in 1:nrow(merged_peaks_in_chr)){
      if(i == 1){
        distance <- abs( merged_peaks_in_chr[i,"pos"] - merged_peaks_in_chr[down[i],"pos"] )
        next
      } 
      if(i == nrow(merged_peaks_in_chr)){
        distance <- abs( merged_peaks_in_chr[i,"pos"] - merged_peaks_in_chr[up[i],"pos"] )
        next
      }
      distance <- min(abs( merged_peaks_in_chr[i,"pos"] - merged_peaks_in_chr[down[i],"pos"] ),
                      abs( merged_peaks_in_chr[i,"pos"] - merged_peaks_in_chr[up[i],"pos"] ))
      
      g <- merged_peaks_in_chr[i,"group"]
      distances_list[[g]] <- c(distances_list[[g]], distance )
    }
  }
  
  distances_list
}