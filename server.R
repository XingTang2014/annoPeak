
shinyServer(function(input, output, session) {
  
	library(RColorBrewer)
  library(GenomicRanges)
  library(VennDiagram)
  library(ggplot2)
  library(reshape2)
  library(ChIPpeakAnno)
  library(xtable)
  library(biomaRt)
  library(GSEABase)
  library(GO.db)
  library(GOstats)
  library(gplots)

  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  
	source("target_gene_prediction.R")
	source("pair_peak.R")	
  source("general.R")	
  source("peak_overlaps.R")
  source("GSEA.R")
  
  
  options(stringsAsFactors = FALSE)
  
  
  ########################################################################
  # prepare
  ########################################################################
  	# *** Read in data matrix ***
	dataL <- reactive({
	  data <- list()
		if(input$dataInput==1){
			if(input$sampleData_A==1){
				temp <- read.table("RbKO_crypts_anti-E2F3.bed", sep="\t", header=FALSE)
				data[["sampleData_A"]] <- temp[,1:min(4,ncol(temp))]
			} 
		  if(input$sampleData_B==1){
				temp <- read.table("RbKO_crypts_anti-Myc.bed", sep="\t", header=FALSE)
				data[["sampleData_B"]] <- temp[,1:min(4,ncol(temp))]
			}
		} else if(input$dataInput==2){
		  for( i in 1:10){
		    upload <- paste("upload", i, sep="")
		    peakSet <- paste("peakSet", i, sep="")
		    # Avoid error message while file is not uploaded yet
		    if (is.null(input[[upload]]))  {next}
		    #if (is.null(input[[i]]))  {return(NULL)}
		    inFile <- input[[upload]]
		    # Get the separator
		    mySep<-switch(input$fileSepDF, '1'=",",'2'="\t",'3'=";", '4'="") #list("Comma"=1,"Tab"=2,"Semicolon"=3)
		    if(file.info(inFile$datapath)$size<=10485800){
		      temp <- read.table(inFile$datapath, sep=mySep, header=FALSE)
		      data[[ input[[peakSet]] ]] <- temp[,1:min(ncol(temp), 4)]
		    } else print("File is bigger than 10MB and will not be uploaded.")
		  }
		} 
	  for(i in names(data)){
	    colnames(data[[i]]) <- c("space", "start", "end", "names")[1:ncol(data[[i]])]
	  }
		return(data)

	})
	
	# *** The plot dimensions ***
	heightSize <- reactive ({ input$myHeight })
	widthSize <- reactive ({ input$myWidth })

	## *** Data in table ***
	tabelize_input_bed <- function() {
	  tables <- list() # create a list to hold all tables
	  data <- dataL()
	  for (variable in names(data)) { # go through all possible values of variables
	    table <- data[[variable]]
	    table <- table[1:min(5, nrow(table)),]
	    tables[[as.character(variable)]] <- 
	      # save table into slot in created list
	      # print table as HTML with additional formatting options
	      print(xtable(table, caption=paste("Peak sets:", variable)),
	            type="html",
	            html.table.attributes='class="data table table-bordered table-condensed"',
	            caption.placement="top")
	  }
	  return(lapply(tables, paste)) # return HTML tables pasted together
	}
	
	output$tables <- renderUI({
	  out <- unlist(tabelize_input_bed()) 
	  # additional options to make rendering possible
	  return(div(HTML(out),class="shiny-html-output"))
	})
	
	
	# *** Target gene prediction ***
	geneTargetList <- reactive({
	  geneTargets <- list()
	  data <- dataL()
	  promoter.left <- as.numeric(input$promoter.left)
	  promoter.right <- as.numeric(input$promoter.right)
	  enhancer.left <- as.numeric(input$enhancer.left)
	  for(i in names(data)){
	    geneTargets[[i]] <- target_gene_prediction(data[[i]], input$genomeBuild,
	                                               promoter.left, promoter.right, enhancer.left)
	  }
	  geneTargets
	})
	output$genetargets <- renderTable({
	  targets <- geneTargetList()
	  return(targets[[1]])
	})
	
	
	########################################################################
	# peak length distribution
	########################################################################
	generatePeakLen <- function(){
	  par(mar=c(5,8,8,2)) # c(bottom, left, top, right)
	  
	  myColours<-gsub("\\s","", strsplit(input$myColours,",")[[1]])
	  myColours<-gsub("0x","#", myColours)
	  
	  
	  if(as.numeric(input$subregion) != 0){
	    data <- select_subregion_peaks()
	  } else {
	    data <- dataL()  
	  }
	  peakLen <- vector()
	  sampleLab <- vector()
	  for(i in names(data)){
	    d <- data[[i]]
	    sampleLab <- c(sampleLab, rep(i, nrow(d)))
	    peakLen <- c(peakLen, as.numeric(d[,3]) - as.numeric(d[,2]))
	  }

	  toPlot <- data.frame(length=peakLen, peakSet=sampleLab)
	  if(input$changePeakSetsLabel){
	    myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	  } else {
	    myLabels <- levels(as.factor(toPlot$peakSet))
	  }
	 
	  print(ggplot(toPlot, aes(x=length, fill=peakSet)) + 
	          geom_histogram(position = "dodge") +
	          scale_color_manual(values=myColours, labels=myLabels) +
	          labs(y="Peak Count", x="Peak Length") +
	          theme_bw(base_size = 18 * input$fontSizes, base_family=input$Font) + 
	          theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
	                panel.grid.minor = element_blank(), 
	                axis.line = element_line(colour = "black")))
	}
	# *** Histgram (using 'generatePeakLen'-function) ***
	output$peakLengthDensity <- renderPlot({
	  generatePeakLen()
	  	}, height = heightSize, width = widthSize)
	
	## *** Download EPS file ***
	output$downloadPlotEPS_length <- downloadHandler(
	  filename <- function() { paste('PeakLenDistribution.eps') },
	  content <- function(file) {
	    postscript(file, horizontal = FALSE, onefile = FALSE, paper = "special", width = input$myWidth/72, height = input$myHeight/72)
	    ## ---------------
	    generatePeakLen()
	    ## ---------------
	    dev.off()
	  },
	  contentType = 'application/postscript'
	)
	## *** Download PDF file ***
	output$downloadPlotPDF_length <- downloadHandler(
	  filename <- function() { paste('PeakLenDistribution.pdf') },
	  content <- function(file) {
	    pdf(file, width = input$myWidth/72, height = input$myHeight/72)
	    ## ---------------
	    generatePeakLen()
	    ## ---------------
	    dev.off()
	  },
	  contentType = 'application/pdf' # MIME type of the image
	)
	## *** Download SVG file ***
	output$downloadPlotSVG_length <- downloadHandler(
	  filename <- function() { paste('PeakLenDistribution.svg') },
	  content <- function(file) {
	    svg(file, width = input$myWidth/72, height = input$myHeight/72)
	    ## ---------------
	    generatePeakLen()
	    ## ---------------
	    dev.off()
	  },
	  contentType = 'image/svg'
	)
	
	######################################
	# select subregions
	select_subregion_peaks <- reactive({
	  regionCode <- c("enhancer", "promoter", "inside", "upstream", "downstream")
	  subregion <- regionCode[as.numeric(input$subregion)]
	  #"enhancer"=2, "promoter"=3, "genebody"=4,  "upstream"=5, "downstream"=6
	  targets <- geneTargetList()
	  lapply(targets, function(x){
	    subset(x, gene_region == subregion)[,c("peak_chr", "peak_start", "peak_end")]
	  })
	})
	
	########################################################################
	# peak region wise distribution
	########################################################################
		# *** Generate the histogram ***
	generateHist <- function(){
	  par(mar=c(5,8,8,2)) # c(bottom, left, top, right)
	  
	  myColours<-gsub("\\s","", strsplit(input$myColours,",")[[1]])
	  myColours<-gsub("0x","#", myColours)
	  
	  toPlot <- peak_distribution(geneTargetList())
  
	  if(input$changePeakSetsLabel){
	    myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	  } else {
	    myLabels <- levels(toPlot$Peakset)
	  }
	  
	  
	  if(input$histYaxis == 1){
	    toPlot[,"yaxis"] <- as.numeric(sub("%", "", toPlot[,"Percentage"])) 
	    print(ggplot(toPlot, aes(Region, yaxis, fill=Peakset)) + 
	            geom_bar(stat="identity", position = position_dodge(width = .8), width = 0.7) + 
	            geom_text(aes(label = PeakCounts), position = position_dodge(width = .8), vjust = -0.5) +
	            geom_text(aes(label = Percentage), position = position_dodge(width = .8), vjust = -2, colour="red") +
	            ylim(c(0, max(toPlot$yaxis) * 1.1)) + 
	            scale_fill_manual(values=myColours, labels=myLabels) +
	            ylab("Peak proportion (%)") +
	            theme_bw(base_size = 18 * input$fontSizes, base_family=input$Font) + 
	            theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
	                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
	  }
	  if(input$histYaxis == 2){
	    print(ggplot(toPlot, aes(Region, PeakCounts, fill=Peakset)) + 
	            geom_bar(stat="identity", position = position_dodge(width = .8), width = 0.7) + 
	            geom_text(aes(label = PeakCounts), position = position_dodge(width = .8), vjust = -0.5) +
	            geom_text(aes(label = Percentage), position = position_dodge(width = .8), vjust = -2, colour="red") +
	            ylim(c(0, max(toPlot$PeakCounts) * 1.1)) + 
	            ylab("Peak Counts") +
	            scale_fill_manual(values=myColours, labels=myLabels) +
	            theme_bw(base_size = 18 * input$fontSizes, base_family=input$Font) + 
	            theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
	                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
	  }
	  
	}
	
	
	# *** Output peak region wise distribution statistics in table below plot ***
	output$histStatsTable <- renderTable({
	  
	  # Create a Progress object
	  progress <- shiny::Progress$new(session)
	  progress$set(message = "Predicting target genes", value = 0)
	  # Close the progress when this reactive exits (even if there's an error)
	  on.exit(progress$close())
	  
	  targets <- geneTargetList()
	  if(input$changePeakSetsLabel){
	    myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	    names(targets) <- myLabels
	  }
	  toPlot <- peak_distribution(targets)
	  return(toPlot)
	  
	  })
	

	# *** Histgram (using 'generateHist'-function) ***
	output$hist <- renderPlot({
	  generateHist()
	}, height = heightSize, width = widthSize)
	
	## *** Download EPS file ***
	output$downloadPlotEPS_hist <- downloadHandler(
	  filename <- function() { paste('PeakDistribution.eps') },
	  content <- function(file) {
	    postscript(file, horizontal = FALSE, onefile = FALSE, paper = "special", width = input$myWidth/72, height = input$myHeight/72)
	    ## ---------------
	    generateHist()
	    ## ---------------
	    dev.off()
	  },
	  contentType = 'application/postscript'
	)
	## *** Download PDF file ***
	output$downloadPlotPDF_hist <- downloadHandler(
	  filename <- function() { paste('PeakDistribution.pdf') },
	  content <- function(file) {
	    pdf(file, width = input$myWidth/72, height = input$myHeight/72)
	    ## ---------------
	    generateHist()
	    ## ---------------
	    dev.off()
	  },
	  contentType = 'application/pdf' # MIME type of the image
	)
	## *** Download SVG file ***
	output$downloadPlotSVG_hist <- downloadHandler(
	  filename <- function() { paste('PeakDistribution.svg') },
	  content <- function(file) {
	    svg(file, width = input$myWidth/72, height = input$myHeight/72)
	    ## ---------------
	    generateHist()
	    ## ---------------
	    dev.off()
	  },
	  contentType = 'image/svg'
	)
	
	# *** Download peak region wise distribution data in csv format ***
	output$downloadHistData <- downloadHandler(
	  filename = function() { "PeakDistributionData.csv" },
	  content = function(file) {
	    write.csv(peak_distribution(geneTargetList()), file, row.names=FALSE)
	  }) ###
	# *** Download target genes data in csv format ***
	output$downloadTargetGenes <- downloadHandler(
	  filename = function() { "targetGenesData.csv" },
	  content = function(file) {
	    write.csv(listToDataframe(geneTargetList()), file, row.names=TRUE)
	  }) ###

	
	########################################################################
	# pair peak
	########################################################################
		# *** Select peak set to do pair peak distances plot ***
	peakSets <- reactive({
	  
	  data <- dataL()  
	  variable <- names(data)
	  if(input$changePeakSetsLabel){
	    variable <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	  }
	  peakSetsL <- list()
	  for(i in variable){
	    peakSetsL[[i]] <- i
	  }
	  print(peakSetsL)
	  return(peakSetsL)
	})
	output$pairPeakSet1 <- renderUI({
	  selectInput(inputId = "pairPeakSet1", label = h5("Peak set 1:"), peakSets(), selected=unlist(peakSets()[1]))
	})
	output$pairPeakSet2 <- renderUI({
	  selectInput(inputId = "pairPeakSet2", label = h5("Peak set 2:"), peakSets(),selected=unlist(peakSets()[2]))
	})

	# *** functions to plot pair peak distances distribution ***
	plotPairPeakPlot <- function(){
	  # Create a Progress object
	  progress <- shiny::Progress$new(session)
	  progress$set(message = "Calculating closest peak distances", value = 0)
	  # Close the progress when this reactive exits (even if there's an error)
	  on.exit(progress$close())
	  
	  if(as.numeric(input$subregion) != 0){
	    data <- select_subregion_peaks()
	    print(lapply(data, nrow))
	  } else {
	    data <- dataL()  
	  }
	  
	  if(input$changePeakSetsLabel){
	    names(data) <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	  }
	  distances_list <- pair_peak_distance(data[[input$pairPeakSet1]],
	                                       data[[input$pairPeakSet2]],
	                                       input$pairPeakSet1, input$pairPeakSet2)
	  
	  if( class(distances_list) != "character"){
	    par(mar=c(5,8,8,2)) # c(bottom, left, top, right)
	    
	    myColours<-gsub("\\s","", strsplit(input$myColours,",")[[1]])
	    myColours<-gsub("0x","#", myColours)
	    
	    toPlot <- matrix(ncol=2,nrow=0)
	    for(i in names(distances_list)){
	      toPlot <- rbind(toPlot, data.frame(distance=distances_list[[i]], variable=rep(i, length(distances_list[[i]]))))
	    }
	    
	    if(input$changePeakSetsLabel){
	      myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	    } else {
	      myLabels <- levels(as.factor(toPlot$variable))
	    }
	    
	    print(ggplot(toPlot, aes(distance, color=variable)) + 
	            geom_density(size=0.65) + 
	            scale_x_log10() + 
	            scale_color_manual(values=myColours, labels=myLabels) + 
	            theme_bw(base_size = 18 * input$fontSizes, base_family=input$Font) + 
	            theme(panel.border = element_blank(), 
	                  panel.grid.major = element_blank(), 
	                  panel.grid.minor = element_blank(), 
	                  axis.line = element_line(colour = "black")))
	  }
	}
	# *** pair peak plot ***
	output$pairPeakPlot <- renderPlot({
	  plotPairPeakPlot()
	}, height = heightSize, width = widthSize)
	
	## *** Download EPS file ***
	output$downloadPlotEPS_pair <- downloadHandler(
	  filename <- function() { paste('nearestPeakDistances.eps') },
	  content <- function(file) {
	    postscript(file, horizontal = FALSE, onefile = FALSE, paper = "special", width = input$myWidth/72, height = input$myHeight/72)
	    plotPairPeakPlot()
	    dev.off()
	  },
	  contentType = 'application/postscript'
	)
	## *** Download PDF file ***
	output$downloadPlotPDF_pair <- downloadHandler(
	  filename <- function() { paste('nearestPeakDistances.pdf') },
	  content <- function(file) {
	    pdf(file, width = input$myWidth/72, height = input$myHeight/72)
	    plotPairPeakPlot()
	    dev.off()
	  },
	  contentType = 'application/pdf' # MIME type of the image
	)
	## *** Download SVG file ***
	output$downloadPlotSVG_pair <- downloadHandler(
	  filename <- function() { paste('nearestPeakDistances.svg') },
	  content <- function(file) {
	    svg(file, width = input$myWidth/72, height = input$myHeight/72)
	    plotPairPeakPlot()
	    dev.off()
	  },
	  contentType = 'image/svg'
	)
	
	
	########################################################################
	# peak overlaps
	########################################################################
	peakOverlaps <- reactive({
	  
	  if(as.numeric(input$subregion) != 0){
	    data <- select_subregion_peaks()
	  } else {
	    data <- dataL()
	  }
	  res <- generate_peak_overlaps(data, input$minimumDist4overlap)
	  if(input$changePeakSetsLabel){
	    myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	    res <- lapply(res, function(d){colnames(d) <- myLabels; rownames(d) <- myLabels; d })
	  } 
	  return(res)
	})
	output$peakOverlapsTable <- renderTable({
	  peakOverlaps()$counts
	})
	output$peakOverlapsRateTable <- renderTable({
	  peakOverlaps()$rates
	})
	output$downloadPeakOverlapsTable <- downloadHandler(
	  filename = function() { "PeakOverlapsCountsTable.csv" },
	  content = function(file) {
	    write.csv(peakOverlaps()$counts, file, row.names=TRUE)
	  }) ###
	output$downloadPeakOverlapsRatesTable <- downloadHandler(
	  filename = function() { "PeakOverlapsRatesTable.csv" },
	  content = function(file) {
	    write.csv(peakOverlaps()$rates, file, row.names=TRUE)
	  }) ###
	
	
	########################################################################
	# gene overlaps
	########################################################################
	select_subregion_targets <- reactive({
	  regionCode <- c("enhancer", "promoter", "inside", "upstream", "downstream")
	  #"enhancer"=2, "promoter"=3, "genebody"=4,  "upstream"=5, "downstream"=6
	  subregion <- regionCode[as.numeric(input$subregion)]
	  targets <- geneTargetList()
	  targets <- lapply(targets, function(x){
	    subset(x, gene_region == subregion)
	  })
	})
	plotGeneVenn <- function(){
	  if(as.numeric(input$subregion) != 0){
	    targets <- select_subregion_targets()
	  } else {
	    targets <- geneTargetList()
	  }
	  
	  targets <- lapply(targets, function(x){
	    setdiff(as.character(unique(x[,"nearest_gene_symbol"])), NA)
	  })
	  
	  myColours<-gsub("\\s","", strsplit(input$myColours,",")[[1]])
	  myColours<-gsub("0x","#", myColours)
	  
	  if(input$changePeakSetsLabel){
	    myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	    names(targets) <- myLabels
	  } 
	  if(length(targets) <= 5){
	    grid.draw(venn.diagram(targets, filename=NULL, 
	                           col=myColours[1:length(targets)],
	                           cat.fontfamily=input$Font,
	                           fontfamily=input$Font,
	                           cex=1.2 * input$fontSizes, 
	                           cat.cex=1.2 * input$fontSizes,
	                           # c(bottom, left, top, right)
	                           margin=c(0.1,0.1,0.1,0.1)), 
	              recording=FALSE)   
	  } else {
	    textplot("Too many sets of peaks to plot venn diagram.
	             Maximum count of acceptable sets of peaks is 5.", col.data="red")
	  }
	 
	}
	# generate gene venn plot output
	output$geneVennPlot <- renderPlot({
	  plotGeneVenn()
	}, height = heightSize, width = widthSize)
	## *** Download EPS file ***
	output$downloadPlotEPS_geneVenn <- downloadHandler(
	  filename <- function() { paste('targetGeneVenn.eps') },
	  content <- function(file) {
	    postscript(file, horizontal = FALSE, onefile = FALSE, paper = "special", width = input$myWidth/72, height = input$myHeight/72)
	    plotGeneVenn()
	    dev.off()
	  },
	  contentType = 'application/postscript'
	)
	output$downloadGeneMergedTable <- downloadHandler(
	  filename = function() { "GeneMergedTable.csv" },
	  content = function(file) {
	    write.csv(geneOverlaps()$merged_table, file, row.names=TRUE)
	  }) ###
	## *** Download PDF file ***
	output$downloadPlotPDF_geneVenn <- downloadHandler(
	  filename <- function() { paste('targetGeneVenn.pdf') },
	  content <- function(file) {
	    pdf(file, width = input$myWidth/72, height = input$myHeight/72)
	    plotGeneVenn()
	    dev.off()
	  },
	  contentType = 'application/pdf' # MIME type of the image
	)
	## *** Download SVG file ***
	output$downloadPlotSVG_geneVenn <- downloadHandler(
	  filename <- function() { paste('targetGeneVenn.svg') },
	  content <- function(file) {
	    svg(file, width = input$myWidth/72, height = input$myHeight/72)
	    plotGeneVenn()
	    dev.off()
	  },
	  contentType = 'image/svg'
	)
	
	geneOverlaps <- reactive({
	  if(as.numeric(input$subregion) != 0){
	    targets <- select_subregion_targets()
	  } else {
	    targets <- geneTargetList()
	  }
	  
	  targets <- lapply(targets, function(x){
	    setdiff(as.character(unique(x[,"nearest_gene_symbol"])), NA)
	  })
	  res <- generate_gene_overlaps(targets)
	  if(input$changePeakSetsLabel){
	    myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	    res <- lapply(res, function(d){colnames(d) <- myLabels; rownames(d) <- myLabels; d })
	  } 
	  return(res)
	})
	
	output$geneOverlapsTable <- renderTable({
	  geneOverlaps()$counts
	})
	output$geneOverlapsRateTable <- renderTable({
	  geneOverlaps()$rates
	})
	output$downloadGeneOverlapsTable <- downloadHandler(
	  filename = function() { "GeneOverlapsCountsTable.csv" },
	  content = function(file) {
	    write.csv(geneOverlaps()$counts, file, row.names=TRUE)
	  }) ###
	output$downloadGeneOverlapsRatesTable <- downloadHandler(
	  filename = function() { "GeneOverlapsRatesTable.csv" },
	  content = function(file) {
	    write.csv(geneOverlaps()$rates, file, row.names=TRUE)
	  }) ###
	
	########################################################################
	# perform functional enrichment
	########################################################################
	geneSetsL <- reactive({
	  data <- dataL()
	  variable <- names(data)
	  if(input$changePeakSetsLabel){
	    variable <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	  }
	  geneSets <- list()
	  geneSets[["None"]] <- "None"
	  for(i in variable){
	    geneSets[[i]] <- i
	  }
	  print(geneSets)
	  return(geneSets)
	})
	output$GeneSet1 <- renderUI({
	  selectInput(inputId = "GeneSet1", label = h5("Gene set 1:"), geneSetsL(), selected=unlist(geneSetsL()[2]))
	})
	output$GeneSet2 <- renderUI({
	  selectInput(inputId = "GeneSet2", label = h5("Gene set 2:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	output$GeneSet3 <- renderUI({
	  selectInput(inputId = "GeneSet3", label = h5("Gene set 3:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	output$GeneSet4 <- renderUI({
	  selectInput(inputId = "GeneSet4", label = h5("Gene set 4:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	output$GeneSet5 <- renderUI({
	  selectInput(inputId = "GeneSet5", label = h5("Gene set 5:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	output$GeneSet6 <- renderUI({
	  selectInput(inputId = "GeneSet6", label = h5("Gene set 6:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	output$GeneSet7 <- renderUI({
	  selectInput(inputId = "GeneSet7", label = h5("Gene set 7:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	output$GeneSet8 <- renderUI({
	  selectInput(inputId = "GeneSet8", label = h5("Gene set 8:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	output$GeneSet9 <- renderUI({
	  selectInput(inputId = "GeneSet9", label = h5("Gene set 9:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	output$GeneSet10 <- renderUI({
	  selectInput(inputId = "GeneSet10", label = h5("Gene set 10:"), geneSetsL(),selected=unlist(geneSetsL()[1]))
	})
	

	# output GSEA results
	get_GSEA_res <- eventReactive(input$goButton, {
	  
	  # Create a Progress object
	  progress <- shiny::Progress$new(session)
	  progress$set(message = "Predicting enriched gene signature. This may take a while", value = 0)
	  # Close the progress when this reactive exits (even if there's an error)
	  on.exit(progress$close())
	  
	  if(as.numeric(input$subregion) != 0){
	    targets <- select_subregion_targets()
	  } else {
	    targets <- geneTargetList()
	  }
	  
	  if(input$changePeakSetsLabel){
	    myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	    names(targets) <- myLabels
	  }
	  samples <- setdiff(unique(c(input$GeneSet1, input$GeneSet2, input$GeneSet3, input$GeneSet4, input$GeneSet5, input$GeneSet6, input$GeneSet7, input$GeneSet8, input$GeneSet9, input$GeneSet10)), "None")
	  
	  if(length(samples)==0){ return(data.frame())}
	  
	  if(input$modeOfGeneSetsSelection == 1){
	    # intersection
	    genes <- unique(targets[[ samples[1] ]][, "nearest_gene_symbol"] )
	    if( length(samples) > 1){
	      for(i in 2:length(samples)){
	        genes <- intersect(genes, unique(targets[[ samples[i] ]][, "nearest_gene_symbol"] ))
	      }  
	    }
	  }
	  if(input$modeOfGeneSetsSelection == 2){
	    # union
	    genes <- unique(targets[[ samples[1] ]][, "nearest_gene_symbol"] )
	    if( length(samples) > 1){
	      for(i in 2:length(samples)){
	        genes <- union(genes, unique(targets[[ samples[i] ]][, "nearest_gene_symbol"] ))
	      }  
	    }
	  }
	  if(input$modeOfGeneSetsSelection == 3){
	    # subtract
	    genes <- unique(targets[[ samples[1] ]][, "nearest_gene_symbol"] )
	    if( length(samples) > 1){
	      for(i in 2:length(samples)){
	        genes <- setdiff(genes, unique(targets[[ samples[i] ]][, "nearest_gene_symbol"] ))
	      }  
	    }
	  }
	    genes <- setdiff(genes, NA)
	  print(genes)
	    res <- run_GSEA(genes, input$genomeBuild)
	    return(res)
	})
	
	output$GSEATable <- renderTable({
	  toPlot <- get_GSEA_res()
	  if(nrow(toPlot) == 0){
	    return(toPlot)
	  } else {
	    toPlot <- toPlot[1:min(10, nrow(toPlot)), 1:8]
	    rownames(toPlot) <- NULL
	    return(toPlot)
	  }
	})
	output$downloadEnrichedSig <- downloadHandler(
	  filename = function() { "EnrichedSignatures.csv" },
	  content = function(file) {
	    write.csv(get_GSEA_res(), file, row.names=FALSE)
	  }) ###
	
	# output GO results
	get_GO_res <- eventReactive(input$goButton,{
	  
	  # Create a Progress object
	  progress <- shiny::Progress$new(session)
	  progress$set(message = "Predicting enriched GO terms. This may take a while", value = 0)
	  # Close the progress when this reactive exits (even if there's an error)
	  on.exit(progress$close())
	  
	  if(as.numeric(input$subregion) != 0){
	    targets <- select_subregion_targets()
	  } else {
	    targets <- geneTargetList()
	  }
	  
	  if(input$changePeakSetsLabel){
	    myLabels <- gsub("^\\s","", strsplit(input$myLabels,",")[[1]])
	    names(targets) <- myLabels
	  }
	  samples <- setdiff(unique(c(input$GeneSet1, input$GeneSet2, input$GeneSet3, input$GeneSet4, input$GeneSet5, input$GeneSet6, input$GeneSet7, input$GeneSet8, input$GeneSet9, input$GeneSet10)), "None")
	  
	  if(length(samples)==0){ return(data.frame())}
	  
	  if(input$modeOfGeneSetsSelection == 1){
	    # intersection
	    genes <- unique(targets[[ samples[1] ]][, "nearest_gene_symbol"] )
	    if( length(samples) > 1){
	      for(i in 2:length(samples)){
	        genes <- intersect(genes, unique(targets[[ samples[i] ]][, "nearest_gene_symbol"] ))
	      }  
	    }
	  }
	  if(input$modeOfGeneSetsSelection == 2){
	    # union
	    genes <- unique(targets[[ samples[1] ]][, "nearest_gene_symbol"] )
	    if( length(samples) > 1){
	      for(i in 2:length(samples)){
	        genes <- union(genes, unique(targets[[ samples[i] ]][, "nearest_gene_symbol"] ))
	      }  
	    }
	  }
	  if(input$modeOfGeneSetsSelection == 3){
	    # subtract
	    genes <- unique(targets[[ samples[1] ]][, "nearest_gene_symbol"] )
	    if( length(samples) > 1){
	      for(i in 2:length(samples)){
	        genes <- setdiff(genes, unique(targets[[ samples[i] ]][, "nearest_gene_symbol"] ))
	      }  
	    }
	  }
	  
	  genes <- setdiff(genes, NA)
	  res <- run_GO(genes, input$genomeBuild)
	  return(res)
	})
	
	output$GOTable <- renderTable({
	  toPlot <- get_GO_res()
	  if(nrow(toPlot) == 0){
	    return(toPlot)
	  } else {
	    return(toPlot[1:min(10, nrow(toPlot)), 1:8])
	  }
	})
	output$downloadEnrichedGO <- downloadHandler(
	  filename = function() { "EnrichedGO.csv" },
	  content = function(file) {
	    write.csv(get_GO_res(), file, row.names=FALSE)
	  }) ###
	

})


# *** Determine extent of range ***
#myRange <- reactive({
#  if(input$whiskerType==0){myRange<-c(-1.5)} 
#  else if(input$whiskerType==1){myRange<-c(0)} 
#  else if (input$whiskerType==2){myRange<-c(5)}
#  return(myRange)
#})


