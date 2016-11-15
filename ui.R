

shinyUI(pageWithSidebar(


  headerPanel("annoPeakR: a web-tool to annotate, visualize and compare peak sets from ChIP-seq/ChIP-exo",
              tags$head(tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
                        tags$style(type="text/css", "select { max-width: 200px; }"),
                        tags$style(type="text/css", "textarea { max-width: 185px; }"),
                        tags$style(type="text/css", ".jslider { max-width: 200px; }"),
                        tags$style(type='text/css', ".well { max-width: 330px; }"),
                        #tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                        tags$style(type='text/css', ".span4 { max-width: 330px; }")
                        ) 
  ),
  
  sidebarPanel(
    conditionalPanel(condition="input.tabs1=='About'",
                     h4("Introduction")
    ),
    conditionalPanel(condition="input.tabs1=='Data upload'",
                     radioButtons("genomeBuild", "Select Genome Build", list("hg18"="hg18","hg19"="hg19","hg38"="hg38", "mm9"="mm9", "mm10"="mm10"), "mm9"),
                     h4("Enter data"),
                     radioButtons("dataInput", "", list("Load sample data"=1,"Upload file"=2)),
                     conditionalPanel(condition="input.dataInput=='1'",
                                      h5("Load sample data:"),
                                      checkboxInput("sampleData_A", "Load sample data A", TRUE),
                                      checkboxInput("sampleData_B", "Load sample data B", TRUE)
                     ),
                     conditionalPanel(condition="input.dataInput=='2'",
                                      h5("Upload delimited text file: "),
                                      fileInput("upload1", "", multiple = FALSE),
                                      textInput("peakSet1", "Name of Peak Set 1", "upload1"),
                                      fileInput("upload2", "", multiple = FALSE),
                                      textInput("peakSet2", "Name of Peak Set 2", "upload2"),
                                      fileInput("upload3", "", multiple = FALSE),
                                      textInput("peakSet3", "Name of Peak Set 3", "upload3"),
                                      fileInput("upload4", "", multiple = FALSE),
                                      textInput("peakSet4", "Name of Peak Set 4", "upload4"),
                                      fileInput("upload5", "", multiple = FALSE),
                                      textInput("peakSet5", "Name of Peak Set 5", "upload5"),
                                      fileInput("upload6", "", multiple = FALSE),
                                      textInput("peakSet6", "Name of Peak Set 6", "upload6"),
                                      fileInput("upload7", "", multiple = FALSE),
                                      textInput("peakSet7", "Name of Peak Set 7", "upload7"),
                                      fileInput("upload8", "", multiple = FALSE),
                                      textInput("peakSet8", "Name of Peak Set 8", "upload8"),
                                      fileInput("upload9", "", multiple = FALSE),
                                      textInput("peakSet9", "Name of Peak Set 9", "upload9"),
                                      fileInput("upload10", "", multiple = FALSE),
                                      textInput("peakSet10", "Name of Peak Set 10", "upload10"),
                                      radioButtons("fileSepDF", "Delimiter:", list("Tab"=2,"Comma"=1,"Semicolon"=3)),#, "Space"=4))
                                      HTML('<p>Data in <a href="http://en.wikipedia.org/wiki/Delimiter-separated_values">delimited text files </a> can be separated by comma, tab or semicolon. 
                                           For example, Excel data can be exported in .csv (comma separated) or .tab (tab separated) format. </p>')
                                      )
                     ),
    conditionalPanel(condition="input.tabs1=='Data visualization'",
                     
                     h4("Select analysis module"),
                     radioButtons("plotType", "", list("Peak associated gene structures"=0, "Peak size distribition"=4, "Peak to nearest peak distances"=1, "Overlapping peak identifications"=2, "Overlapping peak-associated genes"=3)),
                     conditionalPanel(condition="input.plotType==0",
                                      h4("Change definition of promoter and enhancer"),
                                      checkboxInput("changePromoterDefinition", "Change definition of promoter region", FALSE),
                                      conditionalPanel(condition="input.changePromoterDefinition",
                                                       textInput("promoter.left", "Up stream of TSS (bp)", 5000),
                                                       textInput("promoter.right", "Down stream of TSS (bp)", 2000)
                                      ),
                                      checkboxInput("changeEnhancerDefinition", "Change definition of enhancer region", FALSE),
                                      conditionalPanel(condition="input.changeEnhancerDefinition",
                                                       textInput("enhancer.left", "Up stream of TSS (bp)", 50000),
                                                       HTML('<p>The downstream bound of enhancer region is the upstream bound of promoter region</p>')
                                      ),
                                      radioButtons("histYaxis", "Scale Y axis by", list( "Peak proportion"=1, "Peak count"=2), selected=1)
                                      ),
                     conditionalPanel(condition="input.plotType!=0",
                                      radioButtons("subregion", label=h5("Subregion selection"),
                                                   choices=list("all"=0, "enhancer"=1, "promoter"=2, "genebody"=3, "upstream"=4, "downstream"=5), selected=0)),
                     conditionalPanel(condition="input.plotType==1",
                                      h4("Select two peak sets to compare"),
                                      uiOutput("pairPeakSet1"),
                                      uiOutput("pairPeakSet2")),
                     conditionalPanel(condition="input.plotType==2",
                                      h4("Maximum allowed distances between overlapped peaks (bp)"),
                                      textInput("minimumDist4overlap", "", "1"),
                                      HTML('<p>Default: 1 bp</p>'),
                                      HTML('<p>Note: Distances are calculated between centers of peaks</p>')),
                     conditionalPanel(condition="input.plotType==3",
                                      h4("Perform functional enrichment analysis"),
                                      radioButtons("modeOfGeneSetsSelection", label = h5("Mode of selection"),
                                                   choices = list("Intersection" = 1, "Union" = 2, "Subtract"=3), 
                                                   selected = 1),
                                      conditionalPanel(condition="input.modeOfGeneSetsSelection==3", 
                                                       h5("Subtract from GeneSet1", style = "color:red")),
                                      actionButton("goButton", "Start Functional Enrichment Analysis"),
                                      HTML('<p>We recommend using Fisher\'s exact test only for promoter region.  Using it with any of the other locus definitions may result in biased enrichment results.  Please download the peak coordinate from other regions you are interested in from our application and upload them to <a href="http://bejerano.stanford.edu/great/public/html/">GREAT</a> to do both hypergeometric test and binomial test.</p>'),
                                      uiOutput("GeneSet1"),
                                      uiOutput("GeneSet2"),
                                      uiOutput("GeneSet3"),
                                      uiOutput("GeneSet4"),
                                      uiOutput("GeneSet5"),
                                      uiOutput("GeneSet6"),
                                      uiOutput("GeneSet7"),
                                      uiOutput("GeneSet8"),
                                      uiOutput("GeneSet9"),
                                      uiOutput("GeneSet10")),
                     
                     h4("Plot options"),
                     checkboxInput("changePeakSetsLabel", "Change peak sets label", FALSE),
                     conditionalPanel(condition="input.changePeakSetsLabel",
                                      textInput("myLabels", "", ""),
                                      HTML('<p>Example: Set 1, Set 2, Set 3</p>')
                     ),
                     checkboxInput("changePeakSetsColor", "Change colors", FALSE),
                     conditionalPanel(condition="input.changePeakSetsColor",
                                      textInput("myColours", "", c("deepskyblue, lightcoral, lightsalmon, tan3")),
                                      HTML('<p>Example: deepskyblue, lightcoral, lightsalmon 
                                           Check <a href="http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf">R Colors</a> for more colors</p>')
                                      ),
                     checkboxInput("plotSize", "Adjust plot size", FALSE),
                     conditionalPanel(condition="input.plotSize",
                                      numericInput("myHeight", "Plot height:", value=550),
                                      numericInput("myWidth", "Plot width:", value=750)
                     ),
                     checkboxInput("changeFont", "Change font", FALSE),
                     conditionalPanel(condition="input.changeFont",
                                      textInput("Font", "", "Bookman"),
                                      HTML('<p>Check <a href="http://www.cookbook-r.com/Graphs/Fonts/#table-of-fonts">R Fonts</a> for more fonts</p>')
                     ),
                     sliderInput("fontSizes",
                                 "Change font sizes:",
                                 min = 0.1,  max = 3, value = 1)
                     )	
    ),
  
  mainPanel(
    tabsetPanel(id="tabs1",
      # Welcome tab
      tabPanel("About",
               HTML('<p> <br> This application was developed to help biologists visualize their ChIP-seq/ChIP-exo results and generate new hypothesis or validate their hypothesis. We hope that you find the annoPeakR useful and we welcome suggestions
                    for additional features by our users. We would like to thank everyone who has made constructive suggestions so far. We will document the addition of new features in the News tab.</p>
                    <p>This application allows users to generate multiple types of plots to compare multiple ChIP-seq/ChIP-exo experiments. Peak sets from ChIP-seq/ChIP-exo experiments
                    can be uploaded in a bed format. Five analysis modules are presented as described in our paper. Additional features become available when checking that option. Plots can be labeled, customized (colors, dimensions) and exported as eps, pdf and svg files.  <br> </p>'),
               h4("Software references"),
               HTML('<p>R Development Core Team. <i><a href="http://www.r-project.org/">R</a>:  A Language and Environment for Statistical Computing.</i> R Foundation for Statistical Computing, Vienna (2013) <br>
                    RStudio and Inc. <i><a href="http://www.rstudio.com/shiny/">shiny</a>: Web Application Framework for R.</i> R package version 0.5.0 (2013) <br> 
                    P. Aboyoun, H. Pages and M. Lawrence. <i><a href="http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html">GenomicRanges</a>: Representation and manipulation of genomic intervals.</i> R package version 1.20.5 <br>
                    Hanbo Chen. <i><a href="https://cran.r-project.org/web/packages/VennDiagram/index.html">VennDiagram</a>: Generate high-resolution Venn and Euler plots.</i> R package version 1.6.9 <br>
                    Hadley Wickham. <i><a href="http://docs.ggplot2.org/current/">ggplot2</a>: a plotting system for R.</i> R package version 1.0.1 <br>
                    Erich Neuwirth. <i><a href="http://cran.r-project.org/web/packages/RColorBrewer/index.html">RColorBrewer</a>: ColorBrewer palettes.</i> R package version 1.1-2<br>
                    Hadley Wickham. <i><a href="https://cran.r-project.org/web/packages/reshape2/index.html">reshape2</a>: Flexibly Reshape Data.</i> R package version 1.4.1 <br>
                    Lihua Julie Zhu. <i><a href="http://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html">ChIPpeakAnno</a>: Batch annotation of the peaks identified from either ChIP-seq, ChIP-chip experiments or any experiments resulted in large number of chromosome ranges.</i> R package version 3.2.2 <br>
                    David B. Dahl. <i><a href="https://cran.r-project.org/web/packages/xtable/index.html">xtable</a>: Export tables to LaTeX or HTML.</i> R package version 1.7-4 <br>
                    Steffen Durinck and Wolfgang Huber. <i><a href="http://bioconductor.org/packages/release/bioc/html/biomaRt.html">biomaRt</a>: Interface to BioMart databases.</i> R package version 2.24.0 <br>
                    Martin Morgan, Seth Falcon, and Robert Gentleman. <i><a ref="http://www.bioconductor.org/packages/release/bioc/html/GSEABase.html">GSEABase</a>: Gene set enrichment data structures and methods.</i> R package version 1.30.2 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/GO.db.html">GO.db</a>: A set of annotation maps describing the entire Gene Ontology.</i> R package version 3.1.2 <br>
                    R. Gentleman and S. Falcon. <i><a ref="http://www.bioconductor.org/packages/release/bioc/html/GOstats.html">GOstats</a>: Tools for manipulating GO and microarrays.</i> R package version 2.34.0 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm9.knownGene.html">TxDb.Mmusculus.UCSC.mm9.knownGene</a>: Annotation package for TxDb object(s).</i> R package version 3.1.2 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm10.knownGene.html">TxDb.Mmusculus.UCSC.mm10.knownGene</a>: Annotation package for TxDb object(s).</i> R package version 3.1.2 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg18.knownGene.html">TxDb.Hsapiens.UCSC.hg18.knownGene</a>: Annotation package for TxDb object(s).</i> R package version 3.1.2 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html">TxDb.Hsapiens.UCSC.hg19.knownGene</a>: Annotation package for TxDb object(s).</i> R package version 3.1.2 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html">TxDb.Hsapiens.UCSC.hg38.knownGene</a>: Annotation package for TxDb object(s).</i> R package version 3.1.2 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html">org.Hs.eg.db</a>: Genome wide annotation for Human.</i> R package version 3.1.2 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html">org.Mm.eg.db</a>: Genome wide annotation for Human.</i> R package version 3.1.2 <br>
                    Marc Carlson. <i><a ref="http://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html">GeneSigDB</a>: a manually curated database and resource for analysis of gene expression signatures.</i> Nucl. Acids Res. (2012) 40 (D1): D1060-
                    D1066. </p>'),
               h6("This application was created by ", a("Xing Tang", href="https://scholar.google.com/citations?user=F9arEIMAAAAJ&hl=en"), 
                  " from The Ohio State University. Please send bugs and feature requests to Xing Tang (tangx1986(at)gmail.com) or (tang.811(at)osu.edu). This application uses the ", 
                  a("shiny package from RStudio", href="http://www.rstudio.com/shiny/"), ".")
               ),
      # Data upload tab
      tabPanel("Data upload",
               h5("The uploaded files should be in", a(".bed", href="https://genome.ucsc.edu/FAQ/FAQformat.html#format1"), "format.", style = "color:blue"),
               h5("The first four columns are required. They are \"chromosome\", \"start\", \"end\", \"peak_id\". " , style = "color:blue"),
               uiOutput("tables"), 
               h6("This application was created by the ", a("Xing Tang", href="https://scholar.google.com/citations?user=F9arEIMAAAAJ&hl=en"), 
                  " from The Ohio State University. Please send bugs and feature requests to Xing Tang (tangx1986(at)gmail.com) or (tang.811(at)osu.edu). This application uses the ", 
                  a("shiny package from RStudio", href="http://www.rstudio.com/shiny/"), ".")
      ),
      
      # Boxplot tab
      tabPanel("Data visualization",
               
               conditionalPanel(condition="input.plotType==4", 
                                downloadButton("downloadPlotEPS_length", "Download eps-file"),
                                downloadButton("downloadPlotPDF_length", "Download pdf-file"),
                                downloadButton("downloadPlotSVG_length", "Download svg-file"),
                                plotOutput("peakLengthDensity", height='100%', width='100%')
               ),
               conditionalPanel(condition="input.plotType==0", 
                                downloadButton("downloadPlotEPS_hist", "Download eps-file"),
                                downloadButton("downloadPlotPDF_hist", "Download pdf-file"),
                                downloadButton("downloadPlotSVG_hist", "Download svg-file"),
                                plotOutput("hist", height='100%', width='100%'),
                                downloadButton("downloadTargetGenes", "Download target genes data as .CSV file"),
                                h4("Peak distribution plot statistics"), tableOutput("histStatsTable"),
                                downloadButton("downloadHistData", "Download peak region wise distribution data as .CSV file")
               ),
               conditionalPanel(condition="input.plotType==1", 
                                downloadButton("downloadPlotEPS_pair", "Download eps-file"),
                                downloadButton("downloadPlotPDF_pair", "Download pdf-file"),
                                downloadButton("downloadPlotSVG_pair", "Download svg-file"),
                                h5("Use this module to survey the changes of binding positions between two peak sets.", style = "color:blue"),
                                p("Check \"FAQ\" tab for detailed information.", style = "color:blue"),
                                plotOutput("pairPeakPlot", height='100%', width='100%')
               ),
               conditionalPanel(condition="input.plotType==2", 
                                h4("Pair-wise overlap counts"),
                                tableOutput("peakOverlapsTable"),
                                downloadButton("downloadPeakOverlapsTable", "Download pair-wise overlap counts data as .CSV file"),
                                h4("Pair-wise overlap rates"),
                                tableOutput("peakOverlapsRateTable"),
                                downloadButton("downloadPeakOverlapsRatesTable", "Download pair-wise overlap rates data as .CSV file")
               ),
               conditionalPanel(condition="input.plotType==3", 
                                h4("Show overlaps of gene sets by venn diagram"),
                                downloadButton("downloadPlotEPS_geneVenn", "Download eps-file"),
                                downloadButton("downloadPlotPDF_geneVenn", "Download pdf-file"),
                                downloadButton("downloadPlotSVG_geneVenn", "Download svg-file"),
                                plotOutput("geneVennPlot", height='100%', width='100%'),
                                downloadButton("downloadGeneMergedTable", "Download merged table for target genes as .CSV file"),
                                h4("Pair-wise overlap counts"),
                                tableOutput("geneOverlapsTable"),
                                downloadButton("downloadGeneOverlapsTable", "Download pair-wise overlap counts data as .CSV file"),
                                h4("Pair-wise overlap rates"),
                                tableOutput("geneOverlapsRateTable"),
                                downloadButton("downloadGeneOverlapsRatesTable", "Download pair-wise overlap rates data as .CSV file"),
                                h4("Enriched GO terms"),
                                tableOutput("GOTable"),
                                downloadButton("downloadEnrichedGO", "Download full table of enriched GO terms as .CSV file"),
                                h4("Enriched gene signatures"),
                                tableOutput("GSEATable"),
                                h6("All the gene signatures are from  ", a("GeneSigDB.", href="http://cccb.dfci.harvard.edu/genesigdb/"), "GeneSigDB collects gene signatures from literature. Check the detailed information of each enriched signature by the provided literature ID."),
                                downloadButton("downloadEnrichedSig", "Download full table of enriched gene signatures as .CSV file")
                                ),
               
               h6("This application was created by ", a("Xing Tang", href="https://scholar.google.com/citations?user=F9arEIMAAAAJ&hl=en"), "from the Ohio State University. Please send bugs and feature requests to Xing (tangx1986(at)gmail.com) or (tang.811@osu.edu). This application uses the ", 
                  a("shiny package from RStudio", href="http://www.rstudio.com/shiny/"), ".")
      ), 
      # Figure legend 
      
      # News
      tabPanel("News",
               h5("July 24, 2015"), 
               p("version 1.0.0")
      ),			
      
      # FAQ 
      tabPanel("FAQ",
               h5("Q: Can user do motif analysis with this application?"),
               p("A: Currently we didn't integrate this analysis directly. A simple solution for users is to
                download annotated peaks from \"Peak associated gene structures\" module of our application
                and upload interested peaks to http://embnet.ccg.unam.mx/rsa-tools/ to do motif analysis with
                \"RSAT peak-motif\" algorithm. User can also try other motif discovery tools, such as MEME, HOMER and etc.."),
               h5("Q: How to compare peaks only from a specific region (eg, promoter)?"),
               p("A: If you want to analyze peaks only from promoter region using annoPeakR, 
                  you can download the \"target genes data\" file from \"Peak regionwise distribution\" module.
                 Open it with excel, then extract and save peaks from promoter region into a file in bed format. 
                 Finally upload the saved peaks onto annoPeakR and do the analysis."),
               h5("Q: How to calculate nearest peak distance?"),
               p("A: For each pair of peak sets, we calculate two nearest peak distance vectors.
                 As the example shows below, we calculate nearest peak distance vectors for both set A and set B.
                 For each peak in set A, we find the nearest peak in set B and vice versa.
                 In this example, we name each peak by the peak set ID (A/B) and the position on chromosome (eg, A.100).
                 The distance vector from A to B is (35, 25, 85). 
                 The distance vector from B to A is c(35, 25). 
                 The distance vectors are used to generate the density plot in the \"Nearest peak distances distribution\" module.
                 Curves for \"A to B\" and \"B to A\" are colored differently in the density plot."),
               img(src = "pair_peak_desc.jpg", height = 160, width = 350),
               h5("Q: I have trouble editing the graphic files."), 
               p("A: For EPS files make sure to 'ungroup' all objects so they can be edited independently. 
                 In Adobe Illustrator you will also need to use the 'release compound path' command. For PDF 
                 files you should 'release clipping mask'. SVG import appears to have problems in Adobe Illustrator 
                 and Corel Draw and should be avoided. EPS, PDF and SVG import all work with Inkscape http://www.inkscape.org/.")
                             )
      
               )
      )
               )
    )







