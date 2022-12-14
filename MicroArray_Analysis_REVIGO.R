#!/MicroArray_Analysis_REVIGO.R

#The script takes a list of CEL files from microarray Affymetrix platform and analyse them.
#The script identifies the differential expressed terms with limma.
#Enrichment analysis is carried out with g:Profiler.
#Heatmap plot are created to show the similarity between each group comparison by using REVIGO.

#The main points of the script was guided from here -> https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor
#Additional Theory -> https://data.bits.vib.be/pub/trainingen/MicroarraysIntro/Theory.pdf

#The inputs are:
#1)sampleTarget -> path to the file where the CEL filename and its group are saved (#Filename #Group).
#2)groupComparison -> List of group comparison written as "Group1-Group2".
#3)enrich_org -> Scientific name of the organism to whom enrichment analysis will bbe carried out with g:Profiler (https://biit.cs.ut.ee/gprofiler/page/organism-list)

#The REVIGO inputs are:
#4)cutoff <- Cut-off allowed values: "0.90" "0.70" "0.50" "0.40" 
#5)goSizes <- "0" #("whole UniProt")
#6)measure <- "SIMREL" #Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG" 

#The outputs are:

#Quality Control Directories
#1)QC_Raw_Intensity shows the raw intensities for each sample.
#2)QC_Chip_PseudoImage_weights shows the pseudo-images based on the weights for each sample.
#3)QC_Chip_PseudoImage_residuals shows the pseudo-images based on the residuals for each sample.
#4)QC_Chip_PseudoImage_pos_residuals shows the pseudo-images based on the positive residuals for each sample where only high positive residuals are drawn in red, while negative and near 0 residuals being drawn in white.
#5)QC_Chip_PseudoImage_neg_residuals shows the pseudo-images based on the negative residuals for each sample where only extreme negative residuals are drawn in blue, while positive negative and near 0 residuals being drawn in white.
#6)QC_Chip_PseudoImage_sign_residuals shows the pseudo-images based on the sign of the residuals for each sample where all negative residuals regardless of magnitude are indicated by blue and all positive residuals by red.
#7)QC_MAplots_rawData shows the MA plots for the raw data.
#8)QC_MAplots_norm shows the MA plots for the normalised data.

#Dendogram of the DEG analysis
#9)dendogram.pdf displays the dendogram of the chosen group comparisons based on the identified DEGs.

#Quality Control Plots file
#10)quality_control_plots.pdf file shows histogram of the intensities of the raw data of the samples, boxplot of the raw data intensities after and before the normalisation for each samples, and PCA of the normalised data.

#Enrichment Analysis Plots
#11)The gProfiler_revigo_heatmap_plots.pdf file displays the heatmaps from the enrichment analysis of the chosen group comparisons (GO:BP, GO:CC, GO:MF, KEGG). The description of the heatmap row is based on the main GO term of the cluster; followed by the top 5 words of all the GO terms' description present in the cluster; followed by the total number of GO terms in the cluster.

# CSV files
#12)all_info.csv shows the information for all the probes in the microarray.
#13)rma_norm.csv shows the normalised data of the probe's intensities.
#14)deg_analysis_duplicates.csv shows the DEGs with duplicated genes since multiple probes are used for one single gene.
#15)deg_unique_group.csv shows the unique list of DEGs for each group comparison.
#16)gProfiler_enrich.csv shows all enriched terms for each group comparison.
#17)revigoCluster_GO:BP.csv shows the GO clusters found for GO:BP.
#18)revigo_GO:BP.csv shows the results after the REVIGO analysis for GO:BP analysis.
#19)revigoClusterHeatmap_GO:BP.csv shows the input file for the heatmap generation for the GO:BP analysis.
#20)revigoCluster_GO:CC.csv shows the GO clusters found for GO:CC.
#21)revigo_GO:CC.csv shows the results after the REVIGO analysis for GO:CC analysis.
#22)revigoClusterHeatmap_GO:CC.csv shows the input file for the heatmap generation for the GO:CC analysis.
#23)revigoCluster_GO:MF.csv shows the GO clusters found for GO:MF.
#24)revigo_GO:MF.csv shows the results after the REVIGO analysis for GO:MF analysis.
#25)revigoClusterHeatmap_GO:MF.csv shows the input file for the heatmap generation for the GO:MF analysis.
#26)heatmapKEGG.csv shows the heatmap results for the KEGG enrichment analysis.


#Load Library
library(oligo) 			#Preprocessing tools for oligonucleotide arrays
library(affycoretools)	#Functions for those repetitive analyses with Affymetrix GeneChips
library(hgu133plus2.db)	#Affymetrix Human Genome U133 Plus 2.0 Array annotation data
library(factoextra)		#Multivariate Data Analyses and Elegant Visualization
library(tidyverse)		#Collection of R packages designed for data science
library(gprofiler2)		#Gene list functional enrichment analysis
library(rvest)			#Download, then manipulate, HTML and XML
library(tidytext)		#Text mining tasks and plot generation.
library(scales)			#Graphical scales map data to aesthetics
library(stringi)		#Language processing tools
library(gridExtra)		#Arrange multiple grid-based plots on a page


###### ##### ##### USER's INPUT ###### ##### #####
#Sample metadata filename
sampleTarget <- "metadata.txt"

#Choose the Group Comparisons 
groupComparison <- c("HGP-Control", "UV-Control", "HGPTERT-Control", "ControlTERT-Control")

#Scientific Name of the organism to perform enrichment analysis with g:Profiler (https://biit.cs.ut.ee/gprofiler/page/organism-list)
enrich_org <- "Homo sapiens" 

###REVIGO's Input format
cutoff <- "0.40" #Allowed values: "0.90" "0.70" "0.50" "0.40" 
goSizes <- "0" #("whole UniProt")
measure <- "SIMREL" #Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"

#################################################

###### ##### ##### General Object ###### ##### #####
qc_plots <- list ()   #Quality Control Plots
enrichPlot <- list()  ##Enrichment Analysis Plots
baseurl <- "http://revigo.irb.hr/"  #REVIGO website link
isPValue <- "yes" #Allowed values in REVIGO: "yes"  "no"
whatIsBetter <- "higher" #Allowed values in REVIGO: "higher" "lower" "absolute" "abs_log"

###### ##### ##### START ANALYSIS ###### ##### #####
##### CEL files
data <- read.celfiles(list.celfiles())
n_samples <- length(list.celfiles())

##### Samples Groups
targets <- readTargets (sampleTarget, sep="")
ph <- data@phenoData
targets <- targets[order(targets$Filename),] 
ph@data[ ,2] <- targets$Group
colnames(ph@data)[2] <- "source"
sample <- rownames (ph@data)
ph@data <- cbind(ph@data,sample)

f <- factor(ph@data$source)
design <- model.matrix(~ 0 + f) 
colnames(design) <- levels(f)

# Sort the samples based on their group membership
targets <- targets [order(targets$Group),]



###### ##### ##### QUALITY CONTROL  MICROARRAY FILEs ###### ##### #####
###### Normalisation using RMA
data.rma <- rma(data)
data.matrix <- exprs(data.rma)

###### Raw Intensities of each array
dir.create("QC_Raw_Intensity")
for (i in 1:n_samples)
{
  name = paste("QC_Raw_Intensity/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(data[,i],main=ph@data$sample[i])
  dev.off()
}

###### Chip pseudo-images based on weights
Pset <- fitProbeLevelModel(data)
#Pset <- fitPLM(data,output.param=list(varcov="none"))  #In case of a large dataset

dir.create("QC_Chip_PseudoImage_weights")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_weights/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type='weights', main=ph@data$sample[i])
  dev.off()
}

###### Chip pseudo-images based on residuals -> 4 Types: residuals, positive residuals, negative residuals or sign of residuals

#residuals -> gives a pseudo-image of residuals
dir.create("QC_Chip_PseudoImage_residuals")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_residuals/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type="residuals", main=ph@data$sample[i])
  dev.off()
}

#pos.resids -> only high positive residuals are drawn in red, while negative and near 0 residuals being drawn in white
dir.create("QC_Chip_PseudoImage_pos_residuals")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_pos_residuals/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type="pos.residuals", main=ph@data$sample[i])
  dev.off()
}

#type="neg.resids" only extreme negative residuals are drawn in blue, while positive negative and near 0 residuals being drawn in white
dir.create("QC_Chip_PseudoImage_neg_residuals")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_neg_residuals/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type="neg.residuals", main=ph@data$sample[i])
  dev.off()
}

#type="sign.resids" gives images where all negative residuals regardless of magnitude are indicated by blue and all positive residuals by red
dir.create("QC_Chip_PseudoImage_sign_residuals")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_sign_residuals/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type="sign.residuals", main=ph@data$sample[i])
  dev.off()
}


###### Histogram -> use the PM intensities
pmexp <- pm(data)
sampleNames <- vector()
logs <- vector()

for (i in 1:n_samples)
{
  sampleNames <- c(sampleNames, rep(ph@data$sample[i], dim(pmexp)[1]))
  logs <- c(logs, log2(pmexp[,i]))
}

logData <- data.frame(logInt=logs, sampleName=sampleNames)
logData$sampleName <- factor(logData$sampleName, levels = targets$Filename)

qc_plots[["Histogram_Raw"]] <- logData %>% 
  ggplot(aes(logInt, colour = sampleName)) + 
  geom_density() +
  ggtitle("Histogram")


###### Boxplots on the normalised data
sampleNames = vector()
normlogs = vector()

for (i in 1:n_samples)
{
  sampleNames = c(sampleNames,rep(ph@data$sample[i],dim(data.matrix)[1]))
  normlogs = c(normlogs,data.matrix[,i])
}

norphata <- data.frame(norm_logInt=normlogs,sampleName=sampleNames)
norphata$sampleName <- factor(norphata$sampleName, levels = targets$Filename)


qc_plots[["Boxplot_After_Normalization"]] <- norphata %>% 
  ggplot(aes(sampleName,norm_logInt)) + geom_boxplot() + ylim(2,16) + ggtitle("after normalization") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8, face = "bold"))

qc_plots[["Boxplot_Before_Normalization"]] <- logData %>% 
  ggplot(aes(sampleName,logInt)) + geom_boxplot() + ylim(2,16) + ggtitle("before normalization") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8, face = "bold"))


###### MA PLot: raw data
dir.create("QC_MAplots_rawData")
for (i in 1:n_samples)
{
  name = paste("QC_MAplots_rawData/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  MAplot(data,which=i)
  dev.off()
}


###### MA PLot: normalised data
dir.create("QC_MAplots_norm")
for (i in 1:n_samples)
{
  name = paste("QC_MAplots_norm/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  MAplot(data.rma,which=i)
  dev.off()
}

###### PCA
sort_targets <- targets [order(targets$Filename),]
color <- sort_targets$Group
data.PC <- prcomp(t(data.matrix), scale = TRUE, center = TRUE)

qc_plots[["PCA"]] <- fviz_pca_ind(data.PC, 
                                  palette = c("green",  "blue", "red", "darkorange", "black", "purple", "aquamarine", "darkgoldenrod1"),
                                  col.ind = sort_targets$Group,  #color by groups
                                  legend.title = "Groups",
                                  addEllipses = TRUE, # Concentration ellipses
                                  ellipse.type = "confidence",
                                  label = FALSE,
                                  repel = TRUE)



###### ##### ##### MICROARRAY ANALYSIS ###### ##### #####

##### Normalisation methods: rma, mas5, germa, justPlier, expresso
data.rma <- rma(data)
data.annot <- annotateEset(data.rma, hgu133plus2.db)  #ProbeID Annotation
data.fit <- lmFit(data.annot,design)

##### Choose the group comparisons
contrast.matrix <- makeContrasts(contrasts=groupComparison, levels=design) 
data.fit.con <- contrasts.fit(data.fit,contrast.matrix)
data.fit.eb <- eBayes(data.fit.con)
data_final <- as.data.frame(data.fit.eb)

##### Adjusting for multiple testing and defining DE genes
DEresults <- decideTests (data.fit.eb, method='global', adjust.method="BH", p.value=0.000005,lfc=1)
DEresults <- rownames_to_column(as.data.frame (DEresults), "genes.PROBEID")
all_info <- merge (x = data_final, y = DEresults, by = "genes.PROBEID", all.y=TRUE)
write.csv(all_info, file = "all_info.csv")

##### Save the normalised data #####
norm <- as.data.frame(data.rma)
tnorm <- as.data.frame (t(norm))
tnorm <- rownames_to_column(tnorm, "PROBEID")
annot <- data.annot@featureData@data
tnorm_annot <- merge(x=tnorm, y=annot, by="PROBEID", all.y=TRUE)
write.csv(tnorm, file = "rma_norm.csv")

##### Create the right DF frame from the "all_info" dataframe
colnames(all_info) <- make.names(colnames(all_info))
all_info$DEG <- rowSums(abs(all_info[(ncol(all_info)-ncol(contrast.matrix)+1) : ncol(all_info)]))

resultDF <- all_info %>% filter(genes.ENTREZID != "NA") %>%    #Remove rows where the ENTREZID is "NA"
                        filter(DEG > 0) %>%                   #Filter for rows where the gene is DE
                        select(starts_with(c("genes", "p.value", "coefficients")), make.names(colnames(contrast.matrix))) #Select columns
                               
write.csv(resultDF, file = "deg_analysis_duplicates.csv")                     

##### Get the unique list of DEGs for each group #####
deg_unique_group <- list()
uniqueDEG <- c()

for (i in make.names(colnames(contrast.matrix))) 
{
  uniDEG <- resultDF %>% filter(eval(parse(text = i)) != 0)
  uniDEG <- unique(uniDEG$genes.ENTREZID)
  deg_unique_group[[i]] <- uniDEG
  uniqueDEG <- append(uniqueDEG, uniDEG)
}

deg_unique_group_DF <- as.data.frame(t(plyr::ldply(deg_unique_group, rbind)))
write.csv(deg_unique_group_DF, file = "deg_unique_group.csv")

##### Create clustering
uniqueDEG <- unique(uniqueDEG)
degClusterDF <- data.frame(matrix(0, ncol=length(colnames(contrast.matrix)), nrow=length(uniqueDEG)))
colnames(degClusterDF) <- make.names(colnames(contrast.matrix))
rownames(degClusterDF) <- uniqueDEG

for (i in names(deg_unique_group)) 
{for (l in deg_unique_group[[i]]) degClusterDF[l, i] <- 1}

dist_mat <- dist(t(degClusterDF), method = 'euclidean')

# Save the Cluster Dendogram plot
pdf("dendogram.pdf") 
plot(hclust(dist_mat, method = 'average'))
dev.off() 


##### Enrichment Analysis -> g:Profiler #####
gProfiler_Res <- NULL
gPr_iSpecies <- paste(substr (tolower(enrich_org), 1, 1), str_split (tolower(enrich_org), pattern = " ")[[1]][2], sep = "")

tryCatch (
  {
    gProfiler_Res <- gost(deg_unique_group, organism = gPr_iSpecies, ordered_query = FALSE,
                          multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                          measure_underrepresentation = FALSE, evcodes = FALSE,
                          user_threshold = 0.05, correction_method = "gSCS",
                          domain_scope = "annotated", custom_bg = NULL,
                          numeric_ns = "", sources = NULL)
  },
  error = function(e) 
  {
    message(e)
    return(gProfiler_Res <- NULL)
  }
)

if (is.null(gProfiler_Res) == FALSE) 
{
  enrich_go_kegg <- gProfiler_Res$result[,c(1,3,9,10,11)] %>% filter(source %in% c("GO:BP", "GO:CC", "GO:MF", "KEGG"))
  write.csv(enrich_go_kegg, file = "gProfiler_enrich.csv")
  
  ##### REVIGO Analysis #####
  
  # REVIGO's Input format
  revigo_session <- html_session(baseurl)
  revigo_form <- html_form(revigo_session)[[1]] 
  
  # Loop through the GO terms
  for (go in c("GO:BP", "GO:CC", "GO:MF")) 
  {
    # Number of Cluster
    nCluster <- 40
    
    # Gene Ontology
    filterGO <- enrich_go_kegg %>% filter (source==go)
    goList <- paste(unique(filterGO$term_id), collapse="\n")
    
    filled_form <- set_values(revigo_form, 'goList'= goList, 
                              'cutoff'= cutoff, 'isPValue'= isPValue, 
                              'whatIsBetter'= whatIsBetter, 'goSizes'= goSizes, 'measure'= measure)
    
    # Get the Revigo Summary Table & Save it as a CSV file
    result_page <- submit_form(revigo_session, filled_form, submit='startRevigo')
    GORevigo <- result_page %>% follow_link("Export results to text table (CSV)")
    httr::content(GORevigo$response, as="text") %>% write_lines("revigoSummaryOnline.csv") 
    
    # Get the REVIGO summary result from the saved file
    revigoSummary <- read.csv2("revigoSummaryOnline.csv", header=TRUE, sep=",")
    revigoSummary$frequency <- gsub('%', '', revigoSummary$frequency)
    revigoSummary$frequency <- as.numeric( as.character(revigoSummary$frequency) );
    
    #Get only the GO terms represented in the "semantic space"
    uniqueGO <- revigoSummary %>% filter(plot_X != "null")
    iNum <- c(4,5,6,7,8) # Change the column into numeric: plot_X - plot_Y - plot_size - uniqueness - dispensability
    uniqueGO [ , iNum] <- apply(uniqueGO [ , iNum], 2, function(x) as.numeric(as.character(x)))
    
    # Add "HeadGO" Column which shows which GO REVIGO term has been put as the main GO for that group of terms
    revigoSummary$HeadGO <- NA
    headGO <- "NA"
    
    for (i in 1:nrow(revigoSummary)) 
    {
      if(revigoSummary$plot_X[i]!="null") {headGO <- revigoSummary$term_ID[i]}
      revigoSummary$HeadGO[i] <- headGO
    }
    write.csv(revigoSummary, file = sprintf("revigo_%s.csv", go))   #Save the REVIGO result
    
    
    ####################### GENE ONTOLOGY REVIGO CLUSTER ANALYSIS   ####################### 
    #Find the clusters present in the REVIGO results based on their sematic space
    clusterDF <- uniqueGO %>% remove_rownames %>% column_to_rownames(var="term_ID")
    if (nCluster >= nrow(clusterDF)) {nCluster <- nrow(clusterDF)-1}
    
    set.seed(100)
    clusters <- kmeans(clusterDF[,c(3,4)], nCluster)
    clusterDF$Cluster <- clusters$cluster
    
    # Group elements from the same cluster
    revigoCluster <- data.frame(matrix(ncol=10, nrow=nCluster))
    colnames(revigoCluster) <- c("Cluster", "PlotX", "PlotY", "RevigoGOs", "RevigoRep", "RevigoDescription", "AllGOs", "AllWords", "AllDescription", "FinalDescription")
    
    for (i in sort(unique(clusterDF$Cluster)))
    {
      filterGO_cluster <- clusterDF %>% filter(Cluster==i)
      
      # REVIGO GOs
      revigoCluster$Cluster[i] <- i
      revigoCluster$PlotX[i] <- mean(filterGO_cluster$plot_X)
      revigoCluster$PlotY[i] <- mean(filterGO_cluster$plot_Y)
      revigoCluster$RevigoGOs[i] <- paste (rownames(filterGO_cluster), collapse=" ")
      revigoCluster$RevigoRep[i] <- filterGO_cluster[which.max(filterGO_cluster$representative),][1,1] #Term with the max value of "representative of the cluster
      
      # REVIGO Cluster Description -> Concagenate the HEADs' term description of the cluster - Split - Delete generic terms (i.e. "of", "a", "an")
      wordOcc <- paste(filterGO_cluster$description, collapse=" ") %>% str_split(pattern = " ") 
      wordOcc <- wordOcc[[1]]
      indexWords <- !(wordOcc %in% stop_words$word)
      wordOcc <- as.data.frame(table (wordOcc[indexWords])) %>% arrange(desc (Freq))
      wordOcc <- as.vector(wordOcc$Var1[1:5])
      revigoCluster$RevigoDescription[i] <- paste (wordOcc[!is.na(wordOcc)], collapse = " ")
      
      #All GOs Columns
      filterAllGO <- revigoSummary %>% filter(HeadGO %in% rownames(filterGO_cluster))
      
      revigoCluster$AllGOs[i] <- paste (filterAllGO$term_ID, collapse=" ")
      revigoCluster$AllWords[i] <- paste(filterAllGO$description, collapse=" ")
      
      wordOcc <- paste(filterAllGO$description, collapse=" ") %>% str_split(pattern = " ") 
      wordOcc <- wordOcc[[1]]
      indexWords <- !(wordOcc %in% stop_words$word)
      wordOcc <- as.data.frame(table (wordOcc[indexWords])) %>% arrange(desc (Freq))
      wordOcc <- as.vector(wordOcc$Var1[1:5])
      revigoCluster$AllDescription[i] <- paste (wordOcc[!is.na(wordOcc)], collapse = " ")
      
      #Final Description of the Cluster
      revigoCluster$FinalDescription[i] <- paste(revigoCluster$RevigoRep[i], revigoCluster$AllDescription[i], nrow(filterAllGO), sep = " * ")
    }
    write.csv(revigoCluster, file = sprintf("revigoCluster_%s.csv", go))   #Save the REVIGO_Cluster result
    
    
    #Get the lowest pvalue of the group comparison's GO terms and the cluster combination
    heatmapGO <- data.frame(matrix(ncol = length (unique(filterGO$query)), nrow = nrow(revigoCluster)))
    colnames(heatmapGO) <- unique(filterGO$query)
    rownames(heatmapGO) <- revigoCluster$FinalDescription
    
    for (n in 1:nrow(revigoCluster)) 
    {
      allGO <- str_split(revigoCluster$AllGOs[n], " ") [[1]] #Get all GOs 
      filterGOGroup <- filterGO %>% filter(term_id %in% allGO)  #Get all the rows which have those GOs
      
      uniqueGroup <- unique (filterGOGroup$query)
      for (k in uniqueGroup) 
      {
        pvaluesList <- filterGOGroup %>% filter(query == k)
        heatmapGO[n, grep(k, colnames(heatmapGO))] <- min(pvaluesList$p_value)
      }
    }
    write.csv(heatmapGO, file = sprintf("revigoClusterHeatmap_%s.csv", go))   #Save the REVIGO_Cluster_Heatmap result
    
    
    #Create the Heatmap plot
    GroupName <- c()
    for (i in 1:nrow(heatmapGO)) 
    {
      groupPresence <- colnames(heatmapGO)[sapply(heatmapGO[i,], function(x) any(!is.na(x)))]
      groupPresence <- paste(groupPresence, collapse = " * ")
      GroupName <- append(GroupName, groupPresence)
    }
    
    heatmapGO$GroupName <- GroupName
    heatmapGO$NGroups <- length (unique(filterGO$query)) - rowSums(is.na(heatmapGO))
    heatmapGO <- heatmapGO %>% arrange(NGroups, GroupName) %>% dplyr::select(-c(NGroups, GroupName))
    
    heatmapGO <- tibble::rownames_to_column(heatmapGO, "Description")
    heatmapGO_melt <- reshape2::melt(heatmapGO,id.vars=c("Description"))
    
    heatmapGO_melt$Description <- factor(heatmapGO_melt$Description, levels = heatmapGO$Description)
    
    plotName <- stri_replace_all_fixed (paste(enrich_org, go), pattern = c(" "), replacement = c("_"), vectorize_all = FALSE)
    enrichPlot[[plotName]] <- heatmapGO_melt %>% 
      ggplot(aes(x=variable, y=factor(Description), fill=value)) +
      geom_raster(aes(fill = value)) +
      scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
      geom_text(size=3, aes(label=scientific(value, digits = 3))) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
      theme(axis.text.y = element_text(size=10, face="bold")) +
      coord_fixed(ratio=ncol(heatmapGO)/nrow(heatmapGO))
  }
  
  ##### KEGG Heatmap creation #####
  filterKEGG <- filter(enrich_go_kegg, source=="KEGG")
  heatmapKEGG <- data.frame(matrix(ncol = length (unique(filterKEGG$query)), nrow = length (unique (filterKEGG$term_id))))
  colnames(heatmapKEGG) <- unique(filterKEGG$query)
  rownames(heatmapKEGG) <- unique (filterKEGG$term_name)
  
  for (i in 1:nrow(filterKEGG)) 
  {heatmapKEGG[filterKEGG$term_name[i], filterKEGG$query[i]] <- filterKEGG$p_value[i]}
  write.csv(heatmapKEGG, file = "heatmapKEGG.csv")   #Save the gProfiler_Heatmap result
  
  
  #Create the Heatmap plot
  GroupName <- c()
  for (i in 1:nrow(heatmapKEGG)) 
  {
    groupPresence <- colnames(heatmapKEGG)[sapply(heatmapKEGG[i,], function(x) any(!is.na(x)))]
    groupPresence <- paste(groupPresence, collapse = " * ")
    GroupName <- append(GroupName, groupPresence)
  }
  
  heatmapKEGG$GroupName <- GroupName
  heatmapKEGG$NGroups <- length (unique(filterGO$query)) - rowSums(is.na(heatmapKEGG))
  heatmapKEGG <- heatmapKEGG %>% arrange(NGroups, GroupName) %>% dplyr::select(-c(NGroups, GroupName))
  
  heatmapKEGG <- tibble::rownames_to_column(heatmapKEGG, "Description")
  heatmapKEGG_melt <- reshape2::melt(heatmapKEGG,id.vars=c("Description"))
  heatmapKEGG_melt$Description <- factor(heatmapKEGG_melt$Description, levels = heatmapKEGG$Description)
  heatValues <- heatmapKEGG_melt$value [!heatmapKEGG_melt$value %in% c(NA)]
  
  plotName <- stri_replace_all_fixed (paste(enrich_org, "KEGG"), pattern = c(" "), replacement = c("_"), vectorize_all = FALSE)
  enrichPlot[[plotName]] <- heatmapKEGG_melt %>% 
    ggplot(aes(x=variable, y=Description, fill=value)) +
    geom_raster(aes(fill = value)) +
    scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
    geom_text(size=3, aes(label=scientific(value, digits = 3))) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
    theme(axis.text.y = element_text(size=10, face="bold")) +
    coord_fixed(ratio=ncol(heatmapGO)/nrow(heatmapGO))
}

#Delete the "revigoSummaryOnline.csv" file
if (file.exists("revigoSummaryOnline.csv")) {file.remove("revigoSummaryOnline.csv")}

##### Save Quality Control Plots #####
pdf("quality_control_plots.pdf", onefile=TRUE)
for (i in seq(length(qc_plots))) {grid.arrange (qc_plots[[i]])}
dev.off()

##### Save Enrichment Analysis Plots #####
pdf("gProfiler_revigo_heatmap_plots.pdf", onefile=TRUE)
for (i in seq(length(enrichPlot))) {grid.arrange (enrichPlot[[i]])}
dev.off()
