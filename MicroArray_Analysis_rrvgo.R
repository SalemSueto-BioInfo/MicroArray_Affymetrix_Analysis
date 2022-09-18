#!/MicroArray_Analysis_rrvgo.R

#The script takes a list of CEL files from microarray Affymetrix platform and analyse them.
#The script identifies the differential expressed terms with limma.
#Enrichment analysis is carried out for each group comparison and their similarity is showed in a heatmap.

#The main points of the script was guided from here -> https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor
#Additional Theory -> https://data.bits.vib.be/pub/trainingen/MicroarraysIntro/Theory.pdf

#The inputs are:
#1)sampleTarget -> path to the file where the CEL filename and its group are saved (#Filename #Group).
#2)groupComparison -> List of group comparison written as "Group1-Group2".
#3)enrich_org -> Scientific name of the organism to whom enrichment analysis will bbe carried out with g:Profiler (https://biit.cs.ut.ee/gprofiler/page/organism-list)
#4)orgDB -> Name of the organism database from the Bioconductor (https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb)
#5)simSimMethod -> GO similarity terms methodology name (Allowed values: "Resnik", "Lin", "Rel", "Jiang", "Wang")
#6)thr -> threshold value for the GO reduce terms -> Allowed values: Large (allowed similarity=0.9), Medium (0.7), Small (0.5), Tiny (0.4)
#7)topWords -> #Highest number of words to describe the members of the GO cluster during the reduced step

#The outputs are:

#Quality Control Directories
#1)QC_Chip_PseudoImage_neg_residuals directory -> Pseudo image negative residuals QC plot for each CEL file.
#2)QC_Chip_PseudoImage_pos_residuals directory -> Pseudo image positive residuals QC plot for each CEL file.
#3)QC_Chip_PseudoImage_residuals directory -> Pseudo image residuals QC plot for each CEL file.
#4)QC_Chip_PseudoImage_sign_residuals directory -> Pseudo image sign residuals QC plot for each CEL file.
#5)QC_Chip_PseudoImage_weights directory -> Pseudo image weigths residuals QC plot for each CEL file.
#6)QC_MAplots_norm directory -> MA plots of the normalised data for each CEL file.
#7)QC_MAplots_rawData directory -> MA plots of the raw data for each CEL file.
#8)QC_Raw_Intensity directory -> Intensity plots of the raw data for each CEL file.

#CSV files
#9)all_info.csv -> All the info for the probes and their DEGs analysis.
#10)deg_analysis_duplicates.csv -> DEG analysis with all the gene's duplicates due to the different probes for the same gene.
#11)deg_unique_group.csv -> DEG list for each group comparison.
#12)goReducedTerms.csv -> Reduced GO terms.
#13)gProfiler_enrich.csv -> Enrichment Analysis
#14)heatmapDF.csv -> Heatmap plot dataframe
#15)rma_norm.csv -> RMA normalisation data

#PDF plot file
#16)dendogram.pdf -> Sample cluster dendogram plot 
#17)gProfiler_reducedTerms_heatmap_plots.pdf -> Enrichment analysis heatmap plots
#18)quality_control_plots.pdf -> Quality control plots of the raw and normalised data

#RData plots file
#19)QC_heatmap_plots.RData -> RData where all the plots are saved.


### Load Library ###
library(oligo)          #Preprocessing tools for oligonucleotide arrays
library(affycoretools)  #Functions for those repetitive analyses with Affymetrix GeneChips
library(tidyverse)      #Collection of R packages designed for data science
library(factoextra)     #Multivariate Data Analyses and Elegant Visualization
library(hgu133plus2.db) #Affymetrix Human Genome U133 Plus 2.0 Array annotation data
library(gprofiler2)     #Gene list functional enrichment analysis
library(rrvgo)          #Reduce + Visualize GO
library(shinydashboard) #Create a shiny version of the result from rrvgo
library(tidytext)       #Text mining tasks and plot generation.
library(gridExtra)      #Arrange multiple grid-based plots on a page
library(scales)         #Graphical scales map data to aesthetics

rm(list=ls())


### USER's INPUT ###
sampleTarget <- "metadata.txt" #Sample filename
groupComparison <- c("HGP-Control", "UV-Control", "HGPTERT-Control") #Group Comparisons
enrich_org <- "Homo sapiens" #Scientific name of the organism to perform enrichment analysis with g:Profiler (https://biit.cs.ut.ee/gprofiler/page/organism-list)

#Gene Ontology Reduced Terms settings: rrvgo package
orgDB <- "org.Hs.eg.db" #Allowed values: https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
simSimMethod <- "Rel" #Allowed values: "Resnik", "Lin", "Rel", "Jiang", "Wang"
thr <- 0.4 #Allowed values: Large (allowed similarity=0.9), Medium (0.7), Small (0.5), Tiny (0.4)
topWords <- 10 #Highest number of words to describe the members of the GO cluster during the reduced step
####################


### General Object ###
qc_plots <- list ()   #Quality Control Plots
enrichPlot <- list()  ##Enrichment Analysis Plots

reducedAllTerms <- data.frame(matrix(ncol = 9, nrow = 0)) #Reduced GOTerms Dataframe
heatmapGO <- data.frame(matrix(ncol = 5, nrow = 0)) #Heatmap Dataframe
heatmapPlots <- list() #Heatmap Plot List


### START ANALYSIS ###
#CEL files
data <- read.celfiles(list.celfiles())
n_samples <- length(list.celfiles())

#Samples Groups
targets <- limma::readTargets (sampleTarget, sep="")
ph <- data@phenoData
targets <- targets[order(targets$Filename),] 
ph@data[ ,2] <- targets$Group
colnames(ph@data)[2] <- "source"
sample <- rownames (ph@data)
ph@data <- cbind(ph@data,sample)

f <- factor(ph@data$source)
design <- model.matrix(~ 0 + f) 
colnames(design) <- levels(f)
targets <- targets [order(targets$Group),] #Sort the samples based on their group membership


### QUALITY CONTROL  MICROARRAY FILEs ###
#Normalisation using RMA
data.rma <- rma(data)
data.matrix <- exprs(data.rma)

#Raw Intensities of each array
dir.create("QC_Raw_Intensity")
for (i in 1:n_samples)
{
  name = paste("QC_Raw_Intensity/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(data[,i],main=ph@data$sample[i])
  dev.off()
}

#Chip pseudo-images based on weights
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

#Chip pseudo-images based on residuals -> 4 Types: residuals, positive residuals, negative residuals or sign of residuals

#residuals -> gives a pseudo-image of residuals
dir.create("QC_Chip_PseudoImage_residuals")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_residuals/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type="residuals", main=ph@data$sample[i])
  dev.off()
}

#pos.residuals -> only high positive residuals are drawn in red, while negative and near 0 residuals being drawn in white
dir.create("QC_Chip_PseudoImage_pos_residuals")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_pos_residuals/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type="pos.residuals", main=ph@data$sample[i])
  dev.off()
}

#neg.residuals -> only extreme negative residuals are drawn in blue, while positive negative and near 0 residuals being drawn in white
dir.create("QC_Chip_PseudoImage_neg_residuals")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_neg_residuals/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type="neg.residuals", main=ph@data$sample[i])
  dev.off()
}

#sign.residuals -> gives images where all negative residuals regardless of magnitude are indicated by blue and all positive residuals by red
dir.create("QC_Chip_PseudoImage_sign_residuals")
for (i in 1:n_samples)
{
  name = paste("QC_Chip_PseudoImage_sign_residuals/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  image(Pset, which=i, type="sign.residuals", main=ph@data$sample[i])
  dev.off()
}

#Histogram -> use the PM intensities
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

#Boxplots on the normalised data
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

#MA PLot: raw data
dir.create("QC_MAplots_rawData")
for (i in 1:n_samples)
{
  name = paste("QC_MAplots_rawData/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  MAplot(data,which=i)
  dev.off()
}


#MA PLot: normalised data
dir.create("QC_MAplots_norm")
for (i in 1:n_samples)
{
  name = paste("QC_MAplots_norm/", ph@data$sample[i], ".jpg", sep="")
  jpeg(name)
  MAplot(data.rma,which=i)
  dev.off()
}

#PCA
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



### MICROARRAY ANALYSIS ###

#Normalisation methods: rma, mas5, germa, justPlier, expresso
data.rma <- rma(data)
data.annot <- annotateEset(data.rma, hgu133plus2.db)  #ProbeID Annotation
data.fit <- limma::lmFit(data.annot,design)

#Choose the group comparisons
contrast.matrix <- limma::makeContrasts(contrasts=groupComparison, levels=design) 
data.fit.con <- limma::contrasts.fit(data.fit,contrast.matrix)
data.fit.eb <- limma::eBayes(data.fit.con)
data_final <- as.data.frame(data.fit.eb)

#Adjusting for multiple testing and defining DE genes
DEresults <- limma::decideTests (data.fit.eb, method='global', adjust.method="BH", p.value=0.000005,lfc=1)
DEresults <- rownames_to_column(as.data.frame (DEresults), "genes.PROBEID")
all_info <- merge (x = data_final, y = DEresults, by = "genes.PROBEID", all.y=TRUE)


#Save the normalised data
norm <- as.data.frame(data.rma)
tnorm <- as.data.frame (t(norm))
tnorm <- rownames_to_column(tnorm, "PROBEID")
annot <- data.annot@featureData@data
tnorm_annot <- merge(x=tnorm, y=annot, by="PROBEID", all.y=TRUE)


#Create the right DF frame from the "all_info" dataframe
colnames(all_info) <- make.names(colnames(all_info))
all_info$DEG <- rowSums(abs(all_info[(ncol(all_info)-ncol(contrast.matrix)+1) : ncol(all_info)]))

resultDF <- all_info %>% filter(genes.ENTREZID != "NA") %>%    #Remove rows where the ENTREZID is "NA"
  filter(DEG > 0) %>%                   #Filter for rows where the gene is DE
  dplyr::select(starts_with(c("genes", "p.value", "coefficients")), make.names(colnames(contrast.matrix))) #Select columns

#Get the unique list of DEGs for each group
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

#Create Clustering Plot
uniqueDEG <- unique(uniqueDEG)
degClusterDF <- data.frame(matrix(0, ncol=length(colnames(contrast.matrix)), nrow=length(uniqueDEG)))
colnames(degClusterDF) <- make.names(colnames(contrast.matrix))
rownames(degClusterDF) <- uniqueDEG

for (i in names(deg_unique_group)) 
{for (l in deg_unique_group[[i]]) degClusterDF[l, i] <- 1}

dist_mat <- dist(t(degClusterDF), method = 'euclidean')


### Enrichment Analysis -> g:Profiler ###
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
  
  ### Loop through the GO terms of the enriched analysis terms found ###
  for (iGO in c("GO:BP", "GO:CC", "GO:MF"))
  {
    ### Filter the GO terms for their type ###
    filterGO <- enrich_go_kegg %>% filter (source==iGO)
    
    ### GO terms Semantic Similarity calculation ###
    simMatrix <- calculateSimMatrix (filterGO$term_id, orgdb=orgDB, ont=strsplit(iGO, split = "GO:")[[1]][2], method=simSimMethod)
    scores <- setNames (-log10(filterGO$p_value), filterGO$term_id)
    reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=thr, orgdb=orgDB)
    reducedTerms$goType <- iGO
    reducedAllTerms <- rbind(reducedAllTerms, reducedTerms)
    
    ### Loop through the Clusters ###
    for (iCluster in unique(reducedTerms$cluster)) 
    {
      #Cluster
      clusterX <- reducedTerms %>% 
        filter (cluster==iCluster) %>% #Get all the members/rows for the present cluster.
        arrange(desc(score), desc(parentSimScore))    #Sort by the "score" and by the "parentSimScore" columns.
      
      clusterXParent <- clusterX[1,]   #Get the first row which has the max value for both the "score" and the "parentSimScore" column. 
      clusterXMembers <- clusterX[-1,]  #Get all the row/s (non-parent term) except the first one.
      
      #Description for terms which are non-parent
      terms_members <- paste(clusterXMembers$term, collapse = ' ') #Get all the terms from the members non-parent as a string
      terms_members <- strsplit (terms_members, split=" ")[[1]] #Get all the terms from the members non-parent as a vector
      
      terms_members <- terms_members [!terms_members %in% strsplit(clusterXParent$term, split=" ") [[1]]] #Remove the parent terms from the members non-parent terms description
      terms_members <- terms_members [!(terms_members %in% stop_words$word)] #Remove un-informative words such as: of, via, and, etc..
      
      terms_members <- as.data.frame(table (terms_members)) %>% arrange(desc (Freq)) #Order the words of the members based on their frequency
      terms_members <- paste(terms_members$terms_members[1:min(c(topWords, nrow(terms_members)))], collapse = " ") #Get the top words with the higher frequency
      
      #Number of GO terms present in the cluster (with the parent term)
      clusterX_nodes <- nrow(clusterX)
      
      ### Loop through the Comparisons: Get the lowest p-value for each comparison for the present cluster ###
      for (iComparison in unique (filterGO$query)) 
      {
        filterComparison <- filterGO %>% 
          filter(query==iComparison) %>%    #Filter for rows of the present comparison
          filter(term_id %in% clusterX$go) #Filter for the GO:terms present in the cluster
        
        #Insert values in the heatmap dataframe
        des <- paste (sprintf("Cluster %s", iCluster), clusterXParent$term, terms_members, clusterX_nodes, sep = " * ")
        
        if(nrow(filterComparison)>0) {minPvalue <- min(filterComparison$p_value)}
        if(nrow(filterComparison)==0) {minPvalue <- NA}
        
        heatmapGO <- rbind(heatmapGO, c(iCluster, des, iComparison, minPvalue, iGO))
      }
    }
  }
  
  colnames(heatmapGO) <- c("Cluster", "Description", "Comparison", "Pvalue", "Type")
  
  
  ##### KEGG terms #####
  filterKEGG <- filter(enrich_go_kegg, source=="KEGG") #Filter for KEGG terms
  
  heatmapKEGG <- data.frame(matrix(ncol = length (unique(filterKEGG$query)), nrow = length (unique (filterKEGG$term_id))))
  colnames(heatmapKEGG) <- unique(filterKEGG$query)
  rownames(heatmapKEGG) <- unique (filterKEGG$term_name)
  
  #Save info in the heatmap df
  for (i in 1:nrow(filterKEGG)) 
  {heatmapKEGG[filterKEGG$term_name[i], filterKEGG$query[i]] <- filterKEGG$p_value[i]}
  
  heatmapKEGG <- heatmapKEGG %>% rownames_to_column(var="Description")
  heatmapKEGG <- reshape2::melt(heatmapKEGG, id.vars="Description")
  heatmapKEGG$Type <- "KEGG"
  heatmapKEGG$Cluster <- 1
  
  keggClusters <- setNames (c(1:length(unique(heatmapKEGG$Description))), unique(heatmapKEGG$Description))
  
  for (ikeggClusters in keggClusters) 
  {
    iDes <- names(keggClusters) [keggClusters == ikeggClusters]
    heatmapKEGG <- heatmapKEGG %>% mutate(Cluster=replace(Cluster, Description==iDes, ikeggClusters)) %>% as.data.frame()
  }
  
  colnames(heatmapKEGG) <- c("Description", "Comparison", "Pvalue", "Type", "Cluster")
  
  
  ### Combine GO and KEGG terms heatmap dataframe ###
  heatmapGOKEGG <- rbind(heatmapGO, heatmapKEGG)
  heatmapGOKEGG$Pvalue <- as.numeric(as.character(heatmapGOKEGG$Pvalue))
  
  
  ### Create Heatmap plots ###
  for (iType in unique(heatmapGOKEGG$Type)) 
  {
    filterGOKEGG <- heatmapGOKEGG %>% filter(Type==iType) #Filter for GO and KEGG type
    
    #Get the order 
    filterGOKEGG$order <- ""
    
    for (iCluster in unique(filterGOKEGG$Cluster)) 
    {
      ffCluster <- filterGOKEGG %>% filter(Cluster==iCluster) %>% filter(Pvalue>0)
      iOrder <- paste(sort (ffCluster$Comparison), collapse = " * ")
      filterGOKEGG <- filterGOKEGG %>% mutate(order=replace(order, Cluster==iCluster, iOrder)) %>% as.data.frame()
    }
    
    filterGOKEGG$count <- str_count(filterGOKEGG$order, "\\*") + 1
    filterGOKEGG <- arrange(filterGOKEGG, count, order)
    filterGOKEGG$Description <- factor(filterGOKEGG$Description, levels = unique(filterGOKEGG$Description))
    
    #Save the heatmap plot
    heatmapPlots[[sprintf ("%s_%s", orgDB, iType)]] <- filterGOKEGG %>% 
      ggplot(aes(x=Comparison, y=Description, fill=Pvalue)) +
      geom_raster(aes(fill = Pvalue)) +
      scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
      geom_text(size=3, aes(label=scientific(Pvalue, digits = 3))) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
      theme(axis.text.y = element_text(size=10, face="bold")) +
      ggtitle(sprintf ("%s_%s", orgDB, iType))
  }
}


### Save Dataframes ###
write.csv(all_info, file = "all_info.csv")                    #All the info for the probes and their DEGs analysis
write.csv(tnorm, file = "rma_norm.csv")                       #RMA normalisation data
write.csv(resultDF, file = "deg_analysis_duplicates.csv")     #DEG analysis with all the gene's duplicates due to the different probes for the same gene 
write.csv(deg_unique_group_DF, file = "deg_unique_group.csv") #DEG list for each group comparison
write.csv(enrich_go_kegg, file = "gProfiler_enrich.csv")      #Enrichment Analysis
write.csv(reducedAllTerms, file = "goReducedTerms.csv")       #Reduced GO terms
write.csv(heatmapGOKEGG, file = "heatmapDF.csv")              #Heatmap plot dataframe

### Save the Cluster Dendogram plot ###
pdf("dendogram.pdf") 
plot(hclust(dist_mat, method = 'average'))
dev.off() 

##### Save Quality Control Plots #####
pdf("quality_control_plots.pdf", onefile=TRUE)
for (i in seq(length(qc_plots))) {grid.arrange (qc_plots[[i]])}
dev.off()

##### Save Enrichment Analysis Plots #####
pdf("gProfiler_reducedTerms_heatmap_plots.pdf", onefile=TRUE)
for (i in seq(length(heatmapPlots))) {grid.arrange (heatmapPlots[[i]])}
dev.off()

### Save the all the plots in an R data ###
allPlots <- c(qc_plots, heatmapPlots)
save(allPlots, file = "QC_heatmap_plots.RData")
