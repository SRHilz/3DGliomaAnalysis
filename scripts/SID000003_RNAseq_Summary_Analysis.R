## Costello Lab
## 2018.11.01
## Stephanie Hilz

# Load required libraries
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR",ask=FALSE, suppressUpdates=TRUE)
biocLite("statmod",ask=FALSE, suppressUpdates=TRUE)
biocLite("mygene",ask=FALSE, suppressUpdates=TRUE)
biocLite("org.Hs.eg.db",ask=FALSE, suppressUpdates=TRUE)
#biocLite("goseq", ask=FALSE, suppressUpdates=TRUE)
library(goseq)
library(org.Hs.eg.db)
if("dplyr" %in% rownames(installed.packages()) == FALSE) {
  install.packages("dplyr", repos="http://cran.rstudio.com/")
}
if("vegan" %in% rownames(installed.packages()) == FALSE) {
  install.packages("vegan", repos="http://cran.rstudio.com/")
}
if("pheatmap" %in% rownames(installed.packages()) == FALSE) {
  install.packages("pheatmap", repos="http://cran.rstudio.com/")
}

annotateGO <- function(GO.wall,invertedGetgo){
  GO.wall$foreground_genes <- 'NA'
  for (cat in GO.wall$category){
    GO.wall[which(GO.wall$category==cat),]$foreground_genes <- paste(invertedGetgo[[cat]],collapse=';')
  }
  return(GO.wall)
}

invertGetgo <- function(getgolist){
  inversion <- list()
  getgolist <- getgolist[!is.na(names(getgolist))]
  for (gene in names(getgolist)){
    for (goterm in getgolist[[gene]]){
      if (!goterm %in% names(inversion)){
        inversion[[goterm]] = c()
      }
      inversion[[goterm]] <- append(inversion[[goterm]], gene)
    }
  }
  return(inversion)
}

performGO <- function(binaryList, outfile){
  print("Table of input values")
  print(table(binaryList))
  pwf=nullp(binaryList,'hg19',"geneSymbol")
  GO.wall=goseq(pwf,"hg19","geneSymbol")
  print("Top 20 most significant GO terms")
  top <- GO.wall[,c(6,2)]
  colnames(top) <- c("term","pvalue")
  rownames(top) <- NULL
  print(head(top,20))
  getgolist <- getgo(names(binaryList[which(binaryList==1)]), 'hg19','geneSymbol')
  getGeneList <- invertGetgo(getgolist)
  GO.wall.anno <- annotateGO(GO.wall,getGeneList)
  write.csv(GO.wall.anno,outfile, quote=F)
}

library(edgeR)
library(statmod)
library(dplyr)
library(vegan)
library(pheatmap)
library(goseq)
library(geneLenDataBase)
library(org.Hs.eg.db)
library(mygene)

#Functions
ensmblRoot <- function(ensmblID){
  return(strsplit(as.character(ensmblID),'[.]')[[1]][1])
}

# User-defined variables
tag <- 'SID000003_20210524_resubmission'
outfolder <- 'SID000003_RNAseq_Summary_Analysis/'
countsType <- '.symbol.coding.counts'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# set file info
countsInfile <- paste0(dataPath,gsub('SID000003_','',tag),countsType)
CPMOutputFile <- paste0(outputPath,outfolder,tag,gsub('.counts','.CPMs.csv', countsType)) #will move CPM files from output to data folder after checking
PCAFile <- paste0(outputPath,outfolder,tag,'_PCA.pdf')
PCALoadingPath <- paste0(outputPath,outfolder)
PCALoadingGOPath <- paste0(outputPath,outfolder)

# Create counts object with cpms
par(mfrow=c(1,1),mar=rep(4,4))
sampleTable_edgeR<-read.delim(countsInfile,sep='\t')#the data
rownames(sampleTable_edgeR) <- sampleTable_edgeR$symbol
sampleTable_edgeR$symbol <- NULL
y<-DGEList(counts=sampleTable_edgeR)#making a DGEList
dim(y)
keep <- rowSums(cpm(y)>1) >= 1#filters out anything without CPM>N in at least 1 library
y <- y[keep,,keep.lib.sizes=FALSE]
dim(y)
boxplot(log(cpm(y),2))
y<-calcNormFactors(y)#default normalization
boxplot(log(cpm(y),2))
write.csv(cpm(y),CPMOutputFile)
write.table(cpm(y),gsub('.csv','.txt',CPMOutputFile), sep='\t', quote=F, col.names=NA)

## Get sample metadata
# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)

# merge by patient ID
merged <- merge(data, subtypedata, by="Patient")

# color assignment by tumor type
merged$color <- 'black'
merged[merged$IDH_Mut==0,]$color <- 'red'
merged[merged$IDH_Mut==0 & merged$TERT==0,]$color <- 'orange'
merged$sampleID <- paste0(merged$Patient,merged$SampleName)

# color assignment by purity
merged$purity <- merged$PyClone
merged[which(merged$PurityEstUsed=='FACETS'),]$purity <- merged[which(merged$PurityEstUsed=='FACETS'),]$FACETS
merged$color2 <- 'orange'
merged[which(merged$purity < 0.3),]$color2 <- 'black'
merged[is.na(merged$purity),]$color2 <- 'grey'

#PubQuality PCA plot with % explained
par(bg = 'white')
par(mfrow=c(1,1))
par(mar=c(4,4,4,4))
data=read.table(CPMOutputFile,header = TRUE, sep = ',', stringsAsFactors=F)
orderedGeneNames <- data$X
data$X <- NULL
nrm_count_matrix=as.matrix(data)
rownames(nrm_count_matrix) <- orderedGeneNames
log=log10(1+nrm_count_matrix)
tlog=t(log)
pca=prcomp(tlog)
summary(pca)
colors <- merged$color2
names(colors) <- merged$sampleID
colors = colors[colnames(nrm_count_matrix)]
raw <- pca$x[,1:4]
xaxis <- 1
yaxis <- 2
xlab <- paste('PC',toString(xaxis),' (',toString(round(100*summary(pca)$importance[2,xaxis])),'%)')
ylab <- paste('PC',toString(yaxis),' (',toString(round(100*summary(pca)$importance[2,yaxis])),'%)')
xaxisExtra <-(max(raw[,xaxis])-min(raw[,xaxis]))/4
yaxisExtra <-(max(raw[,yaxis])-min(raw[,yaxis]))/4
pdf(PCAFile)
plot(raw[,xaxis], raw[,yaxis],  col=colors, pch=20, 
     xlim=c(min(raw[,xaxis])-xaxisExtra,max(raw[,xaxis])+xaxisExtra),
     ylim=c(min(raw[,yaxis])-yaxisExtra,max(raw[,yaxis])+yaxisExtra),
     xlab=xlab,
     ylab=ylab,
     cex=1,
     cex.lab=1,
     cex.axis=1)
text(raw[,xaxis], raw[,yaxis], names(colors), c(1,1), cex=.5)
dev.off()

# Get top 100 genes for a PC and do gene ontology
PCs <- c("PC1","PC2","PC3","PC4")
for (PC in PCs){
  loading <- names(sort(abs(pca$rotation[,PC]), decreasing=T)[1:100])
  if (PC=="PC1"){
    masterLoading <- as.data.frame(loading)
    masterLoading$new <- TRUE
    colnames(masterLoading) <- c("gene",PC)
  } else {
    toAdd <- loading
    toAdd <- as.data.frame(toAdd)
    toAdd$new <- TRUE
    colnames(toAdd) <- c("gene",PC)
    masterLoading <- merge(masterLoading, toAdd, by="gene", all=TRUE)
  }
  forgo <- rep(0,length(rownames(nrm_count_matrix)))
  names(forgo) <- rownames(nrm_count_matrix)
  forgo[loading] <- 1
  print(table(forgo))
  PCALoadingGOFile <- paste0(PCALoadingGOPath,tag,'PCALoadingGO_',PC,'.csv')
  performGO(forgo, PCALoadingGOFile)
  PCALoadingFile <- paste0(PCALoadingPath,tag,'PCALoadingList_',PC,'.csv')
  write.csv(loading,PCALoadingFile)
}
rownames(masterLoading) <- masterLoading$gene
masterLoading$gene <- NULL
masterLoading[is.na(masterLoading)] <- FALSE
masterLoadingMatrix <- as.matrix(apply(masterLoading, 2, as.numeric))
rownames(masterLoadingMatrix) <- rownames(masterLoading)