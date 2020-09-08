## Costello Lab
## 2019.09.13
## Stephanie Hilz
## Usage: Takes in CPMs, study config files, and cell state gene enrichment lists and performs an array of analyses 
##  to look at how cell states are distributed in glioma

computeZscore <- function(cpmList){
  zscores <- c()
  cpmList <- log(1+cpmList,2)
  mean <- mean(cpmList)
  sd <- sd(cpmList)
  for (i in cpmList){
    z <- (i - mean)/sd
    zscores <- append(zscores, z)
  }
  names(zscores) <- names(cpmList)
  return(zscores)
}

pairwiseDist2pointsLPS <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + (((p2[3]-p1[3]))^2) ))
}

computeMeanZscore <- function(cpmMatrix){
  zscoreMatrix <- cpmMatrix
  for (i in 1:nrow(cpmMatrix)){
    unlisted <- cpmMatrix[i,] %>% unlist
    zscoreMatrix[i,] <- computeZscore(unlisted)
  }
  zScoreMean <- colMeans(zscoreMatrix)
  return(zScoreMean)
}

library(gtools)
library(gplots)
library(ggpmisc)
library(kableExtra)
library(dplyr)

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# User-defined variables and file paths (will also need to make sure VenteicherSignatures folder is in dataPath)
tag <- 'SID000017'
outfolder <- 'SID000017_generate_cancerSEA_dataset/'
CPMFile <- paste0(dataPath,'SID000003_20190913_expanded_gbm.symbol.coding.CPMs.csv')
signaturesPath <- paste0(dataPath,'CancerSEASignatures/')

## read in cpm data for all patients
CPMs <- read.csv(CPMFile)
rownames(CPMs) <- CPMs$X
CPMs$X <- NULL
sampleN <- ncol(CPMs)
sampleNames <- colnames(CPMs)
my_palette <- colorRampPalette(c("black", "pink", "red"))(n = 1000)

# read in relevant metadata
sampleData <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)
sampleData$sampleID <- paste0(sampleData$Patient,sampleData$SampleName)
sampleData <- sampleData[which(sampleData$sampleID %in% sampleNames),]

# restrict sampleNames to those also in sampleData
sampleNames <- sampleNames[sampleNames %in% sampleData$sampleID]

## read in signatures and remove genes for which we don't have CPMs for
signatureFiles <- list.files(path=signaturesPath, pattern="*.txt", full.names=FALSE, recursive=FALSE)
signaturesAll <- list()
for (file in signatureFiles){
  signatureName <- gsub('.txt','',file)
  print(signatureName)
  signature <- read.table(paste0(signaturesPath, file), header = TRUE, 
                                colClasses=c("NULL","character")) %>% unlist
  colnames(signature) <- NULL
  print(length(signature))
  signaturesAll[[signatureName]] <- signature[signature %in% rownames(CPMs)] 
  print(length(signaturesAll[[signatureName]]))
}

# set up structures to collect data in
signatureNames <- names(signaturesAll)
zscores <- setNames(replicate(length(signatureNames),data.frame()),signatureNames)
zscoresdf <- c()
geneExpression <- setNames(replicate(length(signatureNames),data.frame()),signatureNames)
patientLables <- sampleData$Patient
geneCorrelation <- data.frame(correlation = numeric(), signatureMeanGeneExpression= numeric(), signature = character())
# populate data structures
for (sig in signatureNames){
  # array of all genes in signature for which we have CPMs
  signatureGenes <- signaturesAll[sig] %>% unlist 
  signatureGeneExpression <- CPMs[signatureGenes,sampleNames]
  # populate our object that collects CPMs for each gene in signature
  geneExpression[[sig]] <- signatureGeneExpression
  # populate our object that collcts scores for each signature
  zscore <- computeMeanZscore(signatureGeneExpression)
  # determine the mean expression of each gene in that signature
  signatureMeanGeneExpression <- rowMeans(signatureGeneExpression)
  # determine the correlation between each gene's expression and the zscore
  signatureGeneCorrelation <- cor(t(signatureGeneExpression), zscore)
  #populate object that collects this information for correlation and expression (by gene)
  toBind <- cbind(gene = rownames(signatureGeneCorrelation), correlation = signatureGeneCorrelation, meanCPM = data.frame(signatureMeanGeneExpression), signature = sig)
  row.names(toBind) <- c()
  geneCorrelation <-  rbind(geneCorrelation, toBind)
  #populate object that collects this information for the zscore (by signature)
  means <- colMeans(signatureGeneExpression)
  scoreOut <- cbind(patientLables,sampleNames,zscore, means)
  rownames(scoreOut) <- NULL
  zscores[[sig]] <- rbind(zscores[[sig]], scoreOut)
  colnames(zscores[[sig]]) <- c('patient','sampleID','zscore', 'means')
  zscores[[sig]]$patient <- as.character(zscores[[sig]]$patient)
  zscores[[sig]]$sampleID <- as.character(zscores[[sig]]$sampleID)
  zscores[[sig]]$zscore <- as.numeric(as.character(zscores[[sig]]$zscore))
  zscores[[sig]]$means <- as.numeric(as.character(zscores[[sig]]$means))
  zscoresdf <- cbind(zscoresdf, zscore)
}
colnames(zscoresdf) <- signatureNames
write.table(zscoresdf, file=paste0(outputPath,outfolder,'cancerSEA_signature_zscores.txt'), sep='\t', col.names = NA)

# plot correlation vs mean expression for each gene in each signature (used to identify genes that server as markers for a signature)
ggplot(geneCorrelation, aes(x=correlation, y=log(signatureMeanGeneExpression,2))) +
  geom_point(size=.8) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
  labs(y='log2 Mean CPM', x='Correlation (Pearson\'s R)') +
  facet_wrap(~signature)

# compare abundance between tumor types
toPlot <- bind_rows(zscores, .id = "signature")
toPlotAgg <- aggregate(zscore~patient+signature, toPlot, mean)
toPlotAgg$subtype <- NA
toPlotAgg[which(toPlotAgg$patient %in% subsetIDHwtc),]$subtype <- 'IDH-wtc'
toPlotAgg[which(toPlotAgg$patient %in% subsetIDHmut),]$subtype <- 'IDH-mut'
toPlotAgg <- toPlotAgg[which(!is.na(toPlotAgg$subtype)),]
toPlotAgg$subtype <- factor(toPlotAgg$subtype, levels=unique(toPlotAgg$subtype))
ggplot(toPlotAgg, aes(x=signature, y=zscore, fill=subtype)) +
  geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c("black","#11cc42"))+
  labs(y= "Mean gene z-score of log2 CPM")+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))

# add in distance info
zscoresdfWithSpatial <- as.data.frame(zscoresdf)
zscoresdfWithSpatial$sampleID <- rownames(zscoresdfWithSpatial)
zscoresdfWithSpatial <- merge(zscoresdfWithSpatial, sampleData, by="sampleID")
# Exploration of gene expression
sig <- 'Stemness'
m <- as.matrix(zscoresdf)
m <- as.matrix(zscoresdfWithSpatial[which(!is.na(zscoresdfWithSpatial$DistCentroid)),c(signatureNames,'DistCentroid','DistPeriph')])
dist.pear <- function(x) as.dist(1-cor(t(x)))
pdf(paste0(outputPath, outfolder,'allsigs_heatmap_with_distances.pdf'), width=20, height=40)
h <- heatmap.2(m, 
               trace="none", 
               dendrogram="both",
               distfun=dist.pear,
               density.info="none",
               cexCol=1,
               scale='row',
               col = my_palette)
dev.off()

# Exploration of z scores
sig <- 'Angiogenesis'
toPlot <- zscores[[sig]][order(zscores[[sig]]$zscore),]
barplot(toPlot$zscore)
boxplot(toPlot$zscore~toPlot$patient, ylab=paste0(sig, " signature Zscore"))

# Exploration of mean CPMs
sig <- 'Angiogenesis'
toPlot <- colMeans(geneExpression[[sig]])
barplot(toPlot)
## add in distance info to samples and their zscores
mergedSubset <- merged[,c('SampleName','L.Coordinate','P.Coordinate','S.Coordinate','DistCentroid','DistPeriph','Histology','Tumor','purity','sampleID')]

# Create a data structure that also includes distances (right now just have for Astros)
zscoresWithDistances <- setNames(replicate(length(signatureNames),data.frame()),signatureNames)
for (s in signatureNames){
  zscoresWithDistances[[s]] <- merge(zscores[[s]],mergedSubset, by="sampleID")
  zscoresWithDistances[[s]]$peripheralCat <- 0
  zscoresWithDistances[[s]][which(zscoresWithDistances[[s]]$DistPeriph <= 7),]$peripheralCat <- 1
  zscoresWithDistances[[s]]$peripheralCat <- factor(zscoresWithDistances[[s]]$peripheralCat, levels=c(0,1))
  zscoresWithDistances[[s]]$centroidCat <- 0
  zscoresWithDistances[[s]][which(zscoresWithDistances[[s]]$DistCentroid <= 17),]$centroidCat <- 1
  zscoresWithDistances[[s]]$centroidCat <- factor(zscoresWithDistances[[s]]$centroidCat, levels=c(0,1))
  zscoresWithDistances[[s]]$patient <- factor(zscoresWithDistances[[s]]$patient,levels=patientsToUse)
}

## Create zscoresWithDistance that does not contain missing Centroid/Periph data
zscoresWithDistancesCentPeriph <- zscoresWithDistances[[sig]][!is.na(zscoresWithDistances[[sig]]$DistCentroid),]

# plot by correlation (change x to go between peripheral and core, and sig to go between signatures)
colors <- as.character(colorKey[patientsToUse[patientsToUse %in% unique(zscoresWithDistancesCentPeriph$patient)]])
ggplot(zscoresWithDistancesCentPeriph, aes(x=DistCentroid, y=zscore, color=patient)) +
  geom_point() +
  scale_colour_manual(values=colors) +
  theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_smooth(aes(colour=factor(patient)), method = "lm", se=F)+
  labs(y=paste0(sig, ' z-score')) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x.npc = "right", label.y.npc = 0.8, formula = y~x, 
               parse = TRUE, size = 3) + 
  stat_fit_glance(method = 'lm',
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.4, size = 3)

## Plot comparison of score per sample by patient by peripheral or centroid
ggplot(zscoresWithDistancesCentPeriph, aes(x = patient, y = zscore, fill=centroidCat)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Dark2") +
  labs(y=paste0(sig, ' z-score')) +
  theme(axis.text.x = element_text(size=10, angle = 90, hjust = 1), axis.title = element_text(size = 10), axis.text.y = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'black'))

## Calculate spearmans
patients <- patientsToUse[patientsToUse %in% unique(zscoresWithDistancesCentPeriph$patient)]
colors <- colors
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 20 #where on the x axis will display p value
y=max(toReturnSignature$difference)/1.5
x.var = "DistPeriph"
y.var = "zscore"
y.inc = y/15
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- zscoresWithDistancesCentPeriph[which(zscoresWithDistancesCentPeriph$patient == patientID),]
  testResult <- cor.test(patientSubset[,x.var], patientSubset[,y.var], method="spearman")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  rho=round(testResult$estimate,3)
  label <- paste0('p=',p,', rho=',rho)
  dataText <- rbind(dataText, c(p,rho,label,x,y,color,patientID), stringsAsFactors=F)
  y <- y + y.inc
}
colnames(dataText) <- c('p','rho','label','x','y','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'_',x.var,'_vs_',sig,'_zscore_difference_stats.txt'),sep='\t', quote=F, row.names = F)

## Finally, look at relationship between pairwise distances and score differences
# create a data structure that also contains distance to do some correlations with
toReturn <- data.frame(signature=character(), patient=character(), samplePair=character(), difference=numeric(), distance=numeric(), stringsAsFactors=F)
for (p in unique(zscoresWithDistances[[sig]]$patient)){
  zscoresWithDistancesLocalP <- zscoresWithDistances[[sig]][which(zscoresWithDistances[[sig]]$patient == p),]
  print(p)
  comparisons <- combn(zscoresWithDistancesLocalP$sampleID,2)
  p <- as.character(p)
  print(dim(comparisons))
  # Loop through and calculate euclidian spatial distance as well as signature difference for each pair
  for (c in seq_len(dim(comparisons)[2])){
    signature <- sig
    sample_a <- comparisons[1,c]
    sample_b <- comparisons[2,c]
    samplePair <- paste0(sample_a, '-', sample_b)
    difference <- abs(zscoresWithDistancesLocalP[which(zscoresWithDistancesLocalP$sampleID == sample_a),]$zscore - zscoresWithDistancesLocalP[which(zscoresWithDistancesLocalP$sampleID == sample_b),]$zscore)
    Sample1LPS <- zscoresWithDistancesLocalP[which(zscoresWithDistancesLocalP$sampleID == sample_a),c("L.Coordinate","P.Coordinate","S.Coordinate")]
    Sample2LPS <- zscoresWithDistancesLocalP[which(zscoresWithDistancesLocalP$sampleID == sample_b),c("L.Coordinate","P.Coordinate","S.Coordinate")]
    distance <- as.numeric(pairwiseDist2pointsLPS(Sample1LPS,Sample2LPS))
    toReturn <- rbind(toReturn, data.frame(signature=signature, patient=p, samplePair=samplePair, difference=difference, distance=distance, stringsAsFactors=F))
  }
}

colors <- as.character(colorKey[patientsToUse[patientsToUse %in% unique(toReturn$patient)]])
ggplot(toReturn, aes(x=distance, y=difference, color=patient)) +
  geom_point() +
  scale_colour_manual(values=colors) +
  theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_smooth(aes(colour=factor(patient)), method = "lm", se=F)+
  labs(y=paste0(sig, ' z-score pairwise difference')) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x.npc = "right", label.y.npc = 1.0, formula = y~x, 
               parse = TRUE, size = 3) + 
  stat_fit_glance(method = 'lm',
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.4, size = 3)


## Calculate spearmans
df <- toReturn
patients <- patientsToUse[patientsToUse %in% unique(df$patient)]
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 20 #where on the x axis will display p value
x.var = "distance"
y.var = "difference"
y=max(df[,y.var])/1.5
y.inc = y/15
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- df[which(df$patient == patientID),]
  testResult <- cor.test(patientSubset[,x.var], patientSubset[,y.var], method="spearman")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  rho=round(testResult$estimate,3)
  label <- paste0('p=',p,', rho=',rho)
  dataText <- rbind(dataText, c(p,rho,label,x,y,color,patientID), stringsAsFactors=F)
  y <- y + y.inc
}
colnames(dataText) <- c('p','rho','label','x','y','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'_spatial_distance_vs_',sig,'_zscore_difference_stats.txt'),sep='\t', quote=F, row.names = F)
