## Costello Lab
## 2018.11.01
## Stephanie Hilz

library(statmod)
library(dplyr)
library(vegan)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(ggplot2)

#Functions
pairwiseDist2pointsLPS <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + (((p2[3]-p1[3]))^2) ))
}

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# User-defined variables
tag <- 'SID000004'
outfolder <- 'SID000004_rnaseq_expression_distances_intra_inter/'
CPMFile <- paste0(dataPath,'SID000003_20190913_expanded_gbm.symbol.coding.CPMs.csv')
purityCutoff <- .7 #only used for subset of analyses

## Get sample metadata
# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T, stringsAsFactors = F)

# create color key by patient from subtype folder
colorKey <- subtypedata$Color
names(colorKey) <- subtypedata$Patient

# merge by patient ID
merged <- merge(data, subtypedata, by="Patient")

# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF
print("Warning: forcing purity NA to .1 - this is not appropriate if an NA indicates no estimate rather than a low FALCON-X value")
merged[is.na(merged$purity),]$purity <- .1

# assign sampleID
merged$sampleID <- paste0(merged$Patient,merged$SampleName)

# subset of merge
mergedSubset <- merged[,c('Patient','SampleName','L.Coordinate','P.Coordinate','S.Coordinate','DistCentroid','DistPeriph','IDH_Mut','Tumor','purity')]
mergedSubset$sampleID <- paste0(mergedSubset$Patient,mergedSubset$SampleName)

# Creation of a distance matrix
data=read.table(CPMFile,header = TRUE, sep = ',', stringsAsFactors=F)
orderedGeneNames <- data$X
data$X <- NULL
nrm_count_matrix=as.matrix(data)
rownames(nrm_count_matrix) <- orderedGeneNames
loggedNrmCountMatrix <- log(nrm_count_matrix+ 0.01, 2)
correlation <- cor(loggedNrmCountMatrix, method="pearson")
dissimilarity <- (1 - correlation)/2

# Separation in to inter intra pairs
dissimilarityToMelt <- dissimilarity
dissimilarityToMelt[lower.tri(dissimilarityToMelt, diag=T)] <- NA #do this to ensure do not get duplicate comparisons or self-self comparisons
dissimilarityPairs <- melt(dissimilarityToMelt)
colnames(dissimilarityPairs) <- c('sample1','sample2','dissimilarity')
dissimilarityPairs <- dissimilarityPairs[!is.na(dissimilarityPairs$dissimilarity),]

# Add in patient info
mergedSubset1 <- mergedSubset
colnames(mergedSubset1)[which(colnames(mergedSubset1) == 'sampleID')] <- 'sample'
colnames(mergedSubset1) <- paste0(colnames(mergedSubset1), '1')
toAnalyze <- merge(dissimilarityPairs, mergedSubset1, by='sample1')
colnames(toAnalyze)[which(colnames(toAnalyze) == 'Patient')] <- 'Patient1'
mergedSubset2 <- mergedSubset
colnames(mergedSubset2)[which(colnames(mergedSubset2) == 'sampleID')] <- 'sample'
colnames(mergedSubset2) <- paste0(colnames(mergedSubset2), '2')
toAnalyze <- merge(toAnalyze, mergedSubset2, by='sample2')

# add in columns for purity difference and purity fold difference
toAnalyze$purityDifference <- abs(toAnalyze$purity1-toAnalyze$purity2)
toAnalyze$purityFold <- toAnalyze$purity1/toAnalyze$purity2
toAnalyze[which(toAnalyze$purityFold < 1),]$purityFold <- 1/toAnalyze[which(toAnalyze$purityFold < 1),]$purityFold

# patients to drop if any
#drop <- c("P260","P452")
#toAnalyze <- toAnalyze[-which(toAnalyze$Patient1 %in% drop | toAnalyze$Patient2 %in% drop),]

# Designate inter intra
toAnalyze$relationship <- 'inter'
toAnalyze[which(toAnalyze$Patient1 == toAnalyze$Patient2),]$relationship <- 'intra'

# Designate inter-intra patient combinations
toAnalyze$patientCombination <- paste0(toAnalyze$Patient1,'-',toAnalyze$Patient2)

# Designate inter molecular combination categories
toAnalyze$interType <- NA
toAnalyze[which(toAnalyze$relationship == 'inter'),]$interType <- paste0(toAnalyze[which(toAnalyze$relationship == 'inter'),]$IDH_Mut1, '-', toAnalyze[which(toAnalyze$relationship == 'inter'),]$IDH_Mut2)

# Save file with data
write.table(toAnalyze, file=paste0(outputPath,outfolder, 'SID000004_pairwise_metrics_expression_dissimilarity.txt'), sep = '\t', row.names = F, quote=F)

# Calculate means for each intra- or inter-patient comparison
toAnalyzeMeans <- data.frame(Patient1=character(),
                             Patient2=character(),
                             IDH_Mut1=character(),
                             IDH_Mut2=character(),
                             dissimilarityMean=numeric(),
                             purityMean1=numeric(),
                             purityMean2=numeric(),
                             purityDifference=numeric(),
                             purityFold=numeric())
for (c in unique(toAnalyze$patientCombination)){
  print(c)
  toAnalyzeC <-  toAnalyze[which(toAnalyze$patientCombination==c),]
  print(dim(toAnalyzeC))
  Patient1 <- toAnalyzeC$Patient1[1]
  Patient2 <- toAnalyzeC$Patient2[1]
  IDH_Mut1 <- toAnalyzeC$IDH_Mut1[1]
  IDH_Mut2 <- toAnalyzeC$IDH_Mut2[1]
  purityDifference <- mean(toAnalyzeC$purityDifference)
  purityFold <- mean(toAnalyzeC$purityFold)
  dissimilarityMean <- mean(toAnalyzeC$dissimilarity)
  purityMean1 <- mean(mergedSubset[mergedSubset$Patient==Patient1,]$purity)
  purityMean2 <- mean(mergedSubset[mergedSubset$Patient==Patient2,]$purity)
  toAnalyzeMeans <- rbind(toAnalyzeMeans,data.frame(Patient1=Patient1,
                                     Patient2=Patient2,
                                     IDH_Mut1=IDH_Mut1,
                                     IDH_Mut2=IDH_Mut2,
                                     dissimilarityMean=dissimilarityMean,
                                     purityMean1=purityMean1,
                                     purityMean2=purityMean2,
                                     purityDifference=purityDifference,
                                     purityFold=purityFold
                                     ))
}
toAnalyzeMeans$patientRelationship <- "inter"
toAnalyzeMeans[which(toAnalyzeMeans$Patient1 == toAnalyzeMeans$Patient2),]$patientRelationship <- 'intra'
toAnalyzeMeans$subtypeRelationship <- "between"
toAnalyzeMeans[which(toAnalyzeMeans$IDH_Mut1 == toAnalyzeMeans$IDH_Mut2),]$subtypeRelationship <- 'within'
toAnalyzeMeans$IDH_Mut1 <- factor(toAnalyzeMeans$IDH_Mut1, levels=c(1,0))
toAnalyzeMeans$IDH_Mut2 <- factor(toAnalyzeMeans$IDH_Mut2, levels=c(1,0))

# Reformat into matrix (including adding back in lower tri just to get it to plot; won't affect original toAnalyzeMeans which does not ) to look at heatmap
toAnalyzeMeansFlipFlop <- toAnalyzeMeans[which(!toAnalyzeMeans$Patient1 == toAnalyzeMeans$Patient2),]
colnames(toAnalyzeMeansFlipFlop)[1] <- 'Patient2'
colnames(toAnalyzeMeansFlipFlop)[2] <- 'Patient1'
toAnalyzeMeansFlipFlop <- toAnalyzeMeansFlipFlop[,c(2,1,seq(3,11,1))]
toAnalyzeMeansToCast <- rbind(toAnalyzeMeans, toAnalyzeMeansFlipFlop)
meanDissimilarityMatrix <- acast(toAnalyzeMeansToCast[,c('Patient1','Patient2','dissimilarityMean')], Patient1 ~ Patient2, value.var="dissimilarityMean")
meanDissimilarityMatrix <- meanDissimilarityMatrix[,match(rownames(meanDissimilarityMatrix),colnames(meanDissimilarityMatrix))]
patientOrder <- as.character(colnames(meanDissimilarityMatrix))
annotationData <- subtypedata[,c('IDH_Mut','X1p19q','X7gain10loss','Tumor')]
rownames(annotationData) <- subtypedata$Patient
annotationData <- annotationData[patientOrder,]
annotationData$Tumor <- factor(annotationData$Tumor)
annotationData$IDH_Mut <- factor(annotationData$IDH_Mut)
annotationData$X1p19q <- factor(annotationData$X1p19q)
annotationData$X7gain10loss <- factor(annotationData$X7gain10loss)
pheatmap(meanDissimilarityMatrix,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "BrBG")))(100),
         show_rownames = TRUE, show_colnames = TRUE,
         breaks = seq(min(meanDissimilarityMatrix), max(meanDissimilarityMatrix), length = 100),
         annotation_col = annotationData,
         annotation_colors = list(IDH_Mut = c('1'='black','0'='white'),
                                  X1p19q = c('1'='blue','0'='white'),
                                  X7gain10loss = c('1'='red','0'='white'),
                                  Tumor = c('Recurrence1'='purple','Primary'='white'))
)


## Now look at means within subtype, and further by if they are inter or intra
toAnalyzeMeansWithinSubtype <- toAnalyzeMeans[which(toAnalyzeMeans$subtypeRelationship == 'within'),]
patientsToDrop <- c('P452') #we drop our unique non-cannonical GBM from this analysis
toAnalyzeMeansWithinSubtype <- toAnalyzeMeansWithinSubtype[which(!(toAnalyzeMeansWithinSubtype$Patient1 %in% patientsToDrop | toAnalyzeMeansWithinSubtype$Patient2 %in% patientsToDrop)),]
ggplot(toAnalyzeMeansWithinSubtype, aes(x=IDH_Mut1, y=dissimilarityMean, fill=patientRelationship)) +
  geom_boxplot(position=position_dodge(0.8)) +
  scale_fill_manual(values=c("gray28","gray72"))+
  labs(list(y = "Mean dissimilarity") )+
  theme(axis.text.x = element_text(size=20, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 20), axis.text.y = element_text(size=20, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))
for (relationship in unique(toAnalyzeMeansWithinSubtype$patientRelationship)){
  print(relationship)
  IDH_Mut <- toAnalyzeMeansWithinSubtype[which(toAnalyzeMeansWithinSubtype$IDH_Mut1==1 & toAnalyzeMeansWithinSubtype$patientRelationship == relationship),]$dissimilarityMean
  IDH_WT <- toAnalyzeMeansWithinSubtype[which(toAnalyzeMeansWithinSubtype$IDH_Mut1==0 & toAnalyzeMeansWithinSubtype$patientRelationship == relationship),]$dissimilarityMean
  print('IDH_mut vs IDH_wt')
  print(wilcox.test(IDH_Mut,IDH_WT, alternative='less'))
}
# create a data structure that also contains distance to do some correlations with
toAnalyzeWithinDistance <- toAnalyze[which(toAnalyze$relationship == 'intra' & !is.na(toAnalyze$L.Coordinate1)),]
patientArray <- data.frame(patient=character(), dissimilarity=numeric(), distance=numeric(), purityDifference=numeric(), stringsAsFactors=F)
for (p in unique(toAnalyzeWithinDistance$Patient1)){
  print(p)
  p <- as.character(p)
  toAnalyzeWithinDistanceLocalP <- toAnalyzeWithinDistance[which(toAnalyzeWithinDistance$Patient1==p),]
  tumorType <- unique(as.character(toAnalyzeWithinDistanceLocalP$Tumor1))
  for (c in 1:nrow(toAnalyzeWithinDistanceLocalP)){
    Sample1LPS <- toAnalyzeWithinDistanceLocalP[c,c("L.Coordinate1","P.Coordinate1","S.Coordinate1")]
    Sample2LPS <- toAnalyzeWithinDistanceLocalP[c,c("L.Coordinate2","P.Coordinate2","S.Coordinate2")]
    samples <- paste0(toAnalyzeWithinDistanceLocalP[c,'sample1'],',',toAnalyzeWithinDistanceLocalP[c,'sample2'])
    distance <- pairwiseDist2pointsLPS(Sample1LPS,Sample2LPS)
    purityDifference <- toAnalyzeWithinDistanceLocalP[c,]$purityDifference
    if (toAnalyzeWithinDistanceLocalP[c,'purity1'] >= purityCutoff & toAnalyzeWithinDistanceLocalP[c,'purity2'] >= purityCutoff){
      purityStatus <- 'highCCF'
    } else {
      purityStatus <- 'lowCCF'
    }
    names(distance) <- "distance"
    dissimilarity <- toAnalyzeWithinDistanceLocalP[c,]$dissimilarity
    patientArray <- rbind(patientArray, data.frame(patient=p, samples=samples, tumorType=tumorType, dissimilarity=dissimilarity, distance=distance, purityDifference=purityDifference, purityStatus=purityStatus, stringsAsFactors = F))
  }
}
patientArray <- patientArray[which(!is.na(patientArray$distance)),]
## normalize by distance and plot by patient ## Ran twice - once with P452 dropped for boxplots, since this is a hypermutated patient and we don't want it dumped in with the others, and once with it present for scatter plot where it is separate so okay to be there
patientArray$normedDissimilarity <- patientArray$dissimilarity/patientArray$distance
patientsToDrop <- c('P452') # 
patientArray <- patientArray[which(!(patientArray$patient %in% patientsToDrop)),]
unique(patientArray$patient) %in% patientOrderRec #check
patientOrderRNAseq <- patientOrderRec[patientOrderRec %in% unique(patientArray$patient)]
patientArray$patient <- factor(patientArray$patient, levels=patientOrderRNAseq)
colors <- as.character(colorKey[patientOrderRNAseq])
boxplot(patientArray$normedDissimilarity~patientArray$patient, col=colors, las=2, ylab='Distance-normalized expression dissimilarity (((1-R)/2)/mm)')
patientArrayHighCCF <- patientArray[which(patientArray$purityStatus == 'highCCF'),]
patientOrderRNAseq <- patientOrderRec[patientOrderRec %in% unique(patientArrayHighCCF$patient)]
patientArrayHighCCF$patient <- factor(patientArrayHighCCF$patient, levels=patientOrderRNAseq)
colors <- as.character(colorKey[patientOrderRNAseq])
boxplot(patientArrayHighCCF$normedDissimilarity~patientArrayHighCCF$patient, col=colors, las=2, ylab='Distance-normalized expression dissimilarity (((1-R)/2)/mm)')
patientMeans <- aggregate(patientArrayHighCCF$normedDissimilarity, by=list(patient=patientArrayHighCCF$patient), mean)
colnames(patientMeans) <- c('patient','normed_dissimilarity_means')
patientMeans$tumorType <- 'Primary'
patientMeans[which(patientMeans$patient %in% unique(patientArrayHighCCF[which(patientArrayHighCCF$tumorType=='Recurrence1'),]$patient)),]$tumorType <- 'Recurrent'
wilcox.test(patientMeans$normed_dissimilarity_means~patientMeans$tumorType, alternative='less')
boxplot(patientMeans$normed_dissimilarity_means~patientMeans$tumorType, col='grey', ylab='Mean normalized expression dissimilarity')
## plot for each patient (can change dissimilarity "Spatial distance (mm)" to purityDifference)
ggplot(patientArrayHighCCF, aes(x=distance, y=dissimilarity, color=patient)) +
  geom_point(size=.5) +
  scale_colour_manual(values=colors) +
  labs(list(x = "Spatial distance (mm)", y = "Dissimilarity") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank()) +
  geom_smooth(aes(colour=factor(patient)), method = "lm", se=F, size=.5)

## Do stat test for whether distance is significantly correlated with dissimilarity
patients <- unique(patientArray$patient)
colors <- rainbow(length(patients))
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
y=40
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- patientArray[which(patientArray$patient == patientID),]
  testResult <- cor.test(patientSubset$distance, patientSubset$dissimilarity, method="pearson")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  R=round(testResult$estimate,3)
  lmResult <- lm(patientSubset$dissimilarity~patientSubset$distance)
  m <- round(coef(lmResult)["patientSubset$distance"],2)
  b <- round(coef(lmResult)["(Intercept)"],2)
  label <- paste0('p=',p,', R=',R,', y = ',m,'x + ',b)
  dataText <- rbind(dataText, c(p,R,label,x,y,m,b,color,patientID), stringsAsFactors=F)
  y <- y + 1.5
}
colnames(dataText) <- c('p','R','label','x','y','m','b','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'spatial_distance_vs_dissimilarity_stats.txt'),sep='\t', quote=F, row.names = F)




## Plot mean purity vs mean dissimilarity

# Get list of unique patient pairings by selecting an arbitrary 'direction'
allPatientCombinations = combn(unique(toAnalyzeMeans$Patient1),2)
selectedComboDirections <- c()
for (n in 1:ncol(allPatientCombinations)){
  a <- as.character(allPatientCombinations[1,n])
  b <- as.character(allPatientCombinations[2,n])
  selectedComboDirections <- append(selectedComboDirections, paste0(a,'_',b))
}
withinCombinations <- paste0(unique(toAnalyzeMeans$Patient1),'_',unique(toAnalyzeMeans$Patient1))
selectedComboDirections <- append(selectedComboDirections,withinCombinations)

# Select these unique pairs from the degenerate dataframe
toAnalyzeMeansUnique <- toAnalyzeMeans[paste0(toAnalyzeMeans$Patient1,'_',toAnalyzeMeans$Patient2)%in%selectedComboDirections,]

# Ensure our HM patient isn't included
toAnalyzeMeansUnique <- toAnalyzeMeansUnique[which(!toAnalyzeMeansUnique$Patient1 == 'P452' | toAnalyzeMeansUnique$Patient2 == 'P452'),]

# Make mean purity vs mean dissimilarity scatterplot
ggplot(toAnalyzeMeansUnique[(toAnalyzeMeansUnique$subtypeRelationship=='within'),], aes(x = (purityMean1+purityMean2)/2, y = dissimilarityMean, color = factor(IDH_Mut1) , shape=patientRelationship ))+
  geom_point()+
  geom_smooth(method=lm,se=F) +
  xlab("Mean CCF") +
  ylab("Mean dissimilarity")+
  scale_color_manual(values=c("#000000",'#fd0006'),labels=c("IDH-mut","IDH-wt,\n7g10I"))+
  scale_shape_manual(values=c(17, 19))+
  theme(legend.title = element_blank())+
  theme(legend.text.align =0.5)


# Do stat test for whether dissimilarity is significantly correlated with mean purity
dissimilarityPurityStats = data.frame(patientRelationship = character(),
                                      IDH.Mut =  character(),
                                      p = numeric(),
                                      R = numeric())
for (pr in c('inter', 'intra')){
  for (idh_mut in c(0,1)){
    toAnalyzeMeans_subset = toAnalyzeMeansUnique[(toAnalyzeMeansUnique$subtypeRelationship=='within') &
                                             (toAnalyzeMeansUnique$patientRelationship==pr) &
                                             (toAnalyzeMeansUnique$IDH_Mut1==idh_mut),]
    testResult <- cor.test((toAnalyzeMeans_subset$purityMean1+toAnalyzeMeans_subset$purityMean2)/2, toAnalyzeMeans_subset$dissimilarityMean, method="pearson")
    p=round(testResult$p.value,3)
    R=round(testResult$estimate[['cor']],3)
    dissimilarityPurityStats = rbind(dissimilarityPurityStats, data.frame(patientRelationship = pr,
                                                                          IDH.Mut =  idh_mut,
                                                                          R,
                                                                          p))
  }
}

write.table(dissimilarityPurityStats, file=paste0(outputPath,outfolder,tag,'purity_vs_dissimilarity_stats.txt'),sep='\t', quote=F, row.names = F)

