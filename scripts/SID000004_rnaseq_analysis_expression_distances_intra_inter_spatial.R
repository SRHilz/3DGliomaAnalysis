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
ensmblRoot <- function(ensmblID){
  return(strsplit(as.character(ensmblID),'[.]')[[1]][1])
}

# User-defined variables
tag <- 'SID000004_20190529_first_submission'
outfolder <- 'SID000004_rnaseq_expression_distances_intra_inter/'
CPMFile <- paste0(dataPath,'SID000003_20190529_first_submission.symbol.coding.CPMs.csv')

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

## Get sample metadata
# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T, stringsAsFactors = F)

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
dissimilarityPairs <- melt(dissimilarity)
colnames(dissimilarityPairs) <- c('sample1','sample2','dissimilarity')

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

# Calculate means for each intra- or inter-patient comparison
toAnalyzeMeans <- data.frame(Patient1=character(), 
                             Patient2=character(), 
                             IDH_Mut1=character(),
                             IDH_Mut2=character(),
                             dissimilarityMean=numeric(),
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
  toAnalyzeMeans <- rbind(toAnalyzeMeans,data.frame(Patient1=Patient1, 
                                     Patient2=Patient2,
                                     IDH_Mut1=IDH_Mut1,
                                     IDH_Mut2=IDH_Mut2,
                                     dissimilarityMean=dissimilarityMean,
                                     purityDifference=purityDifference,
                                     purityFold=purityFold
                                     ))
}
toAnalyzeMeans$patientRelationship <- "inter"
toAnalyzeMeans[which(toAnalyzeMeans$Patient1 == toAnalyzeMeans$Patient2),]$patientRelationship <- 'intra'
toAnalyzeMeans$subtypeRelationship <- "between"
toAnalyzeMeans[which(toAnalyzeMeans$IDH_Mut1 == toAnalyzeMeans$IDH_Mut2),]$subtypeRelationship <- 'within'

# Reformat into matrix to look at heatmap
meanDissimilarityMatrix <- acast(toAnalyzeMeans[,c('Patient1','Patient2','dissimilarityMean')], Patient1 ~ Patient2, value.var="dissimilarityMean")
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
ggplot(toAnalyzeMeansWithinSubtype, aes(x=IDH_Mut1, y=dissimilarityMean, fill=patientRelationship)) +
  geom_boxplot(position=position_dodge(0.8)) +
  scale_fill_manual(values=c("gray28","gray72"))+
  labs(list(y = "Mean dissimilarity") )+
  theme(axis.text.x = element_text(size=20, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 20), axis.text.y = element_text(size=20, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))
for (relationship in unique(toAnalyzeMeansWithinSubtype$patientRelationship)){
  print(relationship)
  GBM <- toAnalyzeMeansWithinSubtype[which(toAnalyzeMeansWithinSubtype$Histology1=="GBM"),]$dissimilarityMean
  Astro <- toAnalyzeMeansWithinSubtype[which(toAnalyzeMeansWithinSubtype$Histology1=="Astro"),]$dissimilarityMean
  Oligo <- toAnalyzeMeansWithinSubtype[which(toAnalyzeMeansWithinSubtype$Histology1=="Oligo"),]$dissimilarityMean
  print('GBM vs Astro')
  print(wilcox.test(GBM,Astro))
  print('GBM vs Oligo')
  print(wilcox.test(GBM,Oligo))
  print('Astro vs Oligo')
  print(wilcox.test(Astro,Oligo))
}

## Also look at purity by subtype (can toggle between fold and difference), and further by if they are inter or intra
ggplot(toAnalyzeMeansWithinSubtype, aes(x=Histology1, y=purityFold, fill=patientRelationship)) +
  geom_boxplot(position=position_dodge(0.8)) +
  scale_fill_manual(values=c("gray28","gray72"))+
  labs(list(y = "Fold difference in purity") )+
  theme(axis.text.x = element_text(size=20, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 20), axis.text.y = element_text(size=20, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))

# correlate patient purity with dissimilarity
## plot for each patient
toAnalyzeWithin <- toAnalyze[which(toAnalyze$relationship == 'intra'),]
ggplot(toAnalyzeWithin, aes(purityDistance, dissimilarity)) +
  geom_point(aes(colour=Patient1)) +
  labs(list(y = "Expression dissimilarity", x = "Difference in purity") )+
  theme(axis.text.x = element_text(size=20, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 20), axis.text.y = element_text(size=20, color="black"), panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_smooth(aes(colour=factor(Patient1)), method = "lm", se=F) 

## do stat test for purity with dissimilarity
patients <- unique(toAnalyzeWithin$patient)
colors <- rainbow(length(patients))
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
y=40
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- toAnalyzeWithin[which(toAnalyzeWithin$patient == patientID),]
  testResult <- cor.test(patientSubset$purityDifference, patientSubset$dissimilarity, method="pearson")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  R=round(testResult$estimate,3)
  lmResult <- lm(patientSubset$dissimilarity~patientSubset$purityDifference)
  m <- round(coef(lmResult)["patientSubset$distance"],2)
  b <- round(coef(lmResult)["(Intercept)"],2)
  label <- paste0('p=',p,', R=',R,', y = ',m,'x + ',b)
  dataText <- rbind(dataText, c(p,R,label,x,y,m,b,color,patientID), stringsAsFactors=F)
  y <- y + 1.5
}
colnames(dataText) <- c('p','R','label','x','y','m','b','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)

# create a data structure that also contains distance to do some correlations with
toAnalyzeWithinDistance <- toAnalyze[which(toAnalyze$relationship == 'intra' & !is.na(toAnalyze$L.Coordinate1)),]
patientArray <- data.frame(patient=character(), dissimilarity=numeric(), distance=numeric(), purityDifference=numeric(), stringsAsFactors=F)
for (p in unique(toAnalyzeWithinDistance$Patient1)){
  print(p)
  p <- as.character(p)
  toAnalyzeWithinDistanceLocalP <- toAnalyzeWithinDistance[which(toAnalyzeWithinDistance$Patient1==p),]
  for (c in 1:nrow(toAnalyzeWithinDistanceLocalP)){
    Sample1LPS <- toAnalyzeWithinDistanceLocalP[c,c("L.Coordinate1","P.Coordinate1","S.Coordinate1")]
    Sample2LPS <- toAnalyzeWithinDistanceLocalP[c,c("L.Coordinate2","P.Coordinate2","S.Coordinate2")]
    distance <- pairwiseDist2pointsLPS(Sample1LPS,Sample2LPS)
    purityDifference <- toAnalyzeWithinDistanceLocalP[c,]$purityDifference
    names(distance) <- "distance"
    dissimilarity <- toAnalyzeWithinDistanceLocalP[c,]$dissimilarity 
    patientArray <- rbind(patientArray, data.frame(patient=p, dissimilarity=dissimilarity, distance=distance, purityDifference=purityDifference, stringsAsFactors = F))
  }
}
## plot for each patient (can change dissimilarity "Spatial distance (mm)" to purityDifference)
unique(patientArray$patient) %in% patientsToUse #check
patientArray$patient <- factor(patientArray$patient, levels=patientsToUse)
colors <- as.character(colorKey[patientsToUse])
ggplot(patientArray, aes(x=distance, y=dissimilarity, color=patient)) +
  geom_point() +
  scale_colour_manual(values=colors) +
  labs(list(x = "Spatial distance (mm)", y = "Dissimilarity") )+
  theme(axis.text.x = element_text(size=20, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 20), axis.text.y = element_text(size=20, color="black"), panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_smooth(aes(colour=factor(patient)), method = "lm", se=F) 



## Do stat test for distinct mutations
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

## purity with distance from periph move to purity section
mergedSM <- merged[which(merged$SampleType=="SM" & merged$Patient %in% patientsToUse & !is.na(merged$DistCentroid)),]
patients <- as.character(unique(mergedSM$Patient))
mergedSM$Patient <- factor(mergedSM$Patient, levels=patientsToUse)
colors <- as.character(colorKey[patientsToUse])
ggplot(mergedSM, aes(x=DistPeriph, y=purity, color=Patient)) +
  geom_point() +
  scale_colour_manual(values=colors) +
  labs(list(x = "Dist. Periph", y = "CCF") )+
  theme(axis.text.x = element_text(size=20, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 20), axis.text.y = element_text(size=20, color="black"), panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_smooth(aes(colour=factor(Patient)), method = "lm", se=F) 


## Do stat test for purity vs centroid or periph
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
y=40
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- mergedSM[which(mergedSM$Patient == patientID),]
  testResult <- cor.test(patientSubset$DistPeriph, patientSubset$purity, method="pearson")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  R=round(testResult$estimate,3)
  lmResult <- lm(patientSubset$purity~patientSubset$DistPeriph)
  m <- round(coef(lmResult)["patientSubset$distance"],2)
  b <- round(coef(lmResult)["(Intercept)"],2)
  label <- paste0('p=',p,', R=',R,', y = ',m,'x + ',b)
  dataText <- rbind(dataText, c(p,R,label,x,y,m,b,color,patientID), stringsAsFactors=F)
  y <- y + 1.5
}
colnames(dataText) <- c('p','R','label','x','y','m','b','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)