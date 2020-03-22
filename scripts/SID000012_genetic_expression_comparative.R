## Costello Lab
## 2018.11.01
## Stephanie Hilz

#Libraries
library(ggplot2)

#Functions
pairwiseDist2pointsLPS <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + (((p2[3]-p1[3]))^2) ))
}

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# User-defined variables
tag <- 'SID000012'
outfolder <- 'SID000012_genetic_expression_comparative/'
geneticFile <- paste0(outputPath,'SID000010_mutational_analysis_shared_distinct_by_sample_subset/SID000010_pairwise_metrics_genetic_distance.txt')
expressionFile <- paste0(outputPath,'SID000004_rnaseq_expression_distances_intra_inter/SID000004_pairwise_metrics_expression_dissimilarity.txt')

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)

# create color key by patient from subtype folder
colorKey <- subtypedata$Color
names(colorKey) <- subtypedata$Patient

## Get genetic data
# read in genetic distance data file and 
geneticDifference <- read.table(geneticFile, sep='\t', header = T, stringsAsFactors = F)
geneticDifference$patient <- gsub('Patient','P', geneticDifference$patient)

# read in expression dissimilarity data file and format to match genetic (final set will only include patients with CCF > threshold + RNAseq + WES data)
expressionDifference <- read.table(expressionFile, sep='\t', header = T, stringsAsFactors = F)
expressionDifference <- expressionDifference[which(expressionDifference$relationship == 'intra'),] # remove inter-patient comparisons
expressionDifference <- expressionDifference[which(expressionDifference$Patient1 %in% unique(geneticDifference$patient)),]
nameType1Patients <- unique(geneticDifference[which(grepl('-',geneticDifference$samples)),]$patient)
nameType2Patients <- unique(geneticDifference[which(!grepl('-',geneticDifference$samples)),]$patient)
nameType1Indices <- which(expressionDifference$Patient1 %in% nameType1Patients)
nameType2Indices <- which(expressionDifference$Patient1 %in% nameType2Patients)
expressionDifference[nameType1Indices,]$sample1 <- paste0(expressionDifference[nameType1Indices,]$Tumor1,'-',expressionDifference[nameType1Indices,]$SampleName1)
expressionDifference[nameType1Indices,]$sample2 <- paste0(expressionDifference[nameType1Indices,]$Tumor2,'-',expressionDifference[nameType1Indices,]$SampleName2)
expressionDifference[nameType2Indices,]$sample1 <- paste0(expressionDifference[nameType2Indices,]$Tumor1,expressionDifference[nameType2Indices,]$SampleName1)
expressionDifference[nameType2Indices,]$sample2 <- paste0(expressionDifference[nameType2Indices,]$Tumor2,expressionDifference[nameType2Indices,]$SampleName2)
expressionDifference$samplesOrder1 <- paste0(expressionDifference$sample1,',',expressionDifference$sample2)
expressionDifference$samplesOrder2 <- paste0(expressionDifference$sample2,',',expressionDifference$sample1)
expressionDifference$matchingOrder <- NA
expressionDifference$spatial_distance <- NA
expressionDifference$snvs_different <- NA
expressionDifference$normedDifference <- NA
for (p in unique(expressionDifference$Patient1)){
  print(p)
  sampleCombos <- geneticDifference[which(geneticDifference$patient == p),]$samples
  for (s in sampleCombos){
    print(s)
    matchingComboIndex <- which((expressionDifference$Patient1 == p  & expressionDifference$samplesOrder1 == s) | (expressionDifference$Patient1 == p  & expressionDifference$samplesOrder2 == s) )
    print(matchingComboIndex)
    if(length(matchingComboIndex) > 0){
      expressionDifference[matchingComboIndex,]$matchingOrder <- s
      expressionDifference[matchingComboIndex,]$snvs_different <- geneticDifference[which(geneticDifference$patient == p & geneticDifference$samples == s),]$snvs_different
      expressionDifference[matchingComboIndex,]$spatial_distance <- geneticDifference[which(geneticDifference$patient == p & geneticDifference$samples == s),]$spatial_distance
      expressionDifference[matchingComboIndex,]$normedDifference <- geneticDifference[which(geneticDifference$patient == p & geneticDifference$samples == s),]$normedDifference
    }
  }
}
expressionDifference <- expressionDifference[which(!is.na(expressionDifference$snvs_different)),c('Patient1','sample1','sample2','purity1','purity2','dissimilarity','matchingOrder','spatial_distance','snvs_different','normedDifference')]
expressionDifference$normedDissimilarity <- expressionDifference$dissimilarity/expressionDifference$spatial_distance

## Plot both normed and un-normed
patients <- unique(patientOrderRec[patientOrderRec %in% expressionDifference$Patient1])
colors <- as.character(colorKey[patients])
expressionDifference$Patient1 <- factor(expressionDifference$Patient1, levels=patients)
ggplot(expressionDifference, aes(x=snvs_different, y=dissimilarity, color=Patient1)) +
  geom_point(size=.5) +
  scale_colour_manual(values=colors) +
  labs(list(x = "Genetic difference (SNVs)", y = "Expression dissimilarity") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank()) +
  geom_smooth(aes(colour=factor(Patient1)), method = "lm", se=F, size=.5)
ggplot(expressionDifference, aes(x=normedDifference, y=normedDissimilarity, color=Patient1)) +
  geom_point(size=.5) +
  scale_colour_manual(values=colors) +
  labs(list(x = "Spatial Distance-normed genetic difference (SNVs/mm)", y = "Spatial distance-normed expression dissimilarity (dism./mm") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank()) +
  geom_smooth(aes(colour=factor(Patient1)), method = "lm", se=F, size=.5)

#do stat test for relationship between genetic distance and expression dissimilarity
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
y=40
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- expressionDifference[which(expressionDifference$Patient1 == patientID),]
  testResult <- cor.test(patientSubset$snvs_different, patientSubset$dissimilarity, method="pearson")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  R=round(testResult$estimate,3)
  lmResult <- lm(patientSubset$dissimilarity~patientSubset$snvs_different)
  m <- round(coef(lmResult)["patientSubset$distance"],2)
  b <- round(coef(lmResult)["(Intercept)"],2)
  label <- paste0('p=',p,', R=',R,', y = ',m,'x + ',b)
  dataText <- rbind(dataText, c(p,rho,label,x,y,m,b,color,patientID), stringsAsFactors=F)
  y <- y + 1.5
}
colnames(dataText) <- c('p','R','label','x','y','m','b','color','patient')
dataText$p.adj<- p.adjust(dataText$p, method='BH')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'genetic_difference_vs_expression_dissimilarity_stats.txt'),sep='\t', quote=F, row.names = F)

