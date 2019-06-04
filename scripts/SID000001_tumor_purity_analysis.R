# Created: 2018.05.28
# By: Stephanie R Hilz
# Usage: Assess the purity of samples; for those spatially mapped look at spatial distributions

library(ggplot2)

pairwiseDist2pointsLPS <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + (((p2[3]-p1[3]))^2) ))
}

# user-defined variables
outfolder <- 'SID000001_tumor_purity_analysis/'
tag <- 'SID000001'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)

# create color key by patient from subtype folder
colorKey <- subtypedata$Color
names(colorKey) <- subtypedata$Patient

# merge by patient ID
merged <- merge(data, subtypedata, by="Patient")

# remove samples that are lacking purity info at all (did not have exome done)
discard <- which(is.na(merged$WES_ID))
merged <- merged[-discard,]

# subset P260 by medial and lateral
merged[which(merged$Patient=='P260' & merged$SampleName %in% paste0('v',seq(10))),]$Patient <- 'P260-l'
merged[which(merged$Patient=='P260' & merged$SampleName %in% paste0('v',seq(from=11,to=20,by=1))),]$Patient <- 'P260-m'

# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF

# convert NAs to 0 (so that are counted; if NA means FACETs could not calculate because were so low)
merged[which(is.na(merged$purity)),]$purity <- .1

# set patient level order (Oligo, Astro, GBM, with Primaries always before Recurrences)
merged$Patient <- factor(merged$Patient, levels=patientOrderSplit)

# set colors (+ molsubtype) for plotting (here focused on subtype differences)
merged$molType <- 'IDH-mut_A' # Astro
merged$color <- subtypeColors['IDH-mut_A'] # green from bp Set1
merged[which(merged$IDH_Mut==0),]$molType <- 'IDH-wt' #red Set1
merged[which(merged$IDH_Mut==0),]$color <- subtypeColors['IDH-wt'] #red Set1
merged[which(merged$X1p19q==1),]$molType <- 'IDH-mut_O' #blue Set1
merged[which(merged$X1p19q==1),]$color <- subtypeColors['IDH-mut_O'] #blue Set1
color <- unique(merged[,c('Patient','color')])$color
names(color) <- unique(merged[,c('Patient','color')])$Patient
color <- color[patientOrderSplit] %>% unlist
molType <- unique(merged[,c('Patient','molType')])$molType
names(molType) <- unique(merged[,c('Patient','molType')])$Patient
molType <- molType[patientOrderSplit]

# plot boxplot
boxplot(merged$purity~merged$Patient, col = 'grey72', ylab = 'Estimated CCF', las='2')

# get out grade by patient 
grade <- unique(merged[,c('Patient','Grade')])$Grade
names(grade) <- unique(merged[,c('Patient','Grade')])$Patient
grade <- grade[patientOrderSplit]

# get mean,med, and var by patient
summary <- aggregate(merged[,'purity'], list(merged$Patient), mean)
colnames(summary) <- c("patientID","means")
summary$molType <- molType
summary$molType <- factor(molType, levels=c("IDH-mut_O","IDH-mut_A","IDH-wt"))
summary$grade <- grade
summary$grade <- as.factor(summary$grade)
summary$medians <- aggregate(merged[,'purity'], list(merged$Patient), median)$x
summary$variance <- aggregate(merged[,'purity'], list(merged$Patient), var)$x

# plot medians by molType (i.e. subtype)
boxplot(summary$medians~summary$molType, ylab = 'Estimated Median Purity', las='2', col='grey')
TukeyHSD(aov(summary$medians~summary$molType))

# look at variance for medians by subtype
aggregate(summary[,'medians'], list(summary$molType), var)

## look at difference in purity as it relates to spatial distance
pairwiseWithinPatient <- data.frame(patient=character(),
                                   samples=character(),
                                   spatial_distance=integer(),
                                   purity_difference=integer(),
                                   stringsAsFactors=FALSE) 
merged$Patient <- as.character(merged$Patient)
merged[which(merged$Patient=='P260-l' | merged$Patient=='P260-m'),]$Patient <- 'P260'
patientsSM <- as.character(unique(merged[which(!is.na(merged$L.Coordinate)),]$Patient))
for (patientID in patientsSM){
  print(patientID)
  samples <- as.character(merged[which(merged$Patient == patientID & !is.na(merged$L.Coordinate)),]$SampleName)
  print(samples)
  combinations <- combn(samples, 2) #will just have one column if using all samples
  for (c in 1:ncol(combinations)){
    samplesForCombo <- combinations[,c]
    print(samplesForCombo)
    sample_a <- samplesForCombo[1]
    sample_b <- samplesForCombo[2]
    samplesForComboOut <- paste0(samplesForCombo, collapse=',')
    coordinates <- merged[which(merged$SampleName %in% samplesForCombo & merged$Patient==patientID),c('L.Coordinate','P.Coordinate','S.Coordinate')]
    exomeDistanceMatrix <- as.matrix(dist(coordinates, method = "euclidean"))
    meanDistance <- mean(exomeDistanceMatrix[lower.tri(exomeDistanceMatrix)])
    purity_a <- merged[which(merged$SampleName==sample_a & merged$Patient==patientID),]$purity
    purity_b <- merged[which(merged$SampleName==sample_b & merged$Patient==patientID),]$purity
    purityDifference <- abs(purity_a-purity_b)
    pairwiseWithinPatient <- rbind(pairwiseWithinPatient, data.frame(patientID, samplesForComboOut, meanDistance,purityDifference))
  }
}
patientOrderToPlot <- patientOrderRec[patientOrderRec %in% unique(pairwiseWithinPatient$patientID)]
colors <- as.character(colorKey[patientOrderToPlot])
pairwiseWithinPatient$patientID <- factor(pairwiseWithinPatient$patientID, levels=patientOrderToPlot)
ggplot(pairwiseWithinPatient, aes(x=meanDistance, y=purityDifference, color=patientID)) +
  geom_point(size=.5) +
  scale_colour_manual(values=colors) +
  labs(list(x = "Spatial distance (mm)", y = "CCF Difference") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank()) +
  geom_smooth(aes(colour=factor(patientID)), method = "lm", se=F, size=.5)

## Do stat test for whether distance is significantly correlated with purity
patients <- unique(pairwiseWithinPatient$patientID)
colors <- rainbow(length(patients))
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
y=40
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- pairwiseWithinPatient[which(pairwiseWithinPatient$patientID == patientID),]
  testResult <- cor.test(patientSubset$meanDistance, patientSubset$purityDifference, method="spearman")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  rho=round(testResult$estimate,3)
  lmResult <- lm(patientSubset$purityDifference~patientSubset$meanDistance)
  m <- round(coef(lmResult)["patientSubset$meanDistance"],2)
  b <- round(coef(lmResult)["(Intercept)"],2)
  label <- paste0('p=',p,', rho=',rho,', y = ',m,'x + ',b)
  dataText <- rbind(dataText, c(p,rho,label,x,y,m,b,color,patientID), stringsAsFactors=F)
  y <- y + 1.5
}
colnames(dataText) <- c('p','R','label','x','y','m','b','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'spatial_distance_vs_CCF_difference_stats.txt'),sep='\t', quote=F, row.names = F)

## purity with distance from periph (used for purity section analysis)
mergedToPlot <- merged[which(!is.na(merged$DistCentroid)),]
patientOrderToPlot <- patientOrderRec[patientOrderRec %in% unique(mergedToPlot$Patient)]
mergedToPlot$Patient <- factor(mergedToPlot$Patient, levels=patientOrderToPlot)
colors <- as.character(colorKey[patientOrderToPlot])
ggplot(mergedToPlot, aes(x=DistPeriph, y=purity, color=Patient)) +
  geom_point(size=.5) +
  scale_colour_manual(values=colors) +
  labs(list(x = "Dist. Periph", y = "CCF") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank()) +
  geom_smooth(aes(colour=factor(Patient)), method = "lm",size=.5, se=F) 


## Do stat test for purity vs centroid or periph
patients <- unique(mergedToPlot$Patient)
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
y=40
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- mergedToPlot[which(mergedToPlot$Patient == patientID),]
  testResult <- cor.test(patientSubset$DistPeriph, patientSubset$purity, method="pearson")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  rho=round(testResult$estimate,3)
  lmResult <- lm(patientSubset$purity~patientSubset$DistPeriph)
  m <- round(coef(lmResult)["patientSubset$DistPeriph"],2)
  b <- round(coef(lmResult)["(Intercept)"],2)
  label <- paste0('p=',p,', rho=',rho,', y = ',m,'x + ',b)
  dataText <- rbind(dataText, c(p,rho,label,x,y,m,b,color,patientID), stringsAsFactors=F)
  y <- y + 1.5
}
colnames(dataText) <- c('p','rho','label','x','y','m','b','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'_distance_periph_vs_purity_difference_stats.txt'),sep='\t', quote=F, row.names = F)
