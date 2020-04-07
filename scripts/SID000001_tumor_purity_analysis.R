# Created: 2018.05.28
# By: Stephanie R Hilz
# Usage: Assess the purity of samples; for those spatially mapped look at spatial distributions

library(ggplot2)
library(dplyr)
library(ggpmisc)
library(stats)
library(RColorBrewer)

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

# create 3D atlas cohort figure (before getting into purity)
merged_3DAtlas <- merged[which(merged$PatientStudy=='3DAtlas'),]
orderToUse <- patientOrder[which(patientOrder %in% merged_3DAtlas$Patient)]
samplesCount <- merged_3DAtlas[which(merged_3DAtlas$SampleType == 'SM'),] %>% count(Patient)
samplesCount$Patient <- factor(samplesCount$Patient, levels=rev(orderToUse))
ggplot(data=samplesCount, aes(x=Patient, y=n))+
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
  coord_flip()
range(samplesCount$n)
median(samplesCount$n)
pairwideDistances <- data.frame(patient=character(), samplePair=character(), distance=numeric(), stringsAsFactors=F)
for (p in unique(merged_3DAtlas$Patient)){
  print(p)
  localP <- merged_3DAtlas[which(merged_3DAtlas$Patient==p),]
  allSampleCombinations <- combn(localP$SampleName,2)
  for (c in 1:ncol(allSampleCombinations)){
    a <- allSampleCombinations[1,c]
    b <- allSampleCombinations[2,c]
    a.coords <- localP[which(localP$SampleName==a),c('L.Coordinate','P.Coordinate','S.Coordinate')]
    b.coords <- localP[which(localP$SampleName==b),c('L.Coordinate','P.Coordinate','S.Coordinate')]
    distance <- pairwiseDist2pointsLPS(a.coords,b.coords)
    toBind <- data.frame(c(patient=p, samplePair=paste0(a,'_','b'), distance=distance))
    pairwideDistances <- rbind(pairwideDistances, toBind)
  }
}
colnames(pairwideDistances) <- c('patient','samplePair','distance')
pairwideDistances$patient <- factor(pairwideDistances$patient, levels=orderToUse)
# violin plot of average pairwise distance
ggplot(pairwideDistances, aes(x=pairwideDistances$patient, y=distance)) + 
  geom_violin(width=1.2,fill='grey')+
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
  geom_boxplot(width=0.1)

# create tumor size plot
subtype3DAtlas <- subtypedata[which(subtypedata$PatientStudy=='3DAtlas'),]
toPlot <- subtype3DAtlas$T2Volume
names(toPlot) <- subtype3DAtlas$Patient
toPlot <- toPlot[patientOrder[patientOrder %in% names(toPlot)]]
barplot(toPlot, las=2 )

# subset P260 by medial and lateral
#merged[which(merged$Patient=='P260' & merged$SampleName %in% paste0('v',seq(10))),]$Patient <- 'P260-l'
#merged[which(merged$Patient=='P260' & merged$SampleName %in% paste0('v',seq(from=11,to=20,by=1))),]$Patient <- 'P260-m'

# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF

# convert NAs to 0 (so that are counted; if NA means FACETs could not calculate because were so low)
merged[which(is.na(merged$purity)),]$purity <- .1

# convert purities >1 to 1
merged[which(merged$purity >1),]$purity <- 1

# set colors (+ molsubtype) for plotting (here focused on subtype differences)
merged$molType <- 'IDH-mut_A' # Astro
merged$color <- subtypeColors['IDH-mut_A'] # green from bp Set1
merged[which(merged$IDH_Mut==0),]$molType <- 'IDH-wt' #red Set1
merged[which(merged$IDH_Mut==0),]$color <- subtypeColors['IDH-wt'] #red Set1
merged[which(merged$X1p19q==1),]$molType <- 'IDH-mut_O' #blue Set1
merged[which(merged$X1p19q==1),]$color <- subtypeColors['IDH-mut_O'] #blue Set1
color <- unique(merged[,c('Patient','color')])$color
names(color) <- unique(merged[,c('Patient','color')])$Patient
color <- color[patientOrder] %>% unlist
molType <- unique(merged[,c('Patient','molType')])$molType
names(molType) <- unique(merged[,c('Patient','molType')])$Patient

# get out grade by patient 
grade <- unique(merged[,c('Patient','Grade')])$Grade
names(grade) <- unique(merged[,c('Patient','Grade')])$Patient

# get mean,med, and var by patient
summary <- aggregate(merged[,'purity'], list(merged$Patient), mean)
colnames(summary) <- c("patientID","means")
summary$molType <- molType[summary$patientID]
summary$molType <- factor(molType, levels=c("IDH-mut_O","IDH-mut_A","IDH-wt"))
summary$grade <- grade[summary$patientID]
summary$grade <- as.factor(summary$grade)
summary$grade4 <- FALSE
summary$recurrence <- TRUE
summary[which(summary$patientID %in% merged[which(merged$Tumor=='Primary'),]$Patient),]$recurrence <- FALSE
summary[which(summary$grade == 4),]$grade4 <- TRUE
summary$IDHwtTERT <- FALSE
summary[which(summary$molType=='IDH-wt'),]$IDHwtTERT <- TRUE
summary$medians <- aggregate(merged[,'purity'], list(merged$Patient), median)$x
summary$mad <- aggregate(merged[,'purity'], list(merged$Patient), mad)$x
patientOrderCCFmad <- summary[order(summary$mad),]$patientID

# plot medians by mad colored by IDH-wtTPM
summary$shape <- 1
summary[which(summary$IDHwtTERT==TRUE),]$shape <- 2
summary$color <- 'black'
#summary[which(summary$IDHwtTERT==TRUE),]$color <- '#823119'
summary[which(summary$patientID=='P452'),]$shape <- 0
plot(summary$median, summary$mad, col=as.character(summary$color), pch=as.numeric(summary$shape), cex=2)

# plot variance by molType (i.e. subtype)
summaryNoP452 <- summary[which(!summary$patientID=='P452'),]
boxplot(summaryNoP452$mad~summaryNoP452$IDHwtTERT, ylab = 'Estimated Purity MAD', las='2', col='grey')
wilcox.test(summaryNoP452$mad~summaryNoP452$IDHwtTERT)

# plot boxplot ordered by variance in purity
merged$Patient <- factor(merged$Patient, levels=orderToUse)
boxplot(merged$purity~merged$Patient, col = 'grey72', ylab = 'Estimated CCF', las='2')
mergedSM <- merged[which(merged$SampleType == 'SM'),]
mergednonSM <- merged[which(merged$SampleType == 'non-SM'),]
stripchart(mergedSM$purity~mergedSM$Patient, vertical = TRUE, #separated out sm and non sm as gives control to color them separately
           method = "jitter", add = TRUE, pch = 20, col = c('#000000'), cex=.5)
stripchart(mergednonSM$purity~mergednonSM$Patient, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = c('#000000'), cex=.5)

# set patient level order (Oligo, Astro, GBM, with Primaries always before Recurrences)
merged$Patient <- factor(merged$Patient, levels=patientOrder)

# create colors for radial plots 
par(mar=c(1,1,1,1))
mergedSMNoP302 <- mergedSM[which(!mergedSM$Patient == 'P302'),]
distanceStepSize <- 0.01
for (p in unique(mergedSMNoP302$Patient)){
  print(p)
  patientSubset <- mergedSMNoP302[which(mergedSMNoP302$Patient==p),]
  metricOfInterest <- 'DistVR'
  values <- patientSubset[order(patientSubset[,metricOfInterest]),c(metricOfInterest,'purity')]
  if (metricOfInterest=='DistVR'){
    offset <- min(values[,metricOfInterest])
    values$normedDist <- round((values[,metricOfInterest]-offset)/(max(values[,metricOfInterest])-offset),2)
  } else {
    values$normedDist <- round(values[,metricOfInterest]/max(values[,metricOfInterest]),2)
  }
  # we first make sure each normed distance is uniqe. If it isn't, we take the mean of all purities at that distance
  redundantNormedDist <- names(table(values$normedDist)[which(table(values$normedDist) > 1)])
  for (rd in redundantNormedDist){
    meanPurity <- mean(values[which(values$normedDist == round(as.numeric(rd),2)),]$purity)
    values[which(values$normedDist == round(as.numeric(rd),2)),]$purity <- meanPurity
    values[which(values$normedDist == round(as.numeric(rd),2)),]$DistVR <- 9999 #tag to show is merged
    values <- unique(values)
  }
  # then we step through each distance break
  finalOutput <- data.frame(distance=numeric(), purity=numeric(), stringsAsFactors = F)
  for (d in 1:(nrow(values)-1)){
    print(d)
    currentDistance <- round(values$normedDist[d],2)
    currentPurity <- values$purity[d]
    nextDistance <- round(values$normedDist[d+1],2)
    nextPurity <- values$purity[d+1]
    if (d == 1){ #should only happen once
      if (currentDistance < distanceStepSize){currentDistance <- distanceStepSize} #for rare cases where normedDist is 0, we advance to step size
      distanceSpanStart <- seq(distanceStepSize,currentDistance,distanceStepSize)
      spanBlockStart <- cbind(distanceSpanStart, currentPurity)
      colnames(spanBlockStart) <- c('distance','purity')
      finalOutput <- rbind(finalOutput, spanBlockStart)
    } 
    if (round((currentDistance + distanceStepSize),2) ==  nextDistance){ # if 
      distanceSpanBetween <- nextDistance
    } else {
      distanceSpanBetween <- seq(currentDistance + distanceStepSize, nextDistance, distanceStepSize)
    }
    distanceSpanBetweenSize <- length(distanceSpanBetween)
    purityStepSize <- abs(nextPurity-currentPurity)/distanceSpanBetweenSize
    if (currentPurity > nextPurity){
      puritySpanBetween <- rev(seq(nextPurity, currentPurity, purityStepSize))[2:(distanceSpanBetweenSize+1)]
    } else if (currentPurity < nextPurity){
      puritySpanBetween <- seq(currentPurity, nextPurity, purityStepSize)[2:(distanceSpanBetweenSize+1)]
    }  else { #if they are equal
      puritySpanBetween <- rep(currentPurity, distanceSpanBetweenSize)
    }
    spanBlockBetween <- cbind(distanceSpanBetween,puritySpanBetween)
    colnames(spanBlockBetween) <- c('distance','purity')
    finalOutput <- rbind(finalOutput, spanBlockBetween)
  }
  if (metricOfInterest == 'DistPeriph'){
    finalOutput$distance <- 1.01-finalOutput$distance
    finalOutput <- finalOutput[order(finalOutput$distance),]
  } 
  # add on 0 and 1 to purity temporarily to include them in scale
  originalLength <- nrow(finalOutput)
  tmpAdd <- cbind(c(0,0),c(0,1))
  colnames(tmpAdd) <- c('distance','purity')
  finalOutput <- rbind(finalOutput, tmpAdd)
  # get colors including extras
  toColor <- finalOutput$purity
  colfunc <- colorRampPalette(c("black", "white"))
  mappedColors <- colfunc(length(toColor))[as.numeric(cut(toColor,breaks = length(toColor)))][1:length(toColor)]
  # remove tmp extras and add colors to data structure
  finalOutput <- finalOutput[1:originalLength,]
  finalOutput$purityColors <- mappedColors[1:originalLength]
  # plot!
  jpeg(paste0(outputPath, outfolder, tag, '_gradient_purity_by_', metricOfInterest, '_',p,'.jpeg'))
  plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
  for (i in finalOutput$distance){
    color <- as.character(finalOutput[which(finalOutput$distance==i),]$purityColors)
    abline(h=i, col=color, lwd=4)
  }
  dev.off()
}
#make gradient for colorbar
par(mfrow=c(1,1))
plot(rep(1,10), col=colfunc(10),pch=15, cex=3, xlim=c(1,40), axes=F, ann=F)


# enhancing vs non-enhancing
mergedWithE <- merged[which(!is.na(merged$MREnhancementWithContrast)),]
mergedWithE$MREnhancementWithContrast <- factor(mergedWithE$MREnhancementWithContrast, levels=unique(mergedWithE$MREnhancementWithContrast))
ggplot(data=mergedWithE, aes(y=purity, x=MREnhancementWithContrast, fill=MREnhancementWithContrast)) +
  geom_boxplot(position=position_dodge(1)) +
  facet_wrap(facets=mergedWithE$Patient)

# read in GBM 
gbm <- read.table(paste0(dataPath, '20190903_GBM_histology_radiology_metrics.txt'), header=T, sep='\t')
gbm$uniqueID <- paste0(gbm$Patient, gbm$SampleName)
merged$uniqueID <- paste0(merged$Patient, merged$SampleName)
gbmMerged <- merge(gbm, merged, by="uniqueID")
gbmMerged$Patient.y <- factor(gbmMerged$Patient.y)
gbmMerged$Target.BV.hyperplasia <- factor(gbmMerged$Target.BV.hyperplasia)
gbmMerged$Target.Necrosis <- factor(gbmMerged$Target.Necrosis)
ggplot(gbmMerged, aes(y=purity, x=Percent.necrosis)) +
  geom_boxplot(position=position_dodge(1)) +
  facet_wrap(facets=gbmMerged$Patient.y)
ggplot(gbmMerged, aes(x=Percent.necrosis, y=purity, colour=Patient.y)) +
  geom_point(size=.8) +
  geom_smooth(aes(colour=factor(Patient.y)), method = "lm", se=F)
plot(gbmMerged$purity,gbmMerged$Percent.necrosis)

# pull out mean clonal VAF info for each sample
vafs <- read.table(paste0(outputPath,'/SID000014_mutational_categorization_clonal_subclonal_by_patient/SID000014_vafs.txt'), sep='\t', header=T, stringsAsFactors = F)
snvCategories <- read.table(paste0(outputPath,'/SID000014_mutational_categorization_clonal_subclonal_by_patient/SID000014_0p3_min_purity_2_min_samples_snv_categories.txt'), sep='\t', header=T, stringsAsFactors=F)
toPlot <- data.frame(sampleID=character(),
                                meanTumorWideVAFs=numeric(),
                                patient=character(),
                                purity=numeric(),
                                PurityEstUsed=character(),
                                stringsAsFactors = F)
for (p in unique(snvCategories$patient)){
  print(p)
  tumorWideSNVs <- snvCategories[which(snvCategories$patient == p & snvCategories$category == 'tumor-wide'),]$snvs
  tumorWideVafsForPatient <- vafs[which(vafs$patient == p & vafs$SNVuniqueID %in% tumorWideSNVs),]
  meanTumorWideVafsForPatient <- aggregate(tumorWideVafsForPatient[,'vaf'], by=list(tumorWideVafsForPatient$sampleID), mean)
  colnames(meanTumorWideVafsForPatient) <- c('sample_type','meanTumorWideVAFs')
  meanTumorWideVafsForPatient$patient <- p
  print(meanTumorWideVafsForPatient)
  #merge with purity info
  purityForPatient <- merged[which(merged$Patient == p),c('sample_type','purity','PurityEstUsed')]
  meanTumorWideVafsForPatientMerged <- merge(meanTumorWideVafsForPatient, purityForPatient, by='sample_type')
  toPlot <- rbind(toPlot, meanTumorWideVafsForPatientMerged)
}
toPlot$x2_meanTumorWideVAFs <- 2*toPlot$meanTumorWideVAFs

## plot purity vs 2x meanVAF
toPlot$PurityEstUsed <- as.factor(toPlot$PurityEstUsed)
patients <- patientOrder[patientOrder %in% toPlot$patient]
toPlot$patient <- factor(toPlot$patient, levels=patients)
colors <- patients
toPlot$color <- 'black'
toPlot[which(toPlot$PurityEstUsed == 'FACETS'),]$color <- 'orange'
colors <- toPlot[which(toPlot)]
ggplot(toPlot, aes(x=x2_meanTumorWideVAFs, y=purity)) +
  geom_point(size=.8) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_smooth(method = "lm", se=F, size=.5)+
  labs(y='tumor purity', x='2 x mean tumor-wide VAFs') +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x.npc = "right", label.y.npc = 0.8, formula = y~x, 
               parse = TRUE, size = 3) + 
  stat_fit_glance(method = 'lm',
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.4, size = 3) +
  facet_wrap(~patient)



