# Created: 2018.05.28
# By: Stephanie R Hilz
# Usage: Assess the purity of samples; for those spatially mapped look at spatial distributions

library(ggplot2)
library(dplyr)
library(ggpmisc)
library(stats)
library(RColorBrewer)
library(kableExtra)
library(ggpubr)
library(ComplexHeatmap)
library(tidyr)

RMSE <- function(error) { sqrt(mean(error^2)) }

pairwiseDist2pointsLPS <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + (((p2[3]-p1[3]))^2) ))
}

countSamplesPatient <- function(merged){
  merged_3DAtlas <- merged[which(merged$PatientStudy=='3DAtlas'),]
  orderToUse <- patientOrder[which(patientOrder %in% merged_3DAtlas$Patient)]
  samplesCount <- merged_3DAtlas[which(merged_3DAtlas$SampleType == 'SM'),] %>% count(Patient)
  samplesCount$Patient <- factor(samplesCount$Patient, levels=rev(orderToUse))
  ggplot(data=samplesCount, aes(x=Patient, y=n))+
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
    coord_flip()
  print(paste0('Range of sample numbers: ',range(samplesCount$n)))
  print(paste0('Median of sample numbers: ',median(samplesCount$n)))
}

plotSampleDistances <- function(merged){
  merged_3DAtlas <- merged[which(merged$PatientStudy=='3DAtlas'),]
  orderToUse <- patientOrder[which(patientOrder %in% merged_3DAtlas$Patient)]
  pairwideDistances <- data.frame(patient=character(), samplePair=character(), distance=numeric(), stringsAsFactors=F)
  for (p in unique(merged_3DAtlas$Patient)){
    print(p)
    localP <- merged_3DAtlas[which(merged_3DAtlas$Patient==p),]
    sampn <- length(localP$SampleName)
    allSampleCombinations <- combn(localP$SampleName,2)
    for (c in 1:ncol(allSampleCombinations)){
      a <- allSampleCombinations[1,c]
      b <- allSampleCombinations[2,c]
      a.coords <- localP[which(localP$SampleName==a),c('L.Coordinate','P.Coordinate','S.Coordinate')]
      b.coords <- localP[which(localP$SampleName==b),c('L.Coordinate','P.Coordinate','S.Coordinate')]
      distance <- pairwiseDist2pointsLPS(a.coords,b.coords)
      toBind <- data.frame(c(patient=p, samplePair=paste0(a,'_',b), distance=distance, sampn=sampn))
      pairwideDistances <- rbind(pairwideDistances, toBind)
    }
  }
  colnames(pairwideDistances) <- c('patient','samplePair','distance','sampn')
  pairwideDistances$patient <- factor(pairwideDistances$patient, levels=orderToUse)
  sampn <- pairwideDistances[,c('patient','sampn')]
  sampn <- sampn[!duplicated(sampn),]
  # violin plot of average pairwise distance
  ggplot(pairwideDistances, aes(x=patient, y=distance)) + 
    geom_violin(width=1.2,fill='grey')+
    geom_text(data=sampn, aes(patient, Inf, label=paste0('',sampn)), inherit.aes=TRUE, parse=FALSE, size=3, vjust=2) +
    theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
    geom_boxplot(width=0.1)
}

addIndelToMat <- function(mat){
  # a function to add indels to the mat
  for (p in colnames(mat)){
    mutFile <- file.path(dataPath, 'mutations_avf',paste0(gsub('P','Patient',p),'.R.mutations.avf.txt'))
    mutData <- read.csv(mutFile, sep='\t', header=T)
    mutDataIndels <- mutData[mutData$algorithm=='Pindel' & mutData$decision=='retain',]$X.gene
    if (any(rownames(mat) %in% mutDataIndels)){
      mat[rownames(mat) %in% mutDataIndels,p] <- 'unassigned indel'
    }
  }
  return(mat)
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

# clear out purity estimate used as this is actually made in this script 
merged$PurityEstUsed <- NULL

# set colors (+ molsubtype) for plotting (here focused on subtype differences)
merged$molType <- 'IDH-mut_A' # Astro
merged$color <- subtypeColors['IDH-mut_A'] # green from bp Set1
merged[which(merged$IDH_Mut==0 & merged$TERT==1),]$molType <- 'IDH-wt-TERT' #red Set1
merged[which(merged$IDH_Mut==0 & merged$TERT==1),]$color <- subtypeColors['IDH-wt'] #red Set1
merged[which(merged$X1p19q==1),]$molType <- 'IDH-mut_O' #blue Set1
merged[which(merged$X1p19q==1),]$color <- subtypeColors['IDH-mut_O'] #blue Set1
merged[which(merged$IDH_Mut==0 & merged$TERT==0),]$molType <- 'IDH-wt-noTERT' #red Set1
merged[which(merged$IDH_Mut==0 & merged$TERT==0),]$color <- 'orange'
color <- unique(merged[,c('Patient','color')])$color
names(color) <- unique(merged[,c('Patient','color')])$Patient
color <- color[patientOrder] %>% unlist
molType <- unique(merged[,c('Patient','molType')])$molType
names(molType) <- unique(merged[,c('Patient','molType')])$Patient

# count number of samples per patient and get some stats
countSamplesPatient(merged)

# create violin plot of distances between samples in each patient
plotSampleDistances(merged)

# create tumor size plot
subtype3DAtlas <- subtypedata[which(subtypedata$PatientStudy=='3DAtlas'),]
toPlot <- subtype3DAtlas$T2Volume
names(toPlot) <- subtype3DAtlas$Patient
toPlot <- toPlot[patientOrder[patientOrder %in% names(toPlot)]]
barplot(toPlot, las=2 )

# convert NAs to 0 (so that are counted; if NA means FACETs could not calculate because were so low)
merged[which(is.na(merged$FACETS)),]$FACETS <- .1

# create an object for plotting with purity infomraiton
toPlot <- merged[,c('Patient','sample_type','FACETS','IDH1_VAF','PyClone')]
toPlot <- toPlot[!is.na(toPlot$PyClone),] # removes the few samples that did not have WES and thus PyClone estimates

# add to object IDH1-based purity
toPlot$IDH1_VAF_purity <- 2 * toPlot$IDH1_VAF
toPlot[which(toPlot$IDH1_VAF_purity > 1),]$IDH1_VAF_purity <- 1
toPlot$IDH1_VAF <- NULL

# add 2x mean tumor wide VAF
# pull out mean clonal VAF info for each sample 
toMerge <- c()
vafs <- read.table(paste0(outputPath,'/SID000014_mutational_categorization_clonal_subclonal_by_patient/SID000014_vafs.txt'), sep='\t', header=T, stringsAsFactors = F)
snvCategories <- read.table(paste0(outputPath,'/SID000014_mutational_categorization_clonal_subclonal_by_patient/SID000014_0p3_min_purity_2_min_samples_snv_categories.txt'), sep='\t', header=T, stringsAsFactors=F)
for (p in unique(snvCategories$patient)){
  print(p)
  tumorWideSNVs <- snvCategories[which(snvCategories$patient == p & snvCategories$category == 'tumor-wide'),]$snvs
  tumorWideVafsForPatient <- vafs[which(vafs$patient == p & vafs$SNVuniqueID %in% tumorWideSNVs),]
  meanTumorWideVafsForPatient <- aggregate(tumorWideVafsForPatient[,'vaf'], by=list(tumorWideVafsForPatient$sampleID), mean)
  colnames(meanTumorWideVafsForPatient) <- c('sample_type','meanTumorWideVAFs')
  meanTumorWideVafsForPatient$patient <- p
  print(meanTumorWideVafsForPatient)
  toMerge <- rbind(toMerge, meanTumorWideVafsForPatient)
}
toMerge$x2_meanTumorWideVAFs <- 2*toMerge$meanTumorWideVAFs
toMerge[which(toMerge$x2_meanTumorWideVAFs > 1),]$x2_meanTumorWideVAFs <- 1
toMerge$meanTumorWideVAFs <- NULL

# combine this with the other purity metrics
toMerge$uniqueID <- paste0(toMerge$patient, toMerge$sample_type)
toPlot$uniqueID <- paste0(toPlot$Patient, toPlot$sample_type)
toPlot <- merge(toPlot, toMerge, by='uniqueID')

# detrmine RMSE from 1:1 line for each purity metric
summary.stats <- data.frame(patient=character(), subtype=character(), r2=numeric(), rmse=numeric(), comparison=character(), stringsAsFactors = F)
for (p in unique(toPlot$Patient)){
  print(p)
  toPlotSubset <- toPlot[which(toPlot$Patient == p),]
  subtype <- molType[p]
  #toPlotSubset$PyClone_vs_FACETS_residuals <- residuals(lm(toPlotSubset$FACETS-toPlotSubset$PyClone ~ 0))
  r2 <- summary(lm(FACETS~PyClone, data=toPlotSubset))$r.squared #NA values are ignored
  rmse <-  RMSE(resid(lm(toPlotSubset$PyClone-toPlotSubset$FACETS ~ 0)))
  summary.stats <- rbind(summary.stats, data.frame(cbind(patient=p, subtype=subtype, r2, rmse, comparison='PyClone_vs_FACETS'), stringsAsFactors = F))
  r2 <- summary(lm(x2_meanTumorWideVAFs~PyClone, data=toPlotSubset))$r.squared #NA values are ignored
  rmse <-  RMSE(resid(lm(toPlotSubset$PyClone-toPlotSubset$x2_meanTumorWideVAFs ~ 0)))
  summary.stats <- rbind(summary.stats, data.frame(cbind(patient=p, subtype=subtype, r2=r2, rmse=rmse, comparison='PyClone_vs_x2Mean'), stringsAsFactors = F))
  r2 <- summary(lm(x2_meanTumorWideVAFs~FACETS, data=toPlotSubset))$r.squared #NA values are ignored
  rmse <-  RMSE(resid(lm(toPlotSubset$FACETS-toPlotSubset$x2_meanTumorWideVAFs ~ 0)))
  summary.stats <- rbind(summary.stats, data.frame(cbind(patient=p, subtype=subtype, r2=r2, rmse=rmse, comparison='FACETS_vs_x2Mean'), stringsAsFactors = F))
  r2 <- summary(lm(IDH1_VAF_purity~PyClone, data=toPlotSubset))$r.squared #NA values are ignored
  rmse <-  RMSE(resid(lm(toPlotSubset$PyClone-toPlotSubset$IDH1_VAF_purity ~ 0)))
  summary.stats <- rbind(summary.stats, data.frame(cbind(patient=p, subtype=subtype, r2=r2, rmse=rmse, comparison='PyClone_vs_IDH'), stringsAsFactors = F))
  r2 <- summary(lm(IDH1_VAF_purity~FACETS, data=toPlotSubset))$r.squared #NA values are ignored
  rmse <-  RMSE(resid(lm(toPlotSubset$FACETS-toPlotSubset$IDH1_VAF_purity ~ 0)))
  summary.stats <- rbind(summary.stats, data.frame(cbind(patient=p, subtype=subtype, r2=r2, rmse=rmse, comparison='FACETS_vs_IDH'), stringsAsFactors = F))
  
}
summary.stats$r2 <- as.numeric(summary.stats$r2)
summary.stats$rmse <- as.numeric(summary.stats$rmse)
summary.stats <- summary.stats[complete.cases(summary.stats),]
ggplot(summary.stats, aes(x=r2, color=comparison)) +
  geom_density() +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) 
ggplot(summary.stats, aes(x=rmse, color=comparison)) +
  geom_density() +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) 
ggplot(summary.stats, aes(x=subtype, y=rmse))+
  geom_boxplot()+
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
  facet_wrap(~comparison)
par(mfrow=c(1,1), mar=c(4,4,4,4))

# compute absolute difference between different metrics
toPlot$PurityEstUsed <- NA
toPlot$refPurityUsed <- NA
toPlot$divergentPyClone_vs_FACETS <- abs(toPlot$PyClone-toPlot$FACETS) > .1
agree <- which(!toPlot$divergentPyClone_vs_FACETS)
toPlot[agree,]$PurityEstUsed <- 'PyClone' # for those patients in which the purity agrees within .1 between PyClone and FACETS we use the PyClone
toPlot[agree,]$refPurityUsed <- 'FACETS' # also note the 2nd method 
toPlot$divergentPyClone_vs_x2Mean <- abs(toPlot$PyClone-toPlot$x2_meanTumorWideVAFs) > .1
agree <- which((!toPlot$divergentPyClone_vs_x2Mean) & is.na(toPlot$PurityEstUsed))
toPlot[agree,]$PurityEstUsed <- 'PyClone' # for those patients in which the purity agrees within .1 between PyClone and 2x mean VAF we use the PyClone
toPlot[agree,]$refPurityUsed <- 'x2Mean' 
toPlot$divergentFACETS_vs_x2Mean <- abs(toPlot$FACETS-toPlot$x2_meanTumorWideVAFs) > .1
agree <- which((!toPlot$divergentFACETS_vs_x2Mean) & is.na(toPlot$PurityEstUsed))
toPlot[agree,]$PurityEstUsed <- 'FACETS' # for those patients in which the purity agrees within .1 between PyClone and 2x mean VAF we use the FACETS
toPlot[agree,]$refPurityUsed <- 'x2Mean' 
toPlot$divergentPyClone_vs_IDH <- abs(toPlot$PyClone-toPlot$IDH1_VAF_purity) > .1
agree <- which((!toPlot$divergentPyClone_vs_IDH) & is.na(toPlot$PurityEstUsed))
toPlot[agree,]$PurityEstUsed <- 'PyClone' # for those patients in which the purity agrees within .1 between PyClone and 2x mean VAF we use the FACETS
toPlot[agree,]$refPurityUsed <- 'IDH1_VAF' 
toPlot$divergentFACETS_vs_IDH <- abs(toPlot$FACETS-toPlot$IDH1_VAF_purity) > .1
agree <- which((!toPlot$divergentFACETS_vs_IDH) & is.na(toPlot$PurityEstUsed))
toPlot[agree,]$PurityEstUsed <- 'FACETS' # for those patients in which the purity agrees within .1 between PyClone and 2x mean VAF we use the FACETS
toPlot[agree,]$refPurityUsed <- 'IDH1_VAF' 
# need to manually inspect above to then fill out below
manualPyClonex2 <- c('P300Recurrence1v2','P327Primary-v3','P327Primary-v6',
                     'P375Primary-v3','P507Primary-v5')
toPlot[toPlot$uniqueID %in% manualPyClonex2,]$PurityEstUsed <- 'PyClone'
toPlot[toPlot$uniqueID %in% manualPyClonex2,]$refPurityUsed <- 'x2Mean'
manualFACETSIDH <- c('P482Primary-v5')
toPlot[toPlot$uniqueID %in% manualFACETSIDH,]$PurityEstUsed <- 'FACETS'
toPlot[toPlot$uniqueID %in% manualFACETSIDH,]$refPurityUsed <- 'x2Mean'
# finally, create a label for these
toPlot$puritySets <- paste0(toPlot$PurityEstUsed, '_vs_',toPlot$refPurityUsed)
toPlot$puritySets <- factor(toPlot$puritySets, levels=unique(toPlot$puritySets))
# and set actual values for main purity and ref purity
toPlot$purity <- toPlot$PyClone
toPlot[toPlot$PurityEstUsed=='FACETS',]$purity <- toPlot[toPlot$PurityEstUsed=='FACETS',]$FACETS
toPlot$refPurity <- toPlot$FACETS
toPlot[toPlot$refPurityUsed=='IDH1_VAF',]$refPurity <- toPlot[toPlot$refPurityUsed=='IDH1_VAF',]$IDH1_VAF_purity
toPlot[toPlot$refPurityUsed=='x2Mean',]$refPurity <- toPlot[toPlot$refPurityUsed=='x2Mean',]$x2_meanTumorWideVAFs

# plot primary vs secondary purity
summary.stats.final <- data.frame(Patient=character(), r2=numeric(), rmse=numeric(), n=numeric(), stringsAsFactors = F)
for (p in unique(toPlot$Patient)){
  print(p)
  toPlotSubset <- toPlot[which(toPlot$Patient == p),]
  toPlotSubset <- toPlotSubset[(!is.na(toPlotSubset$refPurity)) &  (!is.na(toPlotSubset$purity)),]
  r2 <- summary(lm(refPurity~purity, data=toPlotSubset))$r.squared #NA values are ignored
  rmse <-  RMSE(resid(lm(toPlotSubset$purity-toPlotSubset$refPurity ~ 0)))
  n <- nrow(toPlotSubset)
  summary.stats.final <- rbind(summary.stats.final, data.frame(cbind(Patient=p, r2, rmse, n), stringsAsFactors = F))
}
summary.stats.final$r2 <- as.numeric(summary.stats.final$r2)
summary.stats.final$rmse <- as.numeric(summary.stats.final$rmse)
summary.stats.final <- summary.stats.final[complete.cases(summary.stats.final),]
ggplot(toPlot, aes(x=purity, y=refPurity, colour=puritySets)) +
  geom_point(size=.8) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
  #geom_smooth(aes(x=purity, y=refPurity), method = "lm", se=F, size=.5, fullrange=T, inherit.aes = F)+
  geom_text(aes(x = -Inf, y = Inf, label=paste0('n = ',n), group = Patient, vjust=3.5, hjust=-.4), data=summary.stats.final, parse=FALSE, size=3, inherit.aes=F)+
  geom_text(aes(x = -Inf, y = Inf, label=paste0('RMSE = ',round(rmse,2)), group = Patient, vjust=1.5, hjust=-.1), data=summary.stats.final, parse=FALSE, size=3, inherit.aes=F)+
  scale_color_brewer(palette="Set1")+
  labs(x='purity', y='secondary purity') +
  geom_abline(slope=1, intercept=0, col='grey') +
  ylim(0,1) +
  xlim(0,1) +
  #stat_poly_eq(aes(x=purity, y=refPurity, label = paste(..rr.label..)),
               #label.x.npc = "right", label.y.npc = 0.1, formula = y~x, 
               #parse = TRUE, size = 2, inherit.aes = F) + 
  #stat_fit_glance(method = 'lm',
  #                geom = 'text',
  #                aes(x=purity, y=refPurity, label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
  #                label.x.npc = 'right', label.y.npc = 0.4, size = 3, inherit.aes = F) +
  facet_wrap(~Patient)

# Add purity estimate to merged (also manually add final to samples spreadsheet so have hard record)
merged$uniqueID <- paste0(merged$Patient, merged$sample_type)
merged <- merge(merged, toPlot[,c('uniqueID','PurityEstUsed')], by='uniqueID')
merged$purity <- merged$PyClone
merged[which(merged$PurityEstUsed=='FACETS'),]$purity <- merged[which(merged$PurityEstUsed=='FACETS'),]$FACETS

# Visualize all IDH-wt tumors
mergedSM <- merged[which(merged$SampleType=='SM'),]
mergeSMIDHwtc <- mergedSM[which(mergedSM$molType=='IDH-wt-TERT'),]
for (p in unique(mergeSMIDHwtc$Patient)){
  print(p)
  print('total samples (SM only')
  total <- nrow(mergeSMIDHwtc[which(mergeSMIDHwtc$Patient == p),])
  print(total)
  print('samples with purity < .6 (SM only')
  belowThresh <- nrow(mergeSMIDHwtc[which(mergeSMIDHwtc$Patient == p & mergeSMIDHwtc$purity < .6),])
  print(belowThresh)
  if (belowThresh/total >= .5){
    print('majority')
  } else {
    print('minority')
  }
}

# Visualize all IDH-mut tumors
mergedSM <- merged[which(merged$SampleType=='SM'),]
mergeSMIDHmut <- mergedSM[which(mergedSM$molType %in% c('IDH-mut_A','IDH-mut_O')),]
for (p in unique(mergeSMIDHmut$Patient)){
  print(p)
  print('total samples (SM only')
  total <- nrow(mergeSMIDHmut[which(mergeSMIDHmut$Patient == p),])
  print(total)
  print('samples with purity < .6 (SM only')
  belowThresh <- nrow(mergeSMIDHmut[which(mergeSMIDHmut$Patient == p & mergeSMIDHmut$purity < .6),])
  print(belowThresh)
  if (belowThresh/total >= .5){
    print('majority')
  } else {
    print('minority')
  }
}

# get out grade by patient 
grade <- unique(merged[,c('Patient','Grade')])$Grade
names(grade) <- unique(merged[,c('Patient','Grade')])$Patient

# get mean,med, and var by patient (only for SM samples)
mergedSM <- merged[which(merged$SampleType=='SM'),]
summary <- aggregate(mergedSM[,'purity'], list(mergedSM$Patient), mean)
colnames(summary) <- c("patientID","means")
summary$molType <- molType[summary$patientID]
summary$molType <- factor(molType, levels=c("IDH-mut_O","IDH-mut_A","IDH-wt-TERT","IDH-wt-noTERT"))
summary$grade <- grade[summary$patientID]
summary$grade <- as.factor(summary$grade)
summary$grade4 <- FALSE
summary$recurrence <- TRUE
summary[which(summary$patientID %in% mergedSM[which(mergedSM$Tumor=='Primary'),]$Patient),]$recurrence <- FALSE
summary[which(summary$grade == 4),]$grade4 <- TRUE
summary$IDHwtTERT <- FALSE
summary[which(summary$molType=='IDH-wt-TERT'),]$IDHwtTERT <- TRUE
summary$medians <- aggregate(mergedSM[,'purity'], list(mergedSM$Patient), median)$x
summary$mad <- aggregate(mergedSM[,'purity'], list(mergedSM$Patient), mad)$x
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
summaryNoP452_Rec <- summaryNoP452[which(summaryNoP452$recurrence==T),]
summaryNoP452_Primary <- summaryNoP452[which(summaryNoP452$recurrence==F),]
wilcox.test(summaryNoP452$median~summaryNoP452$IDHwtTERT, alternative='greater')
wilcox.test(summaryNoP452$mad~summaryNoP452$IDHwtTERT, alternative='less')
ggplot(data = summaryNoP452, aes(y=medians, x=IDHwtTERT, fill = IDHwtTERT)) + 
  geom_boxplot(position="dodge", outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(y='medians', x='subtype') +
  scale_fill_manual(values=c('white','#11cc42')) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
  stat_compare_means(aes(group = IDHwtTERT), label = "p.format", method='wilcox.test', method.args = list(alternative = "less"))
ggplot(data = summaryNoP452, aes(y=mad, x=IDHwtTERT, fill = IDHwtTERT)) + 
  geom_boxplot(position="dodge", outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(y='mad', x='subtypes') +
  scale_fill_manual(values=c('white','#11cc42')) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
  stat_compare_means(aes(group = IDHwtTERT), label = "p.format", method='wilcox.test', method.args = list(alternative = "greater"))

# compare variance of medians (is the variance equal?)
hist(summaryNoP452[which(summaryNoP452$IDHwtTERT==F),]$medians)
fligner.test(medians ~ IDHwtTERT, data = summaryNoP452)

# plot boxplot of purities by patient
orderToUse <- patientOrder
mergedSM <- merged[which(merged$SampleType=='SM'),]
mergedSM$Patient <- factor(mergedSM$Patient, levels=orderToUse)
boxplot(mergedSM$purity~mergedSM$Patient, col = 'white', ylab = 'Estimated CCF', las='2', outline = FALSE, boxlty = 0, outlier.shape = NA, whisklty = 0, staplelty = 0)
stripchart(mergedSM$purity~mergedSM$Patient, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = c('#000000'), cex=.5) 

# set patient level order (Oligo, Astro, GBM, with Primaries always before Recurrences)
merged$Patient <- factor(merged$Patient, levels=patientOrder)

# bring in PyClone data with muts annotated for oncoprint
vafsPerSample <- read.table(mutationDataFile, sep='\t', header=T)

# make oncoprint
vafsPerSample$onco_type <- vafsPerSample$pyclone_category
vafsPerSample[which(vafsPerSample$clonal),]$onco_type <- 'clonal'
vafsPerSample[which(is.na(vafsPerSample$clonal)),]$onco_type <- 'unassigned snv'
toSpread <- vafsPerSample[vafsPerSample$SangerCancerGeneCensus.=='YES',c('SNVuniqueID','patient','onco_type')]
toSpread <- toSpread[!duplicated(toSpread),]
matByLocus <- spread(toSpread, 'patient','onco_type')
matByLocus$gene <- sapply(matByLocus$SNVuniqueID, function(x) strsplit(x,'_')[[1]][1])
mat <- c()
for (g in unique(matByLocus$gene)){
  geneSubset <- matByLocus[matByLocus$gene==g,colnames(matByLocus)[!colnames(matByLocus) %in% c('SNVuniqueID','gene')]]
  nastring <- paste(rep('NA',nrow(geneSubset)), collapse=',')
  collapsed <- sapply(colnames(geneSubset), function(x) paste(geneSubset[,x], collapse=','))
  collapsed <- gsub(nastring,NA, collapsed)
  collapsed <- gsub('NA,','',collapsed)
  collapsed <- gsub(',NA','',collapsed)
  mat <- rbind(mat, c(g, collapsed))
}
rownames(mat) <- mat[,1]
mat <- mat[,-1]
mat <- mat[rowSums(!is.na(mat)) > 2,] #require that mutation is in at least 2 patients
mat <- addIndelToMat(mat) # add indels for any genes that make the cut based on SNVs
# set up alter function for mat
col <- c(clonal="forestgreen",`tumor-wide`="red",shared="blue",private="orange", `unassigned snv`='black', `unassigned indel`='black')
alter_fun <- list(
  background = function(x,y,w,h) grid.rect(x,y,w*0.9,h*0.9,gp=gpar(fill='#e0e0e0',col=NA)),
  clonal = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.9, gp=gpar(fill=col["clonal"], col = NA)),
  `tumor-wide` = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.5, gp=gpar(fill=col["tumor-wide"], col = NA)),
  shared = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.4, gp=gpar(fill=col["shared"], col = NA)),
  private = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.2, gp=gpar(fill=col["private"], col = NA)),
  `unassigned snv` = function(x,y,w,h) grid.points(x,y,pch=16, gp=gpar(col=col["unassigned snv"])),
  `unassigned indel` = function(x,y,w,h) grid.points(x,y,pch=17, gp=gpar(col=col["unassigned indel"]))
)
# create bottom annotation
clustersPerPatient <- vafsPerSample[,c('patient','cluster_id','onco_type')]
clustersPerPatient <- clustersPerPatient[!duplicated(clustersPerPatient),]
clustersPerPatient <- clustersPerPatient[complete.cases(clustersPerPatient),]
numClustersPerPatient <- table(clustersPerPatient$patient)
clustersTypePerPatient <- table(clustersPerPatient$patient,clustersPerPatient$onco_type)
propShared <- clustersTypePerPatient[,'shared']/rowSums(clustersTypePerPatient)
propPrivate <- clustersTypePerPatient[,'private']/rowSums(clustersTypePerPatient)
# compute average number of clones per sample and patient
clonesPerSample <- vafsPerSample[,c('uniqueID2','cluster_id','pyclone_called')]
clonesPerSample <- clonesPerSample[complete.cases(clonesPerSample),]
clonesPerSample <- clonesPerSample[clonesPerSample$pyclone_called,]
clonesPerSample <- clonesPerSample[!duplicated(clonesPerSample),]
numClonesPerSample <- table(clonesPerSample$uniqueID2)
pycloneSamplesUsed <- apply(vafsPerSample[which(vafsPerSample$pyclone_category_used),c('patient','sampleID')],1,paste, collapse='')
clonesPerSampleByPatient <- numClonesPerSample[unique(pycloneSamplesUsed)]
clonesPerSampleByPatient <- data.frame(cbind(clonesPerSampleByPatient,sapply(names(clonesPerSampleByPatient), substr, start=1, stop=4)))
clonesPerSampleByPatient$clonesPerSampleByPatient <- as.numeric(clonesPerSampleByPatient$clonesPerSampleByPatient)
clonesPerSampleByPatient$V2 <- factor(clonesPerSampleByPatient$V2, levels=unique(clonesPerSampleByPatient$V2))
medianClonesPerSampleByPatientDf <- aggregate(clonesPerSampleByPatient~V2, data=clonesPerSampleByPatient, median)
medianClonesPerSampleByPatient <- medianClonesPerSampleByPatientDf$clonesPerSampleByPatient
names(medianClonesPerSampleByPatient) <- medianClonesPerSampleByPatientDf$V2
medianClonesPerSampleByPatient <- medianClonesPerSampleByPatient[colnames(mat)]
# tag recurrent samples
recurrent <- colnames(mat) %in% subsetRecurrent
# define colors
propClone.col.spectrum=colorRamp2(c(0,max(propShared)), c('white','black'))
medClone.col.spectrum=colorRamp2(c(0,5,max(medianClonesPerSampleByPatient)), c('lightyellow','darkorange','red'))
ann_colors <- list(
  propShared=propClone.col.spectrum,
  propPrivate=propClone.col.spectrum,
  medianClonesPerSampleByPatient=medClone.col.spectrum
)

df <- data.frame(cbind(propShared, propPrivate, medianClonesPerSampleByPatient, recurrent))
bottom_annotation <- HeatmapAnnotation(numClusts=anno_barplot(c(numClustersPerPatient)), df=df, col=ann_colors)

oncoColOrder <- match(patientOrder, colnames(mat))
oncoPrint(mat, 
          alter_fun=alter_fun, 
          col=col,
          column_order = oncoColOrder,
          bottom_annotation = bottom_annotation,
          show_column_names = TRUE)

# look at relationship between tp53 and medainClonesPerSample
tp53tumors <- ifelse(names(medianClonesPerSampleByPatient) %in% colnames(mat)[!is.na(mat['TP53',])], 'TP53 mut','TP53 wt')
names(tp53tumors) <- colnames(mat)
tp53tumors['P260'] <- 'TP53 mut' #this one is germline, so not called in MuTect but we know it is there.
subtypedata$TP53Mut <- tp53tumors[match(subtypedata$Patient, names(tp53tumors))]
subtypedata$MedClonesPerSample <- as.numeric(medianClonesPerSampleByPatient[match(subtypedata$Patient, names(medianClonesPerSampleByPatient))])
subtypedata$binary <- ifelse(subtypedata$IDH_Mut==1,'IDH mut','IDH wt TERT')
subtypedata[subtypedata$Patient=='P452',]$binary <- NA
subtypedata$numClones <- as.numeric(numClustersPerPatient[match(subtypedata$Patient, names(numClustersPerPatient))])
p <- ggboxplot(subtypedata, x = "TP53Mut", y = "MedClonesPerSample",
              color="TP53Mut", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "wilcox")

p <- ggboxplot(subtypedata, x = "binary", y = "MedClonesPerSample",
               color="binary", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "wilcox")

p <- ggboxplot(subtypedata, x = "TP53Mut", y = "numClones",
               color="TP53Mut", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "wilcox")

p <- ggboxplot(subtypedata, x = "binary", y = "numClones",
               color="binary", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "wilcox")

p <- ggboxplot(subtypedata, x = "TP53Mut", y = "numClones",
               color="TP53Mut", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "wilcox")
# plot purity by TP53 status; might also want to create a model that controls for purity
medianPurityPerPatient <- aggregate(purity~Patient, data=merged, median)
subtypedata$medianPurity <- medianPurityPerPatient[match(medianPurityPerPatient$Patient, subtypedata$Patient),]$purity
p <- ggboxplot(subtypedata, x = "TP53Mut", y = "medianPurity",
               color="TP53Mut", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "wilcox")


# look at relationship between metrics and purity 
# compute number of samples per patient
samplesPerPatient <- vafsPerSample[,c('patient','sampleID')]
samplesPerPatient <- samplesPerPatient[!duplicated(samplesPerPatient),]
numSamplesPerPatient <- table(samplesPerPatient$patient)
subtypedata$numSamples <- as.numeric(numSamplesPerPatient[match(subtypedata$Patient, names(numSamplesPerPatient))])
subtypedata$numClones <- as.numeric(numClustersPerPatient[match(subtypedata$Patient, names(numClustersPerPatient))])
# plot samples vs clones
ggplot(subtypedata, aes(x=numSamples, y=numClones)) + 
  geom_point() +
  geom_text(label=subtypedata$Patient, size=3, nudge_x=1) +
  geom_smooth(method="lm", color="blue", se=F) +
  scale_y_continuous(breaks=c(2,4,6,8,10,12,14)) +
  stat_cor(method = "pearson", label.x = 3, label.y = 12) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
ggplot(subtypedata, aes(x=T2Volume, y=numClones)) + 
  geom_point() +
  geom_text(label=subtypedata$Patient, size=3, nudge_x=1) +
  geom_smooth(method="lm", color="blue", se=F) +
  scale_y_continuous(breaks=c(2,4,6,8,10,12,14)) +
  stat_cor(method = "pearson", label.x = 3, label.y = 12)+
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
# add number of clones per sample to merged
vafsPerSample$uniqueID2 <- paste0(vafsPerSample$patient, vafsPerSample$sampleID)
merged$uniqueID2 <- paste0(merged$Patient, merged$sample_type)
merged$numClusters <- as.numeric(numClonesPerSample[match(merged$uniqueID2, names(numClonesPerSample))])
# compute mutations per sample - one for calls and one for vafs
mutsPerSample <- vafsPerSample[,c('uniqueID2','SNVuniqueID','mut_called','vaf')]
mutsPerSample <- mutsPerSample[complete.cases(mutsPerSample),]
mutsPerSample <- mutsPerSample[!duplicated(mutsPerSample),]
mutsPerSampleCalled <- mutsPerSample[mutsPerSample$mut_called,]
numMutsPerSampleCalled <- table(mutsPerSampleCalled$uniqueID2)
mutsPerSampleNonzero <- mutsPerSample[mutsPerSample$vaf > 0,]
numMutsPerSampleNonzero <- table(mutsPerSampleNonzero$uniqueID2)
merged$numMutsCalled <- as.numeric(numMutsPerSampleCalled[match(merged$uniqueID2, names(numMutsPerSampleCalled))])
merged$numMutsNonzero <- as.numeric(numMutsPerSampleNonzero[match(merged$uniqueID2, names(numMutsPerSampleNonzero))])
# plot purity vs clusters (I also looked at purity x coverage and coverage alone, and these were not sig associated)
ggplot(merged, aes(x=purity, y=numClusters)) + 
  geom_point() +
  geom_smooth(method="lm", color="blue", se=F) +
  scale_y_continuous(breaks=c(2,4,6,8,10,12,14)) +
  stat_cor(method = "pearson", label.x = .1, label.y = 12) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
# plot purity vs muts called
ggplot(merged[!merged$Patient=='P452',], aes(x=purity, y=numMutsCalled)) + #leaving off HM patient since is such an outlier
  geom_point() +
  geom_smooth(method="lm", color="blue", se=F) +
  stat_cor(method = "pearson", label.x = .1, label.y = 12) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
# plot purity vs muts detected at all
ggplot(merged[!merged$Patient=='P452',], aes(x=purity, y=numMutsNonzero)) + #leaving off HM patient since is such an outlier
  geom_point() +
  geom_smooth(method="lm", color="blue", se=F) +
  stat_cor(method = "pearson", label.x = .1, label.y = 12) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
# plot purity vs muts called with cutoff
ggplot(merged[(!merged$Patient=='P452') & merged$purity >= 0.2,], aes(x=purity, y=numMutsCalled)) + #leaving off HM patient since is such an outlier
  geom_point() +
  geom_smooth(method="lm", color="blue", se=F) +
  stat_cor(method = "pearson", label.x = .1, label.y = 12) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
# plot purity vs muts detected at all with cutoff
ggplot(merged[(!merged$Patient=='P452') & merged$purity >= 0.1,], aes(x=purity, y=numMutsNonzero)) + #leaving off HM patient since is such an outlier
  geom_point() +
  geom_smooth(method="lm", color="blue", se=F) +
  stat_cor(method = "pearson", label.x = .1, label.y = 12) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
# plot purity vs number of clones with cutoff
ggplot(merged[(!merged$Patient=='P452') & merged$purity >= 0.5,], aes(x=purity, y=numClusters)) + #leaving off HM patient since is such an outlier
  geom_point() +
  geom_smooth(method="lm", color="blue", se=F) +
  stat_cor(method = "pearson", label.x = .1, label.y = 12) +
  theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))


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



