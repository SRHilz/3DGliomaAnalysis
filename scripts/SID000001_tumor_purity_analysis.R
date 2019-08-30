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
mergedSM <- merged[which(merged$SampleType == 'SM'),]
mergednonSM <- merged[which(merged$SampleType == 'non-SM'),]
stripchart(mergedSM$purity~mergedSM$Patient, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = c('#363795'), cex=.5)
stripchart(mergednonSM$purity~mergednonSM$Patient, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 17, col = c('#E05828'), cex=.5)

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

# enhancing vs non-enhancing
mergedWithE <- merged[which(!is.na(merged$MREnhancementWithContrast)),]
mergedWithE$MREnhancementWithContrast <- factor(mergedWithE$MREnhancementWithContrast, levels=unique(mergedWithE$MREnhancementWithContrast))
ggplot(data=mergedWithE, aes(x=Patient, y=purity, fill=MREnhancementWithContrast)) +
  geom_boxplot(position=position_dodge(1))

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

