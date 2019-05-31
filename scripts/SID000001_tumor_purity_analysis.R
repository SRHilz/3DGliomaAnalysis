# Created: 2018.05.28
# By: Stephanie R Hilz
# Usage: Assess the purity of samples; for those spatially mapped look at spatial distributions

library(ggplot2)

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)

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
merged[which(is.na(merged$purity)),]$purity <- 0

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
boxplot(merged$purity~merged$Patient, col = color, ylab = 'Estimated CCF', las='2')

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
summary$grade <- as.factor(means$grade)
summary$medians <- aggregate(merged[,'purity'], list(merged$Patient), median)$x
summary$variance <- aggregate(merged[,'purity'], list(merged$Patient), var)$x

# plot medians by molType (i.e. subtype)
boxplot(summary$medians~summary$molType, ylab = 'Estimated Median Purity', las='2', col='grey')
TukeyHSD(aov(summary$medians~summary$molType))

# look at variance for medians by subtype
aggregate(summary[,'medians'], list(summary$molType), var)
