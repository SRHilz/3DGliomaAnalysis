# Created: 2018.08.29
# By: Stephanie R Hilz
# Usage: For a list of patients, calculates the number of clonal, shared, and 
#   private mutations in that patient for a specified number of samples that 
#   meet a specified minimum tumor purity; for patients with >6 high-purity samples,
#   an average of all possible combinations of 6 samples is taken.

library(ggplot2)
library(kableExtra)
library(reshape)
library(reshape2)
library(dplyr)

#user-defined variables
outfolder <- 'SID000007_mutational_analysis_clonal_subclonal_by_patient/'
tag <- 'SID000007'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# specify cutoffs - which patients to look at, how many samples will be used, and how pure each must be
samplesPerPatientToUse <- 6 #can also be NA; in this case all samples fitting the criteria will be used
requireSM <- FALSE

# read in sample data file
merged <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# subset by if has WES data
toUse <- merged[which(!is.na(merged$WES_ID)),]$WES_ID
merged <- merged[which(merged$WES_ID %in% toUse),]

# subset by sample type if desired
if (requireSM){
  merged <- merged[which(merged$SampleType=='SM'),]
}

# bring in mutation categories for each patient
vafsPerSample <- read.table(mutationDataFile, sep='\t', header=T)

# Remove patients with less than the required number of high-CCF samples and reinspect
patientsToUse <- unique(vafsPerSample[which(vafsPerSample$pyclone_category_nsamp >= samplesPerPatientToUse),]$patient)

# Also remove patients not in merged (should be none)
patientsToUse <- patientsToUse[patientsToUse %in% merged$Patient]

# Set final patients to use
patientsToUse <- patientOrder[patientOrder %in% patientsToUse]

# Subset our categories to this number
vafsPerSample_Filtered <- vafsPerSample[vafsPerSample$patient %in% patientsToUse,]

# Create an empty data frame for storing the number of total mutations called per patient that are clonal, shared, or private - will only
#  be done for small subset of patients with comparable samples
categoryPerPatient <- data.frame(patient=character(),
                                           clonal=integer(),
                                           tumor_wide=integer(),
                                           shared=integer(),
                                           private=integer(),
                                           total=integer(),
                                           un=integer(),
                                           stringsAsFactors=FALSE) 

for (p in patientsToUse){
  # subset variants to patient
  vafsPerSample_Filtered_subset <- vafsPerSample_Filtered[vafsPerSample_Filtered$patient==p,]
  
  # pull out all samples that were usable for pyclone
  samplesToUse <- unique(vafsPerSample_Filtered_subset[which(vafsPerSample_Filtered_subset$pyclone_category_used==T),]$sampleID)
  
  # pull out list of clonal variants
  clonal_muts <- unique(vafsPerSample_Filtered_subset[which(vafsPerSample_Filtered_subset$clonal),]$SNVuniqueID)

  # iterate through all possible combinations of samplesPerPatientToUse
  combinations <- combn(samplesToUse, samplesPerPatientToUse)
  toAverage <- c()
  for (c in 1:ncol(combinations)){
    localSamples <- combinations[,c]
    vafsPerSample_Filtered_subset_local <- vafsPerSample_Filtered_subset[vafsPerSample_Filtered_subset$sampleID %in% localSamples,]
    
    # create variant by sample call matrix
    mat <- acast(vafsPerSample_Filtered_subset_local[,c('SNVuniqueID','sampleID','pyclone_called')], SNVuniqueID~sampleID, value.var="pyclone_called")
    
    # drop NAs - these are unclassified by PyClone and we will compute stats on them later
    mat <- mat[complete.cases(mat),]
    
    # compute call numbers
    callNums <- rowSums(mat)
    
    # compute number of tumor-wide and clonal mutations
    tumor_wide_muts_local <- names(callNums[callNums == samplesPerPatientToUse]) # this might seem silly to do here - something tumor-wide in all samples will be tumor-wide in 6, but the concverse is not true, and that is why we do this here - to catch cases where a variant appears tumor-wide but is not
    tumor_wide_muts_local <- tumor_wide_muts[!tumor_wide_muts_local %in% clonal_muts]
    tumor_wide <- length(tumor_wide_muts_local)
    
    # compute number of clonal mutations
    clonal_muts_local <- tumor_wide_muts_local[tumor_wide_muts_local %in% clonal_muts]
    clonal <- length(clonal_muts)
    
    # compute number of shared mutations
    shared_muts_local <- names(callNums[callNums > 1 & callNums < samplesPerPatientToUse])
    shared <- length(shared_muts_local)
    
    # compute number of private mutations
    private_muts_local <- names(callNums[callNums == 1])
    private <- length(private_muts_local)
    
    # append to average
    toAverage <- rbind(toAverage, c(tumor_wide, clonal, shared, private))
  }
  toAverage <- data.frame(toAverage)
  colnames(toAverage) <- c('tumor_wide','clonal','shared','private')
  
  # compute final averaged numbers
  tumor_wide <- mean(toAverage$tumor_wide)
  clonal <- mean(toAverage$clonal)
  shared <- mean(toAverage$shared)
  private <- mean(toAverage$private)
  
  # total called 
  total_muts <- unique(vafsPerSample_Filtered_subset$SNVuniqueID)
  total <- length(total_muts)
  
  # uncategorized
  un_muts <- unique(vafsPerSample_Filtered_subset[is.na(vafsPerSample_Filtered_subset$pyclone_category),]$SNVuniqueID)
  un <- length(un_muts)
  
  # add to patient-level information about type of mutation data frame
  categoryPerPatient_Local <- data.frame(patient=p,clonal=clonal,tumor_wide=tumor_wide,shared=shared,private=private, total=total, un=un)
  categoryPerPatient <- rbind(categoryPerPatient, categoryPerPatient_Local)
}

## Add a total from the downsampled calls for all that pyclone could use
categoryPerPatient$total_ds <- rowSums(categoryPerPatient[,c('clonal','tumor_wide','shared','private')])

## Create boxplot of absolute clonal in IDH-wt vs IDH-mut
IDHwt <- categoryPerPatient[which(categoryPerPatient$patient %in% subsetIDHwtc),]#excluding 452 bc hypermutated
IDHmut <- categoryPerPatient[which(categoryPerPatient$patient %in% subsetIDHmut),]
boxplot(IDHmut$clonal, IDHwt$clonal, names=c('IDH-mut','IDH-wtc'), col='grey', ylab="Absolute clonal mutations per patient")
wilcox.test(IDHwt$clonal, IDHmut$clonal)

# [Consider adding here a model assessing this after controlling for purity + patient-level effects] - have decided not to focus on this for the paper since is not a main topic

## Clonal, shared, and private per patient plot
categoryPerPatient$patient <- factor(categoryPerPatient$patient, levels = patientsToUse)
categoryPerPatient$clonalRatio <- (categoryPerPatient$clonal) /categoryPerPatient$total_ds
categoryPerPatient$tumorwideRatio <- (categoryPerPatient$tumor_wide) /categoryPerPatient$total_ds
categoryPerPatient$sharedRatio <- (categoryPerPatient$shared) /categoryPerPatient$total_ds
categoryPerPatient$privateRatio <- (categoryPerPatient$private) /categoryPerPatient$total_ds
ratio3Table <- categoryPerPatient[,c('patient',paste0(c('clonal','tumorwide','shared','private'),'Ratio'))] 
ratio3toPlot <- melt(ratio3Table)
ratio3toPlot$variable <- gsub('Ratio','',ratio3toPlot$variable)
ratio3toPlot$variable <- factor(ratio3toPlot$variable, levels = c('clonal','tumorwide','shared','private'))
#plot
colors <- c('#000080','#18BCE8','#EA9D15','#9DC105','grey')
ggplot(data = ratio3toPlot, aes(x = patient, y = value, fill = variable)) + 
  labs(x="Patients",y="Proportion of SNVs") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10, colour = 'black'), axis.title = element_text(size = 10), axis.text.y = element_text(size=10, colour = 'black'), panel.background = element_rect(fill = 'white', colour = 'black'))

## Uncategorized vs categorized plot
categoryPerPatient$unRatio <- (categoryPerPatient$un) / categoryPerPatient$total
categoryPerPatient$catRatio <- 1 - categoryPerPatient$unRatio
ratio1Table <- categoryPerPatient[,c('patient',paste0(c('un','cat'),'Ratio'))] 
ratio1toPlot <- melt(ratio1Table)
ratio1toPlot$variable <- gsub('Ratio','',ratio1toPlot$variable)
ratio1toPlot$variable <- factor(ratio1toPlot$variable, levels = c('un','cat'))
#plot
colors <- c('grey','black')
ggplot(data = ratio1toPlot, aes(x = patient, y = value, fill = variable)) + 
  labs(x="Patients",y="Proportion of SNVs") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10, colour = 'black'), axis.title = element_text(size = 10), axis.text.y = element_text(size=10, colour = 'black'), panel.background = element_rect(fill = 'white', colour = 'black'))



## Plot for subclonal vs clonal
categoryPerPatient$subclonal <- categoryPerPatient$tumor_wide + categoryPerPatient$shared + categoryPerPatient$private
categoryPerPatient$subclonalRatio <- (categoryPerPatient$subclonal) /(categoryPerPatient$subclonal + categoryPerPatient$clonal)
ratio2Table <- categoryPerPatient[,c('patient',paste0(c('clonal','subclonal'),'Ratio'))] 
ratio2toPlot <- melt(ratio2Table)
ratio2toPlot$variable <- gsub('Ratio','',ratio2toPlot$variable)
ratio2toPlot$variable <- factor(ratio2toPlot$variable, levels = c('subclonal','clonal'))
# plot
colors <- c('#BF2026','#4379BD')
ggplot(data = ratio2toPlot, aes(x = patient, y = value, fill = variable)) + 
  labs(x="Patients",y="Mutation Clonality") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))


## Barchart subclonal
toPlot <- categoryPerPatient$subclonal
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'Number of subclonal mutations')

par(mfrow=c(1,3))
## Barchart clonal
toPlot <- categoryPerPatient$clonal
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'Number of clonal mutations', col='#000080',ylim=c(0,100))

## Barchart shared
toPlot <- categoryPerPatient$shared
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'Number of shared mutations', col='#EA9D15',ylim=c(0,200))

## Barchart private
toPlot <- categoryPerPatient$private
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'Number of private mutations', col='#9DC105', ylim=c(0,50))

## Barchart meanCCF
toPlot <- categoryPerPatient$meanCCF
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'MeanCCF', col='#069547')


