# Created: 2018.08.29
# By: Stephanie R Hilz
# Usage: For a list of patients, calculates the number of clonal, shared, and 
#   private mutations in that patient for a specified number of samples that 
#   meet a specified minimum tumor purity.

library(ggplot2)
library(kableExtra)
library(reshape)
library(dplyr)

getPatient <- function(x){
  strsplit(x, '_')[[1]][1]
}

getMutsBin <- function(avfile){ #from Tali's phylo tree code, pulls out all covered muts
  muts <- read.delim(avfile, as.is=TRUE)
  if('X.gene' %in% colnames(muts)){
    geneIndex <- grep('X.gene',colnames(muts))
    colnames(muts)[geneIndex] <- 'gene'
  }
  
  sanger <- grepl("_sanger_status", paste(names(muts),collapse=""))
  covered <- "mutect_covered_in_all" %in% names(muts)
  
  ## filter out mutations according to rules:
  # remove indels unless sanger validated
  pindel.idx <- which(muts$algorithm == "Pindel")
  if(sanger) { pindel.idx <- intersect(pindel.idx, which(muts$any_sample_somatic != "SOMATIC")) }
  # remove SNVs not covered in all exomes unless sanger validated
  covered.idx <- c()
  if(covered) { covered.idx <- which(muts$algorithm == "MuTect" & is.na(muts$mutect_covered_in_all)) }
  if(covered & sanger) { covered.idx <- which(muts$algorithm == "MuTect" & is.na(muts$mutect_covered_in_all) & muts$any_sample_somatic != "SOMATIC") }
  # remove SNVs with any reads in normal unless sanger validated
  norm.idx <- norm.idx <- which(as.numeric(muts$Normal_ref_reads) > 0 )
  if(sanger) { norm.idx <- which(as.numeric(muts$Normal_ref_reads) > 0 & muts$any_sample_somatic != "SOMATIC") }
  # remove mutations defined above
  ridx <- c(pindel.idx, covered.idx)#ridx <- c(pindel.idx, covered.idx, norm.idx)
  muts <- muts[-ridx, ]
  
  ## pull out the relevant columns
  if(sanger) { calls <- which(sapply(names(muts), function(x) { regexpr("_called_sanger$", x) != -1} )) ## these columns contain the sanger-corrected mutations calls
  } else { calls <- which(sapply(names(muts), function(x) { regexpr("_called$", x) != -1} )); calls <- setdiff(calls, which(names(muts)=="samples_called")) }
  calls <- calls[order(names(muts)[calls])] ## need to be in order to compare to epigenetic trees
  muts$uniq <- paste(muts$gene, muts$contig, muts$position, muts$ref_allele, muts$alt_allele, sep="_")
  muts.bin <- muts[ , calls]
  rownames(muts.bin) <- muts$uniq
  
  ## (optional) remove variants that did not look good upon advanced variant filtering
  if (file.exists(avfFile)){
    avfData <- read.table(avfFile,header=TRUE, sep="\t", comment.char = '$')
    if('X.gene' %in% colnames(avfData)){
      geneIndex <- grep('X.gene',colnames(avfData))
      colnames(avfData)[geneIndex] <- 'gene'
    }
    avfData$uniq <- paste(avfData$gene, avfData$contig, avfData$position, avfData$ref_allele, avfData$alt_allele, sep="_")
    keep <- avfData[which(avfData$decision=='retain'),]$uniq
    muts.bin <- muts.bin[rownames(muts.bin) %in% keep,]
  }
  
  ## make logical
  for(i in 1:ncol(muts.bin)) { muts.bin[ , i] <- as.logical(muts.bin[ , i]) }
  
  ## only keep mutations for which at least one sample has the mutation (ie remove false negatives from sanger-validation)
  keep = vector()
  for(i in 1:nrow(muts.bin)) {
    keep[i] = FALSE
    for(j in 1:ncol(muts.bin)) {
      if(muts.bin[i,j] == TRUE) keep[i] = TRUE
    }
  }
  muts.bin <- muts.bin[which(keep), ]
  return(muts.bin)
}

#user-defined variables
outfolder <- 'SID000007_mutational_analysis_clonal_subclonal_by_patient/'
tag <- 'SID000007'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# specify cutoffs - which patients to look at, how many samples will be used, and how pure each must be;
#  if there are more samples than samplesPerPatientToUse, will average result from all possible combinations
#  of this number.
purityCutoff <- .7
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

# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF

# Filter based on purity and inspect numbers of samples per patient
merged <- merged[which(merged$purity >= purityCutoff),]
kable(table(merged$Patient)) %>%
        kable_styling(bootstrap_options = c("striped", "hover")) %>%
        row_spec(0, bold=TRUE, background = "#C2BFBA")

# Remove patients with less than the required number of high-CCF samples and reinspect
for (p in unique(merged$Patient)){
  if (nrow(merged[which(merged$Patient==p),]) < samplesPerPatientToUse){
    merged <- merged[which(!merged$Patient==p),]
  }
}
kable(table(merged$Patient)) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold=TRUE, background = "#C2BFBA")

# Set final patients to use
patientsToUse <- patientOrder[patientOrder %in% unique(merged$Patient)]

# Create an empty column for the periph and centroid coding
merged$peripheral <- NA
merged$centroid <- NA
merged[which(merged$DistPeriph <= 7),]$peripheral <- 1
merged[which(merged$DistPeriph > 7),]$peripheral <- 0
merged[which(merged$DistCentroid <= 17),]$centroid <- 1
merged[which(merged$DistCentroid > 17),]$centroid <- 0

# Create an empty data frame for storing the number of total mutations called per patient that are clonal, shared, or private - will only
#  be done for small subset of patients with comparable samples
categoryPerPatient <- data.frame(patient=character(),
                                           clonal=integer(),
                                           shared=integer(),
                                           private=integer(),
                                           distance=integer(),
                                           meanCCF=numeric(),
                                           stringsAsFactors=FALSE) 

for (p in patientsToUse){
  patientID <- as.character(gsub('P','Patient',p))
  print(patientID)
  avfFile <- paste0(dataPath,patientID,'.R.mutations.avf.txt')
  
  # run muts.bin function, which creates the same logical matrix used for phylo trees
  muts.bin <- getMutsBin(avfFile) #muts.bin for all samples in the patient, not yet subsetted
  
  # make the column names so they can string match our sample names
  colnames(muts.bin) <- gsub('[.]','-',colnames(muts.bin))
  colnames(muts.bin) <- gsub('_called','',colnames(muts.bin))
  
  # subset this by samples that fit our standards (those left in merged)
  highPuritySamples <- as.character(merged[which(merged$Patient == p),]$sample_type )
  muts.bin <- muts.bin[,highPuritySamples]  
  
  # set number of samples to use per patient if not previously specified
  if (is.na(samplesPerPatientToUse)){#if we don't specify a number, we use all
    samplesPerPatientToUse <- length(highPuritySamples)
  }
  
  # identify all possible combinations of the specified number of samples for this patient
  combinations <- combn(highPuritySamples, samplesPerPatientToUse) #will just have one column if using all samples
  
  # calculate number of clonal, shared, and private mutations for the patient
  clonalToAverage <- c()
  sharedToAverage <- c()
  privateToAverage <- c()
  for (c in 1:ncol(combinations)){
    samplesForCombo <- combinations[,c]
    #calcualte for particular sample set
    muts.bin.subset <- muts.bin[,samplesForCombo]
    clonalCombo <- length(rownames(muts.bin.subset[which(rowSums(muts.bin.subset)== ncol(muts.bin.subset)),]))
    sharedCombo <- length(rownames(muts.bin.subset[which(rowSums(muts.bin.subset)< ncol(muts.bin.subset) & rowSums(muts.bin.subset) > 1),]))
    privateCombo <- length(rownames(muts.bin.subset[which(rowSums(muts.bin.subset)== 1),]))
    #append to vector
    clonalToAverage <- append(clonalToAverage, clonalCombo)
    sharedToAverage <- append(sharedToAverage, sharedCombo)
    privateToAverage <- append(privateToAverage, privateCombo)
  }
  clonal <- mean(clonalToAverage)
  shared <- mean(sharedToAverage)
  private <- mean(privateToAverage)
  
  # calculate mean CCF
  meanCCF <- mean(merged[which(merged$Patient == p),]$purity)
  
  # finally also calculate average distance among samples (can only do if SM)
  meanDistanceToAverage <- c()
  for (c in 1:ncol(combinations)){
    samplesForCombo <- combinations[,c]
    coordinates <- merged[which(merged$sample_type %in% samplesForCombo & merged$Patient==p),c('L.Coordinate','P.Coordinate','S.Coordinate')]
    exomeDistanceMatrix <- as.matrix(dist(coordinates, method = "euclidean"))
    meanDistance <- mean(exomeDistanceMatrix[lower.tri(exomeDistanceMatrix)], na.rm=TRUE)
    meanDistanceToAverage <- append(meanDistanceToAverage, meanDistance)
    meanDistanceAverage <- mean(meanDistanceToAverage)
  }
    
  # add to patient-level information about type of mutation data frame
  categoryPerPatient_Local <- data.frame(patient=p,clonal=clonal,shared=shared,private=private, distance=meanDistanceAverage, meanCCF=meanCCF)
  categoryPerPatient <- rbind(categoryPerPatient, categoryPerPatient_Local)
}

## Create boxplot of average distance among samples
barplot(categoryPerPatient$distance, names.arg = categoryPerPatient$patient, las=2, ylab="Average pairwisde distance (mm)", col='black')

## Create boxplot of absolute tumor-wide in IDH-wt vs IDH-mut
IDHwt <- categoryPerPatient[which(categoryPerPatient$patient %in% c('P413','P454')),]#excluding 452 bc hypermutated
IDHmut <- categoryPerPatient[which(!categoryPerPatient$patient %in% c('P413','P454','P452')),]
boxplot(IDHmut$clonal, IDHwt$clonal, names=c('IDH-mut','IDH-wt'), col='grey', ylab="Absolute tumor-wide mutations per patient")
wilcox.test(IDHwt$clonal, IDHmut$clonal)

## Clonal, shared, and private per patient plot
categoryPerPatient$patient <- factor(categoryPerPatient$patient, levels = patientsToUse)
categoryPerPatient$clonalRatio <- (categoryPerPatient$clonal) /(rowSums(categoryPerPatient[,2:4]))
categoryPerPatient$sharedRatio <- (categoryPerPatient$shared) /(rowSums(categoryPerPatient[,2:4]))
categoryPerPatient$privateRatio <- (categoryPerPatient$private) /(rowSums(categoryPerPatient[,2:4]))
ratio3Table <- categoryPerPatient[,c('patient',paste0(c('clonal','shared','private'),'Ratio'))] 
ratio3toPlot <- melt(ratio3Table)
ratio3toPlot$variable <- gsub('Ratio','',ratio3toPlot$variable)
ratio3toPlot$variable <- factor(ratio3toPlot$variable, levels = c('private','shared','clonal'))
#plot
colors <- c('#9DC105','#EA9D15','#18BCE8')
ggplot(data = ratio3toPlot, aes(x = patient, y = value, fill = variable)) + 
  labs(x="Patients",y="Proportion of SNVs") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10, colour = 'black'), axis.title = element_text(size = 10), axis.text.y = element_text(size=10, colour = 'black'), panel.background = element_rect(fill = 'white', colour = 'black'))

## Plot for subclonal vs clonal
categoryPerPatient$subclonal <- categoryPerPatient$shared + categoryPerPatient$private
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

par(mfrow=c(1,3))
## Barchart subclonal
toPlot <- categoryPerPatient$subclonal
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'Number of subclonal mutations')

## Barchart clonal
toPlot <- categoryPerPatient$clonal
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'Number of tumor-wide mutations', col='#18BCE8',ylim=c(0,150))

## Barchart shared
toPlot <- categoryPerPatient$shared
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'Number of shared mutations', col='#EA9D15',ylim=c(0,150))

## Barchart private
toPlot <- categoryPerPatient$private
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'Number of private mutations', col='#9DC105', ylim=c(0,150))

## Barchart meanCCF
toPlot <- categoryPerPatient$meanCCF
names(toPlot) <- categoryPerPatient$patient
toPlot <- toPlot[patientsToUse]
par(las=2)
barplot(toPlot, ylab= 'MeanCCF', col='#069547')


