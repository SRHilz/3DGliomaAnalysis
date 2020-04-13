# Created: 2018.08.29
# By: Stephanie R Hilz
# Usage: For a list of patients, calculates the number of different mutations
#  between n pieces, as well as shared. For SM samples only.

library(ggplot2)
library(kableExtra)

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

# specify path where data files are located as well as file names
outfolder <- 'SID000010_mutational_analysis_shared_distinct_by_sample_subset/'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# specify cutoffs - which patients to look at, how many samples will be used, and how pure each must be;
#  if there are more samples than samplesPerPatientToUse, will average result from all possible combinations
#  of this number.
purityCutoff <- .7
patientsToUse <- c('Patient303','Patient327','Patient375','Patient453','Patient482','Patient413','Patient454','Patient475','Patient485','Patient260','Patient276','Patient300','Patient302','Patient450')
samplesPerPatientToUse <- 2 # cannot use NA for this analysis

# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)

# merge by patient ID
merged <- merge(data, subtypedata, by="Patient")

# set up color key
colorKey <- subtypedata$Color
names(colorKey) <- subtypedata$Patient

# subset by if has WES data
toUse <- data[which(!is.na(data$WES_ID)),]$WES_ID
merged <- merged[which(merged$WES_ID %in% toUse),]

# subset by sample type if desired
merged <- merged[which(merged$SampleType=='SM'),]

# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF

# Filter based on purity and inspect numbers of samples per patient
merged <- merged[which(merged$purity >= purityCutoff),]
kable(table(merged$Patient)) %>%
        kable_styling(bootstrap_options = c("striped", "hover")) %>%
        row_spec(0, bold=TRUE, background = "#C2BFBA")

callsPerSampleSubset <- data.frame(patient=character(),
                             samples=character(),
                                spatial_distance=integer(),
                                snvs_shared=integer(),
                                snvs_different=integer(),
                                tumorType=character(),
                                stringsAsFactors=FALSE) 

for (patientID in patientsToUse){
  print(patientID)
  patient_ID <- gsub('Patient','P',patientID)
  tumorType=unique(as.character(merged[which(merged$Patient==patient_ID),]$Tumor))
  
  # files to use
  avfFile <- paste0(dataPath,patientID,'.R.mutations.avf.txt')
  
  # run muts.bin function, which creates the same logical matrix used for phylo trees
  muts.bin <- getMutsBin(avfFile) #muts.bin for all samples in the patient, not yet subsetted
  
  # make the column names so they can string match our sample names
  colnames(muts.bin) <- gsub('[.]','-',colnames(muts.bin))
  colnames(muts.bin) <- gsub('_called','',colnames(muts.bin))
  
  # subset this by samples that fit our standards (those left in merged)
  highPuritySamples <- as.character(merged[which(merged$Patient == patient_ID),]$sample_type )
  muts.bin <- muts.bin[,highPuritySamples]  
  
  # identify all possible combinations of the specified number of samples for this patient
  combinations <- combn(highPuritySamples, samplesPerPatientToUse) #will just have one column if using all samples
  
  # calculate number of clonal, shared, and private mutations for the patient
  clonalToAverage <- c()
  sharedToAverage <- c()
  privateToAverage <- c()
  for (c in 1:ncol(combinations)){
    samplesForCombo <- combinations[,c]
    samplesForComboOut <- paste0(samplesForCombo, collapse=',')
    #calcualte for particular sample set
    muts.bin.subset <- muts.bin[,samplesForCombo]
    shared <- length(rownames(muts.bin.subset[which(rowSums(muts.bin.subset)== ncol(muts.bin.subset)),])) #called in all
    private <- length(rownames(muts.bin.subset[which(rowSums(muts.bin.subset)== 1),])) # called in only one
    #also calculate average distance among samples 
    coordinates <- merged[which(merged$sample_type %in% samplesForCombo & merged$Patient == patient_ID),c('L.Coordinate','P.Coordinate','S.Coordinate')]
    exomeDistanceMatrix <- as.matrix(dist(coordinates, method = "euclidean"))
    meanDistance <- mean(exomeDistanceMatrix[lower.tri(exomeDistanceMatrix)])
    #add to main data structure
    callsPerSampleSubset_Local <- data.frame(patient=patientID,samples=samplesForComboOut,spatial_distance=meanDistance,snvs_shared=shared, snvs_different=private, tumorType=tumorType)
    callsPerSampleSubset <- rbind(callsPerSampleSubset, callsPerSampleSubset_Local)
  }
}

## Normalize difference by distance and output to file
callsPerSampleSubset$normedDifference <- callsPerSampleSubset$snvs_different / callsPerSampleSubset$spatial_distance
colors <- as.character(colorKey[gsub('Patient','P',unique(callsPerSampleSubset$patient))])
boxplot(callsPerSampleSubset$normedDifference~callsPerSampleSubset$patient, col=colors, las=2, ylab='Distance-normalized genetic difference (#SNVs/mm)')
write.table(callsPerSampleSubset, file=paste0(outputPath,outfolder, 'SID000010_pairwise_metrics_genetic_distance.txt'), sep = '\t', row.names = F, quote=F)

## Do stat test for distinct mutations
dataText <- data.frame(p=numeric(), R=numeric(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
y=40
for (i in rev(seq_along(patientsToUse))){
  patientID=patientsToUse[i]
  print(patientID)
  color <- colors[i]
  patientSubset <- callsPerSampleSubset[which(callsPerSampleSubset$patient == patientID),]
  testResult <- cor.test(patientSubset$spatial_distance, patientSubset$snvs_different, method="pearson")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  R=round(testResult$estimate,3)
  lmResult <- lm(patientSubset$snvs_different~patientSubset$spatial_distance)
  m <- round(coef(lmResult)["patientSubset$spatial_distance"],2)
  b <- round(coef(lmResult)["(Intercept)"],2)
  dataText <- rbind(dataText, c(p,R,x,y,m,b,color,patientID), stringsAsFactors=F)
  y <- y + 1.5
}
colnames(dataText) <- c('p','R','x','y','m','b','color','patient')
dataText$adj.p <- NA
for (n in unique(dataText$patient)){
  dataTextNSubset.index <- which(dataText$patient == n & !is.na(dataText$p))
  dataText[dataTextNSubset.index,]$adj.p <- p.adjust(dataText[dataTextNSubset.index,]$p, method = "BH", n = length(dataText[dataTextNSubset.index,]$p))
}
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
dataText$label <- paste0('p=',dataText$adj.p,', R=',dataText$R)
write.table(dataText, file=paste0(outputPath,outfolder,'SID000010_pairwise_genetic_distance_spatial_distance_stats.txt'),sep='\t', quote=F, row.names = F)

## Do stats for Primary vs Recurrent mean snvs_different
patientMedians <- aggregate(callsPerSampleSubset$normedDifference, by=list(patient=callsPerSampleSubset$patient), median)
colnames(patientMedians) <- c('patient','normedDifference_medians')
patientVariance <- aggregate(callsPerSampleSubset$normedDifference, by=list(patient=callsPerSampleSubset$patient), var)
colnames(patientVariance) <- c('patient','normedDifference_variance')
patientMedians<- merge(patientVariance,patientMedians, by='patient')
patientMedians$tumorType <- 'Primary'
patientMedians[which(patientMedians$patient %in% unique(callsPerSampleSubset[which(callsPerSampleSubset$tumorType=='Recurrence1'),]$patient)),]$tumorType <- 'Recurrent'
patientMedians$tumorSubtype <- 'IDH-wt'
patientMedians[as.character(patientMedians$patient) %in% gsub('P','Patient',as.character(subtypedata[subtypedata$IDH_Mut == 1,]$Patient)),]$tumorSubtype <- 'IDH-mut'
wilcox.test(patientMedians$normedDifference_medians~patientMedians$tumorType)
wilcox.test(patientMedians$normedDifference_var~patientMedians$tumorType)
boxplot(patientMedians$normedDifference_medians~patientMedians$tumorType, col='grey', ylab='Median normalized genetic difference')
wilcox.test(patientMedians$normedDifference_medians~patientMedians$tumorSubtype)
boxplot(patientMedians$normedDifference_medians~patientMedians$tumorSubtype, col='grey', ylab='Median normalized genetic difference')

## Do multiple regression on distance vs snvs
fit <- lm(snvs_different ~ spatial_distance + patient, data=callsPerSampleSubset)
summary(fit)

## Plot different
colors <- as.character(colorKey[gsub('Patient','P',levels(callsPerSampleSubset$patient))])
ggplot(callsPerSampleSubset, aes(x=spatial_distance, y=snvs_different, color=patient)) +
  geom_point(size=.5) +
  scale_colour_manual(values=colors) +
  theme(axis.text.x = element_text(size=10, colour='black'), axis.title = element_text(size = 10), axis.text.y = element_text(size=10, colour='black'), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank(), legend.text=element_text(size=10), legend.title=element_text(size=10))+
  geom_smooth(aes(colour=factor(patient)), method = "lm", se=F, size=.5) #+
  #geom_text(
    #data    = dataText,
    #mapping = aes(x = x, y = y, label = label),
    #hjust   = -0.01,
    #colour=dataText$color,
    #size=5
  #)

## Do stat test for shared mutations
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
y=40
for (i in rev(seq_along(patientsToUse))){
  patientID=patientsToUse[i]
  print(patientID)
  color <- colors[i]
  patientSubset <- callsPerSampleSubset[which(callsPerSampleSubset$patient == patientID),]
  testResult <- cor.test(patientSubset$spatial_distance, patientSubset$snvs_shared, method="pearson")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  R=round(testResult$estimate,3)
  lmResult <- lm(patientSubset$snvs_shared~patientSubset$spatial_distance)
  m <- round(coef(lmResult)["patientSubset$spatial_distance"],2)
  b <- round(coef(lmResult)["(Intercept)"],2)
  label <- paste0('p=',p,', R=',R,', y = ',m,'x + ',b)
  dataText <- rbind(dataText, c(p,R,label,x,y,m,b,color,patientID), stringsAsFactors=F)
  y <- y + 2
}
colnames(dataText) <- c('p','R','label','x','y','m','b','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)

## Plot shared
colors <- as.character(colorKey[gsub('Patient','P',levels(callsPerSampleSubset$patient))])
ggplot(callsPerSampleSubset, aes(x=spatial_distance, y=snvs_shared, color=patient)) +
  geom_point(size=.5) +
  scale_colour_manual(values=colors) +
  theme(axis.text.x = element_text(size=10, colour='black'), axis.title = element_text(size = 10), axis.text.y = element_text(size=10, colour='black'), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank(), legend.text=element_text(size=10), legend.title=element_text(size=10))+
  geom_smooth(aes(colour=factor(patient)), method = "lm", se=F, size=.5) +
  geom_text(
    data    = dataText,
    mapping = aes(x = x, y = y, label = label),
    hjust   = -0.01,
    colour=colors,
    size=5
  )




