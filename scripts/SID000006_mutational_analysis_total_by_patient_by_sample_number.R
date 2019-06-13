# Created: 2018.08.29
# By: Stephanie R Hilz
# Usage: For a list of patients, calculates the average total number of mutations
#  identified when 1, 2, 3, 4...n number of samples are used. It furthermore
#  notes this calculation for 2 through n-1 samples for the case where the distance between
#  samples is minimized vs maximized, so that this can be compared to the average.

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

# specify cutoffs - which patients to look at, how many samples will be used, and how pure each must be;
#  if there are more samples than samplesPerPatientToUse, will average result from all possible combinations
#  of this number.
purityCutoff <- .7
patientsToUse <- c('Patient303','Patient327','Patient375','Patient453','Patient413','Patient454') # patients with >1 tumor of a particular type with >=7 high-CCF spatially-mapped samples
outfolder <- 'SID000006_mutational_analysis_total_by_patient_by_sample_number/'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

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

# require that all be SM
merged <- merged[which(merged$SampleType=='SM'),]

# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF

# Filter based on purity and inspect numbers of samples per patient
merged <- merged[which(merged$purity >= purityCutoff),]
kable(table(merged$Patient)) %>%
        kable_styling(bootstrap_options = c("striped", "hover")) %>%
        row_spec(0, bold=TRUE, background = "#C2BFBA")

# Create an empty data frame for storing the number of mutations called in each sample based on n and combination
mutTotalPerSampleSubset <- data.frame(patientID=character(),
                                           n=integer(),
                                           samplesForComboFused=character(),
                                           meanDistance=integer(),
                                           tag=character(),#the tag can be one of four levels - NA, distMin, distMax, or sampleMax
                                           mutsTotal=integer(),
                                           mutsClonal=integer(),
                                           Histology=character(),
                                           Tumor=character(),
                                           stringsAsFactors=FALSE) 

for (patientID in patientsToUse){
  
  print(patientID)
  patient_ID <- gsub('Patient','P',patientID)
  Histology <- unique(as.character(merged[which(merged$Patient == patient_ID),]$Histology))
  Tumor <- unique(as.character(merged[which(merged$Patient == patient_ID),]$Tumor))
  
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
  
  # k is total number of samples for this patient (need to be of correct purity + have exome)
  k <- ncol(muts.bin)
  
  # now for n from 1 to k
  for (n in 1:k){
    
    # define some empty variables we will update later - this is to help us denote which is the min and max distance combos by n
    distanceMin <- c(1000000, NA) #start with a value very large so all will be smaller, NA will be replaced with sample combo
    distanceMax <- c(0, NA) #start with a value very small so all will be larger, NA will be replaced with sample combo
    
    # identify all possible combinations of the specified number of samples for this patient
    combinations <- combn(highPuritySamples, n) #will just have one column if using all samples
    
    # calculate number of total mutations for the patient plus distance between them
    for (c in 1:ncol(combinations)){
      
      # specify sample subset
      samplesForCombo <- combinations[,c]
      samplesForComboFused <- paste(combinations[,c], collapse=',')
      
      # calcualte total number of muts called for this particular sample subset set as well as clonal only for that subset
      muts.bin.subset <- as.data.frame(muts.bin[,samplesForCombo])
      mutsTotal <- nrow(as.data.frame(muts.bin.subset[which(rowSums(muts.bin.subset) > 0),])) # pulls out all muts called in at least one sample of this subset
      mutsClonal <- nrow(as.data.frame(muts.bin.subset[which(rowSums(muts.bin.subset)== ncol(muts.bin.subset)),]))
      
      # calculate average distance between all samples for this particular sample subset set
      if (n > 1){
        coordinates <- merged[which(merged$sample_type %in% samplesForCombo & merged$Patient==patient_ID),c('L.Coordinate','P.Coordinate','S.Coordinate')]
        exomeDistanceMatrix <- as.matrix(dist(coordinates, method = "euclidean"))
        meanDistance <- mean(exomeDistanceMatrix[lower.tri(exomeDistanceMatrix)])
      } else {
        meanDistance <- NA
      }
      
      # determine if this distance is larger or smaller than others for this patient so far
      if (n > 1 & n < k){
        if (meanDistance > distanceMax[1]){
          distanceMax[1] <- meanDistance
          distanceMax[2] <- samplesForComboFused
        }
        if (meanDistance < distanceMin[1]){
          distanceMin[1] <- meanDistance
          distanceMin[2] <- samplesForComboFused
        }
      }
      
      # check if we have hit k to know whether this is the instance we should tag as sampleMax
      if (n==k){
        tag <- 'sampleMax'
      } else {
        tag <- NA
      }
      
      # rbind to final output structure
      names <- c('patientID','n','samplesForComboFused','meanDistance','tag','mutsTotal','mutsClonal','Histology','Tumor')
      mutTotalPerSampleSubset <- rbind.data.frame(mutTotalPerSampleSubset,c(patientID, n, samplesForComboFused, meanDistance, tag, mutsTotal, mutsClonal, Histology, Tumor), stringsAsFactors = F)
      colnames(mutTotalPerSampleSubset) <- names
    }
    
    # only if we are at n for which 1 < n < k is true
    if (n > 1 & n < k){
      print("Finishing up and determining tags")
      # for this n, tag the sample combo that ended up being the sample max
      mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$samplesForComboFused == distanceMax[2] & 
                                    mutTotalPerSampleSubset$patientID == patientID &
                                    mutTotalPerSampleSubset$n == n),]$tag <- 'distMax'
    
      # for this n, tag the sample combo that ended up being the sample min
      mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$samplesForComboFused == distanceMin[2] & 
                                    mutTotalPerSampleSubset$patientID == patientID &
                                    mutTotalPerSampleSubset$n == n),]$tag <- 'distMin'
    }
  }
}

# save data into a supplemental table
write.table(mutTotalPerSampleSubset, file=paste0(outputPath, outfolder, "SID000006_TableS_mutational_analysis_total_by_patient_by_sample_number.txt"), quote=F, row.names=F)

# whip data into the proper types so that we can both do some math on it and plot it
mutTotalPerSampleSubset$mutsTotal <- as.numeric(mutTotalPerSampleSubset$mutsTotal)
mutTotalPerSampleSubset$mutsClonal <- as.numeric(mutTotalPerSampleSubset$mutsClonal)
mutTotalPerSampleSubset$patientID <- factor(mutTotalPerSampleSubset$patientID, levels = patientsToUse)
mutTotalPerSampleSubset$n <- factor(mutTotalPerSampleSubset$n, levels = order(unique(as.integer(mutTotalPerSampleSubset$n))))

# create a new column which will be the mutsTotal-(mutsClonal for the sampleMax of that patient) = mutsSubclonal
mutTotalPerSampleSubset$mutsSubclonal <- as.numeric(NA)
for (patientID in patientsToUse){
  mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID),]$mutsSubclonal <- mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID),]$mutsTotal - mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$tag == 'sampleMax'),]$mutsClonal
}

## Subclonal by number of samples plot
ggplot(mutTotalPerSampleSubset, aes(x=n, y=mutsSubclonal, color=patientID)) +
  geom_point(aes(shape=patientID)) +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values=shapes) +
  theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank(), legend.text=element_text(size=15), legend.title=element_text(size=15))

## Subclonal by number of samples plot for each n relationship with distance
  # this first part makes the data labels with p value and R2
mutTotalPerSampleSubset$meanDistance <- as.numeric(mutTotalPerSampleSubset$meanDistance)
mutTotalPerSampleSubsetUsableN <- mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$n %in% c(2:7)),]
dataText <- data.frame(n=character(), p=numeric(), R=numeric(), x=numeric(), y=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
for (n in unique(mutTotalPerSampleSubsetUsableN$n)){
  print(n)
  y = 65
  for (i in rev(seq_along(patientsToUse))){
    patient=patientsToUse[i]
    print(patient)
    color <- colors[i]
    subset <- mutTotalPerSampleSubsetUsableN[which(mutTotalPerSampleSubsetUsableN$patientID==patient & mutTotalPerSampleSubsetUsableN$n == n),c('meanDistance','mutsSubclonal')]
    subset <- subset[complete.cases(subset),]
    meanDistance =  subset$meanDistance
    mutsSubclonal <- subset$mutsSubclonal
    if (length(meanDistance) > 1){
      testResult <- cor.test(meanDistance, mutsSubclonal, method="pearson")
      p=formatC(testResult$p.value,format = "e", digits = 2)
      R=round(testResult$estimate,3)
    } else {
      p=NA
      R=NA
    }
    dataText <- rbind(dataText, c(n,p,R,x,y,color,patient), stringsAsFactors=F)
    y <- y + 5
  }
}
colnames(dataText) <- c('n','p','R','x','y','color','patient')
dataText$adj.p <- NA
for (n in unique(dataText$n)){
  dataTextNSubset.index <- which(dataText$n == n & !is.na(dataText$p))
  dataText[dataTextNSubset.index,]$adj.p <- p.adjust(dataText[dataTextNSubset.index,]$p, method = "BH", n = length(dataText[dataTextNSubset.index,]$p))
}
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
dataText$label <- paste0('p=',dataText$adj.p,', R=',dataText$R)
write.table(dataText, file=paste0(outputPath,outfolder,'SID000006_mutsSubclonal_distance_n_stats.txt'),sep='\t', quote=F, row.names = F)
# plot
ggplot(mutTotalPerSampleSubsetUsableN, aes(x=meanDistance, y=mutsSubclonal, color=patientID)) +
  geom_point(aes(shape=patientID)) +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values=shapes) +
  theme(axis.text.x = element_text(size=10, colour='black'), axis.title = element_text(size = 10), axis.text.y = element_text(size=10, colour='black'), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank(), legend.text=element_text(size=10), legend.title=element_text(size=10))+
  facet_wrap(vars(n)) +
  geom_smooth(method=lm,se=FALSE, fullrange=TRUE) +
  geom_text(
    data    = dataText,
    mapping = aes(x = x, y = y, label = label),
    hjust   = -0.01,
    colour=dataText$color
  )


# now make plots showing the relationship between p value (from correlation above) and n
dataText$patient <- factor(dataText$patient, levels=patientsToUse)
dataText$log2p <- -log(as.numeric(dataText$p),2)
ggplot(dataText, aes(x=n, y=log2p, color=patient, group=patient)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values=shapes) +
  theme(axis.text.x = element_text(size=20, colour='black'), axis.title = element_text(size = 20), axis.text.y = element_text(size=20, colour='black'), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank(), legend.text=element_text(size=15), legend.title=element_text(size=15))

# now make plots showing the relationship between R (from correlation above) and n
dataText$R <- as.numeric(dataText$R)
ggplot(dataText, aes(x=n, y=R, color=patient, group=patient)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values=shapes) +
  theme(axis.text.x = element_text(size=20, colour='black'), axis.title = element_text(size = 20), axis.text.y = element_text(size=20, colour='black'), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank(), legend.text=element_text(size=15), legend.title=element_text(size=15))

# create a new data structure which will just have the average subclonal for each n per pateint
mutSubclonalMeanPerPatientPerN <- data.frame(patientID=character(),
                                             n=numeric(),
                                             mutsSubclonalMean=numeric(),
                                             mutsSubclonalMin=numeric(),
                                             mutsSubclonalMax=numeric(),
                                             mutsSubclonalSD=numeric(),
                                             stringsAsFactors=FALSE) 
mutTotalPerSampleSubset$n <- as.numeric(as.character(mutTotalPerSampleSubset$n))
for (patientID in patientsToUse){
  print(patientID)
  k <- max(mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID),]$n)
  for (n in 1:k){
    print(n)
    mutsSubclonalMean <- mean(mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$n == n),]$mutsSubclonal)
    mutsSubclonalSD <- sd(mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$n == n),]$mutsSubclonal)
    if (n > 1 & n < k){
      mutsSubclonalMin <- mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$n == n & mutTotalPerSampleSubset$tag == 'distMin'),]$mutsSubclonal
      mutsSubclonalMax <- mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$n == n & mutTotalPerSampleSubset$tag == 'distMax'),]$mutsSubclonal
    } else {
      mutsSubclonalMin <- mutsSubclonalMean
      mutsSubclonalMax <- mutsSubclonalMean
    }
    mutSubclonalMeanPerPatientPerN <- rbind.data.frame(mutSubclonalMeanPerPatientPerN, c(patientID, n, mutsSubclonalMean, mutsSubclonalMin, mutsSubclonalMax,mutsSubclonalSD), stringsAsFactors = F)
    }
}
colnames(mutSubclonalMeanPerPatientPerN) <- c('patientID','n','mutsSubclonalMean','mutsSubclonalMin','mutsSubclonalMax', 'mutsSubclonalSD')
mutTotalPerSampleSubset$n <- factor(as.character(mutTotalPerSampleSubset$n), levels = order(unique(as.integer(mutTotalPerSampleSubset$n))))

## Subclonal mean by number of samples plot
mutSubclonalMeanPerPatientPerN$n <- factor(as.character(mutSubclonalMeanPerPatientPerN$n), levels = order(unique(as.integer(mutSubclonalMeanPerPatientPerN$n))))
mutSubclonalMeanPerPatientPerN$patientID <- factor(mutSubclonalMeanPerPatientPerN$patientID, levels=patientsToUse)
mutSubclonalMeanPerPatientPerN$mutsSubclonalMean <- as.numeric(mutSubclonalMeanPerPatientPerN$mutsSubclonalMean)
mutSubclonalMeanPerPatientPerN$mutsSubclonalMin <- as.numeric(mutSubclonalMeanPerPatientPerN$mutsSubclonalMin)
mutSubclonalMeanPerPatientPerN$mutsSubclonalMax <- as.numeric(mutSubclonalMeanPerPatientPerN$mutsSubclonalMax)
mutSubclonalMeanPerPatientPerN$mutsSubclonalSD <- as.numeric(mutSubclonalMeanPerPatientPerN$mutsSubclonalSD)

ggplot(mutSubclonalMeanPerPatientPerN, aes(x=n, y=mutsSubclonalMean, group=patientID, colour=patientID)) +
  geom_point(aes(shape=patientID, colour=patientID)) +
  geom_line(aes(colour = patientID)) +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values=shapes) +
  geom_errorbar(aes(ymin=mutsSubclonalMean-mutsSubclonalSD, ymax=mutsSubclonalMean+mutsSubclonalSD), width=.2, position=position_dodge(0.05)) +
  theme(axis.text.x = element_text(size=10, color="black"),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank())

# create a new data structure which will just have the average clonal for each n per pateint
mutClonalMeanPerPatientPerN <- data.frame(patientID=character(),
                                             n=numeric(),
                                             mutsClonalMean=numeric(),
                                             mutsClonalMin=numeric(),
                                             mutsClonalMax=numeric(),
                                             mutsClonalSD=numeric(),
                                             stringsAsFactors=FALSE) 
mutTotalPerSampleSubset$n <- as.numeric(as.character(mutTotalPerSampleSubset$n))
for (patientID in patientsToUse){
  print(patientID)
  k <- max(mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID),]$n)
  for (n in 1:k){
    print(n)
    mutsClonalMean <- mean(mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$n == n),]$mutsClonal)
    mutsClonalSD <- sd(mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$n == n),]$mutsClonal)
    if (n > 1 & n < k){
      mutsClonalMin <- mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$n == n & mutTotalPerSampleSubset$tag == 'distMin'),]$mutsClonal
      mutsClonalMax <- mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$patientID == patientID & mutTotalPerSampleSubset$n == n & mutTotalPerSampleSubset$tag == 'distMax'),]$mutsClonal
    } else {
      mutsClonalMin <- mutsClonalMean
      mutsClonalMax <- mutsClonalMean
    }
    mutClonalMeanPerPatientPerN <- rbind.data.frame(mutClonalMeanPerPatientPerN, c(patientID, n, mutsClonalMean, mutsClonalMin, mutsClonalMax, mutsClonalSD), stringsAsFactors = F)
  }
}
colnames(mutClonalMeanPerPatientPerN) <- c('patientID','n','mutsClonalMean','mutsClonalMin','mutsClonalMax', 'mutsClonalSD')
mutTotalPerSampleSubset$n <- factor(as.character(mutTotalPerSampleSubset$n), levels = order(unique(as.integer(mutTotalPerSampleSubset$n))))

## plot to see how many samples are needed to know what is clonal
mutClonalMeanPerPatientPerN$n <- factor(as.character(mutClonalMeanPerPatientPerN$n), levels = order(unique(as.integer(mutClonalMeanPerPatientPerN$n))))
mutClonalMeanPerPatientPerN$patientID <- factor(as.character(mutClonalMeanPerPatientPerN$patientID), levels = patientsToUse)
mutClonalMeanPerPatientPerN$mutsClonalMean <- as.numeric(mutClonalMeanPerPatientPerN$mutsClonalMean)
mutClonalMeanPerPatientPerN$mutsClonalMin <- as.numeric(mutClonalMeanPerPatientPerN$mutsClonalMin)
mutClonalMeanPerPatientPerN$mutsClonalMax <- as.numeric(mutClonalMeanPerPatientPerN$mutsClonalMax)
mutClonalMeanPerPatientPerN$mutsClonalSD <- as.numeric(mutClonalMeanPerPatientPerN$mutsClonalSD)

ggplot(mutClonalMeanPerPatientPerN, aes(x=n, y=mutsClonalMean, group=patientID, colour=patientID)) +
  geom_point(aes(shape=patientID, colour=patientID)) +
  geom_line(aes(colour = patientID)) +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values=shapes) +
  geom_errorbar(aes(ymin=mutsClonalMean-mutsClonalSD, ymax=mutsClonalMean+mutsClonalSD), width=.2, position=position_dodge(0.05)) +
  theme(axis.text.x = element_text(size=10, color="black"),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank())

## Clonal by number of samples plot for each n relationship with distance
# this first part makes the data labels with p value and R2
mutTotalPerSampleSubset$meanDistance <- as.numeric(mutTotalPerSampleSubset$meanDistance)
mutTotalPerSampleSubsetUsableN <- mutTotalPerSampleSubset[which(mutTotalPerSampleSubset$n %in% c(2:7)),]
dataText <- data.frame(n=character(), p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 0 #where on the x axis will display p value
for (n in unique(mutTotalPerSampleSubsetUsableN$n)){
  print(n)
  y = 80
  for (i in rev(seq_along(patientsToUse))){
    patient=patientsToUse[i]
    print(patient)
    color <- colors[i]
    subset <- mutTotalPerSampleSubsetUsableN[which(mutTotalPerSampleSubsetUsableN$patientID==patient & mutTotalPerSampleSubsetUsableN$n == n),c('meanDistance','mutsClonal')]
    subset <- subset[complete.cases(subset),]
    meanDistance =  subset$meanDistance
    mutsClonal <- subset$mutsClonal
    if (length(meanDistance) > 1){
      testResult <- cor.test(meanDistance, mutsClonal, method="pearson")
      p=formatC(testResult$p.value,format = "e", digits = 2)
      R=round(testResult$estimate,3)
    } else {
      p=NA
      R=NA
    }
    label <- paste0('p=',p,', R=',R)
    dataText <- rbind(dataText, c(n,p,R,label,x,y,color,patient), stringsAsFactors=F)
    y <- y + 5
  }
}
colnames(dataText) <- c('n','p','R','label','x','y','color','patient')
dataText$adj.p <- NA
for (n in unique(dataText$n)){
  dataTextNSubset.index <- which(dataText$n == n & !is.na(dataText$p))
  dataText[dataTextNSubset.index,]$adj.p <- p.adjust(dataText[dataTextNSubset.index,]$p, method = "BH", n = length(dataText[dataTextNSubset.index,]$p))
}
dataText$label <- paste0('p=',dataText$adj.p,', R=',dataText$R)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'_mutsClonal_distance_n_stats.txt'),sep='\t', quote=F, row.names = F)
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
# plot
ggplot(mutTotalPerSampleSubsetUsableN, aes(x=meanDistance, y=mutsClonal, color=patientID)) +
  geom_point(aes(shape=patientID)) +
  scale_colour_manual(values=colors) +
  scale_shape_manual(values=shapes) +
  theme(axis.text.x = element_text(size=10, colour='black'), axis.title = element_text(size = 10), axis.text.y = element_text(size=10, colour='black'), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank(), legend.text=element_text(size=10), legend.title=element_text(size=10))+
  facet_wrap(vars(n)) +
  geom_smooth(method=lm,se=FALSE, fullrange=TRUE) +
  geom_text(
    data    = dataText,
    mapping = aes(x = x, y = y, label = label),
    hjust   = -0.01,
    colour=dataText$color
  )

## Calculation of change in number and percent of mutations per patient per N from a specified max sample number
maxN <- 7
mutClonalMeanPerPatientPerNDiffCalc <- mutClonalMeanPerPatientPerN[which(!mutClonalMeanPerPatientPerN$n==1),]
mutClonalMeanPerPatientPerNDiffCalc$diff <- NA
mutClonalMeanPerPatientPerNDiffCalc$diffPercent <- NA
for (patientID in patientsToUse){
  clonalEst <- mutClonalMeanPerPatientPerNDiffCalc[which(mutClonalMeanPerPatientPerNDiffCalc$patientID==patientID & mutClonalMeanPerPatientPerNDiffCalc$n==maxN),]$mutsClonalMean
  mutClonalMeanPerPatientPerNDiffCalc[which(mutClonalMeanPerPatientPerNDiffCalc$patientID==patientID),]$diff <- mutClonalMeanPerPatientPerNDiffCalc[which(mutClonalMeanPerPatientPerNDiffCalc$patientID==patientID),]$mutsClonalMean - clonalEst
  mutClonalMeanPerPatientPerNDiffCalc[which(mutClonalMeanPerPatientPerNDiffCalc$patientID==patientID),]$diffPercent <- 100*mutClonalMeanPerPatientPerNDiffCalc[which(mutClonalMeanPerPatientPerNDiffCalc$patientID==patientID),]$diff/clonalEst
}
mutClonalMeanPerPatientPerNDiffCalc$type <- 'LGG'
mutClonalMeanPerPatientPerNDiffCalc[which(mutClonalMeanPerPatientPerNDiffCalc$patientID %in% HGGSamples),]$type <- 'HGG'
mutClonalMeanPerPatientPerNDiffCalc$type <- factor(mutClonalMeanPerPatientPerNDiffCalc$type, levels=c("LGG","HGG"))
toPlot <- mutClonalMeanPerPatientPerNDiffCalc[which(mutClonalMeanPerPatientPerNDiffCalc$n %in% 1:6),]
ggplot(toPlot, aes(x=n, y=diff, fill=type)) +
  geom_boxplot(position=position_dodge()) +
  scale_fill_manual(values=colorsType) +
  theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'), legend.key=element_blank(), legend.text=element_text(size=15), legend.title=element_text(size=15))

