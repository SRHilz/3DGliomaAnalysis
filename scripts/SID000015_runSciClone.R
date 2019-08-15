# Created: 2018.06.29
# By: Stephanie R Hilz
# Usage: Generates the input for SciClone for a patient then runs it. Requires patients to be sequenced very deeply OR have a large # of mutations (i.e. hypermutated)

#Dependencies to install if not installed already
source("http://bioconductor.org/biocLite.R")
biocLite("IRanges")
install.packages("devtools")

#Libraries to load (must do before next step)
library(devtools)

#Package to import from github using devtools
install_github("genome/bmm")
install_github("genome/sciClone")
library(sciClone)

#Functions
fixcn <- function(cp_input){
  cp_out <- data.frame(cp_input[,c(1,4,5,7)])
  colnames(cp_out) <- c('V1','V2','V3','V5')
  cp_out <- cp_out[which(!is.na(cp_out$V5)),]
  cp_out$V1 <- as.factor(cp_out$V1)
  cp_out$V2 <- as.integer(cp_out$V2)
  cp_out$V3 <- as.integer(cp_out$V3)
  return(cp_out)
}

generateSampleSpecificMutDF <- function(muts, sample){
  mutationsUnique <- paste(muts$gene, muts$contig, muts$position, muts$ref_allele, muts$alt_allele, sep="_")
  toReturn <- muts[,c('contig','position',paste0(sample,'_ref_Q20reads'),paste0(sample,'_alt_Q20reads'))]
  colnames(toReturn) <- c('chr','start','refCount', 'varCount')
  toReturn$VAF <- 100*(toReturn$varCount)/(toReturn$varCount + toReturn$refCount)
  rownames(toReturn) <- mutationsUnique
  return(toReturn)
}

generateSampleSpecifcCNDF <- function(cnFile){
  data <- readRDS(cnFile)
  cnData <- fixcn(data$output)
  return(cnData)
}

generateSampleSpecificLOH <- function(cnFile){
  data <- readRDS(cnFile)
  cutoff <- quantile(data$output$dhMean, na.rm=T)[2]+.1
  dataLOH <- data$output[which(data$output$dhMean > cutoff),c('chromosome','dhStart','dhEnd')]
  return(dataLOH)
}

# specify which patient to run SciClone on
patientID <- 'Patient452'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# file paths
tag <- 'SID000015'
CNDirectory <- paste0(dataPath,'CopyNumberSegmentationFiles/')
avfFile <- paste0(dataPath, patientID, '.R.mutations.avf.txt')
outputDir <- paste0(outputPath,tag,'_runSciClone')

# read in sample data file
sampleData <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# subset by patient
patientSampleData <- sampleData[which(sampleData$Patient == gsub('Patient','P',patientID)),]

# get out samples
samples <- patientSampleData$sample_type

# get CN files
files <- patientSampleData$CN

## Read in mutfile and do a bit of formatting to help make sure it meshes well with the copy number data format
muts <- read.delim(avfFile, as.is=TRUE)
if('X.gene' %in% colnames(muts)){
  geneIndex <- grep('X.gene',colnames(muts))
  colnames(muts)[geneIndex] <- 'gene'
}
colnames(muts) <- gsub('[.]','-',colnames(muts))
muts$contig <- gsub('chr','',muts$contig)
muts[which(muts$contig == 'X'),]$contig <- '23'

## Pull out relevant columns of mutfile
allowedChromosomes <- as.character(c(seq(1,22,1)))
muts <- muts[which(muts$algorithm == "MuTect" & muts$decision=="retain" & muts$contig %in% allowedChromosomes),]
annotation <- muts[,c('contig','position','gene')]

## Iterate through samples, making a sample-specific dataframe for each and appending it to a final list for mutations and copy number; also generates an analysis-wide list of blacklisted LOH regions
cnList <- vector("list", length(samples))
mutsList <- vector("list", length(samples))
regions_to_exclude <- data.frame(chromosome=integer(),
                                 dhStart=integer(), 
                                 dhEnd=integer(), 
                                 stringsAsFactors=FALSE) 
for (i in seq_along(samples)){
  mutsList[[i]] <- generateSampleSpecificMutDF(muts, samples[i])
  CNFile <- paste0(CNDirectory,files[i])
  cnList[[i]] <- generateSampleSpecifcCNDF(CNFile)
  regions_to_exclude <- rbind(regions_to_exclude, generateSampleSpecificLOH(CNFile))
}
regions_to_exclude <- unique(regions_to_exclude)

# Change this to 1 to output an intermediate plot at every iteration (as shown in supplement)
plotIntermediateResults <- 0

sc <- sciClone(vafs=mutsList, 
               sampleNames=samples, 
               annotation=annotation,
               maximumClusters = 8, 
               useSexChrs=FALSE, 
               copyNumberCalls=cnList, 
               copyNumberMargins=0.25, 
               minimumDepth=50, 
               verbose=1, 
               doClusteringAlongMargins=F, 
               regionsToExclude = regions_to_exclude,
               plotIntermediateResults = plotIntermediateResults)


writeClusterTable(sc,paste0(outputDir, "medial.cluster.tsv"))
writeClusterSummaryTable(sc,paste0(outputDir,"medial.cluster.summary.tsv"))
sc.plot2d(sc,paste0(outputDir, '2d_clustering_plots_medial.pdf'), scale=1.8)

connectivity.matrix <- getConnectivityMatrix(sc)

write.table(file=paste0(outputDir, "medial.connectivity.matrix.tsv"), connectivity.matrix, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

