# Created: 2018.08.29
# By: Stephanie R Hilz
# Usage: For a list of patients, tags which mutations are tumor-wide and 
#   subclonal in that patient for a specified number of samples that 
#   meet a specified minimum tumor purity.

## Important note- Patient454 sample8 (v8) was excluded as upon phylogenetic analysis it did not contain mutations found in all other samples

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
outfolder <- 'SID000014_mutational_categorization_clonal_subclonal_by_patient/'
tag <- 'SID000014'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# specify cutoffs - which patients to look at, how many samples will be used, and how pure each must be;
#  if there are more samples than samplesPerPatientToUse, will average result from all possible combinations
#  of this number.
purityCutoff <- .3
samplesPerPatientToUse <- 2 #can also be NA; in this case all samples fitting the criteria will be used
requireSM <- FALSE

# read in sample data file
merged <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# subset by if has WES data
toUse <- data[which(!is.na(data$WES_ID)),]$WES_ID
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

# drop any samples excluded after manual inspection of phylogenetic trees
merged <- merged[which(!(merged$Patient=='P454' & merged$SampleName == 'v8')),]

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
SNVcategoryPerPatient <- data.frame(patient=character(),
                                           SNVuniqueID=character(),
                                           category=integer(),
                                           stringsAsFactors=FALSE) 

vafsPerSample <- data.frame(patient=character(),
                            SNVuniqueID=character(),
                            sampleID=character(),
                            vaf=numeric(),
                            stringsAsFactors=FALSE)

for (p in patientsToUse){
  patientID <- as.character(gsub('P','Patient',p))
  print(patientID)
  avfFile <- paste0(dataPath,patientID,'.R.mutations.avf.txt')
  
  # get vafs
  muts <- read.delim(avfFile, as.is=TRUE)
  if('X.gene' %in% colnames(muts)){
    geneIndex <- grep('X.gene',colnames(muts))
    colnames(muts)[geneIndex] <- 'gene'
  }
  muts <- muts[which(muts$algorithm == "MuTect"),]
  muts$uniq <- paste(muts$gene, muts$contig, muts$position, muts$ref_allele, muts$alt_allele, sep="_")
  rownames(muts) <- muts$uniq
  colnames(muts) <- gsub('[.]','-',colnames(muts))
  alt_counts <- muts[,grepl('_alt_Q20reads', colnames(muts))]
  colnames(alt_counts) <- gsub('_alt_Q20reads','',colnames(alt_counts))
  ref_counts <- muts[,grepl('_ref_Q20reads', colnames(muts))]
  colnames(ref_counts) <- gsub('_ref_Q20reads','',colnames(ref_counts))
  vafs <- alt_counts/(alt_counts+ref_counts)
  vafs$Normal <- NULL
  vafs$SNVuniqueID <- rownames(vafs)
  
  # melt data structure and add in patient info
  vafsPerSample_Local <- melt(vafs, id.vars='SNVuniqueID')
  colnames(vafsPerSample_Local) <- c('SNVuniqueID','sampleID','vaf')
  vafsPerSample_Local$patient <- p
  
  # run muts.bin function, which creates the same logical matrix used for phylo trees
  muts.bin <- getMutsBin(avfFile) #muts.bin for all samples in the patient, not yet subsetted
  
  # make the column names so they can string match our sample names
  colnames(muts.bin) <- gsub('[.]','-',colnames(muts.bin))
  colnames(muts.bin) <- gsub('_called','',colnames(muts.bin))
  
  # identify all samples with high enough purity for assessing what is tumor-wide for this patient
  minPuritySamples <- merged[which(merged$Patient == p),]$sample_type
  
  # determines which mutations are tumor-wide and which are not
  muts.bin.subset <- muts.bin[,minPuritySamples]
  tumorWideSNVs <- rownames(muts.bin.subset[which(rowSums(muts.bin.subset)== ncol(muts.bin.subset)),])
  regionalSNVs <- rownames(muts.bin.subset[which(rowSums(muts.bin.subset) < ncol(muts.bin.subset)),])
  category <- c(rep('tumor-wide',length(tumorWideSNVs)), rep('regional',length(regionalSNVs)))
  SNVcategoryPerPatient_Local <- as.data.frame(cbind(snvs=c(tumorWideSNVs, regionalSNVs), category=category))
  SNVcategoryPerPatient_Local$patient <- p
  
  # add to patient-level information about vafs for each sample to data frame
  vafsPerSample <- rbind(vafsPerSample, vafsPerSample_Local)
  
  # add to patient-level information about type of mutation to data frame
  SNVcategoryPerPatient <- rbind(SNVcategoryPerPatient, SNVcategoryPerPatient_Local)
}

## Output to files
cutofftag <- paste0('_',gsub('[.]','p',purityCutoff))
write.table(SNVcategoryPerPatient, file=paste0(outputPath,outfolder,tag,cutofftag,'_min_purity_',samplesPerPatientToUse,'_min_samples_snv_categories.txt'),sep='\t', quote=F, row.names = F)
write.table(vafsPerSample, file=paste0(outputPath,outfolder,tag,'_vafs.txt'),sep='\t', quote=F, row.names = F)

