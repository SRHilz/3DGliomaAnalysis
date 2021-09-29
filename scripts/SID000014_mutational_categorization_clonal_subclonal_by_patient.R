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
exploratoryPurityCutoff <- .3 # used at the very beginging just to get a sense of what is tumor-wide
samplesPerPatientToUse <- 2 #can also be NA; in this case all samples fitting the criteria will be used
pycloneCellPrevCutoff <- 0
requireSM <- FALSE

# read in sample data file
sampleData <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# subset by if has WES data
toUse <- sampleData[which(!is.na(sampleData$WES_ID)),]$WES_ID
sampleData <- sampleData[which(sampleData$WES_ID %in% toUse),]

# subset by sample type if desired
if (requireSM){
  sampleData <- sampleData[which(sampleData$SampleType=='SM'),]
}

#######################################################################################
### PART1 - Exploratory (only used for evaluating PyClone and FACETS purity)
# set up structure just for exploratory
merged <- sampleData 
 
# specify purity metric to use
merged$purity <- merged$PyClone

# Filter based on purity and inspect numbers of samples per patient
merged <- merged[which(merged$purity >= exploratoryPurityCutoff),]
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
  avfFile <- file.path(dataPath,'mutations_avf',paste0(patientID,'.R.mutations.avf.txt'))
  
  # get vafs
  muts <- read.delim(avfFile, as.is=TRUE)
  muts <- muts[muts$decision=='retain',]
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
  minPuritySamples <- merged[which(merged$Patient == p & merged$purity >= exploratoryPurityCutoff),]$sample_type
  
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
cutofftag <- paste0('_',gsub('[.]','p',exploratoryPurityCutoff))
write.table(SNVcategoryPerPatient, file=paste0(outputPath,outfolder,tag,cutofftag,'_min_purity_',samplesPerPatientToUse,'_min_samples_snv_categories.txt'),sep='\t', quote=F, row.names = F)
write.table(vafsPerSample, file=paste0(outputPath,outfolder,tag,'_vafs.txt'),sep='\t', quote=F, row.names = F)

#######################################################################################

### PART 2: Used for all other mutation classification in the manuscript, which are done using PyClone

# specify purity metric to use
sampleData$purity <- sampleData$PyClone
sampleData[sampleData$PurityEstUsed=='FACETS',]$purity <- sampleData[sampleData$PurityEstUsed=='FACETS',]$FACETS

# Set final patients to use
patientsToUse <- patientOrder

# Read in pyclone clone cellular prevalence data 
pycloneClone <- read.table(file.path(dataPath, 'pycloneiv_100rr_40cl_results_cluster_abundance.tsv'), sep='\t', header=T)

# Read in pyclone variant to clone data
pycloneVar <- read.table(file.path(dataPath, 'pycloneiv_100rr_40cl_results_variant_to_cluster.tsv'), sep='\t', header=T)

# Create an empty data frame for storing the number of total mutations called per patient that are clonal, shared, or private
vafsPerSample <- c()
colNames <- c('patient','sampleID','SNVuniqueID','type','SangerCancerGeneCensus-','vaf','mut_called','mut_category','mut_category_used','mut_category_nsamp','cluster_id','cellular_prevalence','cellular_prevalence_std','upper_95CI','lower_95CI','pyclone_called','clonal','pyclone_category','pyclone_category_used','pyclone_category_nsamp','cluster_assignment_prob')

for (p in patientsToUse){
  patientID <- as.character(gsub('P','Patient',p))
  print(patientID)
  avfFile <- file.path(dataPath,'mutations_avf',paste0(patientID,'.R.mutations.avf.txt'))
  
  # get vafs
  muts <- read.delim(avfFile, as.is=TRUE)
  muts <- muts[muts$decision=='retain',]
  if('X.gene' %in% colnames(muts)){
    geneIndex <- grep('X.gene',colnames(muts))
    colnames(muts)[geneIndex] <- 'gene'
  }
  muts <- muts[which(muts$algorithm == "MuTect"),]
  muts$SNVuniqueID <- paste(muts$gene, gsub('chr','',muts$contig), muts$position, muts$ref_allele, muts$alt_allele, sep="_")
  rownames(muts) <- muts$SNVuniqueID
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
  
  # make muts.bin SNV unique ID compatible downstream
  rownames(muts.bin) <- gsub('_chr','_',rownames(muts.bin))
  
  # make the column names so they can string match our sample names
  colnames(muts.bin) <- gsub('[.]','-',colnames(muts.bin))
  colnames(muts.bin) <- gsub('_called','',colnames(muts.bin))
  
  # identify all samples with high enough purity for assessing what is tumor-wide for this patient
  minPurityCallsSamples <- sampleData[which(sampleData$Patient == p & sampleData$purity >= finalPurityCutoffCalls),]$sample_type
  
  # determines which mutations are tumor-wide and which are not
  muts.bin.subset <- muts.bin[,minPurityCallsSamples]
  tumorWideSNVs <- rownames(muts.bin.subset[which(rowSums(muts.bin.subset)== ncol(muts.bin.subset)),])
  privateSNVs <- rownames(muts.bin.subset[which(rowSums(muts.bin.subset)== 1),])
  sharedSNVs <- rownames(muts.bin.subset[which(rowSums(muts.bin.subset) < ncol(muts.bin.subset) & rowSums(muts.bin.subset) > 1),])
  category <- c(rep('tumor-wide',length(tumorWideSNVs)), rep('shared',length(sharedSNVs)), rep('private',length(privateSNVs)))
  mut_category_nsamp <- ncol(muts.bin.subset)
  SNVcategoryPerPatient_Local <- as.data.frame(cbind(SNVuniqueID=c(tumorWideSNVs, sharedSNVs, privateSNVs), mut_category=category, mut_category_nsamp=mut_category_nsamp))
  
  # melt muts.bin to get call
  muts.bin$SNVuniqueID <- row.names(muts.bin)
  muts.binPerSample_Local <- melt(muts.bin, id.vars='SNVuniqueID')
  colnames(muts.binPerSample_Local) <- c('SNVuniqueID','sampleID','mut_called')
  
  # add distribution of calls to call (mut-level)
  muts.binPerSample_Local <- merge(SNVcategoryPerPatient_Local, muts.binPerSample_Local, by='SNVuniqueID', all=TRUE)
  
  # add call information to vafs
  vafsPerSample_Local$SNVuniqueIDSamp <- paste0(vafsPerSample_Local$SNVuniqueID,'_',vafsPerSample_Local$sampleID)
  muts.binPerSample_Local$SNVuniqueIDSamp <- paste0(muts.binPerSample_Local$SNVuniqueID,'_',muts.binPerSample_Local$sampleID)
  vafsPerSample_Local <- merge(vafsPerSample_Local, muts.binPerSample_Local[,c('SNVuniqueIDSamp','mut_called','mut_category','mut_category_nsamp')], by='SNVuniqueIDSamp', all=T)
  
  # annotate whether sample calls were used to decide to determine mut_category
  vafsPerSample_Local$mut_category_used <- ifelse(vafsPerSample_Local$sampleID %in% minPurityCallsSamples, TRUE, FALSE)
  
  # add sanger cancer annotation
  vafsPerSample_Local <- merge(vafsPerSample_Local, muts[,c('SNVuniqueID','SangerCancerGeneCensus-','type')], by='SNVuniqueID', all=T)

  # bring in pyclone clonal abundance per sample
  pycloneClone_Local <- pycloneClone[pycloneClone$Patient==patientID,]
  
  # determine whether prevalence  @ lower 95CI is over threshold
  pycloneClone_Local$pyclone_called <- pycloneClone_Local$lower_95CI > pycloneCellPrevCutoff
  
  # create call matrix
  pyclone.bin <- spread(pycloneClone_Local[,c('sample_id','cluster_id','pyclone_called')], sample_id, pyclone_called )
  rownames(pyclone.bin) <- pyclone.bin$cluster_id
  pyclone.bin$cluster_id <- NULL
  
  # subset to remove low purity samples if needed
  minPurityNonzeroSamples <- sampleData[which(sampleData$Patient == p & sampleData$purity >= finalPurityCutoffNonzero),]$sample_type
  pyclone.bin.subset <- pyclone.bin[,minPurityNonzeroSamples] # currently using all
  tumorWideClusters <- rownames(pyclone.bin.subset[which(rowSums(pyclone.bin.subset)== ncol(pyclone.bin.subset)),])
  privateClusters <- rownames(pyclone.bin.subset[which(rowSums(pyclone.bin.subset)== 1),])
  sharedClusters <- rownames(pyclone.bin.subset[which(rowSums(pyclone.bin.subset) < ncol(pyclone.bin.subset) & rowSums(pyclone.bin.subset) > 1),])
  category <- c(rep('tumor-wide',length(tumorWideClusters)), rep('shared',length(sharedClusters)), rep('private',length(privateClusters)))
  pyclone_category_nsamp <- ncol(pyclone.bin.subset)
  PyClonePerPatient_Local <- as.data.frame(cbind(cluster_id=c(tumorWideClusters, sharedClusters, privateClusters), pyclone_category=category, pyclone_category_nsamp=pyclone_category_nsamp))

  # add categories to original cluster structure
  pycloneClone_Local <- merge(PyClonePerPatient_Local, pycloneClone_Local, by='cluster_id', all=TRUE)
  
  # bring in pyclone variant to clone mapping
  pycloneVar_Local <- pycloneVar[pycloneVar$Patient==patientID,]
  
  # add variant clone mapping to clonal abundace
  pycloneAll_Local <- merge(pycloneClone_Local[,c('cluster_id','sample_id','cellular_prevalence','cellular_prevalence_std','upper_95CI','lower_95CI','pyclone_called','clonal','pyclone_category','pyclone_category_nsamp')], pycloneVar_Local[,c('mutation_id','cluster_id','cluster_assignment_prob')], by='cluster_id', all=T)
  
  # annotate whether sample calls were used to decide pyclone_category
  pycloneAll_Local$pyclone_category_used <- ifelse(pycloneAll_Local$sample_id %in% minPurityNonzeroSamples, TRUE, FALSE)
  
  # add pyclone to vafs
  pycloneAll_Local$SNVuniqueIDSamp <- paste0(pycloneAll_Local$mutation_id,'_',pycloneAll_Local$sample_id)
  pycloneAll_Local$sample_id <- NULL
  pycloneAll_Local$mutation_id <- NULL
  vafsPerSample_Local <- merge(vafsPerSample_Local, pycloneAll_Local, by='SNVuniqueIDSamp', all.x=T, all.y=F) # we want everything with a vaf, but nothing in PyClone not in VAFs 
  
  # clean up final data structure
  vafsPerSample_Local$SNVuniqueIDSamp <- NULL
  vafsPerSample_Local <- vafsPerSample_Local[,colNames]
  
  # remove any samples not in the final sample set
  sampleData_Local <- sampleData[sampleData$Patient==p,]
  vafsPerSample_Local <- vafsPerSample_Local[vafsPerSample_Local$sampleID %in% sampleData_Local$sample_type,]
  
  # add to patient-level information about vafs for each sample to data frame
  vafsPerSample <- rbind(vafsPerSample, vafsPerSample_Local)
  
}

# Check all looks as expected
sum(is.na(vafsPerSample$patient))
sum(is.na(vafsPerSample$sampleID))
sum(is.na(vafsPerSample$SNVuniqueID))
sum(is.na(vafsPerSample$mut_called)) # some will be NA because not covered in all, is okay
sum(is.na(vafsPerSample$vaf))
sum(is.na(vafsPerSample$pyclone_called)) # many will be NA here because there are variants that couldn't be used by PyClone

## Output to files
cutofftag <- paste0('_',gsub('[.]','p',finalPurityCutoff))
write.table(vafsPerSample, file=paste0(outputPath,outfolder,tag,cutofftag,'_min_final_purity_samples_snv_vafs_categories.txt'),sep='\t', quote=F, row.names = F)



