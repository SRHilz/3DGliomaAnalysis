# Created: 2018.10.23
# By: Stephanie R Hilz
# Usage: For a specific patient, pulls out shared mutations (historically this code did much more and actually id shared clusters, but now all of this is done in SID000014)

library(ggplot2)
library(kableExtra)
library(reshape2)
library(gplots)
library(gdata)
library(gtools)
library(grDevices)
library(vegan)

#user-defined variables
outfolder <- 'SID000020_lobular_vs_intermixed_characteristics_analysis/'
tag <- 'SID000020'
infolder_clusterdiff <- 'SID000011_mutational_analysis_shared_group_expansion_analysis/'

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# specify cutoffs - which patients to look at, how many samples will be used
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

# subset to just shared
vafsPerSample <- vafsPerSample[vafsPerSample$pyclone_category=='shared',]

# Remove patients with less than the required number of high-CCF samples and reinspect
patientsToUse <- unique(vafsPerSample[which(vafsPerSample$pyclone_category_nsamp >= samplesPerPatientToUse),]$patient)

# Also remove patients not in merged (should be none)
patientsToUse <- patientsToUse[patientsToUse %in% merged$Patient]

# Set final patients to use
patientsToUse <- patientOrder[patientOrder %in% patientsToUse]

# Subset our categories to this number
vafsPerSample_Filtered <- vafsPerSample[vafsPerSample$patient %in% patientsToUse,]

# Bring in annotation of clusters as lobular or intermixed
infile_clusterdiff=paste0(outputPath, infolder_clusterdiff, 'clusters_all_patients_post_filtering_final_spatial_diff.txt')
clusterDiff <- read.table(infile_clusterdiff, sep='\t', header=T)

# Make common cluster unique ID and add this annotation to our clusters
vafsPerSample_Filtered$patientCluster <- paste0(vafsPerSample_Filtered$patient,'-',vafsPerSample_Filtered$cluster_id)
vafsPerSample_Filtered <- merge(vafsPerSample_Filtered, clusterDiff[,c('patientCluster','cloneDistribution')], by='patientCluster',all.x=TRUE)

# Create matrix of mutations for oncoprint by clone with patient and distribution information
vafsPerSample_mutsByClone <- vafsPerSample_Filtered[,c('patientCluster','patient','SNVuniqueID','type','SangerCancerGeneCensus.','cluster_id','cloneDistribution')]
vafsPerSample_mutsByClone <- vafsPerSample_mutsByClone[!duplicated(vafsPerSample_mutsByClone),]
vafsPerSample_mutsByClone$gene <- sapply(vafsPerSample_mutsByClone$SNVuniqueID, function(x) strsplit(x, '_')[[1]][1])
vafsPerSample_mutsByClone[is.na(vafsPerSample_mutsByClone$type),]$type <- 'Splicing' # mutations that are splicing get type NA, so we correct this here
# make matrix of mutated genes in each clone
mat <- c()
rowNames <- unique(vafsPerSample_mutsByClone$gene)
colNames <- unique(vafsPerSample_mutsByClone$patientCluster)
for (g in rowNames){
  toBind <- c()
  for (c in colNames){
    muts <- paste(unique(vafsPerSample_mutsByClone[vafsPerSample_mutsByClone$gene==g & vafsPerSample_mutsByClone$patientCluster==c,]$type), collapse=',')
    toBind <- append(toBind, muts)
  }
  mat <- rbind(mat, toBind)
}
rownames(mat) <- rowNames
colnames(mat) <- colNames
mat_noP452 <- mat[,colnames(mat)[!grepl('P452',colnames(mat))]]
genesOfInterest <- rownames(mat_noP452[rowSums(!mat_noP452=='') > 2,]) #we don't include 452 because evertthing is mutated in this patient
mat_subset <- mat[genesOfInterest,]
clonesOfInterest <- colnames(mat_subset[,colSums(!mat_subset=='') > 0])
mat_subset <- mat_subset[,clonesOfInterest]
# set up alter function for mat
col <- c(Missense="orange",Nonsense="red",`Non-stop`="pink",Splicing="magenta",Silent="green",unknown="goldenrod")
alter_fun <- list(
  background = function(x,y,w,h) grid.rect(x,y,w*0.9,h*0.9,gp=gpar(fill='#e0e0e0',col=NA)),
  Missense = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.9, gp=gpar(fill=col["Missense"], col = NA)),
  Nonsense = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.8, gp=gpar(fill=col["Nonsense"], col = NA)),
  `Non-stop` = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.7, gp=gpar(fill=col["Non-stop"], col = NA)),
  Splicing = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.6, gp=gpar(fill=col["Splicing"], col = NA)),
  Silent = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.3, gp=gpar(fill=col["Silent"], col = NA)),
  unknown = function(x,y,w,h) grid.rect(x,y,w*0.9, h*0.3, gp=gpar(fill=col["unknown"], col = NA))
)
# make top annotation df (clone annotation)
top_df <- vafsPerSample_mutsByClone[,c('patientCluster','patient','cloneDistribution')]
top_df <- top_df[!duplicated(top_df),]
rownames(top_df) <- top_df[,'patientCluster']
top_df$patientCluster <- NULL
top_df <- top_df[clonesOfInterest,]
top_annotation <- HeatmapAnnotation(df=top_df)
# make side annotation df (gene annotation)
left_df <- vafsPerSample_mutsByClone[,c('gene','SangerCancerGeneCensus.')]
left_df <- left_df[!duplicated(left_df),]
rownames(left_df) <- left_df[,'gene']
left_df$gene <- NULL
left_annotation <- rowAnnotation(df=left_df[genesOfInterest,])
# print oncoprint
colOrder <- rownames(top_df[order(top_df$cloneDistribution),])
oncoPrint(mat_subset, 
          alter_fun=alter_fun, 
          col=col,
          column_split = top_df$cloneDistribution,
          top_annotation = top_annotation,
          left_annotation = left_annotation,
          show_column_names = FALSE)

# test association between SangerCensus and intermixed/lobular by gene (proportion lobular per gene)
set.seed(1234)
toTest <- vafsPerSample_mutsByClone[,c('gene','patientCluster','patient','SangerCancerGeneCensus.','cloneDistribution')]
toTest <- toTest[!duplicated(toTest),]
toTestFinal <- data.frame(gene=character(), `SangerCancerGeneCensus.`=character(), proportionLobular=numeric(), stringsAsFactors = F)
for (g in unique(toTest$gene)){
  local <- toTest[toTest$gene==g,c('gene','patient','SangerCancerGeneCensus.','cloneDistribution')]
  status <- unique(local$SangerCancerGeneCensus.)
  if (duplicated(local$patient)){
    dim(local)
    print(p)
    for (p in unique(local[duplicated(local$patient),]$patient)){
      print(g)
      local_p_dist <- local[local$patient==p,]$cloneDistribution
      print(local_p_dist)
      local_p_dist_random <- sample(local_p_dist, 1)
      local[local$patient==p,]$cloneDistribution <- local_p_dist_random
    }
    local <- local[!duplicated(local),]
    dim(local)
  }
  proportionLobular <- nrow(local[local$cloneDistribution=='lobular',])/nrow(local)
  toTestFinal <- rbind(toTestFinal, c(g,status,proportionLobular))
}
colnames(toTestFinal) <- c('gene','SangerCancerGeneCensus.','proportionLobular')
toTestFinal$proportionLobular <- as.numeric(toTestFinal$proportionLobular)
ggboxplot(toTestFinal, x = "SangerCancerGeneCensus.", y = "proportionLobular",
               color="SangerCancerGeneCensus.", palette = "jco",
               add = "jitter") +
  stat_compare_means(method = "wilcox") #proportionLobular here is not within a patient but across patients
aggregate(proportionLobular~SangerCancerGeneCensus., data=toTestFinal, mean)

# test association between SangerCensus and intermixed/lobular by patient (proportion of clusters with at least one sanger for lobular vs intermixed)
toTestNoP452 <- toTest[!toTest$patient=='P452',]
toTestNoP452$subtype <- ifelse(toTestNoP452$patient %in% subsetIDHwtc, 'IDHwtc','IDHmut')
toTestFinal <- data.frame(patient=character(), subtype=character(), cloneDistribution=character(),proportionMinOneSanger=character(), stringsAsFactors = F)
for (p in unique(toTestNoP452$patient)){
  local <- toTestNoP452[toTestNoP452$patient==p,]
  subtype <- unique(local$subtype)
  lobular <- local[local$cloneDistribution=='lobular',]
  intermixed <- local[local$cloneDistribution=='intermixed',]
  lobular.total = NA # number of total clusters
  lobular.sanger = 0 # number of clusters with at lease one sanger gene
  intermixed.total = NA # number of total clusters
  intermixed.sanger = 0 # number of clusters with at lease one sanger gene
  for (c in unique(lobular$patientCluster)){
    if (is.na(lobular.total)){
      lobular.total = 0 # initialize it
    }
    lobular.total = lobular.total + 1
    if(any(lobular$patientCluster==c & lobular$SangerCancerGeneCensus.=='YES')){
      lobular.sanger = lobular.sanger + 1
    }
  }
  for (c in unique(intermixed$patientCluster)){
    if (is.na(intermixed.total)){
      intermixed.total = 0 # initialize it
    }
    intermixed.total = intermixed.total + 1
    if(any(intermixed$patientCluster==c & intermixed$SangerCancerGeneCensus.=='YES')){
      intermixed.sanger = intermixed.sanger + 1
    }
  }
  proportionLobularSanger <- lobular.sanger/lobular.total
  proportionIntermixedSanger <- intermixed.sanger/intermixed.total
  toTestFinal <- rbind(toTestFinal, c(p,subtype,'lobular',proportionLobularSanger))
  toTestFinal <- rbind(toTestFinal, c(p,subtype,'intermixed',proportionIntermixedSanger))
}
colnames(toTestFinal) <- c('patient','subtype','cloneDistribution','proportionSanger')
toTestFinal$proportionSanger <- as.numeric(toTestFinal$proportionSanger)
ggboxplot(toTestFinal, x = "cloneDistribution", y = "proportionSanger",
          color="cloneDistribution", palette = "jco",
          add = "jitter", facet.by = "subtype") +
  stat_compare_means(method = "wilcox") #proportionLobular here is not within a patient but across patients
aggregate(proportionLobular~SangerCancerGeneCensus., data=toTestFinal, mean)


# Create df of samples with max lobular and intermixed clones per patient
toCompare <- vafsPerSample_Filtered[,c('patient','patientCluster','sampleID','cellular_prevalence','cloneDistribution')]
toCompare <- toCompare[!duplicated(toCompare),] # remove gene-level info
toCompare_final <- data.frame(patient=character(), sample=character(), cloneDistribution=character(), cellular_prevalence=numeric(), stringsAsFactors = F)
for (p in unique(toCompare$patient)){
  toCompare_Patient <- toCompare[toCompare$patient==p,]
  intermixed_max <- max(toCompare_Patient[toCompare_Patient$cloneDistribution=='intermixed',]$cellular_prevalence)
  if (!is.infinite(intermixed_max)){
    intermixed_sample <- toCompare_Patient[toCompare_Patient$cellular_prevalence==intermixed_max & toCompare_Patient$cloneDistribution=='intermixed',]$sampleID
    toCompare_final <- rbind(toCompare_final, c(p,intermixed_sample,'intermixed',intermixed_max))
  }
  lobular_max <- max(toCompare_Patient[toCompare_Patient$cloneDistribution=='lobular',]$cellular_prevalence)
  if (!is.infinite(lobular_max)){
    lobular_sample <- toCompare_Patient[toCompare_Patient$cellular_prevalence==lobular_max & toCompare_Patient$cloneDistribution=='lobular',]$sampleID
    toCompare_final <- rbind(toCompare_final, c(p,lobular_sample,'lobular',lobular_max))
  }
}
colnames(toCompare_final) <- c('patient','sampleID','cloneDistribution','cellular_prevalence')
toCompare_final <- toCompare_final[toCompare_final$cellular_prevalence >= .5,]
toCompare_final <- toCompare_final[!toCompare_final$patient=='P452',]
toCompare_final$subtype <- ifelse(toCompare_final$patient %in% subsetIDHwtc,'IDH-wtc','IDH-mut')
# write
write.table(toCompare_final, file=paste0(outputPath, outfolder, 'metadata_for_rnaseq_comparison_lobular_vs_intermixed.txt'), quote=F, row.names=F, sep='\t')

  