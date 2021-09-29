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

calcWithinGroupDistancePyclone <- function(df, k, distances){
  #VAFs and distances should already be reshaped so that any replicates are averaged, and any missing samples are removed. Should both contain the same samples with the same names
  df_cluster <- df[df$cluster==k,]
  within <- df_cluster[which(df_cluster$pyclone_called),]$sample
  all <- df_cluster$sample
  n <- length(within)
  distances[lower.tri(distances)] <- NA
  diag(distances) <- NA
  pairwiseDistances <- unmatrix(distances)
  pairwiseDistances <- pairwiseDistances[!is.na(pairwiseDistances)]
  withinDistancesIndex <- c()
  for (i in within){
    for (j in within){
      for (l in 1:length(pairwiseDistances)){
        if (paste(i,j,sep=':')==names(pairwiseDistances)[l] | paste(j,i,sep=':')==names(pairwiseDistances)[l]){
          withinDistancesIndex <- append(withinDistancesIndex,l)
        } 
      }
    }
  }
  withinDistances <- pairwiseDistances[unique(withinDistancesIndex)]
  withinDistancesMean <- mean(withinDistances)
  randomDistancesMeans <- c()
  combinations <- combn(all,n)
  for (c in 1:dim(combinations)[2]){
    if (c %% 100 == 0){
      percent_done <- 100*(c/dim(combinations)[2])
      print(paste0('still on ',unique(df$patient),' cluster ',k,'.  ',percent_done , '% done'))
    }
    combinationIndexes <- c()
    for (i in combinations[,c]){
      for (j in combinations[,c]){
        for (l in 1:length(pairwiseDistances)){
          if (paste(i,j,sep=':')==names(pairwiseDistances)[l] | paste(i,j,sep=':')==names(pairwiseDistances)[l]){
            combinationIndexes <- append(combinationIndexes,l)
          }
        }
      }
    }
    randomDistancesMeans <- append(randomDistancesMeans, mean(pairwiseDistances[unique(combinationIndexes)]))
  }
  #   if(k==5 & vafThreshold==.1){
  #     print(randomDistancesMeans)
  #     print(n)
  #     print(vafThreshold)
  #     hist(randomDistancesMeans,col='blue')
  #   }
  distancesMean <- mean(randomDistancesMeans)
  p <- round(length(randomDistancesMeans[which(randomDistancesMeans <= withinDistancesMean)])/length(randomDistancesMeans),2)
  toReturn <- list(within=withinDistancesMean, random=randomDistancesMeans, p=p, n=n, samples=within)
  return(toReturn)
}

#user-defined variables
outfolder <- 'SID000011_mutational_analysis_shared_group_expansion_analysis/'
tag <- 'SID000011'

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

# Remove patients with less than the required number of high-CCF samples and reinspect
patientsToUse <- unique(vafsPerSample[which(vafsPerSample$pyclone_category_nsamp >= samplesPerPatientToUse),]$patient)

# Also remove patients not in merged (should be none)
patientsToUse <- patientsToUse[patientsToUse %in% merged$Patient]

# Set final patients to use
patientsToUse <- patientOrder[patientOrder %in% patientsToUse]

# Subset our categories to this number
vafsPerSample_Filtered <- vafsPerSample[vafsPerSample$patient %in% patientsToUse,]

toOutput <- data.frame(patient=character(),
                                 cluster=character(),
                                 sample=character(),
                                 pyclone_called=integer(),
                                 stringsAsFactors=FALSE) 

for (p in patientsToUse){
  print(p)
  # subset variants to patient
  vafsPerSample_Filtered_subset <- vafsPerSample_Filtered[vafsPerSample_Filtered$patient==p,]
  
  # subset to shared
  vafsPerSample_Filtered_subset <- vafsPerSample_Filtered_subset[which(vafsPerSample_Filtered_subset$pyclone_category=='shared'),]
  
  # condense to just clone info
  vafsPerSample_Filtered_subset_clones <- vafsPerSample_Filtered_subset[,c('patient','cluster_id','sampleID','pyclone_called','pyclone_category_used')]
  vafsPerSample_Filtered_subset_clones <- vafsPerSample_Filtered_subset_clones[!duplicated(vafsPerSample_Filtered_subset_clones),]
  
  # filter out samples not used due to low purity
  vafsPerSample_Filtered_subset_clones <- vafsPerSample_Filtered_subset_clones[which(vafsPerSample_Filtered_subset_clones$pyclone_category_used),]
  vafsPerSample_Filtered_subset_clones$pyclone_category_used <- NULL
  
  # update colnames
  colnames(vafsPerSample_Filtered_subset_clones) <- c('patient','cluster','sample','pyclone_called')
  
  # create a distance matrix
  coordinates <- merged[which(merged$Patient==p),c('L.Coordinate','P.Coordinate','S.Coordinate')]
  rownames(coordinates) <- merged[which(merged$Patient==p),'sample_type']
  exomeDistanceMatrix <- as.matrix(dist(coordinates, method = "euclidean"))
  
  # Subset distance matrix to correspond with vafs (since vafs come from a file already filtered for purity, should already have low-purity samples removed)
  toUse <- mixedsort(unique(vafsPerSample_Filtered_subset_clones$sample))
  reshapedDistances <- exomeDistanceMatrix[toUse,toUse]
  
  #Main
  outputRandom <- c()
  outputSummary <- c()
  for (k in unique(vafsPerSample_Filtered_subset_clones$cluster)){
    print(k)
    y <- calcWithinGroupDistancePyclone(vafsPerSample_Filtered_subset_clones, k, reshapedDistances)
    #print(y)
    meanRandom <- mean(y$random)
    outputRandom <- rbind(outputRandom,cbind(k,y$random))
    outputSummary <- rbind(outputSummary,c(k,y$within,meanRandom,y$p,y$n, paste(y$samples,collapse=',')))
  }
  colnames(outputRandom)<-c('cluster','meanDistance')
  colnames(outputSummary) <- c('cluster','meanDistanceWithin', 'meanDistanceRandom','p','n','samples')
  outputRandom <- as.data.frame(outputRandom)
  outputRandom$meanDistance <- as.numeric(outputRandom$meanDistance)
  outputRandom$cluster <- factor(outputRandom$cluster, levels = sort(unique(outputRandom$cluster)))
  outputSummary <- as.data.frame(outputSummary)
  outputSummary$meanDistanceWithin <- as.numeric(as.character(outputSummary$meanDistanceWithin))
  outputSummary$meanDistanceRandom <- as.numeric(as.character(outputSummary$meanDistanceRandom))
  outputSummary$cluster <- factor(outputSummary$cluster, levels = sort(unique(outputRandom$cluster)))
  ggplot(data = outputRandom, aes(y=meanDistance, x=cluster)) + 
    labs(y="Average pairwise distance (mm)") + 
    geom_boxplot(position="dodge") + 
    geom_point(data=outputSummary, aes(y=meanDistanceWithin, x=cluster, size=4), position=position_dodge(width=0.8), color='red', shape=17) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
  ggsave(file=paste0(outputPath, outfolder, p,'_expected_vs_observed_distances.pdf'))
  saveRDS(outputSummary, file=paste0(outputPath, outfolder, p,'_outputSummary'))
  write.table(outputSummary, file=paste0(outputPath, outfolder, p,'_outputSummary.txt'), quote=F,sep='\t', row.names=F)
  
  
  # add to patient-level information about type of mutation data frame
  toOutput <- rbind(toOutput, vafsPerSample_Filtered_subset_clones)
}

write.table(toOutput, file=paste0(outputPath, outfolder, 'clusters_all_patients_post_filtering_presence_absence.txt'), quote=F, row.names=F, sep='\t')

# get colors for vaf plotting (note - this part requires the vafs object to already exist from above for the correct patient - these colors can be fed into periphery.R to automate visualization, or entered manually into Slicer)
# note - I typically use #DADEDF (218 222 223) to designate samples that are excluded from the analysis
cluster <- read.table(paste0(outputPath, outfolder, patientID, '_clusters_post_filtering.txt'), sep='\t', header=F)
cluster_vaf_means <- list()
for (c in unique(cluster$V2)){
  print(c)
  snvs <- as.character(cluster[cluster$V2==c,]$V1)
  cluster_vaf_means[[c]] <- colMeans(vafs[snvs,])
}
toColor <- cluster_vaf_means[[1]]#change this index to specify a different cluster
rbPal <- colorRampPalette(c("blue","red"))
mappedColors <- rbPal(length(toColor))[as.numeric(cut(toColor,breaks = length(toColor)))][1:length(toColor)]
names(mappedColors) <- names(toColor)
mappedColors <- mappedColors[toUse]
#for slicer
col2rgb(mappedColors)
#make gradient for colorbar
par(mfrow=c(1,1))
plot(rep(1,10), col=rbPal(10),pch=15, cex=3, xlim=c(1,40), axes=F, ann=F)

## PART 2 - only run this code once you have generated everything above for all patients you use here
# Final aggregate figure 
toPlotFinal <- c()
for (patientID in patientsToUse){
  file=paste0(outputPath, outfolder, patientID,'_outputSummary')
  data <- readRDS(file)
  data$patient <- patientID
  toPlotFinal <- rbind(toPlotFinal,data)
}
# annotate patient colors
toPlotFinal$color <- NA
for (p in unique(toPlotFinal$patient)){
  toPlotFinal[which(toPlotFinal$patient == p),]$color <- as.character(colorKey[patient_ID])
}
toPlotFinal$patientCluster <- paste0(gsub("Patient",'P',toPlotFinal$patient),'-',toPlotFinal$cluster)
toPlotFinal$patient <- factor(toPlotFinal$patient, levels=patientsToUse)
toPlotFinal$type <- 'Other'
toPlotFinal[which(toPlotFinal$patient %in% subsetIDHwtc),]$type <- 'IDHwtTERTp'
# plot by difference 
toPlotFinal$difference <- toPlotFinal$meanDistanceWithin-toPlotFinal$meanDistanceRandom
colors <- as.character(colorKey[gsub('Patient','P',patientOrder)])
# Look at proportion of clusters per patient that are < random
proportionsToPlot <- data.frame(patient=character(), proportion=numeric(), subtype=character(), ordinaltype=character(), stringsAsFactors=F)
for (p in unique(toPlotFinal$patient)){
  toPlotFinalPatient <- toPlotFinal[which(toPlotFinal$patient==p),]
  lessThanThreshold <- nrow(toPlotFinalPatient[which(toPlotFinalPatient$difference < 0),])
  totaln <- nrow(toPlotFinalPatient)
  proportion <- lessThanThreshold/totaln
  subtype <- 'Other'
  if (p %in% subsetIDHwtc){subtype <- 'GBMTERT'}
  ordinaltype <- 'Primary'
  if (p %in% subsetRecurrent){ordinaltype <- 'Recurrent'}
  toBind <- data.frame(patient=p,proportion=proportion, subtype=subtype, ordinaltype=ordinaltype, stringsAsFactors=F)
  proportionsToPlot <- rbind(proportionsToPlot, toBind)
}
patientOrder <- proportionsToPlot[order(proportionsToPlot$proportion, decreasing=TRUE),]$patient
# make barplot order by proportions
toPlotFinal$patient <- factor(toPlotFinal$patient, levels=patientOrder)
toPlotFinal <- toPlotFinal[order(sapply(toPlotFinal$patientCluster, function(x) strsplit(x,'-')[[1]][2])),]
toPlotFinal <- toPlotFinal[order(toPlotFinal$patient),]
toPlotFinal$patientCluster <- factor(toPlotFinal$patientCluster, levels=toPlotFinal$patientCluster)
ggplot(data = toPlotFinal, aes(x=patientCluster, y=difference, fill=patient))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(size=10, angle=90, hjust=1, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
  scale_fill_manual(values=colors)
# then run stat tests after dropping P452
par(mfrow=c(1,4))
proportionsToPlotNoP452 <- proportionsToPlot[which(!proportionsToPlot$patient=='Patient452'),]
wilcox.test(proportionsToPlotNoP452[which(proportionsToPlotNoP452$subtype=='GBMTERT'),]$proportion,proportionsToPlotNoP452[which(proportionsToPlotNoP452$subtype=='Other'),]$proportion)
boxplot(proportionsToPlotNoP452$proportion~proportionsToPlotNoP452$subtype)
stripchart(proportion ~ subtype, vertical = TRUE, data = proportionsToPlotNoP452, 
           method = "jitter", add = TRUE, pch = 20, col = 'black')
wilcox.test(proportionsToPlot[which(proportionsToPlot$ordinaltype=='Primary'),]$proportion,proportionsToPlot[which(proportionsToPlot$ordinaltype=='Recurrent'),]$proportion)
boxplot(proportionsToPlot$proportion~proportionsToPlot$ordinaltype)
stripchart(proportion ~ ordinaltype, vertical = TRUE, data = proportionsToPlotNoP452, 
           method = "jitter", add = TRUE, pch = 20, col = 'black')
proportionsToPlot$TP53 <- ifelse(proportionsToPlot$patient %in% subsetTP53mut, 'TP53mut','TP53wt')
wilcox.test(proportionsToPlot[which(proportionsToPlot$TP53=='TP53mut'),]$proportion,proportionsToPlot[which(proportionsToPlot$TP53=='TP53wt'),]$proportion)
boxplot(proportionsToPlot$proportion~proportionsToPlot$TP53)
stripchart(proportion ~ TP53, vertical = TRUE, data = proportionsToPlotNoP452, 
           method = "jitter", add = TRUE, pch = 20, col = 'black')
proportionsToPlot_IDHmutOnly <- proportionsToPlot[proportionsToPlot$subtype=='Other',]
proportionsToPlot_IDHmutOnly$codel <- ifelse(proportionsToPlot_IDHmutOnly$patient %in% subset1p19q, 'codel','non-codel')
wilcox.test(proportionsToPlot_IDHmutOnly[which(proportionsToPlot_IDHmutOnly$codel=='codel'),]$proportion,proportionsToPlot_IDHmutOnly[which(proportionsToPlot_IDHmutOnly$codel=='non-codel'),]$proportion)
boxplot(proportionsToPlot_IDHmutOnly$proportion~proportionsToPlot_IDHmutOnly$codel)
stripchart(proportion ~ codel, vertical = TRUE, data = proportionsToPlot_IDHmutOnly, 
           method = "jitter", add = TRUE, pch = 20, col = 'black')

# add annotated column and output differences
toPlotFinal$cloneDistribution <- ifelse(toPlotFinal$difference < 0, 'lobular','intermixed')
write.table(toPlotFinal, file=paste0(outputPath, outfolder, 'clusters_all_patients_post_filtering_final_spatial_diff.txt'), quote=F, row.names=F, sep='\t')
