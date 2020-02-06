# Created: 2018.10.23
# By: Stephanie R Hilz
# Usage: For a specific patient, identifies shared mutations, groups them by clustering, 
  # and then does a spatial analysis on clusters

library(ggplot2)
library(kableExtra)
library(reshape2)
library(gplots)
library(gdata)
library(gtools)
library(grDevices)
library(vegan)

makeLongID <- function(gene, contig, position,refBase, altBase){
  longID <-paste(gene,contig,position,refBase,altBase, sep='_')
  return(longID)
}

getInfo <- function(avfile, samples){
  df <- read.delim(avfile, as.is=TRUE)
  if('X.gene' %in% colnames(df)){
    geneIndex <- grep('X.gene',colnames(df))
    colnames(df)[geneIndex] <- 'gene'
  }
  colnames(df) <- gsub('[.]','-',colnames(df))
  subset <- df[which(df$algorithm=='MuTect'),]
  n <- nrow(subset)
  toReturn <- c()
  for (name in samples){
    print(name)
    ref <- subset[,paste0(name,'_ref_Q20reads')]
    alt <- subset[,paste0(name,'_alt_Q20reads')]
    coverage <- ref + alt
    vaf <- alt/coverage
    sampleID <- rep(name,n)
    gene <- makeLongID(subset$gene,subset$contig,subset$position, subset$ref_allele, subset$alt_allele)
    chunk <- cbind(sampleID,gene,coverage,vaf)
    toReturn <- rbind(toReturn,chunk)
  }
  toReturn <- as.data.frame(toReturn, stringsAsFactors=F)
  toReturn$coverage <- as.numeric(toReturn$coverage)
  toReturn$vaf <- as.numeric(toReturn$vaf)
  toReturn$sampleID <- as.factor(toReturn$sampleID)
  colnames(toReturn) <- c('sampleID','uniqueID','coverage','vaf')
  return(toReturn)
}

getCensusGeneAnnotation <- function(geneSymbols,cancerCensusGenesFile){
  ## Develop sanger cancer gene census annotation for amplification sidebar
  head(geneSymbols)
  censusGenes <- read.csv(cancerCensusGenesFile)
  genes <- sapply(geneSymbols, getGene)
  head(genes)
  cancerGeneAnnotationDF <- as.data.frame(genes)
  cancerGeneAnnotationDF$Tier <- 'black'
  tier2Genes <- censusGenes[which(censusGenes$Tier=='2'),]$Gene.Symbol
  tier1Genes <- censusGenes[which(censusGenes$Tier=='1'),]$Gene.Symbol
  if (any(cancerGeneAnnotationDF$genes %in% tier2Genes)){
    cancerGeneAnnotationDF[cancerGeneAnnotationDF$genes %in% tier2Genes,]$Tier <- 'blue'
  }
  if (any(cancerGeneAnnotationDF$genes %in% tier1Genes)){
    cancerGeneAnnotationDF[cancerGeneAnnotationDF$genes %in% tier1Genes,]$Tier <- 'forestgreen'
  }
  return(cancerGeneAnnotationDF$Tier)
}

getGene<- function(x){
  strsplit(x,'_')[[1]][1]
}

plotMeanVafs <- function(meanVAFs, thresholds, VAFs, variants){
  par(mfrow=c(1,1), mar=c(14,5,2,4)+.1)
  plot(log(meanVAFs,10), las=2, col='white', ylim=c(-4,0), ylab="log10 VAF", cex.lab=1.5, cex.axis=1.5, xlab='', xaxt="n")
  lines(log(meanVAFs, 10), las=2, col='black', lty=1, lwd=2)
  for (j in variants){
    points(log(as.numeric(VAFs[j,]), 10), las=2, col='grey', pch=20, cex=1.5)
  }
  axis(1, at=seq(length(meanVAFs)),labels=names(meanVAFs), las=2, cex.axis=1.5)
  colors <- c('darkred','red3','red','deeppink')
  for (i in 1:length(thresholds)){
    abline(h=log(thresholds[i],10), col = colors[i])
    text(1,log(thresholds[i],10),labels = thresholds[i], adj = c(.2,-.2), col =colors[i], cex=1.5)
  }
  legend('bottomright', legend=c('cluster mean','individual'), bty='n', pch=c(NA,20), col=c('black','grey'), lty = c(1,0), lwd = c(2, 0))
}

calcWithinGroupDistance <- function(VAFs, k, clusters, distances, vafThreshold){
  #VAFs and distances should already be reshaped so that any replicates are averaged, and any missing samples are removed. Should both contain the same samples with the same names
  variants <- names(clusters[clusters==k])
  meanVAFs <- colMeans(VAFs[variants,])#will be the average VAF per sample for everything in that cluster
  within <- names(meanVAFs[which(meanVAFs > vafThreshold)])
  all <- names(meanVAFs)
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

## PART 1 - Patient-specific analysis - for one patient at a time (aggregate analysis for final fig at end)
# specify which patient you want to use and file paths
patientID <- 'Patient276'
cancerCensusGenesPath <- paste0(dataPath,'cancer_gene_census.csv')
outfolder <- 'SID000011_mutational_analysis_shared_group_expansion_analysis/'
purityCutoff <- .7

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

# specify high-CCF samples
patient_ID <- gsub('Patient','P',patientID)
highPuritySamples <- as.character(merged[which(merged$Patient == patient_ID),]$sample_type )

# specify mutfile to use
avfFile <- paste0(dataPath,patientID,'.R.mutations.avf.txt')

# run muts.bin function, which creates the same logical matrix used for phylo trees
muts.bin <- getMutsBin(avfFile) #muts.bin for all samples in the patient, not yet subsetted

# make the column names so they can string match our sample names
colnames(muts.bin) <- gsub('[.]','-',colnames(muts.bin))
colnames(muts.bin) <- gsub('_called','',colnames(muts.bin))

# ID shared mutations  for patient
muts.bin.subset <- muts.bin[,highPuritySamples]
sharedNames <- rownames(muts.bin.subset[which(rowSums(muts.bin.subset)< ncol(muts.bin.subset) & rowSums(muts.bin.subset)>1),])

# get mut info
info <- getInfo(avfFile, highPuritySamples)

# subset for shared mutations only
infoShared <- info[which(info$uniqueID %in% sharedNames),]

# create matrix of vafs
vafs <- acast(infoShared, uniqueID~sampleID, value.var="vaf")
colnames(vafs) <- gsub('[.]','-',colnames(vafs))

# plot
toPlot <- vafs
genes <- sapply(rownames(toPlot), getGene)
annotationColors <- getCensusGeneAnnotation(genes, cancerCensusGenesPath)
heatmap.2(toPlot, trace="none", sepcolor='black', colsep=1:ncol(toPlot), rowsep=1:nrow(toPlot),sepwidth=c(0.01,0.01), colRow=annotationColors, margins=c(13,13))

# final create a correlation matrix
correlation <- cor(t(toPlot), method="pearson")
dissimilarity <- (1 - correlation)/2
toPlot3 <- dissimilarity
genes <- sapply(rownames(toPlot3), getGene)
annotationColors <- getCensusGeneAnnotation(genes, cancerCensusGenesPath)
heatmap.2(toPlot3, trace="none", sepcolor='black', colsep=1:ncol(toPlot3), rowsep=1:nrow(toPlot3),sepwidth=c(0.01,0.01), colRow=annotationColors, col=bluered, margins=c(13,13))

# #Identify optimal number for k (both Elbow method and Calinski) - can combine this info with hclust to know where to best cut the tree
kmax <- 4
wss <- (nrow(vafs)-1)*sum(apply(vafs,2,var))
for (i in 2:kmax){
  wss[i] <- sum(kmeans(vafs,centers=i)$withinss)
}
plot(1:kmax, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
fit <- cascadeKM(scale(dissimilarity, center = TRUE,  scale = TRUE), 1, 30, iter = 1000)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)

# Cut the tree
if (patientID=="Patient303" | patientID=="Patient450" | patientID=='Patient260' | patientID=="Patient482"){
  optimal.k=5
}
if (patientID=="Patient327" | patientID=="Patient453" | patientID=="Patient454"){
  optimal.k=4
}
if (patientID=="Patient375" | patientID=="Patient413"){
  optimal.k=3
}
if (patientID=="Patient300"){
  optimal.k=6
}
if (patientID=="Patient452"){
  optimal.k=10
}
if (patientID=='Patient276'){
  optimal.k=2
}
clusters <- cutree(hclust(as.dist(dissimilarity)), k=optimal.k)

# are all on same chrom, or have too low a correlation (hardcoded in for each patient after manual inspection)
banned <- c()
clusterFrequency <- table(clusters)
banned <- names(clusterFrequency[clusterFrequency==1])
clusters <- clusters[!clusters %in% banned]
write.table(clusters[order(clusters)], file=paste0(outputPath, outfolder, patientID, '_clusters_post_filtering.txt'), quote=F, col.names=F, sep='\t')

# create a distance matrix
coordinates <- merged[which(merged$Patient==patient_ID),c('L.Coordinate','P.Coordinate','S.Coordinate')]
rownames(coordinates) <- merged[which(merged$Patient==patient_ID),'sample_type']
exomeDistanceMatrix <- as.matrix(dist(coordinates, method = "euclidean"))

# Subset distance matrix to correspond with vafs (since vafs come from a file already filtered for purity, should already have low-purity samples removed)
toUse <- mixedsort(colnames(vafs)[colnames(vafs) %in% colnames(exomeDistanceMatrix)])
vafs <- vafs[,toUse]
reshapedDistances <- exomeDistanceMatrix[toUse,toUse]

#Main
thresholds <- c(.1,.05,.01,.001)
outputRandom <- c()
outputSummary <- c()
for (i in 1:max(clusters)){
  for(j in 1:length(thresholds)){
    if (i %in% clusters){
      print(i)
      y <- calcWithinGroupDistance(vafs, i, clusters, reshapedDistances, thresholds[j])
      print(y)
      meanRandom <- mean(y$random)
      outputRandom <- rbind(outputRandom,cbind(i,thresholds[j],y$random))
      outputSummary <- rbind(outputSummary,c(i,thresholds[j],y$within,meanRandom,y$p,y$n, paste(y$samples,collapse=',')))
    }
  }
}
colnames(outputRandom)<-c('cluster','threshold','meanDistance')
colnames(outputSummary) <- c('cluster','threshold','meanDistanceWithin', 'meanDistanceRandom','p','n','samples')
outputRandom <- as.data.frame(outputRandom)
outputRandom$meanDistance <- as.numeric(outputRandom$meanDistance)
outputRandom$cluster <- factor(outputRandom$cluster, levels = sort(unique(clusters)))
outputRandom$threshold <- factor(outputRandom$threshold, levels = thresholds)
outputSummary <- as.data.frame(outputSummary)
outputSummary$meanDistanceWithin <- as.numeric(as.character(outputSummary$meanDistanceWithin))
outputSummary$meanDistanceRandom <- as.numeric(as.character(outputSummary$meanDistanceRandom))
outputSummary$cluster <- factor(outputSummary$cluster, levels = sort(unique(clusters)))
outputSummary$threshold <- factor(outputSummary$threshold, levels = thresholds)
thresholdColors <- c('black','gray40','gray63','gray85')
ggplot(data = outputRandom, aes(y=meanDistance, x=cluster, fill = threshold)) + 
  labs(y="Average pairwise distance (mm)") + 
  geom_boxplot(position="dodge") + 
  geom_point(data=outputSummary, aes(y=meanDistanceWithin, x=cluster, size=4), position=position_dodge(width=0.8), color='red', shape=17) +
  scale_fill_manual(values=thresholdColors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
saveRDS(outputSummary, file=paste0(outputPath, outfolder, patientID,'_outputSummary'))
write.table(outputSummary, file=paste0(outputPath, outfolder, patientID,'_outputSummary.txt'), quote=F,sep='\t', row.names=F)

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
patientsToUse <- c("Patient303","Patient327",'Patient375','Patient453','Patient450','Patient260','Patient452','Patient482',"Patient413",'Patient454','Patient276')
IDHwtTERTp <- c("Patient413",'Patient454','Patient276')
Recurrent <- c("Patient300","Patient260",'Patient450','Patient276')
toPlotFinal <- c()
for (patientID in patientsToUse){
  file=paste0(outputPath, outfolder, patientID,'_outputSummary')
  data <- readRDS(file)
  data$patient <- patientID
  toPlotFinal <- rbind(toPlotFinal,data)
}
# remove any where the n number of samples that contain an expansion event == total number of samples; also annotate patient colors while we loop
toExclude <- c()
toPlotFinal$color <- NA
for (p in unique(toPlotFinal$patient)){
  print(p)
  patient_ID <- gsub('Patient','P',p)
  totalSampleN <- nrow(merged[which(merged$Patient==patient_ID & !is.na(merged$WES_ID) & merged$SampleType=='SM'),])
  print(totalSampleN)
  toExclude <- append(toExclude, which(toPlotFinal$patient == p & toPlotFinal$n == totalSampleN))
  toPlotFinal[which(toPlotFinal$patient == p),]$color <- as.character(colorKey[patient_ID])
  print(as.character(colorKey[patient_ID]))
}
toPlotFinal <- toPlotFinal[-toExclude,]
toPlotFinal$patientCluster <- paste0(gsub("Patient",'P',toPlotFinal$patient),'-',toPlotFinal$cluster)
order <- unique(toPlotFinal$patientCluster)
toPlotFinal$patientCluster <- factor(toPlotFinal$patientCluster, levels=order)
toPlotFinal$type <- 'Other'
toPlotFinal[which(toPlotFinal$patient %in% PrimaryGBMs),]$type <- 'PrimaryGBM'
# plot by difference 
toPlotFinal$difference <- toPlotFinal$meanDistanceWithin-toPlotFinal$meanDistanceRandom
ggplot(data = toPlotFinal, aes(x=patientCluster, y=difference, group=threshold))+
  geom_point(aes(colour=threshold),size=10, shape='+')+
  theme(axis.text.x = element_text(size=10, angle=90, hjust=1, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
  geom_hline(yintercept=0,colour='black')+
  scale_colour_manual(values=thresholdColors)
toPlotFinalSingleThresh <- toPlotFinal[which(toPlotFinal$threshold == 0.05),]
toPlotFinalSingleThresh$patient <- factor(toPlotFinalSingleThresh$patient, levels=patientOrder)
toPlotFinalSingleThresh <- toPlotFinalSingleThresh[order(match(toPlotFinalSingleThresh$patient, patientOrder)),]
toPlotFinalSingleThresh$patientCluster <- factor(toPlotFinalSingleThresh$patientCluster, levels=toPlotFinalSingleThresh$patientCluster)
colors <- as.character(colorKey[gsub('Patient','P',patientOrder)])
ggplot(data = toPlotFinalSingleThresh, aes(x=patientCluster, y=difference, fill=patient))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(size=10, angle=90, hjust=1, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
  scale_fill_manual(values=colors)
# Look at proportion of clusters per patient that are < random
proportionsToPlot <- data.frame(patient=character(), proportion=numeric(), subtype=character(), ordinaltype=character(), stringsAsFactors=F)
for (p in unique(toPlotFinalSingleThresh$patient)){
  toPlotFinalSingleThreshPatient <- toPlotFinalSingleThresh[which(toPlotFinalSingleThresh$patient==p),]
  lessThanThreshold <- nrow(toPlotFinalSingleThreshPatient[which(toPlotFinalSingleThreshPatient$difference < 0),])
  totaln <- nrow(toPlotFinalSingleThreshPatient)
  proportion <- lessThanThreshold/totaln
  subtype <- 'Other'
  if (p %in% IDHwtTERTp){subtype <- 'GBMTERT'}
  ordinaltype <- 'Primary'
  if (p %in% Recurrent){ordinaltype <- 'Recurrent'}
  toBind <- data.frame(patient=p,proportion=proportion, subtype=subtype, ordinaltype=ordinaltype, stringsAsFactors=F)
  proportionsToPlot <- rbind(proportionsToPlot, toBind)
}
patientOrder <- proportionsToPlot[order(proportionsToPlot$proportion, decreasing=TRUE),]$patient
proportionsToPlotNoP452 <- proportionsToPlot[which(!proportionsToPlot$patient=='Patient452'),]
ks.test(proportionsToPlotNoP452[which(proportionsToPlotNoP452$subtype=='GBMTERT'),]$proportion,proportionsToPlotNoP452[which(proportionsToPlotNoP452$subtype=='Other'),]$proportion)
boxplot(proportionsToPlotNoP452$proportion~proportionsToPlotNoP452$subtype)
ks.test(proportionsToPlot[which(proportionsToPlot$ordinaltype=='Primary'),]$proportion,proportionsToPlot[which(proportionsToPlot$ordinaltype=='Recurrent'),]$proportion)
boxplot(proportionsToPlot$proportion~proportionsToPlot$ordinaltype)


