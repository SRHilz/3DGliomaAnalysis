## Costello Lab
## 2018.10.15
## Stephanie Hilz
## Usage: Takes results of SID000015_runSciClone.R and plots them

library(ggplot2)
library(kableExtra)
library(reshape2)
library(gplots)
library(gdata)
library(gtools)
library(grDevices)

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

calcWithinGroupDistanceSciClone <- function(sciClone, k, distances, vafThreshold){
  #sciCone object of mean cluster vafs and distances should already be reshaped so that any replicates are averaged, and any missing samples are removed. Should both contain the same samples with the same names
  meanVAFs <- sciClone[k,]#will be the average VAF per sample for everything in that cluster
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



# user-defined cutoffs
purityCutoff <- 0.7

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# data file paths
sciCloneClusterFile <- paste0(outputPath,'SID000015_runSciClone/P452.cluster.summary.tsv.means.txt')
sciCloneGenesClustersFile <- paste0(outputPath,'SID000015_runSciClone/P452.cluster.tsv')
outfolder <- 'SID000016_P452_sciclone_analysis/'
cancerCensusGenesPath <- paste0(dataPath,'cancer_gene_census.csv')

# read in cluster mean data
sciClone <- read.table(sciCloneClusterFile,sep='\t', header=T, stringsAsFactors = F, row.names=1)
colnames(sciClone) <- gsub('Primary.','',colnames(sciClone))

my_palette <- colorRampPalette(c("blue","red"))(n = 299)

heatmap.2(as.matrix(sciClone),
          dendrogram = 'column',
          Rowv=FALSE,
          trace='none',
          col=my_palette
          )

# reformat data for plotting
meltedSciClone <- sciClone
meltedSciClone$cluster <- rownames(meltedSciClone)
meltedSciClone <- melt(meltedSciClone, id.vars="cluster")
colnames(meltedSciClone) <- c('cluster','sample','meanVAF')
meltedSciClone$cluster <- as.factor(meltedSciClone$cluster)

ggplot(meltedSciClone, aes(x = sample, y = meanVAF, fill = cluster)) +
  geom_bar(stat="identity", position=position_dodge()) 

par(mfrow=c(1,1), mar=rep(4,4))
colors <- rainbow(nrow(sciClone))
for (c in 1:nrow(sciClone)){
  col = colors[c]
  toPlot <- sciClone[c,] %>% unlist
  if (c==1){
    plot(toPlot, col='white', xaxt='n', ylim=c(0,1), ylab = 'mean cluster VAF', xlab='')
  }
  lines(toPlot, col = col)
}
axis(1,at=seq(colnames(sciClone)), labels=colnames(sciClone), cex=.9, las=2)
legend('topleft',legend=rownames(sciClone),col=colors, lty=1, bty='n', cex=.8)

## Pull out genes in clusters and annotate with Sanger Info
sciCloneGenesClusters <- read.table(sciCloneGenesClustersFile,sep='\t', header=T, stringsAsFactors = F)
sciCloneGenesClusters$sanger <- getCensusGeneAnnotation(sciCloneGenesClusters$gene)
sciCloneGenesClusters$uniq <- paste0(sciCloneGenesClusters$gene, '_chr', sciCloneGenesClusters$chr, '_',sciCloneGenesClusters$st)

patientID <- 'Patient452'
infile <- paste0(dataPath,patientID,'.R.mutations.avf.txt')
muts <- read.delim(infile, as.is=TRUE)
if('X.gene' %in% colnames(muts)){
  geneIndex <- grep('X.gene',colnames(muts))
  colnames(muts)[geneIndex] <- 'gene'
}
muts$uniq <- paste(muts$gene, muts$contig, muts$position, sep="_")
sciCloneGenesClustersAnno <- merge(sciCloneGenesClusters, muts, by='uniq')
sciCloneGenesClustersAnno[which(sciCloneGenesClustersAnno$SangerCancerGeneCensus.=='YES' & (!sciCloneGenesClustersAnno$type=='Silent')),c('uniq','type','SangerCancerGeneCensus.','cluster'),]
countsNonSilentSanger <- table(sciCloneGenesClustersAnno[which(sciCloneGenesClustersAnno$SangerCancerGeneCensus.=='YES' & (!sciCloneGenesClustersAnno$type=='Silent')),c('uniq','type','SangerCancerGeneCensus.','cluster'),]$cluster)
countsNonSilentAll <- table(sciCloneGenesClustersAnno[which((!sciCloneGenesClustersAnno$type=='Silent')),c('uniq','type','SangerCancerGeneCensus.','cluster'),]$cluster)
countsNonSilentDifference <- countsNonSilentAll - countsNonSilentSanger
toPlot <- as.data.frame(rbind(cbind(countsNonSilentSanger,names(countsNonSilentAll),'sanger'),cbind(countsNonSilentDifference,names(countsNonSilentDifference),'non-sanger')), stringsAsFactors=F)
colnames(toPlot) <- c('counts','cluster','type')
toPlot$counts <- as.numeric(toPlot$counts)
toPlot$cluster <- as.factor(toPlot$cluster)
toPlot$type <- as.factor(toPlot$type)
ggplot(toPlot, aes(x=cluster, y=counts, fill=type))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(size=20, colour='black'), axis.title = element_text(size = 20), axis.text.y = element_text(size=20, colour='black'), panel.background = element_rect(fill = 'white', colour = 'black'))

## look at distance among each clone and see if it is closer than chance - same as SID000011 save do not need in this case to predict clusters, as are already predicted

## PART 1 - Patient-specific analysis - for one patient at a time (aggregate analysis for final fig at end)
# specify which patient you want to use and file paths

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

# subset clusters by these samples
colnames(sciClone) <- gsub('v','Primary-v',colnames(sciClone))
sciCloneToCompare <- sciClone[,highPuritySamples]

# create a distance matrix (pairwise distance in mm between samples)
coordinates <- merged[which(merged$Patient==patient_ID),c('L.Coordinate','P.Coordinate','S.Coordinate')]
rownames(coordinates) <- merged[which(merged$Patient==patient_ID),'sample_type']
exomeDistanceMatrix <- as.matrix(dist(coordinates, method = "euclidean"))

# Subset distance matrix to correspond with vafs (since vafs come from a file already filtered for purity, should already have low-purity samples removed)
toUse <- mixedsort(colnames(sciCloneToCompare)[colnames(sciCloneToCompare) %in% colnames(exomeDistanceMatrix)])
sciCloneToCompare <- sciCloneToCompare[,toUse]
reshapedDistances <- exomeDistanceMatrix[toUse,toUse]

#Main
thresholds <- c(.1,.05,.01,.001)
clusters <- rownames(sciCloneToCompare)
outputRandom <- c()
outputSummary <- c()
for (i in clusters){
  for(j in 1:length(thresholds)){
    y <- calcWithinGroupDistanceSciClone(sciCloneToCompare, i, reshapedDistances, thresholds[j])
    meanRandom <- mean(y$random)
    outputRandom <- rbind(outputRandom,cbind(i,thresholds[j],y$random))
    outputSummary <- rbind(outputSummary,c(i,thresholds[j],y$within,meanRandom,y$p,y$n, paste(y$samples,collapse=',')))
  }
}
colnames(outputRandom)<-c('cluster','threshold','meanDistance')
colnames(outputSummary) <- c('cluster','threshold','meanDistanceWithin', 'meanDistanceRandom','p','n','samples')
outputRandom <- as.data.frame(outputRandom)
outputRandom$meanDistance <- as.numeric(as.character(outputRandom$meanDistance))
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
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10, color='black'), axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))
saveRDS(outputSummary, file=paste0(outputPath, outfolder, patientID,'_outputSummary'))
write.table(outputSummary, file=paste0(outputPath, outfolder, patientID,'_outputSummary.txt'), quote=F,sep='\t', row.names=F)
# simple plot with single threshold
toPlotFinalSingleThresh <- outputSummary[which(outputSummary$threshold == 0.05),]
toPlotFinalSingleThresh <- toPlotFinalSingleThresh[which(!(toPlotFinalSingleThresh$n == 1 | toPlotFinalSingleThresh$n == length(highPuritySamples))),]
toPlotFinalSingleThresh$difference <- toPlotFinalSingleThresh$meanDistanceWithin - toPlotFinalSingleThresh$meanDistanceRandom
ggplot(data = toPlotFinalSingleThresh, aes(x=cluster, y=difference))+
  ylim(-7,5)+
  geom_bar(stat="identity", fill='grey')+
  theme(axis.text.x = element_text(size=10, angle=90, hjust=1, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))



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



