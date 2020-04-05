## Costello Lab
## 2019.09.13
## Stephanie Hilz
## Usage: Takes tumor metadata + other data files, creates a final matrix for clustering. Also has a function for comparing the
##  difference in any metric in the matrix and the difference in spatial distance between the samples

computeZscore <- function(cpmList){
  zscores <- c()
  mean <- mean(cpmList, na.rm=T)
  sd <- sd(cpmList, na.rm=T)
  for (i in cpmList){
    z <- (i - mean)/sd
    zscores <- append(zscores, z)
  }
  names(zscores) <- names(cpmList)
  return(zscores)
}

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

dist.pear <- function(x) as.dist(1-cor(t(x)))

pairwiseDist2pointsLPS <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + (((p2[3]-p1[3]))^2) ))
}

computeMeanZscore <- function(cpmMatrix){
  zscoreMatrix <- cpmMatrix
  for (i in 1:nrow(cpmMatrix)){
    unlisted <- cpmMatrix[i,] %>% unlist
    zscoreMatrix[i,] <- computeZscore(unlisted)
  }
  zScoreMean <- colMeans(zscoreMatrix)
  return(zScoreMean)
}

modelPurity <- function(dataMatrix, toCorrelate, outstatsFileName, outheatmapFileName){
  # for purity, fits a linear model to see its relationship to each cell type/state signature
  #  that controls for repeat measures 
  dataMatrix$Patient <- as.factor(dataMatrix$Patient)
  allPairsPurity <- cbind(rbind(c("purity"),toCorrelate))
  allPairsPurityResults <- data.frame(a=character(), b=character(),estimate=numeric(), pvalue=numeric())
  for (i in 1:ncol(allPairsPurity)){
    print(i)
    a <- allPairsPurity[1,i]
    b <- allPairsPurity[2,i]
    if (!(a == b)){
      # pdf(paste0(outPath, 'correlation_plot',a,'_values_',b,'_within_patient_lm.pdf'))
      # print(ggplot(dataMatrix, aes_string(x=a, y=b, colour="Patient")) +
      #         geom_point()+
      #         geom_smooth(method = "lm", se=F, size=.5)+
      #         stat_poly_eq(aes(label = paste(..rr.label..)),
      #                      label.x.npc = "right", label.y.npc = 1, formula = y~x, 
      #                      parse = TRUE, size = 3) + 
      #         stat_fit_glance(method = 'lm',
      #                         geom = 'text',
      #                         aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
      #                         label.x.npc = 'right', label.y.npc = 0.4, size = 3))
      # dev.off()
      fit <- lmer( dataMatrix[,b] ~ dataMatrix[,a] + (1|dataMatrix[,'Patient']))
      print(summary(fit))
      estimate <- summary(fit)$coefficients['dataMatrix[, a]','Estimate']
      pvalue <- summary(fit)$coefficients["dataMatrix[, a]","Pr(>|t|)"]
      toBind <- data.frame(a,b,estimate,pvalue, stringsAsFactors = F)
      allPairsPurityResults <- rbind.data.frame(allPairsPurityResults, toBind)
    }
  }
  allPairsPurityResults$corr.pvalue <- p.adjust(allPairsPurityResults$pvalue, method='BH')
  allPairsPurityResults$isSignificant <- 0
  allPairsPurityResults[which(allPairsPurityResults$corr.pvalue < 0.05),]$isSignificant <- 1
  write.table(allPairsPurityResults, file=paste0(outFolder,outstatsFileName), sep='\t',quote=F, row.names=F)
  ofInterest <- 'isSignificant'
  colOrder <- allPairsPurityResults[order(allPairsPurityResults$isSignificant),]$b
  mofInterest <- acast(allPairsPurityResults, a~b, value.var=ofInterest)
  mofInterest <- mofInterest[,colOrder]
  mofInterest <- rbind(mofInterest,mofInterest)
  if (ncol(mofInterest) > 1){
    pdf(paste0(outFolder, outheatmapFileName), width=10, height=10)
    print(heatmap.2(mofInterest[order(rowSums(is.na(mofInterest))),order(colSums(is.na(mofInterest)))],
              Rowv=F,
              Colv=F,
              trace='none',
              margins = c(7,7)))
    graphics.off()
    my_palette <- colorRampPalette(c("mediumblue", "grey", "red"))(n = 299)
    ofInterest <- 'estimate'
    allPairsPurityResultsSig <- allPairsPurityResults[which(allPairsPurityResults$isSignificant==1),]
    colOrder <- allPairsPurityResultsSig[order(allPairsPurityResultsSig$estimate),]$b
    mofInterest <- acast(allPairsPurityResultsSig, a~b, value.var=ofInterest)
    mofInterest <- mofInterest[,colOrder]
    mofInterest <- rbind(mofInterest,mofInterest)
    outheatmapEstFileName <- gsub('heatmap','heatmapEstimate',outheatmapFileName)
    break1 <- -.0000001
    break2 <- .0000001
    absolute <- max(abs(mofInterest), na.rm=T)
    pdf(paste0(outFolder, outheatmapEstFileName), width=10, height=10)
    my_palette <- colorRampPalette(c("mediumblue", "grey", "orangered4"))(n = 299)
    col_breaks = c(seq(-absolute,break1,length=100), # for yellow (~1st quartile)
                   seq(break1+.001,break2-.001,length=100), # for black (~within 2nd and 3rd quartiles)
                   seq(break2,absolute,length=100)) # for blue (~4th quartile)
    hm.parameters <- list(mofInterest,cellwidth=10, cellheight=10, fontsize=10, border_col=NA, cluster_rows = F, cluster_cols = F, col= my_palette, breaks = col_breaks, na_col='white')
    do.call("pheatmap", hm.parameters)
    dev.off()
  }
  return(allPairsPurityResults)
}

modelIDHInteractPurity <- function(dataMatrix, toCorrelate, outstatsFileName, outheatmapFileName){ #tests for influence of IDH mut status 
  # for purity, fits a linear model to see its relationship to each cell type/state signature
  #  that controls for repeat measures 
  dataMatrix$Patient <- as.factor(dataMatrix$Patient)
  allPairsPurity <- cbind(rbind(c("purity"),toCorrelate))
  allPairsPurityResults <- data.frame(a=character(), b=character(),estimate=numeric(), pvalue=numeric())
  for (i in 1:ncol(allPairsPurity)){
    print(i)
    a <- allPairsPurity[1,i]
    b <- allPairsPurity[2,i]
    if (!(a == b)){
      fit <- lmer( dataMatrix[,b] ~ dataMatrix[,a]*dataMatrix[,'IDH_Mut'] + (1|dataMatrix[,'Patient']))
      print(summary(fit))
      estimate <- summary(fit)$coefficients['dataMatrix[, a]','Estimate']
      pvalue <- summary(fit)$coefficients["dataMatrix[, a]","Pr(>|t|)"]
      pvalueIDH <- summary(fit)$coefficients['dataMatrix[, "IDH_Mut"]',"Pr(>|t|)"]
      toBind <- data.frame(a,b,estimate,pvalue, pvalueIDH, stringsAsFactors = F)
      allPairsPurityResults <- rbind.data.frame(allPairsPurityResults, toBind)
    }
  }
  allPairsPurityResults$corr.pvalue <- p.adjust(allPairsPurityResults$pvalue, method='BH')
  allPairsPurityResults$corr.pvalueIDH <- p.adjust(allPairsPurityResults$pvalueIDH, method='BH')
  allPairsPurityResults$isSignificant <- 0
  if (any(allPairsPurityResults$corr.pvalueIDH < 0.05)){
    allPairsPurityResults[which(allPairsPurityResults$corr.pvalueIDH < 0.05),]$isSignificant <- 1
  }
  write.table(allPairsPurityResults, file=paste0(outFolder,outstatsFileName), sep='\t',quote=F, row.names=F)
  return(allPairsPurityResults)
}

plotFinal <- function(ResultsIDHMut, ResultsIDHWT,  outstatsFileName){
  merged <- merge(ResultsIDHMut,ResultsIDHWT, by='b')
  mergedSubset <- merged[,c('b','estimate.x','isSignificant.x','estimate.y','isSignificant.y')]
  colnames(mergedSubset) <- c('b','estimate.IDHMut','isSignificant.IDHMut','estimate.IDHWT','isSignificant.IDHWT')
  mofInterest <- mergedSubset[,c('estimate.IDHWT','estimate.IDHMut')]#[which(mergedSubset$isSignificant.IDHMut ==1 | mergedSubset$isSignificant.IDHWT ==1),c('estimate.IDHWT','estimate.IDHMut')]
  rownames(mofInterest) <- mergedSubset$b#[which(mergedSubset$isSignificant.IDHMut ==1 | mergedSubset$isSignificant.IDHWT ==1),]$b
  my_palette <- colorRampPalette(c("mediumblue", "grey", "red"))(n = 299)
  colOrder <- order(rowMeans(mofInterest))
  rowOrder <- c('estimate.IDHMut','estimate.IDHWT')
  mofInterest <- t(mofInterest)[rowOrder,colOrder]
  break1 <- -.001
  break2 <- .001
  absolute <- max(abs(mofInterest), na.rm=T)
  pdf(paste0(outFolder, outstatsFileName), width=10, height=10)
  my_palette <- colorRampPalette(c("blue", "grey", "red"))(n = 299)
  col_breaks = c(seq(-absolute,break1,length=100), # for yellow (~1st quartile)
                 seq(break1-(break1/1000),break2-(break2/1000),length=100), # for black (~within 2nd and 3rd quartiles)
                 seq(break2,absolute,length=100)) # for blue (~4th quartile)
  hm.parameters <- list(mofInterest,cellwidth=10, cellheight=10, fontsize=10, border_col=NA, cluster_rows = F, cluster_cols = F, col= my_palette, breaks = col_breaks, na_col='white')
  do.call("pheatmap", hm.parameters)
  graphics.off()
}

modelDistances <- function(dataMatrix, toCorrelate, allPairsPurityResults, outstatsFileName, outheatmapFileName, purityCorrect){
  # for each distance metric, fit a linear model to see its relationship to each cell type/state signature
  #  that controls for repeat measures (allows for different intercept for each patient)
  allPairs <- cbind(rbind(c("DistCentroid"),toCorrelate),rbind(c("DistPeriph"),toCorrelate),rbind(c("DistVR"),toCorrelate))#combn(toCorrelateAll, 2)
  allPairsResults <- data.frame(a=character(), b=character(),purityCorrection=numeric(), estimate=numeric(), pvalue=numeric())
  for (i in 1:ncol(allPairs)){
    print(i)
    a <- allPairs[1,i]
    b <- allPairs[2,i]
    if (!(a == b)){
      print(b)
      # pdf(paste0(outPath, 'correlation_plot',a,'_values_',b,'_within_patient_lm.pdf'))
      # print(ggplot(mtmp, aes_string(x=a, y=b, colour="Patient")) +
      #         geom_point()+
      #         geom_smooth(method = "lm", se=F, size=.5)+
      #         stat_poly_eq(aes(label = paste(..rr.label..)),
      #                      label.x.npc = "right", label.y.npc = 1, formula = y~x, 
      #                      parse = TRUE, size = 3) + 
      #         stat_fit_glance(method = 'lm',
      #                         geom = 'text',
      #                         aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
      #                         label.x.npc = 'right', label.y.npc = 0.4, size = 3))
      # dev.off()
      purityCorrection <- 0
      if (purityCorrect==0){
        fit <- lmer( dataMatrix[,b] ~ dataMatrix[,a] + (1|dataMatrix[,'Patient']))
      } else {
        if (allPairsPurityResults[which(allPairsPurityResults$b == b),]$isSignificant == 1){ # still only correct if purity even has a role in this
          fit <- lmer( dataMatrix[,b] ~ dataMatrix[,a] + dataMatrix[,'purity'] + (1|dataMatrix[,'Patient']))
          purityCorrection <- 1
        } else {
        fit <- lmer( dataMatrix[,b] ~ dataMatrix[,a] + (1|dataMatrix[,'Patient']))
        }
      }
      print(summary(fit))
      estimate <- summary(fit)$coefficients['dataMatrix[, a]','Estimate']
      pvalue <- summary(fit)$coefficients["dataMatrix[, a]","Pr(>|t|)"]
      toBind <- data.frame(a,b,purityCorrection,estimate,pvalue, stringsAsFactors = F)
      allPairsResults <- rbind.data.frame(allPairsResults, toBind)
    }
  }
  allPairsResults$corr.pvalue <- NA
  for (metric in dataTypes['relTumorDistance'][[1]]){
    print(metric)
    allPairsResults[which(allPairsResults$a == metric),]$corr.pvalue <- p.adjust(allPairsResults[which(allPairsResults$a == metric),]$pvalue, method='BH')
  }
  allPairsResults$isSignificant <- 0
  if (any(allPairsResults$corr.pvalue < 0.05) & nrow(allPairsResults) > 3){
    allPairsResults[which(allPairsResults$corr.pvalue < 0.05),]$isSignificant <- 1
    write.table(allPairsResults, file=paste0(outFolder,outstatsFileName), sep='\t',quote=F, row.names=F)
    ofInterest <- 'isSignificant'
    DistVROnly <- allPairsResults[which(allPairsResults$a=='DistVR'),]
    colOrder <- DistVROnly[order(DistVROnly$isSignificant),]$b
    mofInterest <- acast(allPairsResults, a~b, value.var=ofInterest)
    mofInterest <- mofInterest[,colOrder]
    if (nrow(mofInterest) > 1){
      pdf(paste0(outFolder, outheatmapFileName), width=10, height=10)
      print(heatmap.2(mofInterest[order(rowSums(is.na(mofInterest))),order(colSums(is.na(mofInterest)))],
                    Rowv=F,
                    Colv=F,
                    trace='none',
                    margins = c(7,7)))
      graphics.off()
    }
    my_palette <- colorRampPalette(c("mediumblue", "grey", "red"))(n = 299)
    ofInterest <- 'estimate'
    allPairsResultsSig <- allPairsResults[which(allPairsResults$isSignificant==1 & !(allPairsResults$b %in% c('DistPeriph','DistCentroid','DistVR'))),]
    colOrder <- allPairsResultsSig[order(allPairsResultsSig$estimate),]$b
    mofInterest <- acast(allPairsResultsSig, a~b, value.var=ofInterest)
    mofInterest <- mofInterest[,colOrder]
    mofInterest <- rbind(mofInterest,mofInterest)
    outheatmapEstFileName <- gsub('heatmap','heatmapEstimate',outheatmapFileName)
    break1 <- -.0001
    break2 <- .0001
    absolute <- max(abs(mofInterest), na.rm=T)
    pdf(paste0(outFolder, outheatmapEstFileName), width=10, height=10)
    my_palette <- colorRampPalette(c("mediumblue", "grey", "orangered4"))(n = 299)
    col_breaks = c(seq(-absolute,break1,length=100), # for yellow (~1st quartile)
                   seq(break1-(break1/1000),break2-(break2/1000),length=100), # for black (~within 2nd and 3rd quartiles)
                   seq(break2,absolute,length=100)) # for blue (~4th quartile)
    hm.parameters <- list(mofInterest,cellwidth=10, cellheight=10, fontsize=10, border_col=NA, cluster_rows = F, cluster_cols = F, col= my_palette, breaks = col_breaks, na_col='white')
    do.call("pheatmap", hm.parameters)
    graphics.off()
  } else{
    write.table(allPairsResults, file=paste0(outFolder,outstatsFileName), sep='\t',quote=F, row.names=F)
  }
  return(allPairsResults)
}

modelIDHInteractDistance <- function(dataMatrix, toCorrelate, allPairsPurityResults, outstatsFileName, outheatmapFileName, purityCorrect){ #tests for influence of IDH mut status 
  # for purity, fits a linear model to see its relationship to each cell type/state signature
  #  that controls for repeat measures 
  allPairs <- cbind(rbind(c("DistCentroid"),toCorrelate),rbind(c("DistPeriph"),toCorrelate),rbind(c("DistVR"),toCorrelate))#combn(toCorrelateAll, 2)
  allPairsResults <- data.frame(a=character(), b=character(),purityCorrection=numeric(), estimate=numeric(), pvalue=numeric())
  for (i in 1:ncol(allPairs)){
    print(i)
    a <- allPairs[1,i]
    b <- allPairs[2,i]
    if (!(a == b)){
      print(b)
      purityCorrection <- 0
      if (purityCorrect==0){
        fit <- lmer( dataMatrix[,b] ~ dataMatrix[,a]*dataMatrix[,'IDH_Mut'] + (1|dataMatrix[,'Patient']))
      } else {
        if (allPairsPurityResults[which(allPairsPurityResults$b == b),]$isSignificant == 1){ # still only correct if purity even has a role in this
          fit <- lmer( dataMatrix[,b] ~ dataMatrix[,a]*dataMatrix[,'IDH_Mut'] + dataMatrix[,'purity'] + (1|dataMatrix[,'Patient']))
          purityCorrection <- 1
        } else {
          fit <- lmer( dataMatrix[,b] ~ dataMatrix[,a]*dataMatrix[,'IDH_Mut'] + (1|dataMatrix[,'Patient']))
        }
      }
      print(summary(fit))
      estimate <- summary(fit)$coefficients['dataMatrix[, a]','Estimate']
      pvalue <- summary(fit)$coefficients["dataMatrix[, a]","Pr(>|t|)"]
      pvalueIDH <- summary(fit)$coefficients['dataMatrix[, "IDH_Mut"]',"Pr(>|t|)"]
      toBind <- data.frame(a,b,purityCorrection,estimate,pvalue, pvalueIDH, stringsAsFactors = F)
      allPairsResults <- rbind.data.frame(allPairsResults, toBind)
    }
  }
  allPairsResults$corr.pvalue <- NA
  allPairsResults$corr.pvalueIDH <- NA
  for (metric in dataTypes['relTumorDistance'][[1]]){
    print(metric)
    allPairsResults[which(allPairsResults$a == metric),]$corr.pvalue <- p.adjust(allPairsResults[which(allPairsResults$a == metric),]$pvalue, method='BH')
    allPairsResults[which(allPairsResults$a == metric),]$corr.pvalueIDH <- p.adjust(allPairsResults[which(allPairsResults$a == metric),]$pvalueIDH, method='BH')
  }
  allPairsResults$isSignificant <- 0
  if (any(allPairsResults$corr.pvalueIDH < 0.05)){
    allPairsResults[which(allPairsResults$corr.pvalueIDH < 0.05),]$isSignificant <- 1
  }
  write.table(allPairsResults, file=paste0(outFolder,outstatsFileName), sep='\t',quote=F, row.names=F)
  return(allPairsResults)
}

plotCorrelationByPatient <- function(dataMatrix, a, b){
  pdf(paste0(outFolder, 'correlation_plot',a,'_values_',b,'_within_patient_lm.pdf'))
  print(ggplot(dataMatrix, aes_string(x=a, y=b, colour="Patient")) +
          geom_point()+
          geom_smooth(method = "lm", se=F, size=.5)+
          theme(axis.text.x = element_text(size=10, angle=90, color='black'), axis.title = element_text(size = 10, color='black'), axis.text.y = element_text(size=10, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
          stat_poly_eq(aes(label = paste(..rr.label..)),
                       label.x.npc = "right", label.y.npc = 1, formula = y~x,
                       parse = TRUE, size = 3) +
          stat_fit_glance(method = 'lm',
                          geom = 'text',
                          aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                          label.x.npc = 'right', label.y.npc = 0.4, size = 3))
  dev.off()
}

library(gtools)
library(gplots)
library(ggpmisc)
library(kableExtra)
library(dplyr)
library(reshape2)
library(ggpubr)
library(lmerTest)
library(pheatmap)

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# User-defined variables and file paths (will also need to make sure VenteicherSignatures folder is in dataPath)
tag <- 'SID000018'
outFolder <- paste0(outputPath, 'SID000018_overall_tumor_landscape/')
cancerSEAResultsFile <- paste0(outputPath,'SID000017_generate_cancerSEA_dataset/cancerSEA_signature_zscores.txt') #cell states
oldhamDeconResultsFile <- paste0(dataPath,'filtered_module_eigenes_Bicor-None_signum0.514_minSize8_minMEcor0.85_16000.txt') #brain cell types

# Create an empty list to store type data
dataTypes <- list()

## Get sample metadata
# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)
# remove samples without exome data
data <- data[which(!is.na(data$WES_ID)),]
# remove samples that do not have distance from periphery/centroid metrics
data <- data[which(!is.na(data$DistPeriph)),]
# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)
# merge by patient ID
merged <- merge(data, subtypedata, by="Patient")
# add in column to link patient info to RNAseq sampleID
rownames(merged) <- paste0(merged$Patient,merged$SampleName)
# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF
merged[which(is.na(merged$FACETS)),]$purity <- .1 #here we choose to set anything where FACETS gave a purity of NA to 0
merged$purity <- 100 * merged$purity
# add data types to list
dataTypes[['relTumorDistance']] <- c('DistCentroid','DistPeriph','DistVR')
dataTypes[['Purity']] <- c('purity')

## read in cancerSEA file
dataCancerSEA <- read.table(cancerSEAResultsFile, sep='\t', header = T, row.names=1,stringsAsFactors = F)
dataTypes[['CancerSEA']] <- colnames(dataCancerSEA)

## read in Oldham lab data
dataOldhamDecon <- read.table(oldhamDeconResultsFile, sep='\t', header = T, row.names=1,stringsAsFactors = F)
conversion <- as.data.frame(rbind(c('lightcyan','lymphocytes'),
                                  c('brown4','m1macro'),
                                  c('brown','microgliamacro'),
                                  c('palevioletred2','granulocyte'),
                                  c('yellowgreen','endothelial'),
                                  c('darkseagreen3','astrocytes'),
                                  c('lavenderblush3','astrocytes2'),
                                  c('black','oligodendrocytes'),
                                  c('saddlebrown','ependymal'),
                                  c('pink','neuron')), stringsAsFactors = F)
colnames(conversion) <- c('color','cellType')
colnames(dataOldhamDecon) <- conversion[match(colnames(dataOldhamDecon),conversion$color),]$cellType
dataTypes[['OldhamDecon']] <- colnames(dataOldhamDecon)

# merge these two datasets together
mtmp <- transform(merge(dataCancerSEA,dataOldhamDecon,by=0), row.names=Row.names, Row.names=NULL)

# add in metadata
mtmp <- transform(merge(mtmp, merged, by=0), row.names=Row.names, Row.names=NULL)

# apply conversions (transform purity and distance from periphery/centroid into z scores)
mtmp$zDistCentroid <- computeZscore(mtmp$DistCentroid) 
mtmp$zWithinDistCentroid <- NA #same as above but within a patient
for (p in unique(mtmp$Patient)){ 
  mtmp[which(mtmp$Patient == p),]$zWithinDistCentroid <- computeZscore(mtmp[which(mtmp$Patient == p),]$DistCentroid) 
}
mtmp$zDistPeriph <- computeZscore(mtmp$DistPeriph) 
mtmp$zWithinDistPeriph <- NA #same as above but within a patient
for (p in unique(mtmp$Patient)){ 
  mtmp[which(mtmp$Patient == p),]$zWithinDistPeriph <- computeZscore(mtmp[which(mtmp$Patient == p),]$DistPeriph) 
}
mtmp$zDistVR <- computeZscore(mtmp$DistVR) 
mtmp$zWithinDistVR <- NA #same as above but within a patient
for (p in unique(mtmp$Patient)){ 
  mtmp[which(mtmp$Patient == p),]$zWithinDistVR <- computeZscore(mtmp[which(mtmp$Patient == p),]$DistVR) 
}
mtmp$zPurity <- computeZscore(mtmp$purity) 
dataTypes[['zrelTumorDistance']] <- c('zDistCentroid','zDistPeriph','zDistVR')
dataTypes[['zWithinrelTumorDistance']] <- c('zWithinDistCentroid','zWithinDistPeriph','zWithinDistVR')
dataTypes[['zPurity']] <- 'zPurity'
mtmp$Patient <- as.factor(mtmp$Patient)

# for purity, fits a linear model to see its relationship to each cell type/state signature and corrects group for mult testing
#  that controls for repeat measures 
# running note - for some reason need to manually specify dataMatrix <- mtmp, or won't run.

##### ##### #####
# 1) PURITY
###
# a) total SNVS
toTest <- as.vector(c('totalSNVsCalled'))
# purity with total SNVs
testName <- '_purity_vs_totalSNVs_lmer_with_re_padjBH_'
purityTotalSNVsResults <- modelPurity(mtmp, toTest, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'))

# repeats above purity for idh-mut tumors only
testName <- '_IDHMutonly_purity_vs_totalSNVs_lmer_with_re_padjBH_'
purityTotalSNVsResultsIDHMut <- modelPurity(mtmp[which(mtmp$IDH_Mut==1),], toTest, paste0(tag,testName, 'stats.txt'),paste0(tag,testName,'heatmap.pdf'))

# repeats above purity for idh-wt TERT tumors only
testName <- '_IDHWTTERTonly_purity_vs_totalSNVs_lmer_with_re_padjBH_'
purityTotalSNVsResultsIDHWT <- modelPurity(mtmp[which(mtmp$IDH_Mut==0 & mtmp$TERT==1),], toTest, paste0(tag, testName, 'stats.txt'),paste0(tag,testName, 'heatmap.pdf'))

# finally tests for idh impact
testName <- '_IDHTest_purity_vs_totalSNVs_lmer_with_re_padjBH_'
purityTotalSNVsResultsIDHTest <- modelIDHInteractPurity(mtmp, toTest, paste0(tag,testName, 'stats.txt'),paste0(tag,testName, 'heatmap.pdf'))

###
# b) Cancer signatures
toTest <- as.vector(dataTypes['CancerSEA'] %>% unlist)
# purity with cancer signatures
testName <- '_purity_vs_cancersignatures_lmer_with_re_padjBH_'
purityCancerSignaturesResults <- modelPurity(mtmp, toTest, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'))

# repeats above purity for idh-mut tumors only
testName <- '_IDHMutonly_purity_vs_cancersignatures_lmer_with_re_padjBH_'
purityCancerSignaturesResultsIDHMut <- modelPurity(mtmp[which(mtmp$IDH_Mut==1),], toTest, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'))

# repeats above purity for idh-wt TERT tumors only
testName <- '_IDHWTTERTonly_purity_vs_cancersignatures_lmer_with_re_padjBH_'
purityCancerSignaturesResultsIDHWT <- modelPurity(mtmp[which(mtmp$IDH_Mut==0 & mtmp$TERT==1),], toTest, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'))

# finally tests for idh impact
testName <- '_IDHTest_purity_vs_cancersignatures_lmer_with_re_padjBH_'
purityCancerSignaturesResultsIDHTest <- modelIDHInteractPurity(mtmp, toTest, paste0(tag,testName, 'stats.txt'),paste0(tag,testName, 'heatmap.pdf'))

# final purity plot
testName <- '_summary_purity_vs_cancersignatures_lmer_with_re_padjBH_'
plotFinal(purityCancerSignaturesResultsIDHMut, purityCancerSignaturesResultsIDHWT, paste0(tag,testName,'_composite_heatmap.pdf'))

###
# c) Brain cell types (Oldham)
toTest <- as.vector(dataTypes['OldhamDecon'] %>% unlist)
# purity with cell types
testName <- '_purity_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
purityOldhamCellTypesResults <- modelPurity(mtmp, toTest, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'))

# repeats above purity for idh-mut tumors only
testName <- '_IDHMutonly_purity_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
purityOldhamCellTypesResultsIDHMut <- modelPurity(mtmp[which(mtmp$IDH_Mut==1),], toTest, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'))

# repeats above purity for idh-wt TERT tumors only
testName <- '_IDHWTTERTonly_purity_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
purityOldhamCellTypesResultsIDHWT <- modelPurity(mtmp[which(mtmp$IDH_Mut==0 & mtmp$TERT==1),], toTest, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'))

# finally tests for idh impact
testName <- '_IDHTest_purity_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
purityOldhamCellTypesResultsIDHTest <- modelIDHInteractPurity(mtmp, toTest, paste0(tag,testName, 'stats.txt'),paste0(tag,testName, 'heatmap.pdf'))

# final purity plot
testName <- '_summary_purity_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
plotFinal(purityOldhamCellTypesResultsIDHMut, purityOldhamCellTypesResultsIDHWT, paste0(tag,testName,'_composite_heatmap.pdf'))

# 2) DISTANCE
# for each distance metric, fit a linear model to see its relationship to each cell type/state signature
#  that controls for repeat measures (allows for different intercept for each patient)
###
# a) purity
toTest <- as.vector(c('purity'))
# distance with total SNVs
testName <- '_distance_vs_purity_lmer_with_re_padjBH_'
distancePurityResults <- modelDistances(mtmp, toTest, NA, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'),0)

# repeats above purity for idh-mut tumors only
testName <- '_IDHMutonly_distance_vs_purity_lmer_with_re_padjBH_'
distancePurityResultsIDHMut <- modelDistances(mtmp[which(mtmp$IDH_Mut==1),], toTest, NA, paste0(tag,testName, 'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 0)

# repeats above purity for idh-wt TERT tumors only
testName <- '_IDHWTTERTonly_distance_vs_purity_lmer_with_re_padjBH_'
distancePurityResultsIDHWT <- modelDistances(mtmp[which(mtmp$IDH_Mut==0 & mtmp$TERT==1),], toTest, NA, paste0(tag, testName, 'stats.txt'),paste0(tag,testName, 'heatmap.pdf'), 0)

# finally tests for idh impact
testName <- '_IDHTest_distance_vs_purity_lmer_with_re_padjBH_'
distancePurityResultsIDHTest <- modelIDHInteractDistance(mtmp, toTest, NA, paste0(tag,testName, 'stats.txt'),paste0(tag,testName, 'heatmap.pdf'), 0)

###
# b) Cancer signatures
toTest <- as.vector(dataTypes['CancerSEA'] %>% unlist)
# distance with cancer signatures
testName <- '_distance_vs_cancersignatures_lmer_with_re_padjBH_'
distanceCancerSignaturesResults <- modelDistances(mtmp, toTest, purityCancerSignaturesResults, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 1) #we do correct for purity with these because is a propery of cancer cells we are interested in

# repeats above purity for idh-mut tumors only
testName <- '_IDHMutonly_distance_vs_cancersignatures_lmer_with_re_padjBH_'
distanceCancerSignaturesResultsIDHMut <- modelDistances(mtmp[which(mtmp$IDH_Mut==1),], toTest, purityCancerSignaturesResultsIDHMut, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 1) #we do correct for purity with these because is a propery of cancer cells we are interested in

# repeats above purity for idh-wt TERT tumors only
testName <- '_IDHWTTERTonly_distance_vs_cancersignatures_lmer_with_re_padjBH_'
distanceCancerSignaturesResultsIDHWT <- modelDistances(mtmp[which(mtmp$IDH_Mut==0 & mtmp$TERT==1),], toTest, purityCancerSignaturesResultsIDHWT, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 1) #we do correct for purity with these because is a propery of cancer cells we are interested in

# finally tests for idh impact
testName <- '_IDHTest_distance_vs_cancersignatures_lmer_with_re_padjBH_'
distanceCancerSignaturesResultsIDHTest <- modelIDHInteractDistance(mtmp, toTest, purityCancerSignaturesResults, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 1) #we do correct for purity with these because is a propery of cancer cells we are interested in

###
# c) Brain cell types (Oldham)
toTest <- as.vector(dataTypes['OldhamDecon'] %>% unlist)
# distance with cancer signatures
testName <- '_distance_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
distanceOldhamCellTypesResults <- modelDistances(mtmp, toTest, NA, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 0) 

# repeats above purity for idh-mut tumors only
testName <- '_IDHMutonly_distance_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
distanceOldhamCellTypesResultsIDHMut <- modelDistances(mtmp[which(mtmp$IDH_Mut==1),], toTest, NA, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 0) 

# repeats above purity for idh-wt TERT tumors only
testName <- '_IDHWTTERTonly_distance_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
distanceOldhamCellTypesResultsIDHWT <- modelDistances(mtmp[which(mtmp$IDH_Mut==0 & mtmp$TERT==1),], toTest, NA, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 0) 

# finally tests for idh impact
testName <- '_IDHTest_distance_vs_oldhambraincelltypes_lmer_with_re_padjBH_'
distanceOldhamCellTypesResultsIDHTest <- modelIDHInteractDistance(mtmp, toTest, NA, paste0(tag,testName,'stats.txt'),paste0(tag,testName,'heatmap.pdf'), 0) 
