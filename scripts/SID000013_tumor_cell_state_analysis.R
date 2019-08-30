## Costello Lab
## 2019.09.13
## Stephanie Hilz
## Usage: Takes in CPMs, study config files, and cell state gene enrichment lists and performs an array of analyses 
##  to look at how cell states are distributed in glioma

computeZscore <- function(cpmList){
  zscores <- c()
  mean <- mean(cpmList)
  sd <- sd(cpmList)
  for (i in cpmList){
    z <- (i - mean)/sd
    zscores <- append(zscores, z)
  }
  names(zscores) <- names(cpmList)
  return(zscores)
}

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

library(gtools)
library(gplots)
library(ggpmisc)
library(kableExtra)

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# User-defined variables and file paths (will also need to make sure VenteicherSignatures folder is in dataPath)
tag <- 'SID000013'
outfolder <- 'SID000013_tumor_cell_states/'
CPMFile <- paste0(dataPath,'SID000003_20190529_first_submission.symbol.coding.CPMs.csv')

## read in files
CPMs <- read.csv(CPMFile)
rownames(CPMs) <- CPMs$X
CPMs$X <- NULL
sampleN <- ncol(CPMs)
sampleNames <- colnames(CPMs)
my_palette <- colorRampPalette(c("black", "pink", "red"))(n = 1000)
patientsToUse <- c('P41','P303','P327','P375','P453','P413','P454','P300','P302','P450') #will also use this order; did not include P340 and 455 as only had 1 sample each over purity; P457 had 0
purityCutoff <- 0.7

coreOligo<- scan(file=paste0(dataPath,'VenteicherSignatures/Venteicher_Oligo_yellow.txt'),what=character())
allOligo<- scan(file=paste0(dataPath,'VenteicherSignatures/Venteicher_Oligo_all.txt'),what=character())
coreAstro<- scan(file=paste0(dataPath,'VenteicherSignatures/Venteicher_Astro_green.txt'),what=character())
allAstro<- scan(file=paste0(dataPath,'VenteicherSignatures/Venteicher_Astro_all.txt'),what=character())
stemness <- scan(file=paste0(dataPath,'VenteicherSignatures/Venteicher_stemness_all.txt'),what=character())
coreOligo <- coreOligo[coreOligo %in% rownames(CPMs)]
allOligo <- allOligo[allOligo %in% rownames(CPMs)]
coreAstro <- coreAstro[coreAstro %in% rownames(CPMs)]
allAstro <- allAstro[allAstro %in% rownames(CPMs)]
stemness <- stemness[stemness %in% rownames(CPMs)]
signatures <- list(coreOligo=coreOligo, allOligo=allOligo, coreAstro=coreAstro, allAstro=allAstro, stemness=stemness)
signatureNames <- names(signatures)

## Get sample metadata
# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)
# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)
# merge by patient ID
merged <- merge(data, subtypedata, by="Patient")
# add in column to link patient info to RNAseq sampleID
merged$sampleID <- paste0(merged$Patient,merged$SampleName)
# subset object to just include information relevant to our RNAseq samples
merged <- merged[which(merged$sampleID %in% sampleNames),]
merged$sampleID <- paste0(merged$Patient,merged$SampleName)

# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF

# Filter based on purity and inspect numbers of samples per patient
merged <- merged[which(merged$purity >= purityCutoff),]
kable(table(merged$Patient)) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold=TRUE, background = "#C2BFBA")


# set color
colorKey <- subtypedata$Color
names(colorKey) <- as.character(subtypedata$Patient)

# set up structures to collect data in
zscores <- setNames(replicate(length(signatureNames),data.frame()),signatureNames)
geneExpression <- setNames(replicate(length(signatureNames),data.frame()),signatureNames)

# populate data structures
for (sig in signatureNames){
  signatureGenes <- signatures[sig] %>% unlist
  signatureGeneExpression <- CPMs[signatureGenes,]
  # populate our object that collects CPMs for each gene in signature
  geneExpression[[sig]] <- rbind(geneExpression[[sig]], signatureGeneExpression)
  # populate our object that collcts scores for each signature
  for (patient in patientsToUse){
    print(patient)
    patientSamples <- merged[which(merged$Patient==patient),]$sampleID
    patientSignatureGeneExpression <- signatureGeneExpression[,patientSamples]
    zscore <- computeMeanZscore(patientSignatureGeneExpression)
    means <- colMeans(patientSignatureGeneExpression)
    scoreOut <- cbind(patient,patientSamples,zscore, means)
    rownames(scoreOut) <- NULL
    zscores[[sig]] <- rbind(zscores[[sig]], scoreOut)
  }
  colnames(zscores[[sig]]) <- c('patient','sampleID','zscore', 'means')
  zscores[[sig]]$patient <- as.character(zscores[[sig]]$patient)
  zscores[[sig]]$sampleID <- as.character(zscores[[sig]]$sampleID)
  zscores[[sig]]$zscore <- as.numeric(as.character(zscores[[sig]]$zscore))
  zscores[[sig]]$means <- as.numeric(as.character(zscores[[sig]]$means))
}

# Exploration of gene expression
sig <- 'stemness'
m <- as.matrix(geneExpression[[sig]])
dist.pear <- function(x) as.dist(1-cor(t(x)))
pdf(paste0(sig,'_heatmap.pdf'), width=30, height=7)
h <- heatmap.2(m, 
               trace="none", 
               dendrogram="both",
               distfun=dist.pear,
               density.info="none",
               cexCol=1,
               scale='row',
               col = my_palette)
dev.off()

# Exploration of z scores
sig <- 'coreAstro'
toPlot <- zscores[[sig]][order(zscores[[sig]]$zscore),]
barplot(toPlot$zscore)
boxplot(toPlot$zscore~toPlot$patient, ylab=paste0(sig, " signature Zscore"))

## add in distance info to samples and their zscores
mergedSubset <- merged[,c('SampleName','L.Coordinate','P.Coordinate','S.Coordinate','DistCentroid','DistPeriph','Histology','Tumor','purity','sampleID')]

# Create a data structure that also includes distances (right now just have for Astros)
zscoresWithDistances <- setNames(replicate(length(signatureNames),data.frame()),signatureNames)
for (s in signatureNames){
  zscoresWithDistances[[s]] <- merge(zscores[[s]],mergedSubset, by="sampleID")
  zscoresWithDistances[[s]]$peripheralCat <- 0
  zscoresWithDistances[[s]][which(zscoresWithDistances[[s]]$DistPeriph <= 7),]$peripheralCat <- 1
  zscoresWithDistances[[s]]$peripheralCat <- factor(zscoresWithDistances[[s]]$peripheralCat, levels=c(0,1))
  zscoresWithDistances[[s]]$centroidCat <- 0
  zscoresWithDistances[[s]][which(zscoresWithDistances[[s]]$DistCentroid <= 17),]$centroidCat <- 1
  zscoresWithDistances[[s]]$centroidCat <- factor(zscoresWithDistances[[s]]$centroidCat, levels=c(0,1))
  zscoresWithDistances[[s]]$patient <- factor(zscoresWithDistances[[s]]$patient,levels=patientsToUse)
}

## Create zscoresWithDistance that does not contain missing Centroid/Periph data
zscoresWithDistancesCentPeriph <- zscoresWithDistances[[sig]][!is.na(zscoresWithDistances[[sig]]$DistCentroid),]

# plot by correlation (change x to go between peripheral and core, and sig to go between signatures)
colors <- as.character(colorKey[patientsToUse[patientsToUse %in% unique(zscoresWithDistancesCentPeriph$patient)]])
ggplot(zscoresWithDistancesCentPeriph, aes(x=DistCentroid, y=zscore, color=patient)) +
  geom_point() +
  scale_colour_manual(values=colors) +
  theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_smooth(aes(colour=factor(patient)), method = "lm", se=F)+
  labs(y=paste0(sig, ' z-score')) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x.npc = "right", label.y.npc = 0.8, formula = y~x, 
               parse = TRUE, size = 3) + 
  stat_fit_glance(method = 'lm',
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.4, size = 3)

## Plot comparison of score per sample by patient by peripheral or centroid
ggplot(zscoresWithDistancesCentPeriph, aes(x = patient, y = zscore, fill=centroidCat)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Dark2") +
  labs(y=paste0(sig, ' z-score')) +
  theme(axis.text.x = element_text(size=10, angle = 90, hjust = 1), axis.title = element_text(size = 10), axis.text.y = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'black'))

## Calculate spearmans
patients <- patientsToUse[patientsToUse %in% unique(zscoresWithDistancesCentPeriph$patient)]
colors <- colors
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 20 #where on the x axis will display p value
y=max(toReturnSignature$difference)/1.5
x.var = "DistPeriph"
y.var = "zscore"
y.inc = y/15
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- zscoresWithDistancesCentPeriph[which(zscoresWithDistancesCentPeriph$patient == patientID),]
  testResult <- cor.test(patientSubset[,x.var], patientSubset[,y.var], method="spearman")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  rho=round(testResult$estimate,3)
  label <- paste0('p=',p,', rho=',rho)
  dataText <- rbind(dataText, c(p,rho,label,x,y,color,patientID), stringsAsFactors=F)
  y <- y + y.inc
}
colnames(dataText) <- c('p','rho','label','x','y','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'_',x.var,'_vs_',sig,'_zscore_difference_stats.txt'),sep='\t', quote=F, row.names = F)

## Finally, look at relationship between pairwise distances and score differences
# create a data structure that also contains distance to do some correlations with
toReturn <- data.frame(signature=character(), patient=character(), samplePair=character(), difference=numeric(), distance=numeric(), stringsAsFactors=F)
for (p in unique(zscoresWithDistances[[sig]]$patient)){
  zscoresWithDistancesLocalP <- zscoresWithDistances[[sig]][which(zscoresWithDistances[[sig]]$patient == p),]
  print(p)
  comparisons <- combn(zscoresWithDistancesLocalP$sampleID,2)
  p <- as.character(p)
  print(dim(comparisons))
  # Loop through and calculate euclidian spatial distance as well as signature difference for each pair
  for (c in seq_len(dim(comparisons)[2])){
    signature <- sig
    sample_a <- comparisons[1,c]
    sample_b <- comparisons[2,c]
    samplePair <- paste0(sample_a, '-', sample_b)
    difference <- abs(zscoresWithDistancesLocalP[which(zscoresWithDistancesLocalP$sampleID == sample_a),]$zscore - zscoresWithDistancesLocalP[which(zscoresWithDistancesLocalP$sampleID == sample_b),]$zscore)
    Sample1LPS <- zscoresWithDistancesLocalP[which(zscoresWithDistancesLocalP$sampleID == sample_a),c("L.Coordinate","P.Coordinate","S.Coordinate")]
    Sample2LPS <- zscoresWithDistancesLocalP[which(zscoresWithDistancesLocalP$sampleID == sample_b),c("L.Coordinate","P.Coordinate","S.Coordinate")]
    distance <- as.numeric(pairwiseDist2pointsLPS(Sample1LPS,Sample2LPS))
    toReturn <- rbind(toReturn, data.frame(signature=signature, patient=p, samplePair=samplePair, difference=difference, distance=distance, stringsAsFactors=F))
  }
}

colors <- as.character(colorKey[patientsToUse[patientsToUse %in% unique(toReturn$patient)]])
ggplot(toReturn, aes(x=distance, y=difference, color=patient)) +
  geom_point() +
  scale_colour_manual(values=colors) +
  theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_smooth(aes(colour=factor(patient)), method = "lm", se=F)+
  labs(y=paste0(sig, ' z-score pairwise difference')) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x.npc = "right", label.y.npc = 1.0, formula = y~x, 
               parse = TRUE, size = 3) + 
  stat_fit_glance(method = 'lm',
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.4, size = 3)


## Calculate spearmans
df <- toReturn
patients <- patientsToUse[patientsToUse %in% unique(df$patient)]
dataText <- data.frame(p=numeric(), R=numeric(), label=character(), x=numeric(), y=numeric(),m=numeric(), b=numeric(), color=character(), patient=character(), stringsAsFactors=F)
x = 20 #where on the x axis will display p value
x.var = "distance"
y.var = "difference"
y=max(df[,y.var])/1.5
y.inc = y/15
for (i in rev(seq_along(patients))){
  patientID=as.character(patients[i])
  print(patientID)
  color <- colors[i]
  patientSubset <- df[which(df$patient == patientID),]
  testResult <- cor.test(patientSubset[,x.var], patientSubset[,y.var], method="spearman")
  p=formatC(testResult$p.value,format = "e", digits = 2)
  rho=round(testResult$estimate,3)
  label <- paste0('p=',p,', rho=',rho)
  dataText <- rbind(dataText, c(p,rho,label,x,y,color,patientID), stringsAsFactors=F)
  y <- y + y.inc
}
colnames(dataText) <- c('p','rho','label','x','y','color','patient')
dataText$x <- as.numeric(dataText$x)
dataText$y <- as.numeric(dataText$y)
write.table(dataText, file=paste0(outputPath,outfolder,tag,'_spatial_distance_vs_',sig,'_zscore_difference_stats.txt'),sep='\t', quote=F, row.names = F)
