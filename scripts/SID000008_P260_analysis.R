# Created: 2018.05.28
# By: Stephanie R Hilz
# Usage: Perform all final analyses and generate plots for P260

library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(tidyr)

# user-specified variables
purityCutoff <- 0.7
lateral_samples <- paste0('v',seq(1,10,1))
medial_samples <- paste0('v',seq(11,20,1))
PDGFRA_lateral_id <- 'PDGFRA_chr4_55133815_C_A'
PDGFRA_medial_id <- 'PDGFRA_chr4_55138638_G_T'
IDH_id <- 'IDH1_chr2_209113112_C_T'

# define file paths
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')
cnFile <- paste0(dataPath, 'Patient260.GeneCopyNumber.rds')
cpmFile <- paste0(dataPath, 'SID000003_20190913_expanded_gbm.symbol.coding.CPMs.csv')
ssgseaInput <-  paste0(dataPath,'Patient260_ssGSEA_results.txt')
ciberSortFile <- paste0(dataPath, '20190529_first_submission.symbol.coding.CPMs.CIBERSORT.Output_Job8.txt')
mutfile <- paste0(dataPath, 'Patient260.R.mutations.avf.txt')

  
# read in mutfile and get relevant vafs
muts <- read.delim(mutfile, as.is=TRUE)
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
colnames(alt_counts) <- gsub('Recurrence1-','',colnames(alt_counts))
ref_counts <- muts[,grepl('_ref_Q20reads', colnames(muts))]
colnames(ref_counts) <- gsub('_ref_Q20reads','',colnames(ref_counts))
colnames(ref_counts) <- gsub('Recurrence1-','',colnames(ref_counts))
vafs <- alt_counts/(alt_counts+ref_counts)

# define sample categories and genes of interest

# read in different data types
cpms <- read.csv(cpmFile, row.names=1)
cn <- readRDS(cnFile)
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# define purity and sample categories of interest
data$purity <- data$FACETS
dataPurity <- data[which(data$Patient=="P260"),c('SampleName','purity')]
lateral_high_purity <- dataPurity[which(dataPurity$purity >= purityCutoff & dataPurity$SampleName %in% lateral_samples),]$SampleName
medial_high_purity <- dataPurity[which(dataPurity$purity >= purityCutoff & dataPurity$SampleName %in% medial_samples),]$SampleName
medail_low_purity <- dataPurity[which(dataPurity$purity < purityCutoff & dataPurity$SampleName %in% medial_samples),]$SampleName
samplesToPlot <- c(lateral_high_purity,medial_high_purity)

# take log of cpms
logcpms <- log(cpms+.001, 2)

# make PDGFRA plot
par(mfrow=c(2,1))
#first plot of variants (mutation data)
IDH1 <- vafs[IDH_id,samplesToPlot] %>% unlist
PDGFRA_medial <- vafs[PDGFRA_medial_id,samplesToPlot] %>% unlist
PDGFRA_lateral <- vafs[PDGFRA_lateral_id,samplesToPlot] %>% unlist
plot(IDH1, col='white', xaxt='n', ylim=c(0,1))
lines(IDH1, col = 'black')
lines(PDGFRA_lateral, col = 'blue')
lines(PDGFRA_medial, col = 'red')
axis(1,at=seq(colnames(vafs[,samplesToPlot])), labels=colnames(vafs[,samplesToPlot]), cex=.9, las=2)
#second plot of cpms
PDGFRAcpms <- cpms['PDGFRA',paste0('P260',samplesToPlot)] %>% unlist
plot(PDGFRAcpms, col='white', xaxt='n')
lines(PDGFRAcpms, col = 'black')
axis(1,at=seq(colnames(vafs[,samplesToPlot])), labels=colnames(vafs[,samplesToPlot]), cex=.9, las=2)

# pdgfra boxplot
par(mfrow=c(1,1))
toPlot <- data.frame(PDGFRAcpms[paste0('P260',lateral_high_purity)],'lateral')
colnames(toPlot) <- c('PDGFRA_cpms','aspect')
toBind <- data.frame(PDGFRAcpms[paste0('P260',medial_high_purity)],'medial')
colnames(toBind) <- c('PDGFRA_cpms','aspect')
toPlot <- rbind(toPlot, toBind)
ggplot(data = toPlot, aes(y=PDGFRA_cpms, x=aspect)) + 
  geom_boxplot(position="dodge", outlier.shape=NA) + 
  geom_jitter(shape=16, size=3, position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
  stat_compare_means(aes(group = PDGFRA_cpms), label = "p.format", method='wilcox.test')
wilcox.test(PDGFRAcpms[paste0('P260',lateral_high_purity)],PDGFRAcpms[paste0('P260',medial_high_purity)])

# make PDGFA and C plot
par(mfrow=c(1,1))
#first plot of variants (mutation data)
#second plot of cpms
PDGFAcpms <- cpms['PDGFA',paste0('P260',samplesToPlot)] %>% unlist
PDGFCcpms <- cpms['PDGFC',paste0('P260',samplesToPlot)] %>% unlist
plot(PDGFAcpms, col='white', xaxt='n')
lines(PDGFAcpms, col = 'black', lty=3)
lines(PDGFCcpms, col = 'black', lty=5)
axis(1,at=seq(colnames(vafs[,samplesToPlot])), labels=colnames(vafs[,samplesToPlot]), cex=.9, las=2)
legend('topleft',c('PDGFA','PDGFC'), col=c('black','black'),lty=c(3,5))

# purity boxplot
par(mfrow=c(1,1))
toPlot <- data.frame(dataPurity[which(dataPurity$SampleName %in% lateral_samples),]$purity,'lateral')
colnames(toPlot) <- c('purity','aspect')
toBind <- data.frame(dataPurity[which(dataPurity$SampleName %in% medial_samples),]$purity,'medial')
colnames(toBind) <- c('purity','aspect')
toPlot <- rbind(toPlot, toBind)
ggplot(data = toPlot, aes(y=purity, x=aspect)) + 
  geom_boxplot(position="dodge", outlier.shape=NA) + 
  geom_jitter(shape=16, size=3, position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
wilcox.test(toPlot$purity~toPlot$aspect)


# create PTEN cpm and cn line plot
dev.off()
par(mfrow=c(2,1))
PTENcn <- cn['PTEN',c(paste0('Recurrence1-',lateral_high_purity),paste0('Recurrence1-',medial_high_purity))] %>% unlist
names(PTENcn) <- c(lateral_high_purity, medial_high_purity)
plot(PTENcn, col='white', xaxt='n')
lines(PTENcn, col = 'black')
axis(1,at=seq(names(PTENcn)), labels=names(PTENcn), cex=.9, las=2)
PTENcpm <- cpms['PTEN',paste0('P260',samplesToPlot)] %>% unlist
names(PTENcpm) <- samplesToPlot
plot(PTENcpm, col='white', xaxt='n')
lines(PTENcpm, col = 'black')
axis(1,at=seq(names(PTENcpm)), labels=names(PTENcpm), cex=.9, las=2)

# create PTEN CN boxplot
par(mfrow=c(1,1))
# purity boxplot
par(mfrow=c(1,1))
toPlot <- data.frame(PTENcn[lateral_high_purity],'lateral')
colnames(toPlot) <- c('PTENcn','aspect')
toBind <- data.frame(PTENcn[medial_high_purity],'medial')
colnames(toBind) <- c('PTENcn','aspect')
toPlot <- rbind(toPlot, toBind)
ggplot(data = toPlot, aes(y=PTENcn, x=aspect)) + 
  geom_boxplot(position="dodge", outlier.shape=NA) + 
  geom_jitter(shape=16, size=3, position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
wilcox.test(toPlot$PTENcn~toPlot$aspect)


# create PTEN CPM boxplot
par(mfrow=c(1,1))
# purity boxplot
par(mfrow=c(1,1))
toPlot <- data.frame(PTENcpm[lateral_high_purity],'lateral')
colnames(toPlot) <- c('PTENcpm','aspect')
toBind <- data.frame(PTENcpm[medial_high_purity],'medial')
colnames(toBind) <- c('PTENcpm','aspect')
toPlot <- rbind(toPlot, toBind)
ggplot(data = toPlot, aes(y=PTENcpm, x=aspect)) + 
  geom_boxplot(position="dodge", outlier.shape=NA) + 
  geom_jitter(shape=16, size=3, position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))
wilcox.test(toPlot$PTENcpm~toPlot$aspect)

# create ssgsea pathway plot for MTOR and ERK
ssgsea <- read.table(ssgseaInput, row.names=1, header=T, sep='\t')
rownames(ssgsea) <- gsub('X','v', rownames(ssgsea))
ofInterest <- c('REACTOME_MTORC1_MEDIATED_SIGNALLING',
                'REACTOME_PROLONGED_ERK_ACTIVATION_EVENTS')
ssgsea <- ssgsea[samplesToPlot,]
toPlot <- ssgsea[,colnames(ssgsea) %in% ofInterest]
colnames(toPlot)[which(colnames(toPlot)=='REACTOME_MTORC1_MEDIATED_SIGNALLING')] <- 'MTOR Sig.'
colnames(toPlot)[which(colnames(toPlot)=='REACTOME_PROLONGED_ERK_ACTIVATION_EVENTS')] <- 'ERK Actv.'
#heatmap.2(t(toPlot), trace="none", Colv=F, col=brewer.pal(9,"YlGn"), cexRow=.5)
toPlotGG <- toPlot
toPlotGG$sample <- rownames(toPlotGG)
toPlotGG$aspect <- 'lateral'
toPlotGG[which(toPlotGG$sample %in% medial_high_purity),]$aspect <- 'medial'
toPlotGG <- melt(toPlotGG)
colnames(toPlotGG) <- c('sample','aspect','pathway','score')
ggplot(toPlotGG, aes(x=pathway, y=score, fill=aspect)) +
  geom_boxplot(position=position_dodge(0.8), outlier.shape=NA) +
  geom_point(position=position_jitterdodge()) +
  scale_fill_manual(values=c("blue","red"))+
  labs(list(y = "Enrichment score") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))
testResults <- data.frame(pathway=character(), pvalue=numeric())
for (p in unique(toPlotGG$pathway)){
  print(paste0('test for lateral vs medial for pathway: ' , p))
  toPlotGGPathway <- toPlotGG[which(toPlotGG$pathway == p),]
  scoreLateral <- toPlotGGPathway[which(toPlotGGPathway$aspect == 'lateral'),]$score
  print(scoreLateral)
  scoreMedial <- toPlotGGPathway[which(toPlotGGPathway$aspect == 'medial'),]$score
  print(scoreMedial)
  if (p == 'MTOR Sig.'){
    result <- ks.test(scoreMedial, scoreLateral, alternative = 'less')
    print(result)
    pvalue <- result$p.value
  } else {
    result <- ks.test(scoreMedial, scoreLateral, alternative = 'greater')
    print(result)
    pvalue <- result$p.value
  }
  testResults <- rbind(testResults, data.frame(p, pvalue))
}

# create immune cell plot for immune cell enrichment differences (uses all sampels for which we have RNAseq for)
## Read in Cibersort results  
cibersortResults <- read.table(ciberSortFile, sep='\t', header=T)
## Add in patient tags
cibersortResults$Patient <- substr(cibersortResults$Input.Sample, 1, 4)
cibersortResults[cibersortResults$Patient %in% c("P41V","P41Y","P41P","P41G"),]$Patient <- "P41"
## Pull out P260
cibersortResults260 <- cibersortResults[which(cibersortResults$Patient == 'P260'),]
## Fix sample names
cibersortResults260$Input.Sample <- gsub('P260','',cibersortResults260$Input.Sample)
rownames(cibersortResults260) <- cibersortResults260$Input.Sample
## Remove low-purity samples
cibersortResults260 <- cibersortResults260[samplesToPlot,]
## Pull out interesting columns
ofInterest <- c('Input.Sample','B.cells.naive','T.cells.CD8','T.cells.CD4.memory.resting','T.cells.regulatory..Tregs.','Macrophages.M0')
cibersortResults260 <- cibersortResults260[,ofInterest]
## Reshape and plot
toPlot <- melt(cibersortResults260, id.vars='Input.Sample')
colnames(toPlot) <- c('sample','cellType','score')
toPlot$aspect <- 'lateral'
toPlot[which(toPlot$sample %in% medial_samples),]$aspect <- 'medial'
ggplot(toPlot, aes(x=cellType, y=score, fill=aspect)) +
  geom_boxplot(position=position_dodge(0.5), width=0.5, outlier.shape=NA) +
  geom_point(position=position_jitterdodge()) +
  scale_fill_manual(values=c("blue","red"))+
  labs(list(y = "Enrichment score") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))
testResults <- data.frame(pathway=character(), pvalue=numeric())
for (p in unique(toPlot$cellType)){
  print(paste0('test for lateral vs medial for pathway: ' , p))
  toPlotCellType <- toPlot[which(toPlot$cellType == p),]
  scoreLateral <- toPlotCellType[which(toPlotCellType$aspect == 'lateral'),]$score
  print(scoreLateral)
  scoreMedial <- toPlotCellType[which(toPlotCellType$aspect == 'medial'),]$score
  print(scoreMedial)
  result <- ks.test(scoreMedial, scoreLateral)
  pvalue <- result$p.value
  testResults <- rbind(testResults, data.frame(p, pvalue))
}

testResults$adj.pvalue <- p.adjust(testResults$pvalue, method='BH')

## Pathology analysis - conducted on samples adjacent to those we sequenced; used pathologist-determined % tumor nuclei for "purity" for these 
file <- '/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Projects/2018_Patient260_Case_Study/Data/Pathology_CollectionSpreadsheet.txt'
data <- read.table(file, sep='\t', header = T, stringsAsFactors = F)
data$Color <- paste0('v',data$Color)
samplesWithTumor <- data[which(data$`Target...tumor.nuclei` > 0),]$Color
data$color <- 'blue'
data[which(data$Color %in% medial_samples),]$color <- 'red'
data$aspect <- as.factor(data$aspect)
# BV Hyperplasia chart - show in all samples, since is a feature of microenvironment
dotchart(data$Target.BV.hyperplasia, data$Color, col=data$color, ylab="BV Hyperplasia")
toChi <- rbind(table(data[data$Target.BV.hyperplasia > 0,]$aspect), table(data[data$Target.BV.hyperplasia <= 0,]$aspect))
rownames(toChi) <- c('<=0','>0')
chisq.test(toChi) # not sig different
# Ki67 boxplot - only in those samples with tumor cells
toPlot <- data.frame(data[which(data$Color %in% lateral_samples & data$Color %in% samplesWithTumor),]$Target.Ki67.BTRC.Count,'lateral')
colnames(toPlot) <- c('Ki67','aspect')
toBind <- data.frame(data[which(data$Color %in% medial_samples & data$Color %in% samplesWithTumor),]$Target.Ki67.BTRC.Count,'medial')
colnames(toBind) <- c('Ki67','aspect')
toPlot <- rbind(toPlot, toBind)
ggplot(data = toPlot, aes(y=Ki67, x=aspect, fill=aspect)) + 
  geom_boxplot(position="dodge", outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_compare_means(aes(group = aspect), label = "p.format", method='wilcox.test')
# CD163 boxplot - because not about only tumor cells, used all samples
toPlot <- data.frame(data[which(data$Color %in% lateral_samples),]$CD163,'lateral')
colnames(toPlot) <- c('CD163','aspect')
toBind <- data.frame(data[which(data$Color %in% medial_samples),]$CD163,'medial')
colnames(toBind) <- c('CD163','aspect')
toPlot <- rbind(toPlot, toBind)
ggplot(data = toPlot, aes(y=CD163, x=aspect, fill=aspect)) + 
  geom_boxplot(position="dodge", outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_compare_means(aes(group = aspect), label = "p.format", method='wilcox.test', method.args = list(alternative = "greater"))
# CD163 boxplot - because not about only tumor cells, used all samples
toPlot <- data.frame(data[which(data$Color %in% lateral_samples),]$ps6,'lateral')
colnames(toPlot) <- c('ps6','aspect')
toBind <- data.frame(data[which(data$Color %in% medial_samples),]$ps6,'medial')
colnames(toBind) <- c('ps6','aspect')
toPlot <- rbind(toPlot, toBind)
ggplot(data = toPlot, aes(y=ps6, x=aspect, fill=aspect)) + 
  geom_boxplot(position="dodge", outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_compare_means(aes(group = aspect), label = "p.format", method='wilcox.test', method.args = list(alternative = "greater"))
ggplot(data, aes(x=aspect, y=ps6, fill=aspect)) + 
  geom_boxplot(fill="white")+
  geom_dotplot(binaxis='y', stackdir='center')+
  scale_fill_manual(values=c("blue", "red"))+
  theme(axis.text.x = element_text(size=12), axis.title = element_text(size = 12), axis.text.y = element_text(size=12), panel.background = element_rect(fill = 'white', colour = 'black'), legend.background=element_blank())+
  stat_compare_means(aes(group = aspect), label = "p.format", method='wilcox.test', method.args = list(alternative = "greater"))


