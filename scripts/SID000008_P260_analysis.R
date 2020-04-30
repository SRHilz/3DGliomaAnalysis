# Created: 2018.05.28
# By: Stephanie R Hilz
# Usage: Perform all final analyses and generate plots for P260

library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

# define file paths
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')
path <- '/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Projects/2016_LoglioExomeAnalysis/Data/'
cnFile <- paste0(dataPath, 'Patient260.GeneCopyNumber.rds')
cpmFile <- paste0(dataPath, 'SID000003_20190529_first_submission.symbol.coding.CPMs.csv')
ssgseaInput <-  paste0(dataPath,'Patient260_ssGSEA_results.txt')
ciberSortFile <- paste0(dataPath, '20190529_first_submission.symbol.coding.CPMs.CIBERSORT.Output_Job8.txt')
mutfile <- paste0(dataPath, 'Patient260.R.mutations.avf.txt')
PDGFRA_lateral_id <- 'PDGFRA_chr4_55133815_C_A'
PDGFRA_medial_id <- 'PDGFRA_chr4_55138638_G_T'
IDH_id <- 'IDH1_chr2_209113112_C_T'
  
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
lateral_high_purity <- c('v1','v2','v3','v4','v5','v6','v7','v8','v9','v10')
medial_high_purity <- c('v12','v13','v14','v18','v19')
medail_low_purity <- c('v11','v15','v16','v17','v20')
samplesToPlot <- c(lateral_high_purity,medial_high_purity)

# read in different data types
cpms <- read.csv(cpmFile, row.names=1)
cn <- readRDS(cnFile)
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

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
boxplot(PDGFRAcpms[paste0('P260',lateral_high_purity)],PDGFRAcpms[paste0('P260',medial_high_purity)], col=c('blue','red'))
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

# pull out purity info for P260
data$purity <- 2*data$IDH1_VAF
dataPurity <- data[which(data$Patient=="P260"),c('SampleName','purity')]
purity <- dataPurity[which(dataPurity$SampleName %in% lateral_high_purity | dataPurity$SampleName %in% medial_high_purity),]$purity
names(purity) <- c(lateral_high_purity, medial_high_purity)

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

# create PTEN boxplot
par(mfrow=c(1,1))
boxplot(PTENcn[lateral_high_purity],PTENcn[medial_high_purity], col=c('blue','red'))
boxplot(PTENcpm[lateral_high_purity],PTENcpm[medial_high_purity], col=c('blue','red'))
wilcox.test(PTENcpm[lateral_high_purity],PTENcpm[medial_high_purity])

# create ssgsea pathway plot for MTOR and ERK
ssgsea <- read.table(ssgseaInput, row.names=1, header=T, sep='\t')
rownames(ssgsea) <- gsub('X','v', rownames(ssgsea))
ofInterest <- c('REACTOME_MTORC1_MEDIATED_SIGNALLING',
                'REACTOME_PROLONGED_ERK_ACTIVATION_EVENTS')
ssgsea <- ssgsea[samplesToPlot,]
toPlot <- ssgsea[,colnames(ssgsea) %in% ofInterest]
colnames(toPlot)[which(colnames(toPlot)=='REACTOME_MTORC1_MEDIATED_SIGNALLING')] <- 'MTOR Sig.'
colnames(toPlot)[which(colnames(toPlot)=='REACTOME_PROLONGED_ERK_ACTIVATION_EVENTS')] <- 'ERK Actv.'
heatmap.2(t(toPlot), trace="none", Colv=F, col=brewer.pal(9,"YlGn"), cexRow=.5)
toPlotGG <- toPlot
toPlotGG$sample <- rownames(toPlotGG)
toPlotGG$aspect <- 'lateral'
toPlotGG[which(toPlotGG$sample %in% medial_high_purity),]$aspect <- 'medial'
toPlotGG <- melt(toPlotGG)
colnames(toPlotGG) <- c('sample','aspect','pathway','score')
ggplot(toPlotGG, aes(x=pathway, y=score, fill=aspect)) +
  geom_boxplot(position=position_dodge(0.8)) +
  scale_fill_manual(values=c("blue","red"))+
  labs(list(y = "Enrichment score") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))
for (p in unique(toPlotGG$pathway)){
  print(paste0('test for lateral vs medial for pathway: ' , p))
  toPlotGGPathway <- toPlotGG[which(toPlotGG$pathway == p),]
  scoreLateral <- toPlotGGPathway[which(toPlotGGPathway$aspect == 'lateral'),]$score
  print(scoreLateral)
  scoreMedial <- toPlotGGPathway[which(toPlotGGPathway$aspect == 'medial'),]$score
  print(scoreMedial)
  if (p == 'MTOR Sig.'){
    print(ks.test(scoreMedial, scoreLateral, alternative = 'less'))
  } else {
    print(ks.test(scoreLateral, scoreMedial, alternative = 'less'))
  }
}

# create immune cell plot for immune cell enrichment differences
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
toPlot[which(toPlot$sample %in% medial_high_purity),]$aspect <- 'medial'
ggplot(toPlot, aes(x=cellType, y=score, fill=aspect)) +
  geom_boxplot(position=position_dodge(0.5), width=0.5) +
  scale_fill_manual(values=c("blue","red"))+
  labs(list(y = "Enrichment score") )+
  theme(axis.text.x = element_text(size=10, color="black",angle = 90, hjust = 1),axis.title = element_text(size = 10), axis.text.y = element_text(size=10, color="black"), panel.background = element_rect(fill = 'white', colour = 'black'))
for (p in unique(toPlot$cellType)){
  print(paste0('test for lateral vs medial for pathway: ' , p))
  toPlotCellType <- toPlot[which(toPlot$cellType == p),]
  scoreLateral <- toPlotCellType[which(toPlotCellType$aspect == 'lateral'),]$score
  print(scoreLateral)
  scoreMedial <- toPlotCellType[which(toPlotCellType$aspect == 'medial'),]$score
  print(scoreMedial)
  print(ks.test(scoreLateral, scoreMedial))
}

## Pathology analysis
file <- '/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Projects/2018_Patient260_Case_Study/Data/Pathology_CollectionSpreadsheet.txt'
data <- read.table(file, sep='\t', header = T, stringsAsFactors = F)
data$lesion <- c(rep('lateral',10),rep('medial',10))
data$lesion <- as.factor(data$lesion)
dataTumor <- data[which(data$`Target...tumor.nuclei` > 0),]
par(mfrow=c(1,2))
dotchart(dataTumor[seq(18,1,-1),]$Target.BV.hyperplasia, dataTumor[seq(18,1,-1),]$Color, col=c(rep('red',8),rep('blue',10)), ylab="BV Hyperplasia")
boxplot(dataTumor$Target.Ki67.BTRC.Count~dataTumor$lesion, col=c('blue','red'), ylab="Ki67")
t.test(dataTumor[which(dataTumor$lesion=='lateral'),]$Target.Ki67.BTRC.Count,dataTumor[which(dataTumor$lesion=='medial'),]$Target.Ki67.BTRC.Count, alternative="less")
ggplot(data, aes(x=aspect, y=CD163, fill=aspect)) + 
  geom_boxplot(fill="white")+
  geom_dotplot(binaxis='y', stackdir='center')+
  scale_fill_manual(values=c("blue", "red"))+
  theme(axis.text.x = element_text(size=12), axis.title = element_text(size = 12), axis.text.y = element_text(size=12), panel.background = element_rect(fill = 'white', colour = 'black'), legend.background=element_blank())+
  stat_compare_means(aes(group = aspect), label = "p.format", method='wilcox.test', method.args = list(alternative = "greater"))
ggplot(data, aes(x=aspect, y=ps6, fill=aspect)) + 
  geom_boxplot(fill="white")+
  geom_dotplot(binaxis='y', stackdir='center')+
  scale_fill_manual(values=c("blue", "red"))+
  theme(axis.text.x = element_text(size=12), axis.title = element_text(size = 12), axis.text.y = element_text(size=12), panel.background = element_rect(fill = 'white', colour = 'black'), legend.background=element_blank())+
  stat_compare_means(aes(group = aspect), label = "p.format", method='wilcox.test', method.args = list(alternative = "greater"))


