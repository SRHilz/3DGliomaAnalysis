# Created: 2018.05.28
# By: Stephanie R Hilz
# Usage: Perform all final analyses and generate plots for P260

library(reshape2)
library(ggplot2)
library(dplyr)

# define file paths
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')
path <- '/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Projects/2016_LoglioExomeAnalysis/Data/'
cnFile <- paste0(dataPath, 'Patient260.GeneCopyNumber.rds')
cpmFile <- paste0(dataPath, 'SID000003_20190529_first_submission.symbol.coding.CPMs.csv')
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



