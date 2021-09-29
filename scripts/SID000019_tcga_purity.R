library(RColorBrewer)
library(ggrepel)
library(stringr)
library(readr)

## PART 1 - Read in purity data for our patients
# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)

# merge by patient ID
merged <- merge(data, subtypedata, by="Patient")

# subset by if has WES data
toUse <- data[which(!is.na(data$WES_ID)),]$WES_ID
merged <- merged[which(merged$WES_ID %in% toUse),]

# specify purity metric to use
merged$purity <- merged$PyClone
merged[merged$PurityEstUsed=='FACETS',]$purity <- merged[merged$PurityEstUsed=='FACETS',]$FACETS

## PART 2 - read in TCGA data
# Read in TCGA-GBM FACETS output
tcgaFile <- file.path('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Projects/2016_LoglioExomeAnalysis/Analysis_and_Results/tcga_gbm/','TCGA_GBM_FACETS.txt')
facets_tcga <- read.table(tcgaFile,sep='\t', header = T)
vial_col <- c()
# as for other analyses, replaced NA in facets with .1, as this is typical when purity is too low to assess
for (i in 1:nrow(facets_tcga)){
  s <- facets_tcga[i,'sample']
  v <- paste0(substr(s, 22, 23), substr(s, 14, 14))
  vial_col <- c(vial_col, v)
  if (is.na(facets_tcga[i,]$purity)){
    facets_tcga[i,]$purity <- 0.1
  }
}
facets_tcga$case <- paste0(str_extract(facets_tcga$sample, 'TCGA.[A-Z0-9]{2}.[A-Z0-9]{4}'), '-', vial_col)
rownames(facets_tcga) <- facets_tcga$case
facets_tcga$sample <- NULL
facets_tcga$case <- NULL

################# Wilcoxon test ########################
amp_path <- '/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Projects/2016_LoglioExomeAnalysis/Tools/3DGliomaAtlas/git_repo/3DGliomaAtlas/data/metadata/'

densityTCGA <- density(facets_tcga$purity)

patients <- unique(merged[which(merged$IDH_Mut==0 & merged$TERT==1 & merged$Tumor=='Primary'),]$Patient) # we want all classic GBMs

avg_pur_df <- data.frame(patient = c(), avg_pur = c(), med_pur = c(), x.avg = c(), x.med = c()) #x is where falls on TCGA density
# as for other analyses, replaced NA in facets with .1, as this is typical when purity is too low to assess
for (p in patients){
  mergedPatient <- merged[which(merged$Patient == p),]
  avg_pur <- mean(mergedPatient$purity)
  med_pur <- median(mergedPatient$purity)
  x.avg <- densityTCGA$y[which.min(abs(densityTCGA$x - avg_pur))]
  x.med <- densityTCGA$y[which.min(abs(densityTCGA$x - med_pur))]
  avg_pur_df <- rbind(avg_pur_df, data.frame(patient = c(p), avg_pur = c(avg_pur), med_pur = c(med_pur), x.avg = c(x.avg), x.med = c(x.med)))
}

wilcox.test(avg_pur_df$avg_pur, facets_tcga$purity, alternative = 'less')
wilcox.test(avg_pur_df$med_pur, facets_tcga$purity, alternative = 'less')

# comparison matrix 
comparison <- cbind(facets_tcga$purity, 'tcga')
nrow(comparison)
toBind <- cbind(avg_pur_df$med_pur, '3DAtlas')
nrow(toBind)
comparison <- rbind(comparison, toBind)
nrow(comparison)
comparison <- as.data.frame(comparison, stringsAsFactors = F)
colnames(comparison) <- c('purity','dataset')
comparison$purity <- as.numeric(comparison$purity)
comparison$dataset <- factor(comparison$dataset, levels=c('tcga','3DAtlas'))
ggplot(data = comparison, aes(y=purity, x=dataset)) + 
  geom_boxplot(position="dodge", outlier.shape = NA, fill='grey') + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color='black'), axis.title = element_text(size = 12, color='black'), axis.text.y = element_text(size=12, color='black'), panel.background = element_rect(fill = 'white', colour = 'black'))+
  stat_compare_means(aes(group = dataset), label = "p.format", method='wilcox.test', method.args = list(alternative = "less"))

# TCGA purity dist with our patients
ggplot(data = facets_tcga, aes(x=purity)) + 
  stat_density(geom="line") + 
  geom_point(aes(x=med_pur, y=x.avg), data=avg_pur_df, size = 0.1) +
  scale_x_continuous(limits = c(0,1)) + xlab('Purity') +
  scale_y_continuous(limits = c(0,3)) + ylab('Density') +
  geom_text_repel(aes(x=med_pur, y=x.med, label=patient), data=avg_pur_df, hjust=1.5, vjust=-1) +
  theme_classic() + expand_limits(x = 0, y = 0)
