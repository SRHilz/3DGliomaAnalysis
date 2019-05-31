# 2017.09.09
# S. Hilz
# Description: Takes mutation data + CN data and looks for evidence of IDH1/2 loss in high cancer cell fraction (CCF) samples

# Dependencies

if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")
library(pheatmap)
library(grid)
library(kableExtra)

# Functions
grabCTforChr <- function(files, names, chr){ #funciton that grabs all TCN values from all PSCN files for the specified chromosome
  
  # first we automatically pull out some helpful info about our window sizes and chr size (max, window, and start)
  max = 0
  for (f in 1:length(files)){
    data <- readRDS(files[f])
    dataChr <- data$data[which(data$data$chromosome==chr),]
    window <- dataChr[2,]$x - dataChr[1,]$x
    start <- dataChr[1,]$x
    if (max(dataChr$x) > max){
      max <- max(dataChr$x)
    }
  }
  
  # create data structure to add data from each file to
  positions <- seq(start, max, by = window)
  toReturn <- as.data.frame(positions, stringsAsFactors=F)
  colnames(toReturn) <- c('x')
  toReturn$x <- as.character(toReturn$x)
  
  # add data from each file
  for (f in 1:length(files)){
    print(f)
    data <- readRDS(files[f])
    dataChr <- data$data[which(data$data$chromosome==chr),c('x','CT')]
    dataChr$x <- as.character(dataChr$x)
    print(dataChr[which(dataChr$x=='250000'),])
    toReturn <- merge(toReturn, dataChr, by='x', all.x=T)
    colnames(toReturn)[which(colnames(toReturn)=='CT')] <- names[f]
  }
  
  # format output as matrix
  toReturn$x <- as.numeric(toReturn$x)
  toReturn <- toReturn[order(toReturn$x),]
  rownames(toReturn) <- toReturn$x
  toReturn <- as.matrix(toReturn[2:dim(toReturn)[2]])

  return(toReturn)
}

getInfo <- function(df, samples){
  subset <- df[which(df$algorithm=='MuTect'),]
  if('X-gene' %in% colnames(subset)){
    geneIndex <- grep('X.gene',colnames(subset))
    colnames(subset)[geneIndex] <- 'gene'
  }
  n <- nrow(subset)
  toReturn <- c()
  for (name in samples){
    print(name)
    ref <- subset[,paste0(name,'_ref_Q20reads')]
    alt <- subset[,paste0(name,'_alt_Q20reads')]
    called <- as.character(subset[,paste0(name,'_called')])
    coverage <- ref + alt
    vaf <- alt/coverage
    sampleID <- rep(name,n)
    gene <- makeLongID(subset$gene,subset$contig,subset$position, subset$ref_allele, subset$alt_allele)
    chunk <- cbind(sampleID,gene,coverage,vaf, called)
    toReturn <- rbind(toReturn,chunk)
  }
  toReturn <- as.data.frame(toReturn, stringsAsFactors=F)
  toReturn$coverage <- as.numeric(toReturn$coverage)
  toReturn$vaf <- as.numeric(toReturn$vaf)
  toReturn$sampleID <- as.factor(toReturn$sampleID)
  colnames(toReturn) <- c('sampleID','uniqueID','coverage','vaf','called')
  return(toReturn)
}

determineIfClonal <- function(n, info){
  clonal <- c()
  infoCalled <- info[which(info$called == "True"),]
  for (variant in unique(infoCalled$uniqueID)){
    if (nrow(infoCalled[which(infoCalled$uniqueID == variant),]) == n){
      clonal <- append(clonal, variant)
    }
  }
  return(clonal)
}

getLocusCN <- function(start, end, chr, rdsFile){
  CNdata <- readRDS(rdsFile)$output
  cn <- min(CNdata[which((CNdata$chromosome == chr) & ((CNdata$tcnStart >= start & CNdata$tcnStart <= end) | (CNdata$tcnEnd >= start & CNdata$tcnEnd <= end) | (CNdata$tcnStart <= start & CNdata$tcnEnd >= end))),]$tcnMean)
}

makeLongID <- function(gene, contig, position,refBase, altBase){
  longID <-paste(gene,'_',contig,'_',refBase,position,altBase, sep='')
  return(longID)
}

# user-defined variables
purityCutoff <- .7
start = 209101803 #start of the locus where we want to look for CN change
end = 209116275 #end of the locus where we want to look for CN change
chr = 2 #chromosome of the locus where we want to look for CN change

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

# read in sample data file
data <- read.table(sampleDataFile, sep='\t', header = T, stringsAsFactors = F)

# create unique sample id
data$sampleID <- paste0(tumorType, '-', data$SampleName)#might need to remove hyphen for some

# read in patient + tumor data file
subtypedata <- read.table(patientTumorDataFile, sep='\t', header = T)

# merge by patient ID
merged <- merge(data, subtypedata, by="Patient")

# set colors (+ molsubtype) for plotting (here focused on subtype differences)
merged$molType <- 'IDH-mut_A' # Astro
merged[which(merged$IDH_Mut==0),]$molType <- 'IDH-wt' #red Set
merged[which(merged$X1p19q==1),]$molType <- 'IDH-mut_O' #blue Set1

# specify purity metric to use
merged$purity <- merged$FACETS
merged[which(merged$PurityEstUsed == 'IDH'),]$purity <- 2*merged[which(merged$PurityEstUsed == 'IDH'),]$IDH1_VAF

# filter based on purity and inspect numbers of samples per patient (removes low-purity samples)
merged <- merged[which(merged$purity >= purityCutoff),]
kable(table(merged$Patient)) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold=TRUE, background = "#C2BFBA")

# identify tumor-wide mutations for each patient (definition of tumor-wide = present in all pieces of CCF >= purity cutoff )
toAnalyze <- data.frame(patient=character(), molType=character(), tumorType=character(), sampleID=character(), CCF=numeric(), IDHvar=numeric(), clonalMean=numeric())
for (p in unique(merged$Patient)){
  print(p)
  # get basic data about case
  molType <- unique(merged[which(merged$Patient == p),]$molType)
  tumorType <- as.character(unique(merged[which(merged$Patient == p),]$Tumor))
  samples <- merged[which(merged$Patient==p),]$sample_type # because we have already filtered by CCF above, will only be high CCF samples
  # get mutational data for case
  patientID <- gsub('P','Patient',p)
  mutFile <- paste0(dataPath,'/',patientID,'.R.mutations.avf.txt')
  mutData <- read.csv(mutFile, header=T, comment.char ='!', sep='\t')
  colnames(mutData) <- gsub('[.]','-',colnames(mutData))
  info <- getInfo(mutData,samples)
  sampleN <- length(unique(info$sampleID))
  clonal_uniqueID <- determineIfClonal(sampleN,info)
  IDH <- clonal_uniqueID[grepl('IDH', clonal_uniqueID)]
  for (s in unique(info$sampleID)){ # for each sample
    print(s)
    # calculate IDH vaf if has IDH mutation
    if(identical(IDH, character(0))){ # if IDH is not present as a clonal variant in this tumor
      IDHvaf <- NA
    } else{
      IDHvaf <- info[which(info$sampleID == s & info$uniqueID == IDH),]$vaf
    }
    # pull out CCF for sample
    CCF <- merged[which(merged$Patient == p & merged$sample_type == s),]$purity
    # calculate mean of all tumor-wide mutations in sample
    clonalMean <- mean(info[which(info$sampleID == s & info$uniqueID %in% clonal_uniqueID),]$vaf)
    ## Add in part here for CN - makes sense I think to wait until I have for everything
    rdsFile <- paste0(dataPath, '/',merged[which(merged$Patient == p & merged$sample_type == s),]$CN)
    locusCN <- getLocusCN(start, end, chr, rdsFile)
    toBind <- c(p,molType,tumorType,s,CCF,IDHvaf,clonalMean) 
    toAnalyze <- rbind(toAnalyze, toBind, stringsAsFactors=F)
  }
}
colnames(toAnalyze) <- c("patient", "molType", "tumorType", "sampleID", "CCF", "IDHvar", "clonalMean")




info$tumorWide[which(info)]


files <- c('/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00363_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00347_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00348_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00349_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00350_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00351_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00352_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00353_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,Z00354_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,SF10711_9-1-22_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,SF10711_9-1-46_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,SF10711_9-2-85_vs_Z00346,chrs=1-22,PairedPSCBS.rds',
           '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient300/CopyNumber/Patient300,SF10711_9-2-123_vs_Z00346,chrs=1-22,PairedPSCBS.rds')
names <- c('P','1','2','4','5','6','7','8','10','9a','9b','9c','9d')

# Code Main

# Read in copy number data for files defined above for a given chromosome; returns a matrix
#  for which each column is a sample, and each row is the copy number for a region equal in size
#  to the window specified in the PSCN pipeline.
chrMatrix <- grabCTforChr(files, names, 2)

# Provide upper bound (prevents insanely large skew in heatmap colors toward gains, which can be infinitely large)
chrMatrix[chrMatrix > 4] <- 4

# Plot as heatmap
break1 <- 1.8
break2 <- 2.2
par(mfrow=c(1,1))
xnames <- rownames(chrMatrix)
i=0
for (i in seq_along(xnames)){
  if (!i%%10 == 0){
    xnames[i+1] <- ''
  }
}
my_palette <- colorRampPalette(c("mediumblue", "grey", "orangered4"))(n = 299)
col_breaks = c(seq(min(chrMatrix, na.rm=T),break1,length=100), # for yellow (~1st quartile)
               seq(break1+.001,break2-.001,length=100), # for black (~within 2nd and 3rd quartiles)
               seq(break2,max(chrMatrix, na.rm=T),length=100)) # for blue (~4th quartile)
hm.parameters <- list(chrMatrix,cellwidth=4, cellheight=.08, fontsize=1, labels_row=xnames, border_col=NA, cluster_rows = F, cluster_cols = F, col= my_palette, breaks = col_breaks, na_col='white')
do.call("pheatmap", hm.parameters)



# Make top plot
means <- colMeans(chrMatrix,na.rm=T)
par(mfrow=c(1,1),mar=rep(4,4))
plot(seq(means),means, type='l', xaxt='n')
axis(1, at = seq(1, length(means), by = 1), las=2)

# Make top plot specific to a region - here we just take the most altered segment in that region



# Confirm location of centromere matches the official "93.3 Mb" location
plot(names(chrMatrix[,1]),is.na(chrMatrix[,1]), cex=.1)
abline(v=93900000, col='red')

