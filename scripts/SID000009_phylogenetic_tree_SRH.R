# Created: 2017.08.15
# By: Stephanie R Hilz
# Usage: Perform phylogenetic analysis from avf files

library(ape)

makeLongID <- function(gene, contig, position,refBase, altBase){
  longID <-paste(gene,'_',contig,'_',refBase,position,altBase, sep='')
  return(longID)
}

getEdgeMembership <- function(df){
  toReturn <- c()
  for (i in seq_len(nrow(df))){
    mut <- df[i,]$uniqueID
    subset <- df[which(df$uniqueID==mut),]
    groupMembership <- paste(unique(subset$sample_type),collapse=',')
    toReturn <- append(toReturn, groupMembership)
  }
  return(toReturn)
}

getEdgeCount <- function(df){
  toReturn <- c()
  for (i in seq_len(nrow(df))){
    edge <- df[i,]$edge_members
    subset <- df[which(df$edge_members==edge),]
    muts_in_edge <- length(unique(subset$uniqueID))
    toReturn <- append(toReturn, muts_in_edge)
  }
  return(toReturn)
}

## Changes since original version taken from Tali but prior to using git
  # added section that allows filtering based off of IGV manual checking from an IGV annotation file
  # added section that outputs three additional columns - mut_freq, edge_membership, and edge_mut_count - to mutationstree file
  # added section that allows filtering based off of automated exome quality filtering (avf)

## Define which patient to make tree for

# read in config file info
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')

## File I/O setup
patientID <- 'Patient450'
avfFile <- paste0(patientID,'.R.mutations.avf.txt')
avfPath <- paste0(dataPath,avfFile)
outfolder <- 'SID000009_phylogenetic_tree/'
outpdf <- paste0(outputPath, outfolder, sub(".txt", "tree.pdf", sub("Variants","PhylogeneticTrees",avfFile)))
outtxt <- paste0(outputPath,outfolder, sub(".txt", "tree.txt", sub("Variants","PhylogeneticTrees",avfFile)))

## read in data
muts <- read.delim(avfPath, as.is=TRUE)
if('X.gene' %in% colnames(muts)){
  geneIndex <- grep('X.gene',colnames(muts))
    colnames(muts)[geneIndex] <- 'gene'
}

## check if mutation file has sanger status & mutect covered in all columns
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
# finally filter out mutations flagged by avf
muts <- muts[which(muts$decision=='retain'),]

## pull out the relevant columns
if(sanger) { calls <- which(sapply(names(muts), function(x) { regexpr("_called_sanger$", x) != -1} )) ## these columns contain the sanger-corrected mutations calls
} else { calls <- which(sapply(names(muts), function(x) { regexpr("_called$", x) != -1} )); calls <- setdiff(calls, which(names(muts)=="samples_called")) }
calls <- calls[order(names(muts)[calls])]
muts.bin <- muts[ , calls]
muts$uniq <- paste(muts$gene, muts$contig, muts$position, muts$ref_allele, muts$alt_allele, sep="_")
rownames(muts.bin) <- muts$uniq

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

## add a germline "normal" sample
muts.norm <- muts.bin
muts.norm$Normal <- FALSE
names(muts.norm) <- sapply(names(muts.norm), function(x) { strsplit(x, "_called")[[1]][1] } )

## save the list of mutations that are actually used to build the tree
names(muts.bin) <- sapply(names(muts.bin), function(x) { strsplit(x, "_called")[[1]][1] } )
cols.to.keep <- c(1:9)
muts.final <- data.frame(matrix(ncol=11,nrow=0))
names(muts.final) <- c(names(muts)[cols.to.keep], "patient_ID", "sample_type")

for(rowID in 1:nrow(muts.bin)) {
	colIDs <- which(as.logical(muts.bin[rowID, ]))
	for(colID in colIDs) {
		newrow <- nrow(muts.final) + 1
		muts.final[newrow, 1:9] <- muts[which(muts$uniq == rownames(muts.bin)[rowID]), cols.to.keep]
		muts.final[newrow, 11] <- names(muts.bin)[colID]
	}
}
patID <- gsub('Patient','P',patientID)
muts.final$patient_ID <- patID
muts.final$uniqueID <- makeLongID(muts.final$gene,muts.final$contig,muts.final$position, muts.final$ref_allele, muts.final$alt_allele)
muts.final <- transform(muts.final, mut_freq = ave(seq(nrow(muts.final)), uniqueID, FUN=length))
muts.final$edge_members <- getEdgeMembership(muts.final)
muts.final$edge_mut_count <- getEdgeCount(muts.final)
muts.final$uniqueID <- NULL
write.table(muts.final, file=outtxt, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#colnames(muts.norm) <- gsub('Primary.v','',colnames(muts.norm)) #only use if want to simplify sample names for multi-sample cases
## print trees to file
pdf(outpdf)
dd <- dist(t(as.matrix(muts.norm)), method="manhattan")
tt <- fastme.ols(dd)
rr <- root(tt, "Normal")
#plot(rr, "c")
plot(rr, "u")  ## "u" draws branch lengths to scale!
text(0, .8, write.tree(rr, digits=2), cex=0.5, pos=4)
dev.off()

save(dd, file=sub(".pdf",".dd.Rdata", outpdf))

### for HM it's sometimes nice to draw a tree with modified branch lengths:
## print "a" to the screen
## save a modified version back as "a" after shortening the very long branches
## save that and as a short tree
## edit in illustrator to make clear that branches aren't to scale:
a <- write.tree(rr, digits=3)
pdf(sub("tree", "tree.short", outpdf))
plot(read.tree(text=a),"u")
text(0, .8, a, cex=0.5, pos=4)
dev.off()
