sciClone <- read.table('/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient452/SciClone/P452.cluster.summary.tsv.means.txt',sep='\t', header=T, stringsAsFactors = F, row.names=1)
colnames(sciClone) <- gsub('Primary.','',colnames(sciClone))

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

getCensusGeneAnnotation <- function(genes){
  cancerCensusGenesFile <- '/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/GeneLists/cancer_gene_census.csv'
  censusGenes <- read.csv(cancerCensusGenesFile)
  cancerGeneAnnotationDF <- as.data.frame(genes)
  cancerGeneAnnotationDF$Tier <- 0
  tier2Genes <- censusGenes[which(censusGenes$Tier=='2'),]$Gene.Symbol
  tier1Genes <- censusGenes[which(censusGenes$Tier=='1'),]$Gene.Symbol
  if (any(cancerGeneAnnotationDF$genes %in% tier2Genes)){
    cancerGeneAnnotationDF[cancerGeneAnnotationDF$genes %in% tier2Genes,]$Tier <- 2
  }
  if (any(cancerGeneAnnotationDF$genes %in% tier1Genes)){
    cancerGeneAnnotationDF[cancerGeneAnnotationDF$genes %in% tier1Genes,]$Tier <- 1
  }
  return(cancerGeneAnnotationDF$Tier)
}

sciCloneGenesClusters <- read.table('/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/Patient452/SciClone/P452.cluster.tsv',sep='\t', header=T, stringsAsFactors = F)
sciCloneGenesClusters$sanger <- getCensusGeneAnnotation(sciCloneGenesClusters$gene)
sciCloneGenesClusters$uniq <- paste0(sciCloneGenesClusters$gene, '_chr', sciCloneGenesClusters$chr, '_',sciCloneGenesClusters$st)

patientID <- 'Patient452'
infile <- paste('/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/Variants/',patientID,'.R.mutations.txt',sep='')
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

