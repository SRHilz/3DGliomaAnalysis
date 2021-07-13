library(oro.dicom)
library(oro.nifti)
library(misc3d)
library(gtools)
library(rgl)
library(brainR)
library(gtools)
library(dplyr)


plotTemplate <- function(tumorModel){
  print('Processing tumor model')
  dtemp <- dim(tumorModel)
  print('Creating tumor contour and plotting tumor')
  tumor <- contour3d(tumorModel, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = 1, alpha = .2, add = FALSE, draw = TRUE, color = 'yellow')
}

definePeriph <- function(tumorModel){
  tumorMaskPositive <- max(tumorModel)
  peripheryMask <- tumorModel 
  peripheryMask[which(peripheryMask <= tumorMaskPositive)] <- 0
  peripheryXYZ <- c()
  for (z in 1:dim(tumorModel)[3]){
    for (y in 1:dim(tumorModel)[2]){
      for (x in 1:dim(tumorModel)[1]){
        if (tumorModel[x,y,z] == tumorMaskPositive){
          if (tumorModel[x+1,y,z] == 0 ||
              tumorModel[x-1,y,z] == 0 ||
              tumorModel[x,y+1,z] == 0 ||
              tumorModel[x,y-1,z] == 0 ||
              tumorModel[x,y,z+1] == 0 ||
              tumorModel[x,y,z-1] == 0){
            peripheryXYZ <- rbind(peripheryXYZ, c(x,y,z))
            peripheryMask[x,y,z] <- 1
          }
        }
      }
    }
  }
  peripheryXYZ <- as.data.frame(peripheryXYZ)
  colnames(peripheryXYZ) <- c('x','y','z')
  output <- list('XYZ' = peripheryXYZ, 'mask' = peripheryMask )
  return(output)
}

plot3dSamples <- function(sampleCoordinates, colors){
  for (name in rownames(sampleCoordinates)){
    sampleColors <- colors[name]
    x <- sampleCoordinates[name,]
    print(paste0('Processing sample',name))
    points3d(x=sampleCoordinates[,1], y=sampleCoordinates[,2], z=sampleCoordinates[,3], level = 1, size = 7, color=sampleColors)
    text3d(x=sampleCoordinates[,1], y=sampleCoordinates[,2], z=sampleCoordinates[,3], texts = rownames(sampleCoordinates), cex=1, adj=-.3)
  }
}

pairwiseDist2pointsXYZ <- function(p1, p2, adj){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((adj$x*(p2[1]-p1[1]))^2) + ((adj$y*(p2[2]-p1[2]))^2) + ((adj$z*(p2[3]-p1[3]))^2) ))
}

calcSamplePairwiseDistances <- function(samplePoints, adj){
  toReturn <- c()
  for(i in 1:dim(sampleCenters)[1]){
    forSample <- c()
    for(j in 1:dim(sampleCenters)[1]){
      forSample <- append(forSample, pairwiseDist2pointsXYZ(sampleCenters[i,],sampleCenters[j,], adj))
    }
    toReturn <- rbind(toReturn, forSample)
  }
  return(toReturn)
}

calcMinDistances <- function(samplePoints, peripheryPoints, adj){
  minDistAllSamples <- c()
  minPointAllSamples <- c()
  for (s in 1:dim(samplePoints)[1]){
    print(s)
    distances <- c()
    for (p in 1:dim(peripheryPoints)[1]){
      d <- pairwiseDist2pointsXYZ(samplePoints[s,], peripheryPoints[p,], adj)
      distances <- append(distances, d)
    }
    distances <- distances %>% unlist
    minDistAllSamples <- append(minDistAllSamples, min(distances))
    minPointAllSamples <- rbind(minPointAllSamples, peripheryPoints[which(distances==min(distances))[1],])#there may be >1 point that is the min, this will just take the first
  }
  minPointAllSamples <- as.data.frame(minPointAllSamples)
  colnames(minPointAllSamples) <- c('x','y','z')
  output <- list('distances' = minDistAllSamples, 'XYZ' = minPointAllSamples )
  return(output)
}

plotDistanceSegments <- function(periphCoord, sampleCoord){
  x <- c()
  y <- c()
  z <- c()
  for (i in 1:dim(periphCoord)[1]){
    x <- cbind(x, c(sampleCoord[i,1], periphCoord[i,1]))
    y <- cbind(y, c(sampleCoord[i,2], periphCoord[i,2]))
    z <- cbind(z, c(sampleCoord[i,3], periphCoord[i,3]))
  }
  segments3d(x=x, y=y, z=z, col='red')
  points3d(x=x[2,], y=y[2,], z=z[2,], cex=2, col='red')
}

plotDistanceSegmentsCentroid <- function(centroid, sampleCoord){
  x <- c()
  y <- c()
  z <- c()
  for (i in 1:dim(sampleCoord)[1]){
    x <- cbind(x, c(sampleCoord[i,1], centroid[1]))
    y <- cbind(y, c(sampleCoord[i,2], centroid[2]))
    z <- cbind(z, c(sampleCoord[i,3], centroid[3]))
  }
  segments3d(x=x, y=y, z=z, col='red')
}

defineCentroid <- function(peripheryPoints){
  x <- mean(peripheryPoints$x)
  y <- mean(peripheryPoints$y)
  z <- mean(peripheryPoints$z)
  return(c(x,y,z))
}

calcDistancesCentroid <- function(samplePoints, centroid, adj){
  distAllSamples <- c()
  for (s in 1:dim(samplePoints)[1]){
    print(s)
    d <- pairwiseDist2pointsXYZ(samplePoints[s,], centroid, adj)
    distAllSamples <- append(distAllSamples, d)
  }
  return(distAllSamples)
}

# Specify patient and load the config file for that patient. Config file contains paths to imaging files + ordering of samples + sample names + colors
patientID <- 'Patient507'
sf <- "sf12446"
tumorPrefix <- "Primary-v"
source('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/scripts/studyConfig.R')
modelsPath <- paste0(dataPath, '/3Dmodels/',patientID,'/',sf)

# Retreive sample model files
sampleModelFiles <- mixedsort(list.files(modelsPath, pattern="coordinates_samples.rds"))

# Read in sample models
sampleCenters <- lapply(paste0(modelsPath, '/', sampleModelFiles), readRDS)[[1]]#this should be a list of all subdirectories in the main sample dir, each containing dicoms for a single data point

# Create names object
names <- rownames(sampleCenters)

# Also generate a longer names object that will match the exome sample identifier
sampleID_long <- paste0(tumorPrefix,names)

# Set up default colors
colors <- rep('blue', length(names))
names(colors) <- names

# Read in tumor model for patient
tumorModel <- readRDS(paste0(modelsPath, '/tumor_t2.rds'))

# bring in adjustemnts (i.e. spacing for pixels) for patient
adj <- readRDS(paste0(modelsPath, '/adj.rds'))

# Plot background of brain and tumor
plotTemplate(tumorModel) 

# Define periphery
tumorPeriphery <- definePeriph(tumorModel)

# View periphery
#contour3d(tumorPeriphery$mask, x=1:dim(tumorPeriphery$mask)[1], y=1:dim(tumorPeriphery$mask)[2], level=1, z=1:dim(tumorPeriphery$mask)[3], alpha = 1, add = TRUE, draw = TRUE, color = 'red')

# Plot samples
plot3dSamples(sampleCenters, colors)

# Calc sample pairwise distances - used only to check, should always use official coordinates
samplePairwiseDistances <- calcSamplePairwiseDistances(sampleCenters, adj)

# Calculate avg distance from periphery
print(mean(samplePairwiseDistances))

# Calc min distance
minDistances <- calcMinDistances(sampleCenters, tumorPeriphery$XYZ, adj)

# Define centroid
tumorCentroid <- defineCentroid(tumorPeriphery$XYZ)

# View centroid
points3d(x=tumorCentroid[1], y=tumorCentroid[2], z=tumorCentroid[3], level=1, size=13, alpha = 1, add = TRUE, draw = TRUE, color = 'red')

# Calc distances from centroid
distancesCentroid <- calcDistancesCentroid(sampleCenters, tumorCentroid, adj)

# Plot segment between sample and closest tumor edge
plotDistanceSegments(minDistances$XYZ, sampleCenters)

# Plot segment between sample and centroid
plotDistanceSegmentsCentroid(tumorCentroid, sampleCenters)

# Calculate avg distance from periphery
print(mean(minDistances$distances))

# Create data frame for distance from periphery
distanceFromPeriphery <- data.frame(sampleID_long=sampleID_long, distance=minDistances$distances)

# Save distance from tumor periphery as R object
save(distanceFromPeriphery,file=paste('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/Coordinates/DistanceFromPeriphery.rda',sep=''))

# Save distance from tumor periphery as txt file
write.table(distanceFromPeriphery, file=paste0('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/Coordinates/DistanceFromPeriphery.txt'),sep='\t', row.names=F, quote=F)

# Calculate avg distance from centroid
print(mean(distancesCentroid))

# Create data frame for distance from centroid
distanceFromCentroid <- data.frame(sampleID_long=sampleID_long, distance=distancesCentroid)

# Save distance from centroid as R object
save(distanceFromCentroid,file=paste('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/Coordinates/DistanceFromCentroid.rda',sep=''))

# Save distance from tumor centoid as txt file
write.table(distanceFromCentroid, file=paste0('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/Coordinates/DistanceFromCentroid.txt'),sep='\t', row.names=F, quote=F)

# Save model as html
writeWebGL_split(dir=getwd(), filename = paste(patientID,'.html',sep=''), width=500,writeIt=TRUE)

# Save model as movie!
movie3d(spin3d(), duration=7, movie = patientID, dir=paste('/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/Imaging/output',sep=''))
