library(oro.dicom)
library(oro.nifti)
library(misc3d)
library(gtools)
library(rgl)
library(brainR)
library(gtools)
library(dplyr)


plotTemplate <- function(pathBrainDicoms, pathTumorDicoms,showBrain=TRUE){
  if (showBrain==TRUE){
    print(paste('Processing ',pathBrainDicoms,sep=''))
    brainTemplate <- readDICOM(pathBrainDicoms, recursive = FALSE, exclude = NULL, verbose = TRUE)
    brainTemplate <- dicom2nifti(brainTemplate,datatype=4, mode="integer")
    brainTemplate@.Data <- correctOrientation(brainTemplate@.Data)
    dtemp <- dim(brainTemplate)
    print('Creating brain contour and plotting brain.')
    brainTemplate[brainTemplate<1000] = 0
    brainTemplate[brainTemplate>=1000] = 1
    brain <- contour3d(brainTemplate, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = 1, alpha = .05, color='grey', draw = TRUE)
    print(paste('Processing ', pathTumorDicoms, sep=''))
    tumorTemplate <- readDICOM(pathTumorDicoms, recursive = TRUE, exclude = NULL, verbose = TRUE)
    tumorTemplate <- dicom2nifti(tumorTemplate,datatype=4, mode="integer")
    tumorTemplate@.Data <- correctOrientation(tumorTemplate@.Data)
    print('Creating tumor contour and plotting tumor')
    tumor <- contour3d(tumorTemplate, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = tumorTemplate@cal_max, alpha = .2, add = TRUE, draw = TRUE, color = 'yellow')
  } else {
    print(paste('Processing ', pathTumorDicoms, sep=''))
    tumorTemplate <- readDICOM(pathTumorDicoms, recursive = TRUE, exclude = NULL, verbose = TRUE)
    tumorTemplate <- dicom2nifti(tumorTemplate,datatype=4, mode="integer")
    tumorTemplate@.Data <- correctOrientation(tumorTemplate@.Data)
    dtemp <- dim(tumorTemplate)
    print('Creating tumor contour and plotting tumor')
    tumor <- contour3d(tumorTemplate, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = tumorTemplate@cal_max, alpha = .2, add = FALSE, draw = TRUE, color = 'yellow')
  }
}

addCavity <- function(pathCavity, pathTumorDicoms, trim=TRUE){
  print(paste('Processing ',pathCavity,sep=''))
  cavityTemplate <- readDICOM(pathCavity, recursive = TRUE, exclude = NULL, verbose = TRUE)
  cavityTemplate <- dicom2nifti(cavityTemplate,datatype=4, mode="integer")
  cavityTemplate@.Data <- correctOrientation(cavityTemplate@.Data)
  dtemp <- dim(cavityTemplate)
  if (trim==TRUE){
    print(paste('Processing ', pathTumorDicoms, sep=''))
    tumorTemplate <- readDICOM(pathTumorDicoms, recursive = TRUE, exclude = NULL, verbose = TRUE)
    tumorTemplate <- dicom2nifti(tumorTemplate,datatype=4, mode="integer")
    tumorTemplate@.Data <- correctOrientation(tumorTemplate@.Data)
    dtumor <- dim(tumorTemplate)
    if(dtumor[1]==dtemp[1] & dtumor[2]==dtemp[2] & dtumor[3]==dtemp[3]){
      cavityTemplate[which(tumorTemplate==0)] <- 0
    }else{print('ERROR - Tumor dicom not of same dimensions as cavity dicom')}
  }
  print('Adding cavity.')
  cavity <- contour3d(cavityTemplate, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = cavityTemplate@cal_max, add = TRUE, draw = TRUE, alpha = 0.2, color='blue')
}

definePeriph <- function(pathTumorDicoms){
  print(paste('Processing ', pathTumorDicoms, sep=''))
  tumorTemplate <- readDICOM(pathTumorDicoms, recursive = TRUE, exclude = NULL, verbose = TRUE)
  tumorTemplate <- dicom2nifti(tumorTemplate,datatype=4, mode="integer")
  tumorTemplate@.Data <- correctOrientation(tumorTemplate@.Data)
  tumorNifti <- tumorTemplate@.Data
  tumorMaskPositive <- max(tumorTemplate@.Data)
  peripheryMask <- tumorNifti 
  peripheryMask[which(peripheryMask <= tumorMaskPositive)] <- 0
  peripheryXYZ <- c()
  for (z in 1:dim(tumorNifti)[3]){
    for (y in 1:dim(tumorNifti)[2]){
      for (x in 1:dim(tumorNifti)[1]){
        if (tumorNifti[x,y,z] == tumorMaskPositive){
          if (tumorNifti[x+1,y,z] == 0 ||
              tumorNifti[x-1,y,z] == 0 ||
              tumorNifti[x,y+1,z] == 0 ||
              tumorNifti[x,y-1,z] == 0 ||
              tumorNifti[x,y,z+1] == 0 ||
              tumorNifti[x,y,z-1] == 0){
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

plot3dSamples <- function(orderedSubdirectories, colors, names=NULL){
  for (i in 1:length(orderedSubdirectories)){
    print(paste('Processing ',orderedSubdirectories[i],sep=''))
    if (is.null(names)){
      splitname <- strsplit(orderedSubdirectories[i],'/') %>% unlist
      name <- splitname[length(splitname)]
      name <- strsplit(name,'\\.')[1] %>% unlist
    } else {
      name <- names[i]
    }
    sampleDicom <- readDICOM(orderedSubdirectories[i], recursive = FALSE, exclude = NULL, verbose = TRUE)
    sampleNifti <- dicom2nifti(sampleDicom,datatype=4, mode="integer")
    sampleNifti@.Data <- correctOrientation(sampleNifti@.Data)
    dtemp <- dim(sampleNifti)
    sample <- contour3d(sampleNifti, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = sampleNifti@cal_max, alpha = 1, color=colors[i], add = TRUE, draw = TRUE)
    text3d(which(sampleNifti@.Data == 1, arr.ind=TRUE)[1,], texts = name, cex=1, adj=-.3)
  }
}

calcSampleCenters <- function(orderedSubdirectories){
  sampleCenteredXYZ <- c()
  for (i in 1:length(orderedSubdirectories)){
    print(paste('Processing ',orderedSubdirectories[i],sep=''))
    splitname <- strsplit(orderedSubdirectories[i],'/') %>% unlist
    name <- splitname[length(splitname)]
    name <- gsub('sample','',name)
    sampleDicom <- readDICOM(orderedSubdirectories[i], recursive = FALSE, exclude = NULL, verbose = TRUE)
    sampleNifti <- dicom2nifti(sampleDicom,datatype=4, mode="integer")
    sampleNifti@.Data <- correctOrientation(sampleNifti@.Data)
    sampleXYZ <- which(sampleNifti==1, arr.ind=TRUE)
    sampleCenteredX <- max(sampleXYZ[,1])-((max(sampleXYZ[,1])-min(sampleXYZ[,1]))/2)
    sampleCenteredY <- max(sampleXYZ[,2])-((max(sampleXYZ[,2])-min(sampleXYZ[,2]))/2)
    sampleCenteredZ <- max(sampleXYZ[,3])-((max(sampleXYZ[,3])-min(sampleXYZ[,3]))/2)
    sampleCenteredXYZ <- rbind(sampleCenteredXYZ, c(sampleCenteredX,sampleCenteredY,sampleCenteredZ))
  }
  return(sampleCenteredXYZ)
}

pairwiseDist2pointsXYZ <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + ((1.5*(p2[3]-p1[3]))^2) ))
}

calcSamplePairwiseDistances <- function(samplePoints){
  toReturn <- c()
  for(i in 1:dim(sampleCenters)[1]){
    forSample <- c()
    for(j in 1:dim(sampleCenters)[1]){
      forSample <- append(forSample, pairwiseDist2pointsXYZ(sampleCenters[i,],sampleCenters[j,]))
    }
    toReturn <- rbind(toReturn, forSample)
  }
  return(toReturn)
}

calcMinDistances <- function(samplePoints, peripheryPoints){
  minDistAllSamples <- c()
  minPointAllSamples <- c()
  for (s in 1:dim(samplePoints)[1]){
    print(s)
    distances <- c()
    for (p in 1:dim(peripheryPoints)[1]){
      d <- pairwiseDist2pointsXYZ(samplePoints[s,], peripheryPoints[p,])
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

correctOrientation <-function(a){#does an x and then y axis reflection
#on a 3D array [x,y,z] of imaging data
  for (z in seq_len(dim(a)[3])){
    for (y in seq_len(dim(a)[2])){
      a[,y,z] <- rev(a[,y,z])
    }
  }
  return(a)
}

defineCentroid <- function(peripheryPoints){
  x <- mean(peripheryPoints$x)
  y <- mean(peripheryPoints$y)
  z <- mean(peripheryPoints$z)
  return(c(x,y,z))
}

calcDistancesCentroid <- function(samplePoints, centroid){
  distAllSamples <- c()
  for (s in 1:dim(samplePoints)[1]){
    print(s)
    d <- pairwiseDist2pointsXYZ(samplePoints[s,], centroid)
    distAllSamples <- append(distAllSamples, d)
  }
  return(distAllSamples)
}

# Specify patient and load the config file for that patient. Config file contains paths to imaging files + ordering of samples + sample names + colors
patientID <- 'Patient340'
sf <- "sf11055"
source(paste('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/ConfigFiles/configImagingData.R',sep=''))

# Explore a scan (user should edit this as they wish to explore different dicom files)
featureToExplore <- brainDicoms
toExplore <- readDICOM(featureToExplore, recursive = FALSE, exclude = NULL, verbose = TRUE)
toExplore <- dicom2nifti(toExplore,datatype=4, mode="integer")
image(toExplore)
orthographic(toExplore, xyz=c(90,120,50))#x moves saggital R->L, y coronal P->A, z axial I -> S; so 0,0,0 is RPI

# Plot background of brain and tumor
plotTemplate(brainDicoms,tumorDicoms,showBrain=F) 

# Add in resection cavity (optional)
addCavity(cavityDicoms, tumorDicoms, trim=FALSE)

# Bring in subdirectories each containing a set of dicom files for a given sample; required to view samples; IMPORTANT: check that order is same as order of names and colors (either defined here or in config file) - if not, need to manually order sampleDicomsMainDir in config file
subdirectories <- mixedsort(list.dirs(sampleDicomsMainDir,recursive=FALSE))#this should be a list of all subdirectories in the main sample dir, each containing dicoms for a single data point

# Define periphery
tumorPeriphery <- definePeriph(tumorDicoms)

# View periphery
contour3d(tumorPeriphery$mask, x=1:dim(tumorPeriphery$mask)[1], y=1:dim(tumorPeriphery$mask)[2], level=1, z=1:dim(tumorPeriphery$mask)[3], alpha = 1, add = TRUE, draw = TRUE, color = 'red')

# Extract centered sample coordinates
sampleCenters <- calcSampleCenters(subdirectories)

# Calc sample pairwise distances - used only to check, should always use official coordinates
samplePairwiseDistances <- calcSamplePairwiseDistances(sampleCenters)

# Calculate avg distance from periphery
print(mean(samplePairwiseDistances))

# Calc min distance
minDistances <- calcMinDistances(sampleCenters, tumorPeriphery$XYZ)

# Define centroid
tumorCentroid <- defineCentroid(tumorPeriphery$XYZ)

# View centroid
points3d(x=tumorCentroid[1], y=tumorCentroid[2], z=tumorCentroid[3], level=1, size=13, alpha = 1, add = TRUE, draw = TRUE, color = 'red')

# Calc distances from centroid
distancesCentroid <- calcDistancesCentroid(sampleCenters, tumorCentroid)

# Plot samples
plot3dSamples(subdirectories, colors, sampleID_short)

# Plot segment between sample and closest tumor edge
plotDistanceSegments(minDistances$XYZ, sampleCenters)

# Plot segment between sample and centroid
plotDistanceSegmentsCentroid(tumorCentroid, sampleCenters)

# View sample orthographic (will give warning message In min(x...etc) because putting in text, please ignore!)
sample <- 20
orthographic(toExplore, xyz=sampleCenters[sample,], text=c(paste('sample ',as.character(sample),sep='')))#x is sagital, y is coronal, z is axial

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
