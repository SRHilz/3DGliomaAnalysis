# a script that takes a tumor volume and a set of samples, and determines the volume that encompasses the samples,
# allowing for a 20mm radius around samples; uses convex hull package
library(oro.dicom)
library(oro.nifti)
library(misc3d)
library(geometry)
library(gtools)
library(RNifti)
library(brainR)

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
    text3d(which(sampleNifti@.Data == 1, arr.ind=TRUE)[1,], texts = name, cex=1, adj=-.4)
  }
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

calcSampleCenters <- function(orderedSubdirectories, names=NULL){
  sampleCenteredXYZ <- c()
  for (i in 1:length(orderedSubdirectories)){
    print(paste('Processing ',orderedSubdirectories[i],sep=''))
    sampleDicom <- readDICOM(orderedSubdirectories[i], recursive = FALSE, exclude = NULL, verbose = TRUE)
    sampleNifti <- dicom2nifti(sampleDicom,datatype=4, mode="integer")
    sampleNifti@.Data <- correctOrientation(sampleNifti@.Data)
    sampleXYZ <- which(sampleNifti==1, arr.ind=TRUE)
    sampleCenteredX <- floor(max(sampleXYZ[,1])-((max(sampleXYZ[,1])-min(sampleXYZ[,1]))/2))
    sampleCenteredY <- floor(max(sampleXYZ[,2])-((max(sampleXYZ[,2])-min(sampleXYZ[,2]))/2))
    sampleCenteredZ <- floor(max(sampleXYZ[,3])-((max(sampleXYZ[,3])-min(sampleXYZ[,3]))/2))
    sampleCenteredXYZ <- rbind(sampleCenteredXYZ, c(sampleCenteredX,sampleCenteredY,sampleCenteredZ))
  }
  return(sampleCenteredXYZ)
}

create6Pointers <- function(orderedSubdirectories, samples, hullRadiusBuffer){# creates the sample center coordinates + 20 mm spacers for chull input
  converted6Point <- matrix(, ncol=3)
  sampleCenteredXYZ <- c()
  otherSampleCenteredXYZ <- c()
  for (i in 1:length(orderedSubdirectories)){
    print(paste('Processing ',orderedSubdirectories[i],sep=''))
    sampleID <- gsub('sample','',strsplit(subdirectories[i],'/')[[1]][length(strsplit(subdirectories[i],'/')[[1]])])
    print(sampleID)
    sampleDicom <- readDICOM(orderedSubdirectories[i], recursive = FALSE, exclude = NULL, verbose = TRUE)
    sampleNifti <- dicom2nifti(sampleDicom,datatype=4, mode="integer")
    sampleNifti@.Data <- correctOrientation(sampleNifti@.Data)
    sampleXYZ <- which(sampleNifti==1, arr.ind=TRUE)
    sampleCenteredX <- max(sampleXYZ[,1])-((max(sampleXYZ[,1])-min(sampleXYZ[,1]))/2)
    sampleCenteredY <- max(sampleXYZ[,2])-((max(sampleXYZ[,2])-min(sampleXYZ[,2]))/2)
    sampleCenteredZ <- max(sampleXYZ[,3])-((max(sampleXYZ[,3])-min(sampleXYZ[,3]))/2)
    referenceNifit <- retrieveNifti(sampleNifti)
    convertedCentered <- voxelToWorld(c(sampleCenteredX,sampleCenteredY,sampleCenteredZ), referenceNifit)
    if (sampleID %in% samples){
      print('here1')
      sampleCenteredXYZ <- rbind(sampleCenteredXYZ, c(convertedCentered))
      for (d in 1:3){
        expansion <- matrix(rep(convertedCentered,2*hullRadiusBuffer+1),nrow=3)
        expansion[d,] <- expansion[d,] + seq(-hullRadiusBuffer,hullRadiusBuffer,1)
        converted6Point <- rbind(converted6Point, t(expansion) )
      }
    } else {
      print('here2')
      otherSampleCenteredXYZ <- rbind(otherSampleCenteredXYZ, c(convertedCentered))
    }
  }
  output <- list('converted6Point'= converted6Point, 'centers' = sampleCenteredXYZ, 'otherCenters' = otherSampleCenteredXYZ, 'niftiRef' = referenceNifit)
  return(output)
}

convertCartesianToVector <- function(x,y,z,dim1,dim2){
  ((dim1*dim2)*(z-1))+(dim1*(y-1))+x
}

pairwiseDist2pointsXYZ <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + ((1*(p2[3]-p1[3]))^2) ))
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

#Gather input
# Specify patient and load the config file for that patient. Config file contains paths to imaging files + ordering of samples + sample names + colors
patientID <- 'Patient454'
sf <- 'sf11731'
hullRadiusBuffer <- 10 # in mm, the extra buffer you want around the samples for your hull
samples <- c(3,4,6)
source(paste('/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/ConfigFiles/configImagingData.R',sep=''))

# Bring in subdirectories each containing a set of dicom files for a given sample; check that order is same as order of colors you want to plot
subdirectories <- mixedsort(list.dirs(sampleDicomsMainDir,recursive=FALSE))#this should be a list of all subdirectories in the main sample dir, each containing dicoms for a single data point

# Add 20 mm around each sample
output <- create6Pointers(subdirectories, samples, hullRadiusBuffer)
expandedSamples <- output$converted6Point
expandedSamples <- as.matrix(expandedSamples[complete.cases(expandedSamples),])
sampleCentersLPS <- output$centers
otherSampleCentersLPS <- output$otherCenters
niftiRef <- output$niftiRef

# Plot sample centers in LPS
rgl.clear()
for (i in 1:nrow(sampleCentersLPS)){
  points3d(x=sampleCentersLPS[i,1], y=sampleCentersLPS[i,2], z=sampleCentersLPS[i,3], cex=2, col='blue')
  text3d(x=sampleCentersLPS[i,1], y=sampleCentersLPS[i,2], z=sampleCentersLPS[i,3], texts = samples[i], cex=1, adj=-.3)
}
for (i in 1:nrow(otherSampleCentersLPS)){
  points3d(x=otherSampleCentersLPS[i,1], y=otherSampleCentersLPS[i,2], z=otherSampleCentersLPS[i,3], cex=2, col='grey')
}

# Create hull, calculate volume in real-world mm
hull.v <- t(convhulln(expandedSamples,options="FA"))
print(hull.v)
hull <- t(convhulln(expandedSamples))

# visualize LPS-space hull
triangles3d(expandedSamples[hull,1],expandedSamples[hull,2],expandedSamples[hull,3],col="blue",alpha=.2)

# Create a hull object that is instead viewable on top of NIFTI pixels
expandedSamplesVoxels <- worldToVoxel(expandedSamples, niftiRef)
hull.pixels <- t(convhulln(expandedSamplesVoxels))

# View with tumor
# Plot background of brain and tumor
plotTemplate(brainDicoms,tumorDicoms,showBrain=F)
# Plot samples
colors <- rep('grey',length(subdirectories))
colors[samples] <- 'red'
plot3dSamples(subdirectories, colors, sampleID_short)
# add hull
triangles3d(expandedSamplesVoxels[hull.pixels,1],expandedSamplesVoxels[hull.pixels,2],expandedSamplesVoxels[hull.pixels,3],col="blue",alpha=.2)

# Save model as html
writeWebGL_split(dir=getwd(), filename = paste(patientID,'_',paste(samples,collapse='_'), '_chull.html',sep=''), width=500,writeIt=TRUE)

# Save model as movie!
movie3d(spin3d(), duration=7, movie = paste(patientID,paste(samples,collapse=''), '_chull',sep=''))
