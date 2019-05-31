repoPath <- '/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/'
dataPath <- paste0(repoPath, 'data/')
outputPath <- paste0(repoPath, 'output/')
sampleDataFile <- paste0(dataPath,'20190529_SampleData_Tumor_v6.txt')
patientTumorDataFile <- paste0(dataPath, '20190401_PatientTumorData_v5.txt')
patientOrder <- c('P340','P41','P49','P259','P302','P450','P303','P327','P375','P326','P453','P456','P300','P260','P372','P373','P413','P454','P455','P457','P452')
patientOrderSplit <- c('P340','P41','P49','P259','P302','P450','P303','P327','P375','P326','P453','P456','P300','P260-l','P260-m','P372','P373','P413','P454','P455','P457','P452')
subtypeColors <- list('IDH-mut_A'='#4DAF4A', 'IDH-mut_O'='#377EB8', 'IDH-wt'='#E41A1C')
