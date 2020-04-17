repoPath <- '/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/'
dataPath <- paste0(repoPath, 'data/')
outputPath <- paste0(repoPath, 'output/')
sampleDataFile <- paste0(dataPath,'20200410_SampleData_Tumor_v11.txt')
patientTumorDataFile <- paste0(dataPath, '20200410_PatientTumorData_v11.txt')
patientOrder <- c('P303','P327','P340','P375','P453','P482','P260','P300','P302','P450','P452','P372','P373','P413','P454','P455','P457','P475','P485','P276', 'P481')
patientOrderRec <- c('P340','P303','P327','P375','P453','P482','P452','P372','P373','P413','P454','P455','P457','P475','P485','P260','P300','P302','P450','P276','P481')
patientOrderOligo <- c('P340','P302','P450','P303','P327','P375','P453','P482','P260','P300','P452','P372','P373','P413','P454','P455','P457','P475','P485','P276', 'P481')
subtypeColors <- list('IDH-mut_A'='#4DAF4A', 'IDH-mut_O'='#377EB8', 'IDH-wt'='#E41A1C')
