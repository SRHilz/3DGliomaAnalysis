repoPath <- '/Users/shilz/Documents/Professional/Positions/UCSF_Costello/Publications/Hilz2018_IDHSpatioTemporal/Scripts/3DGliomaAnalysis/'
dataPath <- paste0(repoPath, 'data/')
outputPath <- paste0(repoPath, 'output/')
sampleDataFile <- paste0(dataPath,'20210207_SampleData_v12.txt')
patientTumorDataFile <- paste0(dataPath, '20210207_TumorData_v12.txt')
mutationDataFile <- file.path(outputPath, 'SID000014_mutational_categorization_clonal_subclonal_by_patient','SID000014_0_min_final_purity_samples_snv_vafs_categories.txt')
patientOrder <- c('P303','P327','P340','P375','P453','P470','P482','P490','P507','P508','P260','P300','P302','P450','P452','P372','P373','P413','P454','P455','P457','P458','P475','P485','P498','P499','P500','P501','P503','P276', 'P481')
patientOrderRec <- c('P340','P303','P327','P375','P453','P470','P490','P507','P508','P482','P452','P372','P373','P413','P454','P455','P457','P458','P475','P485','P498','P499','P500','P501','P503','P260','P300','P302','P450','P276','P481')
patientOrderOligo <- c('P340','P490','P302','P450','P303','P327','P375','P453','P470','P507','P508','P482','P260','P300','P452','P372','P373','P413','P454','P455','P457','P458','P475','P485','P498','P499','P500','P501','P503','P276', 'P481')
subtypeColors <- list('IDH-mut_A'='#4DAF4A', 'IDH-mut_O'='#377EB8', 'IDH-wt'='#E41A1C')
subsetIDHwtc <- c('P276','P372','P373','P413','P454','P455','P457','P458','P475','P481','P485','P498','P499','P500','P501','P503')
subsetIDHmut <- c('P260','P300','P302','P303','P327','P340','P375','P453','P470','P490','P507','P508','P482','P450')
subsetTP53mut <- c('P303','P327','P375','P453','P470','P507','P508','P260','P300','P413','P455','P475','P485','P499','P500')
subsetRecurrent <- c('P260','P300','P302','P450','P276','P481')
subset1p19q <- c('P340','P470','P490')
finalPurityCutoffCalls <- .2 # determined through exploration and used for all remaining analyses 
finalPurityCutoffNonzero <- .1 # determined through exploration and used for all remaining analyses 