



export_cytofkit_RData <-function(expressionData = combined_data_raw,
                                 dimReducedRes = dimReducedRes,
                                 clusterRes = clusterRes, 
                                 dimReductionMethod = names(dimReducedRes),
                                 visualizationMethods =names(dimReducedRes), #visualizationMethods,
                                 progressionRes = NULL,
                                 projectName = "projectName",
                                 rawFCSdir = raw_fcs_dir,
                                 resultDir = raw_fcs_dir,
                                 RDataName = NULL,
                                 dimRedMarkers = colnames(tSNE_input_data),
                                 sampleNames = unique(File_ID)){
  
  
  if (is.null(RDataName)) RDataName = paste0(projectName,".RData")
  
  for (i in c(1:length(dimReducedRes))){
    row.names(dimReducedRes[[i]])<-row.names(expressionData)}
  for (i in c(1:length(clusterRes))){
    names(clusterRes[[i]])<-row.names(expressionData)}
  
  ## wrap the results
  
  
  analysis_results <- list(expressionData = expressionData,
                           dimReductionMethod = names(dimReducedRes),
                           visualizationMethods =names(dimReducedRes), #visualizationMethods,
                           dimReducedRes = dimReducedRes,
                           clusterRes = clusterRes, 
                           progressionRes = NULL,
                           projectName = projectName,
                           rawFCSdir = rawFCSdir,
                           resultDir = resultDir,
                           dimRedMarkers = dimRedMarkers,
                           sampleNames = sampleNames)
  
  save(analysis_results, file = RDataName)
  
  cat(RDataName,paste0(" has been generated."))
}