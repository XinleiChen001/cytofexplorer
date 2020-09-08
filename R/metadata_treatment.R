
#============function======
#'@para rawname      fcs文件中原始名称
#'@para markernames  combind_data_raw的名称
#'@para add_cols     要在all_markers.csv中加入的列名 
#'  
#'@export


markers2csv<-function(rawname,markername,add_col=c("transform","tSNE","PhenoGraph","tsne_heatmap","expr_para","heatmap")){
  if(!file.exists(paste0(metadata_dir,"/all_markers.csv"))){

    nmarker<-length(markername)
    ncol<-length(add_col)
    
    all_markers_input<-matrix(data=rep("",nmarker*ncol),
                              ncol=ncol,
                              dimnames = list(markername,add_col))
    
    all_markers_input<-data.frame(rawname=c(rawname,"File_ID"),
                                  markers=markername,
                                  all_markers_input)
    write.csv(all_markers_input,paste0(metadata_dir,"/all_markers.csv"),row.names = FALSE)
    cat("all_markers.csv has been output successfully\n")
  }else{
    message("Warning: all_markers.csv is already exist ,stop outputing to avoid overwriting.\n")
  }


}
#=======function======




#============function======
#'@para markernames  combind_data_raw的名称
#'@para add_cols     要在all_markers.csv中加入的列名 
#'  
#'@export


groups2csv<-function(filename,add_col=c("Short_name")){
  if(!file.exists(paste0(metadata_dir,"/all_samples.csv"))){
    
    nfile<-length(filename)
    ncol<-length(add_col)

    all_samples_input<-matrix(data=rep("",nfile*ncol),
                              ncol=ncol,
                              dimnames = list(filename,add_col))
    
    all_samples_input<-data.frame(File_ID=filename,
                                  all_samples_input)
    write.csv(all_samples_input,paste0(metadata_dir,"/all_samples.csv"),row.names = FALSE)
    cat("all_samples.csv has been output successfully\n")
  }else{
    message("Warning: all_samples.csv is already exist ,stop outputing to avoid overwriting.\n")
  }
  
}
