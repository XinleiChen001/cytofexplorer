###########################################################################
#                                                                         #
#                          版本信息      发布时间                         #
#                                                                         #
########################################################################### 

message("\n\nbackup version infor: cytofexplorer2.0\n")
backup_dir<-paste0(projectdir,"/","03_backups")
backup_scripts <- list.files(backup_dir,pattern='.R$', full=TRUE)

for (i in backup_scripts) {
  
  if(sub("^.*/","",i)!="Version_control.R")
    {source(i)
     cat("Loading ", sub("^.*/","",i),"\n")
    }
}
