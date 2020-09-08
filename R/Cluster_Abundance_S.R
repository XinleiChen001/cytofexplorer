

#'@title 生成Cluster Aundance的火山图
#'@description
#'依据stat_by_cluster生成的list数据，生成火山图，图中颜色代表不同的Marker，数字和字母等符号代表不同的cluster；
#'结果采用了与SPADEvizR相同的形式，标注有“cite from SPADEVizR”或者“Adapted from SPADEVizR” ，使用时请注意在文章中引用原始文献。



#'@param cluster_stat       stat_by_cluster()函数生成的list数据
#'@param cond1              major_cond的一个，火山图（Volcano map）靠左一侧的条件，例如"B_Adjacent"
#'@param cond2              major_cond的一个，火山图（Volcano map）靠右一侧的条件，例如"A_Tumor"
#'@param xlimit             一个包含两个数字的向量（vector），例如c(-5,+5),用于手动指定火山图横坐标范围，默认值为NULL，此时会自动设置xlimit
#'@param ylimit             一个包含两个数字的向量（vector），例如c(0,3.5),用于手动指定火山图纵坐标范围，默认值为NULL，此时会自动设置xlimit
#'@param dif.level          一个数字，手动设置两组marker差异的的阈值大小，例如取默认值2时，两组均值相差2倍以上算显著（简单的说就是设置火山图上的两条竖线位置）
#'@param output_dir         输出数据文件夹名称
#'@param volcano.stat.paired  Volcano图两个condition是否是成对样本，默认为NULL，自动取值为stat.paired  
#'如使用全局统计默认参数，以下参数无需设置
#'
#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著），调整该数值可以调整火山图上的横线位置
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#' @return                none
#' @export


draw_abundance_volcano2<-function(cluster_stat,
                                  abundance_metadata,
                                  cond1,
                                  cond2,
                                  xlimit=NULL,
                                  ylimit=NULL ,
                                  dif.level=2,
                                  output_dir=paste0("Aundance_volcano_plots2_",cond1,"_vs_",cond2),
                                  volcano.stat.paired=NULL,
                                  #如使用全局统计默认参数，以下参数无需设置
                                  stat.paired=NULL,
                                  stat.method=NULL,
                                  set.equal.var=NULL,
                                  conf.level=NULL,
                                  p.adjust.method=NULL
                                ){


  
  
  if(0){
    xlimit=NULL
    ylimit=NULL
    dif.level=2
    output_dir=paste0("Aundance_volcano_plots_",cond1,"_vs_",cond2)
    #如使用全局统计默认参数，以下参数无需设置
    stat.paired=NULL
    stat.method=NULL
    set.equal.var=NULL
    conf.level=NULL
    p.adjust.method=NULL
  }
  


  p_para="log_p_value"
  adjust_label=""


  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  wdir=getwd()

  cat(paste0("Output to folder: ",wdir,"/",output_dir,"\n"))

  cluster_percent_summary<-cluster_stat[["cluster_percent_summary"]]
  major_cond             <- as.character(cluster_stat[["metadata"]]$major_cond)



  if(is.null(stat.paired)) stat.paired        <- cluster_stat[["metadata"]]$stat.paired
  if(is.null(stat.method)) stat.method        <- as.character(cluster_stat[["metadata"]]$stat.method)
  if(is.null(conf.level)) conf.level          <- cluster_stat[["metadata"]]$conf.level
  if(is.null(set.equal.var)) set.equal.var    <- as.character(cluster_stat[["metadata"]]$set.equal.var)
  if(is.null(p.adjust.method)) p.adjust.method    <- as.character(cluster_stat[["metadata"]]$p.adjust.method)
  merged_cluster_id <-cluster_stat[["merged_cluster_id"]]
  
  if(is.null(volcano.stat.paired)) volcano.stat.paired<-  stat.paired
  

  if(volcano.stat.paired){
    if_paired<-"paired"}else{
      if_paired<-"unpaired"}

  
  
  #cluster_percent_data<-full_join(cluster_percent_summary,groups,by="File_ID")
  cluster_percent_data<-abundance_metadata[["cluster_percent_data"]]
  cluster_percent_ctrl<-abundance_metadata[["cluster_percent_ctrl"]]
  cluster_names<-intersect(paste0("cluster",merged_cluster_id),colnames(cluster_percent_data))
  
  #计算各个cluster比例的变化
  cluster_percent_data_dif<-cluster_percent_data
  cluster_percent_data_dif[,cluster_names]<-cluster_percent_data_dif[,cluster_names]-cluster_percent_ctrl[,cluster_names]
  cluster_percent_data_difference<-cluster_percent_data
  cluster_percent_data_log_ratio<-cluster_percent_data
  No_zero<-0.01
  cluster_percent_data_log_ratio[,cluster_names]<-log((cluster_percent_data_log_ratio[,cluster_names]+No_zero)/(cluster_percent_ctrl[,cluster_names]+No_zero),base=2)
  #cluster_percent_data<-cluster_percent_data_dif
  
  
  
  
  #火山图
  draw_volcanos<-function(cluster_percent_data,stat_type_tag)
                    {
        
      #计算Cluster Ratio
    
      cond1_cluster<-
        cluster_percent_data %>%
        #filter_(paste0(major_cond,"==cond1"))
           dplyr::filter_at(vars(matches(major_cond)),all_vars(.==cond1)) 
    
      cond2_cluster<-
        cluster_percent_data %>%
        #filter_(paste0(major_cond,"==cond2"))
           dplyr::filter_at(vars(matches(major_cond)),all_vars(.==cond2)) 
        
    
      volcano_stat_result<-combined_stat(cond1_data=cond1_cluster[,colnames(cluster_percent_summary)[-1]],
                                         cond2_data=cond2_cluster[,colnames(cluster_percent_summary)[-1]],
                                         alternative="two.sided",
                                         stat.paired=volcano.stat.paired,
                                         conf.level=conf.level,
                                         stat.method=stat.method,
                                         p.adjust.method=p.adjust.method,
                                         set.equal.var=set.equal.var,
                                         silent=T)
    
      #summarise_each(funs(use_method))
      cond1_data_m<-
        cond1_cluster %>%
        summarise_if(is.numeric,mean,na.rm=F)
      cond2_data_m<-
        cond2_cluster %>%
        summarise_if(is.numeric,mean,na.rm=F)
      #summarise_each(funs(use_method))
    
      No_zero<-0.01  #add to prevent divided by zero
    
      #difference=log2((cond2_data_m+No_zero)/(cond1_data_m+No_zero))
      difference<-cond2_data_m-cond1_data_m
    
      cluster_volcano_data<-data.frame(cluster=paste0("Cluster",c(1:nrow(volcano_stat_result))),
                                       difference=t(difference),
                                       volcano_stat_result)
    
      cluster_percent_mean<-apply(cluster_percent_summary[,-1],2,mean)
    
      #为了方便绘图，所有小于0.0001的p值等于0.0001
      cluster_volcano_data$p_values[cluster_volcano_data$p_values<0.0001]=0.0001
      cluster_volcano_data$p_adjust[cluster_volcano_data$p_adjust<0.0001]=0.0001
      cluster_volcano_data$log_p_value=-log10(cluster_volcano_data$p_value)
      cluster_volcano_data$log_p_adjust=-log10(cluster_volcano_data$p_adjust)
      cluster_volcano_data$cluster_percent_mean = cluster_percent_mean
      cluster_volcano_data$Significant<-ordered(cluster_volcano_data$p_values<1-conf.level & abs(cluster_volcano_data$difference)>=log2(dif.level),levels=c(FALSE,TRUE))
      cluster_volcano_data$cluster_num<-sub("cluster","",row.names(cluster_volcano_data))
    
      x_threshold<-round(log2(dif.level),digits=1)
    
      if(is.null(xlimit)){
        xuplimit=max(abs(difference)+2)
        xlimit=c(-1*xuplimit,xuplimit)}
    
      if(is.null(ylimit)){
        yuplimit=min(max(cluster_volcano_data$log_p_value+1.5),4)
        ylimit=c(0,yuplimit)}
    
      mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                       legend.key = element_rect(fill = "white", colour = "white"), #图标
                       legend.background = (element_rect(colour= "white", fill = "white")))
    
    
      cluster_volcano<-ggplot(cluster_volcano_data)+
        geom_point(aes_string(x="difference",y=p_para,size="cluster_percent_mean",color="Significant"))+
        labs(y="-Log10(p_values)",title="Response of Cluster Abundance")+
        labs(x=paste0("Difference","\n", cond1, " <- enriched -> ", cond2,"\n\n Cite from SPADEVizR"))+
        labs(subtitle=paste0(adjust_label,"p values of ",stat_type_tag," were calculated by ",if_paired," ",stat.method))+
        mytheme+
        expand_limits(x=xlimit,y=ylimit)+
        labs(size= "Percentage %")+
        geom_vline(xintercept = c(x_threshold,-1*x_threshold), linetype="dashed",color = "blue", size=0.55)+
        geom_hline(yintercept=-log10(1-conf.level), linetype="dashed", color = "blue",size=0.55)+
        scale_x_continuous(breaks=c(seq(-1*ceiling(xuplimit),ceiling(xuplimit),by=1),x_threshold,-1*x_threshold),minor_breaks = NULL)+
        scale_y_continuous(breaks=c(seq(0,6,by=1),round(-log10(1-conf.level),digits=1)),minor_breaks = NULL)+
        theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
        theme(panel.background=element_rect(colour='black',size = 1 ))+
        theme(axis.line=element_line(colour='black',linetype = 1,size=0.3))+
        scale_color_manual(values =colorRampPalette(c("grey","red"))(2))+
        geom_text_repel(aes_string(x="difference",y=p_para),label=cluster_volcano_data$cluster_num,box.padding=0.25,point.padding=0.5)
    
      #准备FDR-pValue画图
      FDR_ranks<-factor(c("other(adjust_p>=10%)","adjust_p<10%","adjust_p<5%"),levels=c("other(adjust_p>=10%)","adjust_p<10%","adjust_p<5%"),ordered=TRUE)
      FDR_pValue<-dplyr::arrange(cluster_volcano_data,desc(p_values))
      FDR_pValue$FDR_rank<-FDR_ranks[1]
      FDR_pValue$FDR_rank[FDR_pValue$p_adjust<0.1]=FDR_ranks[2]
      FDR_pValue$FDR_rank[FDR_pValue$p_adjust<0.05]=FDR_ranks[3]
    
      FDR_table_data<-FDR_pValue[(FDR_pValue$p_values<0.05)&(!is.na(FDR_pValue$p_values)),c(1:6,9,10,12)]
    
      round_3<-function (var){return(round(var,3))}
    
      if(nrow(FDR_table_data)>0)
        FDR_table_data[,c(2,3,4)]<-apply(FDR_table_data[,c(2,3,4)],2,round_3)
      
        suppressWarnings(FDR_table<-ggtable(FDR_table_data))
    
      cat(paste0("Output Statistical Report(arranged).csv...\n"))
    
      FDR_pValue_arranged<-dplyr::arrange(FDR_pValue,p_values)
      write.csv(FDR_pValue_arranged,paste0("./",output_dir,"/",stat_type_tag," Volcano.csv"),row.names = TRUE)
    
    
      FDR_pValue$cluster_num<-factor(FDR_pValue$cluster_num,levels=FDR_pValue$cluster_num,ordered=TRUE)
    
      FDR_pValue_plot<-ggplot(FDR_pValue)+
        geom_bar(aes_string(x="cluster_num",y="log_p_value",fill="FDR_rank"),stat="identity",width = 0.8)+
        labs(y="-Log10(p_values)",title="Summary of Adjust p Value")+
        labs(x=paste0("Clusters"))+
        coord_flip()+
        labs(subtitle=paste0("Adjust p values were calculated based on ",p.adjust.method," adjust method"))+
        mytheme+
        expand_limits(y=c(0,2))+
        theme(panel.background=element_rect(colour='white',size = 1 ))+
        theme(axis.line.x=element_line(colour='black',linetype = NULL,size=0.3))+
        theme(axis.ticks.y=element_blank())+
        scale_fill_manual(values = c("lightgreen","navy","red"))+
        scale_y_continuous(breaks=c(seq(0,3,by=0.5)),minor_breaks = NULL)+
        theme(legend.position = c(0.7,0.4),legend.title=element_blank())+
        #theme(axis.text.y=element_blank())+
        theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))
    
      return_figs<-list(cluster_volcano,FDR_pValue_plot)
      return(return_figs)
      
      }
  # 
  # return_figs<- draw_volcanos(cluster_percent_data,stat_type_tag="cluster percentage")
  # cat(paste0("Output Cluster Abundance Volcano plot...\n"))
  # pdf(file=paste0("./",output_dir,"/cluster percentage",cond1," enriched ",cond2,".pdf"), width=15, height=10)
  # #pdf(file=paste0("./",output_dir,"/1.pdf"), width=15, height=10)
  # multiplot(return_figs[[1]],return_figs[[2]],cols = 2)
  # #multiplot(FDR_table)
  # dev.off()
  # 
  return_figs<- draw_volcanos(cluster_percent_data_dif,stat_type_tag="Percentage Change")
  cat(paste0("Output Cluster Abundance Volcano plot...\n"))
  pdf(file=paste0("./",output_dir,"/Percentage Change_",cond1," enriched ",cond2,".pdf"), width=15, height=10)
  multiplot(return_figs[[1]],return_figs[[2]],cols = 2)
  #multiplot(FDR_table)
  dev.off()
  
  return_figs<- draw_volcanos(cluster_percent_data_log_ratio,stat_type_tag="Percentage log2 ratio")
  cat(paste0("Output Cluster Abundance Volcano plot...\n"))
  pdf(file=paste0("./",output_dir,"/Percentage log2 ratio_",cond1," enriched ",cond2,".pdf"), width=15, height=10)
  multiplot(return_figs[[1]],return_figs[[2]],cols = 2)
  #multiplot(FDR_table)
  dev.off()
  
}




