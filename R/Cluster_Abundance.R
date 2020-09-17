
#2019_11_4  在abundance_boxplot里面加入groups_to_show 参数
#2019_11_5  在abundance_boxplot debug
#2019_11_22 将abundance_boxplot升级为abundance_report
#2019_11_26 修复abundance顺序bug
#2019_12_5  修复abundance非成对样本的bug
#2019_12_19 在abundance_report中增加Rowv=T,Colv=F,dendrogram="row",几个heatmap 聚类参数
#2020_01_04 cluster abundance histogram 的cluster模式
#2020_01_20 cluster abundance histogram 自动对齐聚类树枝图
#2020_01_30 对cluster abundance hide_ctrl参数进行扩展（可以同时导出带有ctrl和不带有ctrl的figure）
#2020_01_30 修复非成对样本计算control的一处 bug
#2020_09_04 matches-->one_of
#2020-9-17  ncolor<-nrow(unique(Plot_Data[,color_cond,drop=F]))加上drop=F




#'说明:
#'
#'SPADEViZR工具包可以对SPADE聚类结果进行可视化的统计分析，发现在组成比例上显著性差异的亚群。详细细节请参照其github主页
#'链接：https://github.com/tchitchek-lab/SPADEVizR
#'文献：SPADEVizR: an R package for visualization, analysis and integration of SPADE results. Gautreau G etc.
#'在本package中我们保持了与其分析结果的风格，并加以扩展，使其同时可以适用于与SPADE以外的分析结果的分析，并可以同时输出padjust前后的结果。所有风格与
#'SPADEvizR相近的图均标注有“cite from SPADEVizR”或者“Adapted from SPADEVizR” ，使用时请注意在文章中引用原始文献。
#'

#'Notice:
#'SPADEViZR R tool package provides a visual analysis tool of the SPADE clustering results and finds population that are significantly different in composition. Please refer to its github homepage for details.
#'Link: https://github.com/tchitchek-lab/SPADEVizR
#'Published Article:SPADEVizR: an R package for visualization, analysis and integration of SPADE results. Gautreau G etc.
#'In this package, we adopt the similar style of its analysis results and extend it so that it can be applied to the clustering analysis other than SPADE,
#'and can output the results before and after p-adjust. All figures similar to SPADEvizR are marked with "cite from SPADEVizR" or "Adapted from SPADEVizR".
#'Please pay attention to it, and cite the original article in your reference.
#'




#' @title cluster统计中细胞数量的质控（quality control）
#'
#' @description
#' 对cluster进行统计分析时，过少的细胞数有可能造成统计结果的不准确，本函数统计分析方法分析每个cluster的绝对数量是否符合要求
#' 结果采用了与SPADEvizR相近的风格，标注有“cite from SPADEVizR”或者“Adapted from SPADEVizR” ，使用时请注意在文章中引用原始文献。
#'
#'@param cluster_stat         stat_by_cluster()函数生成的list数据
#'@param count_threshold      数字，指定统计要求的细胞数，默认50
#'@param percent_threshold    数字，指定统计要求的细胞数，默认0.5，即0.5%
#'@param group_to_show        指定要统计的condition（major_cond中的一个或者多个），如不设置直接统计全部；
#'@param output_dir           输出数据文件夹名称
#'
#'如使用统计默认参数，以下参数无需设置
#'
#'@param stat.paired        逻辑变量，TRUE或者FALSE，表示是否是成对检验
#'@param conf.level         显著性水平，默认0.95,（即p<0.05为显著）
#'@param stat.method        一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#'@param set.equal.var      一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#'@param p.adjust.method    对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#'@return none
#'@export


cluster_significant_report<-function(cluster_stat,
                                   count_threshold=50,
                                   percent_threshold=0.5,
                                   group_to_show=NULL,
                                   output_dir="Aundance_significant_report",

                                   #如使用全局统计默认参数，以下参数无需设置
                                   stat.paired=NULL,
                                   stat.method=NULL,
                                   set.equal.var=NULL,
                                   conf.level=NULL,
                                   p.adjust.method=NULL){


  
  if(0){
  group_to_show=NULL
  output_dir="Aundance_significant_report"
  
  #如使用全局统计默认参数，以下参数无需设置
  stat.paired=NULL
  stat.method=NULL
  set.equal.var=NULL
  conf.level=NULL
  p.adjust.method=NULL
  
}
  if(!is.null(p.adjust.method)){

  if(p.adjust.method=="none"){
     adjust_label=""}else{
      adjust_label="adjust " }
    }


  clustern          <- cluster_stat[["metadata"]]$clustern
  groups            <-cluster_stat[["groups"]]
  major_cond       <-as.character(cluster_stat[["metadata"]]$major_cond)
  group_names<-unique(groups[,major_cond])
  group_n<-length(group_names)

  group_list<-as.list(group_names)

  if(is.null(group_to_show)){
    group_to_show=c(1:group_n)
  }
  group_to_show_n<-length(group_to_show)
  if(is.numeric(group_to_show)){
  group_to_show<-group_names[group_to_show]}

  group_list<-as.list(group_to_show)
  if(group_to_show_n>1){
    group_list[[group_to_show_n+1]]<-group_names
   }
  significant_count_stats<-data.frame(cluster_num=as.character(merged_cluster_id))
  significant_percent_stats<-data.frame(cluster_num=as.character(merged_cluster_id))
  
  significant_count_stats$cluster_num<-as.character(  significant_count_stats$cluster_num)
  significant_percent_stats$cluster_num<-as.character(significant_percent_stats$cluster_num)


  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  wdir=getwd()

  cat(paste0("Output to folder: ",wdir,"/",output_dir,"\n"))


  pdf(file=paste0("./",output_dir,"/","Cluster Abundence(cell count and percentage).pdf"), width=15, height=10)

  for(list_i in c(length(group_list):1)){

    #list_i=3
    significant_plot<-cluster_significant_plot(cluster_stat,
                             count_threshold=count_threshold,
                             percent_threshold=percent_threshold,
                             group_to_show=group_list[[list_i]],

                             #如使用全局统计默认参数，以下参数无需设置
                             stat.paired=stat.paired,
                             stat.method=stat.method,
                             set.equal.var=set.equal.var,
                             conf.level=conf.level,
                             p.adjust.method=p.adjust.method)

   multiplot(significant_plot[[1]],significant_plot[[2]],cols = 2)  #layout=matrix(c(1,2),cols = 2)
   cluster_count_stat <- significant_plot[[3]]
   cluster_percent_stat <- significant_plot[[4]]

   significant_count_stat<-cluster_count_stat[c("significant","cluster_num")]
   
   significant_count_stats<-full_join(significant_count_stats,significant_count_stat,by="cluster_num")
   significant_count_stats_colname<-colnames(significant_count_stats)
   significant_count_stats_colname[ncol(significant_count_stats)]<-paste(group_list[[list_i]],collapse="_")
   colnames(significant_count_stats) <-significant_count_stats_colname

   significant_percent_stat<-cluster_percent_stat[c("significant","cluster_num")]
   significant_percent_stats<-full_join(significant_percent_stats,significant_percent_stat,by="cluster_num")
   significant_percent_stats_colname<-colnames(significant_percent_stats)
   significant_percent_stats_colname[ncol(significant_percent_stats)]<-paste(group_list[[list_i]],collapse="_")
   colnames(significant_percent_stats) <-significant_percent_stats_colname

   }




  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))

  count_melt<-melt(as.matrix(significant_count_stats[,-1]))
  count_melt$Significant<-as.factor(count_melt$value)
  count_summary_plot<-ggplot(count_melt, aes(y=Var1, x=Var2, fill=Significant))+
    xlab('Tissue_Type')+
    ylab("Clusters")+
    geom_tile(color="white", size=0.1)+
    scale_fill_manual(values =colorRampPalette(c("grey","red"))(2))+
    labs(title="Abundance Significant Summary\nCell Count")+
    scale_y_reverse(breaks=c(seq(1,clustern,by=1)),minor_breaks = NULL)+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
    #theme(panel.background=element_rect(colour='white',size = 1 ))+
    theme(axis.line=element_line(colour='white',linetype = 1,size=0.3))+
    #scale_y_reverse()+
    mytheme


  percent_melt<-melt(as.matrix(significant_percent_stats[,-1]))
  percent_melt$Significant<-as.factor(percent_melt$value)
  percent_summary_plot<-ggplot(percent_melt, aes(y=Var1, x=Var2, fill=Significant))+
    xlab('Tissue_Type')+
    ylab("Clusters")+
    geom_tile(color="white", size=0.1)+
    scale_fill_manual(values =colorRampPalette(c("grey","red"))(2))+
    labs(title="Abundance Significant Summary\nPercentage")+
    scale_y_reverse(breaks=c(seq(1,clustern,by=1)),minor_breaks = NULL)+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
    #theme(panel.background=element_rect(colour='white',size = 1 ))+
    theme(axis.line=element_line(colour='white',linetype = 1,size=0.3))+
    mytheme

  multiplot(count_summary_plot,percent_summary_plot,cols = 2)
  dev.off()
  cat("Fig export finished. ")

}





#' @title cluster统计中细胞数量的质控（quality control）
#'
#' @description
#' 对cluster进行统计分析时，过少的细胞数有可能造成统计结果的不准确，本函数统计分析方法分析每个cluster的绝对数量是否符合要求
#' 结果采用了与SPADEvizR相近的风格，标注有“cite from SPADEVizR”或者“Adapted from SPADEVizR” ，使用时请注意在文章中引用原始文献。

#'
#'@param cluster_stat         stat_by_cluster()函数生成的list数据
#'@param count_threshold      数字，指定统计要求的细胞数，默认50
#'@param percent_threshold    数字，指定统计要求的细胞数，默认0.5，即0.5%
#'@param group_to_show        指定要统计的condition（major_cond中的一个或者多个），如不设置直接统计全部；
#'@param output_dir           输出数据文件夹名称
#'
#'如使用统计默认参数，以下参数无需设置
#'
#'@param stat.paired        逻辑变量，TRUE或者FALSE，表示是否是成对检验
#'@param stat.method        一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#'@param set.equal.var      一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#'@param conf.level         显著性水平，默认0.95,（即p<0.05为显著）
#'@param p.adjust.method    对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#'@return none
#'@export

cluster_significant_plot<-function(cluster_stat,
                                     count_threshold=50,
                                     percent_threshold=0.5,
                                     group_to_show=NULL,
                                     output_dir="Aundance_significant_report",

                                     #如使用全局统计默认参数，以下参数无需设置
                                     stat.paired=NULL,
                                     stat.method=NULL,
                                     set.equal.var=NULL,
                                     conf.level=NULL,
                                     p.adjust.method=NULL
                                    ){



  clustern          <- cluster_stat[["metadata"]]$clustern
  cluster_cell_count<- cluster_stat[["cluster_cell_count"]]
  groups            <-cluster_stat[["groups"]]
  all_markers       <-cluster_stat[["all_markers"]]
  major_cond<-as.character(cluster_stat[["metadata"]]$major_cond)
  cluster_count_summary<-cluster_stat[["cluster_count_summary"]]
  cluster_percent_summary<-cluster_stat[["cluster_percent_summary"]]
  groups<-cluster_stat[["groups"]]
  

  

  if(is.null(stat.paired)) stat.paired        <- cluster_stat[["metadata"]]$stat.paired
  if(is.null(stat.method)) stat.method        <- as.character(cluster_stat[["metadata"]]$stat.method)
  if(is.null(conf.level)) conf.level          <- cluster_stat[["metadata"]]$conf.level
  if(is.null(set.equal.var)) set.equal.var    <- as.character(cluster_stat[["metadata"]]$set.equal.var)
  if(is.null(p.adjust.method)) p.adjust.method    <- as.character(cluster_stat[["metadata"]]$p.adjust.method)


  if(stat.paired){
    if_paired<-"paired"}else{
      if_paired<-"unpaired"}


  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }



  if(!is.null(p.adjust.method)){
  if(p.adjust.method=="none"){
    adjust_label=""}else{
      adjust_label="adjust " }
  }

  # if(is.numeric(group_to_show)){
  # select_groups<-paste0("\"",unique(groups[,major_cond])[group_to_show],"\"")}else{
  # select_groups<-paste0("\"",group_to_show,"\"")
  # }
  # 
  # select_groups<-paste(select_groups,collapse =",")
  # select_groups_short<-gsub(",","_",select_groups)
  # select_groups_short<-gsub("\"","",select_groups_short)
  # 
  # 
  
  if(is.numeric(group_to_show)){
    select_groups<-paste0(unique(groups[,major_cond])[group_to_show])}else{
      select_groups<-as.character(group_to_show)
    }
  
  select_groups2<-paste(select_groups,collapse =",")
  select_groups_short<-gsub(",","_",select_groups2)
  select_groups_short<-gsub("\"","",select_groups_short)


  
  cluster_count_summary<-cluster_stat[["cluster_count_summary"]]
  
  
  cluster_count_summary<-cluster_count_summary %>%
                           full_join(groups,by="File_ID") %>%
                               #  filter_(paste0(major_cond,"%in%c(",select_groups,")"))%>%
                              dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% c(select_groups))) %>%
                                 dplyr::select(num_range("cluster",1:clustern))

  #filter(groups,Tissue_Type%in%c("Biopsy","PBMC" ))

  #cluster cell count statisic
  cluster_count_mean<-apply(cluster_count_summary,2,mean)
  cluster_count_pvalue<-combined_stat(cond1_data=log(cluster_count_summary+1,base=10),
                                      cond2_data=NULL,
                                      alternative="greater",
                                      mu=log(count_threshold+1,base=10),
                                      stat.method=stat.method,
                                      stat.paired=F,
                                      conf.level=conf.level,
                                      p.adjust.method=p.adjust.method,
                                      set.equal.var=set.equal.var,
                                      silent=T)

  cluster_count_stat<-data.frame( cluster_count_pvalue,cluster_count_mean)
  cluster_count_stat$log_p_adjust<--log(cluster_count_stat$p_adjust,base=10)
  cluster_count_stat$log_count_mean<-log(cluster_count_stat$cluster_count_mean+1,base=10)
  cluster_count_stat$significant<-ordered(cluster_count_stat$p_adjust<1-conf.level & cluster_count_stat$cluster_count_mean>count_threshold,levels=c(TRUE,FALSE))
  cluster_count_stat$cluster_num<-as.character(sub("cluster","",row.names(cluster_count_stat)))
  #cluster_count_stat$cluster_count_mean[cluster_count_stat$cluster_count_mean==0]<-1



  #cluster cell percentage statisic

  cluster_percent_summary<-cluster_percent_summary %>%
                                full_join(groups,by="File_ID") %>%
                                       #filter_(paste0(major_cond,"%in%c(",select_groups,")"))%>%
                                       dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% c(select_groups))) %>%
                                          dplyr::select(num_range("cluster",1:clustern))



  cluster_percent_mean<-apply(cluster_percent_summary,2,mean)

  cluster_percent_pvalue<-combined_stat(cond1_data=cluster_percent_summary,
                                        cond2_data=NULL,
                                        alternative="greater",
                                        mu=percent_threshold,

                                        stat.method=stat.method,
                                        stat.paired=F,
                                        conf.level=conf.level,
                                        p.adjust.method=p.adjust.method,
                                        set.equal.var=set.equal.var,
                                        silent=T)


  cluster_percent_stat<-data.frame(cluster_percent_pvalue,cluster_percent_mean)
  cluster_percent_stat$log_pvalue<--log(cluster_percent_stat$p_values,base=10)
  cluster_percent_stat$log_p_adjust<--log(cluster_percent_stat$p_adjust,base=10)
  cluster_percent_stat$cluster_num<-as.character(sub("cluster","",row.names(cluster_percent_stat)))
  cluster_percent_stat$significant<-ordered(cluster_percent_stat$p_adjust<1-conf.level & cluster_percent_stat$cluster_percent_mean>percent_threshold,levels=c(TRUE,FALSE))

  
 
  #Output abundance statisitc result

  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))

  #cat(paste0("Outputing Cluster Relative Abundence(percentage).pdf\n"))

  #pdf(file=paste0("Cluster Abundence(cell count and percentage).pdf"), width=15, height=10)

  xlim_r=max(cluster_percent_stat$log_p_adjust)*1.2
  ylim_r=max(cluster_percent_stat$cluster_percent_mean*1.2)

  relative_abundance_plot<-ggplot(cluster_percent_stat)+
    geom_point(aes_string(x="log_p_adjust",y="cluster_percent_mean",size="cluster_percent_mean",colour="significant"))+
    scale_size("Percentage %",trans= "sqrt")+
    labs(x=paste0("-Log10(",adjust_label,"_p_values)\n\n Cite from SPADEVizR "),y="Average percentage %",title=paste0("Cluster Relative Abundence(percentage)\nGroups: ",select_groups))+
    labs(subtitle=paste0(adjust_label,"p values of cluster percentage were calculated by ",if_paired," ",stat.method))+
    mytheme+
    expand_limits(x=c(0,xlim_r),y=c(0,ylim_r))+
    geom_vline(xintercept = -log((1-conf.level),base = 10), linetype="dashed",color = "blue", size=0.8)+
    geom_hline(yintercept=percent_threshold, linetype="dashed", color = "blue",size=0.8)+
    theme(panel.grid.major=element_line(colour='grey80',linetype = 1,size=0.3))+
    scale_x_continuous(breaks=c(seq(0,xlim_r,by=2),round(-log(1-conf.level,base=10),digits=1)),minor_breaks = NULL)+
    scale_y_continuous(breaks=c(seq(0,ylim_r,by=5),percent_threshold),minor_breaks = NULL)+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
    theme(panel.background=element_rect(colour='black',size = 1 ))+
    theme(axis.line=element_line(colour='black',linetype = 1,size=0.3))+
    scale_color_manual(values =colorRampPalette(c("red","grey"))(2))+
    geom_text_repel(aes_string(x="log_p_adjust",y="cluster_percent_mean"),label=cluster_percent_stat$cluster_num,box.padding=0.25,point.padding=0.5)

   #输出p值
   write.csv(cluster_percent_stat,paste0("./",output_dir,"/",select_groups_short,"_percent_stat.csv"),row.names = FALSE)

  #ggtable(cluster_percent_stat[,c(1,2,3,4,5,8)])


  xlim_a<-max(cluster_count_stat$log_p_adjust)*1.2
  ylim_a=max(cluster_count_stat$cluster_count_mean*1.5)
  absolute_abundance_plot<-ggplot(cluster_count_stat)+
    geom_point(aes_string(x="log_p_adjust",y="cluster_count_mean",size="cluster_count_mean",colour="significant"))+
    scale_size("Cell count",trans= "log10",limits = c(1,ylim_a))+
    labs(x=paste0("-Log10(",adjust_label,"_p_values)\n\n Cite from SPADEVizR "),y="Average Cell Count",title=paste0("Cluster Absolute Cell Count\nGroups: ",select_groups))+
    labs(subtitle=paste0(adjust_label,"p values of log tranformed cluster cell count were calculated by ",if_paired," ",stat.method))+
    mytheme+
    expand_limits(x=c(0,xlim_a),y=c(1,ylim_a))+
    geom_vline(xintercept = -log((1-conf.level),base = 10), linetype="dashed",color = "blue", size=0.8)+
    geom_hline(yintercept=count_threshold, linetype="dashed", color = "blue",size=0.8)+
    theme(panel.grid.major=element_line(colour='grey80',linetype = 1,size=0.3))+
    scale_x_continuous(breaks=c(seq(0,xlim_a,by=2),round(-log(1-conf.level,base=10),digits=1)),minor_breaks = NULL)+
    scale_y_log10(breaks=c(10^seq(0,ceiling(log10(ylim_a)),by=1),count_threshold),minor_breaks = NULL)+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
    theme(panel.background=element_rect(colour='black',size = 1 ))+
    theme(axis.line=element_line(colour='black',linetype = 1,size=0.3))+
    scale_color_manual(values =colorRampPalette(c("red","grey"))(2))+
    geom_text_repel(aes_string(x="log_p_adjust",y="cluster_count_mean"),label=cluster_count_stat$cluster_num,box.padding=0.25,point.padding=0.5)

  write.csv(cluster_percent_stat,paste0("./",output_dir,"/",select_groups_short,"_count_stat.csv"),row.names = FALSE)

  significant_plot<-list()
  significant_plot[[1]]<-absolute_abundance_plot
  significant_plot[[2]]<-relative_abundance_plot
  significant_plot[[3]]<-cluster_count_stat
  significant_plot[[4]]<-cluster_percent_stat

  return(significant_plot)


}




#' @title 画各个Cluster Abundance的boxplot，Heatmap
#'
#' @description
#' 对cluster进行统计分析时，过少的细胞数有可能造成统计结果的不准确，本函数统计分析方法分析每个cluster的绝对数量是否符合要求
#'
#'@param cluster_stat  stat_by_cluster()函数生成的list数据
#'@param cluster_id    指定输出cluster的id，例如：如果要输出前五个cluster，cluster_id=c(1:5); 如果要输出全部cluster，保持其默认值NULL即可。
#'@param color_cond    指定依据哪个条件进行着色
#'@param subgroup_cond     heatmap_ctrl对应的组所在的列名
#'@param group_colorset    组别的配色方案
#'@param output_dir        输出数据文件夹名称，默认"Abundance_boxplot"
#'@param heatmap_ctrl  指定heatmap、boxplot的control组名称，例如“PBMC”   

#'@param heatmap1_color    指定Cluster Percentage heatmap配色方案
#'@param heatmap2_color    指定1.2_Change of Cluster Percentage和1.3_log2_Percentage ratio heatmap配色方案
#'@param Rowv,Colv         逻辑变量，分别heatmap设定行和列是否聚类
#'@param dendrogram        显示heatmap行或者列的树形图,"both","row","column","none"
#'@param hide_ctrl         在图片中隐藏对照组（仅用于多个ctrl的情况）

#boxplot parameters
#'@param boxplot_line_width       boxplot线宽
#'@param comparisons              list变量，用来指定需要显示p值的组别，例如list(c("PBMC","Biopsy")，c("PBMC","Tumor"))就是要分别显示PBMC与Biopsy和Tumor两组之间的p值；如果comparisons=NULL，则显示各组与对照组的p值
#'@param comparisons.stat.paired  指定组别是否为成对样本   
#'@param comparisons.stat.method  指定组别的p值计算方法
#'@param show_pvalues  逻辑变量，TRUE或者FALSE，决定是否在图上显示p Value


#barplot参数
#'@param barplot_clustered           逻辑变量，决定丰度直方图是否依照各个cluster比例聚类，默认为TRUE
#'@param barplot_direction           取值为"v"或者"h",分别决定直方图的方向为竖直或者水平
#'@param barplot_annotation          逻辑变量，是否注释组别，默认为TRUE
#'@param cluster_colorset             barplot中cluster的颜色方案



#'如使用统计默认参数，以下参数无需设置
#'@param group_seq          设置各组顺序，格式实例： group_seq=c("PBMC","Biopsy")，注：括号内为各组名称；如不设置，默认与总体参数保持一致
#'@param groups_to_show     想要展示部分组时，指定展示组的组名，默认为NULL，展示所有组
#'@param stat.paired        逻辑变量，TRUE或者FALSE，表示是否是成对检验(如有多个control，则只subgroup之内是否是成对样本)
#'@param stat.method        一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#'@param set.equal.var      一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#'@param conf.level         显著性水平，默认0.95,（即p<0.05为显著）
#'@param p.adjust.method    对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#'
#'
#'@return a list containing ggplot2 objects  
#'@export


draw_abundance_report<-function(cluster_stat,
                                cluster_id=NULL,
                                color_cond = NULL,
                                subgroup_cond=NULL,
                                group_colorset =brewer_color_sets,
                                output_dir=paste0("Abundance_Report"),
                                
                                #heatmap parameters
                                heatmap_ctrl,
                                heatmap1_color =colorRampPalette(c("black","yellow")),
                                heatmap2_color =colorRampPalette(c("blue","white","red")),
                                Rowv=T,Colv=F,dendrogram="row",
                                hide_ctrl=NULL,
                                
                               #boxplot parameters
                               boxplot_line_width=1,
                               comparisons=NULL,
                               comparisons.stat.paired=NULL,
                               comparisons.stat.method=NULL,
                               show_pvalues=TRUE,
                              

                               #barplot参数
                               barplot_clustered=TRUE,
                               barplot_direction="v",
                               barplot_annotation=TRUE,
                               cluster_colorset =dif_seq_rainbow,
                              
                               
                               #如使用全局统计默认参数，以下参数无需设置
                               group_seq=NULL,
                               groups_to_show=NULL,
                               stat.paired=NULL,
                               stat.method=NULL,
                               set.equal.var=NULL,
                               conf.level=NULL,
                               p.adjust.method=NULL,
                              
                               #shinnyapp parameters
                               use_shiny_para=FALSE,
                               outputfig =TRUE){
  
  
  
  if(0){
    cluster_id=NULL
    show_pvalues=TRUE
    use_shiny_para=FALSE
    outputfig =TRUE
    group_colorset =brewer_color_sets
    color_cond=NULL
    heatmap_ctrl=NULL
    boxplot_line_width=1
    comparisons.stat.paired=NULL
    comparisons.stat.method=NULL
    subgroup_cond=NULL
    comparisons=NULL
    output_dir=paste0("Abundance_Report")
    heatmap1_color =colorRampPalette(c("black","yellow"))
    heatmap2_color =colorRampPalette(c("blue","white","red"))
    Rowv=T
    Colv=F
    dendrogram="row"
    barplot_clustered=TRUE
    barplot_direction="h"
    #如使用全局统计默认参数，以下参数无需设置
    group_seq=NULL
    groups_to_show=NULL
    stat.paired=NULL
    stat.method=NULL
    set.equal.var=NULL
    conf.level=NULL
    p.adjust.method=NULL
  }
  
  
  cat(paste0("Outputing Abundence Report\n"))
  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }
  if (!dir.exists(paste0("./",output_dir,"/csvs"))) {
    dir.create(paste0("./",output_dir,"/csvs"))
  }
  

  col_color_index  <- cluster_stat[["col_color_index"]]
  merged_cluster_id<-cluster_stat[["merged_cluster_id"]]
  groups           <- cluster_stat[["groups"]]
  major_cond       <-as.character(cluster_stat[["metadata"]]$major_cond)
  #heatmap_ctrl            <- as.character(cluster_stat[["metadata"]]$heatmap_ctrl)
  cluster_count_summary   <-cluster_stat[["cluster_count_summary"]]
  cluster_percent_summary <-cluster_stat[["cluster_percent_summary"]]
  clustern <-cluster_stat[["metadata"]]$clustern
  
  if(is.null(color_cond)) color_cond<-major_cond
  
  if(is.null(groups_to_show)){
    groups_to_show=unique(groups[,major_cond,drop=T])  
  }   
  
  if (!all(heatmap_ctrl %in% groups_to_show))
  { message("heatmap_ctrl should be included in groups_to_show\n")
    return(NULL)
  }  
  
  

  if (use_shiny_para==FALSE){

    if(is.null(stat.paired)) stat.paired        <- cluster_stat[["metadata"]]$stat.paired
    if(is.null(stat.method)) stat.method        <- as.character(cluster_stat[["metadata"]]$stat.method)
    if(is.null(conf.level)) conf.level          <- cluster_stat[["metadata"]]$conf.level
    if(is.null(set.equal.var)) set.equal.var    <- as.character(cluster_stat[["metadata"]]$set.equal.var)
    if(is.null(p.adjust.method)) p.adjust.method    <- as.character(cluster_stat[["metadata"]]$p.adjust.method) 


    if (!is.null(group_seq)){
      
      major_cond_groups<-unique(groups[,major_cond])
      
      if (length(major_cond_groups)!=length(group_seq)){
      
      group_seq<-c(group_seq,major_cond_groups[!(major_cond_groups %in% group_seq)])}
      
      groups[,major_cond]<-factor(groups[,major_cond,drop=T],
                                 levels=group_seq,
                                 ordered=T)
    }else { group_seq   <- as.character(cluster_stat[["group_seq"]])}

    axis_font_size   =8
    label_font_size  =7
    title_font_size  =8
    
    group_table<-unique(groups[,major_cond])
    label_font_angle=0
    #label_font_angle=45
    
    # hjust =1
    # vjust = 0.5
    hjust =0.5
    vjust =0
    file_format="PDF"
    
  }else{
    shinny_para<-readRDS("shinny_para")
    eval(parse(text=paste0("colorset<-",shinny_para$colorset)))
    boxplot_ncol  =shinny_para$boxplot_ncol
    boxplot_width =shinny_para$boxplot_width
    boxplot_height=shinny_para$boxplot_height
    axis_font_size   =shinny_para$axis_font_size
    label_font_size  =shinny_para$label_font_size
    title_font_size  =shinny_para$title_font_size
    label_font_angle=shinny_para$label_font_angle
    hjust           =shinny_para$hjust
    vjust           =shinny_para$vjust
    line_width     =shinny_para$line_width
    show_pvalues  =shinny_para$show_pvalues
    file_format   =shinny_para$file_format
    export        =shinny_para$export
    stat.method   =shinny_para$stat.method
    stat.paired   =shinny_para$stat.paired
    conf.level    =shinny_para$conf.level
    p.adjust.method=shinny_para$p.adjust.method
    set.equal.var     =shinny_para$set.equal.var
    group_table   =shinny_para$group_table
    cluster_table =shinny_para$cluster_table
    cluster_id<-as.numeric(cluster_table)
  }


  if(is.null(cluster_id[1])){
    cluster_id=merged_cluster_id
  }
  
  cluster_names<-paste0("cluster",cluster_id)
  cluster_percent_data<-cluster_percent_summary%>%
                         dplyr::full_join(groups,by="File_ID")%>%
                         dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% groups_to_show))%>%
                         arrange_at("Short_name") %>%
                         arrange_at(major_cond)
  

 
  # show_group_id=which(cluster_percent_data[,major_cond,drop=T] %in% group_table)
  # cluster_percent_data<-cluster_percent_data[show_group_id,]

  
  if(stat.method=="auto") {
    stat.method="t.test"
    cat("Warning: Outputing fiures in auto mode, set stat.method to t-test\n" )

  }else { stat.method=sub("-",".",stat.method)}
  if(is.null(comparisons.stat.paired)) comparisons.stat.paired<-stat.paired
  if(is.null(comparisons.stat.method)) {comparisons.stat.method<-stat.method}else{comparisons.stat.method=sub("-",".",comparisons.stat.method)}
  

  #cluster差值
  if (is.null(subgroup_cond))  {
    cluster_percent_data$whole<-"Yes"
    subgroup_cond<-"whole"
  } else if(subgroup_cond=="whole"){
    cluster_percent_data$whole<-"Yes"
    }
  
  #unique(cluster_percent_data$CompareB4_to_9)
  #计算control matrix

  cluster_percent_ctrl<-cluster_percent_data
# ==============
#   if(stat.paired){
#     for(this_subgroup_cond in unique(cluster_percent_data[,subgroup_cond,drop=T])){
#       #this_subgroup_cond=unique(cluster_percent_data[,subgroup_cond,drop=T])[1]
#       
#       
#      single_subgroup_data<-dplyr::filter_at(cluster_percent_data,vars(one_of(subgroup_cond)),all_vars(.==this_subgroup_cond))
#      single_abundance_ctrl<-dplyr::filter_at(single_subgroup_data,vars(one_of(major_cond)),all_vars(.%in% heatmap_ctrl))
#          for(this_major_cond in unique(single_subgroup_data[,major_cond,drop=T])){
#            #this_major_cond<-unique(single_subgroup_data[,major_cond,drop=T])[2]
#            this_major_cond_id<-which((cluster_percent_ctrl[,major_cond,drop=T])==this_major_cond)
#            cluster_percent_ctrl[this_major_cond_id,cluster_names]<-single_abundance_ctrl[,cluster_names]
#          }
#      #cluster_percent_ctrl<-rbind.data.frame(cluster_percent_ctrl,single_subgroup_data)
#      }  
#   }else{
#     for(this_subgroup_cond in unique(cluster_percent_data[,subgroup_cond,drop=T])){
#       #this_subgroup_cond=unique(cluster_percent_data[,subgroup_cond,drop=T])[1]
#       single_subgroup_data<-dplyr::filter_at(cluster_percent_data,vars(one_of(subgroup_cond)),all_vars(.==this_subgroup_cond))
#       single_abundance_ctrl<-dplyr::filter_at(single_subgroup_data,vars(one_of(major_cond)),all_vars(.%in% heatmap_ctrl))
#       single_abundance_ctrl_one_row<-apply(single_abundance_ctrl[,cluster_names],2,median)
#       
#       for (row_i in 1:nrow(single_abundance_ctrl)) {
#         single_abundance_ctrl[row_i,cluster_names]<-single_abundance_ctrl_one_row
#       }
#       for(this_major_cond in unique(single_subgroup_data[,major_cond,drop=T])){
#         #this_major_cond<-unique(single_subgroup_data[,major_cond,drop=T])[2]
#         this_major_cond_id<-which((cluster_percent_ctrl[,major_cond,drop=T])==this_major_cond)
#         cluster_percent_ctrl[this_major_cond_id,cluster_names]<-single_abundance_ctrl[,cluster_names]
#       }
#     }
# 
#    }
#====   

  
  for(this_subgroup_cond in unique(cluster_percent_data[,subgroup_cond,drop=T])){
                #this_subgroup_cond=unique(cluster_percent_data[,subgroup_cond,drop=T])[1]
                single_subgroup_data<-dplyr::filter_at(cluster_percent_data,vars(one_of(subgroup_cond)),all_vars(.==this_subgroup_cond))
                single_abundance_ctrl<-dplyr::filter_at(single_subgroup_data,vars(one_of(major_cond)),all_vars(.%in% heatmap_ctrl))
                
                
                if(stat.paired==F){
                          single_abundance_ctrl_one_row<-apply(single_abundance_ctrl[,cluster_names],2,median)
                          for (row_i in 1:nrow(single_abundance_ctrl)) {
                               single_abundance_ctrl[row_i,cluster_names]<-as.list(single_abundance_ctrl_one_row)#4.0
                           }
                        }
                for(this_major_cond in unique(single_subgroup_data[,major_cond,drop=T])){
                          #this_major_cond<-unique(single_subgroup_data[,major_cond,drop=T])[2]
                          this_major_cond_id<-which((cluster_percent_ctrl[,major_cond,drop=T])==this_major_cond)
                          cluster_percent_ctrl[this_major_cond_id,cluster_names]<-single_abundance_ctrl[,cluster_names]
                        }
     }

  
  
      #计算各个cluster比例的变化
      cluster_percent_data_dif<-cluster_percent_data
      cluster_percent_data_dif[,cluster_names]<-cluster_percent_data_dif[,cluster_names]-cluster_percent_ctrl[,cluster_names]
      cluster_percent_data_log_ratio<-cluster_percent_data
      No_zero<-0.01
      cluster_percent_data_log_ratio[,cluster_names]<-log((cluster_percent_data_log_ratio[,cluster_names]+No_zero)/(cluster_percent_ctrl[,cluster_names]+No_zero),base=2)

      
abundance_report<-function(hide_ctrl){      
        
      file_tag=""    
      if(hide_ctrl==T){
             cluster_percent_data<-dplyr::filter_at(cluster_percent_data,vars(one_of(major_cond)),all_vars(!(. %in% heatmap_ctrl)))
             cluster_percent_data_dif<-dplyr::filter_at(cluster_percent_data_dif,vars(one_of(major_cond)),all_vars(!(. %in% heatmap_ctrl)))
             cluster_percent_data_log_ratio<-dplyr::filter_at(cluster_percent_data_log_ratio,vars(one_of(major_cond)),all_vars(!(. %in% heatmap_ctrl)))
             file_tag<-"(HideCtrls)"
      }
      
  
      #=====================输出各个亚群的abundance boxplot================================
      
          #计算输出图片尺寸
          if(is.null(comparisons)) comparisons_n=2 else
            comparisons_n<-length(comparisons)+2
          
          boxplot_n<-length(cluster_names)
          major_cod_ID<-colnames(cluster_percent_data)==major_cond
          cond_n<-length(unique(cluster_percent_data[,major_cod_ID,drop=T]))
          singlewidth<-cond_n*0.5+1
          singleheight<-4+(comparisons_n-2)*0.2
          
          boxplot_ncol<-ceiling(sqrt(boxplot_n*4*singleheight/3/singlewidth))
          boxplot_nrow<-ceiling(boxplot_n/boxplot_ncol)
        
          boxplot_width<-boxplot_ncol*singlewidth*72
          boxplot_height<-boxplot_nrow*singleheight*72
          
          
          if (use_shiny_para){
          boxplot_ncol  =shinny_para$boxplot_ncol
          boxplot_width =shinny_para$boxplot_width  
          boxplot_width =shinny_para$boxplot_width
          boxplot_height=shinny_para$boxplot_height}
        
      
         #定义cluster_boxplot函数
         mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), #坐标系及坐标轴
                           legend.key = element_rect(fill = "white", colour = "white"), #图标
                           legend.background = (element_rect(colour= "white", fill = "white")))
        #cluster_name=cluster_names[1]  
        
         cluster_boxplot<-function(cluster_name,show_pvalues=T,Plot_Data,ylab){
                          ydata<-Plot_Data %>%dplyr::select(cluster_name)
                          yrange<-max(ydata)-min(ydata)
                          ncolor<-nrow(unique(Plot_Data[,color_cond,drop=F]))
                          #ncolor<-nrow(unique(Plot_Data[,color_cond]))
                          
                          cluster_boxplot<-ggplot(Plot_Data,aes_string(x=major_cond,y=cluster_name,fill=color_cond))+
                                            geom_boxplot(outlier.shape= NA,lwd=boxplot_line_width)+
                                            geom_jitter(shape=16, position=position_jitter(0.2))+
                                            mytheme+
                                            labs(title=cluster_name)+
                                            #scale_colour_manual(values=group_colorset(ncolor))+
                                            scale_fill_manual(values=group_colorset(ncolor))+
                            
                                            labs(y=ylab)+
                                            scale_y_continuous(limits=c((min(ydata)-0.1*yrange),(max(ydata)+0.1*yrange*comparisons_n)))+
                                            theme(legend.position = "none")+
                                            theme(axis.text= element_text(angle=label_font_angle,hjust =hjust,vjust = vjust,size=label_font_size))+
                                            theme(axis.title = element_text(size=axis_font_size))+
                                            theme(title =element_text(size=title_font_size,face="bold"))
                          if(show_pvalues & is.null(comparisons)){
                                            cluster_boxplot<-cluster_boxplot+
                                            stat_compare_means(method = stat.method,
                                                               paired=stat.paired,
                                                               ref.group=heatmap_ctrl,
                                                               aes(label = paste0("p = ", ..p.format..)))
                                            }
                                          
                          if(show_pvalues & !is.null(comparisons)){
                                            cluster_boxplot<-cluster_boxplot+
                                            stat_compare_means(method = comparisons.stat.method,
                                                               paired=comparisons.stat.paired,
                                                               comparisons=comparisons,
                                                               aes(label = paste0("p = ", ..p.format..)))
                                            }
                          
          
                          return(cluster_boxplot)
                          }
      
           cluster_boxplot_list<-lapply(cluster_names,cluster_boxplot,show_pvalues=show_pvalues,Plot_Data=cluster_percent_data,ylab="Percentage %")
        
        if(length(heatmap_ctrl)>1){
          
            cluster_boxplot_list_dif<-lapply(cluster_names,cluster_boxplot,show_pvalues=show_pvalues,Plot_Data=cluster_percent_data_dif,ylab="Change of Percentage %")
            cluster_boxplot_list_log_ratio<-lapply(cluster_names,cluster_boxplot,show_pvalues=show_pvalues,Plot_Data=cluster_percent_data_log_ratio,ylab="log2(Percentage Ratio)")
         }
        # if(outputfig)
        #   
        # {    
          wdir=getwd()
          cat(paste0("Output to folder: ",wdir,"/",output_dir,"\n"))
          cat("Outputing cluster abundance boxplot",file_tag,"...\n")
          write.csv(cluster_percent_data,paste0("./",output_dir,"/csvs","/","2.0_cluster_Percentage with group",file_tag,".csv"))
          
        # if(file_format=="PDF"){
          pdf(file=paste0("./",output_dir,"/","2.Cluster Abundance Boxplot",file_tag,".pdf"), width=boxplot_width/72, height=boxplot_height/72)
              multiplot(plotlist=cluster_boxplot_list,cols = boxplot_ncol)
              if(length(heatmap_ctrl)>1){
              multiplot(plotlist=cluster_boxplot_list_dif,cols = boxplot_ncol)
              multiplot(plotlist=cluster_boxplot_list_log_ratio,cols = boxplot_ncol)}
          dev.off()
      
        
      #=====================输出各个亚群的abundance heatmap==================================================================
            
            #定义heatmap函数
            
            draw_heatmap<-function(xdata,
                                  heatmap_color,
                                  key.title){
              
                        #cat("Output heatmap\n")
                        cluster_heatmap_data<-xdata %>%
                                              arrange_at("Short_name") %>%
                                              arrange_at(major_cond)%>%
                                              data.frame()
                        color_cond_names<-cluster_heatmap_data[,color_cond,drop=T]
                        col_color_index<-as.numeric(factor(color_cond_names,levels = unique(color_cond_names),order=T))
                        col_color<-group_colorset(length(unique(col_color_index)))[col_color_index]
                        rownames(cluster_heatmap_data)<-cluster_heatmap_data$Short_name
                        cluster_heatmap_data<-as.matrix(cluster_heatmap_data[,cluster_names])
                        cluster_heatmap_data<-t(cluster_heatmap_data)
                        
                        write.csv(cluster_heatmap_data,paste0("./",output_dir,"/csvs","/",key.title,file_tag,".csv"))
                        x.inv<-try(heatmap.2(cluster_heatmap_data,
                                             Rowv=Rowv,Colv=Colv,dendrogram=dendrogram,
                                             scale="none",
                                             key=T,keysize =0.8,key.title = paste0(key.title),
                                             trace="none",
                                             col=heatmap_color,
                                             density.info="none",
                                             srtCol=45,
                                             #margin=c(25,25),
                                             colCol = col_color,
                                             ColSideColors = col_color,
                                             main=paste0(key.title,"Heatmap",file_tag,"\n"),
                                             na.color="grey"),silent = TRUE)
                        
                        if ('try-error' %in% class(x.inv)) {
                          message("Warning- Could not output",key.title,file_tag,"\n")
                          paste(cluster_heatmap_data)
                        }
                        
                      }
           
            cat("Outputing cluster heatmap",file_tag,"...\n")
            #par(mar=c(7.5,7.5,7.5,7.5)+10) 
            pdf(file=paste0("./",output_dir,"/","1.Cluster Abundance Heatmap",file_tag,".pdf"), width=max(nrow(cluster_percent_data)/5+4,8), height=max(length(cluster_names)/4+3,8))
           
            draw_heatmap(cluster_percent_data,heatmap1_color,key.title="1.1_Cluster Percentage")
            
            if(length(heatmap_ctrl)>1){
            draw_heatmap(cluster_percent_data_dif,heatmap2_color,key.title="1.2_Change of Cluster Percentage")
            draw_heatmap(cluster_percent_data_log_ratio,heatmap2_color,key.title="1.3_log2_Percentage ratio")}
      
            dev.off()
      
      #=====================输出abundance barplot=====================================================================================
      
            cluster_barplot_data<-cluster_percent_data
            cluster_barplot_data_col<-colnames(cluster_barplot_data)
            Short_name_seq<-as.character(unique(cluster_barplot_data$Short_name))
            cluster_barplot_data$Short_name<-factor( cluster_barplot_data$Short_name,levels=Short_name_seq,ordered =T)     
            cluster_barplot_data<-cluster_barplot_data %>%
                                  arrange_at("Short_name") %>%
                                  arrange_at(major_cond)%>%
                                  data.frame()
            color_cond_names<-cluster_barplot_data[,color_cond,drop=T]
            col_color_index<-as.numeric(factor(color_cond_names,levels = unique(color_cond_names),order=T))
            cluster_barplot_data$col_color_index<-col_color_index
      
            
            #画树状图
            if(barplot_clustered==T){      
                        
                                    
                                    row.names(cluster_barplot_data)<- as.character(cluster_barplot_data$Short_name)
                                    model <- hclust(dist(cluster_barplot_data[,cluster_names]), "complete")
                                    dhc <- as.dendrogram(model)
                                    
                                    cluster_barplot_data$Short_name<-factor( cluster_barplot_data$Short_name,levels=labels(dhc),ordered =T)     
                                    cluster_barplot_data<-cluster_barplot_data %>%arrange_at("Short_name") #根据聚类结果重排文件顺序
                              
                                    ddata <- dendro_data(dhc, type = "rectangle")
                                    filenum<-nrow(cluster_barplot_data)
                                    p <-  ggplot(segment(ddata)) + 
                                          geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
                                          theme(plot.margin=unit(c(2,0,0,2), unit = "cm"))+
                                          scale_x_continuous(limits = c(1,filenum+0.01),breaks = NULL,minor_breaks = NULL,expand = c(0.5/filenum,0.5/filenum))+
                                          #coord_flip() + 
                                          scale_y_continuous(expand = c(0, 0))+
                                          mytheme+
                                          theme(panel.background=element_rect(fill='transparent', color='white'))+
                                          theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
      
                                  }
             
             #画barplot
             col_color<-group_colorset(length(unique(cluster_barplot_data$col_color_index)))[cluster_barplot_data$col_color_index]
             id.vars.names<-cluster_barplot_data_col[!(cluster_barplot_data_col %in% cluster_names)]
             #measure.vars.names<-cluster_barplot_data_col[cluster_barplot_data_col %in% cluster_names]
             cluster_barplot_data_melt<-melt(data=cluster_barplot_data[,-ncol(cluster_barplot_data)],
                                            id=id.vars.names,
                                            variable.name = "cluster",
                                            value.name = "cluster_percentage"
                                             )
             
             #cluster_barplot_data_melt$Short_name<-factor(cluster_barplot_data_melt$Short_name,levels =levels(cluster_barplot_data$Short_name),ordered = T )
                  #画出纵向barplot
             cluster_histogram<-ggplot(data=cluster_barplot_data_melt,aes_string(x="Short_name",fill="cluster",y="cluster_percentage")) +
                                geom_bar(position="stack", stat="identity")+
                                #labs(subtitle=paste0("Cluster Percentage"))+
                                mytheme+
                                scale_fill_manual(values=cluster_colorset(clustern))+
                                #scale_x_continuous(limits = c(0,55))+
                                scale_x_discrete(expand = c(0,0))+
                                scale_y_continuous(limits = c(0,100.5),expand = c(0,0))+
                                #theme(axis.text.x= element_text(angle=label_font_angle,hjust =hjust,vjust = vjust+1.5,size=label_font_size,color = col_color))+
                                theme(axis.text.x= element_text(angle=label_font_angle,hjust =hjust,vjust = vjust+1.5,size=label_font_size,color = col_color))+
               
                                #theme(axis.text.x= element_text(angle=-90,size=label_font_size,color = col_color))+
                                xlab("")+
               
                                theme(axis.ticks = element_blank())+
                                theme(plot.margin=unit(c(0,0,2,2), unit = "cm"))+
                                #theme(legend.position = "bottom")+
                                theme(panel.background=element_rect(fill='transparent', color='white'))
             if(barplot_annotation==TRUE){
                                annotation_data<-cluster_barplot_data
                                annotation_data$col_color_index<-as.factor(annotation_data$col_color_index)
                                annotation_bar<-ggplot(data=annotation_data,aes_string(x="Short_name",fill="col_color_index",y=1))+
                                                geom_bar(position="stack", stat="identity")+
                                                mytheme+
                                                scale_fill_manual(values=(group_colorset(length(unique(col_color_index)))))+
                                                theme(legend.position = "none")+
                                                scale_x_discrete(expand = c(0,0))+
                                                scale_y_continuous(expand = c(0,0))+
                                                theme(axis.ticks = element_blank(),axis.text=element_blank(),axis.title = element_blank())+
                                                theme(plot.margin=unit(c(0,0,0,0), unit = "cm"))
                                               
                                annotation_length=0.1
                               }else{
                                annotation_bar<-ggplot()+mytheme
                                annotation_length=0
                                 }
                            
            histo_legend<-as_ggplot(get_legend(cluster_histogram))
            cluster_histogram2=cluster_histogram+theme(legend.position = "none")
              
            cat("Outputing cluster barplot",file_tag,"...\n")
            
              if(barplot_direction=="v"){
                              pdf(file=paste0("./",output_dir,"/","Cluster Abundance Histogram.pdf"), width=max(nrow(cluster_barplot_data)/5+4,10), height=max(length(cluster_names)/4+3,8))
                              #pdf(file=paste0("./",output_dir,"/","Cluster Abundance Histogram.pdf"), width=max(nrow(cluster_barplot_data)/5+4,8), height=max(length(cluster_names)/4+3,8))
                
                                merge_plot1<-plot_grid(p,
                                                       annotation_bar,
                                                       cluster_histogram2,
                                                       align="v",
                                                       ncol=1,
                                                       rel_heights=c(1,annotation_length,4))
                                merge_plot1<-plot_grid(merge_plot1,
                                                       histo_legend,
                                                       align="h",
                                                       ncol=2,
                                                       rel_widths=c(4,1))
                                multiplot(merge_plot1)
                              dev.off()
                              }
            
            
              if(barplot_direction=="h"){
                            #输出横向barplot图
                            suppressMessages(
                                      p<-p+coord_flip()+
                                           scale_y_reverse(expand = c(0, 0))+
                                           theme(plot.margin=unit(c(2,0,0,2), unit = "cm")))
                            
                            cluster_histogram<-ggplot(data=cluster_barplot_data_melt,aes_string(x="Short_name",fill="cluster",y="cluster_percentage")) +
                                                coord_flip()+
                                                geom_bar(position="stack", stat="identity")+
                                                #labs(subtitle=paste0("Cluster Percentage"))+
                                                mytheme+scale_fill_manual(values=cluster_colorset(clustern))+
                                                scale_y_continuous(limits = c(0,100.5),expand = c(0, 0),position = "right")+
                                                scale_x_discrete(position = "top",expand = c(0, 0))+
                                                theme(axis.text.y= element_text(size=label_font_size,hjust=3,color = col_color))+
                                                theme(axis.ticks = element_blank())+
                                                theme(plot.margin=unit(c(2,2,0,0), unit = "cm"))+
                                                theme(legend.position = "bottom")+
                                                theme(panel.background=element_rect(fill='transparent', color='white'))
                   
                            histo_legend<-as_ggplot(get_legend(cluster_histogram))
                            cluster_histogram2=cluster_histogram+theme(legend.position = "none")
      
                            if(barplot_annotation==TRUE){
                                              suppressWarnings(
                                                   annotation_bar<-annotation_bar+
                                                                   coord_flip()+
                                                                   theme(plot.margin=unit(c(0,0.1,0,0), unit = "cm")))
                                                   annotation_length=0.1
                            }else{
                                              annotation_bar<-ggplot()+mytheme
                                              annotation_length=0
                                                 }
              
                            pdf(file=paste0("./",output_dir,"/","3.Cluster Abundance Histogram",file_tag,".pdf"), height=max(nrow(cluster_barplot_data)/5+5,8), width = max(length(cluster_names)/4+3,8))
                         
                                merge_plot1<-plot_grid(p,
                                                       annotation_bar,
                                                       cluster_histogram2,
                                                       align="h",nrow = 1,
                                                       rel_widths=c(1,annotation_length,4))
                                merge_plot1<-plot_grid(merge_plot1,
                                                       histo_legend,
                                                       align="v",ncol=1,
                                                       rel_heights = c(5,1))  
                                multiplot(merge_plot1)
                                dev.off()
                  }
         
            
      #================   
  
}
   if(length(heatmap_ctrl)==1){hide_ctrl=F}

   if(is.null(hide_ctrl)){
     
          if(Colv==T){  abundance_report(hide_ctrl = F)
                        abundance_report(hide_ctrl = T)
          }else{
                        abundance_report(hide_ctrl = F) 
           }     
   }else{
          if(hide_ctrl==T)  {message(paste0("hide_ctrl is set to TRUE,heatmap_ctrl groups were hided in all output figures\n"))}
          abundance_report(hide_ctrl)
   }
   
                
   message("Fig output finished\n")
      

    return_list<-list()
    return_list[["cluster_percent_data"]]<-cluster_percent_data
    return_list[["cluster_percent_ctrl"]]<-cluster_percent_ctrl
    
    return(return_list)} 





#'@title 生成Cluster Aundance的火山图
#'@description
#'依据stat_by_cluster生成的list数据，生成火山图，图中颜色代表不同的Marker，数字和字母等符号代表不同的cluster；
#'结果采用了与SPADEvizR相同的形式，标注有“cite from SPADEVizR”或者“Adapted from SPADEVizR” ，使用时请注意在文章中引用原始文献。


#'
#'@param cluster_stat       stat_by_cluster()函数生成的list数据
#'@param cond1              major_cond的一个，火山图（Volcano map）靠左一侧的条件，例如"B_Adjacent"
#'@param cond2              major_cond的一个，火山图（Volcano map）靠右一侧的条件，例如"A_Tumor"

#'@param xlimit             一个包含两个数字的向量（vector），例如c(-5,+5),用于手动指定火山图横坐标范围，默认值为NULL，此时会自动设置xlimit
#'@param ylimit             一个包含两个数字的向量（vector），例如c(0,3.5),用于手动指定火山图纵坐标范围，默认值为NULL，此时会自动设置xlimit
#'@param dif.level          一个数字，手动设置两组marker差异的的阈值大小，例如取默认值2时，两组均值相差2倍以上算显著（简单的说就是设置火山图上的两条竖线位置）
#'@param output_dir         输出数据文件夹名称
#'
#'如使用全局统计默认参数，以下参数无需设置
#'
#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著），调整该数值可以调整火山图上的横线位置
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#' @return                none
#' @export


draw_abundance_volcano<-function(cluster_stat,
                            cond1,
                            cond2,
                            xlimit=NULL,
                            ylimit=NULL ,
                            dif.level=2,
                            output_dir=paste0("Aundance_volcano_plots_",cond1,"_vs_",cond2),
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



  if(stat.paired){
    if_paired<-"paired"}else{
      if_paired<-"unpaired"}

  cluster_boxplot_data<-full_join(cluster_percent_summary,groups,by="File_ID")

    
    
  #计算Cluster Ratio

  cond1_cluster<-
    cluster_boxplot_data %>%
    #filter_(paste0(major_cond,"==cond1"))
       dplyr::filter_at(vars(one_of(major_cond)),all_vars(.==cond1)) 

  cond2_cluster<-
    cluster_boxplot_data %>%
    #filter_(paste0(major_cond,"==cond2"))
       dplyr::filter_at(vars(one_of(major_cond)),all_vars(.==cond2)) 
    
  


  volcano_stat_result<-combined_stat(cond1_data=cond1_cluster[,colnames(cluster_percent_summary)[-1]],
                                     cond2_data=cond2_cluster[,colnames(cluster_percent_summary)[-1]],
                                     alternative="two.sided",
                                     stat.paired=stat.paired,
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

  log_ratio=log2((cond2_data_m+No_zero)/(cond1_data_m+No_zero))

  cluster_volcano_data<-data.frame(cluster=paste0("Cluster",c(1:nrow(volcano_stat_result))),
                                   log_ratio=t(log_ratio),
                                   volcano_stat_result)

  cluster_percent_mean<-apply(cluster_percent_summary[,-1],2,mean)

  #为了方便绘图，所有小于0.0001的p值等于0.0001
  cluster_volcano_data$p_values[cluster_volcano_data$p_values<0.0001]=0.0001
  cluster_volcano_data$p_adjust[cluster_volcano_data$p_adjust<0.0001]=0.0001
  cluster_volcano_data$log_p_value=-log10(cluster_volcano_data$p_value)
  cluster_volcano_data$log_p_adjust=-log10(cluster_volcano_data$p_adjust)
  cluster_volcano_data$cluster_percent_mean = cluster_percent_mean
  cluster_volcano_data$Significant<-ordered(cluster_volcano_data$p_values<1-conf.level & abs(cluster_volcano_data$log_ratio)>=log2(dif.level),levels=c(FALSE,TRUE))
  cluster_volcano_data$cluster_num<-sub("cluster","",row.names(cluster_volcano_data))

  x_threshold<-round(log2(dif.level),digits=1)

  if(is.null(xlimit)){
    xuplimit=max(abs(log_ratio)+2)
    xlimit=c(-1*xuplimit,xuplimit)}

  if(is.null(ylimit)){
    yuplimit=min(max(cluster_volcano_data$log_p_value+1.5),4)
    ylimit=c(0,yuplimit)}

  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))



  cluster_volcano<-ggplot(cluster_volcano_data)+
    geom_point(aes_string(x="log_ratio",y=p_para,size="cluster_percent_mean",color="Significant"))+
    labs(y="-Log10(p_values)",title="Differential Abundant Clusters")+
    labs(x=paste0("Log2_Ratio","\n", cond1, " <- enriched -> ", cond2,"\n\n Cite from SPADEVizR"))+
    labs(subtitle=paste0(adjust_label,"p values of cluster percentage were calculated by ",if_paired," ",stat.method))+
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
    geom_text_repel(aes_string(x="log_ratio",y=p_para),label=cluster_volcano_data$cluster_num,box.padding=0.25,point.padding=0.5)

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
  write.csv(FDR_pValue_arranged,paste0("./",output_dir,"/Cluster Abundance Statistical Report(arranged).csv"),row.names = TRUE)


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

  cat(paste0("Output Cluster Abundance Volcano plot...\n"))
  pdf(file=paste0("./",output_dir,"/Cluster Abundance Volcano_",cond1," enriched ",cond2,".pdf"), width=15, height=10)
  multiplot(cluster_volcano,FDR_pValue_plot,cols = 2)
  #multiplot(FDR_table)
  dev.off()

}




