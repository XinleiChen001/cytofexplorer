
#2019_12_20 对cluster_merge的功能进行完善,支持将所有cluster合并为一个进行分析(即对整个文件进行统计)
#2019_12_20 对clusterlabel进行改进,保持合并后cluster与label的对应关系



#'@title 生成marker表达的火山图
#'@description
#'依据stat_by_cluster生成的list数据，生成火山图，图中颜色代表不同的Marker，数字和字母等符号代表不同的cluster；
#'
#'@param cluster_stat       stat_by_cluster()函数生成的list数据
#'@param cond1              major_cond的一个，火山图（Volcano map）靠左一侧的条件
#'@param cond2              major_cond的一个，火山图（Volcano map）靠右一侧的条件
#'@param Usecolorset        火山图采用的配色方案 默认dif_seq_rainbow, 其他取值：rainbow，primary.color或者其他函数
#'@param show_adjust_p      一个逻辑变量， 取值TRUE时火山图显示校正后的p值，取值FALSE时显示原始p值
#'@param cluster_id         指定输出cluster的id，例如：如果要输出前五个cluster，cluster_id=c(1:5); 如果要输出全部cluster，保持其默认值NULL即可。
#'@param xlimit             一个包含两个数字的向量（vector），例如c(-5,+5),用于手动指定火山图横坐标范围，默认值为NULL，此时会自动设置xlimit
#'@param ylimit             一个包含两个数字的向量（vector），例如c(0,3.5),用于手动指定火山图纵坐标范围，默认值为NULL，此时会自动设置xlimit
#'@param dif.level          一个数字，手动设置两组marker差异的的阈值大小，例如取默认值2时，两组均值相差2倍以上算显著（简单的说就是设置火山图上的两条竖线位置）
#'@param output_dir         输出数据文件夹名称

#'如使用全局统计默认参数，以下参数无需设置
#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著），调整该数值可以调整火山图上的横线位置
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#' @return                none
#' @export


draw_expr_volcano2<-function(cluster_stat,
                             expr_metadata,
                             cond1,
                             cond2,
                            Usecolorset=dif_seq_rainbow,
                            show_adjust_p=FALSE,
                            cluster_id=NULL,
                            xlimit=NULL,
                            ylimit=NULL,
                            dif.level=2,
                            # subgroup_cond = NULL,
                            # heatmap_ctrl=NULL,
                            # groups_to_show=NULL,
                            output_dir=paste0("Expr_volcano_plots2_",cond1,"_vs_",cond2),
                            volcano.stat.paired=NULL,
                            #如使用全局统计默认参数，以下参数无需设置
                            stat.paired=NULL,
                            stat.method=NULL,
                            set.equal.var=NULL,
                            conf.level=NULL,
                            p.adjust.method=NULL
                            ){
  

  
  if(0){
    

    Usecolorset=dif_seq_rainbow
    show_adjust_p=FALSE
    cluster_id=NULL
    xlimit=NULL
    ylimit=NULL
    dif.level=2
    # subgroup_cond = NULL,
    # heatmap_ctrl=NULL,
    # groups_to_show=NULL,
    output_dir=paste0("Expr_volcano_plots2_",cond1,"_vs_",cond2)
    volcano.stat.paired=NULL
    #如使用全局统计默认参数，以下参数无需设置
    stat.paired=NULL
    stat.method=NULL
    set.equal.var=NULL
    conf.level=NULL
    p.adjust.method=NULL
    
  }
  
  
  #expr_metadata=expr_data
  
    
# 
# 
#   Usecolorset=primary.colors
#   cond1="Tumor (A)"
#   cond2="Adjacent(B)"
# 
# 
#   Usecolorset=dif_seq_rainbow
#   show_adjust_p=FALSE
#   cluster_id=NULL
#   xlimit=NULL
#   ylimit=NULL
#   dif.level=2
#   output_dir=paste0("Expr_volcano_plots_",cond1,"_vs_",cond2)
# 
#   #如使用全局统计默认参数，以下参数无需设置
#   stat.paired=NULL
#   stat.method=NULL
#   set.equal.var=NULL
#   conf.level=NULL
#   p.adjust.method=NULL

  
  clustern          <- cluster_stat[["metadata"]]$clustern
  summerise_method          <- cluster_stat[["metadata"]]$summerise_method
  major_cond             <- cluster_stat[["metadata"]]$major_cond
  all_markers             <- cluster_stat[["all_markers"]]
  merged_cluster_id      <-cluster_stat[["merged_cluster_id"]]
  
  # if(is.null(groups_to_show)){
  #   groups_to_show=unique(groups[,major_cond,drop=T])  
  # } 
  
  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }
  
  wdir=getwd()
  
  cat(paste0("Output to folder: ",wdir,"/",output_dir,"\n"))
  
  ht_markers<-as.character(subset(all_markers,expr_para==1)$markers)
  
  major_cond<-as.character(major_cond)
  shape_index<-c(49:57,65:90,97:122)  #数字，大写字母，小写字母  shape
  shape_char<-intToUtf8(shape_index,multiple=T)
  
  
  cluster_expr_raw  <-cluster_stat[["cluster_expr_raw"]]
  if(is.null(stat.paired)) stat.paired        <- cluster_stat[["metadata"]]$stat.paired
  if(is.null(stat.method)) stat.method        <- as.character(cluster_stat[["metadata"]]$stat.method)
  if(is.null(conf.level)) conf.level          <- cluster_stat[["metadata"]]$conf.level
  if(is.null(set.equal.var)) set.equal.var    <- as.character(cluster_stat[["metadata"]]$set.equal.var)
  if(is.null(p.adjust.method)) p.adjust.method    <- as.character(cluster_stat[["metadata"]]$p.adjust.method)
  
  if(is.null(cluster_id)){
    cluster_id=merged_cluster_id
  }
  
  major_cond       <-as.character(cluster_stat[["metadata"]]$major_cond)
  cluster_name       <-as.character(cluster_stat[["metadata"]]$cluster_name)
  
  
  if(is.null(volcano.stat.paired)) volcano.stat.paired<-  stat.paired
  
  if(summerise_method=="mean"){
    use_method<-mean}else
      if(summerise_method=="median")
      { use_method<-median
      }
  
  
  if(show_adjust_p==TRUE){
    p_para<-"log_p_adjust"
    adjust_label="Adjust "}else{
      p_para<-"log_p_values"
      adjust_label=""
    }
  ht_markers<-as.character(subset(all_markers,expr_para==1)$markers)
  trans_markers<-as.character(subset(all_markers,transform==1)$markers)
  used_markers<-union(ht_markers,trans_markers)
  
  if(stat.method=="auto") stat.method="t-test or wilcox-test\n(Statistical method is selected automaticly)"
  
  if(volcano.stat.paired){
    if_paired<-"paired"}else{
      if_paired<-"unpaired"}
  
  
  cond_name<-major_cond
  cond_name<-as.character(cond_name)
  #画火山图
  volcano_data_all<-list()
  
  
  
  expr_ctrl_merge<-expr_metadata$expr_ctrl_merge %>%
    arrange_at("Short_name") %>%
    arrange_at(cond_name) 
  
  expr_trans_merge<-expr_metadata$expr_trans_merge  %>%
    arrange_at("Short_name") %>%
    arrange_at(cond_name) 
  
  expr_trans_ratio_merge<-expr_trans_merge
  expr_trans_ratio_merge[,ht_markers]<-expr_trans_merge[,ht_markers]-expr_ctrl_merge[,ht_markers]
  #View(expr_trans_merge)
  if(is.null(cluster_id)){
    cluster_id=unique(expr_trans_merge[,cluster_name,drop=T])
  }
  tail(expr_trans_merge$metacluster)
  for(cluster_i in cluster_id){
    
    #cluster_i=7
    # 
    # expr_trans<-cluster_expr_trans[[cluster_i]] %>%
    #   dplyr::filter_at(vars(matches(major_cond)),all_vars(.%in% groups_to_show))
    # 
    # expr_raw<-cluster_expr_raw[[cluster_i]] %>%
    #   dplyr::filter_at(vars(matches(major_cond)),all_vars(.%in% groups_to_show))
    # 
    # if (is.null(subgroup_cond))  {
    #   expr_raw$whole<-"Yes"
    #   expr_trans$whole<-"Yes"
    #   subgroup_cond<-"whole"
    # } else if(subgroup_cond=="whole"){
    #   expr_raw$whole<-"Yes"
    #   expr_trans$whole<-"Yes"
    #   
    # }
    # 
    # expr_ctrl=NULL
    # 
    # if(length(heatmap_ctrl)!=nrow(unique(expr_raw[,subgroup_cond]))) message("Error:Number of Heatmap_ctrl is not equal to number of subgroup_cond")
    # 
    # 
    # for(subgroup_cond_i in t(unique(expr_raw[,subgroup_cond]))){
      
    #   #subgroup_cond_i<-"Term"
    #   expr_subgroup_cond<-as.data.frame(matrix(0,ncol(expr_trans),nrow=0))
    #   colnames(expr_subgroup_cond)<-colnames(expr_trans)
    #   
    #   expr_trans_cond<-dplyr::filter_(expr_trans,paste0(subgroup_cond,"==subgroup_cond_i"))
    #   #expr_trans_cond<-dplyr::filter_at(expr_trans,vars(matches(subgroup_cond)),all_vars(.==subgroup_cond_i))
    #   
    #   
    #   heatmap_ctrl_i<-intersect(unique(expr_trans_cond[,cond_name,drop=T]),heatmap_ctrl)
    #   
    #   #计算paired samples
    #   if(stat.paired)
    #   {for(conditions in t(unique(expr_trans_cond[,cond_name]))){
    #     
    #     
    #     
    #     # expr_ctrl_add<-dplyr::filter_(expr_trans_cond,paste0(cond_name,"==conditions"))
    #     
    #     expr_ctrl_add<-dplyr::filter_at(expr_trans_cond,vars(matches(cond_name)),all_vars(.==conditions))
    #     
    #     #single_subgroup_cond<-dplyr::filter_(expr_trans_cond,paste0(cond_name,"==heatmap_ctrl_i"))
    #     single_subgroup_cond<-dplyr::filter_at(expr_trans_cond,vars(matches(cond_name)),all_vars(.==heatmap_ctrl_i))
    #     
    #     not_groups_para_id<-which(!(colnames(expr_ctrl_add)%in% colnames(groups)))
    #     expr_ctrl_add[,not_groups_para_id]<-single_subgroup_cond[,not_groups_para_id]
    #     expr_subgroup_cond<-rbind.data.frame(expr_subgroup_cond,expr_ctrl_add)}
    #   } else{
    #     #unpaired_control_cond<-dplyr::filter_(expr_trans_cond,paste0(cond_name,"==heatmap_ctrl_i"))
    #     unpaired_control_cond<-dplyr::filter_at(expr_trans_cond,vars(matches(cond_name)),all_vars(.==heatmap_ctrl_i))
    #     
    #     single_subgroup_cond<-unpaired_control_cond[1,]
    #     single_subgroup_cond[,used_markers]<-summarise_if(unpaired_control_cond[,used_markers],is.numeric,funs(use_method))
    #     expr_subgroup_cond<-expr_trans_cond
    #     for(row_i in c(1:nrow(expr_subgroup_cond))){
    #       expr_subgroup_cond[row_i,used_markers]<-single_subgroup_cond[,used_markers]
    #     }
    #   }
    #   
    #   if(is.null(expr_ctrl))expr_ctrl<-expr_subgroup_cond
    #   else{ expr_ctrl<-rbind.data.frame(expr_ctrl,expr_subgroup_cond) }
    #   
    # }
    # #str(expr_subgroup_cond)
    # 
    # expr_ctrl<-
    #   expr_ctrl %>%
    #   arrange_at("Short_name") %>%
    #   arrange_at(cond_name)
    # 
    # # str(expr_ctrl)
    # # str(expr_trans)
    # 
    # ht_expr_trans_ratio<-expr_trans[,ht_markers]-expr_ctrl[,ht_markers]
    # rownames(ht_expr_trans_ratio)<-expr_trans$Short_name
    # 
    # expr_ratio_boxplot<-ht_expr_trans_ratio
    # expr_ratio_boxplot$Short_name<-rownames(ht_expr_trans_ratio)
    # expr_ratio_boxplot<-left_join(expr_ratio_boxplot,groups,by="Short_name")
    # rownames(expr_ratio_boxplot)<-expr_ratio_boxplot$Short_name
    # 

    #View(expr_trans_ratio_merge)

    
    cond1_data<-
      expr_trans_ratio_merge %>%
      #dplyr::filter_(paste0(major_cond,"==cond1"))%>%
      dplyr::filter_at(vars(matches(major_cond)),all_vars(.== cond1)) %>%
      dplyr::filter_at(vars(matches(cluster_name)),all_vars(.== cluster_i)) %>%
      dplyr::select(ht_markers)
      
      cond2_data<-
        expr_trans_ratio_merge %>%
        #dplyr::filter_(paste0(major_cond,"==cond1"))%>%
        dplyr::filter_at(vars(matches(major_cond)),all_vars(.== cond2)) %>%
        dplyr::filter_at(vars(matches("metacluster")),all_vars(.== cluster_i)) %>%
        dplyr::select(ht_markers)

    
    #write.csv(expr_ratio_boxplot,"expr_arcsih_ratio_merge.csv")
    
    
    volcano_cluster_data<-combined_stat(cond1_data,
                                        cond2_data,
                                        alternative="two.sided",
                                        mu=0,
                                        transform.method="raw",
                                        stat.paired=volcano.stat.paired,
                                        conf.level=conf.level,
                                        stat.method=stat.method,
                                        set.equal.var=set.equal.var,
                                        p.adjust.method=p.adjust.method,
                                        silent=TRUE)
    
    
    #计算Ratio
    
    cond1_data_m<-
      cond1_data %>%
      summarise_if(is.numeric,use_method,na.rm=F)
    
    #summarise_each(funs(use_method))
    
    cond2_data_m<-
      cond2_data %>%
      summarise_if(is.numeric,use_method,na.rm=F)
    #summarise_each(funs(use_method))
    
    #No_zero<-0.25  #add to prevent divided by zero
    
    #arcsinh_ratio_difference=log2((cond2_data_m+No_zero)/(cond1_data_m+No_zero))
    arcsinh_ratio_difference=cond2_data_m-cond1_data_m
    
    
    volcano_data<-data.frame(marker=colnames(cond1_data),
                             cluster_i=cluster_i,
                             cluster_shape=shape_char[cluster_i],
                             volcano_cluster_data,
                             arcsinh_ratio_difference=t(arcsinh_ratio_difference)
    )
    
    
    if(cluster_i==1){volcano_data_all=volcano_data}else{
      volcano_data_all<-rbind(volcano_data_all,volcano_data)
    }
    
    
  }
  
  
  
  
  cat(paste0("Output Statistical Report.csv...\n"))
  write.csv(volcano_data_all,paste0("./",output_dir,"/","Statistical Report.csv"),row.names = TRUE)
  
  shape_index<-c(49:57,65:90,97:122)  #数字，大写字母，小写字母  shape
  shape_char<-intToUtf8(shape_index,multiple=T)
  
  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))
  
  
  #准备FDR-pValue画图
  FDR_ranks<-factor(c("other(adjust_p>=10%)","adjust_p<10%","adjust_p<5%"),levels=c("other(adjust_p>=10%)","adjust_p<10%","adjust_p<5%"),ordered=TRUE)
  FDR_pValue<-dplyr::arrange(volcano_data_all,p_values)
  FDR_pValue$row_id<-factor(c(1:nrow(FDR_pValue)))
  FDR_pValue$FDR_rank<-FDR_ranks[1]
  FDR_pValue$FDR_rank[FDR_pValue$p_adjust<0.1]=FDR_ranks[2]
  FDR_pValue$FDR_rank[FDR_pValue$p_adjust<0.05]=FDR_ranks[3]
  cat(paste0("Output Statistical Report(arranged).csv...\n"))
  write.csv(FDR_pValue,paste0("./",output_dir,"/","Statistical Report(arranged).csv"),row.names = TRUE)
  
  FDR_table_data<-FDR_pValue[(FDR_pValue$p_values<0.05)&(!is.na(FDR_pValue$p_values)),c(1:6,8,10)]
  round_2<-function (var){return(round(var,3))}
  if(nrow(FDR_table_data)>0)
    FDR_table_data[,c(4,5,7)]<-apply(FDR_table_data[,c(4,5,7)],2,round_2)
  
  FDR_table_plot<-ggtable(FDR_table_data,show_row_name = F)
  
  
  FDR_pValue<-dplyr::arrange(FDR_pValue,desc(p_values))
  FDR_pValue$row_id<-factor(c(1:nrow(FDR_pValue)))
  FDR_pValue$p_values=-log10(FDR_pValue$p_values)
  # FDR_pValue$FDR_rank<-factor(FDR_pValue$FDR_rank)
  
  FDR_pValue<-na.omit(FDR_pValue)
  
  FDR_pValue_plot<-ggplot(FDR_pValue)+
    geom_bar(aes_string(x="row_id",y="p_values",fill="FDR_rank"),stat="identity",width = 0.8)+
    labs(y="-Log10(p_values)",title="Summary of Adjust p Value")+
    labs(x=paste0(""))+
    coord_flip()+
    labs(subtitle=paste0("Adjust p values were separately calculated based on the p values of each cluster "))+
    mytheme+
    #expand_limits(y=c(0,3))+
    theme(panel.background=element_rect(colour='white',size = 1 ))+
    theme(axis.line.x=element_line(colour='black',linetype = NULL,size=0.3))+
    theme(axis.ticks.y=element_blank())+
    scale_fill_manual(values = c("lightgreen","navy","red"))+
    scale_y_continuous(breaks=c(seq(0,3,by=0.5)),minor_breaks = NULL)+
    theme(legend.position = c(0.7,0.4),legend.title=element_blank())+
    theme(axis.text.y=element_blank())+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))
  
  
  
  #为了方便绘图，所有小于0.0001的p值等于0.0001
  volcano_data_all$p_values[volcano_data_all$p_values<0.0001]=0.0001
  volcano_data_all$p_adjust[volcano_data_all$p_adjust<0.0001]=0.0001
  volcano_data_all$log_p_values=-log10(volcano_data_all$p_values)
  volcano_data_all$log_p_adjust=-log10(volcano_data_all$p_adjust)
  
  volcano_data_all$cluster_i
  
  #volcano_data_all$cluster_i<-as.integer(volcano_data_all$cluster_i)
  #clustern<-length(unique(volcano_data_all$cluster_i))
  # shape_list<-data.frame(cluster_i=c(1:clustern),
  #                        shape_index=shape_index[1:clustern])
  
  shape_list<-data.frame(cluster_i=merged_cluster_id,
                         shape_index=shape_index[merged_cluster_id]) 
  volcano_data_all<-left_join(volcano_data_all,shape_list,by="cluster_i")
  
  marker_n<-length(unique(volcano_data_all$marker))

  
  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))
  
  
  
  
  cat(paste0("Output Volcanoplot...\n"))
  
  x_threshold<-round(log2(dif.level),digits=1)
  
  
  if(is.null(xlimit)){
    xuplimit=max(abs(volcano_data_all$arcsinh_ratio_difference)+2,na.rm = TRUE)
    xlimit=c(-1*xuplimit,xuplimit)}
  
  
  Y_ID<-colnames(volcano_data_all)==p_para
  
  if(is.null(ylimit)){
    yuplimit=min(max(volcano_data_all[Y_ID]+1.5,na.rm = TRUE),4)
    ylimit=c(0,yuplimit)}
  
  head(volcano_data_all)
  volcanoplot<-ggplot(volcano_data_all)+
    scale_shape_identity() +
    geom_point(aes_string(x="arcsinh_ratio_difference",y=p_para,colour="marker",shape="shape_index"),size=4)+
    labs(y="-Log10(p_values)",title="Differential response of Markers(Arcsinh ratio comparition) ")+
    labs(x=paste0("Difference","\n", cond1, " <- versus -> ", cond2))+
    labs(subtitle=paste0(adjust_label,"p values of marker expression arcsinh ratio were calculated by ",if_paired," ",stat.method))+
    mytheme+
    expand_limits(x=xlimit,y=ylimit)+
    geom_vline(xintercept = c(x_threshold,-1*x_threshold), linetype="dashed",color = "blue", size=0.55)+
    geom_hline(yintercept=-log10(1-conf.level), linetype="dashed", color = "blue",size=0.55)+
    scale_x_continuous(breaks=c(seq(-5,5,by=1)),minor_breaks = NULL)+
    scale_y_continuous(breaks=c(seq(0,6,by=1),1.3),minor_breaks = NULL)+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
    theme(panel.background=element_rect(colour='black',size = 1 ))+
    theme(axis.line=element_line(colour='black',linetype = 1,size=0.3))+
    scale_color_manual(values = Usecolorset(marker_n))
  
  shape_table<- t(data.frame(Cluster=shape_list$cluster_i,
                             Symbol=intToUtf8(shape_list$shape_index,multiple = T)))
  row.names(shape_table)<-c("Cluster","Symbol")
  
  volcano_note<-ggtable(shape_table,show_row_name=TRUE,show_col_name=F)
  # layout<-cbind(matrix(c(rep(1,20),2,2,0,0),ncol = 2, byrow = T),
  #          matrix(c(rep(3,12),rep(4,10),0,0),ncol = 2, byrow = T))
  
  
  layout<-cbind(matrix(c(rep(1,20),2,2,0,0),ncol = 2, byrow = T),
                matrix(c(rep(3,12),rep(3,10),0,0),ncol = 2, byrow = T))
  
  
  
  
  pdf(file=paste0("./",output_dir,"/","Volcano_",cond1," vs ",cond2,".pdf"), width=25, height=15)
  
  multiplot(volcanoplot,volcano_note,FDR_pValue_plot,layout=layout)
  #multiplot(FDR_table_plot)
  
  dev.off()
  
  cat(paste0("Output Stastic Heatmap ...\n"))
  head(volcano_data_all)
  heatmap_data<-tidyr::spread(volcano_data_all[,c("marker","cluster_i","log_p_values")],"cluster_i","log_p_values")
  rownames(heatmap_data)<-heatmap_data$marker
  heatmap_data$marker<-NULL
  heatmap_data<-as.matrix(heatmap_data)

  
  # pheatmap(heatmap_data)
  # 
  cellnotes<-as.character(round(heatmap_data,2))
  cellnotes[heatmap_data<(-log10(1-conf.level))|is.na(heatmap_data)|is.nan(heatmap_data)]<-""
  cellnotes<-matrix(data=cellnotes,ncol=ncol(heatmap_data))
  heatmap_cms<-"Significant -log p values were marked in the heatmap."
  if (sum(cellnotes!="")==0) heatmap_cms<-"No Significant Response was found."

  write.csv(heatmap_data,paste0("./",output_dir,"/","Stastic heatmap data.csv"),row.names = TRUE)
  pdf(file=paste0("./",output_dir,"/","heatmap_",cond1," vs ",cond2,".pdf"), width=15, height=15)
  
  
  x.inv<-try(heatmap.2(heatmap_data,
            Rowv=F,Colv=F,dendrogram="none",
            scale="none",
            key=T,keysize =1,key.title = "",
            trace="none",
            col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
            srtCol=0,
            margin=c(15,15),
            #colCol = mg_col_color,
            #ColSideColors = mg_col_color,
            main=paste0("-log p values",heatmap_cms),
            na.color="grey",
            cellnote=cellnotes),silent = T)
  if ('try-error' %in% class(x.inv)) {
    message("Warning- Could not output Expression Heatmap\n")
    heatmap_data}
 
  
  
  heatmap_ratio_data<-tidyr::spread(volcano_data_all[,c("marker","cluster_i","arcsinh_ratio_difference")],"cluster_i","arcsinh_ratio_difference")
  rownames(heatmap_ratio_data)<-heatmap_ratio_data$marker
  heatmap_ratio_data$marker<-NULL
  heatmap_ratio_data<-as.matrix(heatmap_ratio_data)
  cellnotes<-as.character(round(heatmap_ratio_data,2))
  cellnotes[heatmap_data<(-log10(1-conf.level))|is.na(heatmap_data)|is.nan(heatmap_data)]<-""
  cellnotes<-matrix(data=cellnotes,ncol=ncol(heatmap_data))
  x.inv<-try(heatmap.2(heatmap_ratio_data,
            Rowv=F,Colv=F,dendrogram="none",
            scale="none",
            key=T,keysize =1,key.title = "",
            trace="none",
            col=colorRampPalette(c("blue","black","yellow"))(100),
            srtCol=0,
            margin=c(15,15),
            #colCol = mg_col_color,
            #ColSideColors = mg_col_color,
            main=paste0("marker_arcsinh_ratio_differences"),
            na.color="grey",
            cellnote=cellnotes),silent = T)
  
  if ('try-error' %in% class(x.inv)) {
    message("Warning- Could not output Expression Heatmap\n")
    heatmap_ratio_data}
  

  
  dev.off()
 
  
  cat(paste0("Output Volcanoplot facet by markers...\n"))
  
  pdf(file=paste0("./",output_dir,"/","Volcano_",cond1," vs ",cond2,"by markers.pdf"),
      width=20, height=20)
  
  
  volcanoplot2<-ggplot(volcano_data_all)+
    scale_shape_identity() +
    geom_point(aes_string(x="arcsinh_ratio_difference",y=p_para,colour="marker",shape="shape_index"),size=4)+
    #geom_point(aes(x=arcsinh_ratio_difference,y=p_adjust,colour=marker_name,shape=shape_index),size=4)+
    labs(y="-Log10(p_values)",title="Differential Expression of Markers")+
    labs(x=paste0("Difference","\n", cond1, " <- versus -> ", cond2))+
    labs(subtitle=paste0(adjust_label,"p values of arcsinh transformed expression were calculated by ",if_paired," ",stat.method))+
    mytheme+
    expand_limits(x=xlimit,y=ylimit)+
    geom_vline(xintercept = c(-1,+1), linetype="dashed",color = "blue", size=0.55)+
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "blue",size=0.55)+
    scale_x_continuous(breaks=c(seq(-5,5,by=1)),minor_breaks = NULL)+
    scale_y_continuous(breaks=c(seq(0,6,by=1),1.3),minor_breaks = NULL)+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
    theme(panel.background=element_rect(colour='black',size = 1 ))+
    theme(axis.line=element_line(colour='black',linetype = 1,size=0.3))+
    scale_color_manual(values = Usecolorset(marker_n))+
    facet_wrap(~marker)
  
  multiplot(volcanoplot2)
  dev.off()
  
  
  cat(paste0("Output Volcanoplot facet by clusters...\n"))
  
  pdf(file=paste0("./",output_dir,"/","Volcano_",cond1," vs ",cond2,"by population.pdf"),
      width=20, height=20)
  
  
  volcanoplot3<-ggplot(volcano_data_all)+
    scale_shape_identity() +
    #geom_point(aes(x=arcsinh_ratio_difference,y=p_adjust,colour=marker_name,shape=shape_index),size=4)+
    geom_point(aes_string(x="arcsinh_ratio_difference",y=p_para,colour="marker",shape="shape_index"),size=4)+
    labs(y="-Log10(p_values)",title="Differential Expression of Markers")+
    labs(x=paste0("Difference","\n", cond1, " <- versus -> ", cond2))+
    labs(subtitle=paste0(adjust_label,"p values of arcsinh transformed expression were calculated by ",if_paired," ",stat.method))+
    mytheme+
    expand_limits(x=xlimit,y=ylimit)+
    geom_vline(xintercept = c(-1,+1), linetype="dashed",color = "blue", size=0.55)+
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "blue",size=0.55)+
    scale_x_continuous(breaks=c(seq(-5,5,by=1)),minor_breaks = NULL)+
    scale_y_continuous(breaks=c(seq(0,6,by=1),1.3),minor_breaks = NULL)+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
    theme(panel.background=element_rect(colour='black',size = 1 ))+
    theme(axis.line=element_line(colour='black',linetype = 1,size=0.3))+
    scale_color_manual(values = Usecolorset(marker_n))+
    facet_wrap(~cluster_i)
  
  multiplot(volcanoplot3)
  dev.off()
  
}




