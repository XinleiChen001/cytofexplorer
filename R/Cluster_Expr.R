



#' @title combined_stat 统计分析函数
#'
#' @description
#' 对进行单组或两组数据进行统计分析，并对结果进行adjust。也可以为数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#'
#' @param cond1_data      一个包含数字的向量（vector），待统计分析数据
#' @param cond2_data      一个包含数字的向量（vector），待统计分析数据，如进行单组分析，cond2_data=NULL
#' @param alternative     一个字符串，反应假设检验的种类,有三种取值："two.sided" (默认), "greater" or "less"
#' @param mu              一个数字，在单样本本统计时表示平均值的真实值，或者如果您正在进行两个样本测试，则表示均值差异（此时一般取值为0）
#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著）
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#' @param  silent         逻辑变量，是否输出中间结果
#' @return 返回一个数据框（data.frame),包含p_value,p_adjust,统计方法(stat.method),方差齐性（set.equal.var）信息
#' @export

combined_stat<-function(cond1_data,
                       cond2_data,
                       alternative="two.sided",
                       mu=0,
                       stat.paired=FALSE,
                       conf.level=0.95,
                       stat.method="t-test",
                       set.equal.var="FALSE",
                       p.adjust.method="BH",
                       silent=F){
s_methods<-NULL
p_values<-NULL
var.equals<-NULL
if(stat.method=="t-test"){
  s_method<-"t-test"}else if(stat.method=="wilcox-test"){
    s_method<-"wilcox-test"}





for(i in 1:ncol(cond1_data)){
  if(silent==F) cat(paste0("Treating with cluster",i,":  "))
  var.equal<-NA

  if(stat.method=="auto"){
    #判断是否是正态分布

    shapiro_result1=shapiro.test(cond1_data[,i,drop=T])$p.value
    if(!is.null(cond2_data)){
    shapiro_result2=shapiro.test(cond2_data[,i,drop=T])$p.value}else{
    shapiro_result2=1
    }

    if(shapiro_result1>0.05 & shapiro_result2>0.05){
      s_method<-"t-test" }else{
        s_method<-"wilcox-test"}
  }
  #t-test
  if(s_method=="t-test"){


    #方差齐性检测
    if(set.equal.var=="auto")
    {var.equal=FALSE
    if(!is.null(cond2_data)){

    pool_data<-c(cond1_data[,i,drop=T],cond2_data[,i,drop=T])
    pool_group<-as.factor(c(rep(1,length(cond1_data[,i,drop=T])),rep(2,length(cond2_data[,i,drop=T]))))
    equal_var<-leveneTest(y =pool_data, group=pool_group)
    if(equal_var$`F value`[1]<0.05) {
      var.equal=FALSE }else{
        var.equal=TRUE
      }
    }
    }else if(set.equal.var=="TRUE"){
             var.equal=TRUE} else if((set.equal.var=="FALSE") ) {
               var.equal=FALSE}

    cond1_with_na<- length(na.omit(cond1_data[,i,drop=T]))==length(cond1_data[,i,drop=T])
    cond2_with_na<- suppressWarnings(length(na.omit(cond2_data[,i,drop=T]))==length(cond2_data[,i,drop=T]))
    if(cond1_with_na & cond2_with_na ){
      pv<-t.test(x=cond1_data[,i,drop=T],
                 y=cond2_data[,i,drop=T],
                 alternative = alternative,
                 paired = stat.paired,
                 var.equal=var.equal,
                 mu=mu,
                 conf.level=conf.level)$p.value

    }else{
      pv<-NA
    }
    #p_values<-cbind(p_values,pv)

  }
  else if(s_method=="wilcox-test"){
    #if(nrow(na.omit(cond1_data))>1 &  nrow(na.omit(cond2_data))>1){
    #if(length(na.omit(cond1_data[,i,drop=T]))==length(cond1_data[,i,drop=T]) &  length(na.omit(cond2_data[,i,drop=T]))==length(cond2_data[,i,drop=T])){
    cond1_with_na<- length(na.omit(cond1_data[,i,drop=T]))==length(cond1_data[,i,drop=T])
    cond2_with_na<- suppressWarnings(length(na.omit(cond2_data[,i,drop=T]))==length(cond2_data[,i,drop=T]))
    if(cond1_with_na & cond2_with_na ){
    pv<-wilcox.test(cond1_data[,i,drop=T],
                      cond2_data[,i,drop=T],
                      alternative = alternative,
                      stat.paired = stat.paired,
                      paired = stat.paired,
                      mu=mu,
                      conf.level=conf.level)$p.value
    }else{

      pv<-NA
    }
    #p_values<-cbind(p_values,pv)

  }

  if(silent==F) cat(paste0(" ",s_method," var.equal=",var.equal,"\n"))
  #p_values<-cbind(p_values,pv)
  #s_methods<-rbind(s_methods,s_method)
  #var.equals<-rbind(var.equals,var.equal)

  p_values<-c(p_values,pv)
  s_methods<-c(s_methods,s_method)
  var.equals<-c(var.equals,var.equal)


}

#adjust p value
p_adjust<-p.adjust(p_values,p.adjust.method,length(p_values))



stat_result<-data.frame(p_values=p_values,
                        p_adjust=p_adjust,
                        s_methods=s_methods,
                        var.equals=var.equals)


return(stat_result)

}



#' @title Statistics of clusters(Abundance and Experssion)
#'
#' @description
#' 统计不同分组中各个cluster中marker表达的差异
#'
#' @param combined_data_raw  在Pepline前面部分读入并合并的表达数据
#' @param all_markers        从 “all_markers.csv”读入的dataframe，包含marker的选择信息
#' @param cluster_name       统计所依据的聚类参数名称，例如"PhenoGraph"、"metacluster"
#' @param summerise_method   统计过程中使用的合计方法：有"mean"或者 "median"两种
#' @param groups             groups 从“all_sample.csv”中读入的dataframe，包含实验分组信息（也就是元数据）
#'
#' Heatmap 参数
#' @param major_cond         在groups所提供的各种分组中，选择一种做为统计差异的依据，例如"Tissue_Type"
#' @param heatmap_ctrl       在major_cond中选择一个在Heatmap中显示做为control的组名，例如"A_Tumor"


#'下面5个参数是为后面统计分析中设置的默认值
#'
#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著）
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#'

#' @return a List containing the data which needed in drawing of heatmap, boxplot and volcano chart
#' @export

stat_by_cluster<- function (combined_data_raw,
                            all_markers,
                            cluster_name="metacluster",
                            summerise_method="median", #"mean" or "median" the static method used in summerise
                            groups,

                            #heatmap 参数
                            major_cond,
                            heatmap_ctrl,

                            #全局统计默认参数设置
                            stat.paired=TRUE,
                            stat.method="t-test",
                            set.equal.var="FALSE",
                            conf.level=0.95,
                            p.adjust.method="BH"

                          ){






  cluster_name="metacluster"
  summerise_method="median" #"mean" or "median" the static method used in summerise





cluster_name<-as.character(cluster_name)
ht_markers<-as.character(subset(all_markers,expr_para==1)$markers)
trans_markers<-as.character(subset(all_markers,transform==1)$markers)
used_markers<-union(ht_markers,trans_markers)

cond_name<-major_cond
cond_name<-as.character(cond_name)

  combined_data_df<- tbl_df(combined_data_raw)

  if(summerise_method=="mean"){
    use_method<-mean}else
      if(summerise_method=="median")
      { use_method<-median
      }

  #generating data for total heatmap
  cluster_expr_matrix<-
    combined_data_df %>%
    group_by_(cluster_name) %>%
    summarise_if(is.numeric,use_method,na.rm=TRUE)



    combined_data_tidy<-
    combined_data_df %>%
    group_by_("File_ID",cluster_name) %>%
        summarise_if(is.numeric,use_method,na.rm=TRUE)



  shape_index<-c(49:57,65:90,97:122)  #数字，大写字母，小写字母  shape
  shape_char<-intToUtf8(shape_index,multiple=T)
  cluster_id<-which(colnames(combined_data_tidy)==cluster_name)
  clustern<-max(unique(combined_data_tidy[cluster_id]))

  cluster_expr_trans<-list()
  cluster_expr_ctrl<-list()
  cluster_stat<-list()
  cluster_cell_count<-list()
  cluster_expr_raw<-list()

  for(cluster_i in 1:clustern){
    #cluster_i=3
    cat(paste0("Start to treatwith cluster",cluster_i,"\n"))

    #统计指定cluster上的表达数据
    expr_raw<-
    combined_data_tidy %>%
          dplyr::filter_(paste0(cluster_name,"==",cluster_i)) %>%
          full_join(groups,by="File_ID") %>%
          arrange_("Short_name") %>%
          arrange_(cond_name)


    #统计指定cluster中含有各个文件的细胞数目
    cell_count<-
    combined_data_df %>%
          dplyr::filter_(paste0(cluster_name,"==",cluster_i)) %>%
          count_("File_ID") %>%
          full_join(groups,by="File_ID") %>%
          arrange_("Short_name") %>%
          arrange_(cond_name)

    expr_raw<-ungroup(expr_raw)


    #计算Arcsin ratio

    expr_trans<-expr_raw
    expr_trans[,trans_markers] <- apply(expr_trans[, trans_markers,drop = FALSE],2,simpleAsinh)
    expr_ctrl<-as.data.frame(matrix(0,ncol(expr_trans),nrow=0))
    colnames(expr_ctrl)<-colnames(expr_trans)


    #paired samples
    if(stat.paired){
    for(conditions in t(unique(expr_raw[,cond_name])))
      expr_ctrl<-rbind.data.frame(expr_ctrl,
                                                dplyr::filter_(expr_trans,paste0(cond_name,"==heatmap_ctrl")))
    }else{
      unpaired_control<-dplyr::filter_(expr_trans,paste0(cond_name,"==heatmap_ctrl"))
      single_ctrl<-unpaired_control[1,]
      single_ctrl[,used_markers]<-summarise_if(unpaired_control[,used_markers],is.numeric,funs(use_method))
      expr_ctrl<-expr_trans
      for(row_i in c(1:nrow(expr_ctrl))){
        expr_ctrl[row_i,used_markers]<-single_ctrl[,used_markers]
      }


      #heatmap_expression_ctrl<-
    }

    col_color_index<-as.numeric(as.factor(t(expr_raw[,cond_name,drop=T])))
    cluster_expr_trans[[cluster_i]]<-expr_trans
    cluster_expr_ctrl[[cluster_i]]<-expr_ctrl
    cluster_cell_count[[cluster_i]]<-cell_count
    cluster_expr_raw[[cluster_i]]<-expr_raw



  }


  #Abundance数据汇总
  cat("Summarising cluster abundance data: Cell Count and Percentage\n ")
  cluster_count_summary<-cluster_cell_count[[1]][,1]
  for(cluster_i in c(1:clustern))
  { cluster_count_summary=full_join(cluster_count_summary,cluster_cell_count[[cluster_i]][,c(1,2)],by="File_ID")
  }

  cluster_count_summary[is.na(cluster_count_summary)]<-0
  colnames(cluster_count_summary)<-c("File_ID",paste0("cluster",c(1:clustern)))

  percentage<-function(var){return(var/sum(var)*100)}
  cluster_percent_summary<-cluster_count_summary
  cluster_percent_summary[,-1]<-t(apply(cluster_percent_summary[,-1,drop=F],1,percentage))





  cluster_stat[["metadata"]]<- data.frame( cluster_name=cluster_name,
                                           clustern=clustern,
                                           summerise_method=summerise_method,
                                           #heatmap 参数
                                           major_cond=major_cond,
                                           heatmap_ctrl=heatmap_ctrl,
                                           stat.paired=stat.paired,
                                           stat.method=stat.method,
                                           p.adjust.method=p.adjust.method,
                                           set.equal.var=set.equal.var,
                                           conf.level=conf.level
                                           )







  cluster_stat[["cluster_expr_matrix"]]<-cluster_expr_matrix
  cluster_stat[["cluster_expr_trans"]]<-cluster_expr_trans
  cluster_stat[["cluster_expr_ctrl"]]<-cluster_expr_ctrl
  #cluster_stat[["volcano_data_all"]]<-volcano_data_all
  cluster_stat[["col_color_index"]]<-col_color_index
  cluster_stat[["cluster_cell_count"]]<-cluster_cell_count
  cluster_stat[["ht_markers"]]<-ht_markers
  cluster_stat[["groups"]]<-groups
  cluster_stat[["all_markers"]]<-all_markers
  cluster_stat[["cluster_count_summary"]]<-cluster_count_summary
  cluster_stat[["cluster_percent_summary"]]<-cluster_percent_summary
  cluster_stat[["cluster_expr_raw"]]<-cluster_expr_raw



  return(cluster_stat)
}


#'@title 生成各个cluster的Heatmap以及boxplot
#'@description
#'依据stat_by_cluster生成的list数据，生成
#'1）heatmap1：指定cluster中，各个样本中的marker的表达情况
#'2）heatmap2：指定cluster中，各个样本中的marker与对照的相对强度（Arcsinh raitio）
#'3）boxplots：统计各marker在该cluster中的表达
#'
#'@param cluster_stat  stat_by_cluster()函数生成的list数据
#'@param cluster_id    指定输出cluster的id，例如：如果要输出前五个cluster，cluster_id=c(1:5); 如果要输出全部cluster，保持其默认值NULL即可。
#'@param output_dir    输出数据文件夹名称


#'如使用统计默认参数，以下参数无需设置
#'
#'@param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#'@param conf.level      显著性水平，默认0.95,（即p<0.05为显著）
#'@param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#'@param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#'@param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#'@return none
#' @export


cluster_expr_report<-function(cluster_stat,
                              cluster_id=NULL,
                              output_dir="Cluster_expr_report(heatmap and boxplot)",

                              #如使用全局统计默认参数，以下参数无需设置
                              stat.paired=FALSE,
                              stat.method="t-test",
                              set.equal.var="FALSE",
                              conf.level=0.95,
                              p.adjust.method="BH"){






  col_color_index  <-  cluster_stat[["col_color_index"]]
  #col_color<-brewer.pal(length(unique(col_color_index)),name="Set1")[col_color_index]
  
  col_color<-brewer_color_sets(length(unique(col_color_index)))[col_color_index]
  cluster_expr_trans<- cluster_stat[["cluster_expr_trans"]]
  cluster_expr_ctrl <- cluster_stat[["cluster_expr_ctrl"]]
  cluster_cell_count<- cluster_stat[["cluster_cell_count"]]
  ht_markers        <- cluster_stat[["ht_markers"]]
  summerise_method  <- cluster_stat[["metadata"]]$summerise_method
  cluster_name      <- cluster_stat[["metadata"]]$cluster_name
  clustern          <- cluster_stat[["metadata"]]$clustern
  major_cond        <- cluster_stat[["metadata"]]$major_cond


  if(is.null(stat.paired)) stat.paired        <- cluster_stat[["metadata"]]$stat.paired
  if(is.null(stat.method)) stat.method        <- as.character(cluster_stat[["metadata"]]$stat.method)
  if(is.null(conf.level)) conf.level          <- cluster_stat[["metadata"]]$conf.level
  if(is.null(set.equal.var)) set.equal.var    <- as.character(cluster_stat[["metadata"]]$set.equal.var)
  if(is.null(p.adjust.method)) p.adjust.method    <- as.character(cluster_stat[["metadata"]]$p.adjust.method)

  heatmap_ctrl    <- as.character(cluster_stat[["metadata"]]$heatmap_ctrl)
  if(is.null(cluster_id)){
    cluster_id=c(1:clustern)
  }


  if(stat.method=="auto") {
     stat.method="t.test"
     cat("Warning: Outputing fiures in auto mode, set stat.method to t-test\n" )

  }else { stat.method=sub("-",".",stat.method)}



  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  wdir=getwd()

  cat(paste0("Output to folder: ",wdir,"/",output_dir,"\n"))



  major_cond<-as.character(major_cond)
  shape_index<-c(49:57,65:90,97:122)  #数字，大写字母，小写字母  shape
  shape_char<-intToUtf8(shape_index,multiple=T)



  #cluster_id<-which(colnames(combined_data_raw)==cluster_name)
  #clustern<-max(unique(combined_data_raw[cluster_id]))



  for(cluster_i in cluster_id){

    #cluster_i=1
    cat(paste0("Output  heatmaps,barplots,boxplots of cluster ",cluster_i,"...\n"))

    pdf(file=paste0("./",output_dir,"/",cluster_name,cluster_i,"(Symbol_", shape_char[cluster_i],").pdf"))

    expr_trans<-cluster_expr_trans[[cluster_i]]
    exr_ctrl<-cluster_expr_ctrl[[cluster_i]]

    #画heatmap
    Population_Symbol<-rawToChar(as.raw(cluster_i+64))
    ht_expr_trans<-dplyr::select(expr_trans,ht_markers)
    ht_expr_trans<-as.matrix(t(ht_expr_trans))
    colnames(ht_expr_trans)<-expr_trans$Short_name
    ht_expr_trans<-apply(ht_expr_trans,2,remove_extremum)

    if (!all(is.na(ht_expr_trans))){

      #Page1 Cluster_i 各个样本表达量
      heatmap.2(ht_expr_trans,
                Rowv=T,Colv=F,dendrogram="row",
                scale="none",
                key=T,keysize =1,key.title = paste0("Arcsinh"),
                trace="none",
                col=colorRampPalette(c("black","yellow")),
                srtCol=45,
                margin=c(15,15),
                colCol = col_color,
                ColSideColors = col_color,
                main=paste0("Cluster ",cluster_i," (",shape_char[cluster_i],")"),
                na.color="grey")
    }



    ht_expr_trans_ratio<-expr_trans[,ht_markers]-exr_ctrl[,ht_markers]
    ht_expr_trans_ratio<-as.matrix(t(ht_expr_trans_ratio))
    colnames(ht_expr_trans_ratio)<-expr_trans$Short_name
    ht_expr_trans_ratio<-apply(ht_expr_trans_ratio,2,remove_extremum)

    if (!all(is.na(ht_expr_trans_ratio))){

      #Page1 各个样本表达量
      heatmap.2(ht_expr_trans_ratio,
                Rowv=T,Colv=F,dendrogram="row",
                scale="none",
                key=T,keysize =1,key.title = paste0("Arcsinh Ratio"),
                trace="none",
                col=bluered,
                srtCol=45,
                margin=c(15,15),
                colCol = col_color,
                ColSideColors = col_color,
                main=paste0("Cluster ",cluster_i," (",shape_char[cluster_i],")\n","arcsinh ratio to control"),
                na.color="grey")


    }

    #各组的细胞数量

    mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                     legend.key = element_rect(fill = "white", colour = "white"), #图标
                     legend.background = (element_rect(colour= "white", fill = "white")))

    cell_count<-cluster_cell_count[[cluster_i]]
    cell_count$Short_name<-factor(cell_count$Short_name,levels=cell_count$Short_name,ordered=TRUE)
    cell_count$n<-cell_count$n+1

    count_map<-ggplot(cell_count,aes_string(x="Short_name",y="n",fill=major_cond))+
      geom_bar(stat = "identity")+
      geom_hline(yintercept = 51, linetype="dashed",color = "white", size=0.85)+
      labs(y="number of cells Log(10)",title=paste0("Cell number of population", cluster_i," from each sample"))+
      theme(axis.text= element_text(angle=45,hjust =1,vjust = 0.5))+
      scale_y_log10(limits=c(1,10*max(cell_count$n)),minor_breaks =50)+
      #scale_fill_brewer(palette = "Set1")
      scale_fill_manual(values=brewer_color_sets(length(unique(col_color_index))))

    multiplot(count_map)
    dev.off()


    boxplot_n<-length(ht_markers)
    major_cod_ID<-colnames(expr_trans)==major_cond
    cond_n<-length(unique(expr_trans[,major_cod_ID,drop=T]))
    singlewidth<-cond_n*0.5+1
    singleheight<-4

    boxplot_ncol<-ceiling(sqrt(boxplot_n*4*singleheight/3/singlewidth))
    boxplot_nrow<-ceiling(boxplot_n/boxplot_ncol)

    boxplot_width<-boxplot_ncol*singlewidth
    boxplot_height<-boxplot_nrow*singleheight

    #各个marker的boxplot

    marker_boxplot<-function(boxplot_marker){
      ggplot(expr_trans,aes_string(x=major_cond,y=boxplot_marker,fill=major_cond))+
        geom_boxplot(outlier.shape= NA)+
        #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
        geom_jitter(shape=16, position=position_jitter(0.2))+
        mytheme+
        labs(title=boxplot_marker)+
        #scale_fill_brewer(palette = "Set1")+
        scale_fill_manual(values=brewer_color_sets(length(unique(col_color_index))))+
        labs(y="Expression Level")+
        theme(legend.position = "none")+
        theme(axis.text= element_text(angle=45,hjust =1,vjust = 0.5,size=7))+
        theme(axis.title = element_text(size=8))+
        theme(title =element_text(size=8,face="bold") )+
        stat_compare_means(method = stat.method,paired=stat.paired,ref.group=heatmap_ctrl,aes(label = paste0("p = ", ..p.format..)))
        #stat_compare_means(method = stat.method,paired=stat.paired,ref.group=heatmap_ctrl,aes(label = paste0(" \n", ..p.signif..)),hide.ns=T)
    }

    all_boxplot_list<-lapply(ht_markers,marker_boxplot)

    pdf(file=paste0("./",output_dir,"/",cluster_name,cluster_i,"(Symbol_", shape_char[cluster_i]," Boxplots).pdf"),
        width=boxplot_width,
        height=boxplot_height)
        multiplot(plotlist=all_boxplot_list,cols = boxplot_ncol)
    dev.off()

  }

}




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
#'
#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著），调整该数值可以调整火山图上的横线位置
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#' @return                none
#' @export


draw_expr_volcano<-function(cluster_stat,
                            cond1,
                            cond2,
                            Usecolorset=dif_seq_rainbow,
                            show_adjust_p=FALSE,
                            cluster_id=NULL,
                            xlimit=NULL,
                            ylimit=NULL,
                            dif.level=2,
                            output_dir=paste0("Expr_volcano_plots_",cond1,"_vs_",cond2),

                            #如使用全局统计默认参数，以下参数无需设置
                            stat.paired=NULL,
                            stat.method=NULL,
                            set.equal.var=NULL,
                            conf.level=NULL,
                            p.adjust.method=NULL
                            ){




  clustern          <- cluster_stat[["metadata"]]$clustern
  summerise_method          <- cluster_stat[["metadata"]]$summerise_method
  major_cond             <- cluster_stat[["metadata"]]$major_cond
  all_markers             <- cluster_stat[["all_markers"]]



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
    cluster_id=c(1:clustern)
  }

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


  if(stat.method=="auto") stat.method="t-test or wilcox-test\n(Statistical method is selected automaticly)"

  if(stat.paired){
    if_paired<-"paired"}else{
      if_paired<-"unpaired"}


  cond_name<-major_cond
  cond_name<-as.character(cond_name)
  #画火山图
  volcano_data_all<-list()

  for(cluster_i in cluster_id){

    expr_raw<-cluster_expr_raw[[cluster_i]]

  #使用Acsinh transformed values 计算Pvalue

  cond1_data<-
    expr_raw %>%
    dplyr::filter_(paste0(cond_name,"==cond1"))%>%
    dplyr::select(ht_markers)

  cond2_data<-
    expr_raw %>%
    dplyr::filter_(paste0(cond_name,"==cond2"))%>%
    dplyr::select(ht_markers)


  volcano_cluster_data<-combined_stat(cond1_data,
                                      cond2_data,
                                      alternative="two.sided",
                                      mu=0,
                                      stat.paired=stat.paired,
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

  No_zero<-0.2  #add to prevent divided by zero

  log_ratio=log2((cond2_data_m+No_zero)/(cond1_data_m+No_zero))



  volcano_data<-data.frame(marker=colnames(cond1_data),
                           cluster_i=cluster_i,
                           cluster_shape=shape_char[cluster_i],
                           volcano_cluster_data,
                           log_ratio=t(log_ratio)
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
  FDR_ranks<-factor(c("other(FDR>=10%)","FDR<10%","FDR<5%"),levels=c("other(FDR>=10%)","FDR<10%","FDR<5%"),ordered=TRUE)
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
    labs(y="-Log10(p_values)",title="Summary of False Discovery Rate(FDR)")+
    labs(x=paste0(""))+
    coord_flip()+
    labs(subtitle=paste0("FDRs were separately calculated based on the p values of each cluster "))+
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



  #volcano_data_all$cluster_i<-as.integer(volcano_data_all$cluster_i)
  #clustern<-length(unique(volcano_data_all$cluster_i))
  shape_list<-data.frame(cluster_i=c(1:clustern),
                         shape_index=shape_index[1:clustern])
  volcano_data_all<-full_join(volcano_data_all,shape_list,by="cluster_i")

  marker_n<-length(unique(volcano_data_all$marker))

  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))




   cat(paste0("Output Volcanoplot...\n"))

   x_threshold<-round(log2(dif.level),digits=1)


   if(is.null(xlimit)){
     xuplimit=max(abs(volcano_data_all$log_ratio)+2,na.rm = TRUE)
     xlimit=c(-1*xuplimit,xuplimit)}


   Y_ID<-colnames(volcano_data_all)==p_para

   if(is.null(ylimit)){
     yuplimit=min(max(volcano_data_all[Y_ID]+1.5,na.rm = TRUE),4)
     ylimit=c(0,yuplimit)}

head(volcano_data_all)
  volcanoplot<-ggplot(volcano_data_all)+
    scale_shape_identity() +
    geom_point(aes_string(x="log_ratio",y=p_para,colour="marker",shape="shape_index"),size=4)+
    labs(y="-Log10(p_values)",title="Differential Expression of Markers")+
    labs(x=paste0("Log2_Ratio","\n", cond1, " <- versus -> ", cond2))+
    labs(subtitle=paste0(adjust_label,"p values of arcsinh transformed expression were calculated by ",if_paired," ",stat.method))+
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




  pdf(file=paste0("./",output_dir,"/","Volcano_",cond1," vs ",cond2,".pdf"), width=20, height=15)

  multiplot(volcanoplot,volcano_note,FDR_pValue_plot,layout=layout)
  multiplot(FDR_table_plot)

  dev.off()

  cat(paste0("Output Volcanoplot facet by markers...\n"))

  pdf(file=paste0("./",output_dir,"/","Volcano_",cond1," vs ",cond2,"by markers.pdf"),
      width=20, height=20)


  volcanoplot2<-ggplot(volcano_data_all)+
    scale_shape_identity() +
    geom_point(aes_string(x="log_ratio",y=p_para,colour="marker",shape="shape_index"),size=4)+
    #geom_point(aes(x=log_ratio,y=p_adjust,colour=marker_name,shape=shape_index),size=4)+
    labs(y="-Log10(p_values)",title="Differential Expression of Markers")+
    labs(x=paste0("Log2_Ratio","\n", cond1, " <- versus -> ", cond2))+
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
    #geom_point(aes(x=log_ratio,y=p_adjust,colour=marker_name,shape=shape_index),size=4)+
    geom_point(aes_string(x="log_ratio",y=p_para,colour="marker",shape="shape_index"),size=4)+
    labs(y="-Log10(p_values)",title="Differential Expression of Markers")+
    labs(x=paste0("Log2_Ratio","\n", cond1, " <- versus -> ", cond2))+
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






#' @title Arcsinh 转化函数

#' @description：
#' 这是一个数据转化函数，由Cytofkit::cytofAsinh() 函数修改而来，根据需要进行了简化，主要去掉了原函数关于扣除背景和随机化的部分。
#'
#' @param value the data need to be transformed,a vector of numeric values.
#' @param cofactor Cofactor for asinh transformation, default 5 for mass cytometry data.
#' @return a vector of numeric values containing transformed data
#' @export

simpleAsinh <- function(value, cofactor = 5) {
  value <- value / cofactor
  value <- asinh(value)
  return(value)
}




#' @title 去除极值
#' @description：
#' 去除极值，将大于0.95分位数的数值赋值为0.95分位数，小于0.05分为数的数值赋值为0.05分位数
#' @param value 待处理的数据
#' @return a vector of numeric values containing treated data
#' @export
#'
remove_extremum<-function(value){
  quantile_95<-quantile(value,0.95,na.rm=TRUE)
  quantile_05<-quantile(value,0.05,na.rm=TRUE)
  value[value>quantile_95]<-quantile_95
  value[value<quantile_05]<-quantile_05
  return(value)
}






#' @title 简易表格
#' @description：
#' 将输入dataframe的内容化成表格
#' @param data 待输入的dataframe
#' @param show_row_name TRUE/FALSE, 决定是否在输出的表格中显示行名
#' @param show_col_name TRUE/FALSE, 决定是否在输出的表格中显示列名
#' @return 返回一个ggplot对象
#' @export

ggtable<-function(data,show_row_name=TRUE,show_col_name=TRUE){

        shape_table_chr<-apply(data,2,as.character)
        if(show_col_name)
        shape_table_chr<-rbind(as.character(colnames(data)),
                               shape_table_chr)
        if(show_row_name){
        shape_table_chr<-cbind(as.character(row.names(data)),
                               shape_table_chr)

        }
        table_data<-as.data.frame(matrix(ncol=3,nrow=0))
        colnames(table_data)<-c("text","x","y")

        table_row_number=1
        row_height=1
        col_width=1

        for(row_i in c(1:nrow(shape_table_chr))){

          for(col_i in c(1:ncol(shape_table_chr))){

            table_data[table_row_number,]$text=shape_table_chr[row_i,col_i]
            table_data[table_row_number,]$y=row_height*row_i
            table_data[table_row_number,]$x=col_width*col_i

            table_row_number=table_row_number+1
          }
        }

        table_theme<-theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.2), #坐标系及坐标轴
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border= element_blank(),
                           axis.title=element_blank(),
                           axis.ticks=element_blank(),
                           axis.text = element_blank(),
                           legend.key = element_rect(fill = "white", colour = "white"), #图标
                           axis.line=element_line(colour='white',linetype = 1,size=0)
                           #plot.margin=unit(c(2,2,2,2), unit = "cm")
        )
        ggplot()+
          geom_text(data = table_data,label=table_data$text,aes(x=x,y=y))+
          scale_y_reverse()+
          table_theme+
          expand_limits(x=c(min(table_data$x-0.5),max(table_data$x+0.5)),y=c(min(table_data$y-0.5),max(table_data$y+0.5)))
       }
