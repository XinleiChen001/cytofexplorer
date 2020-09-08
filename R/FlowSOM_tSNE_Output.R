


#修改说明：
#2019——10-28  增加tsne_heatmap按照major_cond分组显示的功能
#2019_11_10   combined_data_sampled开始使用Raw Data,对density,tsne heatmap的代码进行适应性调整； 
#             tsne heatmap增加raw,0_to_max,min_to_max等bar展示方法
#             将phenograph_heatmap挪进此步骤
#2019_11_11   为tsne散点图增加等高线背景
#2019-12-20   修复分组sample的bug
#2020-01-10   为draw_tsne_figs增加cluster_barplot_data参数
#2020-01-20   增加cluster_marker_preview函数
#2020-04-18   修复cluster_marker_perview bug: 用clustername替代PhenoGraph
#2020-4-18    修复长宽比问题bug
#2020-5-9     在draw_tsne_figs和tsne heatmap两个函数中增加网路的功能


#' @title dif_seq_rainbow
#'
#' @description
#' 生成一系列差异较大的颜色，适合于在散点图上显示不同population的细胞

#'
#' @param n                需要产生的颜色种类数量
#' @return                 返回指定数量的颜色序列
#' @export

dif_seq_rainbow<-function(n){

  library("colorspace")


  div=11
  depth= n %/% div+2
  if (n %% div!=0) depth=depth+1

  base_color<-c("red","#FF6B00","gold","greenyellow","green","cyan","#007FFF","blue","#4A00E9","darkviolet","#C90069")
  color_matrix<-matrix(0,ncol=depth,nrow=0)

  for(col_i in base_color){

    seq_col<-colorRampPalette(c("white",col_i,"black"))(depth)

    color_matrix<-rbind(color_matrix,seq_col)
  }

  color_matrix<-color_matrix[,c(-1,-1*depth)]

  color_matrix_order<-data.frame(n=c(1:(depth-2)),t(color_matrix))
  color_matrix_order$n<-abs(color_matrix_order$n-mean(color_matrix_order$n))
  color_matrix_order2<-color_matrix_order[order(color_matrix_order[,"n"]),]
  color_matrix_order3<-t(color_matrix_order2[,-1])

  set.seed(456)
  rdm<-function(var){sample(var,length(var))}

  #rdm(1)

  color_matrix_order4<-apply(color_matrix_order3,2,rdm)

  color_vector<-as.vector(color_matrix_order4)[c(1:n)]
  #if(n>30) cat("NOTE:You are using dif_seq_hcl to produce more than 30 colors, some of them may be difficult to distinguish\n")
  return(color_vector)
}






#' @title 生成颜色系列
#'
#' @description
#' RColorBrewer Set1,Set2,Set3合在一起使用，生成包含29种差异较大的颜色的系列

#'
#' @param n                需要产生的颜色种类数量
#' @return                 返回指定数量的颜色序列
#' @export


brewer_color_sets<-function(n)
  
{
  library(RColorBrewer)
  
  brewer_sets<-c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"))[-c(6,19)]  #剔除两个黄色
  
  #colorset<-as.data.frame(matrix(ncol=1,nrow = n))
  
  if(n>0) return(brewer_sets[1:n])
  
}





#' @title Elbow Test
#'
#' @description
#' Phenograph 的k值是指在计算k 近邻网络(knn)的所取得k值，它的意义是网络中每个细胞的最近邻居数。K值越大，最后得到的cluster越少，反之越大。
#' 因此过大的k值会引起“underclustering”，导致cluster偏少，过小的k值则会导致“overclustering”,导致cluster偏高。
#' 选择合适的k值就是在两者之间寻求平衡。比较直观的方法就是在k-cluster图中，找到曲线的 “elbow piont”。
#'
#'
#' @param x                需要进行测试的表达数据矩阵
#' @param k_from           自然数，测试中k的取值范围的起始点
#' @param k_to             自然数，测试中k的取值范围的结束点
#' @param k_step           自然数，k取值的步长
#' @return 返回值一个矩阵，记录每个k值对应的cluster数目
#' @export


PG_elbow<-function(x,k_from=5,k_to=100,k_step=5){

  test_range<-seq(from=k_from,to=k_to,by=k_step)
  PhenoGraph_elbow<-c(0,0)

  for(k in test_range){

    cluster_PhenoGraph <-as.numeric(membership(Rphenograph(data = x,k)))
    if(k==5){PhenoGraph_elbow<-c(k,max(cluster_PhenoGraph))
    }else{PhenoGraph_elbow<-rbind(PhenoGraph_elbow,c(k,max(cluster_PhenoGraph)))
    }
  }

  colnames(PhenoGraph_elbow)<-c("PhenoGraph k","cluser number")
  plot(PhenoGraph_elbow,type="b")
  return(PhenoGraph_elbow)
}







#' @title metacluster_multi
#'
#' @description
#'利用FlowSOM或其他方法进行聚类以后，可以对产生的cluster进一步聚类，这些新产生的类被称为“metacluster”。
#'metacluster_multi可以利用FlowSOM自带的metacluster函数，以及PhenoGraph算法产生metacluster。
#'#'
#' @param data             包含所有cluster的表达矩阵，dataframe或者矩阵格式
#' @param method           要使用的聚类方法名称，取值有“PhenoGraph”，“metaClustering_consensus”,“metaClustering_hclust”, “metaClustering_kmeans”
#' @param nCluster         使用“metaClustering_consensus”和“metaClustering_kmeans”两种方法时，用其指定产生metaCluster的数量
#' @param max              使用“metaClustering_hclust”时，确定一个metacluster的上限，函数将在此下寻找合适的metacluster数量
#' @param k                使用“PhenoGraph"时，设置knn网络的k值，默认值为15
#'
#' @return                 一个数组（Numeric array）,包含对应每个cluster的metacluster数值
#' @export


metacluster_multi<-function(data,
                            method="metaClustering_hclust",
                            nCluster=NULL,
                            max=50,
                            k=15){
  cat(paste0("Metaclustering with ",method))
  if (method!="PhenoGraph"){
    metaClusters <-suppressMessages(MetaClustering(data=data,
                                                   method=method,
                                                   nClus=nCluster,
                                                   max=max))
  }else{
    metaClusters <-as.numeric(membership(Rphenograph(data = data,k=k)))
  }

  return(metaClusters)
}



#' @title equal_sample
#'
#' @description
#' 从各个Meta cluster或文件中抽取定量细胞进行可视化
#' 可视化前的Sampling有两个功能，第一，减少后续用于降维分析的细胞总数，一般推荐在10万以内，以保证降维效果。
#' 第二， 解决原始数据中文件或者cluster存在细胞数目不均衡的情况。增加对小文件或cluster的可见性。

#' @param x                需要进行Sampleing表达数据矩阵,除了表达数据，还包括两个通道：“File_ID”,以及聚类产生的cluster通道。
#' @param sample_type      Sample的方式，#sample_type="by_cluster",从每个cluster抽取相同细胞；sample_type="by_file"，从每个样本抽取相同细胞
#' @param sample_method    默认"ceil",   #两种取值："ceil", "all"，取"all"时，采取所有细胞
#' @param sample_num         每个样本或者每个cluster抽取的细胞数量；
#' @param cluster_name     cluster通道的名称
#' @return 返回值一个矩阵，包含所有sample取得的表达数据
#' @export

equal_sample<-function( x,
                        sample_type="by_file", #sample_type="by_cluster",从每个cluster抽取相同细胞；sample_type="by_file"，从每个样本(文件)抽取相同细胞; "by_cond",从每个major condition抽取相同细胞
                        sample_method="ceil",   #四种取值："ceil", "all"
                        sample_num=1000, #每个样本或者每个cluster抽取的细胞数量；
                        cluster_name="metaCluster",
                        groups=NULL,
                        sample_cond=NULL
)
{
  input_col_name<-colnames(x)  #记录x的原始列明
#x=combined_data_analysed
x$row_name<-rownames(x)
if(sample_type=="by_cond") {
  #sample_cond="Tissue_Type"
  sample_per_cond<-unique(groups[,as.character(sample_cond)])
  groups$File_ID<-as.character(groups$File_ID)
  x$File_ID<-as.character(x$File_ID)
  x<-dplyr::left_join(x,groups,by="File_ID")
  
}


if(sample_type=="by_cluster") sample_col=cluster_name
if(sample_type=="by_file") sample_col="File_ID"
if(sample_type=="by_cond") sample_col=sample_cond

combined_data_sampled<-as.data.frame(matrix(nrow=0,ncol=ncol(x)))
for(meta_i in unique(x[,sample_col])){
  #meta_i=unique(x[,sample_col][1])
  cell_n=length(which(x[,sample_col]==meta_i))
  if (sample_method=="ceil"){
    cell_n=min(cell_n,sample_num)
  }
  
  data_sampled<-x%>%
    #dplyr::filter_(paste0(sample_col,"==meta_i"))%>%
    dplyr::filter_at(vars(matches(sample_col)),all_vars(.==meta_i)) %>%
    sample_n(cell_n)
  data_sampled
  combined_data_sampled<-rbind(data_sampled,combined_data_sampled)
}



rownames(combined_data_sampled)<-combined_data_sampled$row_name
combined_data_sampled$row_name=NULL

clustern<-length(unique(x[,sample_col]))

#cat(paste0("Total sampled events：",nrow(combined_data_sampled)))

layout(matrix(c(1,3,2,4),ncol=2))
smr_data_total<-x %>% count_(cluster_name)   #smr is short for summerise

plot1<-barplot(smr_data_total$n,
               names.arg=smr_data_total[,1,drop=T],
               main="Events of clusters(Total)")

smr_data_sampled<-combined_data_sampled %>% count_(cluster_name)

plot2<-barplot(smr_data_sampled$n,
               names.arg=smr_data_sampled[,1,drop=T],
               main=paste0("Events of clusters (Sampled ",sample_type,")"))


smr_data_total<-x %>% count_("File_ID")

plot3<-barplot(smr_data_total$n,
               #names.arg=smr_data_sampled[,1,drop=T],
               main="Events of Files (Total)")

smr_data_sampled<-combined_data_sampled %>% count_("File_ID")

plot4<-barplot(smr_data_sampled$n,
               #names.arg=smr_data_sampled[,1,drop=T],
               #axis(side=1,lwd=0,las=1),
               main=paste0("Events of Files (Sampled ",sample_type,")"))


layout(matrix(c(1),ncol=1))

return(combined_data_sampled[,input_col_name])

}





#' @title draw_tsne_figs
#' @description
#' 自动生成系列tSNE图；
#'
#' @param combined_data_plot     数据框或者矩阵，附带有(meta)cluster聚类、降维结果(t_sne_1,t_sne_2)的表达矩阵
#' @param edge_data              画网络图的时候，存储网络连接数据
#' @param groups                 实验组的Metadata信息,其中必需包含列名：”Short_name“
#' @param cluster_color          tSNE图的亚群的配色方案，默认dif_seq_rainbow
#' @param cluster_name           选择groups中一个列名做为展现差异的major_cond
#' @param output_dir             输出数据文件夹名称
#' @param major_cond             选择groups中一个列名做为展现差异的major_cond
#' @param cluster_barplot_data   包含cluster通道信息的Vector，用来生成cluster的丰度直方图
#' @param reduction_dm1          combined_data_plot降维产生的维度，默认"tsne_1"
#' @param reduction_dm2          combined_data_plot降维产生的维度，默认"tsne_2"
#' @param dot_size               散点图点的大小，默认是4
#' @param cluster_labels_size    控制cluster label的大小，默认30
#' @param show_axis              是否显示坐标轴，默认是TRUE
#' @param show_cluster_labels    决定是否显示cluster label，TRUE 显示Label，FALSE 不显示
#' @param show_contour           布尔型变量，决定是否显示等高线背景
#' @param contour_line_size      设置等高线的宽度，默认0.5
#' @param contour_bins           设置等高线的层数，默认25
#' @param contour_colour         设置等高线的颜色，默认"grey"
#' @param group_seq              设置各组顺序，格式实例： group_seq=c("PBMC","Biopsy")，注：括号内为各组名称；
#' @param groups_to_show         指定要显示的condition（major_cond中的一个或者多个），如不设置直接统计全部；
#' @param edge_line_size         连线的宽度      
#' @param edge_layer             连线的图层位置，"back"和"front"，决定连线图层在散点图的后面或前面
#' @param output_format          输出数据文件格式："tiff"和“pdf”,默认"pdf"

#' @export


draw_tsne_figs<-function(
  combined_data_plot,
  edge_data=NULL,
  groups,
  cluster_color=dif_seq_rainbow,
  cluster_name="metacluster",
  output_dir="cluster_tsne_plots",
  major_cond="Tissue_Type",
  cluster_barplot_data=NULL,
  reduction_dm1="tsne_1",
  reduction_dm2="tsne_2",
  dot_size=4,
  cluster_labels_size=15,
  show_axis=TRUE,
  show_cluster_labels=TRUE,
  show_contour=FALSE,
  edge_colour="grey",
  contour_line_size=0.5,
  contour_bins=25,
  contour_colour="grey",
  group_seq=NULL,
  groups_to_show=NULL,
  edge_line_size=3,
  edge_layer="front",
  output_format="pdf"){
  
  
  if(0){
    combined_data_plot<-combined_rawdata_plot
    cluster_color=dif_seq_rainbow
    cluster_name="metacluster"
    output_dir="cluster_tsne_plots"
    major_cond="Tissue_Type"
    reduction_dm1="tsne_1"
    reduction_dm2="tsne_2"
    cluster_barplot_data=NULL
    dot_size=4
    text_size=30
    show_contour=FALSE
    contour_line_size=0.5
    contour_bins=25
    contour_colour="grey"
    group_seq=NULL
    output_format="pdf"
    
  }  
  
  
  
  library(colorRamps)
  library(Rmisc)
  
  if (!is.null(group_seq)){
    groups[,major_cond]<-factor(groups[,major_cond,drop=T],
                                levels=group_seq,
                                ordered=T)
  }
  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }
  
  
  if(is.null(groups_to_show)){
    groups_to_show=unique(groups[,major_cond,drop=T])  
  } 
  
  
  rdm_name<-sub("[0-9]","",reduction_dm1)
  combined_data_plot$File_ID<-as.character(combined_data_plot$File_ID)
  combined_data_plot<-full_join(combined_data_plot,groups,by="File_ID")
  combined_data_plot2<-combined_data_plot
  #%>%
  #dplyr::filter(Merge_in_fig=="Yes")
  
  #str(combined_data_plot2)
  cluster_index<-colnames(combined_data_plot)==cluster_name
  combined_data_plot[,cluster_index]<-as.factor(combined_data_plot[,cluster_index])
  
  major_cond_index<-colnames(combined_data_plot)==major_cond
  major_cond_n<-length(unique(combined_data_plot[,major_cond_index,drop=TRUE]))
  
  
  
  combined_data_plot<-combined_data_plot %>%
    group_by_at(cluster_name)%>%
    dplyr::filter_at(vars(matches(major_cond)),all_vars(.%in% groups_to_show))
  

  centers<-eval(parse(text=paste0("summarise(combined_data_plot,",reduction_dm1,"=median(",reduction_dm1,"),",reduction_dm2,"=median(",reduction_dm2,"))")))
  
  
  for (file_id in groups$File_ID)
  {
    file_center<-centers
    file_center$File_ID<-file_id
    
    if (file_id==groups[1,"File_ID"]){
      file_centers<-file_center
    }else{
      file_centers<-rbind(file_centers,file_center)
    }
  }
  
  file_centers<-full_join(file_centers,groups,by="File_ID")
  #str(file_centers)
  
  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))
  head(combined_data_plot)
  
  combined_data_plot<-combined_data_plot %>%
    dplyr::arrange(Short_name) %>%
    dplyr::arrange_(major_cond)%>%
    as.data.frame()
  
  
  #生成cluster abundance的barplot
  if(!is.null(cluster_barplot_data)){
  cluster_data<-data.frame(cluster=cluster_barplot_data)
  cluster_abandance_stat<-cluster_data %>%
                            dplyr::group_by(cluster)%>%
                               dplyr::summarise(num=n())
  
  cluster_num<-length(unique(combined_data_plot[,cluster_name]))
  if(nrow(cluster_abandance_stat)!=cluster_num){
    message("Warning:cluster number in cluster_barplot_data is not equal to combined_data_plot.\n")
  }
      
  cluster_abandance_stat$percentage=cluster_abandance_stat$num/sum(cluster_abandance_stat$num)*100
  cluster_abandance_stat$fill=cluster_color(nrow(cluster_abandance_stat))[cluster_abandance_stat$cluster]
  cluster_abandance_stat$cluster<-factor(cluster_abandance_stat$cluster,levels=cluster_abandance_stat$cluster,ordered = T)
  cluster_abundance_barplot<-ggplot(cluster_abandance_stat)+
                               geom_bar(aes(x=cluster,y=percentage,fill=cluster),stat = "identity")+
                               scale_fill_manual(values = cluster_abandance_stat$fill)+
                                 mytheme+
                                 coord_flip()+#横向
                                 scale_x_discrete(limits=rev(levels(cluster_abandance_stat$cluster)))
  
  if(output_format=="pdf"){
    pdf(file=paste0("./",output_dir,"/",rdm_name,"_",cluster_name," barplot.pdf"),width = 300/72,height = 20/72*nrow(cluster_abandance_stat))}else if(output_format=="tiff"){
      tiff(filename = paste0("./",output_dir,"/","tSNE-",cluster_name," barplot.tiff"),width=1300,height=1200)} else
      { message(paste0("Error-Output format ",output_format," is not supported\n"))
        return(NULL)}
  
  multiplot(cluster_abundance_barplot)
  dev.off()
  }
  draw_single_tsne_fig<-function(dotcolor_para=NULL,
                                 dotcolor_n=length(unique(combined_data_plot[,dotcolor_para])),
                                 facet_para=NULL,
                                 facet_ncol=NULL){
  

            #产生facet图背景的数据：
            #登高线图
            combined_data_density<-combined_data_plot
            
            if((show_contour==TRUE)& (!is.null(facet_para))){
              combined_data_density<-as.data.frame(matrix(nrow=0,ncol=ncol(combined_data_plot)))
              colnames(combined_data_density)<-colnames(combined_data_plot)
              for(major_cond_i in unique(combined_data_plot[,facet_para,drop=T])){
                combined_data_density_i<-combined_data_plot
                combined_data_density_i[,facet_para]<-major_cond_i
                combined_data_density<-rbind(combined_data_density,
                                             combined_data_density_i)   
              }
            }
            #edge
            combined_edge_data<-edge_data
            if((!is.null(edge_data))& (!is.null(facet_para))){
                      combined_edge_data<-as.data.frame(matrix(nrow=0,ncol=ncol(edge_data)))
                      colnames(combined_edge_data)<-colnames(edge_data)
                      for(major_cond_i in unique(combined_data_plot[,facet_para,drop=T])){
                        combined_edge_data_i<-edge_data
                        combined_edge_data_i[,facet_para]<-major_cond_i
                        combined_edge_data<-rbind(combined_edge_data,
                                                   combined_edge_data_i)   
                      }
                      
            }
            


            plot1<-ggplot(data=combined_data_plot,aes_string(x=reduction_dm1,y=reduction_dm2))
            
            if(show_contour==TRUE){
              
              plot1<-plot1+ geom_density_2d(data=combined_data_density,colour=contour_colour,size=contour_line_size,bins=contour_bins)+
                scale_x_continuous(limits=c(min(combined_data_plot[,reduction_dm1])*1.1,max(combined_data_plot[,reduction_dm1])*1.1))+
                scale_y_continuous(limits=c(min(combined_data_plot[,reduction_dm2])*1.1,max(combined_data_plot[,reduction_dm2])*1.1))
            }
            
            if((!is.null(edge_data)) & (edge_layer=="back")){
              plot1<-plot1+geom_line(data=combined_edge_data,aes_string(x=reduction_dm1,y=reduction_dm2,group="edge_id"),color=edge_colour,size=edge_line_size)
            }
          
            plot1<-plot1+
              geom_point(aes_string(colour=dotcolor_para),size=dot_size)+
              mytheme+
              theme(text=element_text(size=30),legend.title=element_blank())+
              #scale_color_manual(values = primary.colors(nrow(centers)+1)[-1])#facet_wrap(~File_ID)
              scale_color_manual(values = cluster_color(dotcolor_n))
            
              
            
            if((!is.null(edge_data)) & (edge_layer=="front")){
              plot1<-plot1+geom_line(data=combined_edge_data,aes_string(x=reduction_dm1,y=reduction_dm2,group="edge_id"),color=edge_colour,size=edge_line_size)
            }
            
            if(show_cluster_labels==T){
              plot1<-plot1+geom_text(data=file_centers,aes_string(x=reduction_dm1,y=reduction_dm2),label=file_centers[,1,drop=TRUE],colour="black",size=cluster_labels_size)
            }
            
            if(!show_axis){
              
              plot1<-plot1+
                     theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.2))+
                     theme(axis.line=element_blank(),
                           axis.ticks = element_blank(),
                           axis.title = element_blank(),
                           axis.text=element_blank())}  
              
            if(!is.null(facet_para)){
               plot1<-plot1+facet_wrap(facets=facet_para,ncol=facet_ncol)}

               
            
             plot1
            # plot1<-plot1+ geom_density_2d(aes_string(x=reduction_dm1,y=reduction_dm2))
            
  }        
  
  #1生成合并tSNE-Phenograph图像
  
  if(output_format=="pdf"){
    pdf(file=paste0("./",output_dir,"/",rdm_name,"_",cluster_name," merge.pdf"),width = 1300/72,height = 1200/72)}else if(output_format=="tiff"){
      tiff(filename = paste0("./",output_dir,"/","tSNE-",cluster_name," merge.tiff"),width=1300,height=1200)} else
      { message(paste0("Error-Output format ",output_format," is not supported\n"))
        return(NULL)}
  
      plot1<-draw_single_tsne_fig(dotcolor_para=cluster_name)
  
      multiplot(plot1)
      dev.off()
  cat(paste0("Exported: ",rdm_name,"_",cluster_name," merge.",output_format,"\n"))
  
  
  
  
  #2生成合并tSNE-文件图像

  if(output_format=="pdf"){
    pdf(file=paste0("./",output_dir,"/",rdm_name,"_","Files"," merge.pdf"),width = 1300/72*(1+nrow(groups)/110),height = 1200/72)}else if(output_format=="tiff"){
    tiff(filename = paste0("./",output_dir,"/","tSNE-","Files"," merge.tiff"),width=1300*(1+nrow(groups)/110),height=1200)} else{
      message(paste0("Error-Output format ",output_format," is not supported\n"))
    return(NULL)}
  
  plot1<-draw_single_tsne_fig(dotcolor_para="Short_name")
  
  multiplot(plot1)
  dev.off()
  cat(paste0("Exported: ",rdm_name,"_","Files"," merge.",output_format,"\n"))
  
  #3生成单个文件的tSNE-Phenograph图像
  nfacet_col=ceiling(sqrt(nrow(groups)))
  nfacet_row=ceiling(nrow(groups)/nfacet_col)
  if(output_format=="pdf"){
    pdf(file=paste0("./",output_dir,"/",rdm_name,"_",cluster_name,"(sample name).pdf"),width = 1200*nfacet_col/72,height = 1200*nfacet_row/72)} else if(output_format=="tiff"){
    tiff(filename = paste0("./",output_dir,"/","tSNE-",cluster_name,"(sample name).tiff"),width = 1200*nfacet_col,height = 1200*nfacet_row/72)}else{
      message(paste0("Error-Output format ",output_format," is not supported\n"))
      return(NULL)}
  
  plot1<-draw_single_tsne_fig(dotcolor_para=cluster_name,facet_para = "Short_name",facet_ncol = ceiling(sqrt(nrow(groups))))
  
  
  multiplot(plot1)
  dev.off()
  
  cat(paste0("Exported: ",rdm_name,"_",cluster_name,"(sample name).",output_format,"\n"))
  
  
  
  #4#生成单个分组的tSNE-Phenograph图像
  
  
  if(output_format=="pdf"){
    pdf(file=paste0("./",output_dir,"/",rdm_name,"_",cluster_name,"(",major_cond,").pdf"),width=1200*major_cond_n/72,height=1200/72)}else if(output_format=="tiff"){
    tiff(filename = paste0("./",output_dir,"/",rdm_name,"_",cluster_name,"(",major_cond,").tiff"),width=1200*major_cond_n,height=1200)}  else{
      message(paste0("Error-Output format ",output_format," is not supported\n"))
    return(NULL)}
  
  plot1<-draw_single_tsne_fig(dotcolor_para=cluster_name,facet_para = major_cond,facet_ncol = major_cond_n)
  
  multiplot(plot1)
  dev.off()
  cat(paste0("Exported: ",rdm_name," colored by ",major_cond,".",output_format,"\n"))
  
  
}





#' @title cluster_marker_preview
#' @description
#' 预览各个cluster中marker的表达，以heatmap和直方图的形式展示（密度图）；
#'
#' @param xdata                  数据框或者矩阵，附带有(meta)cluster表达矩阵
#' @param groups                 实验组的Metadata信息,其中必需包含列名：”Short_name“
#' @param cluster_color          tSNE图的亚群的配色方案，默认dif_seq_rainbow
#' @param cluster_name           用过分析的cluster名称
#' @param output_dir             输出数据文件夹名称
#' @param major_cond             选择groups中一个列名做为展现差异的major_cond
#' @param cluster_id             指定显示的cluster
#' @param xlim                   指定x轴范围，例如：c(0,5)，只在free_x=F的时候有效
#' @param free_x                 True或者FALSE，设置是否每张图X-轴取值范围是否自动，建议设置为F
#' @param trans_method           数据转化方式，有四种："CytofAsinh"，"simpleAsinh"；"0_to_Max"，所有Marker的信号强度除以最大值，线性转换，通过除以各通道信号的最大值把数值scale到0~1；"Min_to_Max"，线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异
#' @export

cluster_marker_preview<-function(rawdata,
                                 groups,
                                 all_markers,
                                 cluster_color=dif_seq_rainbow,
                                 expression_color=rev(brewer.pal(n = 7, name ="RdYlBu")),
                                 cluster_name="metacluster",
                                 output_dir="cluster_marker_preview",
                                 major_cond="Tissue_Type",
                                 summerise_method="median",
                                 Rowv=T,
                                 Colv=T,
                                 cluster_id=NULL,
                                 xlim=NULL,
                                 free_x=F,
                                 colored_by="cluster", #"expression"
                                 
                                 trans_method="simpleAsinh"
                                 )

{
  library(colorRamps)
  library(Rmisc)
  library(reshape2)
  #transdata=combined_transdata_analysed

  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }
  
  
  if(trans_method=="simpleAsinh") trans_fun=simpleAsinh
  if(trans_method=="0_to_Max")   trans_fun=normalize_Zero_Max
  if(trans_method=="Min_to_Max") trans_fun=normalize_min_max
  
  if(summerise_method=="mean"){
    use_method<-mean}else
      if(summerise_method=="median")
      { use_method<-median
      }
  #all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
  transform_id=which(all_markers$transform==1)
  transdata<-rawdata
  transdata[,transform_id] <- apply(transdata[,transform_id,drop=F],MARGIN=2,trans_fun) #MAGIN=2 按照列进行Normalize，MAGIN=1按照行进行Normalize
  
  
  transdata$File_ID<-as.character(transdata$File_ID)
  transdata<-dplyr::full_join(transdata,groups,by="File_ID")
  density_plot_marker_id<-which(all_markers[,"density_plot"]==1)
  density_plot_marker<-as.character(all_markers[density_plot_marker_id,]$markers)
  #summerise
  heatmap_data <-transdata%>% 
                    #dplyr::filter_at(vars(one_of(cluster_name)),all_vars(.%in% groups_to_show))%>%
                    group_by_at(c(cluster_name)) %>%
                    summarise_if(is.numeric,use_method,na.rm=TRUE)%>%
                    select_at(vars(one_of(density_plot_marker,cluster_name)))%>%
                    data.frame()

  
  
  #生成树形图及排序
  row.names(heatmap_data)<-as.character(heatmap_data[,cluster_name])
      
  #输出heatmap
  width_fig=ncol(heatmap_data)*0.5+4
  height_fig=nrow(heatmap_data)*0.2+4
      
  library(pheatmap)
   pheatmap(mat=heatmap_data[,density_plot_marker],
           cluster_rows=Rowv,
           cluster_cols=Colv,
           scale = "none",
           col=colorRampPalette(expression_color)(100),
           cellwidth = 10,
           cellheight = 10,
           fontsize = 6,
           filename = paste0("./",output_dir,"/","heatmap",".pdf")
           )
  
      
      
    #row.names(raw_data)<-paste0("cluster",raw_data[,cluster_name])
    #heatmap_data<-heatmap_data[,density_plot_marker]
    model_cluster <- hclust(dist(heatmap_data[,density_plot_marker]), "complete")
    model_marker  <- hclust(dist(t(heatmap_data[,density_plot_marker])), "complete")
    dhc_cluster   <- as.dendrogram(model_cluster)
    dhc_marker    <- as.dendrogram(model_marker)
    melt_htdata<-melt(heatmap_data,id.vars = cluster_name,measure.vars=density_plot_marker,variable.name="Markers",value.name = "Expression")
    melt_htdata[,cluster_name]<-as.character(melt_htdata[,cluster_name])
    melt_htdata[,"Markers"]<-as.character(melt_htdata[,"Markers"])
    
    #melt_htdata$metacluster<-factor(melt_htdata$metacluster,levels =paste0(c(nrow(heatmap_data):1)),ordered = TRUE )
    str(melt_htdata)
    if (is.null(cluster_id)){
     cluster_id=unique(transdata[,cluster_name])} 
 
  transdata<-transdata%>%
             dplyr::filter_at(vars(matches(cluster_name)),all_vars(.%in% cluster_id)) %>%
             dplyr::select(one_of(c(cluster_name,density_plot_marker)))
             transdata[,cluster_name]<-as.character(transdata[,cluster_name])

  #rearrange clusters and markers
  plot_data<-tidyr::gather_(transdata,key="Markers",value="exprs",density_plot_marker)
  plot_data<-dplyr::full_join(plot_data,melt_htdata,by=c(cluster_name,"Markers"))
  
  #排序
  if(Rowv==T){
    plot_data[,cluster_name]<-factor(plot_data[,cluster_name],levels =labels(dhc_cluster),ordered = TRUE )
    melt_htdata[,cluster_name]<-factor(melt_htdata[,cluster_name],levels =labels(dhc_cluster),ordered = TRUE )}
  if(Colv==T){
    plot_data$Markers    <-factor(plot_data$Markers,levels =labels(dhc_marker),ordered = TRUE)
    melt_htdata[,"Markers"]<-factor(melt_htdata[,"Markers"],levels =labels(dhc_marker),ordered = TRUE )
    }
  
  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   #legend.key = element_rect( colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))
  head(plot_data)
  

  #2生成density-plot
  
  col_plot<-length(density_plot_marker)*3.5+3.5
  row_plot<-length(unique(transdata[,cluster_name]))*3.5
  
  pdf(file=paste0("./",output_dir,"/","Density_plots",".pdf"),width = col_plot,height = row_plot)
  
  plot1<-ggplot()
  
  if(  1  ){plot1<- plot1+
            geom_rect(data = melt_htdata,aes(fill=Expression),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha=0.8)
  }
  
  if(colored_by=="expression"){
           plot1<-plot1+
                  geom_density(data=plot_data,aes_string(x= "exprs",fill="Expression"),adjust =1)+
                  mytheme+
                  scale_fill_gradientn(colours =colorRampPalette(expression_color)(100))
                  #scale_fill_manual(values = cluster_color(length(unique(xdata3[,cluster_name]))))#+
                  #facet_wrap(as.formula(paste0(cluster_name," ~Markers") ),scales = "free_y",ncol = length(density_plot_marker))
          
  }else if(colored_by=="cluster"){
          plot1<-plot1+
                  geom_density(data=aes_string(x= "exprs",fill=cluster_name),adjust =1)+
                  mytheme+
                  scale_fill_manual(values = cluster_color(length(unique(transdata[,cluster_name]))))#+
                  #facet_wrap(as.formula(paste0(cluster_name," ~Markers") ),scales = "free_y",ncol = length(density_plot_marker))
  }

  

  if(!is.null(xlim) & free_x==F){
                plot1<-plot1+
                scale_x_continuous(limits=xlim)}
  
  if(free_x==F){
                plot1<-plot1+
                facet_wrap(as.formula(paste0(cluster_name," ~Markers") ),scales = "free_y",ncol = length(density_plot_marker))}
  
  if(free_x==T){
                plot1<-plot1+
                facet_wrap(as.formula(paste0(cluster_name," ~Markers") ),scales = "free",ncol = length(density_plot_marker))}
  
  plot1<-plot1+guides(fill = guide_colorbar(label.theme = element_text(size = 40),
                                            barwidth  = unit(3.5*0.3,"inches"), 
                                            barheight = unit(row_plot*0.8, "inches"), ##图例的高度
                                            raster = T ## 如果为TRUE，则颜色条将呈现为栅格对象。 如果为FALSE，则颜色条呈现为一组矩形
                                            )
                      )
               #theme(legend.margin = margin(t=10,unit="inches"))
  
  multiplot(plot1)
  dev.off()
  cat(paste0("Exported: ","tSNE-","Files"," merge.pdf\n"))
  
}





#' @title draw_density_plots
#' @description
#' 自动生成每个marker的直方图（密度图）；
#'
#' @param combined_data_plot     数据框或者矩阵，附带有(meta)cluster、降维结果(t_sne_1,t_sne_2)的表达矩阵
#' @param groups                 实验组的Metadata信息,其中必需包含列名：”Short_name“
#' @param cluster_color          tSNE图的亚群的配色方案，默认dif_seq_rainbow
#' @param cluster_name           用过分析的cluster名称
#' @param  output_dir            输出数据文件夹名称
#' @param major_cond             选择groups中一个列名做为展现差异的major_cond
#' @param cluster_id             指定显示的cluster
#' @param xlim                   指定x轴范围，例如：c(0,5)，只在free_x=F的时候有效
#' @param free_x                 True或者FALSE，设置是否每张图X-轴取值范围是否自动，建议设置为F
#' @param trans_method           数据转化方式，有四种："CytofAsinh"，"simpleAsinh"；"0_to_Max"，所有Marker的信号强度除以最大值，线性转换，通过除以各通道信号的最大值把数值scale到0~1；"Min_to_Max"，线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异

#' @export

draw_density_plots<-function(combined_data_plot,
                             groups,
                             all_markers,
                             cluster_color=dif_seq_rainbow,
                             cluster_name="metacluster",
                             output_dir="cluster_density_plots",
                             major_cond="Tissue_Type",
                             cluster_id=NULL,
                             xlim=NULL,
                             free_x=F,
                             trans_method="simpleAsinh"
                             )
  
{
  library(colorRamps)
  library(Rmisc)
  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  combined_data_plot$File_ID<-as.character(combined_data_plot$File_ID)
  combined_data_plot<-full_join(combined_data_plot,groups,by="File_ID")
  
  density_plot_marker_id<-which(all_markers[,"density_plot"]==1)
  density_plot_marker<-as.character(all_markers[density_plot_marker_id,"markers"])
  
  
  if(trans_method=="CytofAsinh") trans_fun=CytofAsinh
  if(trans_method=="simpleAsinh") trans_fun=simpleAsinh
  if(trans_method=="0_to_Max") trans_fun=normalize_Zero_Max
  if(trans_method=="Min_to_Max") trans_fun=normalize_min_max
  
  combined_data_plot[, density_plot_marker_id] <- apply(combined_data_plot[, density_plot_marker_id,drop = FALSE],
                                                       2,trans_fun)
  
  if (is.null(cluster_id)){
    cluster_id=unique(combined_data_plot[,cluster_name])} 
  
  combined_data_plot2<-combined_data_plot%>%
    #dplyr::filter(Merge_in_fig=="Yes")%>%
    #dplyr::filter(cluster_name %in% cluster_id)%>%
    dplyr::filter_at(vars(matches(cluster_name)),all_vars(.%in% cluster_id)) %>%
    dplyr::select(one_of(c(cluster_name,density_plot_marker)))
  
  
  
  head(combined_data_plot2)

  
  #str(combined_data_plot2)
  cluster_index<-colnames(combined_data_plot2)==cluster_name
  combined_data_plot2[,cluster_index]<-as.factor(combined_data_plot2[,cluster_index])
  
  
  
  combined_data_plot3<-tidyr::gather_(combined_data_plot2,key="Markers",value="exprs",density_plot_marker)
  
  head(combined_data_plot3)
  

  
  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))
  
  #2生成合并tSNE-文件图像
  
  col_plot<-length(density_plot_marker)*3.5
  row_plot<-length(unique(combined_data_plot3[,cluster_name]))*3.5
  
  pdf(file=paste0("./",output_dir,"/","Density_plots",".pdf"),width = col_plot,height = row_plot)
  
  plot1<-ggplot(combined_data_plot3)+
    geom_density(aes_string(x= "exprs",fill=cluster_name),adjust =1)+
    mytheme+
    scale_fill_manual(values = cluster_color(length(unique(combined_data_plot3[,cluster_name]))))#+
  #facet_wrap(as.formula(paste0(cluster_name," ~Markers") ),scales = "free_y",ncol = length(density_plot_marker))
  
  if(!is.null(xlim) & free_x==F){
    plot1<-plot1+
      scale_x_continuous(limits=xlim)}
  
  if(free_x==F){
    plot1<-plot1+
      facet_wrap(as.formula(paste0(cluster_name," ~Markers") ),scales = "free_y",ncol = length(density_plot_marker))}
  
  if(free_x==T){
    plot1<-plot1+
      facet_wrap(as.formula(paste0(cluster_name," ~Markers") ),scales = "free",ncol = length(density_plot_marker))}
  
  multiplot(plot1)
  dev.off()
  cat(paste0("Exported: ","tSNE-","Files"," merge.pdf\n"))
  
}





#' @title draw_tsne_heatmaps
#' @description
#' 生成各个marker的tSNE；
#'
#' @param combined_data_plot     数据框或者矩阵，附带有(meta)cluster、降维结果(t_sne_1,t_sne_2)的表达矩阵
#' @param edge_data              画网络图的时候，存储网络连接数据
#' @param heatmap_tsne_markers   vector，指定marker的名称,默认NULL，此时将生成所有通道的heatmap
#' @param single_file            是否把所有heatmap合成单个文件，默认FALSE
#' @param groups                 实验组的Metadata信息,其中必需包含列名：”Short_name“
#' @param major_cond             选择groups中一个列名做为展现差异的major_cond
#' @param groups_to_show         vector，指定在图上显示哪些组的数据（组名应属于all_samples的major_cond列）,默认为NULL，显示所有组
#' @param trans_method           数据转化方式，有四种："CytofAsinh"，"simpleAsinh"；"0_to_Max"，所有Marker的信号强度除以最大值，线性转换，通过除以各通道信号的最大值把数值scale到0~1；"Min_to_Max"，线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异
#' @param show_legend            是否显示legend（colorbar）
#' @param legend_type            legend(colorbar)是局部（"local"(默认)）,所有marker用不同的colorbar）还是全局（"global",所有marker共用一个colorbar）
#' @param legend_limit           仅在legend_type是"global"的情况下有效，手动设置colorbar的上下限，取值为一个vector，例如c(0,5) 0为最小值，5为最大值
#' @param show_scale             是否显示X 轴和Y轴的刻度，默认是TRUE
#' @param show_axis              是否显示坐标轴，默认是TRUE
#' @param dot_size               图中点的大小，默认是4
#' @param dot_colours            设置colorbar的颜色，默认为jet color函数
#' @param output_dir             输出数据文件夹名称，默认是TRUE
#' @param output_format          输出数据文件格式："tiff"和“pdf”
#' @param reduction_dm1          combined_data_plot降维产生的维度，默认"tsne_1"
#' @param reduction_dm2          combined_data_plot降维产生的维度，默认"tsne_2"
#' @export


draw_tsne_heatmaps<-function(combined_data_plot,
                             edge_data=NULL,
                             heatmap_tsne_markers=NULL,
                             trans_method="simpleAsinh",
                             single_file=FALSE,
                             groups=groups,
                             major_cond=NULL,
                             groups_to_show=NULL,
                             dot_size=4,
                             dot_colours=colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")),
                             show_legend=TRUE,
                             legend_type="local",
                             legend_limit=NULL,
                             show_scale=TRUE,
                             show_axis=TRUE,
                             
                             show_contour=FALSE,
                             edge_colour="grey",
                             line_size=1,
                             edge_layer="back",
                             
                             contour_line_size=0.5,
                             contour_bins=25,
                             contour_colour="grey",
                             
                            
                             output_dir="tsne_heatmap",
                             output_format="tiff",
                             reduction_dm1="tsne_1",
                             reduction_dm2="tsne_2"
){
  
  if(0){

    combined_data_plot=combined_rawdata_plot
    edge_data=NULL
    heatmap_tsne_markers=NULL
    trans_method="simpleAsinh"
    single_file=FALSE
    groups=groups
    major_cond=NULL
    groups_to_show=NULL
    dot_size=4
    dot_colours=colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    show_legend=TRUE
    legend_type="local"
    legend_limit=NULL
    show_scale=TRUE
    show_axis=TRUE
    
    show_contour=FALSE
    edge_colour="grey"
    contour_line_size=0.5
    contour_bins=25
    contour_colour="grey"
    
    output_dir="tsne_heatmap"
    output_format="tiff"
    reduction_dm1="tsne_1"
    reduction_dm2="tsne_2"
  }
  
  
  
  
  
  
  if (!dir.exists(paste0("./",output_dir," ",trans_method))) {
    dir.create(paste0("./",output_dir," ",trans_method))
  }
  
  if(trans_method=="CytofAsinh") trans_fun=CytofAsinh
  if(trans_method=="simpleAsinh") trans_fun=simpleAsinh
  if(trans_method=="0_to_Max") trans_fun=normalize_Zero_Max
  if(trans_method=="Min_to_Max") trans_fun=normalize_min_max
  
  
  if(is.null(heatmap_tsne_markers)) {
    File_ID_num<-which(colnames(combined_data_plot)=="File_ID")
    heatmap_tsne_markers<-colnames(combined_data_plot[2:File_ID_num-1])}
  
  combined_data_plot[,heatmap_tsne_markers]<-apply(combined_data_plot[,heatmap_tsne_markers,drop=F],
                                                   MARGIN = 2,
                                                   remove_extremum)
  

  combined_data_plot[,heatmap_tsne_markers]<-apply(combined_data_plot[,heatmap_tsne_markers,drop=F],MARGIN=2,trans_fun) #MAGIN=2 按照列进行Normalize，MAGIN=1按照行进行Normalize

  
  xlim<-c(min(combined_data_plot[,reduction_dm1])*1.1,max(combined_data_plot[,reduction_dm1])*1.1)
  ylim<-c(min(combined_data_plot[,reduction_dm2])*1.1,max(combined_data_plot[,reduction_dm2])*1.1)

  add_cell<-combined_data_plot[1,,drop=F]
  suppressWarnings(add_cell[1,]<-0)
  add_cell[1,c(reduction_dm1,reduction_dm2)]=c(1000,1000)
  add_cell$File_ID<-combined_data_plot[1,"File_ID"]

  

  combined_data_plot<-rbind(combined_data_plot,
                            add_cell)
  
  
  combined_data_plot$File_ID<-as.character(combined_data_plot$File_ID)
  
  combined_data_plot<-dplyr::full_join(combined_data_plot,groups,by="File_ID")
  
  tail(combined_data_plot)

  #根据groups_to_show裁剪数据
  
  if((!is.null(groups_to_show))&(is.null(major_cond))) {
    message("Warning: groups_to_show could only be set when major_cond exist, all groups will be showed\n")
    groups_to_show<-NULL
  }
  
  if(is.null(groups_to_show)){
    groups_to_show=unique(groups[,major_cond,drop=T])  
  } else {
    
    combined_data_plot<-combined_data_plot %>%
      dplyr::filter_at(vars(matches(major_cond)),all_vars(.%in% groups_to_show))
    
  }
  
  
  
  #jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))
  
  
  heatmap_tsne_markers<-as.matrix(heatmap_tsne_markers,drop=F)
  
  if (is.null(legend_limit))
    global_legend_limits<-c(min(combined_data_plot[,heatmap_tsne_markers]),max(combined_data_plot[,heatmap_tsne_markers]))
  else
    global_legend_limits<-legend_limit
  
  
  
  #产生facet图背景的数据：
  #等高线图
  if(show_contour==TRUE){
    combined_data_density<-as.data.frame(matrix(nrow=0,ncol=ncol(combined_data_plot)))
    colnames(combined_data_density)<-colnames(combined_data_plot)
    for(major_cond_i in unique(combined_data_plot[,major_cond])){
      combined_data_density_i<-combined_data_plot
      combined_data_density_i[,major_cond]<-major_cond_i
      combined_data_density<-rbind(combined_data_density,
                                   combined_data_density_i)   
    }
  }
  #edge
  
  
  if(!is.null(edge_data)){
  edge_data$edge_id<-as.factor(edge_data$edge_id) 
    if(!is.null(major_cond)){
        edge_data_major_cond<-as.data.frame(matrix(nrow=0,ncol=ncol(edge_data)))
        colnames(edge_data_major_cond)<-colnames(edge_data)
        for(major_cond_i in unique(combined_data_plot[,major_cond])){
          edge_data_major_cond_i<-edge_data
          edge_data_major_cond_i[,major_cond]<-major_cond_i
          edge_data_major_cond<-rbind(edge_data_major_cond,
                                      edge_data_major_cond_i)
        }
    }else{
      
      edge_data_major_cond<-edge_data
    }
    
  }
  
  single_heatmap_tsne<-function(marker_name){
    
    #marker_name="CD28"
    
    single_heatmap<-ggplot()
    
    if(show_contour==TRUE){
      single_heatmap<-single_heatmap+ geom_density_2d(data=combined_data_density,aes_string(x=reduction_dm1,y=reduction_dm2),colour=contour_colour,size=contour_line_size,bins=contour_bins)+
        scale_x_continuous(limits=c(min(combined_data_plot[,reduction_dm1])*1.15,max(combined_data_plot[,reduction_dm1])*1.15))+
        scale_y_continuous(limits=c(min(combined_data_plot[,reduction_dm2])*1.15,max(combined_data_plot[,reduction_dm2])*1.15))
    }
    
    if((!is.null(edge_data)) & (edge_layer=="back")){
      single_heatmap<-single_heatmap+geom_line(data=edge_data_major_cond,aes_string(x=reduction_dm1,y=reduction_dm2,group="edge_id"),color=edge_colour,size=line_size)
      
    }

    
    single_heatmap<-single_heatmap+
                    geom_point(data=combined_data_plot,aes_string(x=reduction_dm1,y=reduction_dm2,colour=marker_name),size=dot_size,alpha=1)+
                    mytheme+
                    #expand_limits(x=xlim,y=ylim)+
                    scale_y_continuous(limits=ylim)+
                    scale_x_continuous(limits=xlim)+
                    theme(text = element_text(size=50))+
                    theme(legend.title=element_blank())+
                    labs(title=marker_name)+
                    theme(legend.key.height=unit(5,"cm"))+
                    theme(plot.margin=unit(c(1,1,1,1),"cm"))

    if((!is.null(edge_data)) & (edge_layer=="front")){
      single_heatmap<-single_heatmap+geom_line(data=edge_data_major_cond,aes_string(x=reduction_dm1,y=reduction_dm2,group="edge_id"),color=edge_colour,size=line_size)
    }
    
    
    if(legend_type=="local"){
      single_heatmap<-single_heatmap+
        scale_color_gradientn(colours = dot_colours(7))
    }else if(legend_type=="global"){
      single_heatmap<-single_heatmap+
        scale_color_gradientn(colours = dot_colours(7),limits=global_legend_limits)
    }
    
    
    if(show_legend==F){
      single_heatmap<-single_heatmap+
        theme(legend.position = "none")
    }
    
    
    if(show_scale==F){
      single_heatmap<-single_heatmap+
        theme(axis.text = element_blank())
    }
    
    
    if(!is.null(major_cond))
    {single_heatmap<-single_heatmap+
      facet_wrap(facets=major_cond,ncol = length(unique(combined_data_plot[,major_cond])))
    }
    
    if(!show_axis){
      single_heatmap<-single_heatmap+
        theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.2))+
        theme(axis.line=element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              axis.text=element_blank())
    }  
    
    return(single_heatmap)
  }

  
  
  heatmap_tsne_list<-lapply(heatmap_tsne_markers,single_heatmap_tsne)
  
  n_figure<-length(heatmap_tsne_list)
  
  
  if(!is.null(major_cond)) { 
    fig_cond_n<- length(unique(combined_data_plot[,major_cond,drop=T]))
    
  }else fig_cond_n=1
  
  
  if(single_file==FALSE){
    
    for(figure_i in 1:n_figure){
      
      cat(paste0("Start to output heatmap of ",heatmap_tsne_markers[figure_i]),"\n")
      if(output_format=="tiff"){
        jpeg(filename = paste0("./",output_dir," ",trans_method,"/","Heatmap ",heatmap_tsne_markers[figure_i],".tiff"),width=1000*fig_cond_n,height=1000)}
      else if(output_format=="pdf"){
        pdf(file=paste0("./",output_dir," ",trans_method,"/","Heatmap ",heatmap_tsne_markers[figure_i],".pdf"),width=1000*fig_cond_n/72,height=1000/72)}
      else{
        message(paste0("Error-Output format ",output_format," is not supported\n"))
        return(NULL)}
      suppressWarnings( multiplot(heatmap_tsne_list[figure_i]))
      dev.off()
    }
  }else{
    fig_row<-ceiling(n_figure^0.5)
    fig_col<-ceiling(n_figure/fig_row)
    
    if(output_format=="tiff"){
      jpeg(filename = paste0("./",output_dir," ",trans_method,"/","tSNE-heatmap.tiff"),height=1000*fig_row,width=1000*fig_col*fig_cond_n)}
    else if(output_format=="pdf"){
      jpeg(filename = paste0("./",output_dir," ",trans_method,"/","tSNE-heatmap.pdf"),height=1000*fig_row/72,width=1000*fig_col*fig_cond_n/72)}  
    else{
      message(paste0("Error-Output format ",output_format," is not supported\n"))
      return(NULL)}
    suppressWarnings(multiplot(plotlist=heatmap_tsne_list,cols = fig_col))
    
    dev.off()
  }
}



#' @title cluster_merge
#' @description
#' 手动合并cluster
#'
#' @param cluster_result   Vector变量，带合并处理cluster   
#' @param merge_list       list变量，指定要合并的cluster名称,例如：list(c(1,2,3),c(6:9)) 将1,2,3合并，将6到9合并；
#' @param cluster_rename   TRUE/FALSE，决定是否对合并后的cluster重新命名

#' @export

cluster_merge<-function(cluster_result,
                        merge_list,
                        cluster_rename=T  ){
                merge_result<-cluster_result
                original_cluster_table<-sort(unique(cluster_result))
                merge_cluster_table<-original_cluster_table
                for (i in c(1:length(merge_list))) {
                  merge_result[merge_result %in% merge_list[[i]]]<-merge_list[[i]][1]
                  merge_cluster_table[merge_cluster_table %in% merge_list[[i]]]<-merge_list[[i]][1]
                }
                if(cluster_rename==T){
                  merge_result<-as.numeric(as.factor(merge_result))
                  merge_cluster_table<-as.numeric(as.factor(merge_cluster_table))
                }
                assign_result<-data.frame(Original=original_cluster_table,
                                          Change_to="------>",
                                          Merge=merge_cluster_table)
                
                message("Cluster merge finished:\n")
                print.data.frame(assign_result,row.names = F)
                return(merge_result)
              }





#' @title One_SENSE_report
#' @description
#' 进行OneSense分析，并进行画出两个坐标轴方向的Heatmap；
#'
#' @param combined_data_sampled     数据框或者矩阵，待处理的表达数据
#' @param perplexity     tsne降维参数，困惑度
#' @param max_iter       tsne降维参数，迭代次数
#' @param all_markers    带有One_SENSE_F和One_SENSE_P的通道选择marker
#' @param output_dir     输出的目录名称
#' @export

One_SENSE_report<-function(combined_data_sampled,

                          #tSNE参数设置：
                           perplexity=30,  #困惑度
                           max_iter=1000,   #迭代次数
                           all_markers=all_markers,
                           output_dir="cluster_onesense_plots"

                           ){

  #One_SENSE_F分析

  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  One_SENSE_F_id=which(all_markers$One_SENSE_F==1)
  One_SENSE_F_name=colnames(combined_data_sampled[,One_SENSE_F_id])
  One_SENSE_F_input_data=combined_data_sampled[,One_SENSE_F_name]

  One_SENSE_F_result <- Rtsne(One_SENSE_F_input_data,
                              initial_dims = ncol(One_SENSE_F_input_data),
                              dims = 1,
                              check_duplicates = FALSE,
                              perplexity=perplexity,
                              max_iter=max_iter)$Y

  #One_SENSE_P分析
  One_SENSE_P_id=which(all_markers$One_SENSE_P==1)
  One_SENSE_P_name=colnames(combined_data_sampled[,One_SENSE_P_id])
  One_SENSE_P_input_data=combined_data_sampled[,One_SENSE_P_name]

  One_SENSE_P_result <- Rtsne(One_SENSE_P_input_data,
                              initial_dims = ncol(One_SENSE_P_input_data),
                              dims = 1,
                              check_duplicates = FALSE,
                              perplexity=perplexity,
                              max_iter=max_iter)$Y



  combined_data_plot <- data.frame(combined_data_sampled,
                                   One_SENSE_F=One_SENSE_F_result,
                                   One_SENSE_P=One_SENSE_P_result)

  One_SENSE_result<-data.frame(One_SENSE_F=One_SENSE_F_result,
                                  One_SENSE_P=One_SENSE_P_result)





  #生成Heatmap



  #统计分析

  One_SENSE_P_Group<-cut(combined_data_plot$One_SENSE_P,50,labels = c(1:50))
  One_SENSE_F_Group<-cut(combined_data_plot$One_SENSE_F,50,labels = c(1:50))



  combined_data_plot<-data.frame(combined_data_plot,
                                     P_Group=One_SENSE_P_Group,
                                     F_Group=One_SENSE_F_Group)




  statics_P_Groups<-aggregate(combined_data_plot[,One_SENSE_P_id],list(combined_data_plot$P_Group),median)[,-1]
  statics_F_Groups<-aggregate(combined_data_plot[,One_SENSE_F_id],list(combined_data_plot$F_Group),median)[,-1]
  head(statics_P_Groups)
  head(statics_F_Groups)

  #标准化：将所有数值Normalize到0~1
  statics_P_Groups=BBmisc::normalize(statics_P_Groups, method = "range", range = c(0, 1), margin = 2)
  statics_P_Groups=as.matrix(statics_P_Groups)

  statics_F_Groups=BBmisc::normalize(statics_F_Groups, method = "range", range = c(0, 1), margin = 2)
  statics_F_Groups=t(as.matrix(statics_F_Groups))


  #绘图
  cat(paste0("Start to output heatmap of P_Heatmap.pdf","\n"))
  #jpeg(filename = paste0("./",output_dir,"/","Heatmap ",heatmap_tsne_markers[figure_i],".jpeg"),width=1000,height=1000,quality = 600)
  pdf(file=paste0("./",output_dir,"/","P_Heatmap.pdf"),width = 8,height = 7)

  heatmap3(x=statics_P_Groups,
           method="average",
           Rowv = NA,
           Colv=NULL,
           balanceColor=FALSE,
           scale="none",
           symm = FALSE,
           showColDendro =F,
           col=colorRampPalette(c("black","blue","green","yellow","Red"),space="rgb",interpolate="linear")(100),
           #col=colorRampPalette(c("black","#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100),

           lasRow=3
  )

  dev.off()

  #绘图
  cat(paste0("Start to output heatmap of F_Heatmap.pdf","\n"))
  #jpeg(filename = paste0("./",output_dir"Heatmap ",heatmap_tsne_markers[figure_i],".jpeg"),width=1000,height=1000,quality = 600)
  pdf(file=paste0("./",output_dir,"/","F_Heatmap.pdf"),width = 8,height = 7)

  heatmap3(x=statics_F_Groups,
           method="average",
           Rowv = NULL,
           Colv=NA,
           balanceColor=FALSE,
           scale="none",
           symm = FALSE,
           showColDendro =T,
           showRowDendro = T,
           col=colorRampPalette(c("black","blue","green","yellow","Red"),space="rgb",interpolate="linear")(100),
           lasRow=2
           )

  dev.off()

  return(One_SENSE_result)

}




