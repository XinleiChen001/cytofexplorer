

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

brewer_sets<-c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"))

colorset<-as.data.frame(matrix(ncol=1,nrow = n))

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
    metaClusters <-suppressMessages(MetaClustering(data=PhenoGraph_input_data,
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
#' @param sample_type      Sample的方式，#sample_type=1,从每个cluster抽取相同细胞；sample_type=2，从每个样本抽取相同细胞
#' @param sample_n         每个样本或者每个cluster抽取的细胞数量；
#' @param cluster_name     cluster通道的名称
#' @return 返回值一个矩阵，包含所有sample取得的表达数据
#' @export

equal_sample<-function( x,
                        sample_type=1, #sample_type=1,从每个cluster抽取相同细胞；sample_type=2，从每个样本抽取相同细胞
                        sample_n=1000, #每个样本或者每个cluster抽取的细胞数量；
                        cluster_name="metaCluster"
)
{
  x=combined_data_analysed

  if(sample_type==1) sample_col=cluster_name
  if(sample_type==2) sample_col="File_ID"

  combined_data_sampled<-as.data.frame(matrix(nrow=0,ncol=ncol(x)))
  for(meta_i in unique(x[,sample_col])){
    #meta_i=8
    cell_n=length(which(x[,sample_col]==meta_i))
    if (cell_n>1000) {
      cell_n=1000
    }
    data_sampled<-x%>%
      dplyr::filter_(paste0(sample_col,"==meta_i"))%>%
      sample_n(cell_n)
    combined_data_sampled<-rbind(data_sampled,combined_data_sampled)
  }

  clustern<-length(unique(x[,sample_col]))

  #cat(paste0("Total sampled events：",nrow(combined_data_sampled)))

  layout(matrix(c(1,3,2,4),ncol=2))

  smr_data_total<-x %>% count_(cluster_name)   #smr is short for summerise

  plot1<-barplot(smr_data_total$n,
                 names.arg=smr_data_total[,1,drop=T],
                 main="Events of metaclusters(Total)")

  smr_data_sampled<-combined_data_sampled %>% count_(cluster_name)

  plot2<-barplot(smr_data_sampled$n,
                 names.arg=smr_data_sampled[,1,drop=T],
                 main="Events of metaclusters (Sampled)")


  smr_data_total<-x %>% count_("File_ID")

  plot3<-barplot(smr_data_total$n,
                 #names.arg=smr_data_sampled[,1,drop=T],
                 main="Events of Files (Total)")

  smr_data_sampled<-combined_data_sampled %>% count_("File_ID")

  plot4<-barplot(smr_data_sampled$n,
                 #names.arg=smr_data_sampled[,1,drop=T],
                 #axis(side=1,lwd=0,las=1),
                 main="Events of Files (Sampled)")


  layout(matrix(c(1),ncol=1))


  return(combined_data_sampled)

}




#' @title draw_tsne_maps
#' @description
#' 自动生成系列tSNE图；
#'
#' @param combined_data_plot     数据框或者矩阵，附带有(meta)cluster、降维结果(t_sne_1,t_sne_2)的表达矩阵
#' @param groups                 实验组的Metadata信息,其中必需包含列名：”Short_name“
#' @param cluster_color          tSNE图的亚群的配色方案，默认dif_seq_rainbow
#' @param cluster_name           选择groups中一个列名做为展现差异的major_cond
#'@param  output_dir              输出数据文件夹名称
#' @param reduction_dm1          combined_data_plot降维产生的维度，默认"tsne_1"
#' @param reduction_dm2          combined_data_plot降维产生的维度，默认"tsne_2"
#' @export


draw_tsne_figs<-function(
  combined_data_plot,
  groups,
  cluster_color=dif_seq_rainbow,
  cluster_name="metacluster",
  output_dir="cluster_tsne_plots",
  major_cond="Tissue_Type",
  reduction_dm1="tsne_1",
  reduction_dm2="tsne_2"){

  library(colorRamps)
  library(Rmisc)


  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }




  combined_data_plot<-full_join(combined_data_plot,groups,by="File_ID")

  cluster_index<-colnames(combined_data_plot)==cluster_name
  combined_data_plot[,cluster_index]<-as.factor(combined_data_plot[,cluster_index])

  major_cond_index<-colnames(combined_data_plot)==major_cond
  major_cond_n<-length(unique(combined_data_plot[,major_cond_index,drop=TRUE]))


  combined_data_group<-combined_data_plot %>%
    group_by_(cluster_name)
  centers<-eval(parse(text=paste0("summarise(combined_data_group,",reduction_dm1,"=median(",reduction_dm1,"),",reduction_dm2,"=median(",reduction_dm2,"))")))



  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))


  #生成合并tSNE-Phenograph图像
  pdf(file=paste0("./",output_dir,"/","tSNE-",cluster_name," merge.pdf"),width = 8,height = 7)
  #jpeg(filename = paste0("tSNE-Phenograph merge.jpeg"),width=1000,height=1000,quality = 600)
  plot1<-ggplot(combined_data_plot)+
    geom_point(aes_string(x=reduction_dm1,y=reduction_dm2,colour=cluster_name),size=1)+
    mytheme+
    geom_text(data=centers,aes_string(x=reduction_dm1,y=reduction_dm2),label=centers[,1,drop=TRUE],colour="black",size=5)+
    #scale_color_manual(values = primary.colors(nrow(centers)+1)[-1])#facet_wrap(~File_ID)
    scale_color_manual(values = cluster_color(nrow(centers)))
  multiplot(plot1)
  dev.off()
  cat(paste0("Exported: ","tSNE-",cluster_name," merge.pdf\n"))


  #生成单个文件的tSNE-Phenograph图像
  jpeg(filename = paste0("./",output_dir,"/","tSNE-",cluster_name,"(sample name).jpeg"),width=1200*major_cond_n,height=1200*ceiling(nrow(groups)/major_cond_n),quality = 100)
  plot2<-ggplot(combined_data_plot)+
    geom_point(aes_string(x=reduction_dm1,y=reduction_dm2,colour=cluster_name),size=4)+
    mytheme+
    theme(text = element_text(size=40))+
    #theme(legend.key.height=unit(10,"cm"))+
    facet_wrap(~Short_name,ncol=major_cond_n)+
    scale_color_manual(values = cluster_color(nrow(centers)))
  multiplot(plot2)
  dev.off()
  cat(paste0("Exported: ","tSNE-",cluster_name,"(sample name).jpeg\n"))



  #生成单个分组的tSNE-Phenograph图像

  tiff(filename = paste0("./",output_dir,"/","tSNE-",cluster_name,"(",major_cond,").tif"),width=1200*major_cond_n,height=1200)
  plot3<-ggplot(combined_data_plot)+
    geom_point(aes_string(x=reduction_dm1,y=reduction_dm2,colour=cluster_name),size=4)+
    mytheme+
    theme(text = element_text(size=40))+
    #theme(legend.key.height=unit(10,"cm"))+
    facet_wrap(facets=major_cond,ncol=major_cond_n)+
    #scale_color_manual(values = primary.colors(nrow(centers)+1)[-1])
    scale_color_manual(values = cluster_color(nrow(centers)))
  multiplot(plot3)
  dev.off()
  cat(paste0("Exported: ","tSNE-",cluster_name,"(",major_cond,").tif\n"))


  #生成两组之间merge的tSNE-Phenograph图像
  tiff(filename = paste0("./",output_dir,"/","tSNE colored by ",major_cond,".tiff"),width=1200,height=1200)
  plot4<-ggplot(combined_data_plot)+
    theme(text = element_text(size=30))+
    #theme(legend.key.height=unit(10,"cm"))+
    geom_point(aes_string(x=reduction_dm1,y=reduction_dm2,colour=major_cond),size=4,alpha=0.6)+
    mytheme+
    scale_color_brewer(palette = "Set1")
  multiplot(plot4)
  dev.off()
  cat(paste0("Exported: ","tSNE colored by ",major_cond,".tiff\n"))


}






#' @title draw_tsne_heatmaps
#' @description
#' 生成各个marker的tSNE；
#'
#' @param combined_data_plot     数据框或者矩阵，附带有(meta)cluster、降维结果(t_sne_1,t_sne_2)的表达矩阵
#' @param heatmap_tsne_markers   vector，指定marker的名称,默认NULL，此时将生成所有通道的heatmap
#' @param single_file            是否把所有heatmap合成单个文件，默认FALSE
#'@param output_dir    输出数据文件夹名称
#' @param reduction_dm1          combined_data_plot降维产生的维度，默认"tsne_1"
#' @param reduction_dm2          combined_data_plot降维产生的维度，默认"tsne_2"
#' @export



draw_tsne_heatmaps<-function(combined_data_plot=combined_data_plot,
                             heatmap_tsne_markers=NULL,
                             single_file=FALSE,
                             output_dir="tsne_heatmap",
                             reduction_dm1="tsne_1",
                             reduction_dm2="tsne_2"
){


  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  if(is.null(heatmap_tsne_markers)) {
    File_ID_num<-which(colnames(combined_data_plot)=="File_ID")
    heatmap_tsne_markers<-colnames(combined_data_plot[2:File_ID_num-1])}

  combined_data_plot[,heatmap_tsne_markers]<-apply(combined_data_plot[,heatmap_tsne_markers],
                                                   MARGIN = 2,
                                                   remove_extremum)



  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))


  heatmap_tsne_markers<-as.matrix(heatmap_tsne_markers)

  single_heatmap_tsne<-function(marker_name){
    single_heatmap<-ggplot(combined_data_plot)+
      geom_point(aes_string(x=reduction_dm1,y=reduction_dm2,colour=marker_name),size=4,alpha=1)+
      mytheme+scale_color_gradientn(colours = jet.colors(7))+
      theme(text = element_text(size=50))+
      theme(legend.title=element_blank())+
      labs(title=marker_name)+
      theme(legend.key.height=unit(5,"cm"))+
      theme(plot.margin=unit(c(1,1,1,1),"cm"))


    return(single_heatmap)
  }

  heatmap_tsne_list<-lapply(heatmap_tsne_markers,single_heatmap_tsne)

  n_figure<-length(heatmap_tsne_list)

  if(single_file==FALSE){

    for(figure_i in 1:n_figure){

      cat(paste0("Start to output heatmap of ",heatmap_tsne_markers[figure_i]),"\n")
      jpeg(filename = paste0("./",output_dir,"/","Heatmap ",heatmap_tsne_markers[figure_i],".jpeg"),width=1000,height=1000,quality = 600)
      #pdf(file=paste0("Heatmap ",heatmap_tsne[figure_i],".pdf"),width = 8,height = 7)
      #    width=8,
      #    height = 6)
      #multiplot(plotlist=heatmap_tsne_list,cols = ceiling(n_figure^0.5))
      multiplot(heatmap_tsne_list[figure_i])
      dev.off()
    }
  }else{
    fig_row<-ceiling(n_figure^0.5)
    fig_col<-ceiling(n_figure/fig_row)
    jpeg(filename = paste0("./",output_dir,"/","tSNE-heatmap.jpeg"),width=1000*fig_row,height=1000*fig_col,quality = 600)
    multiplot(plotlist=heatmap_tsne_list,cols = fig_col)
    dev.off()
  }
}

