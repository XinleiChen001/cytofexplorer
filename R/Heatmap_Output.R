

#' @title draw_expr_heatmap
#'
#' @description
#' 生成cluster的Heatmap
#'
#'@param xdata                data.frame,或者matrix，cluster expression matrix
#'@param trans_method         数据转化方式，有三种："CytofAsinh"，Arcsinh转换所有表达数据；"0_to_Max"，所有Marker的信号强度除以最大值，线性转换，通过除以各通道信号的最大值把数值scale到0~1；"Min_to_Max"，线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异
#'@param Rowv,Colv            逻辑变量，分别设定行和列是否聚类
#'@param dendrogram           显示行或者列的树形图,"both","row","column","none"
#'@param color_style          heatmap颜色模式，1：黑-黄；2：白：红；3: jet; 其他数字：手动模式，使用colorkeys参数设置
#'@param colorkeys            手动设置heatmap颜色，默认c("black","yellow")
#'
#'@return none
#'@export

draw_expr_heatmap<-function(xdata,
                        trans_method="CytofAsinh",
                        Rowv=T,
                        Colv=T,
                        dendrogram="both",
                        output_dir=paste0("Cluster_Expression_Heatmap"),
                        color_style=1,
                        colorkeys=c("black","yellow")
){
  

 
  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }
  write.csv(xdata,paste0("./",output_dir,"/Cluster_expression_heatmap.csv"),row.names = FALSE)
  
  if(trans_method=="CytofAsinh") trans_fun=simpleAsinh
  if(trans_method=="0_to_Max") trans_fun=normalize_Zero_Max
  if(trans_method=="Min_to_Max") trans_fun=normalize_min_max

  if(!is.null(color_style)){
  if(color_style==1) colorkeys=c("black","yellow")
  if(color_style==2) colorkeys=c("white","Red")
  if(color_style==3) colorkeys=c( "#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
  }


  #数据转换
  
  transformed_data <- apply(xdata,MARGIN=2,trans_fun) #MAGIN=2 按照列进行Normalize，MAGIN=1按照行进行Normalize

  rownames(transformed_data)<-rownames(xdata)
  ### 画Heatmap
  
  width_fig=ncol(xdata)*0.5+4
  height_fig=nrow(xdata)*0.2+4

  library(pheatmap)
  pheatmap(mat=transformed_data,
           cluster_rows=Rowv,
           cluster_cols = Colv,
           scale = "none",
           col=colorRampPalette(colorkeys)(100),
           #legend_labels=trans_method,
           cellwidth = 10,
           cellheight = 8,
           fontsize = 6,
           filename = paste0("./",output_dir,"/Cluster_expression_heatmap.pdf")
           )
  

}



### 数据预处理
#定义转化函数



## 所有Marker的信号强度除以最大值
##线性转换，通过除以各通道信号的最大值把数值scale到0~1

### 数据预处理
#定义转化函数

normalize_Zero_Max<-function(value)
{
  if(!all(value==0)){
    value<-value/max(value)
  }
  return(value)
}




## Min to Max
## 线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异
### 数据预处理

#定义转化函数
normalize_min_max<-function(value)
{
  if(!all(value==0)){
    value<-(value-min(value))/(max(value)-min(value))
  }
  return(value)
}




