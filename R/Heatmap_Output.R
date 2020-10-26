
#draw_expr_heatmaps  
#2020_09_04 matches-->one_of
#2020_10_24 增加max_b参数
#2020—10——25 修复barplot防反的bug


#' @title draw_expr_heatmaps
#'
#' @description
#' 生成cluster的Heatmap
#'
#'@param xdata                包含原始表达数据及FileID的data frame，取值combinded_data_raw
#'@param cluster_stat         stat_by_cluster的结果
#'@param all_markers          从 “all_markers.csv”读入的dataframe，包含marker的选择信息
#'@param groups               groups 从“all_sample.csv”中读入的dataframe，包含实验分组信息（也就是元数据)
#'@param major_cond           在groups所提供的各种分组中，选择一种做为统计差异的依据，例如"Tissue_Type"
#'@param summerise_method     统计过程中使用的合计方法：有"mean"或者 "median"两种
#'@param cluster_name         统计所依据的聚类参数名称，例如"PhenoGraph"、"metacluster"
#'@param groups_to_show       想要展示部分组时，指定展示组的组名，默认为NULL，展示所有组
#'@param cluster_id           指定输出cluster的id，例如：如果要输出前五个cluster，cluster_id=c(1:5); 如果要输出全部cluster，保持其默认值NULL即可。
#'@param trans_method         数据转化方式，有五种：
#'                            "CytofAsinh"，CytofAsinh转换所有表达数据；
#'                            "simpleAsinh",simpleAsinh转换所有表达数据,（公式为asinh(x/5)）,这是默认设置；
#'                            "logicle",logicle转换所有表达数据,logical在0附近近似线性，远离0时，则表现类似对数，同时保证了不同转化方法的平滑过度。；
#'                            "0_to_Max"，所有Marker的信号强度除以最大值，线性转换，通过除以各通道信号的最大值把数值scale到0~1；
#'                            "Min_to_Max"，线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异；
#'                            "simpleAsinh_0_to_Max" 先Arcsinh转换，然后除以各通道信号的最大值把数值scale到0~1范围内；
#'                            "simpleAsinh_Min_to_Max" 先Arcsinh转换，然后最小值转换成0，最大值转换成1，最大限度展现population的表达差异；
#'                            注意：为方便观察，density plot x轴数据固定使用simpleAsinh，不受该参数影响
#'@param max_b               "0_to_Max"和"simpleAsinh_0_to_Max“设定Max数值的最小值，以避免过于突出显示一些较弱通道的信号。                           
#'@param Rowv,Colv            逻辑变量，分别设定行和列是否聚类
#'@param output_dir           输出数据文件夹名称

#'@param use_marker_col        字符串,取值为all_markers的一个列名，用来手动指定all_markers的列作为marker选择的依据，默认，"heatmap"
#'@param density_plot_fill_by  决定density plot填充颜色所依据的参数：取值有 "markers","clusters"，"groups"，"expression"，"none"(注意，离散参数名后面有s)
#'@param density_plot_color_by 决定density plot线条颜色所依据的参数：取值有 "markers","clusters"，"groups"，"expression"，"none"(注意，离散参数名后面有s)
#'@param density_bkgd_fill_by  决定density plot背景颜色所依据的参数，取值有两种方式，一种保持与density_plot_fill_by相同的值，一种直接输入R能识别的颜色名称例如："grey"
#'@param facet_by_groups       逻辑变量，TRUE时，将groups_to_show中的组别输出独立的heatmap，FALSE时，将所有组数据输出成单个的Heatmap

#'@param expression_color_set  调色板或任意其他颜色函数，指定expression通道的颜色设置
#'@param marker_color_set      调色板或任意其他颜色函数，指定markers通道的颜色设置
#'@param group_color_set       调色板或任意其他颜色函数，指定groups通道的颜色设置
#'@param cluster_color_set     调色板或任意其他颜色函数，指定clusters通道的颜色设置


#'@param show_abundance_barplot 逻辑变量，决定是否在heatmap上显示abundent barplot
#'@param output_density_plot    逻辑变量，决定是否显示density plot


#'@param density_line_size     density plot参数，设置线条粗细，默认2
#'@param density_fill_alpha    density plot参数，设置图片透明度，默认0.8（取值0~1之间，数值越大透明度越低）
#'@param density_plot_free_x   density plot参数，布尔变量，决定是否每幅小图可以根据其包含数据自行决定x标轴范围，默认为F
#'@param density_plot_xlim     density plot参数，设置x轴的范围，取值为一个vector，例如c(0,5) 就是指定x轴范围是0~5，注意，该参数仅在density_plot_free_x=F的时候有效
#'@param strip_label_size      density plot参数，设置每个density plot小图上面的标签（包含cluster和marker名称）的字体大小

#'
#'@return none
#'@export

draw_expr_heatmaps<-function(xdata,
                             cluster_stat=NULL,

                             #以下5个参数均可从cluster_stat读取
                             all_markers=NULL,
                             groups=NULL,
                             major_cond=NULL,
                             
                             summerise_method =NULL,
                             cluster_name=NULL,
                             groups_to_show=NULL,
                             cluster_id=NULL,
                             trans_method="simpleAsinh",
                             max_b=0,
             
                             #marker和cluster排序
                             Rowv=T,
                             Colv=T,
                             output_dir=paste0("Cluster_Expression_Heatmap"),
                             use_marker_col="heatmap",
                             
                             #Density图片axe设置
                             density_plot_fill_by="expression",
                             density_plot_color_by="none",
                             density_bkgd_fill_by="none",  #"white",same with density
                             facet_by_groups=F,
                             
                             #各个参数颜色设置
                             expression_color_set=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu"))),
                             marker_color_set=rainbow,
                             group_color_set=brewer_color_sets,
                             cluster_color_set=dif_seq_rainbow,
                             

                             show_abundance_barplot=TRUE,
                             output_density_plot=TRUE,
                             density_line_size=2,
                             density_fill_alpha=0.8,
                             density_plot_free_x=F,
                             density_plot_xlim=NULL,
                             strip_label_size=15){
  
  if(0){

    xdata=combined_data_raw
   
    
    #以下5个参数均可从cluster_stat读取
    all_markers=NULL
    groups=NULL
    major_cond=NULL
    
    summerise_method =NULL
    cluster_name=NULL
    groups_to_show=NULL
    
    cluster_id=NULL
    trans_method="simpleAsinh"
    
    #marker和cluster排序
    Rowv=T
    Colv=T
    #dendrogram="both"
    output_dir=paste0("Cluster_Expression_Heatmap")
    use_marker_col="heatmap"
    
    #Density图片axe设置
    density_plot_fill_by="expression"
    density_plot_color_by="none"
    density_bkgd_fill_by="none"  #"white",same with density
    facet_by_groups=F
    
    #各个参数颜色设置
    expression_color_set=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))
    marker_color_set=rainbow
    group_color_set=brewer_color_sets
    cluster_color_set=dif_seq_rainbow
    
    
    show_abundance_barplot=TRUE
    
    output_density_plot=TRUE
    density_line_size=2
    density_fill_alpha=0.8
    density_plot_free_x=F    
    density_plot_xlim=NULL
    strip_label_size=15
    
  }
  
  merged_cluster_id<-cluster_stat[["merged_cluster_id"]]
  if(is.null(groups))groups<-cluster_stat[["groups"]]
  if(is.null(major_cond))   major_cond <- as.character(cluster_stat[["metadata"]]$major_cond)
  if(is.null(summerise_method))    summerise_method  <-as.character(cluster_stat[["metadata"]]$summerise_method)
  if(is.null(cluster_name))   cluster_name      <- as.character(cluster_stat[["metadata"]]$cluster_name)
  if(is.null(all_markers))   all_markers       <-cluster_stat[["all_markers"]]
                            heatmap_ID        <-all_markers%>%
                                                dplyr::filter_at(vars(one_of(use_marker_col)),all_vars(.== 1))
                            heatmap_ID        <-heatmap_ID[,"markers",drop=T] 
                            # 
                            # legend_breaks=NA
                            # legend_labels=NA
  
  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  #if(is.null(output_fig_name)) output_fig_name<-"Cluster_expression_heatmap"

                            
if(trans_method=="CytofAsinh") trans_fun=cytofAsinh
if(trans_method=="simpleAsinh") trans_fun=simpleAsinh
if(trans_method=="logicle") trans_fun=logicle
if(trans_method=="0_to_Max") trans_fun=normalize_Zero_Max
if(trans_method=="Min_to_Max") trans_fun=normalize_min_max
if(trans_method=="simpleAsinh_0_to_Max") trans_fun=simpleAsinh_Zero_Max
if(trans_method=="simpleAsinh_Min_to_Max") trans_fun=simpleAsinh_Min_to_Max

                            
                            
  if(summerise_method=="mean"){
          use_method<-mean}else
            if(summerise_method=="median")
            { use_method<-median
            }
  
  if(is.null(groups_to_show)){
              groups_to_show=unique(groups[,major_cond,drop=T])  
            } 
  groups_to_show_list<-list()
  if(facet_by_groups==F)
            { groups_to_show_list[[1]]<-groups_to_show}else{
            groups_to_show_list<-as.list(groups_to_show)
           }
  # if (is.null(cluster_id)){
  #   cluster_id=unique(xdata[,cluster_name])} 
  if(is.null(cluster_id[1])){
    cluster_id=merged_cluster_id
  }
  if(density_bkgd_fill_by=="none"){
    density_bkgd_fill_by<-"white"}
  
  if(density_bkgd_fill_by!=density_plot_fill_by){
    x.inv<-try(colorRampPalette(density_bkgd_fill_by)(1),silent = TRUE)
    
    if ('try-error' %in% class(x.inv)) {
      message("Warning-  Wrong color name detected, density_background_color is set to white\n")
      density_bkgd_fill_by<-"white"}
  }
  if(density_bkgd_fill_by=="groups"){
     density_bkgd_fill_by<-"white"
     message("Warning-  density_background_color cannot set to groups ,set value to white\n")}

 
  
  #根据所有组整合数据，设置cluster和marker顺序
  #xdata<-combined_data_raw
  datacolname<-colnames(xdata)
  datacolname[1:(nrow(all_markers)-1)]<-as.character(all_markers$markers)[1:(nrow(all_markers)-1)]
  colnames(xdata)<-datacolname
  event_data_raw<-xdata %>%
                  dplyr::full_join(groups,by="File_ID")%>%
                  dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% groups_to_show))%>%
                  dplyr::filter_at(vars(one_of(cluster_name)),all_vars(.%in% cluster_id))
  

  
  #Generate density plot data (to be faceted or not)
  
  event_data_trans<-event_data_raw
  event_data_trans[,heatmap_ID]<-apply(event_data_trans[,heatmap_ID,drop=F],MARGIN=2,simpleAsinh) #MAGIN=2 按照列进行Normalize，MAGIN=1按照行进行Normalize
  dstplot_data<-event_data_trans%>%   #dstplot_data: short for density plot data
                dplyr::filter_at(vars(one_of(cluster_name)),all_vars(.%in% cluster_id)) %>%
                dplyr::select(one_of(c(major_cond,cluster_name,heatmap_ID)))
  #dstplot_data[,cluster_name]<-as.character(dstplot_data[,cluster_name])
  
  dstplot_data_gather<-tidyr::gather_(dstplot_data,key="markers",value="exprs",heatmap_ID)
  colnames(dstplot_data_gather)<-sub(cluster_name,"clusters",colnames(dstplot_data_gather))
  colnames(dstplot_data_gather)<-sub(major_cond,"groups",colnames(dstplot_data_gather))
  
  #Generate heatmap data not to be faceted
  
  cluster_data_raw <-event_data_raw%>% 
                    group_by_at(c(cluster_name)) %>%
                    summarise_if(is.numeric,use_method,na.rm=TRUE)%>%
                    select_at(vars(one_of(heatmap_ID,cluster_name)))%>%
                    data.frame()
  
  row.names(cluster_data_raw)<-as.character(cluster_data_raw[,cluster_name])
  
  #row.names(cluster_data_raw)<-paste0("cluster",cluster_data_raw[,cluster_name])
  cluster_data_trans     <-cluster_data_raw

  if(trans_method %in% c("0_to_Max","simpleAsinh_0_to_Max")){
       cluster_data_trans[,heatmap_ID]    <- apply(cluster_data_raw[,heatmap_ID],MARGIN=2,trans_fun,b=max_b) #MAGIN=2 按照列进行Normalize，MAGIN=1按照行进行Normalize
  }else{
       cluster_data_trans[,heatmap_ID]    <- apply(cluster_data_raw[,heatmap_ID],MARGIN=2,trans_fun) #MAGIN=2 按照列进行Normalize，MAGIN=1按照行进行Normalize
  }
  
    heatmap_data<-cluster_data_trans[,c(cluster_name,heatmap_ID)]
  
  melt_htdata<-melt(heatmap_data,id.vars = cluster_name,measure.vars=heatmap_ID,variable.name="markers",value.name = "expression")
  colnames(melt_htdata)<-sub(cluster_name,"clusters",colnames(melt_htdata))
  max_expr<-max(melt_htdata$expression)
  
  model_cluster <- hclust(dist(heatmap_data[heatmap_ID]), "complete")
  model_marker  <- hclust(dist(t(heatmap_data[heatmap_ID])), "complete")
  dhc_cluster   <- as.dendrogram(model_cluster)
  dhc_marker    <- as.dendrogram(model_marker)
  
  event_data_trans
  
  
  abundance_boxplot_data<-event_data_raw
  
  
  if(facet_by_groups==T){
   #Generate heatmap data to be faceted
            cluster_data_raw_facet <-event_data_raw%>% 
                                      group_by_at(c(cluster_name,major_cond)) %>%
                                      summarise_if(is.numeric,use_method,na.rm=TRUE)%>%
                                      select_at(vars(one_of(heatmap_ID,cluster_name,major_cond)))%>%
                                      data.frame()
            #row.names(cluster_data_raw_facet)<-as.character(cluster_data_raw_facet[,cluster_name])
            
            #row.names(cluster_data_raw)<-paste0("cluster",cluster_data_raw[,cluster_name])
            cluster_data_trans_facet     <-cluster_data_raw_facet
            
            if(trans_method %in% c("0_to_Max","simpleAsinh_0_to_Max")){
                  cluster_data_trans_facet[,heatmap_ID]    <- apply(cluster_data_raw_facet[,heatmap_ID],MARGIN=2,trans_fun,b=max_b) #MAGIN=2 按照列进行Normalize，MAGIN=1按照行进行Normalize
            }else{
                  cluster_data_trans_facet[,heatmap_ID]    <- apply(cluster_data_raw_facet[,heatmap_ID],MARGIN=2,trans_fun) #MAGIN=2 按照列进行Normalize，MAGIN=1按照行进行Normalize
            }

            
            heatmap_data_facet<-cluster_data_trans_facet[,c(cluster_name,heatmap_ID,major_cond)]
            melt_htdata_facet<-melt(heatmap_data_facet,id.vars = c(cluster_name,major_cond),measure.vars=heatmap_ID,variable.name="markers",value.name = "expression")
            colnames(melt_htdata_facet)<-sub(cluster_name,"clusters",colnames(melt_htdata_facet))
            colnames(melt_htdata_facet)<-sub(major_cond,"groups",colnames(melt_htdata_facet))
            melt_htdata<-melt_htdata_facet
            dstplot_data_gather<-dplyr::full_join(dstplot_data_gather,melt_htdata_facet,by=c("clusters","markers","groups"))
            groups_to_show_list<-as.list(groups_to_show)
            #max_expr<-max(dstplot_data_gather$expression)
            
  }else{
    
    dstplot_data_gather<-dplyr::full_join(dstplot_data_gather,melt_htdata,by=c("clusters","markers"))
    dstplot_data_gather$groups="merge"
    groups_to_show_list<-list()
    #groups_to_show_list[[1]]<-c("merge",unique(as.character(event_data_raw[,major_cond])))
    groups_to_show_list[[1]]<-c("merge")
    abundance_boxplot_data[,major_cond]="merge"
    melt_htdata$groups<-"merge"
    
  }


  head(dstplot_data_gather)
  head(melt_htdata)

  #排序
  #设置默认顺序
  
  dstplot_data_gather[,"clusters"]<-factor(dstplot_data_gather[,"clusters"],levels =rev(cluster_id),ordered = TRUE )
  melt_htdata[,"clusters"]<-factor(melt_htdata[,"clusters"],levels =rev(cluster_id),ordered = TRUE )
  dstplot_data_gather$markers    <-factor(dstplot_data_gather$markers,levels =heatmap_ID,ordered = TRUE)
  melt_htdata[,"markers"]        <-factor(melt_htdata[,"markers"],    levels =heatmap_ID,ordered = TRUE )
  
  #按照聚类结果排序

  if(Rowv==T){
              dstplot_data_gather[,"clusters"]<-factor(dstplot_data_gather[,"clusters"],levels =labels(dhc_cluster),ordered = TRUE )
              melt_htdata[,"clusters"]<-factor(melt_htdata[,"clusters"],levels =labels(dhc_cluster),ordered = TRUE )}
  if(Colv==T){
              dstplot_data_gather$markers    <-factor(dstplot_data_gather$markers,levels =labels(dhc_marker),ordered = TRUE)
              melt_htdata[,"markers"]<-factor(melt_htdata[,"markers"],levels =labels(dhc_marker),ordered = TRUE )
  }

  head(melt_htdata)
  
  
  #作图
  for (this_group in groups_to_show_list) {
    
      #this_group="merge"  
              this_dstplot_data_gather<-dplyr::filter_at(dstplot_data_gather,vars(one_of("groups")),all_vars(.%in% as.character(this_group)))
              this_melt_htdata<-dplyr::filter_at(melt_htdata,vars(one_of("groups")),all_vars(.%in% this_group))
    
              # str(this_melt_htdata)
              
              mytheme <- theme(panel.background = element_rect(fill = "transparent", colour = "NA", size = 0.2), #坐标系及坐标轴
                               legend.key = element_rect(fill = "white", colour = "white"), #图标
                               #legend.key = element_rect( colour = "white"), #图标,
                               panel.grid = element_blank()) 
              

              
              #生成boxplot
              
                
                #生成cluster abundance的barplot(not facet)
                cluster_abandance_stat<-abundance_boxplot_data %>%
                                        dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% as.character(this_group)))%>%
                                        group_by_at(c(cluster_name)) %>%
                                        dplyr::summarise(num=n())
                
                # cluster_num<-length(unique(event_data_raw[,cluster_name]))
                # if(nrow(event_data_raw)!=cluster_num){
                #   message("Warning:cluster number in cluster_barplot_data is not equal to combined_data_plot.\n")
                # }
                
                cluster_abandance_stat$percentage=cluster_abandance_stat$num/sum(cluster_abandance_stat$num)*100
                cluster_abandance_stat<-data.frame(cluster_abandance_stat)
                #cluster_abandance_stat$fill=cluster_color(nrow(cluster_abandance_stat))[cluster_abandance_stat$cluster]
                cluster_abandance_stat[,cluster_name]<-as.character(cluster_abandance_stat[,cluster_name])
                cluster_abandance_stat[,cluster_name]<-factor(cluster_abandance_stat[,cluster_name],levels=(labels(dhc_cluster)),ordered = T)
                cluster_abundance_barplot<-ggplot(cluster_abandance_stat)+
                geom_bar(aes_string(x=cluster_name,y="percentage"),stat = "identity")+
                #scale_fill_manual(values = cluster_abandance_stat$fill)+
                mytheme+
                xlab(NULL)+
                scale_x_discrete(breaks=NULL)+ 
                coord_flip()+#横向
                #ylab("Percentage")+
                theme(axis.text.x = element_text(size = 20,color="black"))+
                theme(axis.title.x = element_text(size = 20,color="black"))+
                
                #scale_x_discrete(limits=rev(levels(cluster_abandance_stat[,cluster_name])))+
                theme(plot.margin=unit(c(0,0,0,0), unit = "cm"))
              

              
                            
                #导出对应的heatmap
                
                gblank<-ggplot()+
                        mytheme+
                        #scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
                        theme(plot.margin=unit(c(0,0,0,0), unit = "cm"))
                      
                ddata_cluster <- dendro_data(dhc_cluster, type = "rectangle")
                ddata_marker <- dendro_data(dhc_marker, type = "rectangle")
                
                if(facet_by_groups==F){
                          
                          ddplot_cluster <-  ggplot(segment(ddata_cluster)) + 
                                              geom_segment(aes(x = x, y = y, xend = xend, yend = yend),size=1) + 
                                              theme(plot.margin=unit(c(0,0,0,0.2), unit = "cm"))+
                                              #scale_x_continuous(limits = c(1,filenum+0.01),breaks = NULL,minor_breaks = NULL,expand = c(0.5/filenum,0.5/filenum))+
                                              scale_x_continuous(expand = c(0.5/length(labels(dhc_cluster)), 0.5/length(labels(dhc_cluster))))+
                                              coord_flip() + 
                                              scale_y_reverse(expand=c(0,0))+
                                              mytheme+
                                              theme(panel.background=element_rect(fill='transparent', color='white'))+
                                              theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
                                            
                          ddplot_marker <-  ggplot(segment(ddata_marker)) + 
                                              geom_segment(aes(x = x, y = y, xend = xend, yend = yend),size=1) + 
                                              theme(plot.margin=unit(c(0.2,0,0,0), unit = "cm"))+
                                              #scale_x_continuous(limits = c(1,filenum+0.01),breaks = NULL,minor_breaks = NULL,expand = c(0.5/filenum,0.5/filenum))+
                                              scale_x_continuous(expand = c(0.5/length(labels(dhc_marker)),0.5/length(labels(dhc_marker))))+
                                              #coord_flip() + 
                                              scale_y_continuous(expand = c(0, 0))+
                                              mytheme+
                                              theme(panel.background=element_rect(fill='transparent', color='white'))+
                                              theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())    
                        }else{
                          ddplot_cluster<-gblank
                          ddplot_marker<-gblank
                               }
                if(Rowv==F){ ddplot_cluster<-gblank }else{
                              #this_melt_htdata[,"clusters"]<-factor(this_melt_htdata[,"clusters"],levels =rev(labels(dhc_cluster)),ordered = TRUE )
                             }
                  
                if(Colv==F){ ddplot_marker<-gblank }
                gheatmap = ggplot(this_melt_htdata, aes_string(x="markers", y="clusters", fill="expression"))+
                            theme(plot.margin=unit(c(0,0,0,0), unit = "cm"))+
                            geom_tile(color="grey60", size=1)+
                            #mytheme+
                            scale_x_discrete(expand = c(0,0))+
                            theme(axis.text.x = element_text(angle=-90,hjust =0,vjust = 0.5))+
                            theme(axis.ticks = element_blank())+
                            theme(axis.text  = element_text(size = 20,color = "black"))+
                            theme(axis.title = element_blank())+
                            theme(legend.title =  element_blank(),legend.text = element_text(size = 20))+
                            guides(fill = guide_colorbar(barwidth  = unit(30/72,"inches"),
                                                         barheight = unit((length(labels(dhc_cluster)))*30/72*0.7, "inches"), ##图例的高度
                                                         raster = T ))+ ## 如果为TRUE，则颜色条将呈现为栅格对象。 如果为FALSE，则颜色条呈现为一组矩形
                            scale_y_discrete(position = "right",expand = c(0,0))+
                            scale_fill_gradientn(colours =expression_color_set(100),limits=c(0,max_expr))
                
                gheatmap2<-gheatmap+theme(legend.position = "none")
                gheatmap_legend<-as_ggplot(get_legend(gheatmap))
                gheatmap<-gheatmap+
                          theme(axis.text=element_text(size=50))
                gheatmap_x_axis<-as_ggplot(get_x_axis(gheatmap))
                gheatmap_y_axis<-as_ggplot(get_y_axis(gheatmap,position ="right"))

                
                # cluster_abundance_barplot_x_axis<-as_ggplot(get_x_axis(cluster_abundance_barplot))
            
                
                
                nmarkers<-length(labels(dhc_marker))
                nclusters<-length(labels(dhc_cluster))
                # plot_width<-(nmarkers+1)*(1+0.25+9/nmarkers)*30/72
                # plot_height<-(nclusters+3)*(1+5/nclusters)*30/72
                plot_width<-(nmarkers+15)*30/72
                plot_height<-(nclusters+5)*30/72

                
                # install.packages("egg")
                library(egg)
                pdf(file=paste0("./",output_dir,"/",this_group,"_heatmap.pdf"),width=plot_width,height=plot_height)
                
                if(show_abundance_barplot){
                        plot2<-ggarrange(gblank,ddplot_marker,gblank,gblank,
                                         ddplot_cluster,gheatmap2,cluster_abundance_barplot,gheatmap_legend,
                                         ncol = 4,
                                         widths = c(5/nmarkers,1,6/nmarkers,4/nmarkers),
                                         heights = c(5/nclusters,1))
                }else{
                        plot2<-ggarrange(gblank,ddplot_marker,gblank,
                                         ddplot_cluster,gheatmap2,gheatmap_legend,
                                         ncol = 3,
                                         widths = c(5/nmarkers,1,4/nmarkers),
                                         heights = c(5/nclusters,1)) 
                        
                       }

                dev.off()

                cat("Exported: Heatmap(s)\n")
              
            #2生成density-plot
            
            if(output_density_plot){
                
              
            
              #定义函数
              if_continuous<-function(x){ 
                
                                          if(substring(x,nchar(x),nchar(x))!="s"){
                                            
                                          output<-data.frame(continuous=T,
                                                             scale_fun=paste0("gradientn(colours=",x,"_color_set(100),limits=c(0,max_expr))"),
                                                             aes_name=x,
                                                             ncolor=100)}else{
                                           aes_name=substring(x,1,nchar(x)-1)
                                           ncolor=length(unique(this_dstplot_data_gather[,x]))
                                           output<-data.frame(continuous=F,
                                                              scale_fun=paste0("manual(values=",aes_name,"_color_set(",ncolor,"))"),
                                                              aes_name=aes_name,
                                                              ncolor=ncolor)
                                                }
                                          return(output)
                                        }                                                      
              
           
              this_dstplot_data_gather[,"clusters"]<-factor(this_dstplot_data_gather[,"clusters"],levels =rev(levels(this_melt_htdata[,"clusters"])),ordered = TRUE )
              this_melt_htdata[,"clusters"]<-factor(this_melt_htdata[,"clusters"],levels =rev(levels(this_melt_htdata[,"clusters"])),ordered = TRUE )
              
               plot1<-ggplot()+mytheme
               if(density_plot_fill_by!=density_bkgd_fill_by ){
                      plot1<-plot1+theme(panel.background=element_rect(fill=density_bkgd_fill_by))}else{
            
                      plot1<-plot1+geom_rect(data = this_melt_htdata,aes_string(fill=density_bkgd_fill_by),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha=0.8)
               }
               plot1<-plot1+
                      #geom_rect(data = this_melt_htdata,aes_string(fill="expression"),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha=0.8)+
                      geom_density(data=this_dstplot_data_gather,aes_string(x= "exprs"),adjust =1,size=density_line_size,alpha=density_fill_alpha)
                      #设置坐标轴
                      if(density_plot_fill_by!="none") {plot1<-plot1+aes_string(fill=density_plot_fill_by)
                                                        fill_formula<-paste0("plot1<-plot1+scale_fill_",if_continuous(density_plot_fill_by)$scale_fun)
                                                        eval(parse(text = fill_formula))      
                                                       }
                      if(density_plot_color_by!="none"){plot1<-plot1+aes_string(color=density_plot_color_by)
                                                        colorsets<-if_continuous(density_plot_color_by)
                                                        color_formula<-paste0("plot1<-plot1+scale_color_",if_continuous(density_plot_color_by)$scale_fun)
                                                        eval(parse(text = color_formula)) 
                      }
               
            
               
               if(!is.null(density_plot_xlim) & density_plot_free_x==F){
                 plot1<-plot1+
                   scale_x_continuous(limits=density_plot_xlim)}
               
               if(density_plot_free_x==F){
                 plot1<-plot1+
                   facet_wrap(as.formula(paste0("clusters"," ~markers") ),scales = "free_y",ncol = length(heatmap_ID))}
               
               if(density_plot_free_x==T){
                 plot1<-plot1+
                   facet_wrap(as.formula(paste0("clusters"," ~markers") ),scales = "free",ncol = length(heatmap_ID))}
               
               plot1<-plot1+theme(legend.title = element_text(size=50),
                                  legend.text = element_text(size=50),
                                  legend.key.width  = unit(3.5*0.3,'inches'),
                                  legend.key.height = unit(3.5*0.5,'inches')
                                  )+
                            theme(strip.text = element_text(size = strip_label_size))
              
               plot1_legend<-as_ggplot(get_legend(plot1))
               plot1<-plot1+theme(legend.position = "none")
           
               
               #输出facet结果
               col_plot<-(nmarkers+7)*3.5
               row_plot<-(nclusters+4.5)*3.5
               
               
               ddplot_marker<-ddplot_marker+
                              geom_segment(data=segment(ddplot_marker),aes(x = x, y = y, xend = xend, yend = yend),size=3)+
                              theme(plot.margin=unit(c(1,0,0,0), unit = "cm"))
                 

               ddplot_cluster<-ddplot_cluster+
                              geom_segment(data=segment(ddplot_cluster),aes(x = x, y = y, xend = xend, yend = yend),size=3)+ 
                              theme(plot.margin=unit(c(0,0,0,1), unit = "cm"))
                 
               if(Rowv==F || facet_by_groups==T){ ddplot_cluster<-gblank }
               if(Colv==F || facet_by_groups==T){ ddplot_marker<-gblank }

               plot3<-ggarrange(gblank,ddplot_marker,gblank,gblank,
                                ddplot_cluster,plot1,gheatmap_y_axis,plot1_legend,
                                gblank,gheatmap_x_axis,gblank,gblank,ncol = 4,
                                widths = c(3/nmarkers,1,0.5/nmarkers,3/nmarkers),
                                heights = c(2.5/nclusters,1,2/nclusters))
               dev.off()
               pdf(file=paste0("./",output_dir,"/",this_group,"_Density_plots.pdf"),width = col_plot,height = row_plot)
               multiplot(plot3)


               dev.off()

               cat(paste0("Exported: ","Density plot",this_group,".pdf\n"))
            }

        }
  
}



#' @title draw_expr_heatmap
#'
#' @description
#' 生成cluster的Heatmap
#'
#'@param xdata                data.frame,或者matrix，cluster expression matrix
#'@param trans_method         数据转化方式，有五种：
#'                            "CytofAsinh"，CytofAsinh转换所有表达数据；
#'                            "simpleAsinh",simpleAsinh转换所有表达数据,这是默认设置；
#'                            "0_to_Max"，所有Marker的信号强度除以最大值，线性转换，通过除以各通道信号的最大值把数值scale到0~1；
#'                            "Min_to_Max"，线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异；
#'                            "simpleAsinh_0_to_Max" 先Arcsinh转换，然后除以各通道信号的最大值把数值scale到0~1范围内；
#'                            "simpleAsinh_Min_to_Max" 先Arcsinh转换，然后最小值转换成0，最大值转换成1，最大限度展现population的表达差异；
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
  
  if(trans_method=="CytofAsinh") trans_fun=cytofAsinh
  if(trans_method=="simpleAsinh") trans_fun=simpleAsinh
  if(trans_method=="0_to_Max") trans_fun=normalize_Zero_Max
  if(trans_method=="Min_to_Max") trans_fun=normalize_min_max
  if(trans_method=="simpleAsinh_0_to_Max") trans_fun=simpleAsinh_Zero_Max
  if(trans_method=="simpleAsinh_Min_to_Max") trans_fun=simpleAsinh_Min_to_Max

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
           filename = paste0("./",output_dir,"/Cluster_expression_heatmap (",trans_method,").pdf")
           )
  

}



### 数据预处理
#定义转化函数



## 所有Marker的信号强度除以最大值
##线性转换，通过除以各通道信号的最大值把数值scale到0~1

### 数据预处理
#定义转化函数

normalize_Zero_Max<-function(value,b=0)
{
  if(!all(value==0)){
    value<-value/max(value,b)
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
  }#else{message("数据全部是0，无法进行0 to max 转化\n")}
  return(value)
}


## 所有Marker的取Arcsinh后，除以最大值
## Arcsinh转换，然后除以各通道信号的最大值把数值scale到0~1范围内

### 数据预处理
#定义转化函数

simpleAsinh_Zero_Max<-function(value, cofactor = 5,b=0) {
  value <- value / cofactor
  value <- asinh(value)
  if(!all(value==0)){
    value<-value/max(value,b)
  }#else{message("数据全部是0，无法进行0 to max 转化\n")}
  return(value)
}



## Arcsinh转换，然后进行线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异

### 数据预处理
#定义转化函数

simpleAsinh_Min_to_Max<-function(value, cofactor = 5) {
  value <- value / cofactor
  value <- asinh(value)
  if(!all(value==0)){
    value<-(value-min(value))/(max(value)-min(value))
  }#else{message("数据全部是0，无法进行0 to max 转化\n")}
  return(value)
}
