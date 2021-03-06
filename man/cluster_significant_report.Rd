% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cluster_Abundance.R
\name{cluster_significant_report}
\alias{cluster_significant_report}
\title{cluster统计中细胞数量的质控（quality control）}
\usage{
cluster_significant_report(
  cluster_stat,
  count_threshold = 50,
  percent_threshold = 0.5,
  group_to_show = NULL,
  output_dir = "Aundance_significant_report",
  stat.paired = NULL,
  stat.method = NULL,
  set.equal.var = NULL,
  conf.level = NULL,
  p.adjust.method = NULL
)
}
\arguments{
\item{cluster_stat}{stat_by_cluster()函数生成的list数据}

\item{count_threshold}{数字，指定统计要求的细胞数，默认50}

\item{percent_threshold}{数字，指定统计要求的细胞数，默认0.5，即0.5%}

\item{group_to_show}{指定要统计的condition（major_cond中的一个或者多个），如不设置直接统计全部；}

\item{output_dir}{输出数据文件夹名称

如使用统计默认参数，以下参数无需设置}

\item{stat.paired}{逻辑变量，TRUE或者FALSE，表示是否是成对检验}

\item{stat.method}{一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。}

\item{set.equal.var}{一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；}

\item{conf.level}{显著性水平，默认0.95,（即p<0.05为显著）}

\item{p.adjust.method}{对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"}
}
\value{
none
}
\description{
对cluster进行统计分析时，过少的细胞数有可能造成统计结果的不准确，本函数统计分析方法分析每个cluster的绝对数量是否符合要求
结果采用了与SPADEvizR相近的风格，标注有“cite from SPADEVizR”或者“Adapted from SPADEVizR” ，使用时请注意在文章中引用原始文献。
}
\details{
说明:

SPADEViZR工具包可以对SPADE聚类结果进行可视化的统计分析，发现在组成比例上显著性差异的亚群。详细细节请参照其github主页
链接：https://github.com/tchitchek-lab/SPADEVizR
文献：SPADEVizR: an R package for visualization, analysis and integration of SPADE results. Gautreau G etc.
在本package中我们保持了与其分析结果的风格，并加以扩展，使其同时可以适用于与SPADE以外的分析结果的分析，并可以同时输出padjust前后的结果。所有风格与
SPADEvizR相近的图均标注有“cite from SPADEVizR”或者“Adapted from SPADEVizR” ，使用时请注意在文章中引用原始文献。

Notice:
SPADEViZR R tool package provides a visual analysis tool of the SPADE clustering results and finds population that are significantly different in composition. Please refer to its github homepage for details.
Link: https://github.com/tchitchek-lab/SPADEVizR
Published Article:SPADEVizR: an R package for visualization, analysis and integration of SPADE results. Gautreau G etc.
In this package, we adopt the similar style of its analysis results and extend it so that it can be applied to the clustering analysis other than SPADE,
and can output the results before and after p-adjust. All figures similar to SPADEvizR are marked with "cite from SPADEVizR" or "Adapted from SPADEVizR".
Please pay attention to it, and cite the original article in your reference.
}
