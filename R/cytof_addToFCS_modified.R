
#'Notice: The following functions are modified or adapted from cytofkit, please cite original articles if they were used in your pepline.
#'For details, please enter citation("cytofkit") in console lines.
#'
#'

#' @title cytof_addToFCS_modified
#'
#' @description
#' 本函数修改自cytofkit::cytof_addToFCS,根据数据分析的需要对原函数进行了简化和修改。功能是将降维和聚类分析的数据整合到FCS文件中。
#' 请注意引用Cytofkit的原始文献，详情运行 citation("cytofkit") 查看。
#'
#' @param data             一个datafram，包含想要写入FCS文件的各个参数
#' @param rawFCSdir        原始FCS文件所在的目录
#' @param origSampNames    原始的FCS文件名称
#' @param analyzedFCSdir   生成FCS文件存放的目录名
#' @param newHeader        一个字符串，做为新生成文件名的前缀，方便文件识别
#'
#'
#' @return                 返回指定数量的颜色序列
#' @export



cytof_addToFCS_modified <- function (data, rawFCSdir, origSampNames = NULL, analyzedFCSdir,newHeader = "")
{


  if (!dir.exists(analyzedFCSdir)) {
    dir.create(analyzedFCSdir)
  }


if (!is.null(data)) {
  colnames(data) <- colnames(data)
  to_add<-data
}


  addColNames <- colnames(to_add)
  sample <- unique(sub("_[0-9]*$", "", row.names(to_add)))
  for (i in 1:length(sample)) {


    if (!is.null(origSampNames)) {
      fn <- paste0(rawFCSdir, .Platform$file.sep, origSampNames[i],
                   ".fcs")
    }    else {
      fn <- paste0(rawFCSdir, .Platform$file.sep, sample[i],
                   ".fcs")
    }
    if (!file.exists(fn)) {
      message(paste("Cannot find raw FCS file:", fn))
      return(NULL)
    }
#   cat("Save to file:", fn, "\n")
    cat("Save to file: ",paste0(analyzedFCSdir,"/",newHeader, sample[i], ".fcs"),"\n")

    fcs <- read.FCS(fn, transformation = FALSE)
    pattern <- paste(sample[i], "_", sep = "")
    to_add_i <- as.data.frame(to_add[grep(pattern, row.names(to_add),
                                          fixed = TRUE), ,drop=F])       #增加drop以适应增加单个通道的需要


    m <- regexpr("_[0-9]*$", row.names(to_add_i))
    cellNo_i <- as.integer(substring(regmatches(row.names(to_add_i),
                                                m), 2))
    sub_exprs <- fcs@exprs[cellNo_i, ]
    params <- parameters(fcs)
    pd <- pData(params)


    keyval <- keyword(fcs)
    for (j in 1:length(addColNames)) {
      if (addColNames[j] %in% colnames(fcs@exprs)) {
        addColNames[j] <- paste(addColNames[j], "_new",
                                sep = "")
      }
      addColName <- addColNames[j]
      channel_number <- nrow(pd) + 1
      channel_id <- paste("$P", channel_number, sep = "")
      channel_name <- addColName
      minRange <- ceiling(min(to_add_i[[j]]))
      maxRange <- ceiling(max(to_add_i[[j]]))
      channel_range <- maxRange - minRange+1                               #修改一个计算错误
      #plist <- matrix(c(channel_name, "<NA>", channel_range,
      #                  minRange, maxRange))
      plist <- matrix(c(channel_name, channel_name, channel_range,
                        minRange, maxRange))


      rownames(plist) <- c("name", "desc", "range", "minRange",
                           "maxRange")
      colnames(plist) <- c(channel_id)
      pd <- rbind(pd, t(plist))
      out_col_names <- colnames(sub_exprs)
      sub_exprs <- cbind(sub_exprs, to_add_i[[j]])
      colnames(sub_exprs) <- c(out_col_names, addColName)
      keyval <- lapply(keyval, function(x) {
        if (class(x) == "character") {
          gsub("\\", "", x, fixed = TRUE)
        }  else {
          x
        }
      })
      keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
      keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
      keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
      keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
      keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name  #added to show shortname of parameters
      keyval[[paste("P", channel_number, "BS", sep = "")]] <- 0
      keyval[[paste("P", channel_number, "MS", sep = "")]] <- 0
      keyval[[paste("P", channel_number, "DISPLAY", sep = "")]] <- "LIN"
    }

    pData(params) <- pd

    out_frame <- flowFrame(exprs = sub_exprs, parameters = params,description = keyval)
    suppressWarnings(write.FCS(out_frame, paste0(analyzedFCSdir,"/",newHeader, sample[i], ".fcs")))


    }
}





#' @title cytofAsinh 定义数据转化函数
#' @description
#' 本函数来自cytofkit::cytofAsinh，用于Cytof数据的转化，实现去背景（减去1然后将负值随机化），和数据Arcsinh转化。
#' @param value      一个vector，待转化的数据
#' @param cofactor   一个数字做为辅助因子，vector中的数据要先除以cofactor，然后再进行Arcsin转化
#' @export

cytofAsinh <- function(value, cofactor = 5) {
  value <- value-1
  loID <- which(value < 0)
  if(length(loID) > 0) {   value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)}
  value <- value / cofactor
  value <- asinh(value)
  return(value)
}
