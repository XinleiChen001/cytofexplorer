
#'The following functions are modified from cytofkit, please cite cytofik article if they were used in your pepline.

#'@title cytof_addToFCS_modified （modified from cytofkit::cytof_addToFCS）
#'
#'@description 
#'在已有FCS文件中增加额外的通道
#' @param data The new data matrix to be added in.
#' @param rawFCSdir The directory containing the original fcs files.
#' @param origSampNames Vector of original names of samples, if samples were renamed.
#' @param analyzedFCSdir The directory to store the new fcs files.
#' @param transformed_cols The column name of the dimension transformed data in \code{data}.
#' @param cluster_cols The column name of the cluster data in \code{data}.
#' @param clusterIDs Table of cell cluster IDs for each clustering method
#' @param specifySampleNames Used only if sample names differ from those in raw fcs.
#' @param inLgclTrans If \verb{TRUE}, apply the inverse lgcl transformation to the the cluster data before saving
#' @param newHeader  生成文件名前缀
#' 
#' @export
#'@return none

cytof_addToFCS_modified <- function (data, rawFCSdir, origSampNames = NULL, analyzedFCSdir, 
          transformed_cols = NULL, cluster_cols = NULL, 
          clusterIDs = data, specifySampleNames = NULL, inLgclTrans = FALSE,newHeader = "cytofkit_") 
{



  lgcl <- logicleTransform(w = 0.1, t = 4000, m = 4.5, a = 0)
  ilgcl <- inverseLogicleTransform(trans = lgcl)
  if (!dir.exists(analyzedFCSdir)) {
    dir.create(analyzedFCSdir)
  }
  if (!is.null(transformed_cols)) {
    transformed <- data[, transformed_cols, drop = FALSE]
    N_transformed <- apply(transformed, 2, function(x) ((x - 
                                                           min(x))/(max(x) - min(x))) * 4.4)
    R_N_transformed <- apply(N_transformed, 2, ilgcl)
    R_N_transformed_l <- apply(transformed, 2, function(x) (x - 
                                                              min(x)) + 0.1)
    colnames(R_N_transformed_l) <- paste0(colnames(R_N_transformed_l), 
                                          "_linear", sep = "")
    R_N_transformed <- cbind(R_N_transformed, R_N_transformed_l)
    row.names(R_N_transformed) <- row.names(data)
  }
  if (!is.null(specifySampleNames)) {
    originalSample <- specifySampleNames
  }
  if (!is.null(cluster_cols)) {
    if (inLgclTrans) {
      for (i in 1:length(cluster_cols)) {
        cCol <- data[, cluster_cols[i]]
        clust_cor_1 <- as.numeric(cCol)%%10
        clust_cor_2 <- floor(as.numeric(cCol)/10)
        clust_cor_1 <- clust_cor_1 + runif(length(clust_cor_1), 
                                           0, 0.2)
        clust_cor_2 <- clust_cor_2 + runif(length(clust_cor_2), 
                                           0, 0.2)
        cluster_cor12 <- cbind(clust_cor_1, clust_cor_2)
        N_clust_cor <- apply(cluster_cor12, 2, function(x) ((x - 
                                                               min(x))/(max(x) - min(x))) * 4.4)
        clust_cor <- apply(N_clust_cor, 2, ilgcl)
        colnames(clust_cor) <- paste(cluster_cols[i], 
                                     c("cor_1", "cor_2"), sep = "_")
        if (i == 1) {
          R_N_clust_cor <- clust_cor
        }        else {
          R_N_clust_cor <- cbind(R_N_clust_cor, clust_cor)
        }
      }
    }    else {
      R_N_clust_cor <- data[, cluster_cols, drop = FALSE]
    }
    row.names(R_N_clust_cor) <- row.names(data)
  }
  if ((!is.null(transformed_cols)) && (!is.null(cluster_cols))) {
    to_add <- cbind(R_N_transformed, R_N_clust_cor)
  }  else if (!is.null(transformed_cols)) {
    to_add <- R_N_transformed
  }  else if (!is.null(cluster_cols)) {
    to_add <- R_N_clust_cor
  }
  if (!is.null(clusterIDs)) {

    
    colnames(clusterIDs) <- colnames(clusterIDs)

    to_add<-clusterIDs
    
  }
  
  #########
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
                                          fixed = TRUE), ,drop=F])  #增加drop以适应增加单个通道的需要
    
    
    m <- regexpr("_[0-9]*$", row.names(to_add_i))
    cellNo_i <- as.integer(substring(regmatches(row.names(to_add_i), 
                                                m), 2))
    sub_exprs <- fcs@exprs[cellNo_i, ]
    params <- parameters(fcs)
    pd <- pData(params)
    
    #ori_channel_number <- nrow(pd) + 1 #record original channel_numer
    
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
      channel_range <- maxRange - minRange+1  # make a modification :+1
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
      keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name  ###added to show shortname of parameters
      keyval[[paste("P", channel_number, "BS", sep = "")]] <- 0
      keyval[[paste("P", channel_number, "MS", sep = "")]] <- 0
      keyval[[paste("P", channel_number, "DISPLAY", sep = "")]] <- "LIN"
    }
   
    pData(params) <- pd
    out_frame <- flowFrame(exprs = sub_exprs, parameters = params,description = keyval)
    suppressWarnings(write.FCS(out_frame, paste0(analyzedFCSdir,"/",newHeader, sample[i], ".fcs")))

    
    }
}
