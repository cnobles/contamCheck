load_intSiteCaller_data <- function(sampleName, dataType, dataDir){
  content <- list.files(path = paste(dataDir, sampleName, sep = "/"))
  isThere <- any(grepl(dataType, content))
  
  if(isThere){
    file <- grep(dataType, content, value = TRUE)
    filePath <- paste(dataDir, sampleName, file, sep = "/")
    load(filePath)
  }
  
  if(dataType == "allSites" & isThere){
    allSites
  }else if(dataType == "multihitData" & isThere){
    multihitData
  }else if(dataType == "keys" & isThere){
    keys
  }else if(dataType == "sites.final" & isThere){
    sites.final
  }else if(isThere){
    stop("dataType not supported by this function. Check dataType.")
  }
}

assign_groups <- function(sites, start = TRUE){
  origin_order <- names(sites)
  sites_flank <- flank(sites, -1, start = start)
  sites_red <- reduce(sites_flank, min.gapwidth = 5L, with.revmap = TRUE)
  revmap <- sites_red$revmap
  groups <- as.character(Rle(
    values = sapply(revmap, "[", 1), 
    lengths = sapply(revmap, length)
  ))
  sites <- sites[unlist(revmap)]
  sites$clusID <- groups
  sites <- sites[origin_order]
  as.numeric(sites$clusID)
}

group_sites <- function(sites){
  position_grp <- assign_groups(sites, start = TRUE)
  breakpoint_grp <- assign_groups(sites, start = FALSE)
  sites$clusID <- paste(sites$primerID, 
                        position_grp, 
                        breakpoint_grp,
                        sep = ":")
  sites
}

bsub <- function(queue="normal", cpus=1, maxmem=NULL, wait=NULL, jobName=NULL, logFile=NULL, command=NULL){
  stopifnot(!is.null(maxmem))
  stopifnot(!is.null(command))
  
  cmd <- paste0("bsub -q ", queue, " -n ", as.character(cpus), " -M ", maxmem)

  if(!is.null(wait)){
    LSF.VERSION <- system2("bsub", "-V", stdout=TRUE, stderr=TRUE)[1]
    if( grepl("openlava", LSF.VERSION, ignore.case=TRUE) ) {
      wait <- sub("done", "ended", wait)
    }
    cmd <- paste0(cmd, " -w \"", wait, "\"")
  }
  
  if(!is.null(jobName)) cmd <- paste0(cmd, " -J \"", jobName, "\"")
  if(!is.null(logFile)) cmd <- paste0(cmd, " -o ", logFile)
  
  cmd <- paste0(cmd, " ", command)
  message(cmd)
  system(cmd)
}
