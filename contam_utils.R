load_intSiteCaller_data <- function(sampleName, dataType, dataDir){
  content <- list.files(path = paste(dataDir, sampleName, sep = "/"))
  isThere <- any(grepl(dataType, content))
  
  if(isThere){
    file <- grep(dataType, content, value = TRUE)
    filePath <- paste(dataDir, sampleName, file, sep = "/")
    load(filePath)
  }else{
    message(paste0(
      "Could not find dataType in primary analysis directory for\n", 
      sampleName))
  }
  
  if(dataType == "allSites" & isThere){
    allSites
  }else if(dataType == "multihitData" & isThere){
    multihitData
  }else if(dataType == "keys" & isThere){
    keys
  }else if(dataType == "sites.final" & isThere){
    sites.final
  }else if(dataType == "primerIDData" & isThere){
    primerIDs
  }else if(isThere){
    stop("dataType not supported by this function. Check dataType.")
  }
}

assign_groups <- function(sites, start = TRUE){
  if(is.null(names(sites))) names(sites) <- seq(1:length(sites))
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

find_primerID_crossover <- function(sites){
  if(!any(grepl("posID", names(mcols(sites))))) sites$posID <- generate_posID(sites)
  sites.df <- as.data.frame(sites)
  IDs <- table(sites.df$primerID)
  two_or_more <- sites.df %>%
    filter(primerID %in% names(IDs[IDs >= 2])) %>%
    select(., patient, specimen, sampleName, primerID, posID) %>%
    group_by(primerID, posID) %>% 
    summarize(
      patients = n_distinct(patient),
      specimens = n_distinct(specimen), 
      sampleNames = n_distinct(sampleName)) %>%
    filter(sampleNames >= 2) %>%
    ungroup(.)
  
  sites_w_shared_id <- sites[sites$primerID %in% two_or_more$primerID]
  sites_w_shared_id <- sites_w_shared_id[
    sites_w_shared_id$posID %in% two_or_more$posID
  ]
  
  ids_shared_by_reps <- two_or_more %>%
    filter(., patients == 1) %>%
    filter(., specimens == 1)
  
  ids_shared_by_spec <- two_or_more %>%
    filter(., patients == 1) %>%
    filter(., specimens > 1)
  
  ids_shared_by_pats <- two_or_more %>%
    filter(., patients >1)
  
  sites_shared_by_reps <- sites[sites$primerID %in% ids_shared_by_reps$primerID]
  sites_shared_by_spec <- sites[sites$primerID %in% ids_shared_by_spec$primerID]
  sites_shared_by_pats <- sites[sites$primerID %in% ids_shared_by_pats$primerID]
  
  list_of_ids <- list(ids_shared_by_reps, ids_shared_by_spec, ids_shared_by_pats)
  id_names <- c("ids_shared_by_reps", "ids_shared_by_spec", "ids_shared_by_pats")
  names(list_of_ids) <- id_names
  list_of_sites <- list(sites_shared_by_reps, sites_shared_by_spec, sites_shared_by_pats)
  site_names <- c("sites_shared_by_reps", "sites_shared_by_spec", "sites_shared_by_pats")
  names(list_of_sites) <- site_names
  output <- list(list_of_ids, list_of_sites)
  names(output) <- c("primerIDs", "sites_sharing_primerIDs")
  output
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
