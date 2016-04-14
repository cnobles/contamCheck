options(stringsAsFactors = FALSE)

#Check and load packages
rPackages <- c("RMySQL", "devtools", "argparse", "plyr", "dplyr")
stopifnot(all(sapply(rPackages, require, 
                     character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))
source_url("https://raw.githubusercontent.com/cnobles/cloneTracker/master/cloneTracker.SOURCE_ME.R")

setArguments <- function(){
  parser <- ArgumentParser(description = "Determine all integration site present across patients.")
  parser$add_argument("-s", "--specimen_database", default = "hiv_specimen.database", 
                      help = "Group to use for specimen data.")
  parser$add_argument("-i", "--intsites_database", default = "hiv_intsites.database", 
                      help = "Reference genome to use in intSiteCaller.")
  
  arguments <- parser$parse_args()
  arguments
}

arguments <- setArguments()

specimenDatabase <- arguments$specimen_database
intSitesDatabase <- arguments$intsites_database

if(specimenDatabase == "hiv_specimen.database"){
  specimenTable <- "nobles.hivsp"
}else{
  specimenTable <- "specimen_management.gtsp"
}

##Set up working directory and load sampleInfo
dataDir <- getwd()
splitDir <- unlist(strsplit(dataDir, split = "/"))
runName <- splitDir[length(splitDir)]
sampleInfoFile <- grep("sampleInfo.tsv", list.files(path = dataDir), value = TRUE)
sampleInfo <- read.csv(paste(dataDir,sampleInfoFile, sep = "/"), sep = "\t")
sampleInfo$specimen <- sapply(strsplit(
  as.character(sampleInfo$alias), 
  "-"), "[", 1)

if(any(grepl("UninfectedControl", sampleInfo$specimen))){
  uncPos <- grep("UninfectedControl", sampleInfo$specimen)
  sampleInfo[uncPos, "specimen"] <- "UNC"
  rm(uncPos)
}

if(any(grepl("NoTemplateControl", sampleInfo$specimen))){
  ntcPos <- grep("NoTemplateControl", sampleInfo$specimen)
  sampleInfo[ntcPos, "specimen"] <- "NTC"
  rm(ntcPos)
}

if(nrow(sampleInfo) > 0){
  message(paste0("Loaded information for ", 
                 nrow(sampleInfo), 
                 " samples from ", 
                 sampleInfoFile, 
                 "."))
}else{
  message(paste0("No data loaded from ", sampleInfoFile, "."))
}

#Get patient information from specimen management database
specimen_query <- unique(sampleInfo[!sampleInfo$specimen == "UNC" &
                                      !sampleInfo$specimen == "NTC", "specimen"])

message("Querying patient information for:")
message(list(specimen_query))
               
junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
dbConn <- dbConnect(MySQL(), group = specimenDatabase)  
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1) 

if(specimenDatabase == "hiv_specimen.database"){
  string <- paste(specimen_query, collapse="' OR parentAlias = '")
  query_request <- paste0(" WHERE parentAlias = '", string, "'")
  sql <- paste0("SELECT parentAlias, patient FROM ", specimenTable,
                query_request)
}else{
  string <- paste(specimen_query, collapse="' OR SpecimenAccNum = '")
  query_request <- paste0(" WHERE SpecimenAccNum = '", string, "'")
  sql <- paste0("SELECT SpecimenAccNum, Patient FROM ", specimenTable,
                query_request)
}

patientInfo <- suppressWarnings( dbGetQuery(dbConn, sql) ) 
dbDisconnect(dbConn) #Disconnect from specimen_management after query
rm(junk, dbConn, string, query_request)

if(specimenDatabase == "hiv_specimen.database"){
  isThere <- unique(patientInfo$parentAlias)
}else{
  isThere <- unique(patientInfo$SpecimenAccNum)
}
stopifnot( (NA %in% match(specimen_query, isThere)) == FALSE)
message("Received patient information for:")
message(list(isThere))

if(specimenDatabase == "hiv_specimen.database"){
  if(any(grepl("UNC", sampleInfo$specimen))){
    patientInfo <- rbind(patientInfo, data.frame(
      "parentAlias" = "UNC", "patient" = "UNC")
      )
  }
  if(any(grepl("NTC", sampleInfo$specimen))){
    patientInfo <- rbind(patientInfo, data.frame(
      "parentAlias" = "NTC", "patient" = "NTC")
      )
  }
  sampleInfo$patient <- patientInfo[match(sampleInfo$specimen, patientInfo$parentAlias), "patient"]
}else{
  if(any(grepl("UNC", sampleInfo$specimen))){
    patientInfo <- rbind(patientInfo, data.frame(
      "SpecimenAccNum" = "UNC", "Patient" = "UNC")
    )
  }
  if(any(grepl("NTC", sampleInfo$specimen))){
    patientInfo <- rbind(patientInfo, data.frame(
      "SpecimenAccNum" = "NTC", "Patient" = "NTC")
    )
  }
  sampleInfo$patient <- patientInfo[match(sampleInfo$specimen, patientInfo$SpecimenAccNum), "Patient"]
}

#Load integration site data from intSiteCaller output
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

message("Loading intSiteCaller information.")

allSites <- lapply(sampleInfo$alias, function(sampleName){
  load_intSiteCaller_data(sampleName, "allSites", dataDir)
})
names(allSites) <- sampleInfo$alias

if( length(allSites) == nrow(sampleInfo)){
  message("intSiteCaller data loaded.")
}

#Process data for contamination
message(paste0("Standardizing ", 
               sum(sapply(allSites, length)), 
               " unique reads."))
allSites <- allSites[sapply(allSites, class) == "GRanges"]
std.allSites <- standardize_intsites(unlist(GRangesList(allSites)), 
                                     standardize_breakpoints = TRUE)
derep.allSites <- GRangesList(lapply(
  split(std.allSites, std.allSites$sampleName), function(x){
  remove_repeats(x, c("sampleName"))
}))
message(paste0("Identified ", 
               sum(sapply(derep.allSites, length)), 
               " total sonic fragments."))

uniq.allSites <- unlist(derep.allSites)
uniq.allSites$specimen <- sampleInfo[
  match(uniq.allSites$sampleName, sampleInfo$alias), "specimen"]
uniq.allSites$patient <- sampleInfo[
  match(uniq.allSites$sampleName, sampleInfo$alias), "patient"]
uniq.allSites <- split(uniq.allSites, uniq.allSites$patient)
message("Checking for contamination between:")
message(list(names(uniq.allSites)))
possible.contam <- track_clones(uniq.allSites, gap = 0L, track.origin = FALSE)

possible.contam.gr <- unlist(possible.contam)

#Save data
message("Saving data.")
save(possible.contam, 
     file = paste0(dataDir, "/", runName, ".possible.contam.RData"))

if(specimenDatabase == "hiv_specimen.database"){
  stats <- data.frame(
    "amount" = c(length(unique(
      sampleInfo[grep("HIVSP", sampleInfo$specimen), "patient"]
    )),
    length(unique(
      sampleInfo[grep("HIVSP", sampleInfo$specimen, invert = TRUE), "patient"]
    )),
    length(unique(
      flank(granges(unlist(uniq.allSites)), -1, start = TRUE)
    )),
    length(unique(
      flank(granges(possible.contam.gr), -1, start = TRUE)
    )),
    length(unique(granges(possible.contam.gr))),
    length(possible.contam.gr),
    length(unique(possible.contam.gr$patient)),
    length(unique(
      possible.contam.gr[grep("UNC", possible.contam.gr$patient)]$posid
    )),
    length(unique(
      possible.contam.gr[grep("NTC", possible.contam.gr$patient)]$posid
    ))),
    "stat" = c("patientCount", "controlCount", "sitesConsidered", "uniqSitesXOver",
               "uniqSonicFragsXOver", "totalSonicFragsXOver", "subjectsInvolved",
               "UninfectedXOver", "NoTemplateXOver")
  )
}else{
  stats <- data.frame(
    "amount" = c(length(unique(
        sampleInfo[grep("GTSP", sampleInfo$specimen), "patient"]
      )),
      length(unique(
        sampleInfo[grep("GTSP", sampleInfo$specimen, invert = TRUE), "patient"]
      )),
      length(unique(
        flank(granges(unlist(uniq.allSites)), -1, start = TRUE)
      )),
      length(unique(
        flank(granges(possible.contam.gr), -1, start = TRUE)
      )),
      length(unique(granges(possible.contam.gr))),
      length(possible.contam.gr),
      length(unique(possible.contam.gr$patient)),
      length(unique(
        possible.contam.gr[grep("UNC", possible.contam.gr$patient)]$posid
      )),
      length(unique(
        possible.contam.gr[grep("NTC", possible.contam.gr$patient)]$posid
    ))),
    "stat" = c("patientCount", "controlCount", "sitesConsidered", "uniqSitesXOver",
               "uniqSonicFragsXOver", "totalSonicFragsXOver", "subjectsInvolved",
               "UninfectedXOver", "NoTemplateXOver")
  )
}
write.table(stats, file = paste0(dataDir, "/", runName, ".contam.stats.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
message("Check complete.")


