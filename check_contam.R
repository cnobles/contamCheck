# Get commandline arguments ----------------------------------------------------
# -s + specimen database for information about each patient in the dataset
# -p to use primerID in checking methods
# -bsub to use the LSF queuing system
# -qsub to use the SGE queuing system (or qsub alternative)
# not supplying either -bsub or -qsub will process the script serially
args <- commandArgs(trailingOnly = FALSE)

codeDir <- dirname(sub("--file=","",grep("--file=", args, value = TRUE)))
specimenDatabase <- args[ grep("-d", args) + 1 ]
usePrimerID <- ifelse(any(grepl("-p", args)) == 1, TRUE, FALSE)
lsfhpc <- ifelse(any(grepl("-bsub", args)) == 1, TRUE, FALSE)
sgehpc <- ifelse(any(grepl("-qsub", args)) == 1, TRUE, FALSE)

if(specimenDatabase == "hiv_specimen.database"){
  specimenTable <- "nobles.hivsp"
}else{
  specimenTable <- "specimen_management.gtsp"
}

dataDir <- getwd()
splitDir <- unlist(strsplit(dataDir, split = "/"))
runName <- splitDir[length(splitDir)]

## Change runName if name is too long, such as from a MiSeq run
if(nchar(runName) > 25){
  flowcell <- unlist(strsplit(runName, split = "-"))
  flowcell <- flowcell[length(flowcell)]
  runName <- flowcell
}

# Print arguments to verify to the user that the correct input was given
message(paste0("Run Identifier: ", runName))
message(paste0("Data Directory: ", dataDir))
message(paste0("Code Directory: ", codeDir))
message(paste0("Specimen Database: ", specimenDatabase))
message(paste0("Use PrimerID: ", ifelse(usePrimerID, "Yes", "No")))
message(paste0("Use LSF HPC (bsub): ", ifelse(lsfhpc, "Yes", "No")))
message(paste0("Use SGE HPC (qsub): ", ifelse(sgehpc, "Yes", "No")))

if(all(lsfhpc, sgehpc)){stop("Select either 'qsub' or 'bsub'.")}

# Check and load dependances ---------------------------------------------------
dependancies <- c("RMySQL", "devtools", "plyr", "dplyr", 
                  "Biostrings", "GenomicRanges", "igraph")

null <- sapply(dependancies, function(x){
  suppressPackageStartupMessages(
    try(library(x, character.only = TRUE), silent = TRUE))
})

dependancies_present <- sapply(dependancies, function(package){
  package <- paste0("package:", package)
  logic <- package %in% search()
})

if(!any(dependancies_present)){
  Unloaded_Packages <- data.frame(package=as.character(dependancies), 
                                  loaded=dependancies_present)
  write.table(Unloaded_Packages, 
              file = "Unloaded_Packages.tsv", 
              quote = FALSE, 
              row.names = FALSE)
  stop("Load required packages. Check Unloaded_Packages.tsv for missing
       dependancies.")
}else{
  remove(dependancies, dependancies_present)
  message("Required packages loaded.")
}

source_url("https://raw.githubusercontent.com/cnobles/cloneTracker/master/cloneTracker.SOURCE_ME.R")
source(file.path(codeDir, "contam_utils.R"))

# Load sampleInfo --------------------------------------------------------------
sampleInfoFile <- grep("sampleInfo.tsv", list.files(path = dataDir), value = TRUE)
sampleInfo <- read.csv(paste(dataDir,sampleInfoFile, sep = "/"), sep = "\t")
sampleInfo$specimen <- sapply(strsplit(as.character(sampleInfo$alias), "-"), "[", 1)

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
  message(paste0("Loaded information for ", nrow(sampleInfo), " samples from ", 
                 sampleInfoFile, "."))
}else{
  message(paste0("No data loaded from ", sampleInfoFile, "."))
}

# Get patient information from specimen management database --------------------
specimen_query <- unique(sampleInfo[!sampleInfo$specimen == "UNC" &
                                      !sampleInfo$specimen == "NTC", "specimen"])

message("Querying patient information for:\n", paste(specimen_query, collapse = ", "))

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
message("Received patient information for:\n", paste(isThere, collapse = ", "))

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
  sampleInfo$patient <- patientInfo[
    match(sampleInfo$specimen, patientInfo$parentAlias), 
    "patient"]
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
  sampleInfo$patient <- patientInfo[
    match(sampleInfo$specimen, patientInfo$SpecimenAccNum), 
    "Patient"]
}

# Load integration site data from intSiteCaller output -------------------------
#!!!!!! Need to add primerID section into this under the primerID conditional !!!!!!#
message("Loading intSiteCaller information.")

allSites <- lapply(sampleInfo$alias, 
                   load_intSiteCaller_data,
                   dataType = "allSites", 
                   dataDir = dataDir)

names(allSites) <- sampleInfo$alias
if( length(allSites) == nrow(sampleInfo)) message("Unique sites data loaded.")
allSites <- allSites[sapply(allSites, length) > 0]

if(usePrimerID){
  primerIDs <- lapply(sampleInfo$alias,
                      load_intSiteCaller_data, 
                      dataType = "primerIDData", 
                      dataDir = dataDir)

  names(primerIDs) <- sampleInfo$alias
  if( length(primerIDs) == nrow(sampleInfo)) message("PrimerID data loaded.")
  primerIDs <- primerIDs[sapply(primerIDs, length) > 0]
}

#Process data for contamination
message(paste0("Standardizing ", sum(sapply(allSites, length)), " unique reads."))
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

if(usePrimerID){
  #Merge primerIDs with std.allSites based on read names
  primerIDs <- do.call(c, lapply(1:length(primerIDs), function(i){primerIDs[[i]]}))
  
  std.allSites$specimen <- sampleInfo[
    match(std.allSites$sampleName, sampleInfo$alias), "specimen"]
  std.allSites$patient <- sampleInfo[
    match(std.allSites$sampleName, sampleInfo$alias), "patient"]
  
  std.allSites$primerID <- primerIDs[
    paste0(std.allSites$sampleName, "%", std.allSites$ID)]
  
  message("PrimerIDs merged to alignment information.")
  
  possible.contam.primerID <- find_primerID_crossover(std.allSites)
  
  possible.contam <- list(possible.contam, possible.contam.primerID)
  names(possible.contam) <- c("sites_crossing_over", "primerIDs_crossing_over")
}

##### Save data #####
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


