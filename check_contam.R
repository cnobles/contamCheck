options(stringsAsFactors = FALSE)

#Check and load packages
rPackages <- c("RMySQL", "devtools")
stopifnot(all(sapply(rPackages, require, 
                     character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))
source_url("https://raw.githubusercontent.com/cnobles/cloneTracker/master/cloneTracker.SOURCE_ME.R")

##Set up working directory and load sampleInfo
dataDir <- getwd()
splitDir <- unlist(strsplit(dataDir, split = "/"))
runName <- splitDir[length(splitDir)]
sampleInfoFile <- grep("sampleInfo", list.files(path = dataDir), value = TRUE)
sampleInfo <- read.csv(paste(dataDir,sampleInfoFile, sep = "/"), sep = "\t")
sampleInfo$GTSP <- sapply(strsplit(
  as.character(sampleInfo$alias), 
  "-"), "[", 1)

if(any(grepl("UninfectedControl", sampleInfo$GTSP))){
  uncPos <- grep("UninfectedControl", sampleInfo$GTSP)
  sampleInfo[uncPos, "GTSP"] <- "UNC"
  rm(uncPos)
}

if(any(grepl("NoTemplateControl", sampleInfo$GTSP))){
  ntcPos <- grep("NoTemplateControl", sampleInfo$GTSP)
  sampleInfo[ntcPos, "GTSP"] <- "NTC"
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
GTSP_query <- grep("GTSP", unique(sampleInfo$GTSP), value = TRUE)

message("Querying patient information for:")
message(list(GTSP_query))
               
junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
dbConn <- dbConnect(MySQL(), group="specimen_management")  
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1) 

string <- paste(GTSP_query, collapse="' OR SpecimenAccNum = '")
query_request <- paste0("WHERE SpecimenAccNum = '", string, "'")
sql <- paste0("SELECT SpecimenAccNum,Patient FROM specimen_management.gtsp ",
              query_request) 
patientInfo <- suppressWarnings( dbGetQuery(dbConn, sql) ) 
dbDisconnect(dbConn) #Disconnect from specimen_management after query
rm(junk, dbConn, string, query_request)

isThere <- unique(patientInfo$SpecimenAccNum)
stopifnot( (NA %in% match(GTSP_query, isThere)) == FALSE)
message("Received patient information for:")
message(list(isThere))

if(any(grepl("UNC", sampleInfo$GTSP))){
  patientInfo <- rbind(patientInfo, data.frame(
    "SpecimenAccNum" = "UNC", "Patient" = "UNC")
    )
}
if(any(grepl("NTC", sampleInfo$GTSP))){
  patientInfo <- rbind(patientInfo, data.frame(
    "SpecimenAccNum" = "NTC", "Patient" = "NTC")
    )
}

sampleInfo$patient <- patientInfo[match(sampleInfo$GTSP, 
                                        patientInfo$SpecimenAccNum), "Patient"]

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
uniq.allSites$GTSP <- sampleInfo[
  match(uniq.allSites$sampleName, sampleInfo$alias), "GTSP"]
uniq.allSites$patient <- sampleInfo[
  match(uniq.allSites$sampleName, sampleInfo$alias), "patient"]
uniq.allSites <- split(uniq.allSites, uniq.allSites$patient)
message("Checking for contamination between:")
message(list(names(uniq.allSites)))
possible.contam <- track_clones(uniq.allSites, gap = 5L, track.origin = FALSE)

#Post crossover check removal of redundant reads
if(length(possible.contam) > 0){
  possible.contam <- unlist(possible.contam)
  possible.contam <- split(possible.contam, possible.contam$sampleName)
  possible.contam <- unlist(unique(possible.contam))
  possible.contam <- split(possible.contam, possible.contam$posid)
}
possible.contam.gr <- unlist(possible.contam)

#Save data
message("Saving data.")
save(possible.contam, 
     file = paste0(dataDir, "/", runName, ".possible.contam.RData"))
stats <- data.frame(
  "amount" = c(length(unique(
      sampleInfo[grep("GTSP", sampleInfo$GTSP), "patient"]
    )),
    length(unique(
      sampleInfo[grep("GTSP", sampleInfo$GTSP, invert = TRUE), "patient"]
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
write.table(stats, file = paste0(dataDir, "/", runName, ".contam.stats.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
message("Check complete.")


