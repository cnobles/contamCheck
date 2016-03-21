library(ggplot2)

#Test data sets

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

density_plot_score_output <- function(score.df){
  scoreTypes <- names(score.df[,2:ncol(score.df)])
  df <- data.frame(
    "sequence" = c(sapply(1:length(scoreTypes), function(i) 
      rep(scoreTypes[i], nrow(score.df)))),
    "score" = c(sapply(1:length(scoreTypes), function(i)
      score.df[,i+1]))
  )
  ggplot(df, aes(score, fill = sequence, colour = sequence)) + geom_density(alpha = 0.1)
}

combine_score_sets <- function(score.list, score.type = "total.score"){
  scoreSets <- names(score.list)
  data.frame(
    "set" = unlist(sapply(1:length(score.list), function(i)
      rep(scoreSets[i], nrow(score.list[[i]])))),
    "score" = unlist(sapply(1:length(scoreSets), function(i)
      score.list[[i]][, score.type]))
  )
}

plot_scores <- function(score.list, score.type = "total.score"){
  df <- combine_score_sets(
    score.list = score.list,
    score.type = score.type)
  ggplot(df, aes(score, fill = set, colour = set)) + 
    geom_density(alpha = 0.1, adjust = 2/3)
}

t.test_scores <- function(score.list, score.type = "total.score"){
  df <- combine_score_sets(
    score.list = score.list,
    score.type = score.type)
  pairwise.t.test(df$score, df$set, p.adjust.method = "holm")
}

#simulated Sites
simSitesInfo <- read.csv(
  "../intSiteData/multihit_clus_branch/simSites/sampleInfo.tsv", 
  sep = "\t")
simSites <- unlist(GRangesList(lapply(simSitesInfo$alias, function(sampleName){
  load_intSiteCaller_data(
    sampleName,
    dataType = "allSites",
    dataDir = "../intSiteData/multihit_clus_branch/simSites")
})))
simSeqs <- get_upstream_seqs(simSites, upstream = 35, genome)
simSeqs <- DNAStringSet(grep("N", simSeqs, value = TRUE, invert = TRUE))
simScores <- test_confidence_scores(simSeqs, ltrbit, extPrimer)

#HAP-1 sites 
hapSitesInfo <- read.csv("../intSiteData/older_branches/blocking_oligo/sampleInfo.tsv", sep = "\t")
hapSites <- unlist(GRangesList(lapply(hapSitesInfo[1:12, "alias"], function(sampleName){
  load_intSiteCaller_data(
    sampleName,
    dataType = "allSites",
    dataDir = "../intSiteData/older_branches/blocking_oligo"
  )
})))
hapSeqs <- get_upstream_seqs(hapSites, upstream = 35, genome)
hapSeqs <- DNAStringSet(grep("N", hapSeqs, value = TRUE, invert = TRUE))
hapScores <- test_confidence_scores(hapSeqs, ltrbit, extPrimer)

#hiv set 0 - 0 mismatch
hiv0SiteInfo <- read.csv("../intSiteData/zero_mismatch/HIVLAT_Set_0_std/sampleInfo.tsv", sep = "\t")
hiv0Sites <- lapply(hiv0SiteInfo[, "alias"], function(sampleName){
  load_intSiteCaller_data(
    sampleName,
    dataType = "allSites",
    dataDir = "../intSiteData/zero_mismatch/HIVLAT_Set_0_std"
  )
})
hiv0Sites <- unlist(GRangesList(hiv0Sites[sapply(hiv0Sites, class) == "GRanges"]))
hiv0Seqs <- get_upstream_seqs(hiv0Sites, upstream = 35, genome)
hiv0Seqs <- DNAStringSet(grep("N", hiv0Seqs, value = TRUE, invert = TRUE))
hiv0Scores <- test_confidence_scores(hiv0Seqs, ltrbit, extPrimer)

#hiv set 1 - 0 mismatch
hiv1SiteInfo <- read.csv("../intSiteData/zero_mismatch/HIVLAT_Set_1_std/sampleInfo.tsv", sep = "\t")
hiv1Sites <- lapply(hiv1SiteInfo[, "alias"], function(sampleName){
  load_intSiteCaller_data(
    sampleName,
    dataType = "allSites",
    dataDir = "../intSiteData/zero_mismatch/HIVLAT_Set_1_std"
  )
})
hiv1Sites <- unlist(GRangesList(hiv1Sites[sapply(hiv1Sites, class) == "GRanges"]))
hiv1Seqs <- get_upstream_seqs(hiv1Sites, upstream = 35, genome)
hiv1Seqs <- DNAStringSet(grep("N", hiv1Seqs, value = TRUE, invert = TRUE))
hiv1Scores <- test_confidence_scores(hiv1Seqs, ltrbit, extPrimer)

#hiv set 1 contamSites
load("../intSiteData/zero_mismatch/HIVLAT_Set_1_std/HIVLAT_Set_1_std.possilble.contam.RData")
hiv1ContamSites <- unlist(possible.contam)
load("../intSiteData/zero_mismatch/HIVLAT_Set_1_alt3LTR/HIVLAT_Set_1_alt3LTR.possilble.contam.RData")
hiv1ContamSites <- c(hiv1ContamSites, unlist(possible.contam))
hiv1ContamSeqs <- get_upstream_seqs(hiv1ContamSites, upstream = 35, genome)
hiv1ContamSeqs <- DNAStringSet(grep("N", hiv1ContamSeqs, value = TRUE, invert = TRUE))
hiv1ContamScores <- test_confidence_scores(hiv1ContamSeqs, ltrbit, extPrimer)

#Partial sets
sets.list <- list(hapScores, hiv1Scores, hiv1ContamScores)
names(sets.list) <- c("hap", "hiv_set_1", "hiv_contam")

sets.list <- list(simScores, hapScores, hiv0Scores, hiv1Scores, hiv1ContamScores)
names(sets.list) <- c("sim", "hap", "hiv_set_0", "hiv_set_1", "hiv_contam")

plot_scores(sets.list, score.type = "total.score")
plot_scores(sets.list, score.type = "ltrbit.score")
plot_scores(sets.list, score.type = "primer.score")
