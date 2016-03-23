library(GenomicRanges)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)
library(devtools)
source_url("https://raw.githubusercontent.com/cnobles/cloneTracker/master/cloneTracker.SOURCE_ME.R")
genome <- BSgenome.Hsapiens.UCSC.hg18

submat <- nucleotideSubstitutionMatrix(
  match = 3, mismatch = 1, baseOnly = TRUE)

primer.full <- DNAString("CAGACCCTTTTAGTCAGTGTGGAAAATC")
primer <- DNAString("GAAAATC")
ltrbit <- DNAString("TCTAGCA")
primerltr <- c(primer, ltrbit)
extPrimer <- primer.full[15:28]
add_primer <- primer.full[15:21]
extPriLTR <- c(extPrimer, ltrbit)

get_upstream_seqs <- function(sites, upstream, genome){
  sites.fl <- flank(sites, width = -1, start = TRUE)
  sites.rd <- reduce(sites.fl, min.gapwidth=0)
  sites.up <- flank(sites.rd, width = upstream, start = TRUE)
  start(sites.up) <- ifelse(strand(sites.up) == "+", start(sites.up), start(sites.up)-1)
  end(sites.up) <- ifelse(strand(sites.up) == "+", end(sites.up)+1, end(sites.up))
  seqs <- getSeq(genome, sites.up)
  names(seqs) <- generate_posID(sites.rd)
  seqs
}

aln_upstream_seqs <- function(pattern, seqs, submat, 
                              gapOpen = 2, gapExt = 1, return = "scores"){
  aln <- pairwiseAlignment(
    pattern = seqs,
    subject = pattern,
    substitutionMatrix = submat,
    gapOpening = gapOpen,
    gapExtension = gapExt,
    type = "overlap")
  
  score.df <- data.frame(
    "posid" = names(seqs),
    "score" = score(aln),
    "aln_length" = width(aln@subject@range),
    "mismatches"= sapply(aln@subject@mismatch, length),
    "indels" = sapply(aln@subject@indel, length),
    "seq" = as.character(seqs)
  )
  if(return == "scores"){
    score.df
  }else if(return == "aln"){
    aln
  }else{
    message("Choose to return either scores or alignments (aln).")
  }
}

test_confidence_scores <- function(seqs, ltrbit, primer){
  score.ltrbit <- aln_upstream_seqs(ltrbit, DNAStringSet(
    seqs, start = unique(width(seqs))-length(ltrbit), 
    end = unique(width(seqs))),
    submat = submat, return = "scores")
  score.primer <- aln_upstream_seqs(primer, DNAStringSet(
    seqs, start = unique(width(seqs))-(length(ltrbit)+length(primer)), 
    end = unique(width(seqs))-length(ltrbit)),
    submat = submat, return = "scores")
  df <- data.frame(
    "posid" = names(seqs),
    "ltrbit.score" = score.ltrbit$score,
    "primer.score" = score.primer$score#,
  )
  df$total.score <- df$ltrbit.score + df$primer.score 
  df
}

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

load_and_process_caller_data <- function(dataDir){
  siteInfo <- read.csv(paste0(dataDir, "/sampleInfo.tsv"), sep = "\t")
  sites <- lapply(siteInfo[, "alias"], function(sampleName){
    load_intSiteCaller_data(sampleName, "allSites", dataDir)
  })
  sites <- unlist(GRangesList(sites[sapply(sites, class) == "GRanges"]))
  sites <- remove_repeats(sites)
  sites <- standardize_intsites(sites, standardize_breakpoints = FALSE)
  seqs <- get_upstream_seqs(sites, upstream = 35, genome)
  seqs <- DNAStringSet(grep("N", seqs, value = TRUE, invert = TRUE))
  test_confidence_scores(seqs, ltrbit, extPrimer)
}

load_and_process_possible_contam <- function(dataDir){
  list <- list.files(dataDir)
  contam.file <- grep("contam.RData", list, value = TRUE)
  if(length(contam.file) == 0){stop("No possible.contam file found.")}
  load(paste0(dataDir, "/", contam.file))
  contamSites <- unlist(possible.contam)
  contamSites <- remove_repeats(contamSites)
  contamSites <- standardize_intsites(contamSites, standardize_breakpoints = TRUE)
  contamSeqs <- get_upstream_seqs(contamSites, upstream = 35, genome)
  contamSeqs <- DNAStringSet(grep("N", contamSeqs, value = TRUE, invert = TRUE))
  test_confidence_scores(contamSeqs, ltrbit, extPrimer)
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
  pairwise.t.test(df$score, df$set, p.adjust.method = "bonferroni")
}
