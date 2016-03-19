library(GenomicRanges)
library(ShortRead)
library(Biostrings)
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
exdPrimer <- primer.full[15:28]
add_primer <- primer.full[15:21]
extPriLTR <- c(exdPrimer, ltrbit)

get_upstream_seqs <- function(sites, upstream, genome){
  sites.fl <- flank(sites, width = -1, start = TRUE)
  sites.rd <- reduce(sites.fl, min.gapwidth=0)
  sites.up <- flank(sites.rd, width = upstream, start = TRUE)
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
  #score.addprimer <- aln_upstream_seqs(add_primer, DNAStringSet(
  #  seqs, start = unique(width(seqs))-21, 
  #  end = unique(width(seqs))-14),
  #  submat = submat, return = "scores")
  df <- data.frame(
    "posid" = names(seqs),
    "ltrbit.score" = score.ltrbit$score,
    "primer.score" = score.primer$score#,
  #  "addprimer.score" = score.addprimer$score
  )
  df$total.score <- df$ltrbit.score + df$primer.score #+ df$addprimer.score
  df
}
