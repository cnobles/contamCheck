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
hapSites <- remove_repeats(hapSites)
hapSeqs <- get_upstream_seqs(hapSites, upstream = 35, genome)
hapSeqs <- DNAStringSet(grep("N", hapSeqs, value = TRUE, invert = TRUE))
hapScores <- test_confidence_scores(hapSeqs, ltrbit, extPrimer)

#hiv set 0 - 0 mismatch
hiv0_0_Scores <- load_and_process_caller_data("~/intSiteData/zero_mismatch/HIVLAT_Set_0_std")
#hiv set 0 - 1 mismatch
hiv0_1_Scores <- load_and_process_caller_data("~/intSiteData/one_mismatch/Sample_Set_0")
#hiv set 0 - 2 mismatch
hiv0_2_Scores <- load_and_process_caller_data("~/intSiteData/multihit_clus_branch/Sample_Set_0")

#hiv set 0 - 1 mismatch crossover
hiv0_CO_1_Scores <- load_and_process_possible_contam("~/intSiteData/one_mismatch/Sample_Set_0")
#hiv set 0 - 2 mismatch crossover
hiv0_CO_2_Scores <- load_and_process_possible_contam("~/intSiteData/multihit_clus_branch/Sample_Set_0")


#hiv set 1 - 0 mismatch
hiv1_0_Scores <- load_and_process_caller_data("~/intSiteData/zero_mismatch/HIVLAT_Set_1_std")
#hiv set 1 - 1 mismatch
hiv1_1_Scores <- load_and_process_caller_data("~/intSiteData/one_mismatch/Set_1")
#hiv set 1 - 2 mismatch
hiv1_2_Scores <- load_and_process_caller_data("~/intSiteData/multihit_clus_branch/HIVLAT_Set_1")

#hiv set 1 - 0 mismatch crossover
hiv1_CO_0_Scores <- load_and_process_possible_contam("~/intSiteData/zero_mismatch/HIVLAT_Set_1_std")
#hiv set 1 - 1 mismatch crossover
hiv1_CO_1_Scores <- load_and_process_possible_contam("~/intSiteData/one_mismatch/Set_1")
#hiv set 1 - 2 mismatch crossover
hiv1_CO_2_Scores <- load_and_process_possible_contam("~/intSiteData/multihit_clus_branch/HIVLAT_Set_1")

#Partial sets looking at the upstream sequence complementartity to the primer sequence
cont.sets <- list(simScores, hapScores)
names(cont.sets) <- c("simSites", "Lib4-HAP1")
cont.plot <- plot_scores(cont.sets, score.type = "primer.score")
cont.plot

hiv0.sets <- list(simScores, hiv0_0_Scores, hiv0_1_Scores, hiv0_2_Scores)
names(hiv0.sets) <- c("simSites", "Set_0_0mismatch", "Set_0_1mismatch", "Set_0_2mismatch")
hiv0.plot <- plot_scores(hiv0.sets, score.type = "primer.score")
hiv0.plot

hiv1.sets <- list(simScores, hiv1_0_Scores, hiv1_1_Scores, hiv1_2_Scores)
names(hiv1.sets) <- c("simSites", "Set_1_0mismatch", "Set_1_1mismatch", "Set_1_2mismatch")
hiv1.plot <- plot_scores(hiv1.sets, score.type = "primer.score")
hiv1.plot

hiv0.CO.sets <- list(simScores, hiv0_CO_1_Scores, hiv0_CO_2_Scores)
names(hiv0.CO.sets) <- c("simSites", "Set_0_CO_1mismatch", "Set_0_CO_2mismatch")
hiv0.CO.plot <- plot_scores(hiv0.CO.sets, score.type = "primer.score")
hiv0.CO.plot

hiv1.CO.sets <- list(simScores, hiv1_CO_0_Scores, hiv1_CO_1_Scores, hiv1_CO_2_Scores)
names(hiv1.CO.sets) <- c("simSites", "Set_1_CO_0mismatch", "Set_1_CO_1mismatch", "Set_1_CO_2mismatch")
hiv1.CO.plot <- plot_scores(hiv1.CO.sets, score.type = "primer.score")
hiv1.CO.plot

hiv.sets <- list(simScores, hapScores, hiv0_0_Scores, hiv1_0_Scores)
names(hiv.sets) <- c("simSites", "Lib4-HAP1", "Set_0_0mismatch", "Set_1_0mismatch")
hiv.plot <- plot_scores(hiv.sets, score.type = "primer.score")
hiv.stats <- t.test_scores(hiv.sets, score.type = "primer.score")
hiv.plot
hiv.stats

#Partial sets looking at the upstream sequence complementarity to the LTRbit.
ltr.cont.plot <- plot_scores(cont.sets, score.type = "ltrbit.score")
ltr.cont.plot
ltr.hiv1.CO.plot <- plot_scores(hiv1.CO.sets, score.type = "ltrbit.score")
ltr.hiv1.CO.plot
ltr.hiv.plot <- plot_scores(hiv.sets, score.type = "ltrbit.score")
ltr.hiv.plot



#Building data set for trainging profile hidden Markov models
ks.test(simScores$primer.score, round(rnorm(5000, mean = simMean, sd = sqrt(simVar))))

hiv0_2_sampleInfo <- read.csv("~/intSiteData/multihit_clus_branch/Sample_Set_0/sampleInfo.tsv",
                              sep = "\t")
hiv0_2_sites <- lapply(hiv0_2_sampleInfo$alias, function(sampleName){
  load_intSiteCaller_data(sampleName = sampleName, 
                          dataType = "allSites", 
                          dataDir = "~/intSiteData/multihit_clus_branch/Sample_Set_0")
})
hiv0_2_sites <- unlist(GRangesList(hiv0_2_sites[
  sapply(hiv0_2_sites, class) == "GRanges"]))
hiv0_2_seqs <- get_upstream_seqs(hiv0_2_sites, upstream = 50, genome = genome)

misprimed.sites <- as.character(hiv0_2_Scores[
  !hiv0_2_Scores$posid %in% hiv0_0_Scores$posid, "posid"
  ])
misprimedScores <- hiv0_2_Scores[hiv0_2_Scores$posid %in% misprimed.sites,]
misprimedScores <- misprimedScores[misprimedScores$primer.score > 26,]
misprimed.sites <- misprimedScores$posid


misprimed.seqs <- hiv0_2_seqs[match(misprimed.sites, names(hiv0_2_seqs))]
subset.misprimed.seqs <- sample(misprimed.seqs, 250)
set.misprimed.seqs <- misprimed.seqs[!names(misprimed.seqs) %in% names(subset.misprimed.seqs)]
writeXStringSet(subset.misprimed.seqs, 
                filepath = "~/intSiteData/subset.misprimed.hiv0.fasta", 
                format = "fasta")
writeXStringSet(set.misprimed.seqs, 
                filepath = "~/intSiteData/set.misprimed.hiv0.fasta", 
                format = "fasta")

hiv1_2_sampleInfo <- read.csv("~/intSiteData/multihit_clus_branch/HIVLAT_Set_1/sampleInfo.tsv",
                              sep = "\t")
hiv1_2_sites <- lapply(hiv1_2_sampleInfo$alias, function(sampleName){
  load_intSiteCaller_data(sampleName = sampleName, 
                          dataType = "allSites", 
                          dataDir = "~/intSiteData/multihit_clus_branch/HIVLAT_Set_1")
})
hiv1_2_sites <- unlist(GRangesList(hiv1_2_sites[
  sapply(hiv1_2_sites, class) == "GRanges"]))
hiv1_2_seqs <- get_upstream_seqs(hiv1_2_sites, upstream = 50, genome = genome)
writeXStringSet(hiv1_2_seqs, filepath = "~/intSiteData/set.hiv1.2mismatch.fasta", format = "fasta")

