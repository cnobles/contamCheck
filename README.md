# contamCheck
Check data returned from INSPIIRED::intSiteCaller for overlapping integration site locations in different patients. If returned, it may be possible that the overlapping integration sites are present due to a contamination and further investigation is required.

Use:
```
#Needs to be run in primary analysis directory
Rscript path/to/contamCheck/check_contam.R
Rscript path/to/contamCheck/check_contam.R -s hiv_specimen.database -i hiv_intsites.database
Rscript path/to/contamCheck/check_contam.R --specimen_database hiv_specimen.database --intsites_database hiv_intsites.database
```

The script produces a tsv file ({DIRECTORY_NAME}.contam.stats.tsv) and an RData file ({DIRECTORY_NAME}.possible.contam.RData) which contains the positions of integrations sites and which samples they were found in.

Limitations: This script can currently only be applied to sequencing runs with multiple patients multiplexed together. Also this script only analyzes integration site positions. Identical positions could be found if there are very high numbers of integration sites in each sample, the birthday problem.
