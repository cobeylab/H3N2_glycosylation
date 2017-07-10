### 160 allele frequencies

The code to plot recent changes in 160 allele frequencies is in [link]().
It takes the GISAID frequencies ([link]()) as raw input.

### Serological analyses

All of the statistical analysis and figures involving serology are in the script [`Serological_analyses.R`](Serology/Serological_analyses.R). 
These analyses use as input
* the FRNT data
* subjects' probabilities of imprinting on 160K (see below)

Some of the multivariate regression models test the effects of primary exposure to the 160K allele. 
The [`160K_imprinting.R`](Serology/160K_imprinting.R) script calculates the probability of primary exposure to this allele, and it also generates figures showing 160K allele frequencies from 1968 through 2010 and the imprinting probabilities for the birth years relevant to the study. 
The code takes as input 
* the frequencies of different 160 alleles over this period ([`160_freq.csv`](Serology/160_freq.csv)).
* the frequencies of H3N2 v. H1N1 and B over the same period ([`H3_subtype_frac.csv`](Serology/H3_subtype_frac.csv)).

The `US_only` column defines if data from strictly the U.S. was used, or if frequencies were calculated based on all global samples from that year.
