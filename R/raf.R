#' Vector of RAFs from the simulated genotype matrix X used to generate two traits
#'
#' This is a named vector of length 334, where names are SNP rsIDs. 
#' For a region of 345 SNPs in chromosme 10p-6030000-6220000 (GRCh37/hg19), containing IL2RA, 
#' we generated a population of 100,000 individuals based on the CEU 1000 Genomes Phase 3 data 
#' using HapGen2. We selected a random sample of 2000 from this population and only retained 
#' the 334 variants with MAF > 0.005 in this sample. 
#'
"raf"
