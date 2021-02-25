#' Simulated genotype matrix used to generate two traits
#'
#' The matrix is 2000 x 334, where rows are individuals and columns are SNP rsIDs. 
#' For a region of 345 SNPs in chromosme 10p-6030000-6220000 (GRCh37/hg19), containing IL2RA, 
#' we generated a population of 100,000 individuals based on the CEU 1000 Genomes Phase 3 data 
#' using HapGen2. We selected a random sample of 2000 from this population and only retained 
#' the 334 variants with MAF > 0.005 in this sample.
#'
"X"
