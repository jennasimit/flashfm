#' Details of SNPs from the simulated genotype matrix X used to generate two traits
#'
#' This is a data.frame with 334 rows, each corresponding to a SNP, and with 
#' 4 columns: ID=rs ID; pos=base-pair position; allele0=non-reference allele; allele1=reference allele
#' For a region of 345 SNPs in chromosme 10p-6030000-6220000 (GRCh37/hg19), containing IL2RA, 
#' we generated a population of 100,000 individuals based on the CEU 1000 Genomes Phase 3 data 
#' using HapGen2. We selected a random sample of 2000 from this population and only retained 
#' the 334 variants with MAF > 0.005 in this sample. 
#'
"snpinfo"
