#' List of GWAS data.frames for each of the two simulated traits; format of z files of FINEMAP
#'
#' This is a list of two data.frames - one for each trait. 
#' It has columns "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se".
#' We generated a population of 100,000 individuals based on the CEU 1000 Genomes Phase 3 data 
#' using HapGen2. We selected a random sample of 2000 from this population and only retained 
#' the 334 variants with MAF > 0.005 in this sample. 
#'
"gwas.list"
