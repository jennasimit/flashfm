#' Output from FLASHFMwithFINEMAP for two simulated traits 
#' 
#' Simulated data was generated using the mvrnorm function in R and genotype matrix X. The two traits have correlation 0.4 and 
#' each has two causal variants, of which one, rs61839660 is shared between traits; trait 1 has second 
#' causal variant rs62626317 and trait 2 has second causal variant rs11594656. The traits are from a sample
#' of 2000, with no missing measurements.
#' 
#' fm2 was generated using the following command, where all called objects exist in the flashfm package, 
#' fstub is a file prefix for FINEMAP input/output files, and FMpath is the path to the software FINEMAP
#' fm2 <- FLASHFMwithFINEMAP(gwas.list,corX,raf,ybar,N,fstub="DIRresults/eg",TOdds=1,covY,cpp=.99,NCORES=2,FMpath="/software/finemap_v1.4_x86_64/finemap_v1.4_x86_64") 
#' JAMexpanded.multi(beta,X,snpinfo,ybar,diag(covY),N,chr=10,fstub=fstub,mafthr=0.005,path2plink="/software/plink",r2=0.99,save.path=fstubJ,related=FALSE,y=NULL)
#'
"fm2"
