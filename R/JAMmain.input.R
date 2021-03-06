#' Output from JAMexpanded.multi for two simulated traits 
#' 
#' Simulated data was generated using the mvrnorm function in R and genotype matrix X. The two traits have correlation 0.4 and 
#' each has two causal variants, of which one, rs61839660 is shared between traits; trait 1 has second 
#' causal variant rs62626317 and trait 2 has second causal variant rs11594656. The traits are from a sample
#' of 2000, with no missing measurements.
#' 
#' JAMmain.input was generated using the following command, where all called objects exist in the flashfm package, 
#' fstub and fstubJ are paths with file prefixes for output data, and path2plink is the path to the software plink 
#' JAMexpanded.multi(beta,X,snpinfo,ybar,diag(covY),N,chr=10,fstub=fstub,mafthr=0.005,path2plink="/software/plink",r2=0.99,save.path=fstubJ,related=FALSE,y=NULL)
#'
"JAMmain.input"
