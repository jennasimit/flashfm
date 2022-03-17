#' @title Wrapper to run FINEMAP  (Benner et al. 2016) in R
#' @param GWAS a data.frame with columns: "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"; this is the same as z file of FINEMAP
#' @param ldfile path to file that contains the SNP correlation matrix, with SNPs in the same order as in the GWAS data.frame; this is the same as the ld file of FINEMAP
#' @param N sample size for trait
#' @param fstub file stub for input/output files of FINEMAP, e.g. if fstub="DIRresults/region1", FINEMAP files of the form "DIRresults/region1.z"  will be created
#' @param FMpath file pathway to FINEMAP software e.g. "/software/finemap_v1.4_x86_64/finemap_v1.4_x86_64"
#' @return snpPP a data.frame of top SNP models and their model PPs, as output from FINEMAP
#' @author Jenn Asimit
#' @export
finemap <- function (GWAS, ldfile, N, fstub,FMpath) 
{
    zfile <- paste0(fstub, ".z")
    write.table(GWAS, file = zfile, quote = FALSE, row.names = FALSE, 
        col.names = TRUE)
    message("z file written to ", zfile)
	nsnps <- nrow(GWAS)
    snpfile <- paste0(fstub, ".snp")
    confile <- paste0(fstub, ".config")
    logfile <- paste0(fstub, ".log")
    crfile <- paste0(fstub, ".cred")
    mfile <- paste0(fstub, ".master")
    write.table("z;ld;snp;config;cred;log;n_samples", file = mfile, 
        quote = FALSE, col.names = FALSE, row.names = FALSE, 
        append = FALSE)
    write.table(paste(zfile, ldfile, snpfile, confile, crfile, 
        logfile, N, sep = ";"), file = mfile, quote = FALSE, 
        col.names = FALSE, row.names = FALSE, append = TRUE)
    fline <- paste0(FMpath, " --sss --log --n-configs-top 1000 --n-causal-snps 5 --in-files ", 
        mfile)
    system(fline)
    fmresults <- read.table(confile, header = TRUE, as.is = TRUE, 
        sep = " ")
    PP <- fmresults$prob
    topmods <- fmresults$config
    nmod <- sapply(topmods, length)
    snpPP <- data.frame( snps = topmods, PP = PP, stringsAsFactors = FALSE)  
    return(snpPP)
}



#' @title Wrapper to run single-trait fine-mapping with FINEMAP on each trait, followed by flashfm and then constuct SNP groups for each approach and summarises results
#' @param gwas.list a list containing a data.frame for each trait; each data.frame has columns "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"; this is the same as z file of FINEMAP
#' @param corX genotype correlation matrix (reference or from sample); needs SNP column and row names OR if no names, it takes the names for the raf vector. The SNPs in corX MUST be in the same order as the SNPs in the raf vector. 
#' @param raf named vector of reference allele frequencies; the name of each allele frequency is the SNP ID and MUST be in same SNP order as in corX
#' @param ybar vector of trait means; if traits are transformed to be standard Normal, could set ybar as 0-vector
#' @param N vector of sample sizes for each trait; recommended to give effective sample sizes using GWAS summary statistics in Neff function
#' @param fstub file stub for input/output files of FINEMAP, e.g. if fstub="DIRresults/region1", FINEMAP files of the form "DIRresults/region1.z"  will be created
#' @param TOdds target odds of no sharing to sharing; default is 1
#' @param covY trait covariance matrix (for at most 5 traits and all traits should have a signal in the region, e.g. min p < 1E-6)
#' @param cpp cumulative posterior probability threshold for selecting top models; default cpp=0.99
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @param FMpath file pathway to FINEMAP software e.g. "/software/finemap_v1.4_x86_64/finemap_v1.4_x86_64"
#' @return list with 2 components: mpp.pp, a list with 4 components giving the SNP-level results (mpp.pp$PP,mpp.pp$MPP) and SNP group level results (mpp.pp$MPPg, mpp.pp$PPg); and snpGroups, 
#' a list with 2 components giving the SNP groups construced under single-trait (snpGroups[[1]]) and multi-trait fine-mapping (snpGroups[[2]])
#' @export
#' @author Jenn Asimit
FLASHFMwithFINEMAP <- function(gwas.list, corX, raf, ybar, N, fstub, TOdds = 1, 
    covY, cpp = 0.99, NCORES,FMpath) 
{
    M <- length(ybar)
    Vy <- diag(covY)
    maf <- raf*(raf<.5) + (1-raf)*(raf>=.5)

    if (is.null(colnames(corX))) {
        colnames(corX) <- names(raf)
        rownames(corX) <- names(raf)
        }	
    
    ldfile <- paste0(fstub,".ld")
	write.table(corX,file=ldfile,quote = FALSE, row.names = FALSE, col.names = FALSE)
     message("ld file written to ", ldfile)
    
    fm <- beta <- vector("list", M)
    for (i in 1:M) {
    	fstubT <- paste0(fstub,".",i)
        fm[[i]] <- finemap(gwas.list[[i]], ldfile, N[i], fstubT, FMpath)
        beta[[i]] <- gwas.list[[i]]$beta
        names(beta[[i]]) <- gwas.list[[i]]$rsid
       names(fm) <- names(gwas.list)
    }
    gc(verbose = FALSE)
	covX <-  cor2cov(corX,sd=sqrt(2*raf*(1-raf)))   
	 
	FMmain.input <- flashfm.input(fm,beta1.list=beta,Gmat=covX,Nall=N,ybar.all=ybar,is.snpmat=FALSE,raf=raf)
	ss.stats <- summaryStats(Xmat=FALSE,ybar.all=ybar,main.input=FMmain.input)
	fm.multi <- flashfm(FMmain.input, TOdds=TOdds,covY,ss.stats,cpp=cpp,maxmod=NULL,fastapprox=FALSE)
	snpGroups <- makeSNPgroups2(FMmain.input,fm.multi,is.snpmat=FALSE,r2.minmerge=0.6)
	mpp.pp <- PPsummarise(fm.multi, snpGroups, minPP=0)

    return(list(mpp.pp = mpp.pp, snpGroups = snpGroups))
}

#' @title Construct a credible set for a trait 
#' @param modPP named vector of model PPs
#' @param cred probability for credible set; default is 0.99
#' @return cs vector of SNPs belonging to the credible set
#' @export
#' @author Jenn Asimit
credset <- function(modPP,cred=0.99) {
	tmp <- modPP[order(modPP, decreasing = TRUE)]
    cpp <- cumsum(tmp)
    wh <- which(cpp <= cred)
    if (!length(wh)) wh <- 1
    wh <- c(wh, max(wh) + 1)
    keepmodPP <- tmp[wh]
    mods <- names(keepmodPP)
	cs <- unique(unlist(strsplit(mods,"%")))	
	return(cs)
}
 

#' @title Construct a credible set for each trait and under each of single and multi-trait fine-mapping
#' @param mpp.pp object created in flashfm output, e.g. fm$mpp.pp if fm is output from FLASHFMwithJAM or FLASHFMwithFINEMAP
#' @param cred probability for credible set; default is 0.99
#' @return list with 3 components: list of single-trait fine-mapping credible sets for each trait, list of multi-trait fine-mapping credible sets for each trait, credible set probability
#' @export
#' @author Jenn Asimit
allcredsets <- function(mpp.pp,cred=.99) {
	M <- length(mpp.pp$PP)
	csfm <- csflfm <- vector("list",M)
	for(i in 1:M) {
		csfm[[i]] <- credset(mpp.pp$PP[[i]][,1],cred)
		csflfm[[i]] <- credset(mpp.pp$PP[[i]][,2],cred)
		}
	return(list(fm=csfm,flashfm=csflfm,cred=cred))
}


