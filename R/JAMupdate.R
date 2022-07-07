cor.refdata2 <- function (corX,r2 = 0.99) 
{
    gmat2t <- tagSNP(corX, threshold = sqrt(r2))
    gt <- unlist(gmat2t)
    tg = gt[names(gt) == "tagsnp"]
    refG <- corX[tg, tg]
    return(list(refG=refG,taglist=gmat2t))
}

make.nonpos.def <- function(refG) {
# copied from https://www.r-bloggers.com/2012/10/fixing-non-positive-definite-correlation-matrices-using-r-2/
# using the method of Rebonato and Jackel (2000), as elaborated by Brissette et al. (2007), to fix the correlation matrix. As per the method, replace the negative eigenvalues with 0 (or a small positive number as Brissette et al. 2007 suggest), then normalize the new vector.
# The paper by Rebonato and Jackel, “The most general methodology for creating a valid correlation matrix for risk management and option pricing purposes”, Journal of Risk, Vol 2, No 2, 2000, presents a methodology to create a positive definite matrix out of a non-positive definite matrix. 
origMat <- refG
origEig <- eigen(origMat)

cholStatus <- try(u <- chol(origMat), silent = TRUE)
cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)

# fix the correl matrix

newMat <- origMat

iter <- 0
while (cholError) {

    iter <- iter + 1
#    cat("iteration ", iter, "\n")

    # replace -ve eigen values with small +ve number
    newEig <- eigen(newMat)
    newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)

    # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
    # eig vectors
    newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)

    # normalize modified matrix eqn 6 from Brissette et al 2007
    newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))

    # try chol again
    cholStatus <- try(u <- chol(newMat), silent = TRUE)
    cholError <- ifelse(class(cholStatus)[1] == "try-error", TRUE, FALSE)
}

colnames(newMat) <- rownames(newMat) <- colnames(refG)
refG <- newMat
return(refG)
}


#' @title Expanded version of JAM (a single-trait fine-mapping approach) that first runs on thinned SNPs and then expands models on tag SNPs; this can run independently on multiple traits
#' This version is more stable than JAMexpandedCor.multi, but slower, so is run only if JAMexpandedCor.multi fails
#' @param beta1 list where each component is a named vector of of single SNP effect estimates for a trait; one vector for each trait
#' @param corX genotype correlation matrix (reference or from sample) 
#' @param raf named vector of reference allele frequencies; the name of each allele frequency is the SNP ID and MUST be in same SNP order as in corX
#' @param ybar vector of trait means; if related samples, this should be based on unrelated samples; if traits are transformed to be standard Normal, could set ybar as 0-vector
#' @param Vy vector of trait variances; if related samples, this should be based on unrelated samples; if traits are transformed to be standard Normal, could set Vy as 1-vector
#' @param N vector of sample sizes for each trait; recommended to give effective sample sizes using GWAS summary statistics in Neff function
#' @param r2 r.squared threshold for thinning SNPs before JAM and finding tag SNPs
#' @param save.path path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1"). 
#' @return list with 3 components: SM a list of snpmod objects giving fine-mapping results for each trait; mbeta a list of joint effect estimates for each trait; nsnps number of SNPs
#' @import R2BGLiMS
#' @export
#' @author Jenn Asimit
JAMexpandedCor.multi2 <- function (beta1, corX, raf, ybar, Vy, N, r2 = 0.99, save.path) {
	covX <- cor2cov(corX, sd = sqrt(2 * raf * (1 - raf)))
    Nlist <- makeNlist.rel(Ne = N)
    N <- diag(Nlist$Nqq)
    M <- length(beta1)
    if (is.null(names(beta1))) {
        ts <- paste0("T", 1:M)
        names(ybar) <- ts
    } else {
        ts <- names(ybar)
    }
    snps <- names(beta1[[1]])
    for (i in 2:M) snps <- intersect(snps, names(beta1[[i]]))
    snps <- intersect(snps, names(raf))
    for (i in 1:M) beta1[[i]] <- beta1[[i]][snps]
    if (is.null(colnames(corX))) {
        colnames(corX) <- names(raf)
        rownames(corX) <- names(raf)
    }
    corX <- corX[snps, snps]
    covX <- covX[snps, snps]
    maf <- raf * (raf <= 0.5) + (1 - raf) * (raf > 0.5)
    nsnps <- ncol(corX)
#    refG <- cor.refdata.fn(corX, r2)
     reftags <- cor.refdata2(corX, r2)
    refGt <- reftags$refG
    taglist <- reftags$taglist
    
    refG <- make.nonpos.def(refGt)
    
#    corX2 <- corX^2
#    taglist <- tagSNP(corX2, threshold = r2)
    out <- list(SM = NULL, mbeta = NULL, Nlist = Nlist)
    out$SM <- vector("list", M)
    names(out$SM) <- ts
    out$mbeta <- vector("list", M)
    names(out$mbeta) <- ts
    dd <- vector("list", M)
    SSy <- vector("list", M)
    Sxy <- vector("list", M)
    for (j in 1:M) {
        BETA <- beta1[[j]][colnames(refG)]
        mafs.ref <- maf[colnames(refG)]
        jam.results <- R2BGLiMS::JAM(marginal.betas = BETA, trait.variance=Vy[j],
            cor.ref = refG, mafs.ref = mafs.ref, model.space.priors = list(a = 1, 
                b = length(BETA), Variables = names(BETA)), n = N[j], max.model.dim =10,
            xtx.ridge.term = 0.001, save.path = save.path)
             
        topmods <- R2BGLiMS::TopModels(jam.results, n.top.models = 1000)
 if(is.null(ncol(topmods))) {
         stop("A single model of 10 variants was selected with PP=1. This may mean no convergence because a causal variant 
         is missing from the data and it has no tags in your data.")
         }
        binout <- as.matrix(topmods[, -ncol(topmods)])
        colnames(binout) <- colnames(topmods)[-ncol(topmods)]
        snpmods <- apply(binout, 1, mod.fn)
        nmod <- apply(binout, 1, sum)
        PP <- topmods[, ncol(topmods)]
        snpPP <- data.frame(rank = 1:length(nmod), size = nmod, 
            logPP = log(PP), PP = PP, str = snpmods, snps = snpmods, 
            stringsAsFactors = FALSE)
        snpPP <- snpPP[order(snpPP$PP, decreasing = TRUE), ]
        expmods <- rlist::list.stack(lapply(snpPP$snps, tagexpand.mod, 
            taglist = taglist))
        wh <- which(duplicated(expmods$snps))
        if (length(wh) > 0) {
            expmods <- expmods[-wh, ]
        }
        row.names(expmods) <- expmods$snps
        check <- sapply(strsplit(expmods[, 2], "%"), function(x) length(x) > 
            length(unique(x)))
        if (sum(check) > 0) 
            expmods <- expmods[-which(check), ]
        mbeta <- lapply(expmods[, 2], multibeta, beta1[[j]], 
            covX, N = N[j], ybar = ybar[j], is.snpmat = FALSE, 
            raf = raf)
        names(mbeta) <- expmods[, 2]
        SSy[[j]] <- Vy[j] * (N[j] - 1) + N[j] * ybar[j]^2
        Vx <- diag(covX)
        Mx <- 2 * raf
        Sxy[[j]] <- c(Sxy.hat(beta1 = beta1[[j]], Mx = Mx, N = N[j], 
            Vx = Vx, muY = ybar[j]), `1` = ybar[j] * N[j])
        names(Sxy[[j]])[length(Sxy[[j]])] <- "one"
        lABF <- sapply(expmods$snps, calcABF, mbeta, SSy = SSy[[j]], 
            Sxy = Sxy[[j]], Vy = Vy[j], N = N[j])
        names(lABF) <- expmods$snps
        wh <- which(expmods$snps == "1")
        if (!length(wh)) {
            dd[[j]] <- data.frame(model = c("1", expmods$snps), 
                tag = c(FALSE, expmods$tag), lBF = c(0, lABF), 
                stringsAsFactors = FALSE)
            l1 <- multibeta("1", beta1[[j]], covX, N = N[j], 
                ybar = ybar[j], is.snpmat = FALSE, raf = raf)
            mbeta <- rlist::list.append(mbeta, `1` = l1)
        } else {
            dd[[j]] <- data.frame(model = expmods$snps, tag = expmods$tag, 
                lBF = lABF, stringsAsFactors = FALSE)
        }
        EXcv <- round(median(snpPP$size)) 
        SM <- makesnpmod(dd[[j]], expected = EXcv, nsnps = nsnps)
        out$SM[[j]] <- SM
        out$mbeta[[j]] <- mbeta
        names(out$SM) <- ts
        names(out$mbeta) <- ts
        out$SM[[j]] <- best.models.cpp(out$SM[[j]], maxmod = 1000)[[1]]
        out$SM[[j]] <- PP2snpmod(out$SM[[j]])
    }
    out$Nlist <- Nlist
    out$nsnps <- nsnps
    out$Gmat = covX
    out$beta1.list = beta1
    out$raf = raf
    return(out)
}



JAMcormulti.tries <- function(beta1, corX, raf, ybar,Vy, N, save.path) {
    tryCatch({
      JAMexpandedCor.multi(beta1, corX, raf, ybar, Vy,N, r2=.99, save.path)
    },
    error=function(e) { 
      JAMexpandedCor.multi2(beta1, corX, raf, ybar, Vy,N, r2=.99, save.path)
    }
  )
}


#' @title Wrapper to run single-trait fine-mapping with JAMexpandedCor.multi on each trait, followed by flashfm and then constuct SNP groups for each approach and summarises results
#' @param beta1 list where each component is a named vector of of single SNP effect estimates for a trait; one vector for each trait
#' @param corX genotype correlation matrix (reference or from sample) 
#' @param raf named vector of reference allele frequencies; the name of each allele frequency is the SNP ID and MUST be in same SNP order as in corX
#' @param ybar vector of trait means; if related samples, this should be based on unrelated samples; if traits are transformed to be standard Normal, could set ybar as 0-vector
#' @param N vector of sample sizes for each trait; recommended to give effective sample sizes using GWAS summary statistics in Neff function
#' @param save.path path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1"). 
#' @param TOdds Vector of target odds of no sharing to sharing
#' @param covY trait covariance matrix (for at most 5 traits and all traits should have a signal in the region, e.g. min p < 1E-6)
#' @param cpp cumulative posterior probability threshold for selecting top models
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @return list with 2 components: mpp.pp, a list with 4 components giving the SNP-level results (mpp.pp$PP,mpp.pp$MPP) and SNP group level results (mpp.pp$MPPg, mpp.pp$PPg); and snpGroups, 
#' a list with 2 components giving the SNP groups construced under single-trait (snpGroups[[1]]) and multi-trait fine-mapping (snpGroups[[2]])
#' @export
#' @author Jenn Asimit
FLASHFMwithJAM <- function (beta1, corX, raf, ybar, N, save.path, TOdds = 1, 
    covY, cpp = 0.99, NCORES) 
{
    M <- length(ybar)
    Vy <- diag(covY)
    main.input <- JAMcormulti.tries(beta1, corX, raf, ybar,Vy, N, save.path)
    gc(verbose = FALSE)
    ss.stats <- summaryStats(Xmat = FALSE, ybar.all = ybar, main.input = main.input)
    fm.multi <- flashfm(main.input, TOdds = TOdds, covY, ss.stats, 
        cpp = cpp, maxmod = NULL, fastapprox = FALSE, NCORES = NCORES)
    snpGroups <- makeSNPgroups2(main.input, fm.multi, is.snpmat = FALSE, 
        min.mppi = 0.001, minsnpmppi = 0.001, r2.minmerge = 0.6)
    mpp.pp <- PPsummarise(fm.multi, snpGroups, minPP = 0.01)
    return(list(mpp.pp = mpp.pp, snpGroups = snpGroups))
}

 
#' @title Wrapper to run single-trait fine-mapping with JAMexpandedCor.multi on each trait, followed by flashfm (using fast approximation version) and then constuct SNP groups for each approach and summarises results
#' @param beta1 list where each component is a named vector of of single SNP effect estimates for a trait; one vector for each trait
#' @param corX genotype correlation matrix (reference or from sample)
#' @param raf named vector of reference allele frequencies; the name of each allele frequency is the SNP ID and MUST be in same SNP order as in corX
#' @param ybar vector of trait means; if related samples, this should be based on unrelated samples; if traits are transformed to be standard Normal, could set ybar as 0-vector
#' @param N vector of sample sizes for each trait; recommended to give effective sample sizes using GWAS summary statistics in Neff function
#' @param save.path path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").
#' @param TOdds Vector of target odds of no sharing to sharing
#' @param covY trait covariance matrix (for at most 6 traits and all traits should have a signal in the region, e.g. min p < 1E-6)
#' @param cpp cumulative posterior probability threshold for selecting top models
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @return list with 2 components: mpp.pp, a list with 4 components giving the SNP-level results (mpp.pp$PP,mpp.pp$MPP) and SNP group level results (mpp.pp$MPPg, mpp.pp$PPg); and snpGroups,
#' a list with 2 components giving the SNP groups construced under single-trait (snpGroups[[1]]) and multi-trait fine-mapping (snpGroups[[2]])
#' @export
#' @author Jenn Asimit 
FLASHFMwithJAMhat <- function (beta1, corX, raf, ybar, N, save.path, TOdds = 1, 
    covY, cpp = 0.99, NCORES) 
{
    M <- length(ybar)
    Vy <- diag(covY)
    main.input <- JAMcormulti.tries(beta1, corX, raf, ybar,Vy, N, save.path)		
    gc(verbose = FALSE)
    ss.stats <- summaryStats(Xmat = FALSE, ybar.all = ybar, main.input = main.input)
    fm.multi <- flashfm(main.input, TOdds = TOdds, covY, ss.stats, 
        cpp = cpp, maxmod = NULL, fastapprox = TRUE, NCORES = NCORES)
    snpGroups <- makeSNPgroups2(main.input, fm.multi, is.snpmat = FALSE, 
        min.mppi = 0.001, minsnpmppi = 0.001, r2.minmerge = 0.6)
    mpp.pp <- PPsummarise(fm.multi, snpGroups, minPP = 0.01)
    return(list(mpp.pp = mpp.pp, snpGroups = snpGroups))
}


