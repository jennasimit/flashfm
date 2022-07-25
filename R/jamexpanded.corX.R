# from hscovar package
tagSNP <- function (mat, threshold = 0.8) 
{
    if (max(abs(mat), na.rm = T) > 1 + 1e-06) 
        stop("Correlation or R-squared matrix is expected.")
    mat <- apply(mat, 1, function(x) {
        x[!is.finite(x)] <- 0
        return(x)
    })
    diag(mat) <- 1
    if ((threshold >= 1) | (threshold <= 0)) 
        stop("Threshold out of range")
    p <- nrow(mat)
    bin <- list()
    counter <- 0
    snpset <- 1:p
    coln <- colnames(mat)
    colnames(mat) <- NULL
    repeat {
        counter <- counter + 1
        ls <- lapply(snpset, function(i) {
            id <- abs(mat[i, snpset]) > threshold
            snpset[id]
        })
        m <- which.max(lapply(ls, length))
        if (length(ls[[m]]) == 1) {
            ts <- ls[[m]]
        }
        else {
            candidate <- apply(mat[ls[[m]], ls[[m]]], 1, function(x) {
                all(abs(x) > threshold)
            })
            ts <- ls[[m]][candidate][ceiling(sum(candidate)/2)]
        }
        snps <- NULL
        bin[[counter]] <- list(snps = coln[ls[[m]]], tagsnp = coln[ts])
        snpset <- setdiff(snpset, ls[[m]])
        if ((!length(snpset) > 0) || (counter == p)) 
            break
    }
    z <- length(unlist(rlist::list.select(bin, tp = snps)))
    message(paste(z, "SNPs have been grouped into", counter, 
        "bins"))
    return(bin)
}

#' @title Thin genotype correlation matrix for JAM input
#' @param corX correlation matrix of genotypes scores
#' @param r2 r.squared threshold for thinning SNPs
#' @return correlation matrix of SNPs thinned at r2 and with colinearity removed
#' @author Jenn Asimit
#' @export
cor.refdata.fn <- function(corX,r2=0.99) {
 corX2 <- corX^2
 gmat2t<-tagSNP(corX2, threshold = r2)
 gt <- unlist(gmat2t)
 tg=gt[names(gt)=="tagsnp"]
 refG<-corX[tg,tg]
 return(refG)
 }

#' @title internal function for expanding models by tag SNPs in JAMexpandedCor.multi
#' @param snpmod intial snp model with snps separated by \code{"\%"}
#' @param taglist list of tag snps for each snp
#' @author Jenn Asimit
tagexpand.mod <- function (snpmod, taglist) {
    if (snpmod == "1") {
        df <- data.frame(str = snpmod, snps = snpmod, size = 0, 
            tag = FALSE, stringsAsFactors = FALSE)
    }
    else {
        snps <- unlist(strsplit(snpmod, "%"))
        ns <- length(snps)
        tsnps <- vector("list", ns)
        for (i in 1:ns) {
        	tsnps[[i]] <- unique(unlist(taglist[[grep(snps[i],taglist)]]))
	        }
        emods <- expand.grid(tsnps, stringsAsFactors = FALSE)
        out <- apply(emods, 1, function(x) {
            paste0(x, collapse = "%")
        })
        Imod <- which(out == snpmod)
        istag <- rep(TRUE, length(out))
        istag[Imod] <- FALSE
        df <- data.frame(str = snpmod, snps = out, size = ns, 
            tag = istag, stringsAsFactors = FALSE)
    }
    return(df)
}

##' Function to convert a correlation matrix to a covariance matrix.
##'
##' The correlation matrix to convert can be either symmetric or triangular. The covariance matrix returned is always a symmetric matrix.
##' @title Correlation Matrix to Covariance Matrix Conversion
##' @param cor.mat the correlation matrix to be converted
##' @param sd a vector that contains the standard deviations of the variables in the correlation matrix
##' @param discrepancy a neighborhood of 1, such that numbers on the main diagonal of the correlation matrix will be considered as equal to 1 if they fall in this neighborhood
##' @export
##' @note The correlation matrix input should be a square matrix, and the length of sd should be equal to the number of variables in the correlation matrix (i.e., the number of rows/columns). Sometimes the correlation matrix input may not have exactly 1's on the main diagonal, due to, eg, rounding; discrepancy specifies the allowable discrepancy so that the function still considers the input as a correlation matrix and can proceed (but the function does not change the numbers on the main diagonal). 
##' @author Ken Kelley (University of Notre Dame; (\email{KKelley@@ND.Edu}) and Keke Lai (the \code{MBESS} package), with modifications by Dustin Fife \email{fife.dustin@@gmail.com}.
cor2cov = function (cor.mat, sd, discrepancy = 0.00001) 
{
    if (dim(cor.mat)[1] != dim(cor.mat)[2]) 
        stop("'cor.mat' should be a square matrix")
    n <- sqrt(length(cor.mat))
    if (n != length(sd)) 
        stop("The length of 'sd' should be the same as the number of rows of 'cor.mat'")
    if (length(sd[sd > 0]) != n) 
        stop("The elements in 'sd' shuold all be positive")
    if (isSymmetric(cor.mat)) 
        IS.symmetric <- TRUE
    else IS.symmetric <- FALSE
    p <- dim(cor.mat)[1]
    q <- p * (p - 1)/2
    if (isTRUE(all.equal(cor.mat[lower.tri(cor.mat)], rep(0, 
        q))) || isTRUE(all.equal(cor.mat[upper.tri(cor.mat)], 
        rep(0, q)))) 
        IS.triangular <- TRUE
    else IS.triangular <- FALSE
    if (!IS.symmetric & !IS.triangular) 
        stop("The object 'cor.mat' should be either a symmetric or a triangular matrix")
    cov.mat <- diag(sd) %*% cor.mat %*% diag(sd)
    colnames(cov.mat) <- rownames(cov.mat) <- colnames(cor.mat)
    return(cov.mat)
}


#' @title Expanded version of JAM (a single-trait fine-mapping approach) that first runs on thinned SNPs and then expands models on tag SNPs; this can run independently on multiple traits
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
JAMexpandedCor.multi <- function(beta1, corX, raf, ybar, Vy, N, r2 = 0.99, save.path) 
{
    covX <- cor2cov(corX,sd=sqrt(2*raf*(1-raf)))
    
    Nlist <- makeNlist.rel(Ne = N)
    N <- diag(Nlist$Nqq)
    M <- length(beta1)
    if (is.null(names(beta1))) {
        ts <- paste0("T", 1:M)
        names(ybar) <- ts
    }
    else {
        ts <- names(ybar)
    }
#    snps <- NULL
#    for (i in 1:M) snps <- union(snps, names(beta1[[i]]))
	snps <- names(beta1[[1]])
	for (i in 2:M) snps <- intersect(snps, names(beta1[[i]]))

    snps <- intersect(snps, names(raf))
    for (i in 1:M) beta1[[i]] <- beta1[[i]][snps]
    if (is.null(colnames(corX))) {
        colnames(corX) <- names(raf)
	rownames(corX) <- names(raf)
	}
    corX <- corX[snps, snps] 
    covX <- covX[snps,snps]   
    maf <- raf * (raf <= 0.5) + (1 - raf) * (raf > 0.5)
         
    nsnps <- ncol(corX)
    refG <- cor.refdata.fn(corX, r2)
    corX2 <- corX^2
    taglist <- tagSNP(corX2, threshold = r2)     
    
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
            jam.results <- R2BGLiMS::JAM(marginal.betas = BETA, 
                trait.variance = Vy[j], cor.ref = refG, mafs.ref=mafs.ref, model.space.priors = list(a = 1, 
                  b = length(BETA), Variables = names(BETA)),  max.model.dim = 10,
                n = N[j], xtx.ridge.term = 0.01, save.path = save.path)        
        topmods <- R2BGLiMS::TopModels(jam.results, n.top.models = 1000)
 	if (is.null(ncol(topmods))) {
            stop("A single model of 10 variants was selected with PP=1. This may mean no convergence because a causal variant is missing from the data and it has no tags in your data.")
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
                covX, N = N[j], ybar = ybar[j], is.snpmat = FALSE,raf=raf)
            names(mbeta) <- expmods[, 2]
            SSy[[j]] <- Vy[j] * (N[j] - 1) + N[j] * ybar[j]^2
            Vx <- diag(covX)
            Mx <- 2 * raf
            Sxy[[j]] <- c(Sxy.hat(beta1 = beta1[[j]], Mx = Mx, 
                N = N[j], Vx = Vx, muY = ybar[j]), `1` = ybar[j] * 
                N[j])
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
                  ybar = ybar[j], is.snpmat = FALSE,raf=raf)
                mbeta <- rlist::list.append(mbeta, `1` = l1)
            }
            else {
                dd[[j]] <- data.frame(model = expmods$snps, tag = expmods$tag, 
                  lBF = lABF, stringsAsFactors = FALSE)
            }
            SM <- makesnpmod(dd[[j]], expected = 2, nsnps = nsnps)
        
        
        out$SM[[j]] <- SM
        out$mbeta[[j]] <- mbeta
        names(out$SM) <- ts
        names(out$mbeta) <- ts
        out$SM[[j]] <- best.models.cpp(out$SM[[j]],maxmod=1000)[[1]]
        out$SM[[j]] <- PP2snpmod(out$SM[[j]])
    }
    out$Nlist <- Nlist
    out$nsnps <- nsnps
    out$Gmat = covX
    out$beta1.list = beta1
    out$raf = raf

    return(out)
}


