#' @title Sample size information needed for flashfm
#' @param Nall vector of sample sizes for each trait
#' @param y if avaialable, matrix of trait measurements or indicators of non-NA trait measurements (columns are traits); 
#' used to get joint sample counts; default is NULL and if not provided an approximation is used based on vector of trait sample sizes
#' @param Nsame a single sample size that is the same for all traits
#' @return list of 4 components: N = number of individuals with all traits measured; 
#' Nqq=matrix of all pair-wise counts of number of individuals with both traits in a pair measured;
#' Nq3 = vector of counts of number of individuals with three traits measured; all triples considered; NULL if M < 4
#' Nq4 = vector of counts of number of individuals with four traits measured; all quadruples considered; NULL if M < 5
#' @author Jenn Asimit
#' @export
makeNlist <- function(Nall,y=NULL,Nsame=NULL) {
 
 M <- length(Nall)
 
 if(M==2 | M==3) {
 Nq3 <- NULL
 Nq4 <- NULL
 }
 
 # if have trait measurements or indicators of non-NA measurements
 if(!is.null(y)) {
 Ii <- apply(y,2, function(x) !is.na(x)*1)
 Ny <- apply(Ii,2,sum)
 Nqq <- matrix(0,nrow=M,ncol=M)
 diag(Nqq) <- Ny
 N <- sum(apply(Ii,1,prod))
 

for(i in 1:(M-1)) {
   for(j in (i+1):M) {
   Nqq[i,j] <- Nqq[j,i] <- sum(Ii[,i]*Ii[,j])	# counts for being in both traits in a pair
   } } 


 if(M == 4) {
  Nq3 <- numeric(M)
  for(i in 1:M) Nq3[i] <- sum(apply(Ii[,-i],1,prod))		# counts for all but one trait
  Nq4 <- NULL
 }


 if(M ==5) {
  Nq3 <- numeric(10)
  nc2 <- combn(1:M,2,simplify=TRUE)
  for(i in 1:10) Nq3[i] <- sum(apply(Ii[,-nc2[,i]],1,prod))	# counts for all but not two traits
  names(Nq3) <- apply(nc2,2,function(x) {paste(x,collapse=".")})
  Nq4 <- numeric(5)				
  for(i in 1:M) Nq4[i] <- sum(apply(Ii[,-i],1,prod))			# counts for all but one trait
  }
 } else if(!is.null(Nsame)){
     N <- Nsame
	 Nqq <- matrix(Nsame,nrow=M,ncol=M)
 	Nq3 <- rep(Nsame,choose(M,3))
 	Nq4 <- rep(Nsame,M)  
 	} else {
 	    N <- min(Nall)
 		Nqq <- diag(Nall)
 		for(i in 1:(M-1)){
 		 for(j in (i+1):M) {
 		  Nqq[i,j] <- Nqq[j,i] <- min(Nall[c(i,j)])
 		 }
 		}
 		if(M==4) {
      		Nq3 <- numeric(M)
  			for(i in 1:M) Nq3[i] <- min(Nall[-i])		# counts for all but one trait
  			Nq4 <- NULL
 				}
 		 if(M ==5) {
  			Nq3 <- numeric(10)
  			nc2 <- combn(1:M,2,simplify=TRUE)
		    for(i in 1:10) Nq3[i] <- min(Nall[-nc2[,i]])	# counts for all but not two traits
		    names(Nq3) <- apply(nc2,2,function(x) {paste(x,collapse=".")})
  			Nq4 <- numeric(5)				
  			for(i in 1:M) Nq4[i] <- min(Nall[-i])			# counts for all but one trait
 				 }		
 				 				
 	}

return(list(N=N, Nqq=Nqq, Nq3=Nq3, Nq4=Nq4))
}

##' @title Marginal PP for models sharing information between traits
##' @param STR list of models for traits 1, 2, ..., n, each given in
##'     the form of a character vector, with entries
##'     \code{"snp1\%snp2\%snp3"}. The null model is given by
##'     \code{"1"} OR \code{"0"}.  It is assumed that all elements of
##'     ABF, PP and pr below follow this same order.
##' @param PP list of posterir probabilities for the models in M
##' @param mbeta list of joint beta estimates for each trait
##' @param covY trait covariance matrix
##' @param SSy matrix of trait cross-products
##' @param Sxy matrix with each column being the cross-product between SNPs and a trait
##' @param kappa single value or vector of values to consider for the
##'     sharing scale parameter.  the value of kappa=1 must be
##'     included, and if not will be prepended.
##' @param N number of individiduals with measurements for all traits
##' @param Nqq matrix in which Nqq[i,j] = number of individuals  measured in both trait i and  trait j
##' @param nsnps number of snps in region
##' @param Mx vector of SNP means
##' @param xcovo SNP covariance matrix
##' @param Nq3  vector of counts of number of individuals with three traits measured; all triples considered; NULL if M < 4
##' @param Nq4  vector of counts of number of individuals with four traits measured; all quadruples considered; NULL if M < 5
#' @param fastapprox logical that is TRUE when fast approximation is used that does not include unequal sample size adjustments; default is FALSE
##' @return list of: - single.pp: list of pp for each model in
##'     \code{STR[[i]]} for trait i - shared.pp: list of pp for each model
##'     in \code{STR[[i]]} for trait i
						
marginalpp <- function(STR, PP, mbeta, covY, SSy, Sxy, kappa, N,Nqq,nsnps,Mx,xcovo,Nq3,Nq4,fastapprox) {  
    
    nq <- diag(Nqq)
    n <- length(STR) # number of traits
    if(n<2)
        stop("Need at least 2 traits")
    if( length(STR)!=n || length(PP)!=n )
        stop("STR and PP need to have the same lengths")
   
    if(is.null(names(STR)))
        names(STR) <- paste0("trait",seq_along(STR))
   qt <- names(STR)
    
    ## calculate model sizes 
    SS <- lapply(STR,strsplit,"%")
    usnps <- sort(unique(unlist(SS)))
    nsnpspermodel <- lapply(SS,function(x) sapply(x,length))
    for(i in seq_along(STR)) {
        wh <- which(STR[[i]] %in% c("0","1"))
        nsnpspermodel[[i]][wh] <- 0
    }
    maxsnps <- max(unlist(nsnpspermodel))
    tau <- outer(0:maxsnps,0:maxsnps,calctau,nsnps=nsnps,kappa=kappa)
   
     
    Vy <- diag(covY)

    vr <- Vres.all(Nqq,mbeta,SSy,Sxy)
    
   
    alt.pp <- calcAdjPP(qt=qt,STR=STR,SS=SS,tau=tau,nsnpspermodel=nsnpspermodel,kappa=kappa,PP=PP,beta=mbeta,SSy=SSy,Sxy=Sxy,xcovo=xcovo,Mx=Mx,N=N,allVres=vr,covY=covY,Nqq=Nqq,Nq3=Nq3,Nq4=Nq4,fastapprox)

    
    for(i in seq_along(alt.pp)){
 	names(alt.pp[[i]]) <- STR[[i]]
 	}
    ret <- lapply(seq_along(qt), function(i) {
        data.frame(single.pp=PP[[i]],
                   shared.pp=alt.pp[[i]])})
    names(ret) <- qt
    return(ret)
}


#' @title Marginal PP for models of a set of traits, sharing information between the traits
#' @param main.input List of 3 components: SM=list of snpmod objects for a set of traits; mbeta=list of joint effects for each trait; nsnps= number of SNPs in the region 
#' This could be obtained from flashfm.input or JAMexpanded.multi.
#' @param TOdds Vector of target odds of no sharing to sharing
#' @param covY trait covariance matrix
#' @param ss.stats output from summaryStats; list of 4 components: Mx = mean of SNPs, xcovo = covariance matrix of SNPs, Sxy = matrix of Sxy values (column traits), ybar=vector of trait means 
#' @param cpp cumulative posterior probability threshold for selecting top models; this is ignored when maxmod is spe$
#' @param maxmod maximum number of top models to output; NULL by default
#' @param fastapprox logical that is TRUE when fast approximation is used that does not include unequal sample size adjustments; default is FALSE
#' @return List consisting of PP: marginal PP for models and MPP: marginal PP of SNP inclusion
#' @export
#' @author Jenn Asimit
flashfm <- function(main.input,TOdds,covY,ss.stats,cpp=0.99,maxmod=NULL,fastapprox=FALSE) {
	
	Nlist <- main.input$Nlist
	Nqq <- as.matrix(Nlist$Nqq)
	Nq3 <- as.vector(Nlist$Nq3)
	Nq4 <- as.vector(Nlist$Nq4)
	N <- Nlist$N
	
	nsnps <- main.input$nsnps
	mbeta <- main.input$mbeta
	SM <- main.input$SM
	
	ybar <- ss.stats$ybar
	Sxy <- ss.stats$Sxy
	xcovo <- ss.stats$xcovo
	Mx <- ss.stats$Mx
	
	nd <- M <- length(SM)
	qt <- names(main.input$SM)    	
	kappas <- c()
	for(j in 1:length(TOdds)) kappas <- c(kappas,MFM::calckappa(nsnps=nsnps,p=2/nsnps,ndis=nd,target.odds=TOdds[j]))
    kappas <- round(kappas)
    traits <- paste(qt, collapse = "-")
    bestmod.thr <- vector("list",M)
 for(i in 1:M) bestmod.thr[[i]] <- best.models.cpp(SM[[i]],cpp.thr=cpp,maxmod)   
   
 STR <- lapply(bestmod.thr, "[[", "str") 
 PP <- lapply(bestmod.thr, "[[", "PP")

 names(STR) <- qt
 names(PP) <- qt

	SSy <- covY*(Nqq-1) + Nqq*(ybar %o% ybar) # SSy[i,j] = sum(Yi*Yj)
	
 for(i in 1:M) mbeta[[i]] <- mbeta[[i]][STR[[i]]]

    pp <- vector("list",length=nd) 
       
     for(kappa in kappas) {
     
     ret <- marginalpp(STR, PP, mbeta, covY, SSy, Sxy, kappa, N,Nqq,nsnps,Mx,xcovo,Nq3,Nq4,fastapprox)    
     for(i in 1:nd) pp[[i]] <- cbind(pp[[i]],ret[[i]]$shared.pp)
     } 
      for(i in 1:nd) {
       pp[[i]] <- cbind(ret[[i]]$single.pp,pp[[i]])
       colnames(pp[[i]]) <- paste("pp",c("null",round(TOdds,2)),sep=".")
       rownames(pp[[i]]) <- rownames(ret[[i]])
       }

   
    mpp <- lapply(pp, MPP.fn)
    names(pp) <- qt
    mpp1 <- lapply(mpp, t)
   
    MPP <- mpp1[[1]] 
    for (k in 2:M) MPP <- gtools::smartbind(MPP, mpp1[[k]], fill = 0)
    return(list(PP = pp, MPP = MPP,sharing=c("null",kappas)))
}





logminus <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}

logplus <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) + exp(y-my.max))
  return(my.res)
}



calctau <- function(n1,n2,nsnps,kappa) {
    num <- lchoose(nsnps,n1)
    denom <- logminus(logplus(lchoose(nsnps-n2,n1),lchoose(nsnps,n1) + log(kappa)),
                      lchoose(nsnps-n2,n1) + log(kappa))
    exp(num - denom)
}
    nullfirst <- function(x,wh) {
        c(x[wh],x[-wh])
    }


               


##' Internal function, logsum (copied from coloc package)
##'
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##'
##' ie, you want sum(x), but have x already stored in logs.  log(sum(exp(x))) might fail,
##' but logsum(x) should work.
##' @title logsum
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
##' @author Claudia Giambartolomei
##' @examples
##' x <- 1:10
##' log(sum(x))
##' MFM:::logsum(log(x))
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}


#' Internal function, Vx.hat
#' @param maf vector of minor allele frequencies
Vx.hat <- function(maf) 2*maf*(1-maf)

#' @title estimates cross-product of each SNP with one trait
#' @param beta1 vector of effect estimates from single SNP models
#' @param Mx vector of mean genotype counts; 2*raf
#' @param N sample size
#' @param Vx vector of genotype count variances
#' @param muY trait variance
#' @return vector of cross-products between SNPs and trait
#' @author Jenn Asimit
Sxy.hat <- function(beta1,Mx,N,Vx,muY) {
beta1*(N-1)*Vx + Mx*muY*N  
}

#' @title variance of model residuals for trait T1 at model index imod
#' @param imod model index
#' @param T1 index of trait 
#' @param SSy matrix of trait cross-products
#' @param Sxy matrix with each column being the cross-product between SNPs and a trait
#' @param mbeta list of joint beta estimates for each trait
#' @param Nqq matrix in which Nqq[i,j] = number of individuals  measured in both trait i and  trait j
#' @author Jenn Asimit
Vres.hat <- function(imod,T1,SSy,Sxy,mbeta,Nqq) {
 Syy <- SSy[T1,T1]
 beta1 <- mbeta[[T1]][[imod]]
 snps1 <- rownames(beta1)
 if(length(snps1)==0) snps1 <- names(beta1)
 y1xb1 <- sum(beta1*Sxy[snps1,T1])
 out <- (Syy - y1xb1)/Nqq[T1,T1]
 return(out) 
}

#' @title variance of model residuals for trait T1 at all models that have joint effect estimates 
#' @param Nqq matrix in which Nqq[i,j] = number of individuals  measured in both trait i and  trait j 
#' @param mbeta list of joint beta estimates for each trait
#' @param SSy matrix of trait cross-products
#' @param Sxy matrix with each column being the cross-product between SNPs and a trait
#' @author Jenn Asimit
Vres.all <- function(Nqq,mbeta,SSy,Sxy) { 
 M <- nrow(Nqq) # number of traits
 V <- structure(vector("list",M), names=rownames(SSy))
 nummods <- sapply(mbeta,length)
 for(i in 1:M) {
  V[[i]] <- sapply(1:nummods[i],Vres.hat,i,SSy,Sxy,mbeta,Nqq)
  names(V[[i]]) <- names(mbeta[[i]])
  }
 return(V)
}

#' @title covariance between residuals of a pair of models for  a trait pair
#' @param imod1 index of model from STR for trait T1, to be called from beta
#' @param imod2 index of model from STR for trait T2, to be called from beta
#' @param T1 index of first trait
#' @param T2 index of second trait
#' @param beta list of joint beta estimates for each trait
#' @param SSy matrix of trait cross-products
#' @param Sxy matrix with each column being the cross-product between SNPs and a trait
#' @param xcovo covariance matrix of c("one"=1,X)
#' @param Mx vector of SNP means
#' @param Nqq has N[i,j] = no. with  both trait i and trait j measured
#' @author Jenn Asimit
calcCres12 <-  function(imod1,imod2,T1,T2,beta,SSy,Sxy,xcovo,Mx,Nqq) { # checked and fine

 beta1 <- beta[[T1]][[imod1]]
 beta2 <- beta[[T2]][[imod2]]
 snps1 <- rownames(beta1)
 snps2 <- rownames(beta2)
 if(length(snps1)==0) snps1 <- names(beta1)
 if(length(snps2)==0) snps2 <- names(beta2)
 
 if(length(snps1)==1) { if(snps1=="1") snps1 <- "one"}
 if(length(snps2)==1) { if(snps2=="1") snps2 <- "one"} 
 snps <- union(snps1,snps2)
 SS12 <- SSy[T1,T2]
 y1xb2 <- sum(beta2*Sxy[snps2,T1])/Nqq[T1,T1]
 y2xb1 <- sum(beta1*Sxy[snps1,T2])/Nqq[T2,T2]
 xx <- (xcovo[snps,snps] + outer(Mx[snps],Mx[snps],"*"))  # approx of t(X)%*%X/Nqq[T1,T2]
 b1 <- structure(vector("numeric",length(snps)), names=snps)
 b1[snps1] <- beta1
 b2 <- structure(vector("numeric",length(snps)), names=snps)
 b2[snps2] <- beta2 
 b1xxb2 <-  t(b1)%*%xx%*%b2
 hij <- SS12/Nqq[T1,T2] - y1xb2 - y2xb1 +b1xxb2
 return(hij)
}

vcalcCres12 <- Vectorize(calcCres12,vectorize.args=c("imod1","imod2")) # calcCres12 that accepts vectors for mod1,mod2

#' @title internal function for calcAdjPP for that gives list of covariance matrix of residuals for all trait pairs
#' @param M number of traits
#' @param nummods list where component i is the number of models for trait i
#' @param beta list of joint beta estimates for each trait
#' @param SSy matrix of trait cross-products
#' @param Sxy matrix with each column being the cross-product between SNPs and a trait
#' @param xcovo covariance matrix of c("one"=1,X)
#' @param Mx vector of SNP means
#' @param Nqq has N[i,j] = no. with  both trait i and trait j measured
#' @author Jenn Asimit
allC12 <- function(M,nummods,beta,SSy,Sxy,xcovo,Mx,Nqq) {
  np <- choose(M,2)
  c2 <- combn(1:M,2,simplify=TRUE)
  pnames <- apply(c2,2,paste,collapse=".")
  Cpairs <- structure(vector("list",np),names=pnames)
  for(i in 1:np)  Cpairs[[i]] <- as.matrix(outer(1:nummods[[c2[1,i]]],1:nummods[[c2[2,i]]],vcalcCres12,c2[1,i],c2[2,i],beta,SSy,Sxy,xcovo,Mx,Nqq))
  return(Cpairs)
}



#' @title internal function for calcAdjPP for the 2-trait case
#' @param mod1 vector of model indices for trait 1
#' @param mod2 vector of model indices for trait 2
#' @param T1 index of trait 1
#' @param T2 index of trait 2
#' @param beta list of joint beta estimates for each trait
#' @param SSy matrix of trait cross-products
#' @param Sxy matrix with each column being the cross-product between SNPs and a trait
#' @param xcovo covariance matrix of c("one"=1,X)
#' @param Mx vector of SNP means
#' @param Nqq has N[i,j] = no. with  both trait i and trait j measured
#' @param Vres list of variance residuals
#' @param covY covariance matrix of traits
#' @param nsnpspermodel list of number of SNPs per model for each model in mod1, mod2
#' @return lbf12-lbf1-lbf2
#' @author Jenn Asimit
calcD12 <- function(mod1,mod2,T1,T2,beta,SSy,Sxy,xcovo,Mx,Nqq,Vres,covY,nsnpspermodel) {
 C12 <- outer(mod1,mod2,vcalcCres12,T1,T2,beta,SSy,Sxy,xcovo,Mx,Nqq)
 V1 <- Vres[[T1]][mod1]
 V2 <- Vres[[T2]][mod2]

 wh1 <- which(names(V1)=="1")
 wh2 <- which(names(V2)=="1")
 
 r1 <- C12^2/V1 # divide by V1[k] for row k
 r12 <- t(t(r1)/V2) # divide by V2[k] for col k  
 c12 <- covY[T1,T2]*(Nqq[T1,T2]-1)/Nqq[T1,T2]  # cov MLE 
 v1 <- covY[T1,T1]*(Nqq[T1,T1]-1)/Nqq[T1,T1]
 v2 <- covY[T2,T2]*(Nqq[T2,T2]-1)/Nqq[T2,T2]
 R12 <- c12^2/(v1*v2)
D12 <- -Nqq[T1,T2]*0.5*(log((1-r12)) - log((1-R12))) 
 
 return(D12)
 } 


#' @title internal function for calcAdjPP for a pair of traits
#' @param i model index for trait 1
#' @param j model index for trait 2
#' @param T1 index of trait 1
#' @param T2 index of trait 2
#' @param SS list consisting of lists of model configuration SNPs for each trait
#' @param tau matrix of adjustment terms
#' @param nsnpspermodel list of number of SNPs per model for each model in STR
#' @param kappa single value of sharing parameter kappa
#' @author Jenn Asimit
calcQ12 <- function(i,j,T1,T2,SS,tau,nsnpspermodel,kappa) {
# contributes to Q for 1 | 2 and 2|1
if(SS[[T1]][[i]][1] =="1" | SS[[T2]][[j]][1] == "1") { #at least one is null model -> tau=1 & intersection is empty
 kadj <- 1
 tij <- 1
 } else {
overlap <- 1*(any(SS[[T1]][[i]] %in% SS[[T2]][[j]]))  
kadj <- ifelse(overlap==0,1,kappa) 
tij <- tau[(nsnpspermodel[[T1]][i]+1),(nsnpspermodel[[T2]][j]+1)] # shift array indices by 1 since for numsnps 0 to maxnum
}
adj1 <- tij*kadj
adj2 <- tij*kadj
return(c(adj1,adj2))
}

vcalcQ12 <- Vectorize(calcQ12,vectorize.args=c("i","j"),SIMPLIFY=FALSE) #last arg is so that have single element output and can apply outer



#' @title internal function for calcAdjPP that gives constant term for delta
#' @param covY trait covariance matrix
#' @param Nqq has N[i,j] = no. with  both trait i and trait j measured
#' @author Jenn Asimit
calcDcon <- function(covY,Nqq) {

M <- nrow(covY)
V <- diag(covY)
Dij <- diag(M)
 
 np <- choose(M,2)
 c2 <- combn(1:M,2,simplify=TRUE)
 for(i in 1:np) {
   Dij[c2[1,i],c2[2,i]] <- covY[c2[1,i],c2[2,i]]*(Nqq[c2[1,i],c2[2,i]]-1)/Nqq[c2[1,i],c2[2,i]]/(V[c2[2,i]]*(Nqq[c2[2,i],c2[2,i]]-1)/Nqq[c2[2,i],c2[2,i]]  )
   Dij[c2[2,i],c2[1,i]] <- covY[c2[1,i],c2[2,i]]*(Nqq[c2[1,i],c2[2,i]]-1)/Nqq[c2[1,i],c2[2,i]]/(V[c2[1,i]]*(Nqq[c2[1,i],c2[1,i]]-1)/Nqq[c2[1,i],c2[1,i]]  )
   }
 
 log(det(Dij))

}



#' @title Calculates trait-adjusted posterior probabilities for all traits at sharing parameter kappa
#' @param qt vector of trait names
#' @param STR list consisting of vectors of model configurations for each trait
#' @param SS list consisting of lists of model configuration SNPs for each trait
#' @param tau matrix of adjustment terms
#' @param nsnpspermodel list of number of SNPs per model for each model in STR
#' @param kappa single value of sharing parameter kappa
#' @param PP list consisting of vectors of posterior probabilities for the model configurations for each trait
#' @param beta list of joint effect estimates for models in STR; multi.beta output
#' @param SSy matrix of trait cross-products
#' @param Sxy matrix with each column being the cross-product between SNPs and a trait
#' @param xcovo SNP covariance matrix
#' @param Mx vector of SNP means
#' @param N number of individuals with measurements for all traits
#' @param allVres list of variance residuals
#' @param covY covariance matrix of traits
#' @param Nqq matrix of all pair-wise counts of number of individuals with both traits in a pair measured;
#' @param Nq3  vector of counts of number of individuals with three traits measured; all triples considered; NULL if M < 4
#' @param Nq4  vector of counts of number of individuals with four traits measured; all quadruples considered; NULL if M < 5
#' @param fastapprox logical that is TRUE when fast approximation is used that does not include unequal sample size adjustments; default is FALSE
#' @return list of trait-adjusted posterior probabilities for each trait at sharing parameter kappa
#' @author Jenn Asimit
calcAdjPP <- function(qt,STR,SS,tau,nsnpspermodel,kappa,PP,beta,SSy,Sxy,xcovo,Mx,N,allVres,covY,Nqq,Nq3,Nq4,fastapprox) {
 
    M <- length(qt)
    np <- choose(M,2)
    c2 <- combn(1:M,2,simplify=TRUE)
    c2names <- apply(c2,2, function(cc) return(paste0("Q",paste(cc,collapse=".Q"))))
	Q <- structure(vector("list",np),names=c2names)
	nummods <- sapply(STR,length)

	for(i in 1:np) { # for each qt pair Q[[i]] is a matrix where Q[[i]][j,k] is a list with 
					# two components adjPP12[modj for T1,modk for T2], adjPP21[modj for T1,modk for T2] where 1=c2[1,i], 2=c2[2,i]
	    
     Q[[i]] <- outer(1:nummods[c2[1,i]],1:nummods[c2[2,i]],vcalcQ12,T1=c2[1,i],T2=c2[2,i],SS,tau,nsnpspermodel,kappa)
      

		}
	
	
	qns <- unlist(strsplit(names(Q),"[.]"))
	PPadj <- structure(vector("list",M),names=qt)
	if(M==2) { 
	  i=1
	  delta <- calcD12(1:nummods[c2[1,i]],1:nummods[c2[2,i]],T1=c2[1,i],T2=c2[2,i],beta=beta,SSy=SSy,Sxy=Sxy,xcovo=xcovo,Mx=Mx,Nqq=Nqq,Vres=allVres,covY=covY,nsnpspermodel)
      tmp <- Q[[i]]
      q12 <- apply(tmp,2,function(x) unlist(lapply(x,"[[",1))) 
 	  q21 <- apply(tmp,2,function(x) unlist(lapply(x,"[[",2)))  	    
 	  
 	  # need sum(delta*PP) = 1
 	  pd1 <- delta + matrix(log(PP[[2]]),nrow=nummods[1],ncol=nummods[2],byrow=TRUE)
 	  pd1 <- t(apply(pd1,1,function(x) x-logsum(x)))
 	  q1 <- log(q12)+pd1; q1 <- apply(q1,1,logsum); q1 <- exp(q1-logsum(q1))
 	  
 	  pd2 <- t(delta) + matrix(log(PP[[1]]),nrow=nummods[2],ncol=nummods[1],byrow=TRUE)
 	  pd2 <- t(apply(pd2,1,function(x) x-logsum(x)))
 	  q2 <- t(log(q21))+pd2; q2 <- apply(q2,1,logsum); q2 <- exp(q2-logsum(q2))
 	
	  PPadj[[1]] <- PP[[1]]*q1/sum(PP[[1]]*q1)
	  PPadj[[2]] <- PP[[2]]*q2/sum(PP[[2]]*q2)
	  
	} else {
		Dcon <- calcDcon(covY,Nqq) 
		Cij <- allC12(M,nummods,beta,SSy,Sxy,xcovo,Mx,Nqq)
		lPP <- lapply(PP,log)
		
		for(i in 1:M) { # for each trait
	    qn  <- paste0("Q",i)
	 	ind <- grep(qn,qns,fixed=TRUE)
	 	whO <- ind[which(ind %% 2 == 1)] # odd indices so first list component 	 
	 	whE <- ind[which(ind %% 2 == 0)]
	 	keep <- NULL
	 	if(length(whO)>0) {
	 	  tmp <- Q[(whO+1)/2]
	 	  keep <- lapply(tmp,function(x) apply(x,2,function(y) unlist(lapply(y,"[[",1))) )
	 	  nk <- names(keep)
	 	  names(keep) <- unlist(strsplit(nk,"[.]"))[c(FALSE,TRUE)]
	 		}
	 	if(length(whE)>0) {
	 	 tmp <- Q[whE/2]
	 	 if(!is.null(keep)) {
	 	 keep2 <- lapply(tmp,function(x) apply(x,2,function(y) unlist(lapply(y,"[[",2))) )
	 	 nk <- names(keep2)
	 	 names(keep2) <- unlist(strsplit(nk,"[.]"))[c(TRUE,FALSE)]
	 	 keep2 <- lapply(keep2,t)
	 	 keep <- append(keep,keep2)	# 2nd component in list pair
	 	 				} else {
	 	 				keep <- lapply(tmp,function(x) apply(x,2,function(y) unlist(lapply(y,"[[",1))) )
	 	 				nk <- names(keep)
	 	 				names(keep) <- unlist(strsplit(nk,"[.]"))[c(TRUE,FALSE)]
	 	 				keep <- lapply(keep,t)
	 	 				}
	 	} # row k of keep corresponds to model k of trait i
	 	# qp  PP adjustment for trait i, sum over all models to get weighted PP wrt to trait i and multiplying by D for each model 
	 	keep <- keep[sort(names(keep),decreasing=FALSE)]
	 	keep <- lapply(keep,log)
	 	keep <- lapply(keep,as.matrix)

	 	
	 	if(M==3) {
	 	
	 	  if(length(STR[[i]])>1) {
	 	  Vind <- c(i,setdiff(1:3,i))
	 	  Cind <- combn(Vind,2,simplify=TRUE)
	 	  pnames <- apply(Cind,2,function(x) {x <- sort(x); paste(x,collapse=".")})
	 	  Ldcon12 <- numeric(M)
	 	  for(l in 1:M) Ldcon12[l] <- calcDcon(covY[Cind[,l],Cind[,l]],Nqq[Cind[,l],Cind[,l]])
	 	  ccind <- apply(Cind,2,sort)
	 	  c2ind <- Ctrans(Vind,ccind)
	 	  
		 Nsame <-1; 
		 if(var(diag(Nqq))==0 | fastapprox) Nsame <- 0 
		 PPadj[[i]] <- ppadjT3(N, nummods[Vind], allVres[Vind], Cij[pnames], Dcon, keep,Nqq[Vind,Vind],Ldcon12,c2ind-1,lPP[Vind],Nsame)
	 	 	} else { PPadj[[i]] <- 1 }
	 	 }
	 	 
	 	
	 	
	 	if(M==4) {
	 	 if(length(STR[[i]])>1) {
	 	  Vind <- c(i,setdiff(1:4,i))
	 	  Cind <- combn(Vind,2,simplify=TRUE)
	 	  pnames <- apply(Cind,2,function(x) {x <- sort(x); paste(x,collapse=".")})
	 	  npair <- length(pnames)
	 	  CijI <- vector("list",npair)
	 	  for(l in 1:npair) {
	 	    ij <- Cind[,l]
	 	    if(ij[1]>ij[2]) {CijI[[l]] <- t(Cij[[paste(sort(ij),collapse=".")]])
	 	    } else { CijI[[l]] <- Cij[[paste(ij,collapse=".")]] }
	 	    names(CijI)[l] <- paste(ij,collapse=".")
	 	  }
	 	  
	 	  Ldcon12 <- numeric(npair)
	 	  for(l in 1:npair) Ldcon12[l] <- calcDcon(covY[Cind[,l],Cind[,l]],Nqq[Cind[,l],Cind[,l]])
	 	  Ldcon123 <- numeric(4)
	 	  CijI.3 <-vector("list",4)
	 	  for(l in 1:4) {
	 	   Ldcon123[l] <- calcDcon(covY[Vind,Vind][-l,-l],Nqq[Vind,Vind][-l,-l])
	 	   krm <- grep(Vind[l],names(CijI))
	 	   CijI.3[[l]] <- CijI[-krm]
	 	 }
	 	 
	 	 Nsame <-1; 
		 if(var(diag(Nqq))==0 | fastapprox) Nsame <- 0 
		
		 PPadj[[i]] <- ppadjT4(N, nummods[Vind], allVres[Vind], CijI, Dcon, keep,Nqq[Vind,Vind],Ldcon12,Nq3[pnames],Ldcon123,CijI.3,lPP[Vind],Nsame)	
	 	 	} else { PPadj[[i]] <- 1 }
	 	 }
	 	 
	 	  	 if(M==5) {
	 	 if(length(STR[[i]])>1) {
	 	  Vind <- c(i,setdiff(1:5,i))
	 	  Cind <- combn(Vind,2,simplify=TRUE)
	 	  pnames <- apply(Cind,2,function(x) {x <- sort(x); paste(x,collapse=".")})
	 	  npair <- length(pnames)
	 	  CijI <- vector("list",npair)
	 	  for(l in 1:npair) {
	 	    ij <- Cind[,l]
	 	    if(ij[1]>ij[2]) {CijI[[l]] <- t(Cij[[paste(sort(ij),collapse=".")]])
	 	    } else { CijI[[l]] <- Cij[[paste(ij,collapse=".")]] }
	 	    names(CijI)[l] <- paste(ij,collapse=".")
	 	  }
	 	  
	 	  Ldcon12 <- numeric(npair)
	 	  for(l in 1:npair) Ldcon12[l] <- calcDcon(covY[Cind[,l],Cind[,l]],Nqq[Cind[,l],Cind[,l]]) # dcon for each pair
	 	  Ldcon123 <- numeric(npair)
	 	  CijI.3 <-vector("list",npair)
	 	  for(l in 1:npair) {
	 	   Ldcon123[l] <- calcDcon(covY[Vind,Vind][-Cind[,l],-Cind[,l]],Nqq[Vind,Vind][-Cind[,l],-Cind[,l]]) # dcon for each triple (leave out a pair)
	 	   krm <- union(grep(Vind[Cind[1,l]],names(CijI)),grep(Vind[Cind[2,l]],names(CijI)))
	 	   CijI.3[[l]] <- CijI[-krm]
	 	    }
	 	    Ldcon1234 <- numeric(5)
	 	    CijI.4 <- vector("list",5)
	 	   for(l in 1:5){
	 	   Ldcon1234[l]  <- calcDcon(covY[Vind,Vind][-l,-l],Nqq[Vind,Vind][-l,-l]) # dcon for quadruples (leave one out)
	 	   krm <- grep(Vind[l],names(CijI))
	 	   CijI.4[[l]] <- CijI[-krm]
	 	   }
	 	  
		Nsame <-1; 
		if(var(diag(Nqq))==0 | fastapprox) Nsame <- 0
		 PPadj[[i]] <- ppadjT5(N, nummods[Vind], allVres[Vind], CijI, Dcon, keep,Nqq[Vind,Vind],Ldcon12,Nq3[pnames],Ldcon123,CijI.3,Nq4[Vind],Ldcon1234,CijI.4,lPP[Vind],Nsame)	
	 	 	} else { PPadj[[i]] <- 1 }
	 	 }
	 	 
	 	 
	 
	 	
	 		 	
	 	}
     }
	return(PPadj) 
	 
}





Ctrans <- function(Vind,Cind) {
 M <- length(Vind)
 ccind <- Cind
 for(j in 1:M) {
  ccind[1,][which(Cind[1,]==Vind[j])] <- j
  ccind[2,][which(Cind[2,]==Vind[j])] <- j
  }
  return(ccind)
}

# from GUESSFM by Chris Wallace
makestr <- function(x) {paste(sort(unique(x)),collapse="%")}

#' @title Best models from a snpmpd object by cpp or maximum number of models - modification of best.models from GUESSFM by Chris Wallace
#' @param d snpmod object
#' @param cpp.thr cumulative posterior probability threshold for selecting top models; this is ignored when maxmod is specified
#' @param maxmod maximum number of top models to output; NULL by default
#' @return data.frame of top SNP models ordered by posterior probability (PP); if maxmod is specified the PPs are re-scaled to sum to 1.
#' @export
best.models.cpp <- function (d, cpp.thr = .99,maxmod=NULL) 
{
#    if (is.list(d)) 
#        return(lapply(d, best.models, pp.thr = pp.thr, cpp.thr = cpp.thr))
#    if (!is(d, "snpmod")) 
#        stop("expected d to be a snpmod object, found ", class(d))
    
        if(is.null(maxmod)){
        cpp <- cumsum(d@models$PP)
        wh <- which(cpp <= cpp.thr)
        if (!length(wh))  wh <- 1
        wh <- c(wh, max(wh) + 1)
                
        d@models <- d@models[wh, ]
        d@models$PP <- d@models$PP/sum(d@models$PP)
        }
        
        if(!is.null(maxmod)){
        d@models <- d@models[order(d@models$PP, decreasing = TRUE) , ][1:min(maxmod,nrow(d@models)),]
        message("Before adjustment, CPP of top",maxmod,"models is", sum(d@models$PP))
        d@models$PP <- d@models$PP/sum(d@models$PP)
        }
           
    out <- cbind(d@models, snps = unlist(lapply(strsplit(d@models$str, "%"), makestr)))
     
    return(out)
}

