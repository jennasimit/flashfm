#' @title Approximate effective sample size from GWAS of related samples
#' @param raf vector of reference (or minor) allele frequencies for all SNPs in the GWAS
#' @param seB vector of standard errors of SNP effect estimates, in same order as SNPs in raf
#' @param Vy trait variance; default is 1 under assumption trait is transformed to a standard Normal distribution
#' @return Estimate of effective sample size 
#' @author Jenn Asimit
#' @export
 Neff <- function(raf,seB,Vy=1) {
  keep <- which(!is.na(seB))
  raf <- raf[keep]
  seB <- seB[keep]
  Vb <- seB^2
  Vx <- 2*raf*(1-raf)
  Nhat <- Vy/(Vx*Vb)
  Ne <- round(median(Nhat))
  return(Ne)
  }

  

#' @title Sample size information needed for flashfm, when samples are related and have effective sample sizes for each trait
#' @param Ne vector or trait effective sample sizes
#' @param y if available, matrix of trait measurements or indicators of non-NA trait measurements (all samples, not just unrelated);columns are traits; 
#' used to get joint sample counts; default is NULL and if not provided an approximation is used based on vector of trait sample sizes
#' @param Nsame a single sample size that is the same for all traits
#' @return list of 4 components: N = number of individuals with all traits measured; 
#' Nqq=matrix of all pair-wise counts of number of individuals with both traits in a pair measured;
#' Nq3 = vector of counts of number of individuals with three traits measured; all triples considered; NULL if M < 4
#' Nq4 = vector of counts of number of individuals with four traits measured; all quadruples considered; NULL if M < 5
#' @author Jenn Asimit
#' @export
makeNlist.rel <- function(Ne,y=NULL,Nsame=NULL){

 M <- length(Ne)

if(M==2 | M==3 | M==6) {
 Nq3 <- NULL
 Nq4 <- NULL
 }
 
 # if have trait measurements or indicators of non-NA measurements
 if(!is.null(y)) {
  Ii <- apply(y,2, function(x) !is.na(x)*1)
  Ny <- apply(Ii,2,sum)
  N <- sum(apply(Ii,1,prod))
  prop <- N/Ny
  N <- round(mean(prop*Ne)) 
  
  Nqq <- matrix(0,nrow=M,ncol=M)
  diag(Nqq) <- round(N/prop)
  
 for(i in 1:(M-1)) {
   for(j in (i+1):M) {
#   Nqq[i,j] <- Nqq[j,i] <- sum(Ii[,i]*Ii[,j])	# counts for being in both traits in a pair
	n12 <- sum(Ii[,i]*Ii[,j])
    prop <- n12/Ny[c(i,j)]
    Nqq[i,j] <- Nqq[j,i] <- round(mean(prop*diag(Nqq)[c(i,j)]))
   } } 


 if(M == 4) {
  Nq3 <- numeric(M)
  for(i in 1:M) {
#   Nq3[i] <- sum(apply(Ii[,-i],1,prod))		# counts for all but one trait
	n123 <- sum(apply(Ii[,-i],1,prod))
	prop <- n123/Ny[-i]
	Nq3[i] <- round(mean(prop*diag(Nqq)[-i]))
   }
  Nq4 <- NULL
 }

 if(M ==5) {
  Nq3 <- numeric(10)
  nc2 <- combn(1:M,2,simplify=TRUE)
  for(i in 1:10) {
  Nq3[i] <- sum(apply(Ii[,-nc2[,i]],1,prod))	# counts for all but not two traits
  prop <- Nq3[i]/Ny[-nc2[,i]]
  Nq3[i] <- round(mean(prop*diag(Nqq)[-nc2[,i]]))
  names(Nq3) <- apply(nc2,2,function(x) {paste(x,collapse=".")})
  }
  Nq4 <- numeric(5)				
  for(i in 1:M) {
   Nq4[i] <- sum(apply(Ii[,-i],1,prod))			# counts for all but one trait
   prop <- Nq4[i]/Ny[-i]
   Nq4[i] <- round(mean(prop*diag(Nqq)[-i]))
   }
  }

 } else if(!is.null(Nsame)){
     N <- Nsame
	 Nqq <- matrix(Nsame,nrow=M,ncol=M)
 	Nq3 <- rep(Nsame,choose(M,3))
 	Nq4 <- rep(Nsame,M)  
 	} else {
 	    N <- min(Ne)
 		Nqq <- diag(Ne)
 		for(i in 1:(M-1)){
 		 for(j in (i+1):M) {
 		  Nqq[i,j] <- Nqq[j,i] <- min(Ne[c(i,j)])
 		 }
 		}
 		if(M==4) {
      		Nq3 <- numeric(M)
  			for(i in 1:M) Nq3[i] <- min(Ne[-i])		# counts for all but one trait
  			Nq4 <- NULL
 				}
 		 if(M ==5) {
  			Nq3 <- numeric(10)
  			nc2 <- combn(1:M,2,simplify=TRUE)
		    for(i in 1:10) Nq3[i] <- min(Ne[-nc2[,i]])	# counts for all but not two traits
		    names(Nq3) <- apply(nc2,2,function(x) {paste(x,collapse=".")})
  			Nq4 <- numeric(5)				
  			for(i in 1:M) Nq4[i] <- min(Ne[-i])			# counts for all but one trait
 				 }		
 				 				
 	}




return(list(N=N, Nqq=Nqq, Nq3=Nq3, Nq4=Nq4))

}