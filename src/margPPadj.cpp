#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
String concat(CharacterVector x) {
//e.g.  concat(c(1,".",2))

  // initialize string output stream
  std::ostringstream ossOut;

  int nChar = x.size();

// concatenate input
  for (int i = 0; i < nChar; i++)
    ossOut << x[i];

  return ossOut.str();  
}


// [[Rcpp::export]]
NumericVector calcDij(const int i, const int j, 
					const NumericMatrix& Cr12, 
					const NumericVector& Vr1,
			   		const NumericVector& Vr2) {
	// need to have smaller trait index listed first
//	if(i > j) { 
//		int tmp = i;
//		i = j;
//		j = tmp;
//	}
	
	NumericVector Dij(2);
	Dij[0] = Cr12(i,j)/Vr2(j);
	Dij[1] = Cr12(i,j)/Vr1(i);
//	Dij[0] = Cr12(i,j);
//	Dij[1] = Cr12(i,j);
	return(Dij);
}	





// [[Rcpp::export]]
NumericMatrix pairs(const int M) {
// returns matrix of 2 rows and M*(M-1)/2 columns; each column is a pair of indices
	Function f("combn");
//	NumericVector inds = indices(M); 
	IntegerVector inds = seq(0,M-1);
	return f(_["x"]=inds,_["m"]=2);
}

// [[Rcpp::export]]
NumericVector indices(const int M) {
	Function f("seq");
	return f(_["from"]=0,_["to"]=M-1,_["by"]=1);
}

// for log scale x, returns x-logsum(x) i.e. so sum(exp(x))=1
// Rccp::export
NumericVector logsum1(const NumericVector& x) {
  const int n=x.size();
  double lsum=0.0;
  double maxX = max(x);
  NumericVector y(n);
  for(int i=0; i<n; i++) 
    lsum += exp(x(i)-maxX);
    double lsumx = maxX + log(lsum); //logsum(x) =log(sum(x))
  for(int i=0; i<n; i++)
    y(i)=x(i)-lsumx;
  return(y);
}

// re-scale log scale x so sum(exp(x))=1 then add K
// Rccp::export
NumericVector logsum1K(const NumericVector& x, const NumericVector& K) {
  const int n=x.size();
  double lsum=0.0;
  double maxX = max(x);
  NumericVector y(n);
  for(int i=0; i<n; i++) 
    lsum += exp(x(i)-maxX);
    double lsumx = maxX + log(lsum); //logsum(x) =log(sum(x))
  for(int i=0; i<n; i++)
    y(i)=x(i)-lsumx + K(i);
  return(y);
}


// [[Rcpp::export]]
NumericMatrix Mlogsum1(const NumericMatrix& x) {
  const int nrows = x.nrow();
  const int ncols = x.ncol();
  double lsum=0.0;
  double maxX = max(x);
  NumericMatrix y(nrows,ncols);
  for(int i=0; i<nrows; i++) {
  	for(int j=0; j<ncols; j++) {
    lsum += exp(x(i,j)-maxX);
    }
  }  
    double lsumx = maxX + log(lsum); //logsum(x) =log(sum(x))
  	for(int i=0; i<nrows; i++) {
  		for(int j=0; j<ncols; j++) {
   		 y(i,j)=x(i,j)-lsumx;
   		 }
   }		 
  return(y);
}



double logsum(const NumericVector& x) {
  const int n=x.size();
  double lsum=0.0;
  double maxX = max(x);
  for(int i=0; i<n; i++) 
    lsum += exp(x(i)-maxX);
  double lsumx = maxX + log(lsum); //logsum(x) =log(sum(x))
  return(lsumx);
}




// [[Rcpp::export]]
double Mlogsum(const NumericMatrix& x) {
  const int nrows = x.nrow();
  const int ncols = x.ncol();
  double lsum=0.0;
  double maxX = max(x);
  for(int i=0; i<nrows; i++) {
  	for(int j=0; j<ncols; j++) {
     lsum += exp(x(i,j)-maxX);
     } 
    }
  double lsumx = maxX + log(lsum); //logsum(x) =log(sum(x))
  return(lsumx);
}







// [[Rcpp::export]]
double calcDelta(const IntegerVector& imods, // vector of model indices for M traits, in same order as traits
						const int N, const int M, const double dcon,
						const List Vr, const List Cr,
						const NumericMatrix& c2) {
// N is the number measured for all M traits						
	
	int T1, T2;
	int mod1, mod2;
	
	NumericVector dd(2);
	double val;
	double sign;
	
	mat Dij(M,M,fill::eye); // identity matrix size M
	
	int np = M*(M-1)/2; // number of pairs
//	NumericMatrix c2 = pairs(M);

	for(int k=0; k<np ; ++k) {
		T1 = c2(0,k);
		T2 = c2(1,k);
		
		mod1 = imods(T1);
		mod2 = imods(T2);
	
		NumericMatrix Cr12 = Cr[k]; // Cij matrices are listed in same order as pairs		  
		NumericVector V1 = Vr[T1];
		NumericVector V2 = Vr[T2];
		dd = calcDij(mod1, mod2, Cr12, V1, V2);
		
		Dij(T1,T2) = dd(0);
		Dij(T2,T1) = dd(1);
		}
		log_det(val, sign, Dij);
		double lDD = val*sign;
//		double out = exp(-N*0.5*(lDD-dcon));	
		double out = -N*0.5*(lDD-dcon);
		return out;
	}


// [[Rcpp::export]]
double calcQD3(const int mod1, const int N,
					const NumericVector& nummods,
					const List Vr, const List Cr, const double dcon,
					const List keep, const NumericMatrix& Nqq,
					const NumericVector& Ldcon12,
					const NumericMatrix& c2,
					const List lPP,const int Nsame) {		
	
	const int T1 = 0;
	const int M = 3;
	
	
	IntegerVector imods(M), im2(M-1);
	
//	const NumericVector& tmp = indices(M);
//	NumericVector notT1 = tmp[ tmp != T1 ];
//	int T2 = notT1(0);
//	int T3 = notT1(1);
	const int T2 = 1;
	const int T3 = 2;
	const int mT2 = nummods(T2);
	const int mT3 = nummods(T3);

	const int Nc = mT2*mT3;
	NumericVector QD(Nc);
	NumericVector Kterm(Nc);	
	const NumericVector& lPP2 = lPP[1];
	const NumericVector& lPP3 = lPP[2];

	
	
	const NumericMatrix& K1 = keep[0];
	const NumericMatrix& K2 = keep[1]; // keep is in order and on log scale
	
	double  ddd,nn;

	const NumericMatrix& c1p = pairs(2);
	
	const int np = M*(M-1)/2; // number of pairs
			
		// calculate all Dijk*PPj*PPk (log sum) for given i and adjust to sum to 1
		imods(T1) = mod1;
		int i = 0;  // each i is for a pair (j,k) models for T2,T3
		for(int j = 0; j < mT2; j++) {
			imods(T2) = j;
			for(int k = 0; k < mT3; k++) {
				imods(T3) = k;				
				QD(i) = calcDelta(imods, N, M, dcon, Vr, Cr, c2) + lPP2(j) + lPP3(k); //need sum(deltajk*PPj*PPk)=1
				Kterm(i) = K1(mod1,j)+K2(mod1,k);
	
		if(Nsame > 0) {
	// calculate d12 adjustment term between all pairs			
				
		for(int l=0; l<np ; ++l) {
			int Tt1 = c2(0,l);
			int Tt2 = c2(1,l);
			int mod1 = imods(Tt1);
			int mod2 = imods(Tt2);
			int N12 = Nqq(Tt1,Tt2);	
			
			if(N12 > N) {
				double nnum = N12-N;
				nn=nnum/N12; 
				List Vr2 = List(2);
				Vr2[0] = Vr[Tt1];
				Vr2[1] = Vr[Tt2];
				IntegerVector imods2(2);
				imods2(0)= mod1;
				imods2(1)= mod2;
				List Cr2l = List(1);
				Cr2l[0] = Cr[l];	// Cij matrices are listed in same order as pairs	
				double dcon12 = Ldcon12(l);
				ddd=calcDelta(imods2, N12, 2, dcon12, Vr2, Cr2l,c1p);
				QD(i) += nn*ddd;
				}
			}
		}		
				i += 1;				
			}
		}



	QD = logsum1K(QD,Kterm); // sum(exp(qd))=1 then add K 
			
	double out = logsum(QD);				

return(out);			
}


// [[Rcpp::export]]
NumericVector ppadjT3(const int N,
					const NumericVector& nummods,
					const List Vr, const List Cr, const double dcon,
					const List keep, const NumericMatrix& Nqq,
					const NumericVector& Ldcon12,
					const NumericMatrix& c2,
					const List lPP, const int Nsame) {	 	
					
	const int mT1 = nummods(0);
	const NumericVector& lPP1 = lPP[0];
	NumericVector QD(mT1);
	NumericVector lppadj(mT1);
	
	for(int i=0; i<mT1; i++) {
		QD(i) = calcQD3(i, N,nummods,Vr, Cr, dcon,keep,  Nqq, Ldcon12, c2, lPP, Nsame);
	}

//	NumericVector lppadj = logsum1K(QD,lPP1);
	 	 
//	 const NumericVector& lppadj1 = logsum1(lppadj);
//	 const NumericVector& ppadj = exp(lppadj1);

	const double lsum = logsum(QD);

	for(int i = 0; i<mT1; i++) {
		lppadj(i) = QD(i)-lsum + lPP1(i);	
	 	}
	 	 
	 const NumericVector& lppadj1 = logsum1(lppadj);
	 const NumericVector& ppadj = exp(lppadj1);
	 					
	 					
return(ppadj);					
}					



// [[Rcpp::export]]
double calcQD4(const int mod1, const int N,
					const NumericVector& nummods,
					const List Vr, const List Cr, const double dcon,
					const List keep, const NumericMatrix& Nqq, 
					const NumericVector& Ldcon12,
					const NumericVector& Nq3,
					const NumericVector& Ldcon123, const List Cr3,
					const List lPP, const int Nsame) {		
	
	const int T1 = 0;
	const int M = 4;
	
	
	IntegerVector imods(M), im2(M-1);
	
	const NumericVector& tmp = indices(M);
	const int T2 = 1;
	const int T3 = 2;
	const int T4 = 3;
	const int mT2 = nummods(T2);
	const int mT3 = nummods(T3);
	const int mT4 = nummods(T4);



	const int Nc = mT2*mT3*mT4;
	NumericVector QD(Nc);
	NumericVector Kterm(Nc);
	
	const NumericVector& lPP2 = lPP[1];
	const NumericVector& lPP3 = lPP[2];
	const NumericVector& lPP4 = lPP[3];
	
	
	
	const NumericMatrix& K1 = keep[0];
	const NumericMatrix& K2 = keep[1]; // keep is in order and on log scale
	const NumericMatrix& K3 = keep[2];
	
	double ddd,nn;

	
	const int np = M*(M-1)/2; // number of pairs
	const NumericMatrix& c2 = pairs(4);		
	const NumericMatrix& c23 = pairs(3);
	const NumericMatrix& c1p = pairs(2);

		// calculate all Dijkl*PPj*PPk*PPl (log sum) for given i and adjust to sum to 1
		imods(T1) = mod1;
		int i = 0;  // each i is for a triple (j,k,l) models for T2,T3, T4

		for(int j = 0; j < mT2; j++) {
			imods(T2) = j;
			for(int k = 0; k < mT3; k++) {
				imods(T3) = k;	
				for(int l = 0; l < mT4; l++) {
					imods(T4) = l;					
				QD(i) = calcDelta(imods, N, M, dcon, Vr, Cr, c2) + lPP2(j) + lPP3(k) +lPP4(l);
				Kterm(i) = K1(mod1,j)+K2(mod1,k)+K3(mod1,l);


	
	// calculate d12 adjustment term between all pairs			
		if(Nsame>0) {
				
		for(int l=0; l<np ; ++l) {
			int Tt1 = c2(0,l);
			int Tt2 = c2(1,l);
			int mod1 = imods(Tt1);
			int mod2 = imods(Tt2);
			NumericVector notT12 = tmp[ (tmp != Tt1) & (tmp != Tt2)];
			int nTt1 = notT12(0);
			int nTt2 = notT12(1);
			int N12 = Nqq(Tt1,Tt2) - Nq3(nTt1) - Nq3(nTt2) + N ;	
			
			if(N12 > 0) {
				double nnum=N12;
				nn=nnum/Nqq(Tt1,Tt2); 
				List Vr2 = List(2);
				Vr2[0] = Vr[Tt1];
				Vr2[1] = Vr[Tt2];
				
				IntegerVector imods2(2);
				imods2(0)= mod1;
				imods2(1)= mod2;
				List Cr2l = Cr[l];	// Cij matrices are listed in same order as pairs	
				double dcon12 = Ldcon12(l);
				int nq12 = Nqq(Tt1,Tt2);
				ddd=calcDelta(imods2, nq12, 2, dcon12, Vr2, Cr2l,c1p);
				QD(i) += nn*ddd;
				}
			}
	
	
	// calculate d123 adjustment term between all triples
		
		
		for(int l=0; l<M; ++l) {
			NumericVector notT1 = tmp[ tmp != l];
			int nTt1 = notT1(0);
			int nTt2 = notT1(1);
			int nTt3 = notT1(2);
	
			int N123 = Nq3(l) -  N ;	
			
			if(N123 > 0) {
				double nnum=N123;
				nn=nnum/Nq3(l);
				List Vr3 = List(3);
				Vr3[0]=Vr[nTt1];
				Vr3[1]=Vr[nTt2];
				Vr3[2]=Vr[nTt3];
				
				IntegerVector imods3(3);
				imods3(0)=imods(nTt1);
				imods3(1)=imods(nTt2);
				imods3(2)=imods(nTt3);
				List Cr3l = Cr3[l];
				double dcon123 = Ldcon123(l);
				int nq123 = Nq3(l);
				ddd=calcDelta(imods3,nq123 , 3, dcon123, Vr3, Cr3l,c23);
			 	QD(i) += nn*ddd;
			}
		
		
		}
		} // loop above run if different N among traits						
		i += 1;	// increment for each (j,k,l) triple	
			} // l
		}  // k  
		}  // j


  
 	QD = logsum1K(QD,Kterm); // sum(exp(qd))=1 then add K 

			
	double out = logsum(QD);			


return(out);			
}


// [[Rcpp::export]]
NumericVector ppadjT4(const int N,
					const NumericVector& nummods,
					const List Vr, const List Cr, const double dcon,
					const List keep, const NumericMatrix& Nqq, 
					const NumericVector& Ldcon12,
					const NumericVector& Nq3,
					const NumericVector& Ldcon123, const List Cr3,
					const List lPP, const int Nsame) {	 	
					
	const int mT1 = nummods(0);
	const NumericVector& lPP1 = lPP[0];
	NumericVector QD(mT1);
//	NumericVector lppadj(mT1);
	
	for(int i=0; i<mT1; i++) {
		QD(i) = calcQD4(i, N,nummods,Vr, Cr, dcon,keep,  Nqq, Ldcon12, Nq3,
					Ldcon123,  Cr3, lPP, Nsame);
	}
	
//	const double lsum = logsum(QD);

//	for(int i = 0; i<mT1; i++) {
//		lppadj(i) = QD(i)-lsum + lPP1(i);	
//	 	}

	NumericVector lppadj = logsum1K(QD,lPP1);
	 	 
	 const NumericVector& lppadj1 = logsum1(lppadj);
	 const NumericVector& ppadj = exp(lppadj1);
	 					
return(ppadj);					
}					



// [[Rcpp::export]]
double calcQD5(const int mod1, const int N,
					const NumericVector& nummods,
					const List Vr, const List Cr, const double dcon,
					const List keep, const NumericMatrix& Nqq,
					const NumericVector& Ldcon12,
					const NumericVector& Nq3,
					const NumericVector& Ldcon123, const List Cr3,
					const NumericVector& Nq4, const NumericVector& Ldcon1234, 
					const List Cr4, const List lPP, const int Nsame) {		
	
	const int T1 = 0;
	const int M = 5;
	
	IntegerVector imods(M), im2(M-1);
	
	
	const int T2 = 1;
	const int T3 = 2;
	const int T4 = 3;
	const int T5 = 4;
	const int mT2 = nummods(T2);
	const int mT3 = nummods(T3);
	const int mT4 = nummods(T4);
	const int mT5 = nummods(T5);

	const int Nc = mT2*mT3*mT4*mT5;
	NumericVector QD(Nc);
	NumericVector Kterm(Nc);

	
	
	const NumericVector& lPP2 = lPP[1];
	const NumericVector& lPP3 = lPP[2];
	const NumericVector& lPP4 = lPP[3];
	const NumericVector& lPP5 = lPP[4];

	int Tt1, Tt2;
	
	const NumericMatrix& K1 = keep[0];
	const NumericMatrix& K2 = keep[1]; // keep is in order and on log scale
	const NumericMatrix& K3 = keep[2];
	const NumericMatrix& K4 = keep[3];
	
//	double delta, d12, ddd,nn, d123, d1234;
	double nn, ddd;

//	NumericVector dd(2);
//	double val;
//	double sign;
	
//	mat Dij(2,2,fill::eye); // identity matrix size M
	const int np = M*(M-1)/2; // number of pairs
	const NumericMatrix& c2 = pairs(M);		
	const NumericVector& tmp = indices(M);
	const NumericMatrix& c234 = pairs(3);
	const NumericMatrix& c2345 = pairs(4);
	const NumericMatrix& c1p = pairs(2);
	
		imods(T1) = mod1;
		int i = 0;  // each i is for a (j,k,l,r) models for T2,T3,T4,T5
		
		for(int j = 0; j < mT2; j++) {
			imods(T2) = j;
			for(int k = 0; k < mT3; k++) {
				imods(T3) = k;	
				for(int l = 0; l < mT4; l++) {
					imods(T4) = l;
					for(int r = 0; r < mT5; r++) {	
					imods(T5) = r;		
						
				QD(i) = calcDelta(imods, N, M, dcon, Vr, Cr, c2) + lPP2(j) + lPP3(k) +lPP4(l) +lPP5(r);
				Kterm(i) = K1(mod1,j)+K2(mod1,k)+K3(mod1,l)+K4(mod1,r);

	// calculate d12 adjustment term between all pairs	if have missing trait measurements		
			
		if(Nsame >0){
				
		for(int l=0; l<np ; ++l) {
			Tt1 = c2(0,l);
			Tt2 = c2(1,l);
			int mod1 = imods(Tt1);
			int mod2 = imods(Tt2);
			
			NumericVector notT12 = tmp[ (tmp != Tt1) & (tmp != Tt2)];
			int nTt1 = notT12(0);
			int nTt2 = notT12(1);
			int nTt3 = notT12(2);
			int N12 = Nqq(Tt1,Tt2) - Nq3(nTt1) - Nq3(nTt2) -Nq3(nTt3) + Nq4(nTt1) + Nq4(nTt2) + Nq4(nTt3)-N ;	
			
			if(N12 > 0) {
				double nnum=N12;
				nn=nnum/Nqq(Tt1,Tt2); 
				List Vr2 = List(2);
				Vr2[0] = Vr[Tt1];
				Vr2[1] = Vr[Tt2];
				IntegerVector imods2(2);
				imods2(0)= mod1;
				imods2(1)= mod2;
				List Cr2l = Cr[l];	// Cij matrices are listed in same order as pairs	
				double dcon12 = Ldcon12(l);
				int nq12 = Nqq(Tt1,Tt2);
				ddd=calcDelta(imods2, nq12, 2, dcon12, Vr2, Cr2l,c1p);
				QD(i) += nn*ddd;
				}
			}


	// calculate d123 adjustment term between all triples
		
		for(int l=0; l<np ; ++l) {
			Tt1 = c2(0,l);
			Tt2 = c2(1,l);
			
			NumericVector notT12 = tmp[ (tmp != Tt1) & (tmp != Tt2)]; // the 3 that are not in the pair are the focus
			int nTt1 = notT12(0);
			int nTt2 = notT12(1);
			int nTt3 = notT12(2);
			int N123 = Nq3(l) - Nq4(Tt1) - Nq4(Tt2) + N ;	
			
			if(N123 > 0) {
				double nnum=N123;
				nn=nnum/Nq3(l); 
				List Vr3 = List(3);
				Vr3[0]=Vr[nTt1];
				Vr3[1]=Vr[nTt2];
				Vr3[2]=Vr[nTt3];
				
				IntegerVector imods3(3);
				imods3(0)=imods(nTt1);
				imods3(1)=imods(nTt2);
				imods3(2)=imods(nTt3);
				List Cr3l = Cr3[l];
				double dcon123 = Ldcon123(l);
				int nq123 = Nq3(l);
				ddd=calcDelta(imods3,nq123 , 3, dcon123, Vr3, Cr3l,c234);		
				QD(i) += nn*ddd;
				}
			}

	
	// calculate d1234 adjustment term between all quadruples
		
		
		for(int l=0; l<M; ++l) {
			NumericVector tmp = indices(M);
			NumericVector notT1 = tmp[ tmp != l];
			int nTt1 = notT1(0);
			int nTt2 = notT1(1);
			int nTt3 = notT1(2);
			int nTt4 = notT1(3);
			int N1234 = Nq4(l) -  N ;	
			
			if(N1234 > 0) {
				double nnum=N1234;
				nn=nnum/Nq4(l);
				List Vr4 = List(4);
				Vr4[0]=Vr[nTt1];
				Vr4[1]=Vr[nTt2];
				Vr4[2]=Vr[nTt3];
				Vr4[3]=Vr[nTt4];
				
				IntegerVector imods4(4);
				imods4(0)=imods(nTt1);
				imods4(1)=imods(nTt2);
				imods4(2)=imods(nTt3);
				imods4(3)=imods(nTt4);
				List Cr4l = Cr4[l];
				double dcon1234 = Ldcon1234(l);
				int nq1234 = Nq4(l);
				ddd=calcDelta(imods4,nq1234 , 4, dcon1234, Vr4, Cr4l,c2345);
			 	QD(i) += nn*ddd;
			}			
		}
		
	}
			 
			i += 1;			
					
			}
			
			
		} 
		
		}
		
		}

		

	QD = logsum1K(QD,Kterm); // sum(exp(qd))=1 then add K 

			
	double out = logsum(QD);			

						

return(out);			
}

// [[Rcpp::export]]
NumericVector ppadjT5(const int N, const NumericVector& nummods,
					const List Vr, const List Cr, const double dcon,
					const List keep, const NumericMatrix& Nqq,
					const NumericVector& Ldcon12,
					const NumericVector& Nq3,
					const NumericVector& Ldcon123, const List Cr3,
					const NumericVector& Nq4, const NumericVector& Ldcon1234, 
					const List Cr4, const List lPP, const int Nsame) {		
					
	const int mT1 = nummods(0);
	const NumericVector& lPP1 = lPP[0];
	NumericVector QD(mT1);
//	NumericVector lppadj(mT1);
	
	for(int i=0; i<mT1; i++) {
		QD(i) = calcQD5(i, N,nummods,Vr, Cr, dcon,keep,  Nqq, Ldcon12, Nq3,
					Ldcon123,  Cr3, Nq4,Ldcon1234, Cr4, lPP, Nsame);
	}
	
	NumericVector lppadj = logsum1K(QD,lPP1);
	 	 
	 const NumericVector& lppadj1 = logsum1(lppadj);
	 const NumericVector& ppadj = exp(lppadj1);
	 					
return(ppadj);					
}					







// [[Rcpp::export]]
double calcQD6fast(const int mod1, const int N,
					const NumericVector& nummods,
					const List Vr, const List Cr, const double dcon,
					const List keep, const List lPP) {		
	
	const int T1 = 0;
	const int M = 6;
	
	IntegerVector imods(M), im2(M-1);
	
	
	const int T2 = 1;
	const int T3 = 2;
	const int T4 = 3;
	const int T5 = 4;
	const int T6 = 5;
	
	const int mT2 = nummods(T2);
	const int mT3 = nummods(T3);
	const int mT4 = nummods(T4);
	const int mT5 = nummods(T5);
	const int mT6 = nummods(T6);

	const int Nc = mT2*mT3*mT4*mT5*mT6;
	NumericVector QD(Nc);
	NumericVector Kterm(Nc);

	
	
	const NumericVector& lPP2 = lPP[1];
	const NumericVector& lPP3 = lPP[2];
	const NumericVector& lPP4 = lPP[3];
	const NumericVector& lPP5 = lPP[4];
	const NumericVector& lPP6 = lPP[5];

//	int Tt1, Tt2;
	
	const NumericMatrix& K1 = keep[0];
	const NumericMatrix& K2 = keep[1]; // keep is in order and on log scale
	const NumericMatrix& K3 = keep[2];
	const NumericMatrix& K4 = keep[3];
	const NumericMatrix& K5 = keep[4];
	

//	const int np = M*(M-1)/2; // number of pairs
	const NumericMatrix& c2 = pairs(M);		
	
	
		imods(T1) = mod1;
		int i = 0;  // each i is for a (j,k,l,r,s) models for T2,T3,T4,T5,T6
		
		for(int j = 0; j < mT2; j++) {
			imods(T2) = j;
			for(int k = 0; k < mT3; k++) {
				imods(T3) = k;	
				for(int l = 0; l < mT4; l++) {
					imods(T4) = l;
					for(int r = 0; r < mT5; r++) {	
					imods(T5) = r;
						for(int s = 0; s < mT6; s++) {	
						imods(T6) = s;				
						
				QD(i) = calcDelta(imods, N, M, dcon, Vr, Cr, c2) + lPP2(j) + lPP3(k) +lPP4(l) +lPP5(r) +lPP6(s);
				Kterm(i) = K1(mod1,j)+K2(mod1,k)+K3(mod1,l)+K4(mod1,r)+K5(mod1,s);
		
	}
			 
			i += 1;			
					
			}
			
			
		} 
		
		}
		
		}

		

	QD = logsum1K(QD,Kterm); // sum(exp(qd))=1 then add K 

			
	double out = logsum(QD);			

						

return(out);			
}

// [[Rcpp::export]]
NumericVector ppadjT6fast(const int N, const NumericVector& nummods,
					const List Vr, const List Cr, const double dcon,
					const List keep, const List lPP) {		
					
	const int mT1 = nummods(0);
	const NumericVector& lPP1 = lPP[0];
	NumericVector QD(mT1);
//	NumericVector lppadj(mT1);
	
	for(int i=0; i<mT1; i++) {
		QD(i) = calcQD6fast(i, N,nummods,Vr, Cr, dcon,keep, lPP);
	}
	
	NumericVector lppadj = logsum1K(QD,lPP1);
	 	 
	 const NumericVector& lppadj1 = logsum1(lppadj);
	 const NumericVector& ppadj = exp(lppadj1);
	 					
return(ppadj);					
}					

