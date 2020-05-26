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
NumericVector calcDij(int i, int j, 
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
NumericMatrix pairs(int M) {
// returns matrix of 2 rows and M*(M-1)/2 columns; each column is a pair of indices
	Function f("combn");
//	NumericVector inds = indices(M); 
	IntegerVector inds = seq(0,M-1);
	return f(_["x"]=inds,_["m"]=2);
}

// [[Rcpp::export]]
NumericVector indices(int M) {
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
double calcDelta(const IntegerVector& imods, // vector of model indices for M traits, in same order as traits
						int N, int M, double dcon,
						List Vr, List Cr,
						NumericMatrix c2) {
// N is the number measured for all M traits						
	int T1;
	int T2;
	int mod1;
	int mod2;
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
double calcQD3(int mod1, int N,
					NumericVector nummods,
					List Vr, List Cr,double dcon,
					List keep, NumericMatrix Nqq,
					NumericVector Ldcon12,
					NumericMatrix c2,
					List lPP,int Nsame) {		
	
	int T1 = 0;
	int M = 3;
	
	
	IntegerVector imods(M), im2(M-1);
	
	NumericVector tmp = indices(M);
//	NumericVector notT1 = tmp[ tmp != T1 ];
//	int T2 = notT1(0);
//	int T3 = notT1(1);
	int T2 = 1;
	int T3 = 2;
	int mT2 = nummods(T2);
	int mT3 = nummods(T3);

	int Nc = mT2*mT3;
	NumericVector QD(Nc);
	NumericVector Kterm(Nc);
	NumericVector lPP2 = lPP[1];
	NumericVector lPP3 = lPP[2];

	int Tt1, Tt2;
	
	NumericMatrix K1 = keep[0];
	NumericMatrix K2 = keep[1]; // keep is in order and on log scale
	
	double  ddd,nn;

	NumericVector dd(2);
	double val;
	double sign;
	
	mat Dij(M-1,M-1,fill::eye); // identity matrix size M
			
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
		int np = M*(M-1)/2; // number of pairs
			
		for(int l=0; l<np ; ++l) {
			Tt1 = c2(0,l);
			Tt2 = c2(1,l);
			int mod1 = imods(Tt1);
			int mod2 = imods(Tt2);
			int N12 = Nqq(Tt1,Tt2);	
			
			if(N12 > N) {
				nn=(N12-N)/N12; 
				NumericMatrix Cr12 = Cr[l]; // Cij matrices are listed in same order as pairs		  
				NumericVector V1 = Vr[Tt1];
				NumericVector V2 = Vr[Tt2];
				dd = calcDij(mod1, mod2, Cr12, V1, V2);
				Dij(0,1) = dd(0);
				Dij(1,0) = dd(1);
				log_det(val, sign, Dij);
				double lDD = val*sign;
				double dcon12 = Ldcon12(l);
				ddd = -N12*0.5*(lDD-dcon12);		
				QD(i) += nn*ddd;
				}
			}
		}		
				i += 1;				
			}
		}

		
		QD = logsum1(QD); // sum(exp(qd))=1
		
	for(i=0; i< Nc; i++) QD(i) += Kterm(i);
	
	double out = logsum(QD);			

return(out);			
}


// [[Rcpp::export]]
NumericVector ppadjT3(int N,
					NumericVector nummods,
					List Vr, List Cr,double dcon,
					List keep, NumericMatrix Nqq,
					NumericVector Ldcon12,
					NumericMatrix c2,
					List lPP, int Nsame) {	 	
					
	int mT1 = nummods(0);
	NumericVector lPP1 = lPP[0];
	NumericVector QD(mT1);
	NumericVector lppadj(mT1);
	
	for(int i=0; i<mT1; i++) {
		QD(i) = calcQD3(i, N,nummods,Vr, Cr, dcon,keep,  Nqq, Ldcon12, c2, lPP, Nsame);
	}

	
	double lsum = logsum(QD);

	for(int i = 0; i<mT1; i++) {
		lppadj(i) = QD(i)-lsum + lPP1(i);	
	 	}
	 	 
	 NumericVector lppadj1 = logsum1(lppadj);
	 NumericVector ppadj = exp(lppadj1);
	 					
return(ppadj);					
}					



// [[Rcpp::export]]
double calcQD4(int mod1, int N,
					NumericVector nummods,
					List Vr, List Cr,double dcon,
					List keep, NumericMatrix Nqq, 
					NumericVector Ldcon12,
					NumericVector Nq3,
					NumericVector Ldcon123, List Cr3,
					List lPP, int Nsame) {		
	
	int T1 = 0;
	int M = 4;
	
	
	IntegerVector imods(M), im2(M-1);
	
	NumericVector tmp = indices(M);
	int T2 = 1;
	int T3 = 2;
	int T4 = 3;
	int mT2 = nummods(T2);
	int mT3 = nummods(T3);
	int mT4 = nummods(T4);

	int Nc = mT2*mT3*mT4;
	NumericVector QD(Nc);
	NumericVector Kterm(Nc);
	NumericVector lPP2 = lPP[1];
	NumericVector lPP3 = lPP[2];
	NumericVector lPP4 = lPP[3];
	
	int Tt1, Tt2;
	
	NumericMatrix K1 = keep[0];
	NumericMatrix K2 = keep[1]; // keep is in order and on log scale
	NumericMatrix K3 = keep[2];
	
	double ddd,nn;

	NumericVector dd(2);
	double val;
	double sign;
	
	mat Dij(2,2,fill::eye); // identity matrix size M
	int np = M*(M-1)/2; // number of pairs
	NumericMatrix c2 = pairs(4);		

		// calculate all Dijkl*PPj*PPk*PPl (log sum) for given i and adjust to sum to 1
		imods(T1) = mod1;
		int i = 0;  // each i is for a triple (j,k,l) models for T2,T3,T4

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
			Tt1 = c2(0,l);
			Tt2 = c2(1,l);
			int mod1 = imods(Tt1);
			int mod2 = imods(Tt2);
			NumericVector tmp = indices(M);
			NumericVector notT12 = tmp[ (tmp != Tt1) & (tmp != Tt2)];
			int nTt1 = notT12(0);
			int nTt2 = notT12(1);
			int N12 = Nqq(Tt1,Tt2) - Nq3(nTt1) - Nq3(nTt2) + N ;	
			
			if(N12 > 0) {
				nn=N12/Nqq(Tt1,Tt2); 
				NumericMatrix Cr12 = Cr[l]; // Cij matrices are listed in same order as pairs		  
				NumericVector V1 = Vr[Tt1];
				NumericVector V2 = Vr[Tt2];
				dd = calcDij(mod1, mod2, Cr12, V1, V2);
				Dij(0,1) = dd(0);
				Dij(1,0) = dd(1);
				log_det(val, sign, Dij);
				double lDD = val*sign;
				double dcon12 = Ldcon12(l);
				ddd = -Nqq(Tt1,Tt2)*0.5*(lDD-dcon12);		
				QD(i) += nn*ddd;
				}
			}
	
	
	// calculate d123 adjustment term between all triples
		
		
		for(int l=0; l<M; ++l) {
			NumericVector tmp = indices(M);
			NumericVector notT1 = tmp[ tmp != l];
			int nTt1 = notT1(0);
			int nTt2 = notT1(1);
			int nTt3 = notT1(2);
	
			int N123 = Nq3(l) -  N ;	
			
			if(N123 > 0) {
				nn=N123/Nq3(l);
				List Vr3 = List(3);
				Vr3[0]=Vr[nTt1];
				Vr3[1]=Vr[nTt2];
				Vr3[2]=Vr[nTt3];
				NumericMatrix c23 = pairs(3);
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
				i += 1;				
			}
		}
		}

		QD = logsum1(QD); // sum(exp(qd))=1


	for(i=0; i< Nc; i++) QD(i) += Kterm(i);
	
	double out = logsum(QD);			

return(out);			
}


// [[Rcpp::export]]
NumericVector ppadjT4(int N,
					NumericVector nummods,
					List Vr, List Cr,double dcon,
					List keep, NumericMatrix Nqq, 
					NumericVector Ldcon12,
					NumericVector Nq3,
					NumericVector Ldcon123, List Cr3,
					List lPP, int Nsame) {	 	
					
	int mT1 = nummods(0);
	NumericVector lPP1 = lPP[0];
	NumericVector QD(mT1);
	NumericVector lppadj(mT1);
	
	for(int i=0; i<mT1; i++) {
		QD(i) = calcQD4(i, N,nummods,Vr, Cr, dcon,keep,  Nqq, Ldcon12, Nq3,
					Ldcon123,  Cr3, lPP, Nsame);
	}
	
	double lsum = logsum(QD);

	for(int i = 0; i<mT1; i++) {
		lppadj(i) = QD(i)-lsum + lPP1(i);	
	 	}
	 	 
	 NumericVector lppadj1 = logsum1(lppadj);
	 NumericVector ppadj = exp(lppadj1);
	 					
return(ppadj);					
}					



// [[Rcpp::export]]
double calcQD5(int mod1, int N,
					NumericVector nummods,
					List Vr, List Cr,double dcon,
					List keep, NumericMatrix Nqq,
					NumericVector Ldcon12,
					NumericVector Nq3,
					NumericVector Ldcon123, List Cr3,
					NumericVector Nq4, NumericVector Ldcon1234, 
					List Cr4,List lPP, int Nsame) {		
	
	int T1 = 0;
	int M = 5;
	
	IntegerVector imods(M), im2(M-1);
	
	NumericVector tmp = indices(M);
	int T2 = 1;
	int T3 = 2;
	int T4 = 3;
	int T5 = 4;
	int mT2 = nummods(T2);
	int mT3 = nummods(T3);
	int mT4 = nummods(T4);
	int mT5 = nummods(T5);

	int Nc = mT2*mT3*mT4*mT5;
	NumericVector QD(Nc);
	NumericVector Kterm(Nc);
	NumericVector lPP2 = lPP[1];
	NumericVector lPP3 = lPP[2];
	NumericVector lPP4 = lPP[3];
	NumericVector lPP5 = lPP[4];

	int Tt1, Tt2;
	
	NumericMatrix K1 = keep[0];
	NumericMatrix K2 = keep[1]; // keep is in order and on log scale
	NumericMatrix K3 = keep[2];
	NumericMatrix K4 = keep[3];
	
//	double delta, d12, ddd,nn, d123, d1234;
	double nn, ddd;

	NumericVector dd(2);
	double val;
	double sign;
	
	mat Dij(2,2,fill::eye); // identity matrix size M
	int np = M*(M-1)/2; // number of pairs
	NumericMatrix c2 = pairs(M);		

		imods(T1) = mod1;
		int i = 0;  // each i is for a quadruple (j,k,l,r) models for T2,T3,T4,T5
		
		for(int j = 0; j < mT2; j++) {
			imods(T2) = j;
			for(int k = 0; k < mT3; k++) {
				imods(T3) = k;	
				for(int l = 0; l < mT4; l++) {
					imods(T4) = l;
					for(int r = 0; r < mT5; r++) {	
					imods(T5) = r;		
					QD(i) = j;		
				QD(i) = calcDelta(imods, N, M, dcon, Vr, Cr, c2) + lPP2(j) + lPP3(k) +lPP4(l) +lPP5(r);
				Kterm(i) = K1(mod1,j)+K2(mod1,k)+K3(mod1,l)+K4(mod1,r);

	// calculate d12 adjustment term between all pairs	if have missing trait measurements		
			
		if(Nsame >0){
				
		for(int l=0; l<np ; ++l) {
			Tt1 = c2(0,l);
			Tt2 = c2(1,l);
			int mod1 = imods(Tt1);
			int mod2 = imods(Tt2);
			NumericVector tmp = indices(M);
			NumericVector notT12 = tmp[ (tmp != Tt1) & (tmp != Tt2)];
			int nTt1 = notT12(0);
			int nTt2 = notT12(1);
			int nTt3 = notT12(2);
			int N12 = Nqq(Tt1,Tt2) - Nq3(nTt1) - Nq3(nTt2) -Nq3(nTt3) + Nq4(nTt1) + Nq4(nTt2) + Nq4(nTt3)-N ;	
			
			if(N12 > 0) {
				nn=N12/Nqq(Tt1,Tt2); 
				NumericMatrix Cr12 = Cr[l]; // Cij matrices are listed in same order as pairs		  
				NumericVector V1 = Vr[Tt1];
				NumericVector V2 = Vr[Tt2];
				dd = calcDij(mod1, mod2, Cr12, V1, V2);
				Dij(0,1) = dd(0);
				Dij(1,0) = dd(1);
				log_det(val, sign, Dij);
				double lDD = val*sign;
				double dcon12 = Ldcon12(l);
				QD(i) += nn*(-Nqq(Tt1,Tt2)*0.5*(lDD-dcon12));		
//				d12 += nn*ddd;
				}
			}

	// calculate d123 adjustment term between all triples
		
		for(int l=0; l<np ; ++l) {
			Tt1 = c2(0,l);
			Tt2 = c2(1,l);
			NumericVector tmp = indices(M);
			NumericVector notT12 = tmp[ (tmp != Tt1) & (tmp != Tt2)]; // the 3 that are not in the pair are the focus
			int nTt1 = notT12(0);
			int nTt2 = notT12(1);
			int nTt3 = notT12(2);
			int N123 = Nq3(l) - Nq4(Tt1) - Nq4(Tt2) + N ;	
			
			if(N123 > 0) {
				nn=N123/Nq3(l); 
				List Vr3 = List(3);
				Vr3[0]=Vr[nTt1];
				Vr3[1]=Vr[nTt2];
				Vr3[2]=Vr[nTt3];
				NumericMatrix c234 = pairs(3);
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
				nn=N1234/Nq4(l);
				List Vr4 = List(4);
				Vr4[0]=Vr[nTt1];
				Vr4[1]=Vr[nTt2];
				Vr4[2]=Vr[nTt3];
				Vr4[3]=Vr[nTt4];
				NumericMatrix c234 = pairs(4);
				IntegerVector imods4(4);
				imods4(0)=imods(nTt1);
				imods4(1)=imods(nTt2);
				imods4(2)=imods(nTt3);
				imods4(3)=imods(nTt4);
				List Cr4l = Cr4[l];
				double dcon1234 = Ldcon1234(l);
				int nq1234 = Nq4(l);
				ddd=calcDelta(imods4,nq1234 , 4, dcon1234, Vr4, Cr4l,c234);
			 	QD(i) += nn*ddd;
			}			
		}
		
	}
			 
				i += 1;		
					
			}
			
		} 
		
		}
		
		}

		QD = logsum1(QD); // sum(exp(qd))=1


	for(i=0; i< Nc; i++) QD(i) += Kterm(i);
	
	double out = logsum(QD);	
						

return(out);			
}

// [[Rcpp::export]]
NumericVector ppadjT5(int N, NumericVector nummods,
					List Vr, List Cr,double dcon,
					List keep, NumericMatrix Nqq,
					NumericVector Ldcon12,
					NumericVector Nq3,
					NumericVector Ldcon123, List Cr3,
					NumericVector Nq4, NumericVector Ldcon1234, 
					List Cr4,List lPP, int Nsame) {		
					
	int mT1 = nummods(0);
	NumericVector lPP1 = lPP[0];
	NumericVector QD(mT1);
	NumericVector lppadj(mT1);
	
	for(int i=0; i<mT1; i++) {
		QD(i) = calcQD5(i, N,nummods,Vr, Cr, dcon,keep,  Nqq, Ldcon12, Nq3,
					Ldcon123,  Cr3, Nq4,Ldcon1234, Cr4, lPP, Nsame);
	}
	
	double lsum = logsum(QD);
//	NumericVector qp(mT1); 
	for(int i = 0; i<mT1; i++) {
		lppadj(i) = QD(i)-lsum + lPP1(i);	
//	 	lppadj(i) += lPP1(i) + qp(i);
	 	}
	 	 
	 NumericVector lppadj1 = logsum1(lppadj);
	 NumericVector ppadj = exp(lppadj1);
	 					
return(ppadj);					
}					
