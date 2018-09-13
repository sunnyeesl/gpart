#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include <string>     // to use std::string
#include <vector>
#include <math.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;  // to avoid typing std::



// [[Rcpp::export]]
IntegerMatrix VCFtogeno(CharacterMatrix vcf){
  int nSNP = vcf.nrow();
  int nSample = vcf.ncol();

  vector<int> genotypes;
  for(int i=0; i < nSNP; i++) {
    for(int j = 0; j<nSample; j++){
      if ( vcf(i,j)[0] == '.' ) genotypes.push_back(-9);
      else  // assume that allele is 0 to 9
        genotypes.push_back((vcf(i,j)[0]- '0') + (vcf(i,j)[2] - '0'));
    }
  }
  IntegerMatrix mat(nSample, nSNP);
  for(int i=0, k=0; i < nSNP; ++i) {
    for(int j=0; j < nSample; ++j, ++k) {
      if ( genotypes[k] < 0 )
        mat(j,i) = NumericVector::get_na();
      else
        mat(j,i) = genotypes[k];
    }
  }
  return mat;
}


// [[Rcpp::export]]
NumericVector pairCubeX(IntegerVector b1, IntegerVector  b2) {
  NumericMatrix counts012(3,3);
  //
  int n1111, n1112, n1122, n1211, n1212, n1222, n2211, n2212, n2222;
  int n;
  int n11;
  // int n12, n21, n22;
  double p, q, a, b, c, dee;
  double minhap, maxhap;
  double xN, d2, yN, yN2, h2;
  double alpha = -1, beta = -1, gamma = -1;
  NumericVector roots(3);
  //
  int b1size = b1.size();
  for (int i = 0; i < b1size; i++) {
    if(b1[i]!=NA_INTEGER && b2[i]!=NA_INTEGER){
      counts012(b1[i],b2[i])++;
    }
  }
  // return counts012;
  n1111 = counts012(0,0);
  n1112 = counts012(0,1);
  n1122 = counts012(0,2);
  n1211 = counts012(1,0);
  n1212 = counts012(1,1);
  n1222 = counts012(1,2);
  n2211 = counts012(2,0);
  n2212 = counts012(2,1);
  n2222 = counts012(2,2);
  n = (n1111 + n1112 + n1122 + n1211 + n1212 + n1222 + n2211 + n2212 + n2222);
  p = ((n1111+n1112+n1122)*2.0+(n1211+n1212+n1222))/(2.0 * n);
  q = ((n1111+n1211+n2211)*2.0+(n1112+n1212+n2212))/(2.0 * n);
  n11 = (2.0*n1111 + n1112 + n1211) ;
  // n12 = (2.0*n1122 + n1112 + n1222);
  // n21 = (2.0*n2211 + n2212 + n1211);
  // n22 = (2.0*n2222 + n2212 + n1222) ;
  dee = -n11*p*q;
  c = -n11*(1.0 - 2.0*p - 2.0*q) - n1212*(1.0 - p - q) + 2.0*n*p*q;
  b = 2.0*n*(1.0 - 2.0*p - 2.0*q) - 2.0*n11 - n1212;
  a = 4.0 * n;

  minhap = n11 / (2.0 * n);
  maxhap = (n11 + n1212) / (2.0 * n);

  xN = -b/(3.0*a);
  d2 = (pow(b,2)-3.0*a*c)/(9*pow(a,2));
  yN = a * pow(xN,3) + b * pow(xN,2) + c * xN + dee;
  yN2 = pow(yN,2);
  h2 = 4 * pow(a,2) * pow(d2,3);

  if (yN2 > h2){
    double number1, number2;
    number1 = 0.0;
    number2 = 0.0;
    if ((1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))) < 0){
      number1 = -pow(-(1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))),1.0/3.0);
    } else {
      number1 = pow((1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))),1.0/3.0);
    }
    if ((1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))) < 0){
      number2 = -pow(-(1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))),1.0/3.0);
    }else {
      number2 = pow((1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))),1.0/3.0);
    }
    alpha = xN + number1 + number2;
    beta = -1;
    gamma = -1;
  } else if (yN2 == h2){
    double delta;
    delta = pow((yN/2.0*a),(1.0/3.0));
    alpha = xN + delta;
    beta = xN + delta;
    gamma = xN - 2.0*delta;
  } else if (yN2 < h2){
    double h, theta, delta;
    h = pow(h2, 0.5);
    theta = ((acos(-yN/h))/3.0);
    delta = pow(d2,0.5);
    alpha = xN + 2.0 * delta * cos(theta);
    beta = xN + 2.0 * delta * cos(2.0 * M_PI/3.0 + theta);
    gamma = xN + 2.0 * delta * cos(4.0 * M_PI/3.0 + theta);
  }

  roots[0] = alpha;
  roots[1] = beta;
  roots[2] = gamma;

  NumericVector total_loglik(3);
  NumericVector total_Dprime(3);
  NumericVector total_D(3);
  double f11, f12, f21, f22;
  double D, D_max, Dp;
  for (int i = 0; i<3; i++){
    if(roots[i]>=minhap-0.00001 && roots[i]<=maxhap+0.00001){
      f11 = roots[i];
      f12 = p - f11;
      f21 = q - f11;
      f22 = 1 - (f11+f12+f21);
      if(f12 <= 0) f12 = abs(f12)+1e-7;
      if(f22 <= 0) f22 = abs(f22)+1e-7;
      if(f21 <= 0) f21 = abs(f21)+1e-7;
      D = f11*f22-f12*f21;
      total_D[i] = D;
      if(D>0){
        D_max =  min(NumericVector::create((1-p)*q, p*(1-q)));
      }else{
        D_max =  max(NumericVector::create(-p*q,-(1-p)*(1-q)));
      }
      Dp = D/D_max;
      total_Dprime[i] =  Dp;
      double loglik;
      loglik = ((2*n1111+n1112+n1211)*log10(f11)+(2*n1122 + n1222 + n1112)*log10(f12)+
        (2*n2211 + n1211 + n2212)*log10(f21)+(2*n2222 + n2212 + n1222)*log10(f22) +
        (n1212)*log10(f11*f22+f12*f21));
      total_loglik[i] = loglik;
    }else{
      total_Dprime[i] =  NA_REAL;
      total_loglik[i] = NA_REAL;
    }
  }

  // return roots;
  // return total_Dprime;

  NumericVector real_loglik(3);
  NumericVector real_Dp(3);
  NumericVector real_D(3);
  int noNAnum =0;
  for(int i=0;i<3;i++){
    if(total_loglik[i] != NA_REAL){
      real_loglik[noNAnum] = total_loglik[i];
      real_Dp[noNAnum] = total_Dprime[i];
      real_D[noNAnum] = total_D[i];
      noNAnum++;
    }
  }
  double maxlog, fDp, fD;
  if(noNAnum>0){
    maxlog = real_loglik[0];
    fDp = real_Dp[0];
    fD = real_D[0];
    if(noNAnum>1){
      for(int i=1;i<noNAnum;i++){
        double Nmaxlog, nDp, nD;
        Nmaxlog = real_loglik[i];
        nDp = real_Dp[i];
        nD = real_D[i];
        if(maxlog < Nmaxlog){
          maxlog = Nmaxlog;
          fDp = nDp;
          fD = nD;
        }
      }
    }
  }else{
    maxlog = NA_REAL;
    fDp = NA_REAL;
    fD = NA_REAL;
  }
  NumericVector res(3);
  res[0] = maxlog;
  res[1] = fDp;
  res[2] = fD;
  return res;
}


NumericVector pairCubeX(IntegerVector b1, IntegerVector  b2);
// [[Rcpp::export]]
NumericMatrix matCubeX(IntegerMatrix geno){
  int ncol = geno.ncol();
  NumericMatrix resMat(ncol, ncol);
  IntegerVector b1, b2;
  NumericVector Dpres;
  double Dp;
  for(int i=0; i<ncol; i++){
    for(int j=0; j<ncol; j++){
      if(i<j){
        b1 = geno(_,i);
        b2 = geno(_,j);
        Dpres = pairCubeX(b1, b2);
        Dp = Dpres[1];
        resMat(i,j) = Dp;
        resMat(j,i) = Dp;
      }
    }
  }
  return resMat;
}



NumericVector pairCubeX(IntegerVector b1, IntegerVector  b2);
// [[Rcpp::export]]
NumericMatrix matCubeX2(IntegerMatrix geno1, IntegerMatrix geno2){
  int ncol1 = geno1.ncol();
  int ncol2 = geno2.ncol();
  NumericMatrix resMat(ncol1, ncol2);
  IntegerVector b1, b2;
  NumericVector Dpres;
  double Dp;

  for(int i=0; i<ncol1; i++){
    for(int j=0; j<ncol2; j++){
      b1 = geno1(_,i);
      b2 = geno2(_,j);
      Dpres = pairCubeX(b1, b2);
      Dp = Dpres[1];
      resMat(i,j) = Dp;
    }
  }
  return resMat;
}

// [[Rcpp::export]]
int estiD(double pA, double pB, int n11,int n12,int n21,int n22, int n1212){

  double pa,pb,eDmin,eDmax;
  double tf11, tf12, tf21, tf22;
  int minD, maxD, pmin, pmax;
  pa = 1-pA;
  pb = 1-pB;
  eDmin = max(NumericVector::create(-pA*pB, -pa*pb));
  eDmax = min(NumericVector::create(pA * pb, pB * pa));
  pmin = pA*pB +eDmin;
  pmax = pA*pB +eDmax;
  // printf("%f, %f \n", eDmin, eDmax);
  double Dnowlik, Dmaxlog= NA_REAL;
  int Dmaxindex = NA_INTEGER;
  bool start = true;
  minD = (int) floor((pmin*1000));
  maxD = (int) floor((pmax*1000));
  // printf("%d %d \n", minD, maxD);
  for(int k = minD; k <= maxD;k++){
    tf11 = k*0.001+pA*pB;
    tf12 = pA - tf11;
    tf21 = pB - tf11;
    tf22 = 1 - (tf11+tf12+tf21);
    // printf("%f %f %f %f \n", tf11, tf12, tf21, tf22);
    if(tf11>0 && tf12>0 && tf21>0 && tf22>0){
      // printf("%d ",k);
      Dnowlik = (n11*log10(tf11)+n12*log10(tf12)+n21*log10(tf21)+n22*log10(tf22) +(n1212)*log10(tf11*tf22+tf12*tf21));
      // printf("%f\n ",Dnowlik);
      if(start == true){
        // printf("start");
        Dmaxlog = Dnowlik;
        Dmaxindex = k;
        // printf("%f \n", Dmaxlog);
        // printf("%f \n", Dnowlik);
        // printf("%d\n",k);
        // printf("%d \n", Dmaxindex);
        start = false;
      }else{
        if(Dnowlik > Dmaxlog){
          Dmaxlog = Dnowlik;
          Dmaxindex = k;
          // printf("Dmaxlog %f \n", Dmaxlog);
          // printf("Dnowlik %f \n", Dnowlik);
          // printf("%d\n",k);
          // printf("%d \n", Dmaxindex);
        }
      }

    }
  }
  return Dmaxindex;
}



NumericVector pairCubeX(IntegerVector b1, IntegerVector  b2);
int estiD(double pA, double pB, int n11,int n12,int n21,int n22, int n1212);
// [[Rcpp::export]]
NumericVector CIDp(IntegerVector b1, IntegerVector b2){
  NumericVector cubeXres(3);
  double eD, eDp, Dp;
  NumericMatrix counts012(3,3);
  double lowDp = NA_REAL, upperDp = NA_REAL;
  //
  int n1111, n1112, n1122, n1211, n1212, n1222, n2211, n2212, n2222;
  int n;
  int n11;
  int n12, n21, n22;
  double p, q;
  // bool Dpyes = false;

  int b1size = b1.size();
  for (int i = 0; i < b1size; i++) {
    if(b1[i]!=NA_INTEGER && b2[i]!=NA_INTEGER){
      counts012(b1[i],b2[i])++;
    }
  }
  // return counts012;
  n1111 = counts012(0,0);
  n1112 = counts012(0,1);
  n1122 = counts012(0,2);
  n1211 = counts012(1,0);
  n1212 = counts012(1,1);
  n1222 = counts012(1,2);
  n2211 = counts012(2,0);
  n2212 = counts012(2,1);
  n2222 = counts012(2,2);
  n = (n1111 + n1112 + n1122 + n1211 + n1212 + n1222 + n2211 + n2212 + n2222);
  p = ((n1111+n1112+n1122)*2.0+(n1211+n1212+n1222))/(2.0 * n);
  q = ((n1111+n1211+n2211)*2.0+(n1112+n1212+n2212))/(2.0 * n);
  n11 = (2.0*n1111 + n1112 + n1211) ;
  n12 = (2.0*n1122 + n1112 + n1222);
  n21 = (2.0*n2211 + n2212 + n1211);
  n22 = (2.0*n2222 + n2212 + n1222);

  double f11, f12, f21, f22;
  cubeXres = pairCubeX(b1,b2);
  Dp = cubeXres[1];
  if(cubeXres[1] == NA){ //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    // estimate D
    int Dmaxindex = estiD(p, q, n11, n12, n21, n22, n1212);
    // printf("Dmaxindex %d", Dmaxindex);
    // Dp = eDp;
    eD = Dmaxindex*0.001;
  }else{
    // Dp = eDp;
    eD = cubeXres[2];
  }
  // printf("eD : %f", eD);
  double Dmax=NA_REAL;
  // printf("eD : %f", eD);
  // printf("p %f \n", p);
  // printf("q %f \n", q);
  if(eD > 0){
    Dmax =  min(NumericVector::create((1-p)*q, p*(1-q)));
    // printf("1");
  }else if(eD ==0){
    // printf("NONONONO \n");
  }else if(eD <0){
    Dmax = max(NumericVector::create(-p*q, -(1-p)*(1-q)));
    // printf("3");
  }

  // if(Dmax==0) printf("Dmax is zero\n");
  // printf("Dmax : %f", Dmax);
  if(Dp != NA_REAL && Dp < 0.9){
    double cutDp;
    for(int i=1; i<91; i++){
      cutDp = 0.01*i;
      if(Dp >= cutDp-0.01 && Dp < cutDp) {
        break;
      }
    }
    // t1 = cutDp;
    double loglikDp, loglik90;
    // printf("Dmax : %f \n", Dmax);
    f11 = cutDp * Dmax +p*q;
    f12 = p - f11;
    f21 = q - f11;
    f22 = 1 - (f11+f12+f21);

    loglikDp = (n11*log10(f11)+n12*log10(f12)+n21*log10(f21)+n22*log10(f22) +(n1212)*log10(f11*f22+f12*f21));
    // printf("loglikDp : %.19f \n", loglikDp);
    f11 = 0.9 * Dmax +p*q;
    f12 = p - f11;
    f21 = q - f11;
    f22 = 1 - (f11+f12+f21);
    loglik90 = (n11*log10(f11)+n12*log10(f12)+n21*log10(f21)+n22*log10(f22) +(n1212)*log10(f11*f22+f12*f21));
    // printf("loglik90 : %f\n", loglik90);
    if(pow(10.0, loglik90) < pow(10.0,loglikDp)/220){
      // printf("delete pair\n");
      lowDp = -3;
      upperDp = -3;

    }
  }
  // printf("lowDp %f \n", lowDp);
  if(lowDp != -3){
    // printf("true");
    NumericVector logliks(101);
    NumericVector sumlogliks(101);
    double maxlog=0, nowlik, totalsum;
    int maxindex=-1;
    for(int k=0; k<101;k++){
      // printf("true");
      f11 = 0.01*k*Dmax +p*q;
      f12 = p - f11;
      f21 = q - f11;
      f22 = 1 - (f11+f12+f21);
      // printf("Dmax %f \n", Dmax);
      // printf("F11 %f \n", f11);
      // printf("f11 %.19f f12 %.19f f21 %.19f f22 %.19f \n", f11, f12, f21, f22);
      if(f11 <= 0) {f11 = 0.0000000001;}
      if(f12 <= 0) {f12 = 0.0000000001;}
      if(f21 <= 0) {f21 = 0.0000000001;}
      if(f22 <= 0) {f22 = 0.0000000001;}
      nowlik = (n11*log10(f11)+n12*log10(f12)+n21*log10(f21)+n22*log10(f22) +(n1212)*log10(f11*f22+f12*f21));
      // printf("nowlik : %.19f \n", nowlik);
      logliks[k] =nowlik;
      if(k ==0){
        maxlog = nowlik;
        maxindex = k;
        sumlogliks[k] = pow(10.0, nowlik);
      }else{
        if(nowlik>maxlog){
          maxlog = nowlik;
          maxindex = k;
          // printf("%f %d \n", nowlik,maxindex);
        }
        sumlogliks[k] = sumlogliks[(k-1)] + pow(10.0, nowlik);
      }
    }

    eDp = 0.01*maxindex;
    // printf("maxindex %d \n", maxindex);
    // printf("eDp %f \n", eDp);
    // printf("Dp %f \n", Dp);
    // Rf_PrintValue(logliks);
    // printf("max index %d \n", maxindex);
    if(cubeXres[1] == NA){ //--------------------------------------------------------------------
      Dp = eDp;
    }
    totalsum = sumlogliks[100];
    for(int k = 0; k < 101; k++){
      if((sumlogliks[k]/totalsum) <= 0.05){
        lowDp = 0.01*(k+1);
      }else if((sumlogliks[k]/totalsum) > 0.95){
        upperDp = 0.01*k;
        break;
      }
    }
  }
  NumericVector Dpres(3);
  Dp = floor((Dp*1000)+0.5f)/1000;
  Dpres[0] = Dp;
  // Dpres[1] = lowDp;
  // Dpres[2] = upperDp;
  Dpres[1] = lowDp;
  Dpres[2] = upperDp;
  return Dpres;
}


NumericVector CIDp(IntegerVector b1, IntegerVector b2);
// [[Rcpp::export]]
NumericMatrix genoDp(IntegerMatrix geno, bool strLD = true, double lower = 0.7, double upper = 0.98){
  int ncol = geno.ncol();
  NumericMatrix resMat(ncol, ncol);
  IntegerVector b1, b2;
  NumericVector Dpres;
  double Dp, lowerb, upperb;

  for(int i=0; i<ncol; i++){
    for(int j=0; j<ncol; j++){
      if(i<j){
        b1 = geno(_,i);
        b2 = geno(_,j);
        Dpres = CIDp(b1, b2);
        Dp = Dpres[0];
        lowerb = Dpres[1];
        upperb = Dpres[2];
        if(strLD == true){
          if(lowerb>= lower && upperb>= upper){
            resMat(i,j) = 1;
            resMat(j,i) = 1;
          }else{
            resMat(i,j) = 0;
            resMat(j,i) = 0;
          }
        }else{
          resMat(i,j) = Dp;
          resMat(j,i) = Dp;
        }
      }
    }
  }
  return resMat;
}


NumericVector CIDp(IntegerVector b1, IntegerVector b2);
// [[Rcpp::export]]
NumericMatrix genoDp2(IntegerMatrix geno1, IntegerMatrix geno2, bool strLD = true, double lower = 0.7, double upper = 0.98){
  int ncol1 = geno1.ncol();
  int ncol2 = geno2.ncol();
  NumericMatrix resMat(ncol1, ncol2);
  IntegerVector b1, b2;
  NumericVector Dpres;
  double Dp, lowerb, upperb;

  for(int i=0; i<ncol1; i++){
    for(int j=0; j<ncol2; j++){
      b1 = geno1(_,i);
      b2 = geno2(_,j);
      Dpres = CIDp(b1, b2);
      Dp = Dpres[0];
      lowerb = Dpres[1];
      upperb = Dpres[2];
      if(strLD == true){
        if(lowerb>= lower && upperb>= upper){
          resMat(i,j) = 1;
          // resMat(j,i) = 1;
        }else{
          resMat(i,j) = 0;
          // resMat(j,i) = 0;
        }
      }else{
        resMat(i,j) = Dp;
        // resMat(j,i) = Dp;
      }
    }
  }
  return resMat;
}



NumericVector pairCubeX(IntegerVector b1, IntegerVector  b2);
// [[Rcpp::export]]
NumericMatrix genoCubeDp(IntegerMatrix geno){
  int ncol = geno.ncol();
  NumericMatrix resMat(ncol, ncol);
  IntegerVector b1, b2;
  NumericVector Dpres;
  double Dp;

  for(int i=0; i<(ncol-1); i++){
    for(int j=i; j<ncol; j++){
      b1 = geno(_,i);
      b2 = geno(_,j);
      Dpres = pairCubeX(b1, b2);
      Dp = Dpres[1];
      resMat(i,j) = Dp;
      resMat(j,i) = Dp;
    }
  }
  return resMat;
}
