#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include <string>     // to use std::string
#include <vector>
#include <math.h>
#include <algorithm>
#include <R.h>

using namespace Rcpp;
using namespace std;  // to avoid typing std::

//
#define EPSILON 0.000000000931322574615478515625
#define MINV(aa, bb) (((aa) > (bb))? (bb) : (aa))
#define RECIP_2_53 0.00000000000000011102230246251565404236316680908203125
#define SMALLISH_EPSILON 0.00000000002910383045673370361328125

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
  eDmin = max(NumericVector::create(-pA * pB, -pa * pb));
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
uint32_t CIDp_strLD(IntegerVector b1, IntegerVector b2, double lower, double upper);
// [[Rcpp::export]]
NumericMatrix genoDp(IntegerMatrix geno, bool strLD = true, double lower = 0.7, double upper = 0.98){
  int ncol = geno.ncol();
  NumericMatrix resMat(ncol, ncol);
  IntegerVector b1, b2;
  NumericVector Dpres;
  uint32_t Dpclass;
  double Dp;

  for(int i=0; i<ncol; i++){
    for(int j=0; j<ncol; j++){
      if(i<j){
        b1 = geno(_,i);
        b2 = geno(_,j);
        if(strLD == true){
          Dpclass = CIDp_strLD(b1, b2,lower = lower, upper=upper);
          if(Dpclass >= 4){
            resMat(i,j) = 1;
            resMat(j,i) = 1;
          }else{
            resMat(i,j) = 0;
            resMat(j,i) = 0;
          }
        }else{
          Dpres = CIDp(b1, b2);
          Dp = Dpres[0];
          resMat(i,j) = Dp;
          resMat(j,i) = Dp;
        }
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

NumericVector CIDp(IntegerVector b1, IntegerVector b2);
uint32_t CIDp_strLD(IntegerVector b1, IntegerVector b2, double lower, double upper);
// [[Rcpp::export]]
NumericMatrix genoDp2(IntegerMatrix geno1, IntegerMatrix geno2, bool strLD = true, double lower = 0.7, double upper = 0.98){
  int ncol1 = geno1.ncol();
  int ncol2 = geno2.ncol();
  NumericMatrix resMat(ncol1, ncol2);
  IntegerVector b1, b2;
  NumericVector Dpres;
  double Dp;
  uint32_t Dpclass;
  for(int i=0; i<ncol1; i++){
    for(int j=0; j<ncol2; j++){
      b1 = geno1(_,i);
      b2 = geno2(_,j);
      if(strLD == true){
        Dpclass = CIDp_strLD(b1, b2,lower = lower, upper=upper);
        if(Dpclass >= 4){
          resMat(i,j) = 1;
        }else{
          resMat(i,j) = 0;
        }
      }else{
        Dpres = CIDp(b1, b2);
        Dp = Dpres[0];
        resMat(i,j) = Dp;
      }
    }
  }
  return resMat;
}


// [[Rcpp::export]]
double calc_lnlike(double known11, double known12, double known21, double known22, double center_ct_d, double freq11, double freq12, double freq21, double freq22, double half_hethet_share, double freq11_incr) {
  double lnlike;
  freq11 += freq11_incr;
  freq22 += freq11_incr;
  freq12 += half_hethet_share - freq11_incr;
  freq21 += half_hethet_share - freq11_incr;
  lnlike = center_ct_d * log(freq11 * freq22 + freq12 * freq21);
  if (known11 != 0.0) {
    lnlike += known11 * log(freq11);
  }
  if (known12 != 0.0) {
    lnlike += known12 * log(freq12);
  }
  if (known21 != 0.0) {
    lnlike += known21 * log(freq21);
  }
  if (known22 != 0.0) {
    lnlike += known22 * log(freq22);
  }
  return lnlike;
}


// [[Rcpp::export]]
uint32_t cubic_real_roots(double coef_a, double coef_b, double coef_c, NumericVector solutions) {
  // Analytically finds all real roots of x^3 + ax^2 + bx + c, saving them in
  // solutions[] (sorted from smallest to largest), and returning the count.
  // Multiple roots are only returned/counted once.
  // Additional research into numerical stability may be in order here.
  double a2 = coef_a * coef_a;
  double qq = (a2 - 3 * coef_b) * (1.0 / 9.0);
  double rr = (2 * a2 * coef_a - 9 * coef_a * coef_b + 27 * coef_c) * (1.0 / 54.0);
  double r2 = rr * rr;
  double q3 = qq * qq * qq;
  double adiv3 = coef_a * (1.0 / 3.0);
  double sq;
  double dxx;
  if (r2 < q3) {
    // three real roots
    sq = sqrt(qq);
    dxx = acos(rr / (qq * sq)) * (1.0 / 3.0);
    sq *= -2;
    solutions[0] = sq * cos(dxx) - adiv3;
    solutions[1] = sq * cos(dxx + (2.0 * PI / 3.0)) - adiv3;
    solutions[2] = sq * cos(dxx - (2.0 * PI / 3.0)) - adiv3;
    // now sort and check for within-epsilon equality
    if (solutions[0] > solutions[1]) {
      dxx = solutions[0];
      solutions[0] = solutions[1];
      if (dxx > solutions[2]) {
        solutions[1] = solutions[2];
        solutions[2] = dxx;
      } else {
        solutions[1] = dxx;
      }
      if (solutions[0] > solutions[1]) {
        dxx = solutions[0];
        solutions[0] = solutions[1];
        solutions[1] = dxx;
      }
    } else if (solutions[1] > solutions[2]) {
      dxx = solutions[1];
      solutions[1] = solutions[2];
      solutions[2] = dxx;
    }
    if (solutions[1] - solutions[0] < EPSILON) {
      solutions[1] = solutions[2];
      return (solutions[1] - solutions[0] < EPSILON)? 1 : 2;
    }
    return (solutions[2] - solutions[1] < EPSILON)? 2 : 3;
  }
  dxx = -pow(fabs(rr) + sqrt(r2 - q3), 1.0 / 3.0);
  if (dxx == 0.0) {
    solutions[0] = -adiv3;
    return 1;
  }
  if (rr < 0.0) {
    dxx = -dxx;
  }
  sq = qq / dxx;
  solutions[0] = dxx + sq - adiv3;
  // use of regular epsilon here has actually burned us
  if (fabs(dxx - sq) >= (EPSILON * 8)) {
    return 1;
  }
  if (dxx >= 0.0) {
    solutions[1] = solutions[0];
    solutions[0] = -dxx - adiv3;
  } else {
    solutions[1] = -dxx - adiv3;
  }
  return 2;
}


// [[Rcpp::export]]
double calc_lnlike_quantile(double known11, double known12, double known21, double known22, double unknown_dh, double freqx1, double freq1x, double freq2x, double freq11_expected, double denom, int32_t quantile) {
  // almost identical to calc_lnlike, but we can skip the equal-to-zero checks
  // when quantile isn't 100
  double tmp11 = quantile * denom + freq11_expected;
  double tmp12 = freq1x - tmp11;
  double tmp21 = freqx1 - tmp11;
  double tmp22 = freq2x - tmp21;
  if (quantile == 100) {
    // One of these values will be ~zero, and we want to ensure its logarithm
    // is treated as a very negative number instead of nan.  May as well do it
    // the same way as Haploview.
    if (tmp11 < 1e-10) {
      tmp11 = 1e-10;
    }
    if (tmp12 < 1e-10) {
      tmp12 = 1e-10;
    }
    if (tmp21 < 1e-10) {
      tmp21 = 1e-10;
    }
    if (tmp22 < 1e-10) {
      tmp22 = 1e-10;
    }
  }
  return known11 * log(tmp11) + known12 * log(tmp12) + known21 * log(tmp21) + known22 * log(tmp22) + unknown_dh * log(tmp11 * tmp22 + tmp12 * tmp21);
}


uint32_t cubic_real_roots(double coef_a, double coef_b, double coef_c, NumericVector solutions);
double calc_lnlike(double known11, double known12, double known21, double known22, double center_ct_d, double freq11, double freq12, double freq21, double freq22, double half_hethet_share, double freq11_incr);
// [[Rcpp::export]]
NumericVector em_phase_hethet(double known11, double known12, double known21, double known22, uint32_t center_ct, uint32_t onside_sol_ct_ptr) {
  //, double freq1x_ptr, double freq2x_ptr, double freqx1_ptr, double freqx2_ptr, double freq11_ptr, uint32_t onside_sol_ct_ptr

  // Returns 1 if at least one SNP is monomorphic over all valid observations;
  // returns 0 otherwise, and fills all frequencies using the maximum
  // likelihood solution to the cubic equation.
  // (We're discontinuing most use of EM phasing since better algorithms have
  // been developed, but the two marker case is mathematically clean and fast
  // enough that it'll probably remain useful as an input for some of those
  // better algorithms...)
  double center_ct_d = (int32_t)center_ct;
  double twice_tot = known11 + known12 + known21 + known22 + 2 * center_ct_d;
  uint32_t sol_start_idx = 0;
  uint32_t sol_end_idx = 1;
  NumericVector solutions(3);
  double twice_tot_recip;
  double half_hethet_share;
  double freq11;
  double freq12;
  double freq21;
  double freq22;
  double prod_1122;
  double prod_1221;
  double incr_1122;
  double best_sol;
  double best_lnlike;
  double cur_lnlike;
  double freq1x;
  double freq2x;
  double freqx1;
  double freqx2;
  double lbound;
  double dxx;
  uint32_t cur_sol_idx;
  NumericVector final_sol(7); //adding
  // shouldn't have to worry about subtractive cancellation problems here
  if (twice_tot == 0.0) {
    return 1;
  }
  twice_tot_recip = 1.0 / twice_tot;
  freq11 = known11 * twice_tot_recip;
  freq12 = known12 * twice_tot_recip;
  freq21 = known21 * twice_tot_recip;
  freq22 = known22 * twice_tot_recip;
  prod_1122 = freq11 * freq22;
  prod_1221 = freq12 * freq21;
  half_hethet_share = center_ct_d * twice_tot_recip;
  // the following four values should all be guaranteed nonzero except in the
  // NAN case
  freq1x = freq11 + freq12 + half_hethet_share;
  freq2x = 1.0 - freq1x;
  freqx1 = freq11 + freq21 + half_hethet_share;
  freqx2 = 1.0 - freqx1;
  if (center_ct) {
    if ((prod_1122 != 0.0) || (prod_1221 != 0.0)) {
      sol_end_idx = cubic_real_roots(0.5 * (freq11 + freq22 - freq12 - freq21 - 3 * half_hethet_share), 0.5 * (prod_1122 + prod_1221 + half_hethet_share * (freq12 + freq21 - freq11 - freq22 + half_hethet_share)), -0.5 * half_hethet_share * prod_1122, solutions);
      while (sol_end_idx && (solutions[sol_end_idx - 1] > half_hethet_share + SMALLISH_EPSILON)) {
        sol_end_idx--;
      }
      while ((sol_start_idx < sol_end_idx) && (solutions[sol_start_idx] < -SMALLISH_EPSILON)) {
        sol_start_idx++;
      }
      if (sol_start_idx == sol_end_idx) {
        // Lost a planet Master Obi-Wan has.  How embarrassing...
        // lost root must be a double root at one of the boundary points, just
        // check their likelihoods
        sol_start_idx = 0;
        sol_end_idx = 2;
        solutions[0] = 0;
        solutions[1] = half_hethet_share;
      } else {
        if (solutions[sol_start_idx] < 0) {
          solutions[sol_start_idx] = 0;
        }
        if (solutions[sol_end_idx - 1] > half_hethet_share) {
          solutions[sol_end_idx - 1] = half_hethet_share;
        }
      }
    } else {
      solutions[0] = 0;
      // bugfix (6 Oct 2017): need to use all nonzero values here
      const double nonzero_freq_xx = freq11 + freq22;
      const double nonzero_freq_xy = freq12 + freq21;
      if ((nonzero_freq_xx + SMALLISH_EPSILON < half_hethet_share + nonzero_freq_xy) && (nonzero_freq_xy + SMALLISH_EPSILON < half_hethet_share + nonzero_freq_xx)) {
        sol_end_idx = 3;
        solutions[1] = (half_hethet_share + nonzero_freq_xy - nonzero_freq_xx) * 0.5;
        solutions[2] = half_hethet_share;
      } else {
        sol_end_idx = 2;
        solutions[1] = half_hethet_share;
      }
    }
    best_sol = solutions[sol_start_idx];
    if (sol_end_idx > sol_start_idx + 1) {
      // select largest log likelihood
      best_lnlike = calc_lnlike(known11, known12, known21, known22, center_ct_d, freq11, freq12, freq21, freq22, half_hethet_share, best_sol);
      cur_sol_idx = sol_start_idx + 1;
      do {
        incr_1122 = solutions[cur_sol_idx];
        cur_lnlike = calc_lnlike(known11, known12, known21, known22, center_ct_d, freq11, freq12, freq21, freq22, half_hethet_share, incr_1122);
        if (cur_lnlike > best_lnlike) {
          cur_lnlike = best_lnlike;
          best_sol = incr_1122;
        }
      } while (++cur_sol_idx < sol_end_idx);
    }
    if (onside_sol_ct_ptr && (sol_end_idx > sol_start_idx + 1)) {
      if (freqx1 * freq1x >= freq11) {
        dxx = freq1x * freqx1 - freq11;
        if (dxx > half_hethet_share) {
          dxx = half_hethet_share;
        }
      } else {
        dxx = 0.0;
      }
      // okay to NOT count suboptimal boundary points because they don't permit
      // direction changes within the main interval
      // this should exactly match haploview_blocks_classify()'s D sign check
      if ((freq11 + best_sol) - freqx1 * freq1x >= 0.0) {
        if (best_sol > dxx + SMALLISH_EPSILON) {
          lbound = dxx + SMALLISH_EPSILON;
        } else {
          lbound = dxx;
        }
        if (best_sol < half_hethet_share - SMALLISH_EPSILON) {
          half_hethet_share -= SMALLISH_EPSILON;
        }
      } else {
        if (best_sol > SMALLISH_EPSILON) {
          lbound = SMALLISH_EPSILON;
        } else {
          lbound = 0.0;
        }
        if (best_sol < dxx - SMALLISH_EPSILON) {
          half_hethet_share = dxx - SMALLISH_EPSILON;
        } else {
          half_hethet_share = dxx;
        }
      }
      for (cur_sol_idx = sol_start_idx; cur_sol_idx < sol_end_idx; cur_sol_idx++) {
        if (solutions[cur_sol_idx] < lbound) {
          sol_start_idx++;
        }
        if (solutions[cur_sol_idx] > half_hethet_share) {
          break;
        }
      }
      if (cur_sol_idx >= sol_start_idx + 2) {
        onside_sol_ct_ptr = cur_sol_idx - sol_start_idx;
      }
    }
    freq11 += best_sol;
  } else if ((prod_1122 == 0.0) && (prod_1221 == 0.0)) {
    final_sol[0] = 1; //true
    return final_sol;
  }
  final_sol[0] = 0; //false
  final_sol[1] = freq1x;
  final_sol[2] = freq2x;
  final_sol[3] = freqx1;
  final_sol[4] = freqx2;
  final_sol[5] = freq11;
  final_sol[6] = onside_sol_ct_ptr;
  return final_sol;
}




double calc_lnlike_quantile(double known11, double known12, double known21, double known22, double unknown_dh, double freqx1, double freq1x, double freq2x, double freq11_expected, double denom, int32_t quantile);
NumericVector em_phase_hethet(double known11, double known12, double known21, double known22, uint32_t center_ct, uint32_t onside_sol_ct_ptr);
// [[Rcpp::export]]
uint32_t haploview_blocks_classify(NumericVector counts, uint32_t lowci_max, uint32_t lowci_min, uint32_t recomb_highci, uint32_t strong_highci, uint32_t strong_lowci, uint32_t strong_lowci_outer, uint32_t is_x, double recomb_fast_ln_thresh) {
  // See comments in the middle of haploview_blocks().  The key insight is that
  // we only need to classify the D' confidence intervals into a few types, and
  // this almost never requires evaluation of all 101 log likelihoods.

  // Note that lowCI and highCI are *one-sided* 95% confidence bounds, i.e.
  // together, they form a 90% confidence interval.
  double known11 = (double)(2 * counts[0] + counts[1] + counts[3]);
  double known12 = (double)(2 * counts[2] + counts[1] + counts[5]);
  double known21 = (double)(2 * counts[6] + counts[3] + counts[7]);
  double known22 = (double)(2 * counts[8] + counts[5] + counts[7]);
  double total_prob = 0.0;
  double lnsurf_highstrong_thresh = 0.0;
  uint32_t onside_sol_ct = 1;
  double right_sum[83];
  double freq1x;
  double freq2x;
  double freqx1;
  double freqx2;
  double freq11_expected;
  double unknown_dh;
  double denom;
  double lnlike1;
  double lnsurf_highindiff_thresh;
  double dxx;
  double dyy;
  double dzz;
  uint32_t quantile;
  uint32_t center;
  NumericVector final_sol(7);
  // if (is_x) {
  //   known11 -= (double)((int32_t)counts[9]);
  //   known12 -= (double)((int32_t)counts[11]);
  //   known21 -= (double)((int32_t)counts[12]);
  //   known22 -= (double)((int32_t)counts[14]);
  // }
  final_sol = em_phase_hethet(known11, known12, known21, known22, counts[4], onside_sol_ct);
  if (final_sol[0]) {
    return 1;
  }
  final_sol[0] = 0; //false
  freq1x = final_sol[1];
  freq2x = final_sol[2];
  freqx1 = final_sol[3];
  freqx2 = final_sol[4];
  dzz = final_sol[5];
  onside_sol_ct = final_sol[6];
  freq11_expected = freqx1 * freq1x;
  dxx = dzz - freq11_expected;
  if (dxx < 0.0) {
    // D < 0, flip (1,1)<->(1,2) and (2,1)<->(2,2) to make D positive
    dyy = known11;
    known11 = known12;
    known12 = dyy;
    dyy = known21;
    known21 = known22;
    known22 = dyy;
    freq11_expected = freqx2 * freq1x;
    dyy = freqx1;
    freqx1 = freqx2;
    freqx2 = dyy;
    dxx = -dxx;
  }
  dyy = MINV(freqx1 * freq2x, freqx2 * freq1x);
  // this will always be in a term with a 0.01 multiplier from now on, so may
  // as well premultiply.
  denom = 0.01 * dyy;
  unknown_dh = (double)((int32_t)counts[4]);

  // force this to an actual likelihood array entry, so we know for sure
  // total_prob >= 1.0 and can use that inequality for both early exit and
  // determining the "futility threshold" (terms smaller than 2^{-53} / 19 are
  // too small to matter).
  center = (int32_t)(((dxx / dyy) * 100) + 0.5);

  lnlike1 = calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, center);

  // Previously assumed log likelihood was always concave, and used geometric
  // series bounds... then I realized this was NOT a safe assumption to make.
  // See e.g. rs9435793 and rs7531410 in 1000 Genomes phase 1.
  // So, instead, we only use an aggressive approach when onside_sol_ct == 1
  // (fortunately, that is almost always the case).
  if (onside_sol_ct == 1) {
    // It's not actually necessary to keep the entire likelihood array in
    // memory.  This is similar to the HWE and Fisher's exact test
    // calculations: we can get away with tracking a few partial sums, and
    // exploit unimodality, fixed direction on both sides of the center,
    // knowledge of the center's location, and the fact that we only need to
    // classify the CI rather than fully evaluate it.
    //
    // Specifically, we need to determine the following:
    // 1. Is highCI >= 0.98?  Or < 0.90?
    // 2. If highCI >= 0.98, is lowCI >= 0.81?  In [0.71, 0.81)?  Equal to
    //    0.70?  In [0.51, 0.70)?  In [0.01, 0.51)?  Or < 0.01?
    //    (Crucially, if highCI < 0.98, we don't actually need to determine
    //    lowCI at all.)
    // To make this classification with as few relative likelihood evaluations
    // as possible (5 logs, an exp call, 8 multiplies, 9 adds... that's kinda
    // heavy for an inner loop operation), we distinguish the following cases:
    // a. D' >= 0.41.  We first try to quickly rule out highCI >= 0.98 by
    //    inspection of f(0.97).  Then,
    //    * If it's below the futility threshold, jump to case (b).
    //    * Otherwise, sum f(0.98)..f(1.00), and then sum other likelihoods
    //      from f(0.96) on down.
    // b. D' < 0.41.  highCI >= 0.98 is impossible since f(0.41) >= f(0.42) >=
    //    ...; goal is to quickly establish highCI < 0.90.  A large fraction of
    //    the time, this can be accomplished simply by inspecting f(0.89); if
    //    it's less than 1/220, we're done because we know there's a 1
    //    somewhere in the array, and the sum of the likelihoods between
    //    f(0.89) and whatever that 1 entry is is bounded above by 12 * (1/220)
    //    due to fixed direction.  Otherwise, we sum from the top down.
    // This should be good for a ~10x speedup on the larger datasets where it's
    // most wanted.
    if (100 - center < 20 * (100 - strong_highci)) {
      dxx = calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, strong_highci) - lnlike1;
      // ln(2^{-53} / 19) is just under -39.6812
      if ((center > strong_highci) || (dxx > -39.6812)) {
        total_prob = exp(dxx);
        for (quantile = 100; quantile > strong_highci; quantile--) {
          total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
        }
        if (total_prob > (1.0 / 19.0)) {
          // branch 1: highCI might be >= 0.98
          lnsurf_highstrong_thresh = total_prob * 20;
          for (quantile = strong_highci - 1; quantile >= recomb_highci; quantile--) {
            total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
          }
          lnsurf_highindiff_thresh = total_prob * 20;
          while (1) {
            dxx = exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
            total_prob += dxx;
            // see comments on branch 2.  this is more complicated because we
            // still have work to do after resolving whether highCI >= 0.98,
            // but the reasoning is similar.
            if (total_prob >= lnsurf_highstrong_thresh) {
              if (quantile >= center) {
                goto haploview_blocks_classify_no_highstrong_1;
              }
              goto haploview_blocks_classify_no_highstrong_2;
            }
            if ((quantile <= lowci_max) && (quantile >= lowci_min)) {
              // We actually only need the [52..100], [71..100], [72..100], and
              // [82..100] right sums, but saving a few extra values is
              // probably more efficient than making this if-statement more
              // complicated.  [99 - quantile] rather than e.g. [quantile]
              // is used so memory writes go to sequentially increasing rather
              // than decreasing addresses.  (okay, this shouldn't matter since
              // everything should be in L1 cache, but there's negligible
              // opportunity cost)
              right_sum[quantile] = total_prob;
            }
            dxx *= ((int32_t)quantile);
            if (total_prob + dxx < lnsurf_highstrong_thresh) {
              while (1) {
                // Now we want to bound lowCI, optimizing for being able to
                // quickly establish lowCI >= 0.71.
                if (dxx * 19 < total_prob) {
                  // less than 5% remaining on left tail
                  if (quantile >= lowci_max) {
                    return 6;
                  }
                  while (quantile > lowci_min) {
                    quantile--;
                    total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
                    if (quantile <= lowci_max) {
                      right_sum[quantile] = total_prob;
                    }
                  }
                  dyy = right_sum[lowci_min] * (20.0 / 19.0);
                  while (total_prob < dyy) {
                    if ((!quantile) || (dxx <= RECIP_2_53)) {
                      total_prob *= 0.95;
                      if (total_prob >= right_sum[strong_lowci_outer]) {
                        // lowCI < 0.70
                        // -> f(0.00) + f(0.01) + ... + f(0.70) > 0.05 * total
                        return 3;
                      } else if (total_prob < right_sum[lowci_max]) {
                        return 6;
                      } else if ((lowci_max > strong_lowci) && (total_prob < right_sum[strong_lowci])) {
                        return 5;
                      }
                      return 4;
                    }
                    quantile--;
                    dxx = exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
                    total_prob += dxx;
                  }
                  return 2;
                }
                quantile--;
                dxx = exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
                total_prob += dxx;
                if ((quantile <= lowci_max) && (quantile >= lowci_min)) {
                  right_sum[quantile] = total_prob;
                }
                dxx *= ((int32_t)quantile);
              }
            }
            quantile--;
          }
        }
      }
      quantile = strong_highci - 1;
    } else {
      quantile = 100;
    }
    // branch 2: highCI guaranteed less than 0.98.  If D' <= 0.875, try to
    // quickly establish highCI < 0.90.
    dxx = calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, recomb_highci) - lnlike1;
    if ((center < recomb_highci) && (dxx < recomb_fast_ln_thresh)) {
      return 0;
    }
    // okay, we'll sum the whole right tail.  May as well sum from the outside
    // in here for a bit more numerical stability, instead of adding exp(dxx)
    // first.
    do {
      total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
    } while (--quantile > recomb_highci);
    total_prob += exp(dxx);
    lnsurf_highindiff_thresh = total_prob * 20;
    haploview_blocks_classify_no_highstrong_1:
      quantile--;
    if (center < recomb_highci) {
      // if we know there's a 1.0 ahead in the likelihood array, may as well
      // take advantage of that
      lnsurf_highstrong_thresh = lnsurf_highindiff_thresh - 1.0;
      while (quantile > center) {
        total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
        if (total_prob >= lnsurf_highstrong_thresh) {
          return 0;
        }
        quantile--;
      }
      if (!center) {
        return 1;
      }
      total_prob += 1;
      quantile--;
    }
    // likelihoods are now declining, try to exploit that to exit early
    // (it's okay if the first likelihood does not represent a decline)
    while (1) {
      dxx = exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
      total_prob += dxx;
      haploview_blocks_classify_no_highstrong_2:
        if (total_prob >= lnsurf_highindiff_thresh) {
          return 0;
        }
        if (total_prob + ((int32_t)quantile) * dxx < lnsurf_highindiff_thresh) {
          // guaranteed to catch quantile == 0
          return 1;
        }
        quantile--;
    }
  }
  for (quantile = 100; quantile >= recomb_highci; quantile--) {
    total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
    if (quantile == strong_highci) {
      lnsurf_highstrong_thresh = total_prob * 20;
    }
  }
  if (total_prob < (1.0 / 19.0)) {
    return 0;
  }
  lnsurf_highindiff_thresh = total_prob * 20;
  while (1) {
    total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
    if (total_prob >= lnsurf_highindiff_thresh) {
      return 0;
    }
    if (quantile <= lowci_max) {
      if (quantile >= lowci_min) {
        right_sum[quantile] = total_prob;
      } else if (!quantile) {
        break;
      }
    }
    quantile--;
  }
  if (total_prob >= lnsurf_highstrong_thresh) {
    return 1;
  }
  total_prob *= 0.95;
  if (total_prob < right_sum[strong_lowci]) {
    if ((lowci_max > strong_lowci) && (total_prob >= right_sum[lowci_max])) {
      return 5;
    }
    return 6;
  }
  if (total_prob >= right_sum[strong_lowci_outer]) {
    if ((lowci_min < strong_lowci_outer) && (total_prob >= right_sum[lowci_min])) {
      return 2;
    }
    return 3;
  }
  return 4;
}


uint32_t haploview_blocks_classify(uint32_t counts, uint32_t lowci_max, uint32_t lowci_min, uint32_t recomb_highci, uint32_t strong_highci, uint32_t strong_lowci, uint32_t strong_lowci_outer, uint32_t is_x, double recomb_fast_ln_thresh);
// [[Rcpp::export]]
uint32_t CIDp_strLD(IntegerVector b1, IntegerVector b2, double lower, double upper){
  NumericMatrix counts012(3,3);
  NumericVector counts(9);
  int b1size = b1.size();
  for (int i = 0; i < b1size; i++) {
    if(b1[i]!=NA_INTEGER && b2[i]!=NA_INTEGER){
      counts012(b1[i],b2[i])++;
    }
  }
  // return counts012;
  counts[0] = counts012(0,0);
  counts[1] = counts012(0,1);
  counts[2] = counts012(0,2);
  counts[3] = counts012(1,0);
  counts[4] = counts012(1,1);
  counts[5] = counts012(1,2);
  counts[6] = counts012(2,0);
  counts[7] = counts012(2,1);
  counts[8] = counts012(2,2);
  uint32_t strLD;
  uint32_t lowci_max = 82;
  uint32_t lowci_min = 52;
  uint32_t recomb_highci = 89;
  double recomb_fast_ln_thresh = -log((int32_t)((100 - recomb_highci) * 20));
  upper = upper*100-1;
  lower = lower*100-1;
  uint32_t strong_highci = (uint32_t)(upper);
  uint32_t strong_lowci_outer = (uint32_t)(lower);
  uint32_t strong_lowci = strong_lowci_outer+1;
  uint32_t is_x = 0;
  strLD = haploview_blocks_classify(counts, lowci_max, lowci_min, recomb_highci, strong_highci, strong_lowci, strong_lowci_outer, is_x ,recomb_fast_ln_thresh);
  return(strLD);
}
