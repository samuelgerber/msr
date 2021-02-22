#ifndef NULL
#define NULL 0
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "NNMSComplex.h"
#include "NNMSComplex2old.h"

#include "KernelDensity.h"
#include "GaussianKernel.h"

extern "C" {


int storeOneLevel(NNMSComplex<double> &msc, int n, SEXP list, int startIndex){
  using namespace FortranLinalg;
  SEXP crystals;
  PROTECT(crystals = Rf_allocVector(INTSXP, n));
  DenseVector<int> c(n, INTEGER(crystals));
  msc.getPartitions(c);
  //convert to 1 index
  Linalg<int>::Add(c, 1, c);

  SEXP ascending;
  PROTECT(ascending = Rf_allocVector(INTSXP, n));
  DenseVector<int> vascending(n, INTEGER(ascending));
  msc.getAscending(vascending);
  //convert to 1 index
  Linalg<int>::Add(vascending, 1, vascending);

  int nCrystals = msc.getNCrystals();

  SEXP mins;
  PROTECT(mins = Rf_allocVector(INTSXP, nCrystals));
  DenseVector<int> vmins(nCrystals, INTEGER(mins));
  msc.getMin(vmins);
  Linalg<int>::Add(vmins, 1, vmins);

  SEXP maxs;
  PROTECT(maxs = Rf_allocVector(INTSXP, nCrystals));
  DenseVector<int> vmaxs(nCrystals, INTEGER(maxs));
  msc.getMax(vmaxs);
  Linalg<int>::Add(vmaxs, 1, vmaxs);

  SET_VECTOR_ELT(list, startIndex  , crystals);
  SET_VECTOR_ELT(list, startIndex+1, maxs);
  SET_VECTOR_ELT(list, startIndex+2, mins);
  SET_VECTOR_ELT(list, startIndex+3, ascending);
  return 4;
}

int storeOneLevelR2(NNMSComplexR2<double> &msc, int n, SEXP list, int startIndex){
  using namespace FortranLinalg;
  SEXP crystals;
  PROTECT(crystals = Rf_allocVector(INTSXP, n));
  DenseVector<int> c(n, INTEGER(crystals));
  msc.getPartitions(c);
  //convert to 1 index
  Linalg<int>::Add(c, 1, c);

  SEXP ascending;
  PROTECT(ascending = Rf_allocVector(INTSXP, n));
  //DenseVector<int> vascending(n, INTEGER(ascending));
  //msc.getAscedning(vascending);
  //convert to 1 index
  //Linalg<int>::Add(vascending, 1, vascending);

  int nCrystals = msc.getNCrystals();
  SEXP mins;
  PROTECT(mins = Rf_allocVector(INTSXP, nCrystals));
  DenseVector<int> vmins(nCrystals, INTEGER(mins));
  msc.getMin(vmins);
  Linalg<int>::Add(vmins, 1, vmins);

  SEXP maxs;
  PROTECT(maxs = Rf_allocVector(INTSXP, nCrystals));
  DenseVector<int> vmaxs(nCrystals, INTEGER(maxs));
  msc.getMax(vmaxs);
  Linalg<int>::Add(vmaxs, 1, vmaxs);

  int nExt = msc.getNAllExtrema();
  SEXP ps;
  PROTECT(ps = Rf_allocVector(REALSXP, nExt-1));
  DenseVector<double> vps(nExt-1, REAL(ps));
  msc.getPersistence(vps);

  SET_VECTOR_ELT(list, startIndex  , crystals);
  SET_VECTOR_ELT(list, startIndex+1, maxs);
  SET_VECTOR_ELT(list, startIndex+2, mins);
  SET_VECTOR_ELT(list, startIndex+3, ascending);
  return 4;
}


//-----  nn based ms decompositions  ----/

SEXP nnmspartition (SEXP Rm, SEXP Rn, SEXP Ry, SEXP Rx,
    SEXP Rknn, SEXP RpLevel, SEXP Rsmooth, SEXP knnEps) {
     using namespace FortranLinalg;

  int knn = *INTEGER(Rknn);
  int n = *INTEGER(Rn);
  int m = *INTEGER(Rm);
  double *x = REAL(Rx);
  double *y = REAL(Ry);
  double pLevel = *REAL(RpLevel);
  double eps = *REAL(knnEps);
  bool smooth = *LOGICAL(Rsmooth);

  if(knn > n){
    knn = n;
  }
  DenseMatrix<double> X(m, n, x);
  DenseVector<double> Y(n, y);


  NNMSComplex<double> msc(X, Y, knn, smooth, eps);
  msc.mergePersistence(pLevel);

  SEXP list;
  PROTECT( list = Rf_allocVector(VECSXP, 5));

  int nExt = msc.getNAllExtrema();
  SEXP ps;
  PROTECT(ps = Rf_allocVector(REALSXP, nExt-1));
  DenseVector<double> vps(nExt-1, REAL(ps));
  msc.getPersistence(vps);
  SET_VECTOR_ELT(list, 0, ps);
  storeOneLevel(msc, n,list, 1);
  UNPROTECT(6);

  msc.cleanup();

  return list;
}




SEXP nnmsc2(SEXP Rm, SEXP Rn, SEXP Ry, SEXP Rx,
    SEXP Rknn, SEXP nLevels, SEXP Rsmooth, SEXP knnEps) {
     using namespace FortranLinalg;

  int knn = *INTEGER(Rknn);
  int n = *INTEGER(Rn);
  int m = *INTEGER(Rm);
  double *x = REAL(Rx);
  double *y = REAL(Ry);
  int nL = *INTEGER(nLevels);
  double eps = *REAL(knnEps);
  bool smooth = *LOGICAL(Rsmooth);

  if(knn > n){
    knn = n;
  }
  DenseMatrix<double> X(m, n, x);
  DenseVector<double> Y(n, y);


  NNMSComplex<double> msc(X, Y, knn, smooth, eps);

  int nExt = msc.getNAllExtrema();
  SEXP ps;
  PROTECT(ps = Rf_allocVector(REALSXP, nExt-1));
  DenseVector<double> vps(nExt-1, REAL(ps));
  msc.getPersistence(vps);
  if(nL > vps.N()){
    nL = vps.N();
  }


  SEXP list;
  PROTECT( list = Rf_allocVector(VECSXP, 1+4*nL));
  SET_VECTOR_ELT(list, 0, ps);
  int nProtect = 2;
  for(int i=0; i < nL; i++){
    msc.mergePersistence(vps(vps.N()-1-i));
    storeOneLevel(msc, n, list, i*4+1);
    nProtect += 4;
  }

  UNPROTECT(nProtect);

  msc.cleanup();

  return list;
}



//------ graph input ms decompsotion -----//

SEXP graphmsc(SEXP Rm, SEXP Rn, SEXP Rk, SEXP Ry,
    SEXP Rx, SEXP Rknn, SEXP Rknnd, SEXP nLevels, SEXP Rsmooth) {


     using namespace FortranLinalg;


  int k = *INTEGER(Rk);
  int *knn = INTEGER(Rknn);
  double *knnd = REAL(Rknnd);
  int n = *INTEGER(Rn);
  int m = *INTEGER(Rm);
  double *x = REAL(Rx);
  double *y = REAL(Ry);
  int nL = *INTEGER(nLevels);
  bool smooth = *LOGICAL(Rsmooth);

  DenseMatrix<double> KNND(k, n, knnd);
  DenseMatrix<int> KNN(k, n, knn);
  DenseMatrix<double> X(m, n, x);
  DenseVector<double> Y(n, y);


  NNMSComplex<double> msc(X, Y, KNN, KNND, smooth);



  int nExt = msc.getNAllExtrema();
  SEXP ps;
  PROTECT(ps = Rf_allocVector(REALSXP, nExt-1));
  DenseVector<double> vps(nExt-1, REAL(ps));
  msc.getPersistence(vps);
  if(nL > vps.N()){
    nL = vps.N();
  }




  SEXP list;
  PROTECT( list = Rf_allocVector(VECSXP, 1+4*nL));
  SET_VECTOR_ELT(list, 0, ps);
  int nProtect = 2;
  for(int i=0; i < nL; i++){
    msc.mergePersistence(vps(vps.N()-1-i));
    storeOneLevel(msc, n, list, i*4+1);
    nProtect += 4;
  }

  UNPROTECT(nProtect);


  msc.cleanup();

  return list;
}








//merging based on R^2 of linear models
SEXP nnmscR2(SEXP Rm, SEXP Rn, SEXP Ry, SEXP Rx, SEXP Rknn, SEXP Rsmooth, SEXP knnEps) {


  using namespace FortranLinalg;

  int knn = *INTEGER(Rknn);
  int n = *INTEGER(Rn);
  int m = *INTEGER(Rm);
  double *x = REAL(Rx);
  double *y = REAL(Ry);
  double eps = *REAL(knnEps);
  bool smooth = *LOGICAL(Rsmooth);

  if(knn > n){
    knn = n;
  }
  DenseMatrix<double> X(m, n, x);
  DenseVector<double> Y(n, y);


  NNMSComplexR2<double> msc(X, Y, knn, -1, smooth, eps);

  int nExt = msc.getNAllExtrema();
  SEXP ps;
  PROTECT(ps = Rf_allocVector(REALSXP, nExt-1));
  DenseVector<double> vps(nExt-1, REAL(ps));
  msc.getPersistence(vps);

  SEXP list;
  PROTECT( list = Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(list, 0, ps);
  storeOneLevelR2(msc, n, list, 1);
  UNPROTECT(6);

  msc.cleanup();

  return list;
}


SEXP nnmsc2R2(SEXP Rm, SEXP Rn, SEXP Ry, SEXP Rx,
    SEXP Rknn, SEXP nLevels, SEXP Rsmooth, SEXP knnEps) {
  using namespace FortranLinalg;

  int knn = *INTEGER(Rknn);
  int n = *INTEGER(Rn);
  int m = *INTEGER(Rm);
  double *x = REAL(Rx);
  double *y = REAL(Ry);
  int nL = *INTEGER(nLevels);
  double eps = *REAL(knnEps);
  bool smooth = *LOGICAL(Rsmooth);

  if(knn > n){
    knn = n;
  }
  DenseMatrix<double> X(m, n, x);
  DenseVector<double> Y(n, y);


  NNMSComplexR2<double> msc(X, Y, knn, nL, smooth, eps);

  int nExt = msc.getNAllExtrema();
  SEXP ps;
  PROTECT(ps = Rf_allocVector(REALSXP, nExt-1));
  DenseVector<double> vps(nExt-1, REAL(ps));
  msc.getPersistence(vps);
  if(nL > vps.N()){
    nL = vps.N();
  }
  int nEntries = 4;

  SEXP list;
  PROTECT( list = Rf_allocVector(VECSXP, 1+nEntries*nL));
  SET_VECTOR_ELT(list, 0, ps);

  int nProtect = 2;
  for(int i=0; i < nL; i++){
    msc.setLevel(i);
    storeOneLevelR2(msc, n, list, i*4+1);
    nProtect += 4;
  }
  UNPROTECT(nProtect);

  msc.cleanup();

  return list;
}



void gkde(int *m, int *nx, double *x, double *sigma, int *nxeval, double *xeval, double *peval){
  using namespace FortranLinalg;
  DenseMatrix<double> X(*m, *nx, x);
  DenseMatrix<double> Xe(*m, *nxeval, xeval);

  GaussianKernel<double> gaussian(*sigma, *m);
  KernelDensity<double> kde(X, gaussian);
  DenseVector<double> p(*nxeval, peval);
  for(int i=0; i<p.N(); i++){
    p(i) = kde.p(Xe, i)/X.N();
  }
}

}//end extern C
