/*
 * EigenTools.h
 *
 *  Created on: 09-01-2011
 *      Author: Jacek Cyranka
 */

#ifndef _CAPD_JACO_EIGEN_H_
#define _CAPD_JACO_EIGEN_H_

#include "capd2alglib.h"
#include "capd/alglib/hsschur.h"
#include "capd/alglib/hessenberg.h"


namespace capd{
namespace jaco{

template<typename DoubleT, int D_>
class EigenTools{
public:

	typedef DoubleT ScalarType;
	typedef typename capd::vectalg::Vector<DoubleT,0> VectorType;
	typedef typename capd::vectalg::Matrix<DoubleT,0,0> MatrixType;
    typedef typename capd::vectalg::Matrix<DoubleT,0,0> DoubleMatrixType;

    typedef typename capd::vectalg::Vector< capd::intervals::Interval< double > ,0> IntervalVectorType;
    typedef typename capd::vectalg::Matrix< capd::intervals::Interval< double > ,0,0> IntervalMatrixType;

	bool schurDecompose(const MatrixType& A, MatrixType& T, MatrixType& S) const;

	///Solves linear system of equations Ax=b.
	void linSolve(const MatrixType& A, VectorType& x, const VectorType& b) const;

	///Calculates the Schur decomposition of the matrix A.
	bool schurMatrix(const MatrixType& A, MatrixType& T, MatrixType& S, std::ostream& out) const;

    //krawczykRefineMatrix(/*macierz ktora chce odwrocic*/, /*pzryblizona macierz odwrotna np. z gausa na doublach*/, /*macierz */)
    bool krawczykRefineMatrix(DoubleMatrixType& A, const DoubleMatrixType& approxAinv, IntervalMatrixType& Ainv) const;
    bool krawczykRefineColumn( DoubleMatrixType& A, const DoubleMatrixType& Ainv, IntervalVectorType& x, int columnIndex) const;    

};

///==========================================function definitions====================================================



template < class DoubleT, int D_>
 bool EigenTools< DoubleT, D_>::krawczykRefineColumn( DoubleMatrixType& A, const DoubleMatrixType& Ainv, IntervalVectorType& x, int columnIndex) const {
   int dim = A.numberOfRows();
   IntervalVectorType m(dim), r(dim);
   IntervalMatrixType der = IntervalMatrixType(A); /// F(x)=A*x, dF(x)=A
   DoubleMatrixType C=Ainv; /// C=dF(x)^{-1}=A^{-1}
   IntervalMatrixType Id = IntervalMatrixType::Identity(dim);
   capd::vectalg::split(x, m, r);
   IntervalVectorType fm = IntervalMatrixType(A) * m - capd::vectalg::transpose(Id)[columnIndex];
   IntervalVectorType K = m - IntervalMatrixType(C) * fm  + ( Id - IntervalMatrixType(C) * der ) * r;
   if(capd::vectalg::intersection(x, K, x))
     return true;
   else
     return false;
 }

 template < class DoubleT, int D_>
 bool EigenTools< DoubleT, D_>::krawczykRefineMatrix(DoubleMatrixType& A, const DoubleMatrixType& approxAinv, IntervalMatrixType& Ainv) const {
   int i;
   IntervalVectorType x;
   int dim=Ainv.numberOfColumns();
   IntervalMatrixType AinvT=capd::vectalg::transpose(Ainv),
              candidate;
   for(i=0; i<dim; i++){
     x=AinvT[i];
     if(!krawczykRefineColumn( A, approxAinv, x, i ))
       return false;
     AinvT[i]=x;
   }
   candidate=capd::vectalg::transpose(AinvT);
   intersection(Ainv, candidate, Ainv);
   return true;
 }



template<typename DoubleT, int D_>
bool EigenTools<DoubleT, D_>::schurDecompose(const MatrixType& A, MatrixType& T, MatrixType& S) const{
   if(A.numberOfRows() != A.numberOfColumns())
      throw std::invalid_argument("schurDecompose works only for square matrices");
    ap::real_2d_array a;
    int n = A.numberOfRows();
    a.setbounds(0,n-1, 0, n-1) ;
    int i,j,k;
    for(i =0; i < n; i++)
      for(j=0; j< n; j++){
        a(i,j) = A[i][j];
      }
    ap::real_2d_array t;
    ap::real_2d_array s;
    ap::real_1d_array tau;
    ap::real_2d_array h;
    ap::real_2d_array h1;
    ap::real_2d_array s1;
    ap::real_2d_array s2;
    ap::real_2d_array s3;
    h1.setbounds(1, n, 1, n);
    s1.setbounds(1, n, 1, n);
    s2.setbounds(1, n, 1, n);
    s3.setbounds(1, n, 1, n);

    rmatrixhessenberg(a, n, tau);

    rmatrixhessenbergunpackq(a, n, tau, s);

    rmatrixhessenbergunpackh(a, n, h);

    //changing indexing of h and s from [0...n-1] to [1...n]
    for(i=0; i<n; i++)
      for(j=0; j<n; j++){
        h1(i+1, j+1)=h(i,j);
        s1(i+1, j+1)=s(i,j);
      }

    if(! upperhessenbergschurdecomposition(h1, n, s2))
        throw std::runtime_error("algorithm for schur decomposition did not converge!");

    ///multiplying orthogonal s1 by orthogonal s2
    for(i=1; i<=n; i++)
      for(j=1; j<=n; j++){
        s3(i, j)=0;
        for(k=1; k<=n; k++){
          s3(i, j)+=s1(i, k)*s2(k, j);
        }
      }

    for(int i =0; i < n; i++){
      for(int j=0; j< n; j++){
        T[j][i]=h1(j+1,i+1);
        S[j][i]=s3(j+1,i+1);
      }
    }
    return true;
}

template<typename DoubleT, int D_>
void EigenTools<DoubleT, D_>::linSolve(const MatrixType& A, VectorType& x, const VectorType& b) const{
  int dim=x.size();
  MatrixType L(dim, dim),
         U(dim, dim);
  capd::matrixAlgorithms::croutDecomposition(A, L, U); //decomposes matrix into lower and upper triangular matrices
  ScalarType s;
  int i, j;
  VectorType y(dim);
  //Ly=b
  for(i=1; i<=dim; ++i){
    s=b[i-1];
    for(j=1; j<i; ++j){
      s-=L(i, j)*y[j-1];
    }
    y[i-1]=s/L(i, i);
  }
  //Ux=y
  for(i=dim; i>=1; --i){
    s=y[i-1];
    for(j=dim; j>=i+1; --j){
      s-=U(i, j)*x[j-1];
    }
    x[i-1]=s/U(i, i);
  }
}

template<typename DoubleT, int D_>
bool EigenTools<DoubleT, D_>::schurMatrix(const MatrixType& A, MatrixType& T, MatrixType& S, std::ostream& out) const{
  bool r=schurDecompose(A, T, S);
  return r;
}

}}

#endif /* _CAPD_JACO_EIGEN_H_ */
