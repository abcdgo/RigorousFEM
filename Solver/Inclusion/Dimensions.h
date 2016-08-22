#ifndef _CAPD_DIMENSIONS_H_
#define _CAPD_DIMENSIONS_H_

#include "config.h"

namespace capd{
namespace jaco{

///class storing dimensions and providing useful functions.
class Dimensions
{
public:
  int m; ///Galerkin projection dimension
  int M; ///near tail dimension

  Dimensions(int m, int M){
    this->m=m;
    this->M=M;
  }

  ///constructor that uses default dimensions (from the config.h file)
  Dimensions(){
    this->m=_m;
    this->M=_M;
  }

  inline int getm() const{
    return this->m;
  }

  inline int getM() const{
    return this->M;
  }

  bool inTail(int index) const;

  template< typename ScalarType>
  ScalarType norm(const ScalarType& real, const ScalarType& complex) const;
  
  template <typename VectorType>
  double maxDiam(const VectorType& v) const;

};

///==========================================function definitions====================================================
bool Dimensions::inTail(int index) const{
  if(index>this->m || index<-this->m)
    return true;
  return false;
}

template< typename ScalarType>
ScalarType Dimensions::norm(const ScalarType& real, const ScalarType& complex) const{
  ScalarType t=power(real, 2)+power(complex, 2);
  ScalarType result=sqrt(abs(t));
  return ScalarType(result.rightBound());
}

template <typename VectorType>
double Dimensions::maxDiam(const VectorType& v) const{
  int i;
  double max=0;
  for(i=0; i<v.size(); ++i){
    if(diam(v[i]).rightBound()>max) max=diam(v[i]).rightBound();
  }
  return max;
}

}}

#endif 
