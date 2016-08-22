#ifndef _CAPD_D1_REAL_H_
#define _CAPD_D1_REAL_H_

#include "Dimensions.h"

namespace capd{
namespace jaco{

///this class is base for all classes that operate on the Fourier modes of one dimensional
///real functions with periodic bd. conditions (invariant subspace a_k=\overline{a_{-k}}).
///Modes (Galerkin projection, near tail) are always stored in an array, and this class
///provides methods to retrieve from this array indicated modes and vice versa. Moreover it
///assumes that arrays stores modes of positive index only, whereas those of neagative index
///are retrieved by taking the conjugate of corresponding positive indexed mode.
template<class ScalarT, class VectorT=capd::vectalg::Vector<ScalarT, _D> >
class D1Real : public Dimensions
{
public:
  typedef typename capd::jaco::D1Real<ScalarT, VectorT> Class;
  typedef VectorT VectorType;
  typedef ScalarT ScalarType;


  D1Real(int m, int M) : Dimensions(m, M)
  {}

  ///constructor that uses default dimensions
  D1Real() : Dimensions(_m, _M)
  {}

  ///returns index of the first mode from the Galerkin projection that is stored.
  inline int firstMode() const
  {
    return 0;
  }

  ///returns index of the last mode from the Galerkin projection that is stored.
  inline int lastMode() const
  {
    return this->m - 1;
  }

  ///returns index of the first mode from the near tail that is stored.
  inline int firstModeInTail() const
  {
    return this->m;
  }

  ///returns index of the last mode in the near tail that is stored.
  inline int lastModeInTail() const
  {
    return this->M - 1;
  }

  ///returns size of an array needed to store 'modes' modes from a Galerkin projection.
  inline int modes2arraySize(int modes) const
  {
    //return 2*modes;
    return modes;
  }

  ///returns size of an array needed to store m modes from a Galerkin projection.
  inline int modes2arraySize() const
  {
    //return 2*this->m;
    return this->m;
  }

  ///translates index in an array storing modes from a Galerkin projection into index of a corresponding mode,
  ///i.e. returns k: Re{a_k} or Im{a_k} is stored in table[i].
  inline int array2mode(int i) const;
  
  ///translates index in an array storing modes from a Galerkin projection into index of a corresponding mode,
  ///i.e. returns k: a_k is stored in table[i]. Regarding that we store real and imaginary parts separately,
  ///information whether table[i] stores Re{a_k} is passed into the variable re.
  inline int array2mode(int i, bool& real) const;

  ///translates index k of a mode a_k from the Galerkin projection into an index in an array storing modes,
  ///i.e. if re==true returns i: Re{a_k} is stored in table[i],
  ///else returns i: Im{a_k} is stored in table[i].
  inline int mode2array(int k, bool re) const;

  ///translates index k of a mode a_k from the near tail into an index in an array storing modes,
  ///i.e. if re==true returns i: Re{a_k} is stored in table[i],
  ///else returns i: Im{a_k} is stored in table[i].
  int mode2arrayTail(int k, bool re) const;

  ///translates index in an array storing modes from a near tail into index of a corresponding mode,
  ///i.e. returns k: Re{a_k} or Im{a_k} is stored in table[i].
  int array2modeTail(int i) const;

  ///translates index in an array storing modes from a near tail into index of a corresponding mode,
  ///i.e. returns k: a_k is stored in table[i]. Regarding that we store real and imaginary parts separately,
  ///information whether table[i] stores Re{a_k} is passed into the variable re.
  int array2modeTail(int i, bool& re) const;

  ///prints modes from a Galerkin projection stored in interval vector v to the output indicated by out,
  ///along with modes energy.
  void printModes(const VectorType& v, const ScalarType& energy, std::ostream& out) const;

  ///prints modes from a Galerkin projection stored in interval vector v to the standard output, along with
  ///modes energy.
  void printModes(const VectorType& v, const ScalarType& energy) const;

  ///prints modes from a Galerkin projection stored in interval vector v to the output indicated by out.
  void printModes(const VectorType& v, std::ostream& out) const;

  ///prints modes from a Galerkin projection stored in interval vector v to the output indicated by out
  ///with additional information regarding symmetry of all modes with respect to zero.
  void printModesWithSymmetryCheck(const VectorType& v, std::ostream& out) const;

  ///more compact version used for printing on standard output.
  void printModes(const VectorType& v) const;

  //returns value of k-th mode from all-over symmetric Galerkin projection,
  ///if re==true return Re{a_k},
  ///else return Im{a_k}.
  ///Value of a_k for k=1,\dots,m is stored explicitely,
  ///whereas value of -a_k is obtained by taking conjugate of a_k.
  template<typename AVector>
  typename AVector::ScalarType modeVal(int k, bool re, const AVector& vec, const ScalarType& a0) const;

  ///returns value of k-th mode in the near tail,
  ///if re==true return Re{a_k},
  ///else return Im{a_k}.
  ///Value of a_k for k=m+1,\dots,M is stored explicitely,
  ///whereas value of -a_k is obtained by taking conjugate of a_k.
  template<typename AVector>
  typename AVector::ScalarType closeTailVal( int k, bool re, const AVector& vec) const;

  ///calculates sum of all modes supremum,
  ///returns \sum_{k=-m}^{k=m}{\sup{|a_k|}} without a_0 mode.
  template<typename AVector, typename TailType>
  typename AVector::ScalarType sumOfSup(const AVector& modes, const TailType& tail, const ScalarType& a0) const;

  ///calculates energy,
  ///returns \sum_{k=-m}^{k=m}{a_k} without a_0 mode.
  template<typename AVector>
  typename AVector::ScalarType energy(const AVector& modes);
  
  ///calculates energy,
  ///returns \sum_{k=-M}^{k=M}{a_k} without a_0 mode.
  template<typename AVector, typename TailType>
  typename AVector::ScalarType energy(const AVector& modes, const TailType& tail);

  ///function calculating real part of the Fourier series  at x, considering modes from the Galerkin projection only,
  ///returns \sum_{k=-m}^m{a_k e^{ikx}}.
  template<typename AVector>
  typename AVector::ScalarType projectionFourierSeries(const ScalarType& x,const AVector& modes);

  ///function calculating real part of the Fourier series  at x, considering modes from the near tail only,
  ///returns \sum_{k=-M}^{-m}{a_k e^{ikx}}+\sum_{k=m}^{M}{a_k e^{ikx}}.
  template<typename AVector>
  typename AVector::ScalarType closeTailFourierSeries(const ScalarType& x, const AVector& modes);

};

///==========================================function definitions====================================================

template<typename ScalarT, typename VectorT>
inline int D1Real<ScalarT, VectorT>::array2mode(int k) const
{
  return k;
}

template<typename ScalarT, typename VectorT>
inline int D1Real<ScalarT, VectorT>::array2mode(int k, bool& real) const
{
  return k;
}

template<typename ScalarT, typename VectorT>
inline int D1Real<ScalarT, VectorT>::mode2array(int k, bool re) const
{
#if __GENERAL_DEBUG__ //flag here because of efficiency reasons
  if(k<firstMode() || k>lastMode()){
    std::cerr<<"Error in mode2array function. Requested index is out of the projection. Index="<<k<<", m="<<this->m<<".\n";
    throw std::runtime_error("Error in mode2array function. Requested index is out of the projection.\n");
  }
#endif

  return k;
}

template<typename ScalarT, typename VectorT>
int D1Real<ScalarT, VectorT>::mode2arrayTail(int k, bool re) const
{
  //it is under debug flag because of efficiency
#if __GENERAL_DEBUG__
    if(k<=this->m && k>=-this->m){ //we have to check if index k in fact corresponds to mode in Tail
      std::cerr<<"Requested index is out of the tail tail. Index="<<k<<", m="<<this->m<<".\n";
      throw std::runtime_error("Requested index is out of the tail tail.\n");
    }
#endif
  return k;
}

template<typename ScalarT, typename VectorT>
int D1Real<ScalarT, VectorT>::array2modeTail(int k) const
{
  return k;
}

template<typename ScalarT, typename VectorT>
int D1Real<ScalarT, VectorT>::array2modeTail(int k, bool& re) const
{
  return k;
}

template<typename ScalarT, typename VectorT>
void D1Real<ScalarT, VectorT>::printModes(const VectorType& v, const ScalarType& energy, std::ostream& out) const
{
  int i;
  double max=0;
  for(i=0; i<v.size(); i++){
    out<<"a_"<<this->array2mode(i)<<": "<<v[i]<<" diam: "<<diam(v[i]).rightBound()<<"\n";
    if(diam(v[i]).rightBound()>max) max=diam(v[i]).rightBound();
  }
  out<<"Energy="<<energy<<" diam="<<diam(energy).rightBound()<<"\n";
  out<<"Max diam="<<max<<"\n";
}

template<typename ScalarT, typename VectorT>
void D1Real<ScalarT, VectorT>::printModes(const VectorType& v, const ScalarType& energy) const
{
  int i;
  double max=0;
  for(i=0; i<v.size()-1; i+=2){
    std::cout<<"Re(a_"<<this->array2mode(i)<<"): "<<v[i]<<" Im(a_"<<this->array2mode(i)<<"): "<<v[i+1]<<"\n";
    if(diam(v[i]).rightBound()>max) max=diam(v[i]).rightBound();
    if(diam(v[i+1]).rightBound()>max) max=diam(v[i+1]).rightBound();
  }
  std::cout<<"Energy="<<energy<<" diam="<<diam(energy).rightBound()<<"\n";
  std::cout<<"Max diam="<<max<<"\n";
}

template<typename ScalarT, typename VectorT>
void D1Real<ScalarT, VectorT>::printModes(const VectorType& v, std::ostream& out) const
{
  int i;
    for(i=0; i<v.size(); i++){
    out<<"a_"<<this->array2mode(i)<<": "<<v[i]<<" diam: "<<diam(v[i]).rightBound()<<"\n";
  }
}

template<typename ScalarT, typename VectorT>
void D1Real<ScalarT, VectorT>::printModesWithSymmetryCheck(const VectorType& v, std::ostream& out) const
{
  int i;
  for(i=0; i<v.size()-1; i+=2){
    out<<"Re(a_"<<this->array2mode(i)<<"): "<<v[i]<<" diam: "<<diam(v[i]).rightBound()<<", Im(a_"<<this->array2mode(i)<<"): "<<v[i+1]<<" diam: "<<diam(v[i+1]).rightBound()<<" symmetry: "<<v[i+1].leftBound()+v[i+1].rightBound()<<"\n";
  }
}


template<typename ScalarT, typename VectorT>
void D1Real<ScalarT, VectorT>::printModes(const VectorType& v) const
{
  int i;
  for(i=0; i<v.size()-1; i+=2){
    std::cout<<"Re(a_"<<this->array2mode(i)<<"): "<<v[i]<<" Im(a_"<<this->array2mode(i)<<"): "<<v[i+1]<<"\n";
  }
}


template<typename ScalarT, typename VectorT>
template<typename AVector>
typename AVector::ScalarType D1Real<ScalarT, VectorT>::modeVal(int k, bool re, const AVector& vec, const ScalarT& a0) const
{
  if(k==0){ //the constant zero mode
    if(re) return a0;
    else return 0;
  }else{
    if(k>0){ //we store a_k for k>0
      return vec[this->mode2array(k, re)];
    }else{  //we take conjugate of a_k
      //returns conjugate of a k-th mode
      int index=this->mode2array(-k, re);
      if(re)
        return vec[index];
      else
        return -vec[index];
    }
  }
}

template<typename ScalarT, typename VectorT>
template<typename AVector>
typename AVector::ScalarType D1Real<ScalarT, VectorT>::closeTailVal( int k, bool re, const AVector& vec) const
{
  if(k>0){
    return vec[this->mode2arrayTail(k, re)];
  }else{
    //returns conjugate of a k-th mode
    int index=this->mode2arrayTail(-k, re);
    if(index%2==0){
      return vec[index];
    }else{
      return -vec[index];
    }
  }
}

template<typename ScalarT, typename VectorT>
template<typename AVector, typename TailType>
typename AVector::ScalarType D1Real<ScalarT, VectorT>::sumOfSup(const AVector& modes, const TailType& tail, const ScalarType& a0) const
{
  ScalarType e=0;
  int i;
  for(i=1; i<=this->m; ++i){
    e+=2*(this->norm(modes[this->mode2array(i, 0)], modes[this->mode2array(i, 1)]));
  }
  int first=this->firstModeInTail();
  int last=this->lastModeInTail();
  for(i=first; i<=last; ++i){
      e+=2*(this->norm(tail(i, 0), tail(i, 1)));
  }
  //zero mode
  e+=(abs(a0));
  return e;
}

template<typename ScalarT, typename VectorT>
template<typename AVector>
typename AVector::ScalarType D1Real<ScalarT, VectorT>::energy(const AVector& modes)
{
  ScalarType e=0;
  int i;
  for(i=1; i<=this->m; i++){
    e+=2*(power(modes[this->mode2array(i, 0)], 2)+power(modes[this->mode2array(i, 1)], 2));
  }
  return e;
}

template<typename ScalarT, typename VectorT>
template<typename AVector, typename TailType>
typename AVector::ScalarType D1Real<ScalarT, VectorT>::energy(const AVector& modes, const TailType& tail)
{
  ScalarType e=0;
  int i;
  for(i=1; i<=this->m; ++i){
    e+=2*(power(modes[this->mode2array(i, 0)], 2)+power(modes[this->mode2array(i, 1)], 2));
  }

  int first=this->firstModeInTail();
  int last=this->lastModeInTail();
  for(i=first; i<=last; ++i){
    e+=2*(power(tail(i, 0), 2)+power(tail(i, 1), 2));
  }
  return e;
}

template<typename ScalarT, typename VectorT>
template<typename AVector>
typename AVector::ScalarType D1Real<ScalarT, VectorT>::projectionFourierSeries(const ScalarType& x,const AVector& modes)
{
  int i;
  ScalarType result=0;
  int first=1;
  int last=this->m;
  for(i=first; i<=last; ++i){
    result+=ScalarType(2.)*modes[this->mode2array(i, 1)]*ScalarType(capd::intervals::cos(i*x));
    result-=ScalarType(2.)*modes[this->mode2array(i, 0)]*ScalarType(capd::intervals::sin(i*x));
  }

  return result;
}

template<typename ScalarT, typename VectorT>
template<typename AVector>
typename AVector::ScalarType D1Real<ScalarT, VectorT>::closeTailFourierSeries(const ScalarType& x, const AVector& modes)
{
  int i;
  ScalarType result=0;
  int first=this->firstModeInTail();
  int last=this->lastModeInTail();
  for(i=first; i<=last; ++i){
      result+=ScalarType(2.)*modes[this->mode2arrayTail(i, 1)]*capd::intervals::cos(i*x);
      result-=ScalarType(2.)*modes[this->mode2arrayTail(i, 0)]*capd::intervals::sin(i*x);
  }
  return result;
}

}}
#endif
