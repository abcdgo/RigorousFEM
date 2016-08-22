#ifndef _CAPD_JACO_DPDEINCLUSIONCW_H_
#define _CAPD_JACO_DPDEINCLUSIONCW_H_

#include "capd/diffIncl/DiffInclusion.hpp"
#include "capd/diffIncl/DiffInclusionCW.hpp"
#include "capd/diffIncl/InclRect2Set.hpp"
#include "tail.h"
#include <time.h>

namespace capd{
namespace jaco{

///provides functions calculating all sets needed to move an instance of DPDEInclRect2Set.h
///enclosures [W_1], [W_2] and T,  interval set [\Delta], and the selection y_c, steps from first
///to eight of Algorithm 1.
template<typename MapT, typename Tail, typename DynSysT >
class DPDEInclusionCW : public capd::diffIncl::DiffInclusionCW<MapT, DynSysT>{
   
  const char* m_fileName;
  bool m_debug; 
 
  public:
  typedef capd::diffIncl::DiffInclusionCW<MapT, DynSysT> BaseClass;
  typedef typename BaseClass::MultiMapType   MultiMapType;
  typedef typename MultiMapType::TailType    TailType;
  typedef typename BaseClass::MapType        MapType;
  typedef typename BaseClass::FunctionType   FunctionType;
  typedef typename BaseClass::MatrixType     MatrixType;
  typedef typename BaseClass::VectorType     VectorType;
  typedef typename BaseClass::ScalarType     ScalarType;
  typedef typename BaseClass::NormType       NormType;

  DPDEInclusionCW( MultiMapType & diffIncl, 
      int order,
      ScalarType const & step,
      NormType const & norm
  );

  //Inclusion storing report in file "fileName"
  DPDEInclusionCW( MultiMapType & diffIncl, 
      int order,
      ScalarType const & step,
      NormType const & norm, 
      const char* fileName
  );

  ///calculates enclosure, [W_2] in the description of Algorithm 1, using provided enclosure (candidate)
  ///for the tail T.
  template <typename TimeT>
  VectorType diffInclusionEnclosure(const TimeT& t, const VectorType& x, const Tail& T) const;

  ///calculates enclosure, [W_2] in the description of Algorithm 1, using the tail saved in a DPDEMultiMap
  ///instance.
  template <typename TimeT>
  VectorType diffInclusionEnclosure(const TimeT& t, const VectorType& x) const;

  ///calculates enclosures, [W_2] and T: T([0,h])\subset T in the description of Algorithm 1.
  ///. Also the nonlinear part N_k(x+T) calculated for the far modes $|k|>m$ is returned as a PolyBd instance N.
  ///Output is written (depending on the debugging flags provided) to the outputs out and out2, see the config.h file.
  template <typename TimeT>
  VectorType diffInclusionEnclosure(const TimeT& t, const VectorType& x, Tail& tail, Tail& N, std::ostream& out) const;

  ///enclosure as above, but with output written to an output out.
  template <typename TimeT>
  VectorType diffInclusionEnclosure(const TimeT& t, const VectorType& x, std::ostream& out) const;

  ///calculates the enclosure [W_1], \varphi([0,h], [x_0], y_c)\subset [W_1].
  template <typename TimeT>
  VectorType dynamicalSystemEnclosure(const TimeT& t, const VectorType& x, std::ostream& out) const;

  ///alias for the diffInclusionEnclosure.
  template <typename TimeT>
  VectorType enclosure(const TimeT& t, const VectorType & x) const;

  ///calculates the [\Delta] in the step 8 of Algorithm 1.
  ///Input
  ///-x, as in the description of Algorithm 1,
  ///-T, the enclosure (candidate) for the T([0,h]),
  ///Output
  ///-Nt, is N_k(x+T),
  ///-W_2, the enclosure [W_2] in the step 1 of Algorithm 1,
  ///-y_c, y_c in the step 3 of Algorithm 1.
  VectorType perturbations(const ScalarType &time, const VectorType & x, Tail& T, Tail& Nt, VectorType& W_2, VectorType& y_c, std::ostream& out);

  inline void setT0(const Tail& T){
    m_diffIncl.setT0(T);
  }

  inline const TailType& getT0() const{
	  return m_diffIncl.getT0();
  }

  ///calculates T(h), last step in Algorithm 1.
  Tail calculateTh(const VectorType& x, const Tail& T, const Tail& N, const VectorType& W_2, std::ostream& out);

  MultiMapType& getMultiMap() ;

  using BaseClass::getStep;


protected:
  using BaseClass::m_norm;
  using BaseClass::m_diffIncl;
  using BaseClass::getCurrentTime;
  
};

///==========================================function definitions====================================================

template <typename MapT, typename Tail, typename DynSysT>
template <typename TimeT>
inline typename DPDEInclusionCW<MapT, Tail, DynSysT>::VectorType DPDEInclusionCW<MapT, Tail, DynSysT>::enclosure
(const TimeT& t, const DPDEInclusionCW<MapT, Tail, DynSysT>::VectorType& x) const{
    return capd::dynsys::enclosure(this->m_diffIncl,t, x, this->getStep());
}

template <typename MapT, typename TailT, typename DynSysT>
DPDEInclusionCW<MapT, TailT, DynSysT>::DPDEInclusionCW(
    MultiMapType& A_diffIncl,
    int A_order,
    ScalarType const & A_step,
    NormType const & A_norm
)
: BaseClass(A_diffIncl, A_order, A_step, A_norm){
    m_debug=false;
}

template <typename MapT, typename TailT, typename DynSysT>
DPDEInclusionCW<MapT, TailT, DynSysT>::DPDEInclusionCW(
    MultiMapType& A_diffIncl,
    int A_order,
    ScalarType const & A_step,
    NormType const & A_norm,
    const char* nazwaPliku
)
: BaseClass(A_diffIncl, A_order, A_step, A_norm){
    m_fileName=nazwaPliku;
    m_debug=true;
}

template <typename MapT, typename TailT, typename DynSysT>
template <typename TimeT>
typename DPDEInclusionCW<MapT, TailT, DynSysT>::VectorType DPDEInclusionCW<MapT, TailT, DynSysT>::diffInclusionEnclosure(const TimeT& t, const VectorType& x, const TailT& T) const{
    return capd::jaco::enclosure( t, this->m_diffIncl, T, x, this->getStep() );
}

template <typename MapT, typename TailT, typename DynSysT>
template <typename TimeT>
inline typename DPDEInclusionCW<MapT, TailT, DynSysT>::VectorType DPDEInclusionCW<MapT, TailT, DynSysT>::diffInclusionEnclosure(const TimeT& t, const VectorType& x) const{
std::cout << "DPDEInclusionCW dodany time ale nie uzyty 4" << std::endl;
  return capd::jaco::enclosure(m_diffIncl, x, getStep());
}

template <typename MapT, typename TailT, typename DynSysT>
template <typename TimeT>
inline typename DPDEInclusionCW<MapT, TailT, DynSysT>::VectorType DPDEInclusionCW<MapT, TailT, DynSysT>::diffInclusionEnclosure(const TimeT& t, const VectorType& x, std::ostream& out) const{
std::cout << "DPDEInclusionCW dodany time ale nie uzyty 3" << std::endl;
  return capd::jaco::enclosure(m_diffIncl, x, getStep(), out);
}

template <typename MapT, typename TailT, typename DynSysT>
template <typename TimeT>
inline typename DPDEInclusionCW<MapT, TailT, DynSysT>::VectorType DPDEInclusionCW<MapT, TailT, DynSysT>::dynamicalSystemEnclosure(const TimeT& t, const VectorType& x, std::ostream& out) const{
std::cout << "DPDEInclusionCW dodany time ale nie uzyty 2" << std::endl;
  return this->m_dynamicalSystem.enclosure(x, out);
}

template<typename MapT, typename TailT, typename DynSysT>
template< typename TimeT>
typename DPDEInclusionCW<MapT, TailT, DynSysT>::VectorType DPDEInclusionCW<MapT, TailT, DynSysT>::diffInclusionEnclosure
(const TimeT& time, const DPDEInclusionCW<MapT, TailT, DynSysT>::VectorType& x, TailT& T, TailT& N, std::ostream& out) const
{
    MultiMapType& mm=this->m_diffIncl;
    typename MultiMapType::PerturbationType* p=mm.getPerturbation();
    TailT T0=p->getT0();
    VectorType W2=p->template enclosureWithTail< TimeT, typename DPDEInclusionCW<MapT, TailT, DynSysT>::VectorType,
        DPDEInclusionCW<MapT, TailT, DynSysT> >(time, this, x, T0, T, N, this->getStep(), out );
    return W2;
}

template <typename MapT, typename TailT, typename DynSysT>
typename DPDEInclusionCW<MapT, TailT, DynSysT>::VectorType DPDEInclusionCW<MapT, TailT, DynSysT>::perturbations(const ScalarType &time, const VectorType& x, TailT& T, TailT& N, VectorType& W_2, VectorType& y_c, std::ostream& f){
  #if __GENERAL_DEBUG__
    f<<"current vector x:\n";
    m_diffIncl.printModes(x, f);
    f<<"current tail:\n";
    getT0().print(f);
  #endif
  clock_t start, end;
  start=clock();
  ///in blocks {} we collect all the steps of Algorithm 1, with reference to corresponding step number
  {///step 1
//  std::cout << "DPDEInclusionCW l.197 time: " << time << std::endl;
  W_2 = diffInclusionEnclosure( time, x, T, N, f); ///returns enclosures [W_2], T: T([0,h])\subset T, N: N_k(x+T)
  }///end of step 1
  end=clock(); ///measures runtime of one step of the diffInclusionEnclosure function
//  std::cout<<"diffInclEncl time="<<(double)(end-start)/CLOCKS_PER_SEC<<"\n";
  #if __GENERAL_DEBUG__
    f<<"current T:\n";
    T.print(f);
  #endif

  #if __GENERAL_DEBUG__
    f<<"found W_2:\n";
    m_diffIncl.printModes(W_2, f);
  #endif

  VectorType galerkinProjectionError, delta, C, D, result;
  {///step 2
  //std::cout << "DPDEInclusionCW l.214 perturbations(), czas z nikad nie przekazywany tylko uzyty getCurrentTime(): " << time << std::endl;
  galerkinProjectionError=m_diffIncl.perturbations(time, W_2, T);
  }///end of step 2

  {///step 3
  y_c=midVector(galerkinProjectionError);
  }///end of step 3

  {///step 5
  delta=y_c-galerkinProjectionError;
  }///end of step 5

  #if __GENERAL_DEBUG__
    f<<"found mid(delta):\n";
    m_diffIncl.printModes(y_c, f);
    f<<"found delta:\n";
    m_diffIncl.printModesWithSymmetryCheck(delta, f);
  #endif

//==================================COMPONENT WISE ESTIMATES=================================
  {///step 6
  C = rightVector(abs(delta));
  }///end of step 6

  int i, j;
  MatrixType J;
  {///step 7
  J = (m_diffIncl.getVectorField()).derivative(time,W_2); /// calculating Jacobian on W_2\times y_c
  for( i=0; i < J.numberOfRows(); ++i){
    for( j = 0; j < J.numberOfColumns(); ++j){
      if(i != j)
        J[i][j] = (capd::abs(J[i][j])).right();
      else
        J[i][j] = (J[i][j]).right();
    }
  }
  }///end of step 7


  {///step 8
  MatrixType At = getStep() * J;
  MatrixType A = MatrixType::Identity(J.numberOfRows());
  MatrixType Sum = A;
  ScalarType AtNorm = right((*m_norm)(At)),
      AnNorm = right((*m_norm)(A));
  ScalarType n = 2.0;  // n = i + 2
  // remainder = |A| * |At/(N+2)| / (1 - |At/(N+2)|)   (the sum of geometric series from N to infinity)
  ScalarType q = AtNorm/n;

  while(true){
    A = A * At / n;
    Sum += A;
    AnNorm *= q;
    n += ScalarType(1.0);
    q = AtNorm / n;
    if(q < 1){
      // remainder = |A| * |At/(N+2)| / (1 - |At/(N+2)|)   (the sum of geometric series from N to infinity)
      ScalarType remainder = right(AnNorm * AtNorm / (n - AtNorm ));
      if(remainder < this->m_errorTolerance)
        break;
    }
  }
//  std::cout << "\n N = " << n-2 << "\n remainder approx : " << remainder;
  // we recompute remainder because norm of A can be smaller than approximation : AnNorm
  ScalarType remainder = right((*m_norm)(A) * AtNorm / (n  - AtNorm ));
//  std::cout << "\n remainder real   : " << remainder;

  for(i=0; i < J.numberOfRows(); ++i)
    for(int j = 0; j < J.numberOfColumns(); ++j)
      Sum[i][j] += remainder * ScalarType(-1.0, 1.0);
  VectorType D = getStep() * (Sum * C);
  result=VectorType(D.dimension());

  for(int i=0; i< D.dimension(); ++i)
    result[i] = ScalarType(-D[i].rightBound(), D[i].rightBound());

  }///end of step 8

//=different estimates, based on LOGARITHMIC NORM============================================
//  ScalarType C = right((*m_norm)(delta));
//  ScalarType l = right((*m_norm)(J));
//
//  ScalarType D = (l.contains(0.0))? C*getStep() : (C*(exp(l*getStep())-1))/l;
//  VectorType result(x.dimension());
//
//  for(int i=0; i< result.dimension(); ++i)
//    result[i] = ScalarType(-D.rightBound(), D.rightBound());
//==========================================================================================
  #if __GENERAL_DEBUG__
    f<<"found Deltha:\n";
    m_diffIncl.printModes(result, f);
  #endif

  return result;
}

template <typename MapT, typename TailT, typename DynSysT>
TailT DPDEInclusionCW<MapT, TailT, DynSysT>::calculateTh(const VectorType& vec, const TailT& T, const TailT& N, const VectorType& W_2, std::ostream& f){
  //for close TailT we calculate explicitly
  MultiMapType& mm=this->m_diffIncl;
  typename MultiMapType::PerturbationType* p=mm.getPerturbation();
  TailT T0=p->getT0();
  TailT Th=p->template calculateTh<ScalarType, VectorType>(vec, T0, T, N, W_2, this->getStep(), f);
  return Th;
}

template <typename MapT, typename TailT, typename DynSysT>
typename DPDEInclusionCW<MapT, TailT, DynSysT>::MultiMapType& DPDEInclusionCW<MapT, TailT, DynSysT>::getMultiMap(){
    return this->m_diffIncl;
}

}}
#endif
