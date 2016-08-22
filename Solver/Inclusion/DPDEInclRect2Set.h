#ifndef _CAPD_JACO_JACOINCLRECT2SET_H
#define _CAPD_JACO_JACOINCLRECT2SET_H

#include "capd/diffIncl/InclRect2Set.hpp"

namespace capd{
namespace jaco{


template<typename MatrixT>
class DPDEInclRect2Set : public capd::diffIncl::InclRect2Set<MatrixT>{
  
public: 
  typedef MatrixT MatrixType;
  typedef capd::diffIncl::InclRect2Set<MatrixType> BaseSet;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;

  explicit DPDEInclRect2Set(int dimension);
  explicit DPDEInclRect2Set(const VectorType& the_x);
  DPDEInclRect2Set(const VectorType& the_x, const VectorType& the_r0);
  DPDEInclRect2Set(const VectorType& the_x, const MatrixType& the_C, const VectorType& the_r0);
  DPDEInclRect2Set(const VectorType& the_x, const MatrixType& the_C,
               const VectorType& the_r0,
               const VectorType& the_r
  );

 
  template<typename DiffIncl>
  void move( DiffIncl& dynsys, int stepNumber, std::ostream& f);

  template<typename DiffIncl>
  void move( DiffIncl & dynsys, DPDEInclRect2Set& result) const;

  VectorType Deltha;

  using BaseSet::get_x;
  using BaseSet::get_r;
  using BaseSet::get_r0;
  using BaseSet::get_B;
  using BaseSet::get_C;
  using BaseSet::operator VectorType;
  using BaseSet::show;
  using BaseSet::affineTransformation;

protected:
  using BaseSet::m_x;
  using BaseSet::m_r;
  using BaseSet::m_r0;
  using BaseSet::m_B;
  using BaseSet::m_C;
  
};

///==========================================function definitions====================================================

template<typename MatrixType>
inline DPDEInclRect2Set<MatrixType>::DPDEInclRect2Set(int dim)
  :  BaseSet(dim) {
}

template<typename MatrixType>
inline DPDEInclRect2Set<MatrixType>::DPDEInclRect2Set(const VectorType& the_x)
  :  BaseSet(the_x) {
}

template<typename MatrixType>
inline DPDEInclRect2Set<MatrixType>::DPDEInclRect2Set(const VectorType& the_x,const VectorType& the_r0)
  :  BaseSet(the_x, the_r0) {
}

template<typename MatrixType>
inline DPDEInclRect2Set<MatrixType>::DPDEInclRect2Set(
      const VectorType& the_x,
      const MatrixType& the_C,
      const VectorType& the_r0
   )
  : BaseSet(the_x, the_C, the_r0){
}

template<typename MatrixType>
inline DPDEInclRect2Set<MatrixType>::DPDEInclRect2Set(
      const VectorType& the_x,
      const MatrixType &the_C,
      const VectorType& the_r0,
      const VectorType& the_r
   ): BaseSet(the_x, the_C, the_r0, the_r){
}

template<typename MatrixType>
template<typename DiffIncl>
void DPDEInclRect2Set<MatrixType>::move(DiffIncl & diffIncl, int stepNumber, std::ostream& f) {

  //std::cout << "diffIncl.getDynamicalSystem().nextStep=" << diffIncl.getDynamicalSystem().computeNextTimeStep(*this,1.e-1) << "\n";

  ///opening a file for debugging purposes

  typedef typename DiffIncl::TailType TailType;
  typedef typename DiffIncl::MultiMapType MultiMapType;

  VectorType x = VectorType(*this);
  VectorType x0=x;

  ///variables used for storing
  TailType T=TailType(); /// T: T([0, h]) \subset T
  TailType N=TailType(); /// N: N(x+T) \subset N
  VectorType W_2; ///[W_2]
  //f<<"Step number "<<stepNumber<<" in progress.\n";
  VectorType y_c;
  {///steps 1-3, 5-8 of Algorithm 1
  Deltha = diffIncl.perturbations(this->getCurrentTime(), x, T, N, W_2, y_c, f);
  }///end steps 1-3, 5-8 of Algorithm 1
  ///in T we have T([0,h])
  ///in N we have N([W_2], T[0, h])

  {///step  4 of Algorithm 1
//  std::cout << "DPDEInclRect2Set, setYc time: " << this->getCurrentTime() << std::endl;
  diffIncl.getDynamicalSystem().setYc(this->getCurrentTime(),y_c); ///setting y_c of a unperturbed projection
  BaseSet::BaseSet::move(diffIncl .getDynamicalSystem()); ///computation of an unperturbed trajectory
  diffIncl.getDynamicalSystem().eraseYc(); ///y_c=0, in order to make the algorithm produce good enclosures in the
                                           ///next step
  }///end step 4 of Algorithm 1

  {///step 10 of Algorithm 1. Rearrangements
  x = midVector( m_x + Deltha);
  VectorType dr = m_x + Deltha - x;
  MatrixType bt = Transpose(m_B);
  m_r = m_r + bt * dr;
  m_x = x;
  }///end step 10 of Algorithm 1

  TailType Th;
  {///step 11 of Algorithm 1
  Th=diffIncl.calculateTh(x, T, N, W_2, f); ///calculates T(h)
  }///end step 11 of Algorithm 1

  #if __GENERAL_DEBUG__
    f<<"current x: \n"<<m_x<<"\n";
    f<<"current r: \n"<<m_r<<"\n";
    f<<"current r0: \n"<<m_r0<<"\n";
    f<<"calculated T(h):\n";
    Th.print(f);
  #endif
  ///we have T(h) is T(0) for the next step
  diffIncl.setT0(Th);
}

template<typename MatrixType>
template<typename DiffIncl>
void DPDEInclRect2Set<MatrixType>::move(DiffIncl & diffIncl, DPDEInclRect2Set<MatrixType>& result) const {
  std::ofstream f;
#if __GENERAL_DEBUG__
  f.open(__FILE_NAME_DEBUG__, std::ofstream::app);
#else
   f=std::cout;
#endif
   f<<"Wrong move function call in DPDEInclRect2Set class. Use move(diffIncl, stepNumber) instead.\n";
   std::cerr<<"Wrong move function call in DPDEInclRect2Set class. Use move(diffIncl, stepNumber) instead.\n";
   throw std::runtime_error("Wrong move function call in DPDEInclRect2Set class. Use move(diffIncl, stepNumber) instead.\n");
#if __GENERAL_DEBUG__
  f.close();
#endif
}

}}


#endif
