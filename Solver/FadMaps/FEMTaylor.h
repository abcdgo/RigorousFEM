#ifndef _FEMTAYLOR_H_
#define _FEMTAYLOR_H_


#include "../DissipativeEnclosure/dissipative_enclosure.hpp"
#include "capd/dynsys/FadTaylor.hpp"


template<class FadMapT, typename StepControlT = capd::dynsys::IEncFoundStepControl , class Taylor = capd::dynsys::FadTaylor<FadMapT,StepControlT> >
class FEMTaylor : public Taylor
{
public:
  typedef FadMapT MapType;
  typedef Taylor BaseClass;
  typedef StepControlT StepControlType;
  typedef typename BaseClass::ScalarType ScalarType;
  typedef typename BaseClass::FScalar FScalar;
  typedef typename BaseClass::TFScalar TFScalar;
  typedef typename BaseClass::MatrixType MatrixType;
  typedef typename BaseClass::VectorType VectorType;
  typedef typename BaseClass::FVector FVector;
  typedef typename BaseClass::TFVector TFVector;
  typedef typename BaseClass::FunctionType FunctionType;

  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;

  FEMTaylor(FadMapT& f, int _order, ScalarType _step, StepControlType stc = StepControlT() ) : Taylor(f, _order, _step, stc)
  {}

  template <typename SetType>
  void operator()(SetType &set)
  {
    set.move(*this);
  }

  ///it is important that we have to recalculate DAG, not only call vectorField's function
  template <typename TimeT>
  inline void setYc(const TimeT &time,const VectorType& y_c){
    this->m_vectorField.setYc(y_c);
    this->m_out = m_vectorField(time,this->m_in); // recording again DAG
  }
    
  ///we do not have to recalculate DAG, because this function is called after moving a set, it is called in order
  ///to get proper enclosures in the next step.
  inline void eraseYc(){
    this->m_vectorField.eraseYc();
  }


  /// the function finds an enclosure for \varphi([0,step],x), using our provided customized algorithm
  template <typename TimeT>
  VectorType enclosure( const TimeT &time, const VectorType& x) const{
std::cout << "FEMTaylor.h l.47 dostaje czas i przesyla go dalej" << std::endl;
    return capd::jaco::enclosure(time, this->m_vectorField, x, this->m_step);
  }

  /// the function finds an enclosure for \varphi([0,step],x), using our provided customized algorithm
  VectorType enclosure(const VectorType& x, std::ostream& out) const{
    return capd::jaco::enclosure(this->m_vectorField, x, this->m_step, out);
  }

  /// the function finds an enclosure for \varphi([0,step],x), using our provided customized algorithm
  VectorType enclosure(const VectorType &x, int* found) const{
    return capd::jaco::enclosure(this->m_vectorField, x, this->m_step);
  }

  /// the function finds enclosure for Jacobian matrix (variational part)
  MatrixType jacEnclosure(const VectorType& enc, const NormType& norm) const{
    return capd::jaco::jacEnclosure(this->m_vectorField,this->m_step,enc,norm);
  }


  template <class SetType>
  ScalarType computeNextTimeStep(const SetType& x, const ScalarType& maxStep) {
    //return this->m_stepControl.computeNextTimeStep(*this,x,maxStep);
    return BaseClass::computeNextTimeStep(*this,x,maxStep);
  }

  template <class SetType>
  ScalarType getFirstTimeStep(const SetType& x, const ScalarType& maxStep) {
    return this->m_stepControl.getFirstTimeStep(*this,x,maxStep);
  }

//  ///has to be OVERRIDDEN in order to make them use provided enclosure finding algorithm
  void encloseMap(const VectorType& x, const VectorType& xx, VectorType& o_phi, MatrixType& o_jacPhi, VectorType& o_rem);

//  ///Computes bound for remainder term in Taylor series
//  VectorType Remainder(const VectorType& v);

}; 

///==========================================function definitions====================================================

//*
//template<typename FadMapT>
template<class FadMapT, typename StepControlT, class Taylor>
void FEMTaylor<FadMapT,StepControlT,Taylor>::encloseMap(
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      MatrixType& o_jacPhi,
      VectorType& o_rem
  )
{
  o_rem = this->Remainder(xx);
  o_jacPhi = this->JacPhi(xx);
  o_phi = this->Phi(x);
}

#endif //_FEMTAYLOR_H_
