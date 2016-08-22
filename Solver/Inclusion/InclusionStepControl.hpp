#ifndef INCLUSION_STEP_CONTROL_H_
#define INCLUSION_STEP_CONTROL_H_


#include "capd/auxil/logger.h"
#include "capd/auxil/minmax.h"
#include "capd/basicalg/TypeTraits.h"
#include "capd/rounding/DoubleRounding.h"
#include "capd/basicalg/power.h"
#include "capd/vectalg/Norm.hpp"
#include "capd/vectalg/iobject.hpp"

namespace capd {
namespace dynsys {


class InclusionStepControl {
public:
  typedef double Real;
  /**
   * Constructor
   * @param minStep     minimal allowed step size
   * @param stepFactor  fraction of maximal possible time step to be taken as optimal time step (optimal = stepFactor * maximal)
   */
  InclusionStepControl(Real minStep = 1. / 1048576., Real stepFactor = 0.25) :
    m_minTimeStep(minStep), m_stepFactor(stepFactor) {
  }

  template <class DS, class SetType>
  typename DS::ScalarType computeNextTimeStep(DS& dynamicalSystem, typename DS::ScalarType const& t,const SetType& _x, typename DS::ScalarType maxStep) const {

    typedef typename DS::ScalarType ScalarType;
    typedef typename TypeTraits<ScalarType>::Real Real;
    typename DS::VectorType x(_x);
    Real u = capd::min(1.,sqr(0.025*dynamicalSystem.getOrder())) ;
    ScalarType optStep = dynamicalSystem.getStep()/ u * Real(1.5);
    dynamicalSystem.setStep(optStep);
    while(optStep >= m_minTimeStep) {
      try {
        dynamicalSystem.enclosure(t,x);
        // if succeed then break
        break;
      } catch(...) {
        optStep = (optStep*Real(0.8)).leftBound();
        dynamicalSystem.setStep(optStep);
      }
    }
    double result = toDouble(rightBound(optStep));
    result = clearMantissaBits(result * u);
    result = capd::max(result,m_minTimeStep);
    return capd::min(ScalarType(Real(result)),maxStep);
  }

  template <class DS, class SetType>
  typename DS::ScalarType getFirstTimeStep(DS& dynamicalSystem, typename DS::ScalarType const& t, const SetType& x, const typename DS::ScalarType& maxStep) const {
    return this->computeNextTimeStep(dynamicalSystem,t,x,maxStep);
  }

  Real getMinStepAllowed() const
  {
    return m_minTimeStep;
  }

  Real m_minTimeStep;
  Real m_stepFactor;  // what part of maximal possible step take as optimal one
};


class IEncFoundStepControlForInclusion {
public:
  typedef double Real;
  /**
   * Constructor
   * @param minStep     minimal allowed step size
   * @param stepFactor  fraction of maximal possible time step to be taken as optimal time step (optimal = stepFactor * maximal)
   */
  IEncFoundStepControlForInclusion(Real minStep = 1. / 1048576., Real stepFactor = 0.25) :
    m_minTimeStep(minStep), m_stepFactor(stepFactor) {
  }

  template <class DS, class SetType>
  typename DS::ScalarType computeNextTimeStep(DS& dynamicalSystem, typename DS::ScalarType const& t,const SetType& _x, typename DS::ScalarType maxStep) const {
    typedef typename DS::ScalarType ScalarType;
    typedef typename TypeTraits<ScalarType>::Real Real;
    typename DS::VectorType x(_x);

    
    ScalarType oldStep = dynamicalSystem.getStep();
    std::cout << "Deltha: " << _x.Deltha << std::endl;
    std::cout << "Deltha max: " << capd::vectalg::size( _x.Deltha ) << std::endl;
    
    Real factor;
    if( capd::vectalg::size(_x.Deltha) > 1e-05)
    {
        factor = Real(0.5);
    }
    else
    {
        factor = Real(1.5);
    }

    ScalarType optStep = dynamicalSystem.getStep()/ m_stepFactor * factor;

    dynamicalSystem.setStep(optStep);
    while(optStep >= m_minTimeStep) {
      try {
        dynamicalSystem.enclosure(t,x);
        // if succeed then break
        break;
      } catch(...) {
        optStep = (optStep*Real(0.8)).leftBound();
        dynamicalSystem.setStep(optStep);
      }
    }
    
    
    dynamicalSystem.setStep(oldStep);//TOMEK
    double result = toDouble(rightBound(optStep));
    result = clearMantissaBits(result * m_stepFactor);
    result = capd::max(result,m_minTimeStep);
    return capd::min(ScalarType(Real(result)),maxStep);
  }

  template <class DS, class SetType>
  typename DS::ScalarType getFirstTimeStep(DS& dynamicalSystem, typename DS::ScalarType const& t, const SetType& x, const typename DS::ScalarType& maxStep) const {
    return this->computeNextTimeStep(dynamicalSystem,t,x,maxStep);
  }

  Real getMinStepAllowed() const
  {
    return m_minTimeStep;
  }

  Real m_minTimeStep;
  Real m_stepFactor;  // what part of maximal possible step take as optimal one
};

}}

#endif //INCLUSION_STEP_CONTROL_H_
