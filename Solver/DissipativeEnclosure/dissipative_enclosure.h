#ifndef _DISSIPATIVEENCLOSURE_H_
#define _DISSIPATIVEENCLOSURE_H_

#include "config1.h"
#include "SetBlowUpException.h"

namespace capd {
namespace jaco {

///finds maximal element in a vector
template <typename VectorT>
typename VectorT::ScalarType max(const VectorT& v);

///inflates given interval z in both directions by c.
template <typename FadMapT>
typename FadMapT::ScalarType inflate(typename FadMapT::ScalarType const & z, typename FadMapT::ScalarType::BoundType const& c);

///inflates given interval z in both directions by c, another template.
template <typename ScalarT>
ScalarT inflate(ScalarT const & z, typename ScalarT::BoundType const& c);

///the function finds an C^0 enclosure for \varphi([0,step],x), using combination of both techniques
///-C^1 rough-enclosure,
///-enclosure algorithm based on isolation designed for dPDEs, see [Z3],
///vField is a vector field, its type is the template parameter, should be derivative of FadMap class.
///Data is logged into out.
template <typename FadMapT>
typename FadMapT::VectorType enclosure(typename FadMapT::ScalarType const & time, FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step);

///finds an enclosure for the solution \varphi([0,step],x). It is dedicated to the DPDE maps that take into account the tail
///(provided as the input parameter T). It uses OutputStream as data logger.
template <typename FadMapT>
typename FadMapT::VectorType enclosureInputTail(FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::TailType const & T, typename FadMapT::ScalarType const & step);

///version working with PolyBds
template <typename EquationT>
typename EquationT::PolyBdType enclosure(EquationT & vField, typename EquationT::PolyBdType const & x,
    typename EquationT::RealType const & step);


///=================08.2011 Older versions, left for compability reasons===================


///the function finds an C^0 enclosure for \varphi([0,step],x), using combination of both techniques
///-C^1 rough-enclosure,
///-enclosure algorithm based on isolation designed for dPDEs, see [Z3],
///vField is a vector field, its type is the template parameter, should be derivative of FadMap class.
///Data is logged into out.
template <typename FadMapT>
typename FadMapT::VectorType enclosure(FadMapT & vField, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step, std::ostream& out);


///like regular enclosure function, but calculates on provided tail, FadMapT should provide functions
///setT0() and getT0().
template <typename FadMapT, typename Tail>
typename FadMapT::VectorType enclosure( typename FadMapT::ScalarType const & currentTime, FadMapT & vField, const Tail& tail, typename FadMapT::MatrixType::RowVectorType const & x,
    typename FadMapT::ScalarType const & step );

///the function finds enclosure for Jacobian matrix (variational part)
///source- "C^1-Lohner algorithm" by P. Zgliczynski.
///Extract from CAPD library.
template <typename FadMapT, typename NormType>
typename FadMapT::MatrixType jacEnclosure(const FadMapT& vectorField, const typename FadMapT::ScalarType& step,
    const typename FadMapT::VectorType &enc, const NormType &the_norm, typename FadMapT::ScalarType* o_logNormOfDerivative = 0);

}
}

#endif // _DISSIPATIVEENCLOSURE_H_ 
