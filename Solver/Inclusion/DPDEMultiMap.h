#ifndef _CAPD_JACO_JACOMULTIMAP_H_ 
#define _CAPD_JACO_JACOMULTIMAP_H_ 

#include "config.h"

#include "capd/dynsys/FadTaylor.h"
#include "capd/diffIncl/MultiMap.h"

namespace capd{
namespace jaco{

///basically P_mF(x+T), is a composed class, consisting of
///- P_mF(x) being an instance of FMapT,
///- P_mF(x+T)-P_mF(x) being an instance of GMapT.
template<class FMapT, class GMapT>
class DPDEMultiMap : public capd::diffIncl::MultiMap< FMapT, GMapT >, public GMapT::Class{
  public:
  typedef FMapT MapType;
  typedef GMapT PerturbationType;
  typedef typename GMapT::TailType TailType;
  typedef typename GMapT::ScalarType ScalarType;
  typedef typename GMapT::VectorType VectorType;
  typedef typename MapType::MatrixType MatrixType;
  typedef typename FMapT::DoubleType DoubleType;

  DPDEMultiMap( MapType & Field,
      PerturbationType & Perturb
  ) : capd::diffIncl::MultiMap< FMapT, GMapT >(Field, Perturb), GMapT::Class(Perturb){}


  inline const TailType& getT0() const{
     return this->m_g->getT0();
  }

  inline void setT0(const TailType& T){
    this->m_g->setT0(T);
  }         

  ///sets value of the zero mode a_0, being a constant real number.
  void setA0(const ScalarType& a0){
	  this->m_f->setA0(a0);
  }

  ///gets value of the zero mode a_0, being a constant real number.
  ScalarType getA0() const{
	  return this->m_f->getA0();
  }

  ///the linear part of the vector field.
  template <typename TimeT>
  inline VectorType L(const TimeT &t, const VectorType& z){
    return this->m_f->L(t,z);
  }

  //returns P_mN(x+T).
  template <typename TimeT>
  inline VectorType N(const TimeT &t, const VectorType& u){
    return this->m_g->N(t,u)+this->m_f->N(t,u);
  }

  ///returns P_mN(x) and writes to the output out.
  template <typename TimeT>
  inline VectorType N(const TimeT &t, const VectorType& u, std::ostream& out){
    return this->m_g->N(t, u, out)+this->m_f->N(t, u);
  }

  ///returns P_mN(x+T), using the provided T (tail).
  template <typename TimeT>
  inline VectorType N(const TimeT &t, const VectorType& u, const TailType& tail){
    return this->m_g->N(t, u, tail);
  }

  ///returns Q_mN(x+T) and writes to the standard output. It returns a PolyBd instance consisting values of
  ///N_k(x+T) for |k|>=m.
  template <typename TimeT>
  inline TailType n(const TimeT &t, const VectorType& u){
    return this->m_g->n(t, u, std::cout);
  }

  ///returns Q_mN(x+T) and writes to the standard output. It returns a PolyBd instance consisting values of
  ///N_k(x+T) for |k|>=m. Uses the provided T.
  template <typename TimeT>
  inline TailType n(const TimeT &t, const VectorType& u, const TailType& T){
    return this->m_g->n(t, u, T);
  }  

  ///true in case that vector field corresponding to a_k can be considered dissipative.
  inline bool isDissipative(int k){
    return this->m_f->isDissipative(k);
  }

  ///returns \lambda_k
  inline ScalarType lambda(int k){
    return this->m_f->lambda(k);
  }

  ///returns leftBound(\ni) because it is the smallest viscosity constant, and nonrigorous calculations are performed using this
  ///value, which is somehow the worst from the whole interval. Returns -\\ni because of representation.
  inline DoubleType getMiddleNi() const{
	  return rightBound(this->m_f->m_ni);
  }

  ///returns the vector field P_mF(x+T), using the provided tail.
  template <typename TimeT>
  inline VectorType operator()(const TimeT &t, const VectorType & X, const TailType& T)
  {
    VectorType result= (*this->m_f)(t, X) + (*this->m_g)(t, X, T);
    return result;
  }

  ///returns the vector field P_mF(x+T).
  template <typename TimeT>
  VectorType operator()( const TimeT &timeRange, const VectorType & X) const {
    return (*this->m_f)(timeRange,X) + (*this->m_g)(timeRange,X);
  }//*/

  ///returns the Galerkin projection error P_mF(x+T)-P_mF(x), using the provided tail.
  template <typename TimeT>
  VectorType perturbations(const TimeT &time, VectorType const & X, const TailType& T)
  {
    return this->m_g->perturbations(time, X, T);
  }

  ///returns the Galerkin projection error P_mF(x+T)-P_mF(x).
  template <typename TimeT>
  VectorType perturbations(const TimeT &time, VectorType const & X)
  {
    return this->m_g->perturbations(time, X);
  }

  inline PerturbationType* getPerturbation() const{
    return this->m_g;
  }

  ///the Jacobian of the linear part.
  inline MatrixType& jacobianL(){
	  return this->m_f->jacobianL();
  }

  /// derivative  of a MultiMap calculated on a set X\times T
  MatrixType calculateJ(const VectorType & X, const TailType& T) const;

  ///function calculating rigorous value (taking care of the infinite dimension )
  ///of the solution  at point x.
  //VectorType calculateFourierSeries(const VectorType& x, const VectorType& modes) const;

  void debugMsg(std::ostream& out){
    out<<"MultiMap\n"; this->getT0().print(out);
  }


    std::string array2mode(int i) {return "DPDEMultiMap.h array2mode";}
    std::string array2modeIndex(int i) {return "DPDEMultiMap.h array2modeIndex";}
    std::string ni(int i) {return "DPDEMultiMap.h Knights of Ni";}

}; 

///==========================================function definitions====================================================


template<class FMapT, class GMapT>
typename DPDEMultiMap<FMapT, GMapT>::MatrixType DPDEMultiMap<FMapT, GMapT>::calculateJ(const VectorType & X, const TailType& T) const {
  MatrixType result= (*this->m_f)[X] + this->m_g->jacobian(X, T);
  return result;
}


}}

#endif // _CAPD_DYNSYS_FADTAYLOR_H_ 
