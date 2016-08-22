#ifndef _CAPD_ALGORITHMS2_H_
#define _CAPD_ALGORITHMS2_H_

#include "../FadMaps/Inclusion/Data.hpp"

namespace capd {
namespace jaco {

///Class providing algorithms used by maps. All algorithms involving calculation of the convolution sum use the symmetries.
template <class TailT, class PDET>
class Algorithms2 : public PDET {

public:
  typedef typename capd::jaco::Algorithms2<TailT, PDET> Class;
  typedef PDET BaseClass;
  typedef TailT TailType;
  typedef typename BaseClass::ScalarType ScalarType;

  TailType m_T;
  Data<ScalarType> &data;

  Algorithms2(int m, int M, ScalarType ni, ScalarType& a0, ScalarType s , Data<ScalarType> &d) : BaseClass(m, M, ni, a0), data(d)
  {
    m_T = this->initializeTail(_C, s - 1, this->m, this->M);
  }

  Algorithms2(ScalarType ni) : BaseClass(ni, 0) {}

  Algorithms2(int m, int M) : BaseClass(m, M, _ni, 0) {}

  Algorithms2() : BaseClass(_ni, 0) {}

  Algorithms2(const Algorithms2& a2) : PDET(a2.m_ni, a2.m_a0), m_T(a2.m_T), data(a2.data) {
    this->m = a2.m;
    this->M = a2.M;
  }

  ///Tabularizes exponents.
  //void calculateExponents(const ScalarType& h);

  ///Calculates b_k = N_k / -\lambda_k, for k: m<|k|<=M (near tail),
  ///if re==true then return Re{b_k},
  ///else return Im{a_k}.
  ScalarType b(const ScalarType& h, int k, bool re, const TailType& N) const;

  ///Calculates b_k = N_k / -\lambda_k, for k: m<|k|<=M (near tail), explicit index version, where i is index in the modes storage.
  ///Returns either Re{a_k} or Im{a_k}, depending on which value is at index i in the modes storage.
  ///Parameter h is not used.
  ScalarType b(const ScalarType& h, int i, const TailType& N) const { return ScalarType(N[i] / (-this->lambdaTail(i))); }

  ///Calculates b_k = N_k / -\lambda_k for k: m<|k| (tail), result is returned as a PolyBd instance.
  ///Parameter h is not used.
  TailType bt(const ScalarType& h, const TailType& N) const;

  ///Calculates g_k^{+-} = (T0_k^{+-}-b_k^{+-})e^{-\lambda_k h}+b_k^{+-},explicit index version, where i is
  ///index in modes storage.
  ///Returns either Re{a_k} or Im{a_k}, depending on which value is at index i in the modes storage.
  ScalarType g(const ScalarType& h, int i, const TailType& T0, const TailType& N) const;

  ///Calculates g_k^{+-} = (T0_k^{+-}-b_k^{+-})e^{-\lambda_k h}+b_k^{+-}, version used by enclosure algorithms.
  ///Returns either Re{a_k} or Im{a_k}, depending on which value is at index i in the modes storage.
  ScalarType g_encl(const ScalarType& h, int i, const ScalarType& x, const ScalarType& b);

  ///Functions calculating values D_1 and D_2, bounding |FS_1{k}| and |FS_2{k}|, in the following sense
  ///|FS_1{k}| <= D_1 / |k|^s  and  |FS_2{k}| <= D_2 / |k|^s.
  /// Needed to bound far part of |N_k| (for k: |k|>M). See Section 10 in the paper.

  ///Calculates D_1 in |FS_1{k}| <= D_1 / |k|^s , bounding the finite part of convolution sum for k: |k| <= 2M.
  template <typename VectorType>
  ScalarType D1(const VectorType& u, const TailType& tail) const;

  ///Calculates D_2 in |FS_2{k}| <= D_2 / |k|^s , bounding the finite part of convolution sum for k: |k| > 2M.
  template <typename VectorType>
  ScalarType D2(const VectorType& x, const TailType& T) const;

  ///Calculates D_\infty in |IS{k}| <= D_\infty / |k|^s, bounding the infinite part of convolution sum for k: |k| > M.
  template <typename VectorType>
  ScalarType Dinf(const VectorType& x, const TailType& T) const;

  ///Calculates product of two complex modes and the imaginary unit (belonging either to asymmetric Galerkin projection
  ///or a tail) i\cdot a_k\cdot a_{k-k_1} = prod, stores the result in variables
  ///-in re stores Re{prod},
  ///-in im stores Im{prod}.
  template <typename VectorType>
  inline void product(int k, int k_1, typename VectorType::ScalarType& re, typename VectorType::ScalarType& im, const VectorType& u,
      const TailType& tail, double coeff = 1) const;

  ///Calculates product of two complex modes and the imaginary unit (belonging ONLY to a symmetric Galerkin projection)
  ///i\cdot a_k\cdot a_{k-k_1} = prod, stores the result in variables
  ///-in re stores Re{prod},
  ///-in im stores Im{prod}.
  template <typename VectorType>
  inline void product(int k, int k_1, typename VectorType::ScalarType& re, typename VectorType::ScalarType& im, const VectorType& u, double coeff = 1) const;

  ///Calculates P_mF(x).
  template <typename VectorType>
  VectorType vectorField(const VectorType& u) const;

  ///calculates P_mN(x), essentially it is the truncated convolution sum, consisting only modes from the symmetric
  ///Galerkin projection.
  template <typename TimeT, typename VectorType>
  VectorType N(const TimeT &time, const VectorType& z) const;

  ///Calculates P_mN(x+T)-P_mN(x), essentially the Galerkin projection error.
  template <typename TimeT, typename VectorType>
  VectorType N(const TimeT &time, const VectorType& u, const TailType& tail) const;

  template <typename TimeT, typename VectorType>
  VectorType perturbations(const TimeT &time, const VectorType& u, const TailType& tail) const;


  ///Calculates N_k(x+T) for k: |k|>M (corresponding to the mode in the far tail).
  //template <typename VectorType>
  //ScalarType n(int k, int real, const VectorType& u, const TailType& tail, std::ostream& f) const;

  ///Calculates Q_mN(x+T), returns {N_k(x+T), |k|>m} as an instance of PolyBd structure (a tail).
  template <typename TimeT, typename VectorType>
  TailType n(const TimeT &time, const VectorType& u, const TailType& tail) const;

  TailType initializeTail(const ScalarType& c, const ScalarType& s, int m, int M) {
    TailType T(0, c, s, m, M); return T;
  }

  template <typename AScalar>
  AScalar inflate(AScalar const & x, double c) const;


        template <typename TimeT, typename VectorType>
        VectorType nonLinear(const TimeT &time, const VectorType& z) const;

        ScalarType firstIntegral(int j, int k, int i) const;
        ScalarType secondIntegral(int j, int k, int i) const;

        ScalarType lambda( int i ) const
        {
            return data.lambda(i);
        }

        ScalarType lambda_k( int i ) const
        {
            return lambda(i);
        }

        ScalarType lambdaTail( int i ) const
        {
            return lambda(i);
        }

        ScalarType V(int K) const
        {
            return leftBound(this->ni(K));
        }
        ScalarType ni(int K) const
        {
            return ScalarType(1);
        }
        bool isDissipative(int i)
        {
            return data.isDissipative(i);
        }

};

///==========================================function definitions====================================================

template <class TailT, class PDET>
typename Algorithms2<TailT, PDET>::ScalarType Algorithms2<TailT, PDET>::b(const ScalarType& h, int k, bool re, const TailType& N) const {
  if(!N.inTail(k))
    throw std::runtime_error("index that is not in the tail requested in the function b()");
  return ScalarType(N(k, re) / (-this->lambda_k(k)));
}

template <class TailT, class PDET>
TailT Algorithms2<TailT, PDET>::bt(const ScalarType& h, const TailType& N) const {
  TailType bt(N);
  int k;
  int first = this->firstModeInTail();
  int last = this->lastModeInTail();
  //std::cout << "first: " << first << std::endl;
  //std::cout << "last: " << last << std::endl;
  for(k = first; k <= last; ++k) {
    if(N.inTail(k)) {
      bt.set(k, 1, b(h, k, 1, N));
      bt.set(k, 0, b(h, k, 0, N));
    }
  }
  //far tail
  setC(bt, C(N) / V(this->M + 1));

  setS(bt, s(N) + this->m_p);
  return bt;
}

template <class TailT, class PDET>
typename Algorithms2<TailT, PDET>::ScalarType Algorithms2<TailT, PDET>::g(const ScalarType& h, int i, const TailType& T0, const TailType& N) const {
  ScalarType g;
  ScalarType bk = b(h, i, N);

  ScalarType ex = exp(h * this->lambdaTail(i));
  g.setLeftBound((((T0[i].leftBound() - bk.leftBound()) * ex) + bk.leftBound()).leftBound());
  g.setRightBound((((T0[i].rightBound() - bk.rightBound()) * ex) + bk.rightBound()).rightBound());
  return g;
}

template <class TailT, class PDET>
typename Algorithms2<TailT, PDET>::ScalarType Algorithms2<TailT, PDET>::g_encl(const ScalarType& h, int i, const ScalarType& x, const ScalarType& b) {
  ScalarType g;

  ///NOT USING TABULARIZED EXPONENTS (SLIGHTLY SLOWER)
  ScalarType ex=exp(h*this->lambda(i));
  g.setLeftBound((((x.leftBound()-b.leftBound())*ex)+b.leftBound()).leftBound());
  g.setRightBound((((x.rightBound()-b.rightBound())*ex)+b.rightBound()).rightBound());
  return g;
}

template <class TailT, class PDET>
template <typename VectorType>
typename Algorithms2<TailT, PDET>::ScalarType Algorithms2<TailT, PDET>::D1(const VectorType& x, const TailType& T) const {
  ScalarType c = C(T);
  ScalarType es = s(T);
  ScalarType M = this->M;
  ScalarType A = this->template sumOfSup<VectorType, TailType> (x, T, this->m_a0);
  ScalarType result = power(ScalarType(2.), es + 1) * c * (c / ((es - 1.) * power(ScalarType(M), es - 1.)) +
      (c * power(ScalarType(2.), es - 1)) / power(ScalarType(2 * M + 1.), es) + A);
  return result;
}

template <class TailT, class PDET>
template <typename VectorType>
typename Algorithms2<TailT, PDET>::ScalarType Algorithms2<TailT, PDET>::D2(const VectorType& u, const TailType& tail) const {

  std::vector<ScalarType> r_du(this->M);
#if __nDF_DEBUG__
  bool m_d=true;
#endif

  ScalarType sumRe, sumIm;
  int k;
  int first = this->M + 1;
  int last = this->M * 2;
  ScalarType A = this->template sumOfSup<VectorType, TailType> (u, tail, this->m_a0);

#if __nDF_DEBUG__ || __D12inf_DEBUG__
  std::ofstream f;
  f.open("nDF.txt", std::ofstream::app);
#endif
  ScalarType max = 0;
  bool b = true;
  for(k = first; k <= last && b; ++k) {
    if(tail.inTail(k)) {
#if __nDF_DEBUG__
      if(m_d) f<<"k="<<k<<"\n";
      if(m_d) tail.print(f);
#endif
      sumRe = 0;
      sumIm = 0;
      int k_1;

      //we do not distinct a_k: k<0 from a_k: k>0, because norm is the same
      for(k_1 = k - this->M; k_1 < k / 2.; ++k_1) {
        if(!(tail.inFarTail(k_1) && tail.inFarTail(k - k_1))) {//we must to check if AT LEAST one mode of convolution is in the tail
#if __nDF_DEBUG__
          if(m_d) f<<"k_1="<<k_1<<" k-k_1="<<k-k_1<<"\n";
#endif
          if(k_1 > this->m || k_1 < -this->m) {
#if __nDF_DEBUG__
            if(m_d) f<<"tail(k_1): Re="<<tail(k_1, 1)<<" Im="<<tail(k_1, 0)<<"\n";
#endif
          } else {
#if __nDF_DEBUG__
            if(m_d) f<<"modes(k_1): Re="<<((fadbad::F<INTERVAL, 0u>)this->template modeVal<VectorType>(k_1, 1, u, this->m_a0)).x()<<" Im="<<((fadbad::F<INTERVAL, 0u>)this->template modeVal<VectorType>(k_1, 0, u, this->m_a0)).x()<<"\n";
#endif
          }
          if(k - k_1 > this->m || k - k_1 < -this->m) {
#if __nDF_DEBUG__
            if(m_d) f<<"tail(k-k_1): Re="<<tail(k-k_1, 1)<<" Im="<<tail(k-k_1, 0)<<"\n";
#endif
          } else {
#if __nDF_DEBUG__
            if(m_d) f<<"modes(k-k_1): Re="<<((fadbad::F<INTERVAL, 0u>)this->template modeVal<VectorType>(k-k_1, 1, u, this->m_a0)).x()<<" Im="<<((fadbad::F<INTERVAL, 0u>)this->template modeVal<VectorType>(k-k_1, 0, u, this->m_a0)).x()<<"\n";
#endif
          }

          this->template product<VectorType> (k, k_1, sumRe, sumIm, u, tail);
#if __nDF_DEBUG__
          if(m_d) f<<"sumRe="<<((fadbad::F<INTERVAL, 0u>)sumRe).x()<<" sumIm="<<((fadbad::F<INTERVAL, 0u>)sumIm).x()<<"\n";
#endif
        }
      }
      //check if k is even then add additional term
      if(k % 2 == 0) {
        this->template product<VectorType> (k, k / 2, sumRe, sumIm, u, tail, 0.5);
      }

#if __nDF_DEBUG__
      if(m_d) f<<"k="<<k<<" calculated sumRe="<<sumRe<<" calculated sumIm="<<sumIm<<" norm="<<this->norm(sumRe, sumIm)<<" pow="<<POW(INTERVAL(k), s(tail))<<" re n("<<k<<")="<<this->template n<VectorType>(k, 1, u, tail, f)<<" im n("<<k<<")="<<this->template n<VectorType>(k, 0, u, tail, f)<<"\n";
#endif

      r_du[k - this->M - 1] = 2. * (A * C(tail) * power(ScalarType(2.), s(tail)) + this->norm(sumRe, sumIm) * power(ScalarType(k), s(tail)));

#if __D12inf_DEBUG__
      ScalarType dre=this->template n<VectorType>(k, 1, u, tail, f);
      ScalarType dim=this->template n<VectorType>(k, 0, u, tail, f);
      f<<"D1["<<(k)<<"]="<<r_du[k-this->M-1]<<" aaa D1/k^s="<<(abs(this->m_N_coeff)*(r_du[k-this->M-1])/POW(INTERVAL(k), s(tail)-1.))<<" ReN_"<<(k)<<" "<<dre<<" ImN_"<<(k)<<" "<<dim<<" norm="<<(SQRT((dre*dre).rightBound()+(dim*dim).rightBound()))<<"\n";
#endif
      if(r_du[k - this->M - 1] > max)
        max = r_du[k - this->M - 1];
      else {
        b = false;//it means that D_1_k started to decrease
      }
    }
  }
#if __D12inf_DEBUG__
  f<<"A="<<A<<" C(tail)="<<C(tail)<<" POW="<<POW(INTERVAL(2.), s(tail))<<"\n";
#endif

#if __nDF_DEBUG__
  f<<"max="<<max<<"\n";
  f.close();
#endif
  return max;
}

template <class TailT, class PDET>
template <typename VectorType>
typename Algorithms2<TailT, PDET>::ScalarType Algorithms2<TailT, PDET>::Dinf(const VectorType& x, const TailType& T) const {
  ScalarType c = C(T);
  ScalarType es = s(T);
  ScalarType M = this->M;
  ScalarType result = (2. * c * c) / ((es - 1.) * power(ScalarType(M), es - 1.));
  return result;
}

template <class TailT, class PDET>
template <typename VectorType>
inline void Algorithms2<TailT, PDET>::product(int k, int k_1, typename VectorType::ScalarType& re, typename VectorType::ScalarType& im, const VectorType& u,
    const TailType& tail, double coeff) const {
  typedef typename VectorType::ScalarType ScalarType;
  ScalarType ak1Re, akk1Re;
  ScalarType ak1Im, akk1Im;

  if(abs(k_1) <= this->m) { ///the mode a_{k_1} is located in the symmetric Galerkin projection.
    ak1Re = BaseClass::template modeVal<VectorType>(k_1, 1, u, this->m_a0);
    ak1Im = BaseClass::template modeVal<VectorType>(k_1, 0, u, this->m_a0);
  } else { ///on contrary, the mode a_{k_1} is located in the tail.
    ak1Re = tail(k_1, 1);
    ak1Im = tail(k_1, 0);
  }

  if(abs(k - k_1) <= this->m) { ///the mode a_{k-k_1} is located in the symmetric Galerkin projection.
    akk1Re = BaseClass::template modeVal<VectorType>(k - k_1, 1, u, this->m_a0);
    akk1Im = BaseClass::template modeVal<VectorType>(k - k_1, 0, u, this->m_a0);
  } else { ///on contrary, the mode a_{k-k_1} is located in the tail.
    akk1Re = tail(k - k_1, 1);
    akk1Im = tail(k - k_1, 0);
  }

  ///we calculate the i\cdot a_{k_1}\cdot a_{k-k_1} times the given coeff.
  re -= (ak1Re * akk1Im + ak1Im * akk1Re) * coeff;
  im += (ak1Re * akk1Re - ak1Im * akk1Im) * coeff;
}

template <class TailT, class PDET>
template <typename VectorType>
inline void Algorithms2<TailT, PDET>::product(int k, int k_1, typename VectorType::ScalarType& re, typename VectorType::ScalarType& im, const VectorType& u, double coeff) const {
  typedef typename VectorType::ScalarType ScalarType;
  ScalarType ak1Re, akk1Re;
  ScalarType ak1Im, akk1Im;

  if(abs(k_1) <= this->m) {
    ak1Re = BaseClass::template modeVal<VectorType>(k_1, 1, u, this->m_a0);
    ak1Im = BaseClass::template modeVal<VectorType>(k_1, 0, u, this->m_a0);
  } else {
    std::cerr << "mode " << k_1 << " is not located in the projection!\n";
    throw std::runtime_error("Requested mode that is not located in the projection. (convolution)");
  }

  if(abs(k - k_1) <= this->m) {
    akk1Re = BaseClass::template modeVal<VectorType>(k - k_1, 1, u, this->m_a0);
    akk1Im = BaseClass::template modeVal<VectorType>(k - k_1, 0, u, this->m_a0);
  } else {
    std::cerr << "mode " << (k - k_1) << " is not located in the projection!\n";
    throw std::runtime_error("Requested mode that is not located in the projection. (convolution)");
  }

  /// we calculate the product i\cdot a_{k_1}\cdot a_{k-k_1} times the given coeff.
  re -= (ak1Re * akk1Im + ak1Im * akk1Re) * coeff;
  im += (ak1Re * akk1Re - ak1Im * akk1Im) * coeff;
}

template <class TailT, class PDET>
template <typename VectorType>
VectorType Algorithms2<TailT, PDET>::vectorField(const VectorType& u) const {
  typedef BaseClass Class;
  VectorType r_du(Class::modes2arraySize(this->m));
  int k;
  typedef typename VectorType::ScalarType ScalarType;
  ScalarType sumRe, sumIm;
  int k_start, k_end;
  int first = this->firstMode();
  int last = this->lastMode();
  for(k = first; k <= last; k++) {

    sumRe = 0;
    sumIm = 0;

    if(k > 0) {
      for(int k_1 = k - this->m; k_1 < k / 2.; ++k_1) {
        if(abs(k - k_1) <= this->m) {
          this->template product<VectorType> (k, k_1, sumRe, sumIm, u, 2);
        }
      }
      if(k % 2 == 0) {
        this->template product<VectorType> (k, k / 2, sumRe, sumIm, u, 1);
      }
    }
    if(k < 0) { //without symmetry
      k_start = -this->m;
      k_end = this->m + k;
      for(int k_1 = k_start; k_1 <= k_end; k_1++) {
        if(abs(k - k_1) <= this->m) { //check if |k-k_1|<=m
          this->template product<VectorType> (k, k_1, sumRe, sumIm, u);
        }
      }
    }
    r_du[Class::mode2array(k, 1)] = -k * this->m_N_coeff * sumRe + this->lambda_k(k) * u[Class::mode2array(k, 1)];
    r_du[Class::mode2array(k, 0)] = -k * this->m_N_coeff * sumIm + this->lambda_k(k) * u[Class::mode2array(k, 0)];
  }
  return r_du;
}

template <class TailT, class PDET>
template <typename TimeT, typename VectorType>
VectorType Algorithms2<TailT, PDET>::nonLinear(const TimeT &time, const VectorType& in) const
{
    VectorType out( in.dimension() );
    VectorType nonLinearVector( in.dimension() );
    VectorType newIn( in.dimension() );

    VectorType remainsFromLinear( in.dimension() );

    VectorType nonAutonomousPart( in.dimension() );
    data.equation->template assembleNonautonomousPart<TimeT,VectorType>(time, nonAutonomousPart);

    newIn = capd::vectalg::matrixByVector<VectorType>( data.P, in );

    data.equation->assembleNonLinearVector( &nonLinearVector, &newIn);

    out = capd::vectalg::matrixByVector<VectorType>( data.matrixForNonlinearPart,
                                                 capd::vectalg::subtractObjects<VectorType>( capd::vectalg::subtractObjects<VectorType>( data.loadVector, nonLinearVector ), nonAutonomousPart ) );

    remainsFromLinear = capd::vectalg::matrixByVector<VectorType>( data.remainsFromDiagonal, in);
    out += remainsFromLinear;
    return out;
}

template <class TailT, class PDET>
template <typename TimeT, typename VectorType>
VectorType Algorithms2<TailT, PDET>::N(const TimeT &time, const VectorType& in) const
{
    return data.extractHead( nonLinear(time, data.upsizeVector(in) ) );
}

template <class TailT, class PDET>
template <typename TimeT, typename VectorType>
VectorType Algorithms2<TailT, PDET>::N(const TimeT &time,const VectorType& u, const TailType& tail) const
{
    return data.extractHead( nonLinear(time, data.joinHeadWithTail(u,tail) ) );
}//returns non-dissipative

template <class TailT, class PDET>
template <typename TimeT, typename VectorType>
TailT Algorithms2<TailT, PDET>::n(const TimeT &time, const VectorType& u, const TailType& tail) const
{
    return data.extractTail( nonLinear( time, data.joinHeadWithTail(u,tail) ), tail );
}//returns dissipative


template <class TailT, class PDET>
typename Algorithms2<TailT, PDET>::ScalarType Algorithms2<TailT, PDET>::firstIntegral(int j, int k, int i) const
{
    if( (k == i) || (k == i-1) || (k==i+1) )
    {
        if( (i-1==k) && (i-1==j) ) return -1.0/3;
        if( (i-1==k) && (i==j) ) return -1.0/6;
        if( (i==k) && (i-1==j) ) return -1.0/6;
        if( (i==k) && (i+1==j) ) return 1.0/6;
        if( (i+1==k) && (i==j) ) return 1.0/6;
        if( (i+1==k) && (i+1==j) ) return 1.0/3;
        return 0;
    }
    else return 0;
}

template <class TailT, class PDET>
typename Algorithms2<TailT, PDET>::ScalarType Algorithms2<TailT, PDET>::secondIntegral(int j, int k, int i) const
{
    if( (k == i) || (k == i-1) || (k==i+1) )
    {
        if( (i-1==k) && (i-1==j) ) return -1.0/6;
        if( (i-1==k) && (i==j) ) return -1.0/3;
        if( (i==k) && (i-1==j) ) return 1.0/6;
        if( (i==k) && (i+1==j) ) return -1.0/6;
        if( (i+1==k) && (i==j) ) return 1.0/3;
        if( (i+1==k) && (i+1==j) ) return 1.0/6;
        return 0;
    }
    else return 0;
}


//first iteration of optimization - reduction of loops by better nesting
template <class TailT, class PDET>
template <typename TimeT, typename VectorType>
VectorType Algorithms2<TailT, PDET>::perturbations(const TimeT &time,const VectorType& u, const TailType& tail) const
{

    VectorType rtn( data.numberOfModes );
    VectorType headTail = data.joinHeadWithTail(u,tail);

    for( int i = 0; i < data.numberOfModes; i++ )//although we cut head(nonDissipativeModes) from rtn we can not iterate from 0-m because later we multiply the vector by matrix B^-1A^-1
    {
        for( int l = data.numberOfNonDissipativeModes; l < data.numberOfModes; l++ )
        {
            for( int n = 0; n < data.numberOfModes; n++)
            {
                rtn[i] += headTail[l] * headTail[n] * (data.bigSum[i])[l][n];
            }
        }
    }

    rtn *= -1;

    capd::vectalg::Matrix<ScalarType,0,0> linear;
    linear = data.diagonalMatrix + data.remainsFromDiagonal;

    for( int i = 0; i < data.numberOfNonDissipativeModes; i++ )// because we return only head(nonDissipativeModes) and the modes are at the begginging of the vector we may have shorter iteration 0-m in other case we had to write 0-M
    {
        for( int j = data.numberOfNonDissipativeModes; j < data.numberOfModes; j++ )
        {
            rtn[i] += linear[i][j] * headTail[j];
        }
    }
std::cout << "pert: " << data.extractHead( rtn ) << std::endl;
return data.extractHead( rtn );
}


template <class TailT, class PDET>
template <typename Scalar>
Scalar Algorithms2<TailT, PDET>::inflate(Scalar const & x, double c) const {
  typedef typename Scalar::BoundType BoundT;
  BoundT left = x.leftBound();
  BoundT right = x.rightBound();
  BoundT dia = diam(x).leftBound();
  BoundT center = (left + right) / 2;
  BoundT sigma = c * dia / 2;//some small value added to prevent from this function doing nothing because its 0
  Scalar r(center - sigma, center + sigma);
  return r;
}

}}

#endif //_CAPD_ALGORITHMS2_H_
