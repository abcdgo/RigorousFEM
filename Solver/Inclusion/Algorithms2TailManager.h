#ifndef _CAPD_JAC0_SECOND_ALGORITHM_SEVENTH_TAIL_VALIDATION_H_
#define _CAPD_JACO_SECOND_ALGORITHM_SEVENTH_TAIL_VALIDATION_H_

#include "../FadMaps/Inclusion/Data.hpp"
#include "../DissipativeEnclosure/SetBlowUpException.h"
#include "Algorithms2.h"


namespace capd {
namespace jaco {


enum AlogirthmMode {
  timeCritical, magnitudeOfTailCritical
};

///Class provides algorithms for handling of tails with variable dimensions.
///The field m_mode defines which mode of the algorithm is used ,currently available are (defined as enumeration class)
///-timeCritical - algorithm aims at using smallest possible finite part of current tail (M), which makes calculations as fast as possible,
///-magnitudeOfTailCritical - algorithm aims at minimising magnitude of current tail, used when projection is already in the trapping region,
///and we only need to push tail into it.

template <class TailT, class PDET>
class Algorithms2TailManager : public Algorithms2<TailT, PDET> {
public:
  typedef typename capd::jaco::Algorithms2<TailT, PDET> BaseClass;
  typedef typename capd::jaco::Algorithms2TailManager<TailT, PDET> Class;
  typedef TailT TailType;
  typedef typename TailType::DoubleType DoubleType;
  typedef typename BaseClass::ScalarType ScalarType;

  int maxM;
  TailType T0;
  ScalarType step;
  int m_mode;

  Data<ScalarType> &data;

  using BaseClass::m_T;

  Algorithms2TailManager(int m, int M, ScalarType ni, ScalarType a0, ScalarType s, Data<ScalarType> &d, int maxM = 1000 ) :
    BaseClass(m, M, ni, a0, s, d), data(d) {
    maxM = maxM;
    m_T = initializeTail(_C, s, this->m, this->M);
    m_mode = timeCritical;
  }

  Algorithms2TailManager(ScalarType ni, int maxM = 1000) :
    BaseClass(ni) {
    this->maxM = maxM;
    m_mode = timeCritical;
  }

  Algorithms2TailManager(int m, int M, int maxM = 1000) :
    BaseClass(m, M) {
    this->maxM = maxM;
    m_mode = timeCritical;
  }

  Algorithms2TailManager() :
    BaseClass() {
    m_mode = timeCritical;
  }

  Algorithms2TailManager(const Algorithms2TailManager& a2) :
    BaseClass(a2), data(a2.data) {
    m_mode = timeCritical;
  }

  void setMode(int newMode) {
    m_mode = newMode;
  }

  int getMode() const{
    return m_mode;
  }

  ///Calculates g_k^{+-} = (T0_k^{+-}-b_k^{+-})e^{-\lambda_k h}+b_k^{+-}, for k: |k|>m, returns the result as an instance
  ///of PolyBd structure.
  TailType gt_sb(const ScalarType& h, const TailType& T0, const TailType& N) const;

  ///Changes dimension M of the tails T0 and T.
  int changeM(int newM, TailType& T0, TailType& T, bool& mWasChanged, bool guard = false);

  ///Sets a provided M.
  void setM(int newM){ this->M = newM; }

  ///Returns the truncation of a double upto the decimalPlaces after the comma.
  inline DoubleType truncate(DoubleType& d, int decimalPlaces) const;

  ///Changes the decrease speed in the far tail of T, by changing the power in the denominator.
  void changeS(TailType& T, const ScalarType& newS) const;

  ///Calculates a value C(N) that bounds the far part of N(x+T), in the following sense N_k(x+T) <= C(N)/|k|^s.
  template <class VectorType>
  ScalarType calculateCN(const VectorType& u, TailType& tail, std::ostream& f, bool changeT = true) const;

  ///Calculates a value C(g) that bounds the far part of g(x+T), in the following sense g_k(x+T) <= C(g)/|k|^s.
  template <class VectorType>
  ScalarType calculateCG(const ScalarType& h, const VectorType& u, const TailType& T0, TailType& T, std::ostream& f, bool changeT = true) const;

  ///Investigates what is potential dimension of the candidate tail T when we increase/decrease its power s(T) in the denominator
  ///C(T) / |k|^s(T).
  template <class VectorType>
  void investigate(TailType& T, const TailType& T0, const VectorType& W_2, const ScalarType& step, DoubleType& increasedS_dL,
      DoubleType& initialS_dL, DoubleType& decreasedS_dL, std::ostream& f) const;

  ///Changes current M of the tails T and T0.
  bool updateM(TailType& T, TailType& T0, int& L, DoubleType& dL, int previousL, DoubleType previousdL, bool& validated, bool& farTailValidate,
      bool& changedLastStep, bool& mWasChanged, std::ostream& f);

  ///function recalculating L, it is using C(g), because it needs current value, due to the fact that T was updated.
  void recalculateL(const TailType& T, const TailType& gt_b, ScalarType& Cg, int& L, DoubleType& dL);

  template <class VectorType>
  bool updateT(TailType& T, TailType& gt_b, TailType& T0, const VectorType& u, ScalarType h, ScalarType newC, int& L, DoubleType& dL, int& previousL,
      DoubleType& previousdL, bool& validated, bool& farTailValidate, std::ostream& f);

  ///The main function, consisting of two parts
  ///-validation of the near tail,
  ///-validation of the far tail.
  ///It handles tails with variable dimensions, the dimension of a tail T may change during execution of this procedure,
  ///as presented in Section 8 of the paper.
  template <typename TimeT, typename VectorType>
  bool validateT(const TimeT &time, const VectorType& W_2, TailType& T0, TailType& T, const typename VectorType::ScalarType& step, bool& validate,
      bool& farTailValidate, double guessedC, int& previousL, typename TailType::DoubleType& previousdL, bool& changedLastStep, bool& firstStep,
      std::ostream& f);

  int getMaxM() const;

  ///Finds suitable exponent s(T), at the beginning of the Step 1 in the Algorithm 1 presented in the paper.
  template <typename VectorType>
  ScalarType findS(TailType& T, TailType& T0, const TailType& Tprev, const VectorType& W_2, const ScalarType& step, bool& changedLastStep,
      std::ostream& f);

  ///Check if potentialM is sufficiently larger than current M (decide if we change the dimension of a tail or not).
  bool suitableM(int newM, bool decrease) const;

  ///Finds potentially suitable dimension M.
  int potentialM(const DoubleType& basis, bool decrease) const;

  ///Returns starting T.
  template <typename VectorType>
  TailType startT(const VectorType& x, TailType& T0, const TailType& Tprev, const typename VectorType::ScalarType& step, bool& changedLastStep,
      std::ostream& f);

  ///Calculates the enclosures [W_2]: \varphi([0,h], [x_0], T)\subset [W_2] and T: T([0,h])\subset T .Step 1 of the Algorithm 1,
  ///as presented in the paper.
  template <typename TimeT,class VectorType, class DiffIncl>
  VectorType enclosureWithTail(const TimeT &time, const DiffIncl* diffIncl, const VectorType& x, TailType& T0, TailType& T, TailType& N,
      const typename DiffIncl::ScalarType& step, std::ostream& f);

  ///Calculates T(h), Step 11 of Algorithm 1. T00 is not used in this algorithm.
  template <class Scalar, class VectorType>
  TailType calculateTh(const VectorType& x, const TailType& T00, const TailType& T, const TailType& N, const VectorType& W_2, const Scalar& step, std::ostream& f);

  void setD(bool f) { this->d = f; }

};

///==========================================function definitions====================================================

template <class TailT, class PDET>
typename Algorithms2TailManager<TailT, PDET>::TailType Algorithms2TailManager<TailT, PDET>::gt_sb(
    const ScalarType& h, const TailType& T0, const TailType& N) const {
  TailType gt(N);
  int k;
  int first = 0;
  int last = gt.getCloseTailSize();

  for(k = first; k < last; ++k) {
    gt.set(k, g(h, k, T0, N));
  }

  ScalarType sb = s(N) + this->m_p;
  ScalarType Cb = C(N) / this->V(this->M + 1);
  if(sb > s(T0)) {
    ScalarType r = sb - s(T0);
    int maximum = this->maximumPoint(h, r, this->M + 1);
    if(maximum < this->M + 1)
      maximum = this->M + 1;
    setC(gt, (C(T0)) * exp(h * this->lambda_k(maximum)) * power(ScalarType(maximum), r) + Cb);
  } else {
    setC(gt, C(T0) * exp(h * this->lambda_k(this->M + 1)) * power(ScalarType(this->M + 1.), sb - s(T0)) + Cb);
  }
  setS(gt, sb);
  return gt;
}

template <class TailT, class PDET>
int Algorithms2TailManager<TailT, PDET>::changeM(int newM, TailType& T0, TailType& T, bool& mWasChanged, bool guard) {
  mWasChanged = true;
  int val;
  val = T.changeM(newM, guard);

  T0.changeM(val, guard);
  this->M = val;
  return val;
}

template <class TailT, class PDET>
inline typename Algorithms2TailManager<TailT, PDET>::DoubleType Algorithms2TailManager<TailT, PDET>::truncate(
    DoubleType& d, int decimalPlaces) const {
  int i;
  DoubleType factor = 1;
  for(i = 0; i < decimalPlaces; ++i)
    factor *= 10;
  return FLOOR(d * factor) / factor;
}

template <class TailT, class PDET>
void Algorithms2TailManager<TailT, PDET>::changeS(TailType& T, const ScalarType& newS) const {
  ScalarType oldS = s(T);
  setS(T, newS);
  setC(T, C(T) * power(ScalarType(this->M + 1), s(T) - oldS));
}

template <class TailT, class PDET>
template <class VectorType>
typename Algorithms2TailManager<TailT, PDET>::ScalarType Algorithms2TailManager<TailT, PDET>::calculateCN(
    const VectorType& u, TailType& tail, std::ostream& f, bool changeT) const {
  ScalarType potentialC = 0;
  ScalarType D1 = this->template D1<VectorType> (u, tail);
  ScalarType D2 = this->template D2<VectorType> (u, tail);
  ScalarType Dinf = this->template Dinf<VectorType> (u, tail);
  ScalarType D_F = (D1 > D2 ? D1 : D2);
  return abs(this->m_N_coeff) * (D_F + Dinf);
}

template <class TailT, class PDET>
template <class VectorType>
typename Algorithms2TailManager<TailT, PDET>::ScalarType Algorithms2TailManager<TailT, PDET>::calculateCG(
    const ScalarType& h, const VectorType& u, const TailType& T0, TailType& T, std::ostream& f, bool changeT) const {
  ScalarType Cb = this->calculateCN(u, T, f, changeT) / this->V(this->M + 1);
  ScalarType sb = s(T) + this->m_p - 1;
  if(sb >= s(T0)) {
    ScalarType r = sb - s(T0);
    int maximum = this->maximumPoint(h, r, this->M + 1);
    if(maximum < this->M + 1)
      maximum = this->M + 1;
    return (C(T0)) * exp(h * this->lambda_k(maximum)) * power(ScalarType(maximum), r) + Cb;
  } else {
    return C(T0) * exp(h * this->lambda_k(this->M + 1)) * power(ScalarType(this->M + 1.), sb - s(T0)) + Cb;
  }
}

template <class TailT, class PDET>
template <class VectorType>
void Algorithms2TailManager<TailT, PDET>::investigate(TailType& T, const TailType& T0, const VectorType& W_2,
    const ScalarType& step, DoubleType& increasedS_dL, DoubleType& initialS_dL, DoubleType& decreasedS_dL, std::ostream& f) const {
  TailType cpyT(T);
  ScalarType sT = s(cpyT);
  ScalarType Cg = this->calculateCG(step, W_2, T0, cpyT, f);
  ScalarType sg = s(cpyT) + this->m_p - 1.;
  if(C(cpyT) != 0 && Cg != 0)
    initialS_dL = POW(C(cpyT) / Cg, 1. / (s(cpyT) - sg)).rightBound();
  else
    initialS_dL = -1;
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
  f<<"investigating s(T)="<<sT<<" Cg="<<Cg<<" CT="<<C(cpyT)<<" initialS_dL="<<initialS_dL<<"\n";
#endif

  //investigation of T with increased power
  this->changeS(cpyT, sT + 1.);
  Cg = this->calculateCG(step, W_2, T0, cpyT, f);
  sg = s(cpyT) + this->m_p - 1.;
  if(C(cpyT) != 0 && Cg != 0)
    increasedS_dL = POW(C(cpyT) / Cg, 1. / (s(cpyT) - sg)).rightBound();
  else
    increasedS_dL = -1;
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
  f<<"investigating s(T)="<<sT+1<<" Cg="<<Cg<<" CT="<<C(cpyT)<<" increasedS_dL="<<increasedS_dL<<"\n";
#endif

  //investigation of T with decreased power (only if s>sufficientS)
  if(sT > 1 + this->m_r) {
    this->changeS(cpyT, sT - 1.);
    Cg = this->calculateCG(step, W_2, T0, cpyT, f);
    sg = s(cpyT) + this->m_p - 1.;
    if(C(cpyT) != 0 && Cg != 0)
      decreasedS_dL = POW(C(cpyT) / Cg, 1. / (s(cpyT) - sg)).rightBound();
    else
      decreasedS_dL = -1;
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
    f<<"investigating s(T)="<<sT-1<<" Cg="<<Cg<<" CT="<<C(cpyT)<<" decreasedS_dL="<<decreasedS_dL<<"\n";
#endif
  }
}

template <class TailT, class PDET>
bool Algorithms2TailManager<TailT, PDET>::updateM(TailType& T, TailType& T0, int& L, DoubleType& dL, int previousL,
    DoubleType previousdL, bool& validated, bool& farTailValidate, bool& changedLastStep, bool& mWasChanged, std::ostream& f) {
  int newM = -1;
  if(dL > previousdL) {
    int maxM = getMaxM();
    if(L < maxM)
      newM = potentialM(L, false);
    else {
      if(L > __SET_BLOWUP_L_THRESHOLD__) {
        std::cerr << "Blow up! Validate tail function.\n";
        throw capd::jaco::SetBlowUpException<ScalarType>(dL, "validate tail");
      }
      if(L < 2 * this->M) //if L is not much larger than current M we increase M by little
        newM = potentialM(this->M, false)/*static_cast<int>(this->M*__INCREASE_L__+__L_CONST__)*/;
      else
        //otherwise if L is much larger than current M we increase M by large, to make it larger than L
        newM = potentialM(L, false);
    }
  } else {
    if(truncate(dL, __L_TRUNCATION_DP__) == truncate(previousdL, __L_TRUNCATION_DP__)) {
      if(L < this->M)
        newM = potentialM(L, true)/*static_cast<int>(L*__INCREASE_L__+__L_CONST__)*/;
      else if(L > this->M)
        newM = potentialM(L, false);
    }
  }
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
  f<<"newM="<<newM<<"\n";
  f<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
  if(!changedLastStep && newM > 0 && previousdL > 0) {
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
    f<<"current M="<<this->M<<"\n";
#endif

    this->changeM(newM, T0, T, mWasChanged);

#if __SM__
    f<<"M was updated, current M: "<<newM<<"\n";
#endif
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
    f<<"validate of the far tail failed, increasing M, newM="<<newM<<", close tail size="<<T.getCloseTailSize()<<"\n";
    f<<"tail NOT validated.\n";
#endif
    validated = false;
    farTailValidate = false;
    return true;
  } else {
    return false;
  }
}

template <class TailT, class PDET>
void Algorithms2TailManager<TailT, PDET>::recalculateL(const TailType& T, const TailType& gt_b, ScalarType& Cg, int& L,
    DoubleType& dL) {
  if(Cg != 0 && C(T) != 0) {
    L = static_cast<int> (CEIL(POW(C(T) / Cg, 1. / (s(T) - s(gt_b))).rightBound())) + 1;
    dL = POW(C(T) / Cg, 1. / (s(T) - s(gt_b))).rightBound();
  } else {
    L = -1;
    dL = -1;
  }
}

template <class TailT, class PDET>
template <class VectorType>
bool Algorithms2TailManager<TailT, PDET>::updateT(TailType& T, TailType& gt_b, TailType& T0, const VectorType& u, ScalarType h,
    ScalarType newC, int& L, DoubleType& dL, int& previousL, DoubleType& previousdL, bool& validated, bool& farTailValidate, std::ostream& f) {
  bool r = false;
  if(!(newC <= C(T))) {
    r = true;
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
    f<<"changing C(T) to "<<__D_2__*newC<<"\n";
#endif
    setCLarger(T, __D_2__ * newC);
    ScalarType Cg = calculateCG(h, u, T0, T, f, false);
    recalculateL(T, gt_b, Cg, L, dL);//recalculates L using current values
    validated = false;
    farTailValidate = false;
  }
  return r;
}

template <class TailT, class PDET>
template <typename TimeT, typename VectorType>
bool Algorithms2TailManager<TailT, PDET>::validateT(const TimeT &time, const VectorType& W_2, TailType& T0, TailType& T,
    const typename VectorType::ScalarType& step, bool& validate, bool& farTailValidate, double guessedC, int& previousL,
    typename TailType::DoubleType& previousdL, bool& changedLastStep, bool& firstStep, std::ostream& f) {
  bool real, validated = true;
  int first = 0;
  int last = T.getCloseTailSize();
  int i, j, k;
  std::vector<double> inflates(this->modes2arraySize(last - first + 1), 1);
  bool inflate = false;

  bool mWasChanged = false;

  ScalarType bk;
  ScalarType g;
  TailType Nt;
  TailType bt, gt_b;

  ///===NEAR TAIL VALIDATION===
  {
    Nt = this->template n<TimeT,VectorType> (time, W_2, T);
    bt = this-> bt(step, Nt);
    gt_b = this->gt_sb(step, T0, Nt);

#if __ENCLOSURE_WITH_TAIL_DEBUG__
    f<<"Calculated n (Tail):\n";
    Nt.print(f);
#endif

    for(i = first; i < last; ++i) {
      bk = bt[i];
      g = gt_b[i];
#if __ENCLOSURE_WITH_TAIL_DEBUG__
      f<<"Validating i="<<i<<" in tail.\n";
        f<<"lambda("<<i<<")="<<this -> lambda_k(i)<<"\n";      
        f<<"bkRe="<<bk<<"\n";
      f<<"gRe="<<g<<"\n";
#endif
      inflate = false;
      if(T0[i].rightBound() < bk.rightBound()) {
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<"T[i].rightBound() "<<T[i].rightBound()<<"<g.rightBound() "<<g.rightBound()<<" odp="<<(T[i].rightBound()<g.rightBound())<<"\n";
#endif
        if(T[i].rightBound() < g.rightBound()) {
          validate = false;
          validated = false;
#if __ENCLOSURE_WITH_TAIL_DEBUG__
          f<<"validate failed, setting T(i, 1).rightBound="<<(1-__D_G__)*g.rightBound()+__D_G__*bk.rightBound()<<"\n";
#endif
          T.setRightBound(i, (1 - __D_G__) * g.rightBound() + __D_G__ * bk.rightBound());
          inflate = true;
        }
      }
      if(T0[i].leftBound() > bk.leftBound()) {
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<"T(i, 1).leftBound() "<<T[i].leftBound()<<"<gRe.leftBound() "<<g.leftBound()<<" odp="<<(T[i].leftBound()>g.leftBound())<<"\n";
#endif
        if(T[i].leftBound() > g.leftBound()) {
          validate = false;
          validated = false;
#if __ENCLOSURE_WITH_TAIL_DEBUG__
          f<<"validate failed, setting T(i, 1).leftBound="<<(1-__D_G__)*g.leftBound()+__D_G__*bk.leftBound()<<"\n";
#endif
          T.setLeftBound(i, (1 - __D_G__) * g.leftBound() + __D_G__ * bk.leftBound());
          inflate = true;
        }
      }
      if(inflate) {
        k = this->array2modeTail(i, real);
        T.set(i, this->template inflate<ScalarType> (T[i], __INFLATE_C__));
        for(j = -__INFLATE_RADIUS_1__; j <= __INFLATE_RADIUS_2__; j++) {
          if(T.inCloseTail(k + j) && j != 0) {
            //                        T.set(i+j, real, this->template inflate<Scalar>(T(i+j, 1), 1.+(__INFLATE_C__-1.)/POW(abs(j), s(T))));
            if(j < 0)
              inflates[this->mode2arrayTail(k + j, real)] *= 1. + (__INFLATE_C__ - 1.) / power(ScalarType(abs(j)), s(Nt) - 1).rightBound();
            else
              inflates[this->mode2arrayTail(k + j, real)] *= 1. + (__INFLATE_C__ - 1.) / power(ScalarType(abs(j)), s(Nt) - 1).rightBound();
          }
        }
      }
    }
    for(i = first; i < last; ++i) {
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __INFLATE_DEBUG__
      f<<"inflates["<<i<<"]="<<inflates[i]<<"\n";
#endif
      if(inflates[i] != 0)
        T.set(i, this->inflate(T[i], inflates[i]));
    }
  }
  ///===FAR TAIL VALIDATION===
  ScalarType Cb = C(bt);
  ScalarType Cg = C(gt_b);
  int L = -1;
  DoubleType L2 = -1;
  ScalarType gk, max;
  DoubleType dL = -1;

  //recalculateL(T, gt_b, Cg, L, dL);

  if(C(T0) != 0 && Cb != 0 && s(bt) != s(T0))
    L2 = CEIL(POW(Cb / C(T0), 1. / (s(bt) - s(T0))).rightBound());//recalculates L_2 using current values

  if(L2 != L2) {
    f << "Critical error.\nL2 escaped the range. L2=" << L2 << "\n";
    std::cerr << "Critical error.\nL2 escaped the range. L2=" << L2 << "\n";
    throw std::runtime_error("Critical error.\nL2 escaped the range.\n");
  }

#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
  f<<"Validating far tail\n";
  f<<"C(N)="<<C(Nt)<<"\n";
  f<<"C(b)="<<Cb<<"\n";
  f<<"C(g)="<<Cg<<"\n";
  f<<"C(T0)="<<C(T0)<<"\n";
  f<<"C(T)="<<C(T)<<"\n";
  f<<"s(N)="<<s(Nt)<<"\n";
  f<<"s(b)="<<s(bt)<<"\n";
  f<<"s(g)="<<s(gt_b)<<"\n";
  f<<"s(T)="<<s(T)<<"\n";
  f<<"s(T0)="<<s(T0)<<"\n";
  f<<"L="<<L<<"\n";
  f<<"dL="<<dL<<"\n";
  f<<"truncated dL="<<truncate(previousdL, __L_TRUNCATION_DP__)<<"\n";
  f<<"previousL="<<previousL<<"\n";
  if(L>=this->M+1) f<<"T(L,1)="<<T(L, 1)<<"\n";
  f<<"L2="<<L2<<"\n";
  f<<"previousdL="<<previousdL<<"\n";
  f<<"changedLastStep="<<changedLastStep<<"\n";
#endif

  ScalarType newC;
  farTailValidate = true;
  if(s(T0) > s(bt)) {

    if(!(C(T0) > Cb * power(ScalarType(this->M + 1.), s(T0) - s(bt)))) { //C(T0)<=C(b)(M+1)^(s(T0)-s(b))
      newC = Cg * power(ScalarType(this->M + 1.), s(T) - s(gt_b));
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
      f<<"case s(T0)>s(b), T0_{M+1}<=b_{M+1}\n";
      f<<"newC="<<__D_2__*newC<<"\n";
#endif
      if(updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f)) {
#if !__CONSTANT_M__
        updateM(T, T0, L, dL, previousL, previousdL, validated, farTailValidate, changedLastStep, mWasChanged, f);
#else
        f<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
      }
    }
    if(!(C(T0) <= Cb * power(ScalarType(this->M + 1.), s(T0) - s(bt)))) { //C(T0)>C(b)(M+1)^(s(T0)-s(b))
      //b_{M+1}<T(0)_{M+1}
      if(L2 < this->M + 1 && L2 > 0) {
        std::cerr << "Case s(T0)>s(b), T0_{M+1}>b_{M+1}\n";
        std::cerr << "Unexpected L2 " << L2 << "<=M+1 " << this->M + 1 << "!\n";
        throw std::runtime_error("Unexpected L2<=M+1!\n");
      }
      //we check if T_i<T(0)_i for i=M+1,...,L_2
      //first case T0 is decreasing faster than T, it is enough to check for i=M+1
      bool t;
      if(s(T0) >= s(T)) {
        newC = C(T0) * power(ScalarType(this->M + 1.), s(T) - s(T0));
        t = updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f);
      }
      //otherwise we check for i=L_2, T is decreasing faster than T0
      else {
        if(L2 > 0) {
          newC = C(T0) * power(ScalarType(L2), s(T) - s(T0));
          t = updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f);
        } else {
          std::cerr << "Case s(b)>s(T0), T0_{M+1}<b_{M+1}\n";
          std::cerr << "Unexpected values of L2=infty and s(T0)<s(T).\n";
          std::cerr << "with s(T0)<s(T) validation is not possible.\n";
          throw std::runtime_error("Unexpected L2<=M+1!\nWith s(T0)<s(T) validation is not possible.\n");
        }
      }
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
      f<<"case s(T0)>s(b), T0_{M+1}>b_{M+1}\nfirst step\nnewC="<<__D_2__*newC<<"\n";
#endif
      if(L2 > 0) {
        newC = Cg * power(ScalarType(L2), s(T) - s(gt_b));
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
        f<<"case s(T0)>s(b), T0_{M+1}>b_{M+1}\nsecond step\nnewC="<<__D_2__*newC<<"\n";
#endif
        t = (t || updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f));
      }
      if(t) {
#if !__CONSTANT_M__
        updateM(T, T0, L, dL, previousL, previousdL, validated, farTailValidate, changedLastStep, mWasChanged, f);
#else
        f<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
      }
    }
  }
  if(s(bt) == s(T0)) {
    ///we check the condition !(C(T0)<C(b)) instead of C(T0)>=C(b) because both may be the case
    ///(if intervals C(T0) and C(b) overlap
    if(!(C(T0) < Cb)) {//C(T0)>=C(b)
      newC = C(T0) * power(ScalarType(this->M + 1.), s(T) - s(T0));
      updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f);
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
      f<<"case s(b)==s(T0), C(T0)>=C(b)\n";
      f<<"newC="<<__D_2__*newC<<"\n";
#endif
    }

    if(!(C(T0) >= Cb)) {//C(T0)<C(b)
      if(s(gt_b) < s(T)) {
        std::cerr << "Case s(b)==s(T0), C(T0)<C(b)\n";
        std::cerr << "Unexpected s(g).\n";
        std::cerr << "With s(g)<s(T) validation is not possible.\n";
        throw std::runtime_error("Unexpected s(T)!\nWith s(g)<s(T) validation is not possible.\n");
      }
      newC = Cg * power(ScalarType(this->M + 1.), s(T) - s(gt_b));
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
      f<<"case s(b)==s(T0), C(T0)<C(b)\n";
      f<<"newC="<<__D_2__*newC<<"\n";
#endif
      if(updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f)) {
#if !__CONSTANT_M__
        updateM(T, T0, L, dL, previousL, previousdL, validated, farTailValidate, changedLastStep, mWasChanged, f);
#else
        f<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
      }
    }
  }

  if(s(bt) > s(T0)) {
    if(!(C(T0) < Cb * power(ScalarType(this->M + 1.), s(T0) - s(bt)))) {
      newC = C(T0) * power(ScalarType(this->M + 1.), s(T) - s(T0));
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
      f<<"case s(b)>s(T0), T0_{M+1}>=b_{M+1}\n";
      f<<"newC="<<__D_2__*newC<<"\n";
#endif
      updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f);
    }
    if(!(C(T0) >= Cb * power(ScalarType(this->M + 1.), s(T0) - s(bt)))) {
      if(L2 < this->M + 1 && L2 > 0) {//L2<\infty
        std::cerr << "Case s(b)>s(T0), T0_{M+1}<b_{M+1}\n";
        std::cerr << "Unexpected L2 " << L2 << "<=M+1 " << this->M + 1 << "!\n";
        throw std::runtime_error("Unexpected L2<=M+1!\n");
      }
      bool t;
      if(s(gt_b) >= s(T)) {
        newC = Cg * power(ScalarType(this->M + 1.), s(T) - s(gt_b));
        t = updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f);
      }
      //otherwise we check for i=L_2, T is decreasing faster than g
      else {
        if(L2 > 0) {//in fact we check if L2<\infty
          newC = Cg * power(ScalarType(L2), s(T) - s(gt_b));
          t = updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f);
        } else {
          std::cerr << "Case s(b)>s(T0), T0_{M+1}<b_{M+1}\n";
          std::cerr << "Unexpected values of L2=infty and s(g)<s(T).\n";
          std::cerr << "Equation is probably not dissipative, reason: p<=q.\n";
          throw std::runtime_error("Unexpected L2<=M+1!Equation is probably not dissipative.\n");
        }
      }
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
      f<<"case s(b)>s(T0), T0_{M+1}<b_{M+1}\nfirst step\n";
      f<<"newC="<<__D_2__*newC<<"\n";
#endif
      if(L2 > 0) {//L2<\infty
        newC = C(T0) * power(ScalarType(L2), s(T) - s(T0));
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
        f<<"case s(b)>s(T0), T0_{M+1}<b_{M+1}\nsecond step\n";
        f<<"newC="<<__D_2__*newC<<"\n";
#endif
        t = (t || updateT(T, gt_b, T0, W_2, step, newC, L, dL, previousL, previousdL, validated, farTailValidate, f));
      }
      if(t) {
#if !__CONSTANT_M__
        updateM(T, T0, L, dL, previousL, previousdL, validated, farTailValidate, changedLastStep, mWasChanged, f);
#else
        f<<"previousdL="<<previousdL<<" dL="<<dL<<"\n";
#endif
      }
    }
  }
  changedLastStep = mWasChanged;
  previousdL = dL;
  previousL = L;
  ///near and far tail validated
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
  if(validated)
  f<<"tail VALIDATED.\n";
  else
  f<<"tail NOT validated.\n";
#endif
  return validated;
}

template <class TailT, class PDET>
int Algorithms2TailManager<TailT, PDET>::getMaxM() const {
  int maxM;
  if(m_mode == timeCritical)
    maxM = 2 * this->m;
  if(m_mode == magnitudeOfTailCritical)
    maxM = __M_MAX__;
  return maxM;
}

template <class TailT, class PDET>
template <typename VectorType>
typename Algorithms2TailManager<TailT, PDET>::ScalarType Algorithms2TailManager<TailT, PDET>::findS(TailType& T,
    TailType& T0, const TailType& Tprev, const VectorType& W_2, const ScalarType& step, bool& changedLastStep, std::ostream& f) {
  //a loop that aims at finding largest possible exponent s(T)

  DoubleType currentS_dL, decreasedS_dL, increasedS_dL;
  this->investigate(T, T0, W_2, step, increasedS_dL, currentS_dL, decreasedS_dL, f);
  int potentialM = this->potentialM(decreasedS_dL, true);
  int currentM = this->potentialM(currentS_dL, false);
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
  f<<"Current mode="<<m_mode<<"\n";
  f<<"findS beginning.\n current M="<<currentM<<" potentialM="<<potentialM<<"\n currentS="<<s(T)<<" currentS_dL="<<currentS_dL<<" C(Tprev)="<<C(Tprev)<<" C(T)="<<C(T)<<"\n";
#endif

  int maxM = getMaxM();
  while(currentM > maxM && s(T) > this->m_sufficientlyLarge) {
    changeS(T, s(T) - 1.);
    this->investigate(T, T0, W_2, step, increasedS_dL, currentS_dL, decreasedS_dL, f);
    potentialM = this->potentialM(decreasedS_dL, true);
    currentM = this->potentialM(currentS_dL, false);
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
    f<<"entered the loop.\n current M="<<currentM<<" potentialM="<<potentialM<<"\n currentS="<<s(T)<<" currentS_dL="<<currentS_dL<<" C(T)="<<C(T)<<"\n";
#endif
  }
  //setting potentially good M
  if(currentM > this->M) {
    if(suitableM(currentM, false)) {
      this->changeM(currentM, T0, T, changedLastStep);
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
      f<<"found s="<<s(T)<<" currentM="<<currentM<<"> Tprev.M="<<Tprev.M<<" changed to "<<currentM<<"\n";
#endif
    }
  }

  if(currentM < this->M) {
    if(suitableM(currentM, true)) {
      this->changeM(currentM, T0, T, changedLastStep);
#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
      f<<"found s="<<s(T)<<" currentM="<<currentM<<"< Tprev.M="<<Tprev.M<<" changed to "<<currentM<<"\n";
#endif
    }
  }
  return s(T);
}

template <class TailT, class PDET>
bool Algorithms2TailManager<TailT, PDET>::suitableM(int newM, bool decrease) const {
  if(decrease) {
    if(newM > 2 * this->m && newM < (2. - __L_THRESHOLD__) * (this->M + 1.))
      return true;
    else
      return false;
  } else {
    int maxM = getMaxM();
    if(newM > 2 * this->m && newM > __L_THRESHOLD__ * (this->M + 1.) && newM < maxM)
      return true;
    else
      return false;
  }
}

template <class TailT, class PDET>
int Algorithms2TailManager<TailT, PDET>::potentialM(const DoubleType& basis, bool decrease) const {
  if(!decrease)
    return static_cast<int> (basis * __INCREASE_L__ + __L_CONST__);
  else
    return static_cast<int> (basis * __DECREASE_L__ + __L_CONST_DECREASE__);
}

template <class TailT, class PDET>
template <typename VectorType>
typename Algorithms2TailManager<TailT, PDET>::TailType Algorithms2TailManager<TailT, PDET>::startT(
    const VectorType& x, TailType& T0, const TailType& Tprev, const typename VectorType::ScalarType& step, bool& changedLastStep, std::ostream& f) {

  TailType T(T0);
  ScalarType small = ScalarType(-1, 1) * __SMALL__;
  int i;
  int first = 0;
  int last = T.getCloseTailSize();
  for(i = first; i < last; i++) {
    T.set(i, T0[i]);
  }
  if(C(T) > 0)
    this->findS(T, T0, Tprev, x, step, changedLastStep, f);

#if __ENCLOSURE_WITH_TAIL_DEBUG__ || __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__
  f<<"found s(T)="<<s(T)<<"\n";
  f<<"set C(T)="<<C(T)<<"\n";
#endif

  return T;
}

template <class TailT, class PDET>
template <typename TimeT, class VectorType, class DiffIncl>
VectorType Algorithms2TailManager<TailT, PDET>::enclosureWithTail(
    const TimeT &time, const DiffIncl* diffIncl, const VectorType& x, TailType& T0, TailType& T, TailType& N, const typename DiffIncl::ScalarType& step, std::ostream& f) {
  VectorType W_2;
  typedef typename DiffIncl::MultiMapType MultiMapType;
  typedef typename DiffIncl::ScalarType ScalarType;

  int j, i;
  int first = this->firstModeInTail();
  int last = this->lastModeInTail();
  ScalarType bk, g;
  TailType bt, gt, T0h;
  bool validated = false;
  int stepNr = 0;
  bool validate = true;
  bool farTailValidate = false;
  bool changedLastStep = false;
  ScalarType gk;

  if(!(s(T0) >= this->m_sufficientlyLarge)) {
    f << "The assumption regarding dissipative PDE is not satisfied, i.e. s(T0)=" << s(T0) << "<" << this->m_sufficientlyLarge << "=d+p+1.\n";
    std::cerr << "The assumption regarding dissipative PDE is not satisfied, i.e. s(T0)=" << s(T0) << "<" << this->m_sufficientlyLarge
        << "=d+r+1.\n";
    throw std::runtime_error("The assumption regarding dissipative PDE is not satisfied, i.e. s(T0)<d+p+1.\n");
  }

  T = this->template startT<VectorType> (x, T0, m_T, step, changedLastStep, f);

  typename TailType::DoubleType previousdL = -1;
  int previousL = -1;
  bool firstStep = true;
  DoubleType guessedC = 0;
  //validation loop
  while(!validated && stepNr++ < __MAX_VALIDATE_STEPS__) {
    //T=this->template guessT<Scalar, VectorType>(x, T0);
#if __ENCLOSURE_WITH_TAIL_DEBUG__
    f<<"current T:\n";
    T.print(f);
    f<<"!validated\n";
#endif
    W_2 = diffIncl->diffInclusionEnclosure(time,x, T);
    validated = validateT(time, W_2, T0, T, step, validate, farTailValidate, guessedC, previousL, previousdL, changedLastStep, firstStep, f);
  }
//  std::cout << "validating steps nr=" << stepNr << "\n";

  if(!validated) {
    std::cerr << "Failed to validate a tail, number of attempts exceeded threshold value" << __MAX_VALIDATE_STEPS__ << ".\n" <<
        "Probably the set has grown too large.\n For the usage refer the documentation, for the examples refer Section 7 in the paper."<<"\n";
    throw std::runtime_error("Failed to validate a tail, number of attempts exceeded threshold value.\n Probably the set has grown too large.\n For the usage refer the documentation, for the examples refer Section 7 in the paper.");
  }

  first = 0;
  last = T.getCloseTailSize();

  //  ONLY ONE REFINEMENT STEP!!! NEVER CHANGE IT!!!
  for(j = 0; j < 1; j++) {
#if __ENCLOSURE_WITH_TAIL_DEBUG__
    f<<"refinement step: "<<j<<"\n";
#endif

    N = this->template n<TimeT,VectorType> (time, W_2, T);
    bt = this->bt(step, N);
    gt = this->gt_sb(step, T0, N);

    for(i = first; i < last; ++i) {
#if __ENCLOSURE_WITH_TAIL_DEBUG__
      f<<"i= "<<i<<"\n";
#endif
      bk = bt[i];
      g = gt[i];
      if(T0[i].rightBound() >= bk.rightBound()) {
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<" T0.right "<<T0[i].rightBound()<<">=bk.right "<<bk.rightBound()<<"\n";
#endif
        T.setRightBound(i, T0[i].rightBound());
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<" T.right=T0.right "<<T0[i].rightBound()<<"\n";
#endif
      } else {
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<" T0.right "<<T0[i].rightBound()<<"<bk.right "<<bk.rightBound()<<"\n";
#endif
        T.setRightBound(i, g.rightBound());
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<"Re T.right=g.right "<<g.rightBound()<<"\n";
#endif
      }
      if(T0[i].leftBound() <= bk.leftBound()) {
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<"T0.left "<<T0[i].leftBound()<<">=bk.left "<<bk.leftBound()<<"\n";
#endif
        T.setLeftBound(i, T0[i].leftBound());
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<"T.left=T0.left "<<T0[i].leftBound()<<"\n";
#endif
      } else {
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<"T0.left "<<T0[i].leftBound()<<"<bk.left "<<bk.leftBound()<<"\n";
#endif
        T.setLeftBound(i, g.leftBound());
#if __ENCLOSURE_WITH_TAIL_DEBUG__
        f<<"T.left=g.left "<<g.leftBound()<<"\n";
#endif
      }
    }

#if __ENCLOSURE_WITH_TAIL_DEBUG__
    f<<"chosen T:\n";
    T.print(f);
    f<<"and corresponding W_2:\n";
#endif

#if __ENCLOSURE_WITH_TAIL_DEBUG__
    this->printModes(W_2, f);
#endif
  }
  W_2 = diffIncl->diffInclusionEnclosure(time, x, T);
  N = this->template n<TimeT,VectorType> (time, W_2, T);
  ///saving calculated T
//  std::cout << "good c: " << C(T) << "\n";
  ///take intersection of the far tail of T with the far tail of T0
  if(!T0.subsetFar(T)) {
    setCLarger(T, C(T0) * power(ScalarType(this->M + 1.), s(T) - s(T0)));
  }

  ///CHECKS IF IN FACT VALIDATED TAIL IS OK, in the sense satisfies T([0,h])\subset T
#if __VERIFY_TAIL__

  if(!T0.subset(T)) {
    int first=0;
    int last=T.getCloseTailSize();
    bool close=false;
    for(i=first; i<last; ++i)
    {
        if(!T0[i].subset(T[i]))
        {
            close=true;
        
            T.set( i, intervalHull(T0[i], T[i]) );
        }
    }
    if(!close) {
      f<<"C(T0) "<<C(T0)<<" > C(T) "<<C(T)<<" C(T)*pow="<<C(T)*power(ScalarType(this->M+1.), s(T0)-s(T))<<" s(T0) "<<s(T0)<<" s(T) "<<s(T)<<"\n";
      std::cerr<<"C(T0) "<<C(T0)<<" > C(T) "<<" C(T)*pow="<<C(T)*power(ScalarType(this->M+1.), s(T0)-s(T))<<C(T)<<" s(T0) "<<s(T0)<<" s(T) "<<s(T)<<"\n";
    }
  }
  if(!T0.subset(T)) {
    throw std::runtime_error("!T0.subset(T)\n");
  }
  T0h=TailType(T);
  bt=this-> bt(step, N);
  gt=this->gt_sb(step, T0, N);
  
  
#if __ENCLOSURE_WITH_TAIL_DEBUG__
  f<<"Calculated g:\n";
  gt.print(f);
  f<<"Checking if T[0, h]subset T\n";
#endif
  //for(i=first; i<last; ++i) {
  for(i=this->m+1; i<last; ++i) {
    bk=bt[i];
    g=gt[i];
    if(T0[i].rightBound()>=bk.rightBound()) {
      T0h.setRightBound(i, T0[i].rightBound());
    } else {
      T0h.setRightBound(i, g.rightBound());
    }
    if(T0[i].leftBound()<=bk.leftBound()) {
      T0h.setLeftBound(i, T0[i].leftBound());
    } else {
      T0h.setLeftBound(i, g.leftBound());
    }
    //test
    if(!(T0h[i].subset(T[i]))) {
      f<<"tail test failed ! T[0,h] "<<T0h[i]<<" subset T "<<T[i]<<" at coordinate i="<<i<<" re\n";
      if(T0h[i].leftBound()<T[i].leftBound()) {
        f<<"T0h(i, 1).leftBound() "<<T0h[i].leftBound()<<" < T(i, 1).leftBound() "<<T[i].leftBound()<<"diff="<<T0h[i].leftBound()-T[i].leftBound()<<"\n";
      }
      if(T0h(i, 1).rightBound()>T(i, 1).rightBound()) {
        f<<"T0h(i, 1).rightBound() "<<T0h[i].rightBound()<<"> T(i, 1).rightBound() "<<T[i].rightBound()<<"diff="<<T0h[i].rightBound()-T[i].rightBound()<<"\n";
      }
      std::cerr<<"tail test failed ! T[0,h] "<<T0h[i]<<" subset T "<<T[i]<<" at coordinate i="<<i<<" re\n";
      throw std::runtime_error("Tail test failed ! T[0,h] subset T \n");
    }

#if __ENCLOSURE_WITH_TAIL_DEBUG__
    f<<"re T[0,h]_"<<i<<"="<<T0h[i]<<"im T[0,h]_"<<i<<"="<<T0h[i]<<"\n";
#endif
  }
  int L;
  int L2;
  if(C(T0)!=0 && C(bt)!=0 && s(bt)-s(T0)!=0)
  L2=static_cast<int>(POW(C(bt)/C(T0), 1./(s(bt)-s(T0))).rightBound())+1;
  else
  L2=-1;
  if(C(gt)!=0 && C(T)!=0 && s(T)-s(gt)!=0)
  L=static_cast<int>(POW(C(T)/C(gt), 1./(s(T)-s(gt))).rightBound())+1;
  else
  L=-1;

  //validation of the far tali

  if(s(bt)>s(T0)) {
    if(L2>0) { //L2<\infty
      if(L2>this->M+1) {
        if(!(C(T)>=C(gt)*power(ScalarType(this->M+1.), s(T)-s(gt)))) {
          f<<"tail test failed!\nCase s(b)>s(T0), L2>M+1\nT_{M+1}<g_{M+1}\n";
          std::cerr<<"tail test failed!\nCase s(b)>s(T0), L2>M+1\nT_{M+1}<g_{M+1}\n";
          throw std::runtime_error("tail test failed! T_{M+1}<g_{M+1}\n");
        }
        if(!(C(T)>=C(gt)*power(ScalarType(L2), s(T)-s(gt)))) {
          f<<"tail test failed!\nCase s(b)>s(T0), L2>M+1\nT_L2<g_L2\n";
          std::cerr<<"tail test failed!\nCase s(b)>s(T0), L2>M+1\nT_L2<g_L2\n";
          throw std::runtime_error("tail test failed! T_L2<g_L2\n");
        }
        if(!(C(T)>=C(T0)*power(ScalarType(L2), s(T)-s(T0)))) {
          f<<"tail test failed!\nCase s(b)>s(T0), L2>M+1\nC(T) "<<C(T)<<" < C(T0)*POW "<<C(T0)*power(ScalarType(L2), s(T)-s(T0))<<"\n";
          std::cerr<<"tail test failed!\nCase s(b)>s(T0), L2>M+1\nC(T) "<<C(T)<<" < C(T0)*POW "<<C(T0)*power(ScalarType(L2), s(T)-s(T0))<<"\n";
          throw std::runtime_error("tail test failed! C(T) < C(T0)*POW \n");
        }
      } else { //L2<=M+1
        if(s(T)>s(T0)) {
          f<<"Tail test failed!\nCase s(b)>s(T0)\nUnexpected values of L2=infty and s(T0)<s(T).\n";
          std::cerr<<"Tail test failed!\nCase s(b)>s(T0)\nUnexpected values of L2=infty and s(T0)<s(T).\n";
          throw std::runtime_error("tail test failed!\nEquation is probably not dissipative.\n");
        }
        if(!(C(T)>=C(T0)*power(ScalarType(this->M+1.), s(T)-s(T0)))) {
          f<<"tail test failed!\nCase s(b)>s(T0), L2<=M+1\nT_{M+1}<T0_{M+1}\n";
          std::cerr<<"tail test failed!\nCase s(b)>s(T0), L2<=M+1\nT_{M+1}<T0_{M+1}\n";
          throw std::runtime_error("tail test failed!\nCase s(b)>s(T0), L2<=M+1\nT_{M+1}<T0_{M+1}\n");
        }
      }

    } else { //L2=\infty
      if(s(gt)>=s(T)) {
        if(!(C(T)>=C(gt)*power(ScalarType(this->M+1.), s(T)-s(gt)))) {
          f<<"tail test failed!\nCase s(b)>s(T0)\nT_{M+1} < g_{M+1}\n";
          std::cerr<<"tail test failed!\nCase s(b)>s(T0)\nT_{M+1} < g_{M+1}\n";
          throw std::runtime_error("tail test failed! T_{M+1} < g_{M+1}\n");
        }
      } else {//s(gt)<s(T) validation is not possible
        f<<"Tail test failed!\nCase s(b)>s(T0)\nUnexpected values of L2=infty and s(g)<s(T).\nEquation is probably not dissipative, reason: p<=q.\n";
        std::cerr<<"Tail test failed!\nCase s(b)>s(T0)\nUnexpected values of L2=infty and s(g)<s(T).\nEquation is probably not dissipative, reason: p<=q.\n";
        throw std::runtime_error("tail test failed!\nEquation is probably not dissipative.\n");
      }
    }

  }
  if(s(bt)==s(T0)) {
    if(C(T0)>=C(bt)) {
      if(!(C(T)>=C(T0)*power(ScalarType(this->M+1.), s(T)-s(T0)))) {
        f<<"tail test failed !\nCase s(b)==s(T0)\nT_{M+1} < T0_{M+1}\n";
        std::cerr<<"tail test failed !\nCase s(b)==s(T0)\nT_{M+1} < T0_{M+1}\n";
        throw std::runtime_error("Tail test failed ! T_{M+1} < T0_{M+1}\n");
      }
    } else { //C(b)>C(T0)
      if(s(gt)<s(T)) {
        f<<"Tail test failed!\nCase s(b)==s(T0)\nUnexpected s(g).\nWith s(g) "<<s(gt)<<" > s(T) "<<s(T)<<" validation is not possible.\n";
        std::cerr<<"Tail test failed!\nCase s(b)==s(T0)\nUnexpected s(g).\nWith s(g) "<<s(gt)<<" > s(T) "<<s(T)<<" validation is not possible.\n";
        throw std::runtime_error("Tail test failed ! T[0,h] subset T\nUnexpected s(T)!\nWith s(g)>s(T) validation is not possible.\n");
      }
      if(!(C(T)>=C(gt)*power(ScalarType(this->M+1.), s(T)-s(gt)))) {
        f<<"tail test failed !\nCase s(b)==s(T0)\nT_{M+1} "<<C(T)<<" <g_{M+1} "<<C(gt)*power(ScalarType(this->M+1.), s(T)-s(gt))<<" \n";
        std::cerr<<"tail test failed !\nCase s(b)==s(T0)\nT_{M+1} "<<C(T)<<" <g_{M+1} "<<C(gt)*power(ScalarType(this->M+1.), s(T)-s(gt))<<" \n";
        //std::cout<<"M(T) "<<T.getM()<<", M(gt) "<<gt.getM()<<"\n";
        throw std::runtime_error("Tail test failed ! T_{M+1}<g_{M+1}\n");
      }
    }
  }

  if(s(bt)<s(T0)) {
    if(L2>0) { //L2<\infty
      if(L2>this->M+1) {
        if(!(C(T)>=C(T0)*power(ScalarType(this->M+1.), s(T)-s(T0)))) {
          f<<"tail test failed !\nCase s(b)<s(T0), L2>M+1\nT_{M+1}<T0_{M+1}\n";
          std::cerr<<"tail test failed !\nCase s(b)<s(T0), L2>M+1\nT_{M+1}<T0_{M+1}\n";
          throw std::runtime_error("Tail test failed ! T_{M+1}<T0_{M+1}\n");
        }
        if(!(C(T)>=C(T0)*power(ScalarType(L2), s(T)-s(T0)))) {
          f<<"tail test failed !\nCase s(b)<s(T0), L2>M+1\nT_{L2}<T0_{L2}\n";
          std::cerr<<"tail test failed !\nCase s(b)<s(T0), L2>M+1\nT_{L2}<T0_{L2}\n";
          throw std::runtime_error("Tail test failed ! T_{L2}<T0_{L2}\n");
        }
        if(!(C(T)>=C(gt)*power(ScalarType(L2), s(T)-s(gt)))) {
          f<<"tail test failed!\nCase s(b)>s(g), L2>M+1\nT_{M+1}<g_{M+1}\nT_{L2}<g_{L2}\n";
          std::cerr<<"tail test failed!\nCase s(b)>s(g), L2>M+1\nT_{M+1}<g_{M+1}\nT_{L2}<g_{L2}\n";
          throw std::runtime_error("tail test failed! T_{L2}<g_{L2}\n");
        }
      } else { //L2<=M+1
        if(s(T)>s(gt)) {
          f<<"Tail test failed!\nCase s(b)<s(T0)\nUnexpected values of L2=infty and s(g)<s(T).\n";
          std::cerr<<"Tail test failed!\nCase s(b)<s(T0)\nUnexpected values of L2=infty and s(g)<s(T).\n";
          throw std::runtime_error("Tail test failed ! T[0,h] subset T \n.");
        }
        if(!(C(T)>=C(gt)*power(ScalarType(this->M+1.), s(T)-s(gt)))) {
          f<<"tail test failed!\nCase s(b)<s(T0), L2<=M+1\nT_{M+1}<g_{M+1}\n";
          std::cerr<<"tail test failed!\nCase s(b)<s(T0), L2<=M+1\nT_{M+1}<g_{M+1}\n";
          throw std::runtime_error("tail test failed!\nCase s(b)<s(T0), L2<=M+1\nT_{M+1}<g_{M+1}\n");
        }
      }
    } else { //L2==\infty
      if(s(T0)>=s(T)) {
        if(!(C(T)>=C(T0)*power(ScalarType(this->M+1.), s(T)-s(T0)))) {
          f<<"tail test failed !\nCase s(b)<s(T0)\nT_{M+1}<T0_{M+1}\n";
          std::cerr<<"tail test failed !\nCase s(b)<s(T0)\nT_{M+1}<T0_{M+1}\n";
          throw std::runtime_error("Tail test failed ! T_{M+1}<T0_{M+1}\n");
        }
      } else {
        f<<"Tail test failed!\nCase s(b)<s(T0)\nUnexpected values of L2=infty and s(T0)<s(T).\n";
        std::cerr<<"Tail test failed!\nCase s(b)<s(T0)\nUnexpected values of L2=infty and s(T0)<s(T).\n";
        throw std::runtime_error("Tail test failed ! T[0,h] subset T \n.");
      }
    }
  }
#endif

  m_T = T;
  this->T0 = T0;
  return W_2;
}

template <class TailT, class PDET>
template <class Scalar, class VectorType>
typename Algorithms2TailManager<TailT, PDET>::TailType Algorithms2TailManager<TailT, PDET>::calculateTh(
    const VectorType& x, const TailType& T00, const TailType& T, const TailType& N, const VectorType& W_2, const Scalar& step, std::ostream& f) {
  TailType& T0 = this->T0;


  int first = 0;
  int last = T.getCloseTailSize();
  int i;

  Scalar bk;
  Scalar g;
  TailType bt, gt;
  bt = this-> bt(step, N);
  gt = this->gt_sb(step, T0, N);

  TailType Th(gt);
#if __CALCULATE_TH_DEBUG__
  T.print(f);
  f<<"W_2:\n";
  this->printModes(W_2, f);
  f<<"Nt:\n";
  N.print(f);
#endif

  for(i = first; i < last; ++i) {
    g = gt[i];
#if __CALCULATE_TH_DEBUG__
    bk=bt(i, 0);
    f<<"i="<<i<<" gRe="<<gRe<<" gIm="<<gIm<<" lambda(k)="<<this->lambdak(i)<<" b(k)="<<bk<<" ReT0="<<T0(i, 1)<<" ImT0="<<T0(i, 0)<<" \n";
#endif
    Th.set(i, g);
  }

  Scalar far_c = C(gt);

  setC(Th, far_c);
#if __SM__
  if(s(T0) != s(Th)){
    f<<"s was updated, current s: "<<s(Th)<<"\n";
  }
#endif
#if __CALCULATE_TH_DEBUG__
  f<<"lambda="<<this->lambdak(this->M+1)<<" far_c="<<far_c<<"\n";
#endif

  return Th;
}

}
}

#endif
