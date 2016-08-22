#ifndef _CAPD_JACO_TAIL_H
#define _CAPD_JACO_TAIL_H

#include "../DissipativeEnclosure/dissipative_enclosure.hpp"

namespace capd {
namespace jaco {

///All tails are dedicated to store complex valued tails and that are considered for symmetric Galerkin
///projection.
///We use operator()(int index, bool re) when we want to obtain value of i-th mode (index), real or imaginary part boolean re,
///We use operator[](int index) when explicit location in a table is provided by index.

///The PolyBd structure representing the infinite part of a self consistent set \prod_{|k|>m}.
template <typename ScalarT, typename DimensionT, int D>
class PolyBd;

///Finite part of a tail \prod_{m<|k|<=M}.
template <typename ScalarT, typename DimensionT, int D>
class NearTail : public DimensionT {
public:
  typedef typename ScalarT::BoundType DoubleType;
  typedef ScalarT ScalarType;
  typedef capd::vectalg::Vector<ScalarType, D> VectorType;

  int m_size;
  VectorType m_tail;

  NearTail(const NearTail& ct2) :
    DimensionT(ct2.m, ct2.M), m_size(ct2.m_size), m_tail(ct2.m_tail) {
  }

  //the constructor that initialises all modes (re and im) by C/|k|^s*[-1,1]
  NearTail(const ScalarType& c, const ScalarType& s, int m, int M);

  NearTail(int m, int M) :
    DimensionT(m, M), m_size(DimensionT::modes2arraySize(M)), m_tail(m_size) {
  }

  NearTail() {
  }

  inline void resize(int newSize) {
    m_size = newSize;
    m_tail.resize(newSize);
  }

  ///Resizes the near tail and initialises all elements to the initialValue.
  inline void resize(int newSize, ScalarType initialValue);

  void inflate(double c);

  inline bool inCloseTail(int index) const;

  ///Returns a value that is at index i in the internal data storage, used for iterating.
  inline ScalarType operator[](const int& i) const {
    return m_tail[i];
  }

  ///Returns a value that is at index i in the internal data storage, used for iterating.
  inline ScalarType& operator[](const int& i) {
    return m_tail[i];
  }

  ///Returns a_k, in fact Re{a_k} or Im{a_k}
  ///if re==true returns Re{a_k},
  ///else returns Im{a_k}.
  /// If i<0 it returns the conjugate of a_k
  inline ScalarType operator()(const int i, const bool re) const;

  //procedure that sets value of i-th mode in near tail (real or imaginary)
  inline void set(int i, bool re, ScalarType value);

  inline void setLeftBound(int i, bool re, DoubleType value);

  inline void setRightBound(int i, bool re, DoubleType value);

  inline void set(int i, ScalarType value);

  inline void setRightBound(int i, DoubleType value);

  inline void setLeftBound(int i, DoubleType value);

  void setTail(const VectorType& vec) {
    m_tail = vec;
  }

  ///Checks if the near tail is subset of a second near tail.
  bool subset(const NearTail& ct2) const;

  ///Returns the middle of the near tail.
  NearTail mid() const;

  ///Calculates sum of all elements modulus supremum in the near part of the tail.
  inline ScalarType sum(bool re) const {
    return sum(this->m + 1, re);
  }

  ///Calculates sum of elements modulus supremum, starting at an index i in the near part of the tail.
  ScalarType sum(int i, bool re) const;

  operator VectorType&() {
    return m_tail;
  }

  ///Euclidean norm of the near tail.
  inline ScalarType euclNorm() const {
    return diam(m_tail).euclNorm();
  }

  ///Takes intersection of the near tail with an another near tail.
  inline void intersect(const NearTail& nt2);

  friend void operator+=(NearTail& nt1, const NearTail& nt2){
    nt1.m_tail=intervalHull(nt1.m_tail, nt2.m_tail);
  }

  friend class PolyBd<ScalarT, DimensionT, D> ;

};

///infinite part of a tail \prod_{|k|>M}
template <typename ScalarT, typename DimensionT>
class FarTail : public DimensionT {
public:
  typedef ScalarT ScalarType;

  ScalarType m_c;
  ScalarType m_s;

  FarTail(const FarTail& ft2) :
    DimensionT(ft2.m, ft2.M) {
    this->m_c = ft2.m_c;
    this->m_s = ft2.m_s;
  }

  FarTail(const ScalarType& c, const ScalarType& s, int m, int M) :
    DimensionT(m, M) {
    this->m_c = c;
    this->m_s = s;
  }

  FarTail(int m, int M) :
    DimensionT(m, M) {
  }

  FarTail() {
  }

  inline ScalarType value(int i) const {
    double j = (i > 0 ? i : -i);
    return (m_c / power(ScalarType(j), m_s)) * ScalarType(-1, 1);
  }

  ///Returns a value that is at index i in the internal data storage, used for iterating.
  inline ScalarType operator[](const int& i) const {
    return value(this->array2modeTail(i));
  }

  ///if re==true returns Re{a_i},
  ///else returns Im{a_i}.
  inline ScalarType operator()(const int i, const bool re) const {
    return value(i);
  }

  inline void setC(const ScalarType& C) {
    m_c = C;
  }

  inline void setS(const ScalarType& s) {
    this->m_s = s;
  }

  inline ScalarType getC() const {
    return m_c;
  }

  inline ScalarType getS() const {
    return m_s;
  }

  inline int getM() const {
    return this->M;
  }

  inline void setM(int M) {
    this->M = M;
  }

  inline bool inFarTail(int index) const;

  ///Checks if the far tail is subset of a second far tail.
  inline bool subset(const FarTail& ft2) const;

  ///Returns \prod_{|k|>M}{0}.
  inline FarTail mid() const {
    return FarTail(0., this->m_s, this->m, this->M);
  }

  ///Calculates sum of all elements absolute supremum in the tail.
  inline ScalarType sum() const {
    return sum(this->M + 1);
  }

  ///Calculates sum of elements absolute supremum in the tail, starting at an index i.
  inline ScalarType sum(int i) const;

  ///Euclidean norm of the far tail.
  inline ScalarType euclNorm() const {
    return this->getC() / ((2 * this->getS() - 1) * power(ScalarType(this->M + 1), 2 * this->getS() - 1));
  }

  ///Takes intersection of the tail with another tail ft2, assumes homogeneity of the powers s. We do not allow to store a tail
  ///with varying exponents.
  inline void intersect(const FarTail& ft2);

  friend void operator+=(FarTail& ft1, const FarTail& ft2){
    if(ft1.m_s != ft2.m_s){
      std::cerr << "operator+= of the FarTail class works only with tails of equal exponents.\n";
      throw std::runtime_error("operator+= of the FarTail class works only with tails of equal exponents.\n");
    }
    if(ft1.m_c < ft2.m_c)
      ft1.m_c = ft2.m_c;
  }

};

//a class PolyBd for storing a Tail as a sum of near TAIL (n dimensional vector) and FAR TAIL
//C(T) and s(T). NearTail + FarTail.
template <typename ScalarT, class DimensionT, int D>
class PolyBd : public DimensionT {
  typedef NearTail<ScalarT, DimensionT, D> NearTailType;
  typedef FarTail<ScalarT, DimensionT> FarTailType;
  NearTailType m_nearTail;
  FarTailType m_farTail;
public:
  typedef typename ScalarT::BoundType DoubleType;
  typedef typename NearTailType::VectorType VectorType;
  typedef ScalarT ScalarType;

  PolyBd(const ScalarType& c, const ScalarType& s, int m, int M) :
    m_nearTail(c, s, m, M), m_farTail(c, s, m, M) {
    this->m = m;
    this->M = M;
  }

  PolyBd(const ScalarType& c_near, const ScalarType& c_far, const ScalarType& s, int m, int M) :
    m_nearTail(c_near, s, m, M), m_farTail(c_far, s, m, M) {
    this->m = m;
    this->M = M;
  }

  PolyBd(int m, int M) :
    m_nearTail(m, M), m_farTail(m, M) {
    this->m = m;
    this->M = M;
  }

  PolyBd() {
  }

  PolyBd(const PolyBd& pb2) :
    m_nearTail(pb2.m_nearTail), m_farTail(pb2.m_farTail) {
    this->m = pb2.m;
    this->M = pb2.M;
  }

  PolyBd(const NearTailType& cl, const FarTail<ScalarType, DimensionT>& far) :
    m_nearTail(cl), m_farTail(far) {
    this->m = cl.m;
    this->M = cl.M;
  }

  inline void set(int i, bool re, const ScalarType& value) {
    m_nearTail.set(i, re, value);
  }

  inline void setLeftBound(int i, bool re, const DoubleType& value) {
    m_nearTail.setLeftBound(i, re, value);
  }

  inline void setRightBound(int i, bool re, DoubleType value) {
    m_nearTail.setRightBound(i, re, value);
  }

  inline void set(int i, ScalarType value) {
    m_nearTail.set(i, value);
  }

  inline void setLeftBound(int i, DoubleType value) {
    m_nearTail.setLeftBound(i, value);
  }

  inline void setRightBound(int i, DoubleType value) {
    m_nearTail.setRightBound(i, value);
  }

  //operator() gets value of a mode that is located int the tail, of index i, re tells if real or imaginary
  //part is requested. It uses closeTailVal function of Dimensions in order to obtain a value
  //from the closeTail, in this way we can write an algorithm that will be used for that.
  //On the other hand, if a mode is in the farTail it's simple, there's no need of using
  //specific algorithm.
  inline ScalarType operator()(const int i, const bool re) const;

  ///Returns a value that is at index i in the internal data storage, used for iterating.
  inline ScalarType operator[](int i) const;

  ///Subtracts two tails.
  PolyBd operator-(const PolyBd& t2) const;

  //Depreciated! use operator() instead.
  inline void updateCloseTail(int index, bool re, const ScalarType& value) {
    m_nearTail.set(index, re, value);
  }

  //Depreciated! use setC(Tail) instead.
  //Function updating C in farTail (C)/|k|^s,
  inline void updateFarTailC(const ScalarType& c) {
    m_farTail.setC(c);
  }

  //Depreciated! use setS(Tail) instead.
  //Procedure updating s in farTail C/|k|^(s).
  inline void updateFarTailS(const ScalarType& s) {
    m_farTail.setS(s);
  }

  //function evaluating value of infinite sum of conjugate of modes from the far tail
  //2\times\sum_{k_1<-M}a_{k-k_1}a_{k_1}, where a_{k-k_1} and a_{k_1} are from the far tail,
  //see Lemma 6 in the article.
  ScalarType infiniteSum(int k) const;

  ///Computes infinite sum of the tail modes absolute supremum starting at index i either upto \infty for i>0,
  ///or -\infty otherwise.
  ScalarType sum(int i, bool re) const;

  ScalarType sum(bool re) const {
    return sum(this->m + 1, re);
  }

  inline int getm() const {
    return this->m;
  }

  inline int getM() const {
    return this->M;
  }

  inline ScalarType getC() const {
    return m_farTail.getC();
  }

  inline ScalarType getS() const {
    return m_farTail.getS();
  }

  ///Checks if mode a_i belongs to the close tail.
  inline bool inCloseTail(int i) const {
    return m_nearTail.inCloseTail(i);
  }

  ///Checks if mode a_i belongs to the far tail.
  inline bool inFarTail(int i) const {
    return m_farTail.inFarTail(i);
  }

  ///Prints the tail to an output stream, formatted.
  void print(std::ostream& stream) const;

  ///Prints the tail to the output stream, unformatted.
  void printRaw(std::ostream& stream) const;

  ///Checks if the current tail is subset of a second tail.
  inline bool subset(const PolyBd& t2) const;

  inline bool subsetFar(const PolyBd& t2) const;

  inline NearTailType& getCloseTail() {
    return m_nearTail;
  }

  inline int getCloseTailSize() const {
    return m_nearTail.m_size;
  }

  inline FarTail<ScalarType, DimensionT>& getFarTail() const {
    return m_farTail;
  }

  ///The middle of a tail.
  inline PolyBd mid() const {
    return PolyBd(this->m_nearTail.mid(), this->m_farTail.mid());
  }

  ///Takes intersection of the near tail of the current tail with a second tail t2.
  inline void closeTailIntersect(const PolyBd& t2);

  ///Takes intersection of the tail with a second tail t2
  inline void intersect(const PolyBd& t2);

  ///Inflates the tail by c.
  void inflate(double c) {
    m_nearTail.inflate(c);
    m_farTail.m_c = c * m_farTail.m_c;
  }

  ///Euclidean norm or the whole tail.
  ScalarType euclNorm() const {
    return m_nearTail.euclNorm() + m_farTail.euclNorm();
  }

  ///Change dimension of the close tail.
  inline int changeM(int newM, bool guard = false);

  ///Takes interval hull of sum of the tail with a second tail
  friend void operator+=(PolyBd& pb1, const PolyBd& pb2) {
    //function works only if dimension of both tails is the same
    if(pb1.getM() == pb2.getM() && pb1.getm() == pb2.getm()) {
      pb1.m_nearTail += pb2.m_nearTail;
      pb1.m_farTail += pb2.m_farTail;
    } else {
      std::cerr << "operator+= of the PolyBd class works only with tails of equal dimensions.\n";
      throw std::runtime_error("operator+= of the PolyBd class works only with tails of equal dimensions.\n");
    }
  }

  friend std::ostream& operator<<(std::ostream& out, const PolyBd& t) {
    t.print(out);
    return out;
  }

};

///==========================================function definitions====================================================

///=============================================global functions=====================================================

template <typename Tail>
inline typename Tail::ScalarType C(const Tail& pb) {
  return pb.getC();
}

template <typename Tail>
inline typename Tail::ScalarType s(const Tail& pb) {
  return pb.getS();
}

template <typename Tail>
inline void setC(Tail& pb, const typename Tail::ScalarType& c) {
  pb.updateFarTailC(c);
}

///Chooses as C the interval b which is larger than an interval c (forall x in b and y in c, x is > than y).
template <typename Tail>
inline void setCLarger(Tail& pb, const typename Tail::ScalarType& c) {
  typedef typename Tail::ScalarType Scalar;
  Scalar newC = Scalar(c.rightBound()) + Scalar(0., diam(c).rightBound());
  pb.updateFarTailC(newC);
}

template <typename Tail>
inline void setS(Tail& pb, const typename Tail::ScalarType& s) {
  pb.updateFarTailS(s);
}

template <typename TailT>
inline TailT mid(const TailT& pb) {
  return pb.mid();
}

///=================================================methods==========================================================

template <typename ScalarT, typename DimensionT, int D>
NearTail<ScalarT, DimensionT, D>::NearTail(const ScalarType& c, const ScalarType& s, int m, int M) :
  DimensionT(m, M), m_size(DimensionT::modes2arraySize(M)), m_tail(m_size) {
  int i, k;
  int first = 0;
  int last = m_size;
  for(i = first; i < last; ++i) {
    k = this->array2modeTail(i);
    if( k != 0 ) set(i, (c / POW(ScalarType(k), s)) * ScalarType(-1, 1));
  }
}

template <typename ScalarT, typename DimensionT, int D>
inline void NearTail<ScalarT, DimensionT, D>::resize(int newSize, ScalarType initialValue) {
  m_size = newSize;
  m_tail.resize(newSize);
  m_tail = initialValue;
}

template <typename ScalarT, typename DimensionT, int D>
void NearTail<ScalarT, DimensionT, D>::inflate(double c) {
  int i;
  for(i = 0; i < m_size; ++i)
    m_tail[i] = capd::jaco::inflate<ScalarType>(m_tail[i], c);
}

template <typename ScalarT, typename DimensionT, int D>
inline bool NearTail<ScalarT, DimensionT, D>::inCloseTail(int index) const {
  if((index <= this->M && index > this->m) || (index >= -this->M && index < -this->m))
    return true;
  return false;
}

template <typename ScalarT, typename DimensionT, int D>
inline typename NearTail<ScalarT, DimensionT, D>::ScalarType NearTail<ScalarT, DimensionT, D>::operator()(const int i, const bool re) const {
  if(!this->inCloseTail(i)) {
    std::cerr << "You requested value of a mode in the near tail, but index (" << i << ") corresponds either to the projection or the far part.\n";
    throw std::runtime_error("You requested value of a mode in the near tail, but index corresponds either to the projection or the far part.\n");
  }
  return DimensionT::template closeTailVal<VectorType>(i, re, m_tail);
}

template <typename ScalarT, typename DimensionT, int D>
inline void NearTail<ScalarT, DimensionT, D>::set(int i, bool re, ScalarType value) {
    int index = this->mode2arrayTail(i, re);
    m_tail[index] = value;
}

template <typename ScalarT, typename DimensionT, int D>
inline void NearTail<ScalarT, DimensionT, D>::setLeftBound(int i, bool re, DoubleType value) {
    m_tail[this->mode2arrayTail(i, re)].setLeftBound(value);
}

template <typename ScalarT, typename DimensionT, int D>
inline void NearTail<ScalarT, DimensionT, D>::setRightBound(int i, bool re, DoubleType value) {
    m_tail[this->mode2arrayTail(i, re)].setRightBound(value);
}

template <typename ScalarT, typename DimensionT, int D>
inline void NearTail<ScalarT, DimensionT, D>::set(int i, ScalarType value) {
  if(i >= 0 && i < m_size) {
    m_tail[i] = value;
  } else {
    std::cerr << "Error. You requested to set value in a tail, of a mode that isn't in the near tail.\n";
    throw std::runtime_error("Error. You requested to set value in a tail, of a mode that isn't in the near tail.\n");
  }
}

template <typename ScalarT, typename DimensionT, int D>
inline void NearTail<ScalarT, DimensionT, D>::setRightBound(int i, DoubleType value) {
  if(i >= 0 && i < m_size) {
    m_tail[i].setRightBound(value);
  } else {
    std::cerr << "Error. You requested to set value in a tail, of a mode that isn't in the near tail.\n";
    throw std::runtime_error("Error. You requested to set value in a tail, of a mode that isn't in the near tail.\n");
  }
}

template <typename ScalarT, typename DimensionT, int D>
inline void NearTail<ScalarT, DimensionT, D>::setLeftBound(int i, DoubleType value) {
  if(i >= 0 && i < m_size) {
    m_tail[i].setLeftBound(value);
  } else {
    std::cerr << "Error. You requested to set value in a tail, of a mode that isn't in the near tail.\n";
    throw std::runtime_error("Error. You requested to set value in a tail, of a mode that isn't in the near tail.\n");
  }
}

template <typename ScalarT, typename DimensionT, int D>
bool NearTail<ScalarT, DimensionT, D>::subset(const NearTail& ct2) const {
  int i;
  for(i = 0; i <= m_size; ++i)
    if(!this->operator[](i).subset(ct2[i]))
      return false;
  return true;
  return true;
}

template <typename ScalarT, typename DimensionT, int D>
NearTail<ScalarT, DimensionT, D> NearTail<ScalarT, DimensionT, D>::mid() const {
  NearTail center(this->m, this->M);
  int last = m_tail.size();
  int i;
  for(i = 0; i < last; ++i) {
    center.m_tail[i] = capd::intervals::mid(m_tail[i]);
  }
  return center;
}

template <typename ScalarT, typename DimensionT, int D>
typename NearTail<ScalarT, DimensionT, D>::ScalarType NearTail<ScalarT, DimensionT, D>::sum(int i, bool re) const {
  if(!inCloseTail(i)) {
    std::cerr << "Requested near tail sum, starting at i=" << i << ", this index is out of the closeTail.\n";
    throw std::runtime_error("Requested near tail sum starting at index that is out of the closeTail.\n");
  }
  int start = i, end = (i > 0 ? this->M : -this->M), j;
  ScalarType sum = 0;
  for(j = start; j <= end; j++) {
    sum += rightBound(abs(this->operator()(j, re)));
  }
  return sum;
}

template <typename ScalarT, typename DimensionT, int D>
inline void NearTail<ScalarT, DimensionT, D>::intersect(const NearTail<ScalarT, DimensionT, D>& nt2) {
  m_tail=capd::vectalg::intersection(m_tail, nt2.m_tail);
}

template <typename ScalarT, typename DimensionT>
inline bool FarTail<ScalarT, DimensionT>::inFarTail(int index) const {
  if(index > this->M || index < -this->M)
    return true;
  return false;
}

template <typename ScalarT, typename DimensionT>
inline bool FarTail<ScalarT, DimensionT>::subset(const FarTail& ft2) const {
  if(ft2.getS() > this->getS())
    return false;
  int maxM = (this->M > ft2.getM() ? this->M : ft2.getM());
  if(this->operator()(maxM + 1, 1).subset(ft2(maxM + 1, 1)))
    return true;
  return false;
}

template <typename ScalarT, typename DimensionT>
inline typename FarTail<ScalarT, DimensionT>::ScalarType FarTail<ScalarT, DimensionT>::sum(int i) const {
  if(!inFarTail(i)) {
    std::cerr << "Requested to sum modes from the far tail, starting at index i=" << i << " that is out of the far tail.\n";
    throw std::runtime_error("Requested to sum modes from the far tail, starting at index that is out of the far tail.\n");
  }
  return m_c / (POW(ScalarType(i > 0 ? i : -i), m_s - 1) * (m_s - 1));
}

template <typename ScalarT, typename DimensionT>
inline void FarTail<ScalarT, DimensionT>::intersect(const FarTail<ScalarT, DimensionT>& ft2){
  if(m_s!=ft2.m_s){
    std::cerr<< "We do not allow to intersect two tails with different exponents.\n";
    throw std::runtime_error("We do not allow to intersect two tails with different exponents.\n");
  }//otherwise we intersect the far tails
  m_c=min(m_c, ft2.m_c);
}

template <typename ScalarT, class DimensionT, int D>
inline typename PolyBd<ScalarT, DimensionT, D>::ScalarType PolyBd<ScalarT, DimensionT, D>::operator()(const int i, const bool re) const {
  if(inFarTail(i))
    return m_farTail(i, 1);
  else {
    return m_nearTail(i, re);
  }
}

template <typename ScalarT, class DimensionT, int D>
inline typename PolyBd<ScalarT, DimensionT, D>::ScalarType PolyBd<ScalarT, DimensionT, D>::operator[](int i) const {
  if(i < m_nearTail.m_size)
    return m_nearTail[i];
  else
    return m_farTail[i];
}

template <typename ScalarT, class DimensionT, int D>
PolyBd<ScalarT, DimensionT, D> PolyBd<ScalarT, DimensionT, D>::operator-(const PolyBd<ScalarT, DimensionT, D>& t2) const {
  PolyBd t1 = *this;
  if(t1.getm() != t2.getm() || t1.getM() < t2.getM() || s(t1) != s(t2))
    throw std::runtime_error("Cannot subtract given two tails.");
  PolyBd result(t1.getm(), t1.getM());
  int first = 0;
  int last = this->closeTailSize();
  int i;
  for(i = first; i < last; ++i) {
    result.set(i, t1[i] - t2[i]);
  }
  setS(result, s(t1));
  setC(result, abs(C(t1) - C(t2)));
  return result;
}

template <typename ScalarT, class DimensionT, int D>
typename PolyBd<ScalarT, DimensionT, D>::ScalarType PolyBd<ScalarT, DimensionT, D>::infiniteSum(int k) const {
  ScalarType c = this->getC();
  ScalarType s = this->getS();
  ScalarType absk = (k > 0 ? k : -k);
  return 4. * c * c * (1. / (2. * s - 1.)) * power(1. / ScalarType(((this->M + absk) * (this->M))), s - 0.5) * ScalarType(-1, 1);
}

template <typename ScalarT, class DimensionT, int D>
typename PolyBd<ScalarT, DimensionT, D>::ScalarType PolyBd<ScalarT, DimensionT, D>::sum(int i, bool re) const {
  if(m_nearTail.inCloseTail(i))
    return m_nearTail.sum(i, re) + m_farTail.sum();
  else
    return m_farTail.sum(i);
}

template <typename ScalarT, class DimensionT, int D>
void PolyBd<ScalarT, DimensionT, D>::print(std::ostream& stream) const {
  int i, t;
  bool re;
  int first = 0;
  int last = this->getCloseTailSize();
  stream << "near tail (" << this->array2modeTail(first) << " < k < " << this->array2modeTail(last) << "):\n";
  for(i = first; i < last; ++i) {
    t = this->array2modeTail(i, re);
    if(re) {
      stream << "Re(a_" << t << "): " << operator[](i);
    } else {
      stream << ", Im(a_" << t << "): " << operator[](i) << "\n";
    }
  }
  stream << "far tail (k>" << t << "): \n|a_k| <= " << this->getC() << " / |k|^" << this->getS() << "\n";
  stream << "first term in the far tail: " << operator[](i) << "\n";
}

template <typename ScalarT, class DimensionT, int D>
void PolyBd<ScalarT, DimensionT, D>::printRaw(std::ostream& stream) const {
  int i;
  stream << "nearTail:\n";
  for(i = 0; i < m_nearTail.m_size; ++i) {
    stream << operator[](i) << "\n";
  }
  stream << "farTail:\n" << this->getC() << "\n" << this->getS() << "\n";
}

template <typename ScalarT, class DimensionT, int D>
inline bool PolyBd<ScalarT, DimensionT, D>::subset(const PolyBd& t2) const {
  if(s(t2) > this->getS()){
    return false;
  }
  int i;
  int max = (m_nearTail.m_size > t2.m_nearTail.m_size ? m_nearTail.m_size : t2.m_nearTail.m_size);
  for(i = 0; i <= max; ++i)
    if(!this->operator[](i).subset(t2[i])) {
//      std::cout << "!this->operator[](" << i << ") " << this->operator[](i) << " .subset(t2[" << i << "]) " << t2[i] << "\n";
      return false;
    }
  return subsetFar(t2);
}

template <typename ScalarT, class DimensionT, int D>
inline bool PolyBd<ScalarT, DimensionT, D>::subsetFar(const PolyBd<ScalarT, DimensionT, D>& t2) const {
  if(m_farTail.subset(t2.m_farTail))
    return true;
  else
    return false;
}

template <typename ScalarT, class DimensionT, int D>
inline void PolyBd<ScalarT, DimensionT, D>::closeTailIntersect(const PolyBd<ScalarT, DimensionT, D>& t2) {
  int first = 0;
  int last = this->closeTailSize();
  int i;
  ScalarType val;
  //near tail intersection
  for(i = first; i < last; ++i) {
    intersection(t2[i], this->operator[](i), val);
    this->set(i, val);
  }
}

template <typename ScalarT, class DimensionT, int D>
inline void PolyBd<ScalarT, DimensionT, D>::intersect(const PolyBd<ScalarT, DimensionT, D>& t2) {
  m_nearTail.intersect(t2.m_nearTail);
  m_farTail.intersect(t2.m_farTail);
}

template <typename ScalarT, class DimensionT, int D>
inline int PolyBd<ScalarT, DimensionT, D>::changeM(int newM, bool guard) {
  if(this->M > newM) { /// new dimension is smaller than current dimension
    ScalarType norm, ntp, max = 0;
    int i;
    for(i = this->M; i >= newM + 1; --i) {
      norm = this->template norm<ScalarType> (m_nearTail[this->mode2arrayTail(i, 1)], m_nearTail[this->mode2arrayTail(i, 0)]);
      ntp = norm * power(ScalarType(i), m_farTail.m_s);
      if(guard) {
        if(!(ntp <= max)){
          if(max < 2 * m_farTail.m_c)
            max = ntp;
          else {
            i = i + 1;
            break;
          }
        }
      } else {
        if(!(ntp <= max))
          max = ntp;
      }
    }

    newM = i;
    if(newM <= this->M) {
      //we have to estimate rest by newC/k^s
      m_nearTail.resize(this->modes2arraySize(newM - this->m), ScalarType(0.));
      m_nearTail.M = newM;
      m_farTail.setC(max);
      m_farTail.setM(newM);
      this->M = newM;
    }
    return newM;
  }
  if(this->M < newM) { /// new dimension is larger than current dimension
    m_nearTail.resize(this->modes2arraySize(newM - this->m));
    int oldM = this->M;
    m_nearTail.M = newM;
    this->M = newM;
    int i;
    for(i = oldM + 1; i <= newM; ++i) {
      this->set(i, 1, ScalarType(-m_farTail(i, 1).leftBound(), m_farTail(i, 1).rightBound()));
      this->set(i, 0, ScalarType(-m_farTail(i, 0).leftBound(), m_farTail(i, 0).rightBound()));
    }
    m_farTail.setM(newM);
  }
  return newM;
}

}
}
#endif
