#ifndef _PERTURBATION_HPP_
#define _PERTURBATION_HPP_

#include "capd/dynsys/FadMap.h"


template <typename Scalar, int D, class AlgorithmT> 
class Perturbation : public capd::dynsys::FadMap<Scalar,D>, public AlgorithmT
{
    public:

        typedef typename capd::dynsys::FadMap<Scalar,D>::MatrixType MatrixType;
        typedef typename capd::dynsys::FadMap<Scalar,D>::VectorType VectorType;

        typedef Scalar ScalarType;

        typedef AlgorithmT AlgorithmType;
        typedef typename AlgorithmT::TailType TailType;
        typedef typename AlgorithmT::Class Class;


        TailType inclusionVector;


        ScalarType m_sqrt2;

        Perturbation(int m, int M, const ScalarType& R/**<-\nu*/, const ScalarType& c/**<initial c*/, const ScalarType& s/**<initial s*/, Data<Scalar> &d) :
            AlgorithmType(m, M, R, ScalarType(0), s, d)
        {
            inclusionVector = this->initializeTail(c, s, m, M);
            m_sqrt2 = sqrt(ScalarType(2.));
        }

        inline void setT0(const TailType& T){ inclusionVector = T; }

        TailType& getT0(){ return inclusionVector; }

        template <typename TimeT, typename VectorType>
        VectorType operator()( const TimeT& t, const VectorType& in, const TailType& T) const 
        {
            return N(t, in, T);
        }

        template <class TimeT, typename VectorType>
        VectorType operator()(const TimeT& t, const VectorType& in) const
        {
            return N(t, in, inclusionVector);
        }


        /// computes simultaneously value and derivative of the map for a given vector and time
        VectorType operator()(Scalar t, const VectorType& x, MatrixType& o_der) const
        {
            return computeDerivative(*this,t,x,o_der);
        }

        MatrixType derivative(const Scalar& t, const VectorType& u) const
        {
            MatrixType der( u.size(), u.size() );
            computeDerivative(*this,t,u,der);
            return der;
        }

        MatrixType operator[](const VectorType& u) const
        {
            MatrixType der( u.dimension(), u.dimension() );
            computeDerivative(*this,capd::TypeTraits<Scalar>::zero(),u,der);
            return der;
        }

        int dimension() const { return this->data.numberOfDissipativeModes; }


        //FUNCTIONS FOR DISSIPATIVE_ENCLOSURE

        int getDimension() const { return dimension(); }

        Scalar lambda( int position ) const
        {
            return AlgorithmT::lambda( position );
        }

        template <typename TimeT, typename AVector>
        AVector L(const TimeT &time,const AVector& in) const
        {
            AVector out( in.dimension() );
            return out;
        }

        template <typename TimeT, typename AVector>
        AVector N(const TimeT &time, const AVector& in) const
        {
            return N(time, in,inclusionVector);
        }

        template <typename TimeT, typename VectorType>
        VectorType N(const TimeT &time, const VectorType& u, const TailType& tail) const
        {
            return AlgorithmT::N(time,u,tail);
        }//returns non-dissipative mods

        template <typename TimeT,typename VectorType>
        TailType n(const TimeT &time, const VectorType& u, const TailType& tail) const
        {
            std::cout << "Perturbation.hpp n()" << std::endl;
            return AlgorithmT::n(time, u,tail);
        }//returns dissipative mods



        template <typename TimeT, typename AVector>
        AVector perturbations(const TimeT &time, const AVector& in) const
        {
            return perturbations(time, in,inclusionVector);
        }

        template <typename TimeT, typename VectorType>
        VectorType perturbations(const TimeT &time, const VectorType& u, const TailType& tail) const 
        {
            return AlgorithmT::perturbations(time,u,tail);
        }

        ///Calculates g_k^{+-} = (T0_k^{+-}-b_k^{+-})e^{-\lambda_k h}+b_k^{+-}, version used by enclosure algorithms.
        ///Returns either Re{a_k} or Im{a_k}, depending on which value is at index i in the modes storage.
        inline Scalar g_encl(const Scalar& h, int i, const Scalar& x, const Scalar& b)
        {
            Scalar g(0);
            Scalar ex = exp(h * lambda(i));
            g.setLeftBound((((x.leftBound() - b.leftBound()) * ex) + b.leftBound()).leftBound());
            g.setRightBound((((x.rightBound() - b.rightBound()) * ex) + b.rightBound()).rightBound());
            return g;
        }

        //EOF FUNCTIONS FOR DISSIPATIVE_ENCLOSURE
};

#endif // _PERTURBATION_HPP_
