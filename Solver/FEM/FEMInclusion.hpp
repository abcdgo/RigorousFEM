#ifndef _FEMINCLUSION_HPP_
#define _FEMINCLUSION_HPP_

#include "capd/dynsys/FadMap.h"

#include "../Equation/Equation.hpp"


template <typename Scalar, int D, typename Equation, class AlgorithmT>
class FEMInclusion : public capd::dynsys::FadMap<Scalar,D>, public AlgorithmT
{
    public:

        typedef typename capd::dynsys::FadMap<Scalar,D>::MatrixType MatrixType;
        typedef typename capd::dynsys::FadMap<Scalar,D>::VectorType VectorType;

        typedef Scalar ScalarType;

        typedef AlgorithmT AlgorithmType;
        typedef typename AlgorithmT::TailType TailType;
        typedef typename AlgorithmT::Class Class;


        FEM<Scalar> *fem;
        Equation *equation;

        TailType inclusionVector;


        ScalarType m_sqrt2;

        FEMInclusion(int m, int M, const ScalarType& R/**<-\nu*/, const ScalarType& c/**<initial c*/, const ScalarType& s/**<initial s*/) :
            AlgorithmType(m, M, R, ScalarType(0), s)
        {
            inclusionVector = this->initializeTail(c, s, m, M);
            m_sqrt2 = sqrt(ScalarType(2.));
        }

        void initTail()
        {
            this->m_T = this->initializeTail( 0, 0, equation->m, equation->M);
        }

        inline void setT0(const TailType& T)
        {
            inclusionVector = T;
        }

        TailType& getT0()
        {
            return inclusionVector;
        }

        template <typename VectorType>
        VectorType operator()(const VectorType& in, const TailType& T) const
        {
            return N(in, T);
        }


        template <typename VectorType>
        VectorType operator()(const VectorType& in) const
        {
            return N(in, inclusionVector);
        }

        MatrixType operator[](const VectorType& u) const
        {
            MatrixType der( u.dimension(), u.dimension() );
            computeDerivative(*this,u,der);
            return der;
        }


        //FUNCTIONS FOR DISSIPATIVE_ENCLOSURE

        int getDimension() const { return dimension(); }

        Scalar lambda( int position ) const
        {
            return equation->lambda( position );
        }

        template <typename AVector>
        AVector L(const AVector& in) const
        {
            AVector out( in.dimension() );
            return out;
        }

        template <typename AVector>
        AVector N(const AVector& in) const
        {
            return equation->N(in,inclusionVector);
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


        template <typename VectorType>
        VectorType N(const VectorType& u, const TailType& tail) const
        {
            return equation->N(u,tail);
        }
        //returns non-dissipative

        template <typename VectorType>
        TailType n(const VectorType& u, const TailType& tail) const
        {
            return equation->n(u, tail);
        }
        //returns dissipative



        //EOF FUNCTIONS FOR DISSIPATIVE_ENCLOSURE

};

#endif // _FEMINCLUSION_HPP_
