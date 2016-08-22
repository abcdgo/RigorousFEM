#ifndef _HEAD_HPP_
#define _HEAD_HPP_

#include "capd/dynsys/FadMap.h"

template <typename Scalar, int D , class AlgorithmT>
class Head : public capd::dynsys::FadMap<Scalar,D>, public AlgorithmT
{

    public:

        typedef typename capd::dynsys::FadMap<Scalar,D>::MatrixType MatrixType;
        typedef typename capd::dynsys::FadMap<Scalar,D>::VectorType VectorType;
        typedef Scalar ScalarType;

        typedef typename Scalar::BoundType DoubleType;//if scalar has got no bound type or is a basic type (ie. double, int, etc...) then you get error here - needed for Jaco's implementation of dissipative inclusion

        VectorType m_y_c;

        Head(int m, int M, const ScalarType& R/**<-\nu*/, const ScalarType& c/**<initial c*/, const ScalarType& s/**<initial s*/, Data<Scalar> &d) :
            AlgorithmT(m, M, R, ScalarType(0), s, d), m_y_c( modes2arraySize(m) )
        {
        }


        inline int modes2arraySize( int modes ) const
        {
            return modes;
        }

        template <typename TimeT, typename AVector>
        AVector operator()(const TimeT& t, const AVector& in) const
        {
            return L(t,in) + N(t,in);
        }

        // computes simultaneously value and derivative of the map for a given vector and time
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
            MatrixType der( u.size(), u.size() );//der( wym(u),wym(u) )
            computeDerivative(*this,capd::TypeTraits<Scalar>::zero(),u,der);
            return der;
        }


        int dimension() const { return this->data.numberOfNonDissipativeModes;}


        //FUNCTIONS FOR DISSIPATIVE_ENCLOSURE

        int getDimension() const { return dimension(); }

        Scalar lambda( int position ) const
        {
            return AlgorithmT::lambda( position );
        }

        bool isDissipative(int position)
        {
            return AlgorithmT::isDissipative( position );
        }

        template <typename TimeT,typename AVector>
        AVector L(const TimeT &time, const AVector& in) const
        {
            return (-1) * capd::vectalg::matrixByVector<AVector>( this->data.NonDissipativeDiagonalMatrix, in );
        }

        template <typename TimeT, typename AVector>
        AVector N(const TimeT &time, const AVector& in) const
        {
            return AlgorithmT::N(time,in);
        }


        ///sets y_c
        inline void setYc(const VectorType& y_c)
        {
            m_y_c=y_c;
        }

        ///y_c=0
        inline void eraseYc()
        {
            m_y_c=0;
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

        std::string array2modeIndex(int i) {return "InclusionHead array2modeIndex";}
        //EOF FUNCTIONS FOR DISSIPATIVE_ENCLOSURE

};
#endif // _HEAD_HPP_
