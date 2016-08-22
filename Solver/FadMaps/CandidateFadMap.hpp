#ifndef _FEMNONLINEARFADMAP_HPP_
#define _FEMNONLINEARFADMAP_HPP_

#include "capd/dynsys/FadMap.h"

#include "../Equation/Equation.hpp"

template <typename Scalar, int D, typename Equation>
class CandidateFadMap : public capd::dynsys::FadMap<Scalar,D>
{
    public:

        typedef typename capd::dynsys::FadMap<Scalar,D>::MatrixType MatrixType;
        typedef typename capd::dynsys::FadMap<Scalar,D>::VectorType VectorType;

        FEM<Scalar> *fem;
        Equation *equation;


        template <typename TimeT, typename AVector>
        AVector operator()(const TimeT& t, const AVector& in) const
        {
            //Au' + Ku = F
            AVector out( fem->mesh->numberOfFreeNodes() );

            AVector nonLinearVector( fem->mesh->numberOfFreeNodes() );
            equation->assembleNonLinearVector( &nonLinearVector, &in );

            AVector nonAutonomousVector( fem->mesh->numberOfFreeNodes() );
            equation->assembleNonautonomousPart( t, nonAutonomousVector );

            out = capd::vectalg::subtractObjects<AVector>(
                    capd::vectalg::matrixByVector<AVector>( (*fem->inverseOfMassMatrix),
                                                   capd::vectalg::subtractObjects<AVector>( capd::vectalg::subtractObjects<AVector>( (*fem->loadVector), nonLinearVector), nonAutonomousVector ) ),

                    capd::vectalg::matrixByVector<AVector>( (*fem->inverseOfMassMatrixTimesStiffnessMatrix), in )
                    );

            return out;
        }

        MatrixType operator[](const VectorType& u) const
        {
            MatrixType der( u.size(), u.size() );
            computeDerivative(*this,capd::TypeTraits<Scalar>::zero(),u,der);
            return der;
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

        int dimension() const {return fem->mesh->numberOfFreeNodes();}
};

#endif // _FEMNONLINEARFADMAP_HPP_
