#ifndef _DET_HPP_
#define _DET_HPP_

#include "D1Real.h"

template <typename Scalar>
class Det : public capd::jaco::D1Real<Scalar>
{
    public:
        Det( int m, int M, Scalar ni, Scalar a0) : capd::jaco::D1Real<Scalar>(m,M)
        {
            m_ni = ni;
            m_a0 = a0;
        }

        Det( Scalar ni, Scalar a0) : capd::jaco::D1Real<Scalar>()//D1Real will take proper m & M from config.h
        {
            m_ni = ni;
            m_a0 = a0;
        }

        Scalar m_ni; ///-\nu,
        Scalar m_a0; ///constant value of the Re{a0},
        Scalar m_p; ///order of the derivative in the linear part,
        //T m_d; ///space dimension,
        Scalar m_r; ///order of the derivative in the nonlinear part,
        Scalar m_N_coeff; ///coefficient in front of non-linear term m_d u_k = -m_N_coeff(u)^2_x + \nu u_xx,
        Scalar m_sufficientlyLarge; ///d+r+1, minimal order of decrase in the far tail at which solutions exists,
                               ///see Theorem 11 in [Z3].


        int maximumPoint( const Scalar& h, const Scalar& r, int M ) const
        {  return static_cast<int>( round( sqrt( ( r.rightBound()/(-2.*this->m_ni*h) ).rightBound() ) ) );  }

};

#endif // _DET_HPP_
