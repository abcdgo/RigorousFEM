#ifndef _DATA_HPP_
#define _DATA_HPP_


#include "../../FEM/FEM.hpp"
#include "../../Equation/BurgersNWave.hpp"

template <typename Scalar>//forward declaration because my headers are mutually inclusive (FEM and Equation)
class FEM;

template <typename Scalar>
class Data
{
    public:
        typedef Scalar ScalarType;


        std::vector<bool> dissipativeMap;
        int numberOfModes, numberOfDissipativeModes, numberOfNonDissipativeModes;
        double lambdaThreshold;
        std::vector< capd::vectalg::Matrix<Scalar,0,0> > bigSum;

        capd::vectalg::Matrix<Scalar,0,0> P, inverseOfP, diagonalMatrix, remainsFromDiagonal, matrixForNonlinearPart;
        capd::vectalg::Matrix<Scalar,0,0> NonDissipativeP, NonDissipativeDiagonalMatrix, NonDissipativeRemainsFromDiagonal, NonDissipativeMatrixForNonlinearPart;

        capd::vectalg::Vector<Scalar,0> diagonalVector, loadVector;
        capd::vectalg::Vector<Scalar,0> NonDissipativeDiagonalVector, NonDissipativeLoadVector;

        Data( FEM<Scalar> &f, double dissipativeThreshold );
        BurgersNWave<Scalar> *equation;

        void dataForDelta();
        Scalar firstIntegral(int j, int k, int i) const;
        Scalar secondIntegral(int j, int k, int i) const;

        template <typename AVector>
        AVector downsizeVector( const AVector &in ) const
        {
            AVector out( numberOfNonDissipativeModes );
            int j = 0;
            for( int i = 0; i < numberOfModes; i++ )
                if( !dissipativeMap[i] )
                {
                    out[j++] = in[i];
                }
            return out;
        }

        template <typename AVector>
        AVector upsizeVector( const AVector &in ) const
        {
            AVector out( numberOfModes );
            int j = 0;
            for( int i = 0; i < numberOfModes; i++ )
                if( !dissipativeMap[i] )
                {
                    out[i] = in[j++];
                }
                else out[i] = (Scalar)0;
            return out;
        }

        template< typename VectorType, typename TailType>
        VectorType joinHeadWithTail( const VectorType &head, const TailType &tail) const
        {
            VectorType out( numberOfModes );

            int h = 0;
            for( int i = 0; i < numberOfModes; i++ )
                if( dissipativeMap[i] ) out[i] = tail[i];
                else out[i] = head[h++];
            return out;
        }

        template< typename VectorType >
        VectorType extractHead( const VectorType &v ) const
        {
            VectorType out( numberOfNonDissipativeModes );

            int o = 0;
            for( int i = 0; i < numberOfModes; i++ )
                if( !dissipativeMap[i] ) out[o++] = v[i];
            return out;
        }

        template< typename VectorType, typename TailType >
        TailType extractTail( const VectorType &v, const TailType &t ) const
        {
            TailType out;
            out = t;
            for( int i = 0; i < numberOfModes; i++ )
                if( dissipativeMap[i] ) out.set( i, v[i] );
            return out;
        }

        bool isDissipative( int i )
        {
            return dissipativeMap[i];
        }

        Scalar lambda( int i )
        {
            return -diagonalMatrix[i][i];
        }

        bool decideIfDissipative( int i )
        {
            if( lambda(i) < -lambdaThreshold ) {return true;}
            else return false;
        }

};

template <typename Scalar>
Data<Scalar>::Data( FEM<Scalar> &f, double dissipativeThreshold )
{
        P = *(f.P);
        inverseOfP = *(f.inverseOfP);
        diagonalMatrix = *(f.diagonalMatrix);
        remainsFromDiagonal = *(f.remainsFromDiagonal);
        matrixForNonlinearPart = *(f.matrixForNonlinearPart);

        loadVector = *(f.loadVector);
        diagonalVector = *(f.diagonalVector);

        numberOfModes = f.mesh->numberOfFreeNodes();
        numberOfDissipativeModes = 0;
        numberOfNonDissipativeModes = 0;

        lambdaThreshold = dissipativeThreshold;

        for( int i = 0; i < numberOfModes; i++ )
        {
            dissipativeMap.push_back( decideIfDissipative(i) );
            if( dissipativeMap[i] == false ) numberOfNonDissipativeModes++;
            else numberOfDissipativeModes++;
        }

        NonDissipativeDiagonalVector = capd::vectalg::Vector<Scalar,0>( numberOfNonDissipativeModes );
        NonDissipativeLoadVector = capd::vectalg::Vector<Scalar,0>( numberOfNonDissipativeModes );

        int j = 0;
        for( int i = 0; i < numberOfModes; i++ )
            if( !dissipativeMap[i] )
            {
                NonDissipativeLoadVector[j] = loadVector[i];
                NonDissipativeDiagonalVector[j] = diagonalVector[i];
                j++;
            }


        NonDissipativeP = capd::vectalg::Matrix<Scalar,0,0>( numberOfNonDissipativeModes, numberOfNonDissipativeModes );
        NonDissipativeDiagonalMatrix = capd::vectalg::Matrix<Scalar,0,0>( numberOfNonDissipativeModes, numberOfNonDissipativeModes );
        NonDissipativeRemainsFromDiagonal = capd::vectalg::Matrix<Scalar,0,0>( numberOfNonDissipativeModes, numberOfNonDissipativeModes );
        NonDissipativeMatrixForNonlinearPart = capd::vectalg::Matrix<Scalar,0,0>( numberOfNonDissipativeModes, numberOfNonDissipativeModes );

        int x = 0;
        for( int i = 0; i < numberOfModes; i++ )
        {
            int y = 0;
            for( int j = 0; j < numberOfModes; j++ )
            {
                if( !dissipativeMap[i] && !dissipativeMap[j] )
                {
                    NonDissipativeP[x][y] = P[i][j];
                    NonDissipativeDiagonalMatrix[x][y] = diagonalMatrix[i][j];
                    NonDissipativeRemainsFromDiagonal[x][y] = remainsFromDiagonal[i][j];
                    NonDissipativeMatrixForNonlinearPart[x][y] = matrixForNonlinearPart[i][j];
                    y++;
                }
            }
            if( !dissipativeMap[i] ) x++;
        }
    dataForDelta();
}


template <typename Scalar>
Scalar Data<Scalar>::firstIntegral(int j, int k, int i) const
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

template <typename Scalar>
Scalar Data<Scalar>::secondIntegral(int j, int k, int i) const
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


template <typename Scalar>
void Data<Scalar>::dataForDelta()
{

    std::vector< capd::vectalg::Matrix<Scalar,0,0> > tmpBigSum;

    for( int i = 0; i < numberOfModes; i++ )
    {

        tmpBigSum.push_back( capd::vectalg::Matrix<Scalar,0,0>( numberOfModes, numberOfModes ) );
        bigSum.push_back( capd::vectalg::Matrix<Scalar,0,0>( numberOfModes, numberOfModes ) );

        for( int l = numberOfNonDissipativeModes; l < numberOfModes; l++ )
        {
            for( int n = 0; n < numberOfModes; n++)
            {
                for( int k = 0; k < numberOfModes; k++ )
                {
                    for( int j = 0; j < numberOfModes; j++ )
                    {
                        if( n < numberOfNonDissipativeModes ) (tmpBigSum[i])[l][n] += P[j][l] * P[k][n] * firstIntegral(j,k,i);// times first integral
                        else (tmpBigSum[i])[l][n] += P[j][l] * P[k][n] * secondIntegral(j,k,i);//times second intergral
                    }
                }
            }
        }
    }
    
    for( int i = 0; i < numberOfModes; i++ )
        for( int j = 0; j < numberOfModes; j++ )
        {
            capd::vectalg::Vector<Scalar,0> tmpVector = capd::vectalg::Vector<Scalar,0>( numberOfModes );

            for( int m = 0; m < numberOfModes; m++ )
            {
                tmpVector[m] = (tmpBigSum[m])[i][j];
            }

            tmpVector = matrixForNonlinearPart * tmpVector;

            for( int m = 0; m < numberOfModes; m++ )
            {
                (bigSum[m])[i][j] = tmpVector[m];
            }
        }
}

#endif // _DATA_HPP_
