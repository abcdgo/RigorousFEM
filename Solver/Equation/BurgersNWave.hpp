#ifndef _BURGERSNWAVE_HPP_
#define _BURGERSNWAVE_HPP_

#include <time.h>

#include "Equation.hpp"
#include "../Mesh/Node.hpp"

template <typename T>
class BurgersNWave : public Equation<T>
{
    public:

        typedef T ScalarType;
        typedef capd::vectalg::Vector<ScalarType,0> VectorType;

        VectorType NonautonomousVector;

        T functionG( Node<T> *node, double time = 0 ) const;
        T functionF( Node<T> *node, double time = 0 );
        T timeDerivativeOfFunctionG( Node<T> *node, double time = 0 );

        BurgersNWave( double hull, double ni, double nonAutonomousParameter )
        {
            this->initialVectorHull = T(-hull,hull);
            this->ni = ni;
            this->nonAutonomousParameter = nonAutonomousParameter;
        }

        void setInitialVector( capd::vectalg::Vector<T,0> *u_0, std::list< Node<T>* > *listOfNodes );
        void init( capd::vectalg::Vector<T,0> *u_0, std::list< Node<T>* > *listOfNodes );
        void candidate( capd::vectalg::Vector<T,0> *u_0, std::string &inStr );

        template <typename ScalarType>
        void assembleNonLinearVector( capd::vectalg::Vector<ScalarType,0> *nonLinearVector, const capd::vectalg::Vector<ScalarType,0> *previousU) const
        {
            for( typename std::list< Element<T>* >::iterator elementsIterator = this->fem->mesh->listOfElements.begin(); elementsIterator != this->fem->mesh->listOfElements.end(); elementsIterator++ )
            {
                Node<T> *n1,*n2;
                //n1 = n2 = NULL;
                //TODO: ONLY IN 1D - more d's then more n1,n2,n3,...
                for( typename std::list< EdgeWithOrientation<T>*>::iterator edgeIterator = (*elementsIterator)->listOfEdges.begin(); edgeIterator != (*elementsIterator)->listOfEdges.end(); edgeIterator++ )
                {
                    n1 = (*edgeIterator)->edge->nodeStart;
                    n2 = (*edgeIterator)->edge->nodeEnd;
                }

                if( n1->getCoordinate(0) > n2->getCoordinate(0) )
                {
                    Node<T> *swap;
                    swap = n1;
                    n1 = n2;
                    n2 = swap;
                }

                //y = a(i)x + b(i)
                ScalarType a1,/*a2,*/a3,b1,b2,b3;

                T x1,x2;//helpful for translation/normalization
                x1 = 0;
                x2 = n2->getCoordinate(0) - n1->getCoordinate(0);

                ScalarType alpha_1;
                ScalarType alpha_2;

                if( n1->isFree() ) alpha_1 = (*previousU)[ (*(this->fem)->cutCoordinateVector)[ n1->getNumber() ] ];
                else alpha_1 = functionG( n1 );

                if( n2->isFree() ) alpha_2 = (*previousU)[ (*(this->fem)->cutCoordinateVector)[ n2->getNumber() ] ];
                else alpha_2 = functionG( n2 );

                //a1 = ( alpha_2 - alpha_1 ) / ( x2 - x1 );//optimalization line below
                a1 = ( alpha_2 - alpha_1 ) / x2 ;//we do translation therefore x1 always equals zero
                b1 = alpha_1;

                b2 = a1;// if y = a1x+b1 then y'= a1 therefore b2 = a1
                //a2 = 0;

                //let's assume k equals n1->getNumber
                if( n1->isFree() )
                {
                    a3 = -1/x2;
                    b3 = 1;

                    //product of this 3 functions -> Ax^2 + Bx + C
                    ScalarType A = a1 * b2 * a3;
                    ScalarType B = a1 * b2 * b3 + b1 * b2 * a3;
                    ScalarType C = b1 * b2 * b3;

                    (*nonLinearVector)[ (*(this->fem)->cutCoordinateVector)[ n1->getNumber() ] ] += (A/3) * capd::intervals::power(x2, 3) + (B/2) * capd::intervals::power(x2, 2) + C*x2;
                }

                //let's assume k equals n2->getNumber
                if( n2->isFree() )
                {
                    a3 = 1/x2;
                    b3 = 0;

                    //product of this 3 functions -> Ax^2 + Bx + C
                    ScalarType A = a1 * b2 * a3;
                    ScalarType B = /*a1 * b2 * b3 + */ b1 * b2 * a3;

                    (*nonLinearVector)[ (*(this->fem)->cutCoordinateVector)[ n2->getNumber() ] ] += (A/3) * capd::intervals::power(x2, 3) + (B/2) * capd::intervals::power(x2, 2) ;/*+ C*x2;*/
                }
            }
        }


        template <typename TimeT, typename VectorType>
        void assembleNonautonomousPart( const TimeT &time, VectorType &in ) const
        {
            for( int i = 0; i < in.dimension(); i++ )
            {
                in[i] = NonautonomousVector[i] * cos(time);
            }
        }

        //template <typename VectorType>
        void assembleNonautonomousVector()
        {

            VectorType out( this->fem->mesh->numberOfFreeNodes() );

            for( typename std::list< Element<T>* >::iterator elementsIterator = this->fem->mesh->listOfElements.begin(); elementsIterator != this->fem->mesh->listOfElements.end(); elementsIterator++ )
            {
                Node<T> *n1,*n2;
                for( typename std::list< EdgeWithOrientation<T>*>::iterator edgeIterator = (*elementsIterator)->listOfEdges.begin(); edgeIterator != (*elementsIterator)->listOfEdges.end(); edgeIterator++ )
                {
                    n1 = (*edgeIterator)->edge->nodeStart;
                    n2 = (*edgeIterator)->edge->nodeEnd;
                }

                if( n1->getCoordinate(0) > n2->getCoordinate(0) )
                {
                    Node<T> *swap;
                    swap = n1;
                    n1 = n2;
                    n2 = swap;
                }

                ScalarType h;//helpful for translation/normalization
                h = n2->getCoordinate(0) - n1->getCoordinate(0);

                if( n1->isFree() )
                {
                    out[ (*(this->fem)->cutCoordinateVector)[ n1->getNumber() ] ] += capd::intervals::cos( n1->getCoordinate(0) ) + ( capd::intervals::sin( n1->getCoordinate(0) ) - capd::intervals::sin( n2->getCoordinate(0) ) ) / h;
                }

                if( n2->isFree() )
                {
                    out[ (*(this->fem)->cutCoordinateVector)[ n2->getNumber() ] ] += -1 * capd::intervals::cos( n2->getCoordinate(0) ) + ( capd::intervals::sin( n2->getCoordinate(0) ) - capd::intervals::sin( n1->getCoordinate(0) ) ) / h;
                }
            }
            out *= this->nonAutonomousParameter;
            NonautonomousVector = out;
        }


        void assembleNonautonomousVectorType2()
        {

            VectorType out( this->fem->mesh->numberOfFreeNodes() );

            for( typename std::list< Element<T>* >::iterator elementsIterator = this->fem->mesh->listOfElements.begin(); elementsIterator != this->fem->mesh->listOfElements.end(); elementsIterator++ )
            {
                Node<T> *n1,*n2;
                for( typename std::list< EdgeWithOrientation<T>*>::iterator edgeIterator = (*elementsIterator)->listOfEdges.begin(); edgeIterator != (*elementsIterator)->listOfEdges.end(); edgeIterator++ )
                {
                    n1 = (*edgeIterator)->edge->nodeStart;
                    n2 = (*edgeIterator)->edge->nodeEnd;
                }

                if( n1->getCoordinate(0) > n2->getCoordinate(0) )
                {
                    Node<T> *swap;
                    swap = n1;
                    n1 = n2;
                    n2 = swap;
                }

                ScalarType h;//helpful for translation/normalization
                h = n2->getCoordinate(0) - n1->getCoordinate(0);

                if( n1->isFree() )
                {
                    out[ (*(this->fem)->cutCoordinateVector)[ n1->getNumber() ] ] += 3.0 * capd::intervals::cos( n1->getCoordinate(0) ) + ( 5.0 * capd::intervals::sin( n1->getCoordinate(0) ) - capd::intervals::sin( n2->getCoordinate(0) ) - 4.0 * capd::intervals::sin( h/2 + n1->getCoordinate(0) ) ) / h;
                }

                if( n2->isFree() )
                {
                    out[ (*(this->fem)->cutCoordinateVector)[ n2->getNumber() ] ] += -1.0 * capd::intervals::cos( n2->getCoordinate(0) ) - 2.0 * capd::intervals::cos( h/2 + n1->getCoordinate(0) ) + ( capd::intervals::sin( n2->getCoordinate(0) ) - 5.0 * capd::intervals::sin( n1->getCoordinate(0) ) + 4.0 * capd::intervals::sin( h/2 + n1->getCoordinate(0) ) ) / h;
                }
            }
            out *= this->nonAutonomousParameter;
            NonautonomousVector = out;
        }

};

template <typename T>
T BurgersNWave<T>::functionG( Node<T> *node, double time ) const
{
    if( node->isFree() ) return 0;
    return 0;
}

template <typename T>
T BurgersNWave<T>::functionF( Node<T> *node, double time )
{
    return 0;
}

template <typename T>
T BurgersNWave<T>::timeDerivativeOfFunctionG( Node<T> *node, double time )
{
    return 0;
}

//Initial vector with hull
template <typename T>
void BurgersNWave<T>::setInitialVector( capd::vectalg::Vector<T,0> *u_0, std::list< Node<T>* > *listOfNodes )
{
    int a = 1, i = 0;
    for( typename std::list<Node<T>*>::iterator nodeIterator = listOfNodes->begin(); nodeIterator != listOfNodes->end(); nodeIterator++ )
    {
        if( (*nodeIterator)->isFree() )
        {
            T ex = ( (*nodeIterator)->getCoordinate(0) * (*nodeIterator)->getCoordinate(0) ) / -4;
            (*u_0)[i++] = (*nodeIterator)->getCoordinate(0) * ( sqrt(a) * exp(ex) / ( 1 + sqrt(a) * exp(ex) ) ) + this->initialVectorHull;
        }
    }
}

//Initial vector without hull
template <typename T>
void BurgersNWave<T>::init( capd::vectalg::Vector<T,0> *u_0, std::list< Node<T>* > *listOfNodes )
{
    int a = 1, i = 0;
    for( typename std::list<Node<T>*>::iterator nodeIterator = listOfNodes->begin(); nodeIterator != listOfNodes->end(); nodeIterator++ )
    {
        if( (*nodeIterator)->isFree() )
        {
            T ex = ( (*nodeIterator)->getCoordinate(0) * (*nodeIterator)->getCoordinate(0) ) / -4;
            (*u_0)[i++] = (*nodeIterator)->getCoordinate(0) * ( sqrt(a) * exp(ex) / ( 1 + sqrt(a) * exp(ex) ) );
        }
    }
}

//Initial vector made from candidade typed from simulation on doubles
template <typename T>
void BurgersNWave<T>::candidate( capd::vectalg::Vector<T,0> *u_0, std::string &inStr )
{
    std::istringstream in(inStr);
    in >> *u_0;
}

#endif // _BURGERSNWAVE_HPP_



