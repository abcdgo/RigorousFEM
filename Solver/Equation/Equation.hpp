#ifndef _EQUATION_HPP_
#define _EQUATION_HPP_


#include "capd/vectalg/Vector.hpp"

#include "../Mesh/Node.hpp"
#include "../Mesh/Element.hpp"

template <typename T>//forward declaration because my headers are mutually inclusive (FEM and Equation)
class FEM;

template <typename T>
class Equation
{
    public:
        FEM<T> *fem;
        T initialVectorHull;
        double ni;
        double nonAutonomousParameter;

        //Auxilary functions for FEM class
        virtual T functionG( Node<T> *node, double time = 0 ) const = 0;
        virtual T functionF( Node<T> *node, double time = 0 ) = 0;
        virtual T timeDerivativeOfFunctionG( Node<T> *node, double time ) = 0;
        virtual capd::vectalg::Vector<T,0>* gradientOfFunctionG( Element<T> *element, std::list<Node<T>*> &elementNodes, int dimension, double time = 0 );

        //Auxilary functions for Numerics class
        virtual void setInitialVector( capd::vectalg::Vector<T,0> *u_0, std::list< Node<T>* > *listOfNodes = NULL ) = 0;
};


template <typename T>
capd::vectalg::Vector<T,0>* Equation<T>::gradientOfFunctionG( Element<T> *element, std::list<Node<T>*> &elementNodes, int dimension, double time )//compute the gradient G on simplex p.170
{
    capd::vectalg::Vector<T,0> *gradient = new capd::vectalg::Vector<T,0>( dimension );
    capd::vectalg::Vector<T,0> *nodalValues = new capd::vectalg::Vector<T,0>( dimension + 1 );

    int i = 0;
    for( typename std::list<Node<T>*>::iterator it = elementNodes.begin(); it != elementNodes.end(); it++,i++ )
    {
        (*nodalValues)[i] = functionG( *it, time );
    }

    for( int column = 0; column < dimension + 1; column++ )
    {
        for( int row = 0; row < dimension ; row++ )
        {
            (*gradient)[row] += (*nodalValues)[column] * (*element->inverseOfMatrixM)[row+1][column];
        }
    }
    delete nodalValues;
    return gradient;
}

#endif // _EQUATION_HPP_
