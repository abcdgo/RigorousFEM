#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include <iostream>

#include "capd/vectalg/Matrix.hpp"
#include "capd/vectalg/Vector.hpp"

#include "Node.hpp"
#include "Edge.hpp"
#include "Comparators.hpp"


#include <math.h>


template <typename T>
struct EdgeWithOrientation
{
    Edge<T> *edge;
    bool isCounterclockwise;

    EdgeWithOrientation(){edge=NULL;}
    EdgeWithOrientation(const EdgeWithOrientation<T> &e){ edge=e.edge; isCounterclockwise=e.isCounterclockwise;}
};


template <typename T>
class Element
{
    public:
        Element( int dimension ){ matrixM = new capd::vectalg::Matrix<T,0,0>(dimension+1,dimension+1);
                                  inverseOfMatrixM = new capd::vectalg::Matrix<T,0,0>(dimension+1,dimension+1);
                                  this->dimension = dimension;}
        ~Element();
        
        void addEdge( Edge<T> *e, bool orientation = true);//elements will by defined by edges counterclocwise oriented (p.136)
        
        T areaOfSimplex( std::list<Node<T>*> &elementNodes );
        Node<T>* centroidOfSimplex( std::list<Node<T>*> &elementNodes );

        T dotProductOfGradients();
        
        void getElementNodes( std::list<Node<T>*> &elementNodes );//changes input (elementNodes)
        
        void printMe();
        
        int dimension;
        std::list< EdgeWithOrientation<T>* > listOfEdges;//each element is defined by some edges //each triangle is defined by 3 edges
        capd::vectalg::Matrix<T,0,0> *matrixM, *inverseOfMatrixM;
};

template <typename T>
Element<T>::~Element()
{
    delete matrixM;
    delete inverseOfMatrixM;
    matrixM = inverseOfMatrixM = NULL;
    for( typename std::list< EdgeWithOrientation<T>*>::iterator it = listOfEdges.begin(); it != listOfEdges.end(); it++ ) delete *it;
}

template <typename T>
void Element<T>::addEdge(Edge<T> *e, bool orientation)
{
    EdgeWithOrientation<T> *ewo = new EdgeWithOrientation<T>; //TODO: deleting this elements (memory leak)

    ewo->edge = e;
    ewo->isCounterclockwise = orientation;

    listOfEdges.push_back(ewo);
}

template <typename T>
void Element<T>::getElementNodes( std::list<Node<T>*> &elementNodes )
{
    for( typename std::list< EdgeWithOrientation<T>*>::iterator edgeIterator = listOfEdges.begin(); edgeIterator != listOfEdges.end(); edgeIterator++ )
    {
        elementNodes.push_back( (*edgeIterator)->edge->nodeStart );
        elementNodes.push_back( (*edgeIterator)->edge->nodeEnd );
    }
    elementNodes.sort( compareNodeByNumber<T>() );
    elementNodes.unique( isNodeEqual<T>() );
}


template <typename T>
T Element<T>::areaOfSimplex( std::list<Node<T>*> &elementNodes )
{
    capd::vectalg::Matrix<T,0,0> *matrix = new capd::vectalg::Matrix<T,0,0>( dimension, dimension );
    typename std::list<Node<T>*>::iterator it, first;
    
    first = elementNodes.begin();
    int column = 0;
    for( it = elementNodes.begin(), it++; it != elementNodes.end(); it++ )
    {
        for( int row = 0; row < (*it)->getDimension(); row++ )
        {
            (*matrix)[row][column] = (*it)->getCoordinate(row) - (*first)->getCoordinate(row);
        }
        column++;
    }

    T areaOfSimplex = capd::abs( capd::matrixAlgorithms::det(*matrix) ) / 2;
    delete matrix;
    return areaOfSimplex;
}

template <typename T>
Node<T>* Element<T>::centroidOfSimplex( std::list<Node<T>*> &elementNodes )
{
    Node<T> *centroid = new Node<T>( dimension );
    centroid->reverseCounter();
    for( int i = 0; i < dimension; i++ ) centroid->setIthCoordinate(i, 0);
    
    for( typename std::list<Node<T>*>::iterator it = elementNodes.begin(); it != elementNodes.end(); it++ )
    {
        for( int i = 0; i < dimension; i++ ) 
        {
            T tmp = (*it)->getCoordinate(i) / elementNodes.size();
            centroid->setIthCoordinate(i, centroid->getCoordinate(i) + tmp );
        }
    }
    
    return centroid;
}

template <typename T>
T Element<T>::dotProductOfGradients()
{
    T gradientDotProduct = 0;
    
    for( int row = 1; row < dimension + 1; row++ )
    {
        for( int column = 0; column < dimension + 1; column++ )
        {
            for( int index = 0; index < column+1; index++ )
            {
                gradientDotProduct += (*inverseOfMatrixM)[row][column] * (*inverseOfMatrixM)[row][index];
            }
        }
    }
    return gradientDotProduct;
}

template <typename T>
void Element<T>::printMe()
{
    std::cout << "Jestem elementem skladajacym sie z krawedzi:" << std::endl;

    typename std::list< EdgeWithOrientation<T>* >::iterator it;
    
    for( it = listOfEdges.begin(); it != listOfEdges.end(); it++ )
    {
        if( !(*it)->isCounterclockwise ) std::cout << "Na odwrot" << std::endl;
        (*it)->edge->printMe();
    }
}

#endif // _ELEMENT_HPP_
