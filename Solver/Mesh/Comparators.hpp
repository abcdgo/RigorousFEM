#ifndef _COMPARATORS_HPP_
#define _COMPARATORS_HPP_

#include "Node.hpp"
#include "Edge.hpp"

template <typename T>
class compareNodeByNumber
{
    public:
        bool operator() ( Node<T> *n1, Node<T> *n2 )
        {
            return n1->getNumber() < n2->getNumber();
        }
};

template <typename T>
class isNodeEqual
{
    public:
        bool operator() ( Node<T> *n1, Node<T> *n2 )
        {
            return n1->getNumber() == n2->getNumber();
        }
};

template <typename T>
class compareNodeByCoordinates
{
    public:
        bool operator() ( Node<T> *n1, Node<T> *n2 )
        {
            for( int i = 0; i < n1->getDimension(); i++ )
            {
                if( n1->getCoordinate(i) < n2->getCoordinate(i) ) return true;
                if( n1->getCoordinate(i) > n2->getCoordinate(i) ) return false;
            }
            std::cout << "compareNodeByCoordinate OPS\n";
            return true;
        }
};

template <typename T>
class compareEdgeByNumber
{
    public:
        bool operator() ( Edge<T> *e1, Edge<T> *e2 )
        {
            return e1->getNumber() < e2->getNumber();
        }
};

#endif // _COMPARATORS_HPP_
