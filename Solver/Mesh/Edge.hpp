#ifndef _EDGE_HPP_
#define _EDGE_HPP_

#include <iostream>

#include "Node.hpp"

template <typename T>
class Edge
{
    public:
        Edge(Node<T> *nodeStart, Node<T> *nodeEnd, bool isFree = true);
        Edge(const Edge<T> &edge);

        bool isFree();
        bool isConstrained();

        int getNumber(){return edgeNumber;}

        static void resetCounter(){edgeCounter = 1;}
        static void reverseCounter(int i = -1){edgeCounter+=i;}

        void printMe();
    //protected:
    //private:
        Node<T> *nodeStart, *nodeEnd;

        int edgeNumber;
        static int edgeCounter;

        bool isFreeEdge;
};

template <typename T>
int Edge<T>::edgeCounter = 1;

template <typename T>
Edge<T>::Edge(Node<T>* nodeStart, Node<T>* nodeEnd, bool isFree)
{
    this->nodeStart = nodeStart;
    this->nodeEnd = nodeEnd;
    edgeNumber = edgeCounter++;
    isFreeEdge = isFree;
}

template <typename T>
Edge<T>::Edge(const Edge<T>& edge)
{
    nodeStart = edge.nodeStart;
    nodeEnd = edge.nodeEnd;

    isFreeEdge = edge.isFree();
    edgeNumber = edge.getNumber();
}


template <typename T>
bool Edge<T>::isFree()
{
    if( isFreeEdge ) return true;
    else return false;
}

template <typename T>
bool Edge<T>::isConstrained()
{
    if( !isFreeEdge ) return true;
    else return false;
}



template <typename T>
void Edge<T>::printMe()
{
    std::cout << "Krawedz nr: " << edgeNumber << " jest";
    if( isFree() ) std::cout << " wolna i";
    else std::cout << " constrained i";
    std::cout << " laczy wierzcholek: " << std::endl;
    std::cout << "a)";
    nodeStart->printMe();
    std::cout << "z" << std::endl;
    std::cout << "b)";
    nodeEnd->printMe();
    std::cout << std::endl;
}

#endif // _EDGE_HPP_
