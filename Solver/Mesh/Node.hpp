#ifndef _NODE_HPP_
#define _NODE_HPP_

#include <iostream>


template <typename T>
class Node
{
    public:
        
        Node(int dimension, bool isFree = true);
        Node(const Node<T> &node);
        ~Node();

        T getCoordinate(int index);
        void setIthCoordinate(int index, T value);

        int getDimension();

        bool isFree();
        bool isConstrained();

        void setFree(){isFreeNode = true;}
        void setConstrained(){isFreeNode = false;}

        bool isSame(Node<T> &node);

        int getNumber(){return nodeNumber;}
        void setNumber( int num ){ nodeNumber = num; }
        static void reverseCounter(int i = -1){nodeCounter+=i;}

        void printMe();//TODO: it's temporary function - later it should be removed
        
    protected:
    private:
        int dimension;
        T* coordinates;

        bool isFreeNode;

        int nodeNumber;
        static int nodeCounter;

        //Node(){coordinates = NULL;nodeNumber = nodeCounter++;isFreeNode=true;}
        Node();
};


template <typename T>
int Node<T>::nodeCounter = 1;

template <typename T>
Node<T>::Node(int dimension, bool isFree)
{
    coordinates = new T[dimension];
    this->dimension = dimension;
    nodeNumber = nodeCounter++;
    isFreeNode = isFree;
}

template <typename T>
Node<T>::Node(const Node<T> &node)
{
    dimension = node.dimension;
    coordinates = new T[node.dimension];
    isFreeNode = node.isFree();
    nodeNumber = node.getNumber();
    for( int i = 0; i < node.dimension; i++) coordinates[i] = node.coordinates[i];
}

template <typename T>
Node<T>::~Node()
{
    if( coordinates != NULL ) delete[] coordinates;
}

template <typename T>
T Node<T>::getCoordinate(int index)
{
    //TODO:
    //if( coordinates != NULL && index < this->dimension && index >= 0)
    return coordinates[index];
}

template <typename T>
void Node<T>::setIthCoordinate(int index, T value)
{
    //TODO:
    //if( index < dimension && index >= 0 )
    coordinates[index] = value;
}

template <typename T>
int Node<T>::getDimension()
{
    return this->dimension;
}


template <typename T>
bool Node<T>::isFree()
{
    if( isFreeNode ) return true;
    else return false;
}

template <typename T>
bool Node<T>::isConstrained()
{
    if( !isFreeNode ) return true;
    else return false;
}

template <typename T>
bool Node<T>::isSame(Node<T>& node)
{
    if( this->dimension != node.getDimension() ) return false;

    for( int i = 0; i < node.getDimension(); i++ )
    {
        if( getCoordinate(i) != node.getCoordinate(i) ) return false;
    }
    return true;

}

template <typename T>
void Node<T>::printMe()
{
    std::cout << "Jestem " << dimension << " wymiarowym wierzcholkiem o numerze: " << nodeNumber << " o wspl: [";
    for(int i = 0; i < dimension; i++) std::cout << coordinates[i] << " ";
    std::cout << "]";
    if( isFree() ) std::cout << " wolny" << std::endl;
    else std::cout << " constrained" << std::endl;
}

#endif // _NODE_HPP_
