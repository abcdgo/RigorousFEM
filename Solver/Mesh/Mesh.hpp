#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <algorithm>
#include <fstream>

#include "Node.hpp"
#include "Edge.hpp"
#include "Element.hpp"
#include "Comparators.hpp"
#include "MeshExceptions.hpp"

//TODO: it would be nice if Mesh had some iterators (simplex iteratr, node iterator, etc.)



template <typename T>
class Mesh
{
    public:
        //Mesh();
        //~Mesh();

        Mesh<T> refine();

        void make1DMesh( int numberOfNodes, T start, T stop, bool symmetric = true );

        int getDimension(){return dimension;}
        int numberOfFreeNodes();

        void makeSquare();
        void makeLine(T start, T end);
        void makeQube();
        void brutforceMesh();

        void toGnuplot();
        void printFreeNodesCoordinates();

        void clearMemory();
        void clearNodes();
        void clearEdges();
        void clearElements();

        void renumerateNodes();//TODO - optimalization stuff


    //protected:
        std::list< Node<T>* > listOfNodes;               //TODO: clean memory
        std::list< Edge<T>* > listOfEdges;
        std::list< Element<T>* > listOfElements;
        
        int dimension;

    private:
            Node<T>* nodeExist( std::list<Node<T>* > &nodes, Node<T> *node );
            Edge<T>* edgeExist( std::list<Edge<T>* > &edges, Edge<T> *edge );

};

template <typename T>
int Mesh<T>::numberOfFreeNodes()
{
    int number = 0;
    for( typename std::list<Node<T>*>::iterator it = listOfNodes.begin(); it != listOfNodes.end(); it++ ) if( (*it)->isFree() ) number++;
    return number;
}


template <typename T>
Mesh<T> Mesh<T>::refine()
{
    Mesh<T> newMesh;

    newMesh.listOfNodes = this->listOfNodes;
    newMesh.dimension = dimension;

    std::list< Node<T>* > temporaryListOfNewNodes;
    std::list< Edge<T>* > temporaryListOfNewEdges;

    Edge<T>::resetCounter();
    
    for( typename std::list< Element<T>* >::iterator elementsIterator = listOfElements.begin(); elementsIterator != listOfElements.end(); elementsIterator++ )
    {
        std::list< Node<T>* > elementNewNodes;
        std::list< Edge<T>* > elementNewEdges;
        
        for( typename std::list< EdgeWithOrientation<T>*>::iterator edgeIterator = (*elementsIterator)->listOfEdges.begin(); edgeIterator != (*elementsIterator)->listOfEdges.end(); edgeIterator++ )
        {
            Node<T> *node = new Node<T>( (*edgeIterator)->edge->nodeStart->getDimension(), (*edgeIterator)->edge->isFree() );

            //setting up new node coordinates n1(x1,y1,...) n2(x2,y2,...) -> new node 1/2(x1+x2,y1+y2,...)
            for( int i = 0; i < (*edgeIterator)->edge->nodeStart->getDimension(); i++ )
                node->setIthCoordinate(i, ( (*edgeIterator)->edge->nodeStart->getCoordinate(i) + (*edgeIterator)->edge->nodeEnd->getCoordinate(i) )/2 );
            
            Node<T> *tmpNode;
            //check if the node already exist
            if( ( tmpNode = nodeExist(temporaryListOfNewNodes, node) ) )
            {
                node->reverseCounter();
                delete node;
                node = tmpNode;

                for( typename std::list<Edge<T>*>::iterator it = temporaryListOfNewEdges.begin(); it != temporaryListOfNewEdges.end(); it++ )
                {
                    if( ( (*it)->nodeStart == node ) || ( (*it)->nodeEnd == node ) ) elementNewEdges.push_back( *it );
                }
            }
            else
            {
                newMesh.listOfNodes.push_back(node);
                temporaryListOfNewNodes.push_back(node);
                
                //if node doesn't exist then edges don't exist -> create new corresponding edges
                Edge<T>* edge = new Edge<T>( (*edgeIterator)->edge->nodeStart, node, ( (*edgeIterator)->edge->nodeStart->isFree() || node->isFree() ) );
                temporaryListOfNewEdges.push_back(edge);
                newMesh.listOfEdges.push_back(edge);
                elementNewEdges.push_back(edge);
                
                edge = new Edge<T>( (*edgeIterator)->edge->nodeEnd, node, ( (*edgeIterator)->edge->nodeEnd->isFree() || node->isFree() ) );
                temporaryListOfNewEdges.push_back(edge);
                newMesh.listOfEdges.push_back(edge);
                elementNewEdges.push_back(edge);
            }
            elementNewNodes.push_back(node);
        }
        //create new Element
        //elementNewNodes + elementNewEdges = new Element
        //first new Edges and new Elment on new Nodes
        
        elementNewNodes.sort( compareNodeByNumber<T>() );

        //TODO: read next todo
        Element<T> *element = new Element<T>(dimension);
        for( typename std::list<Node<T>*>::iterator it1 = elementNewNodes.begin(); it1 != elementNewNodes.end(); it1++ )
        {
            for( typename std::list<Node<T>*>::iterator it2 = it1; it2 != elementNewNodes.end(); it2++)
            {
                if( it1 == it2 ) continue;
                Edge<T>* edge = new Edge<T>( *it1, *it2 );
                element->addEdge(edge);
                elementNewEdges.push_back(edge);
            }
        }
        //TODO: when dimension = 1 then new element will be empty so we don't want to push it to listOfElements
        if( dimension != 1 ) newMesh.listOfElements.push_back(element);
        else delete element;

        //now we have all edges in elementNewEdges

        while( !elementNewEdges.empty() )
        {
            element = new Element<T>(dimension);
            Edge<T> *tmpEdge = elementNewEdges.front();
            elementNewEdges.pop_front();
            element->addEdge( tmpEdge );
            //iterate through all edges to separate ones starting at the same Node
            for( typename std::list<Edge<T>*>::iterator it = elementNewEdges.begin(); it != elementNewEdges.end(); it++ )
            {
                if( tmpEdge->nodeStart->isSame( *((*it)->nodeStart )) )
                {
                    element->addEdge( *it );
                    elementNewEdges.erase(it++);
                }
            }
            //iterate through separated edges and edges to make element completed
            for( typename std::list<EdgeWithOrientation<T>*>::iterator it1 = element->listOfEdges.begin(); it1 != element->listOfEdges.end(); it1++ )
            {
                for( typename std::list<EdgeWithOrientation<T>*>::iterator it2 = element->listOfEdges.begin(); it2 != element->listOfEdges.end(); it2++ )
                {
                    if( it1 == it2 ) continue;
                    for( typename std::list<Edge<T>*>::reverse_iterator it3 = elementNewEdges.rbegin(); it3 != elementNewEdges.rend(); it3++ )
                    {
                        if(
                            ( (*it1)->edge->nodeEnd->isSame( *((*it3)->nodeStart )) && (*it2)->edge->nodeEnd->isSame( *((*it3)->nodeEnd )) ) ||
                            ( (*it1)->edge->nodeEnd->isSame( *((*it3)->nodeEnd )) && (*it2)->edge->nodeEnd->isSame( *((*it3)->nodeStart )) )
                          )
                        {
                            element->addEdge( *it3 );
                            elementNewEdges.erase( --it3.base() );
                        }
                    }
                }
            }
            newMesh.listOfElements.push_back(element);
        }
        //EOF create new Element
    }

    clearEdges();
    clearElements();
    //if( dimension ==  1 ) renumerateNodes();

    return newMesh;
}

template <typename T>
void Mesh<T>::make1DMesh( int numberOfNodes, T left, T right, bool symmetric )
{
    if( left >= right ) throw BoundaryException<T>( "Mesh.hpp make1DMesh - wrong bounds of mesh.", left, right );

    dimension = 1;

    T meshDivision = capd::abs( (right - left)/(numberOfNodes + 1) );

    Node<T> *node1 = new Node<T>(dimension,false);
    node1->setIthCoordinate( 0, left );
    listOfNodes.push_back( node1 );

    Node<T> *node2;
    Edge<T> *edge;
    Element<T> *element;

    int stepMultipler = 1;

    do{
        node2 = new Node<T>(dimension,true);
        node2->setIthCoordinate( 0, left + stepMultipler * meshDivision );
        listOfNodes.push_back( node2 );

        edge = new Edge<T>(node1,node2);
        listOfEdges.push_back( edge );

        element = new Element<T>(dimension);
        element->addEdge(edge,false);
        listOfElements.push_back(element);

        stepMultipler++;
        node1 = node2;
    }while( numberOfNodes-- );
    node2->setConstrained();
}

//check if the node already exist in nodes
template <typename T>
Node<T>* Mesh<T>::nodeExist(std::list<Node<T>*> &nodes, Node<T> *node)
{
    for( typename std::list<Node<T>*>::iterator nodesIterator = nodes.begin(); nodesIterator != nodes.end(); nodesIterator++ )
    {
        if( (*nodesIterator)->isSame(*node) ){return *nodesIterator;}
    }
    return NULL;
}

//check if the edge already exist in edges
template <typename T>
Edge<T>* Mesh<T>::edgeExist(std::list<Edge<T>*>& edges, Edge<T>* edge)
{
    for( typename std::list<Edge<T>*>::iterator edgesIterator = edges.begin(); edgesIterator != edges.end(); edgesIterator++ )
    {
        if( ( (*edgesIterator)->nodeStart->isSame(*(edge->nodeStart) ) ) &&
            ( (*edgesIterator)->nodeEnd->isSame(*(edge->nodeEnd) ) ) )
            return *edgesIterator;
    }
    return NULL;
}

//
template <typename T>
void Mesh<T>::renumerateNodes()
{

    if( dimension != 1 )
    {
        std::cout << "Mesh.hpp renumerateNodes() -> dimension != 1" << std::endl;
        return;// TODO: temporary code works only for dimension == 1
    }

    listOfNodes.sort( compareNodeByCoordinates<T>() );

    int nodeCount = 1;
    for( typename std::list<Node<T>*>::iterator nodesIterator = listOfNodes.begin(); nodesIterator != listOfNodes.end(); nodesIterator++ )
    {
        (*nodesIterator)->setNumber(nodeCount++);
    }
}

template <typename T>
void Mesh<T>::makeLine(T start, T end)
{
    dimension = 1;
    Node<T> *n1 = new Node<T>(dimension,false), *n2 = new Node<T>(dimension,false);

    n1->setIthCoordinate(0, start);
    n2->setIthCoordinate(0, end);

    listOfNodes.push_back(n1);
    listOfNodes.push_back(n2);

    Edge<T> *e1 = new Edge<T>(n1,n2);

    listOfEdges.push_back(e1);

    Element<T> *el1 = new Element<T>(dimension);
    el1->addEdge(e1,false);

    listOfElements.push_back(el1);
}

template <typename T>
void Mesh<T>::makeSquare()
{
    dimension = 2;
    Node<T> *n1 = new Node<T>(dimension,false), *n2 = new Node<T>(dimension,false), *n3 = new Node<T>(dimension,false), *n4 = new Node<T>(dimension,false);

    n1->setIthCoordinate(0, 0);
    n1->setIthCoordinate(1, 0);

    n2->setIthCoordinate(0, 1);
    n2->setIthCoordinate(1, 0);

    n3->setIthCoordinate(0, 0);
    n3->setIthCoordinate(1, 1);

    n4->setIthCoordinate(0, 1);
    n4->setIthCoordinate(1, 1);
    
    listOfNodes.push_back(n1);
    listOfNodes.push_back(n2);
    listOfNodes.push_back(n3);
    listOfNodes.push_back(n4);

    Edge<T> *e1 = new Edge<T>(n1,n2,false), *e2 = new Edge<T>(n1, n3, false), *e3 = new Edge<T>(n1,n4), *e4 = new Edge<T>(n2,n4,false), *e5 = new Edge<T>(n3,n4,false);

    listOfEdges.push_back(e1);
    listOfEdges.push_back(e2);
    listOfEdges.push_back(e3);
    listOfEdges.push_back(e4);
    listOfEdges.push_back(e5);

    Element<T> *el1 = new Element<T>(dimension);
    el1->addEdge(e2,false);
    el1->addEdge(e3);
    el1->addEdge(e5,false);

    Element<T> *el2 = new Element<T>(dimension);
    el2->addEdge(e1);
    el2->addEdge(e4);
    el2->addEdge(e3,false);

    listOfElements.push_back(el1);
    listOfElements.push_back(el2);
}

template <typename T>
void Mesh<T>::brutforceMesh()
{
    dimension = 2;

    Node<T> *n1 = new Node<T>(dimension,false), *n2 = new Node<T>(dimension,false), *n3 = new Node<T>(dimension,false),
            *n4 = new Node<T>(dimension,false), *n5 = new Node<T>(dimension,false), *n6 = new Node<T>(dimension,false),
            *n7 = new Node<T>(dimension), *n8 = new Node<T>(dimension), *n9 = new Node<T>(dimension),
            *n10 = new Node<T>(dimension,false), *n11 = new Node<T>(dimension,false), *n12 = new Node<T>(dimension),
            *n13 = new Node<T>(dimension), *n14 = new Node<T>(dimension), *n15 = new Node<T>(dimension,false),
            *n16 = new Node<T>(dimension,false), *n17 = new Node<T>(dimension), *n18 = new Node<T>(dimension),
            *n19 = new Node<T>(dimension), *n20 = new Node<T>(dimension,false), *n21 = new Node<T>(dimension,false),
            *n22 = new Node<T>(dimension,false), *n23 = new Node<T>(dimension,false), *n24 = new Node<T>(dimension,false),
            *n25 = new Node<T>(dimension,false);

    n1->setIthCoordinate(0, 0);
    n1->setIthCoordinate(1, 0);

    n2->setIthCoordinate(0, 0.25);
    n2->setIthCoordinate(1, 0);

    n3->setIthCoordinate(0, 0.5);
    n3->setIthCoordinate(1, 0);
    
    n4->setIthCoordinate(0, 0.75);
    n4->setIthCoordinate(1, 0);

    n5->setIthCoordinate(0, 1);
    n5->setIthCoordinate(1, 0);

    n6->setIthCoordinate(0, 0);
    n6->setIthCoordinate(1, 0.25);

    n7->setIthCoordinate(0, 0.25);
    n7->setIthCoordinate(1, 0.25);

    n8->setIthCoordinate(0, 0.5);
    n8->setIthCoordinate(1, 0.25);
    
    n9->setIthCoordinate(0, 0.75);
    n9->setIthCoordinate(1, 0.25);

    n10->setIthCoordinate(0, 1);
    n10->setIthCoordinate(1, 0.25);

    n11->setIthCoordinate(0, 0);
    n11->setIthCoordinate(1, 0.5);

    n12->setIthCoordinate(0, 0.25);
    n12->setIthCoordinate(1, 0.5);

    n13->setIthCoordinate(0, 0.5);
    n13->setIthCoordinate(1, 0.5);
    
    n14->setIthCoordinate(0, 0.75);
    n14->setIthCoordinate(1, 0.5);

    n15->setIthCoordinate(0, 1);
    n15->setIthCoordinate(1, 0.5);

    n16->setIthCoordinate(0, 0);
    n16->setIthCoordinate(1, 0.75);

    n17->setIthCoordinate(0, 0.25);
    n17->setIthCoordinate(1, 0.75);

    n18->setIthCoordinate(0, 0.5);
    n18->setIthCoordinate(1, 0.75);
    
    n19->setIthCoordinate(0, 0.75);
    n19->setIthCoordinate(1, 0.75);

    n20->setIthCoordinate(0, 1);
    n20->setIthCoordinate(1, 0.75);
    
    n21->setIthCoordinate(0, 0);
    n21->setIthCoordinate(1, 1);

    n22->setIthCoordinate(0, 0.25);
    n22->setIthCoordinate(1, 1);

    n23->setIthCoordinate(0, 0.5);
    n23->setIthCoordinate(1, 1);
    
    n24->setIthCoordinate(0, 0.75);
    n24->setIthCoordinate(1, 1);

    n25->setIthCoordinate(0, 1);
    n25->setIthCoordinate(1, 1);

    listOfNodes.push_back(n1);
    listOfNodes.push_back(n2);
    listOfNodes.push_back(n3);
    listOfNodes.push_back(n4);
    listOfNodes.push_back(n5);
    listOfNodes.push_back(n6);
    listOfNodes.push_back(n7);
    listOfNodes.push_back(n8);
    listOfNodes.push_back(n9);
    listOfNodes.push_back(n10);    
    listOfNodes.push_back(n11);
    listOfNodes.push_back(n12);
    listOfNodes.push_back(n13);
    listOfNodes.push_back(n14);
    listOfNodes.push_back(n15);
    listOfNodes.push_back(n16);
    listOfNodes.push_back(n17);
    listOfNodes.push_back(n18);
    listOfNodes.push_back(n19);
    listOfNodes.push_back(n20);
    listOfNodes.push_back(n21);
    listOfNodes.push_back(n22);
    listOfNodes.push_back(n23);
    listOfNodes.push_back(n24);
    listOfNodes.push_back(n25);

    Edge<T> *e1 = new Edge<T>(n1,n6,false), *e2 = new Edge<T>(n1,n7), *e3 = new Edge<T>(n1,n2,false), *e4 = new Edge<T>(n2,n7),
            *e5 = new Edge<T>(n2,n8), *e6 = new Edge<T>(n2,n3,false), *e7 = new Edge<T>(n3,n8), *e8 = new Edge<T>(n3,n9),
            *e9 = new Edge<T>(n3,n4,false), *e10 = new Edge<T>(n4,n9), *e11 = new Edge<T>(n4,n10), *e12 = new Edge<T>(n4,n5,false),
            *e13 = new Edge<T>(n5,n10,false), *e14 = new Edge<T>(n6,n11,false), *e15 = new Edge<T>(n6,n12), *e16 = new Edge<T>(n6,n7),
            *e17 = new Edge<T>(n7,n12), *e18 = new Edge<T>(n7,n13), *e19 = new Edge<T>(n7,n8), *e20 = new Edge<T>(n8,n13),
            *e21 = new Edge<T>(n8,n14), *e22 = new Edge<T>(n8,n9), *e23 = new Edge<T>(n9,n14), *e24 = new Edge<T>(n9,n15),
            *e25 = new Edge<T>(n9,n10), *e26 = new Edge<T>(n10,n15,false), *e27 = new Edge<T>(n11,n16,false), *e28 = new Edge<T>(n11,n17),
            *e29 = new Edge<T>(n11,n12), *e30 = new Edge<T>(n12,n17), *e31 = new Edge<T>(n12,n18), *e32 = new Edge<T>(n12,n13),
            *e33 = new Edge<T>(n13,n18), *e34 = new Edge<T>(n13,n19), *e35 = new Edge<T>(n13,n14), *e36 = new Edge<T>(n14,n19),
            *e37 = new Edge<T>(n14,n20), *e38 = new Edge<T>(n14,n15), *e39 = new Edge<T>(n15,n20,false), *e40 = new Edge<T>(n16,n21,false),
            *e41 = new Edge<T>(n16,n22), *e42 = new Edge<T>(n16,n17), *e43 = new Edge<T>(n17,n22), *e44 = new Edge<T>(n17,n23),
            *e45 = new Edge<T>(n17,n18), *e46 = new Edge<T>(n18,n23), *e47 = new Edge<T>(n18,n24), *e48 = new Edge<T>(n18,n19),
            *e49 = new Edge<T>(n19,n24), *e50 = new Edge<T>(n19,n25), *e51 = new Edge<T>(n19,n20), *e52 = new Edge<T>(n20,n25,false),
            *e53 = new Edge<T>(n21,n22,false), *e54 = new Edge<T>(n22,n23,false), *e55 = new Edge<T>(n23,n24,false), *e56 = new Edge<T>(n24,n25,false);

    listOfEdges.push_back(e1);
    listOfEdges.push_back(e2);
    listOfEdges.push_back(e3);
    listOfEdges.push_back(e4);
    listOfEdges.push_back(e5);
    listOfEdges.push_back(e6);
    listOfEdges.push_back(e7);
    listOfEdges.push_back(e8);
    listOfEdges.push_back(e9);
    listOfEdges.push_back(e10);
    listOfEdges.push_back(e11);
    listOfEdges.push_back(e12);
    listOfEdges.push_back(e13);
    listOfEdges.push_back(e14);
    listOfEdges.push_back(e15);
    listOfEdges.push_back(e16);
    listOfEdges.push_back(e17);
    listOfEdges.push_back(e18);
    listOfEdges.push_back(e19);
    listOfEdges.push_back(e20);
    listOfEdges.push_back(e21);
    listOfEdges.push_back(e22);
    listOfEdges.push_back(e23);
    listOfEdges.push_back(e24);
    listOfEdges.push_back(e25);
    listOfEdges.push_back(e26);
    listOfEdges.push_back(e27);
    listOfEdges.push_back(e28);
    listOfEdges.push_back(e29);
    listOfEdges.push_back(e30);
    listOfEdges.push_back(e41);
    listOfEdges.push_back(e42);
    listOfEdges.push_back(e43);
    listOfEdges.push_back(e44);
    listOfEdges.push_back(e45);
    listOfEdges.push_back(e46);
    listOfEdges.push_back(e47);
    listOfEdges.push_back(e48);
    listOfEdges.push_back(e49);
    listOfEdges.push_back(e50);
    listOfEdges.push_back(e51);
    listOfEdges.push_back(e52);
    listOfEdges.push_back(e53);
    listOfEdges.push_back(e54);
    listOfEdges.push_back(e55);
    listOfEdges.push_back(e56);

    Element<T> *el1 = new Element<T>(dimension);
    el1->addEdge(e1);
    el1->addEdge(e2);
    el1->addEdge(e16);
    
    Element<T> *el2 = new Element<T>(dimension);
    el2->addEdge(e2);
    el2->addEdge(e3);
    el2->addEdge(e4);

    Element<T> *el3 = new Element<T>(dimension);
    el3->addEdge(e4);
    el3->addEdge(e5);
    el3->addEdge(e19);

    Element<T> *el4 = new Element<T>(dimension);
    el4->addEdge(e5);
    el4->addEdge(e6);
    el4->addEdge(e7);

    Element<T> *el5 = new Element<T>(dimension);
    el5->addEdge(e7);
    el5->addEdge(e8);
    el5->addEdge(e22);

    Element<T> *el6 = new Element<T>(dimension);
    el6->addEdge(e8);
    el6->addEdge(e9);
    el6->addEdge(e10);

    Element<T> *el7 = new Element<T>(dimension);
    el7->addEdge(e10);
    el7->addEdge(e11);
    el7->addEdge(e25);

    Element<T> *el8 = new Element<T>(dimension);
    el8->addEdge(e11);
    el8->addEdge(e12);
    el8->addEdge(e13);

    Element<T> *el9 = new Element<T>(dimension);
    el9->addEdge(e14);
    el9->addEdge(e15);
    el9->addEdge(e29);

    Element<T> *el10 = new Element<T>(dimension);
    el10->addEdge(e15);
    el10->addEdge(e16);
    el10->addEdge(e17);

    Element<T> *el11 = new Element<T>(dimension);
    el11->addEdge(e17);
    el11->addEdge(e18);
    el11->addEdge(e32);
    
    Element<T> *el12 = new Element<T>(dimension);
    el12->addEdge(e18);
    el12->addEdge(e19);
    el12->addEdge(e20);

    Element<T> *el13 = new Element<T>(dimension);
    el13->addEdge(e20);
    el13->addEdge(e21);
    el13->addEdge(e35);

    Element<T> *el14 = new Element<T>(dimension);
    el14->addEdge(e21);
    el14->addEdge(e22);
    el14->addEdge(e23);

    Element<T> *el15 = new Element<T>(dimension);
    el15->addEdge(e23);
    el15->addEdge(e24);
    el15->addEdge(e38);

    Element<T> *el16 = new Element<T>(dimension);
    el16->addEdge(e24);
    el16->addEdge(e25);
    el16->addEdge(e26);

    Element<T> *el17 = new Element<T>(dimension);
    el17->addEdge(e27);
    el17->addEdge(e28);
    el17->addEdge(e42);

    Element<T> *el18 = new Element<T>(dimension);
    el18->addEdge(e28);
    el18->addEdge(e29);
    el18->addEdge(e30);

    Element<T> *el19 = new Element<T>(dimension);
    el19->addEdge(e30);
    el19->addEdge(e31);
    el19->addEdge(e45);

    Element<T> *el20 = new Element<T>(dimension);
    el20->addEdge(e31);
    el20->addEdge(e32);
    el20->addEdge(e33);

    Element<T> *el21 = new Element<T>(dimension);
    el21->addEdge(e33);
    el21->addEdge(e34);
    el21->addEdge(e48);
    
    Element<T> *el22 = new Element<T>(dimension);
    el22->addEdge(e34);
    el22->addEdge(e35);
    el22->addEdge(e36);

    Element<T> *el23 = new Element<T>(dimension);
    el23->addEdge(e36);
    el23->addEdge(e37);
    el23->addEdge(e51);

    Element<T> *el24 = new Element<T>(dimension);
    el24->addEdge(e37);
    el24->addEdge(e38);
    el24->addEdge(e39);

    Element<T> *el25 = new Element<T>(dimension);
    el25->addEdge(e40);
    el25->addEdge(e41);
    el25->addEdge(e53);

    Element<T> *el26 = new Element<T>(dimension);
    el26->addEdge(e41);
    el26->addEdge(e42);
    el26->addEdge(e43);

    Element<T> *el27 = new Element<T>(dimension);
    el27->addEdge(e43);
    el27->addEdge(e44);
    el27->addEdge(e54);

    Element<T> *el28 = new Element<T>(dimension);
    el28->addEdge(e44);
    el28->addEdge(e45);
    el28->addEdge(e46);

    Element<T> *el29 = new Element<T>(dimension);
    el29->addEdge(e46);
    el29->addEdge(e47);
    el29->addEdge(e55);

    Element<T> *el30 = new Element<T>(dimension);
    el30->addEdge(e47);
    el30->addEdge(e48);
    el30->addEdge(e49);

    Element<T> *el31 = new Element<T>(dimension);
    el31->addEdge(e49);
    el31->addEdge(e50);
    el31->addEdge(e56);

    Element<T> *el32 = new Element<T>(dimension);
    el32->addEdge(e50);
    el32->addEdge(e51);
    el32->addEdge(e52);
    
    listOfElements.push_back(el1);
    listOfElements.push_back(el2);
    listOfElements.push_back(el3);
    listOfElements.push_back(el4);
    listOfElements.push_back(el5);
    listOfElements.push_back(el6);
    listOfElements.push_back(el7);
    listOfElements.push_back(el8);
    listOfElements.push_back(el9);
    listOfElements.push_back(el10);
    listOfElements.push_back(el11);
    listOfElements.push_back(el12);
    listOfElements.push_back(el13);
    listOfElements.push_back(el14);
    listOfElements.push_back(el15);
    listOfElements.push_back(el16);
    listOfElements.push_back(el17);
    listOfElements.push_back(el18);
    listOfElements.push_back(el19);
    listOfElements.push_back(el20);
    listOfElements.push_back(el21);
    listOfElements.push_back(el22);
    listOfElements.push_back(el23);
    listOfElements.push_back(el24);
    listOfElements.push_back(el25);
    listOfElements.push_back(el26);
    listOfElements.push_back(el27);
    listOfElements.push_back(el28);
    listOfElements.push_back(el29);
    listOfElements.push_back(el30);
    listOfElements.push_back(el31);
    listOfElements.push_back(el32);
}

template <typename T>
void Mesh<T>::makeQube()
{
    dimension = 3;
    Node<T> *n1 = new Node<T>(dimension), *n2 = new Node<T>(dimension), *n3 = new Node<T>(dimension),
            *n4 = new Node<T>(dimension), *n5 = new Node<T>(dimension), *n6 = new Node<T>(dimension),
            *n7 = new Node<T>(dimension), *n8 = new Node<T>(dimension);

    n1->setIthCoordinate(0, 0);
    n1->setIthCoordinate(1, 0);
    n1->setIthCoordinate(2, 0);

    n2->setIthCoordinate(0, 1);
    n2->setIthCoordinate(1, 0);
    n2->setIthCoordinate(2, 0);

    n3->setIthCoordinate(0, 1);
    n3->setIthCoordinate(1, 0);
    n3->setIthCoordinate(2, 1);

    n4->setIthCoordinate(0, 0);
    n4->setIthCoordinate(1, 0);
    n4->setIthCoordinate(2, 1);
    
    n5->setIthCoordinate(0, 0);
    n5->setIthCoordinate(1, 1);
    n5->setIthCoordinate(2, 0);

    n6->setIthCoordinate(0, 1);
    n6->setIthCoordinate(1, 1);
    n6->setIthCoordinate(2, 0);

    n7->setIthCoordinate(0, 1);
    n7->setIthCoordinate(1, 1);
    n7->setIthCoordinate(2, 1);

    n8->setIthCoordinate(0, 0);
    n8->setIthCoordinate(1, 1);
    n8->setIthCoordinate(2, 1);    
    
    listOfNodes.push_back(n1);
    listOfNodes.push_back(n2);
    listOfNodes.push_back(n3);
    listOfNodes.push_back(n4);
    listOfNodes.push_back(n5);
    listOfNodes.push_back(n6);
    listOfNodes.push_back(n7);
    listOfNodes.push_back(n8);

    Edge<T> *e1 = new Edge<T>(n1,n2), *e2 = new Edge<T>(n2,n3), *e3 = new Edge<T>(n3,n4), *e4 = new Edge<T>(n1,n4),
            *e5 = new Edge<T>(n5,n6), *e6 = new Edge<T>(n6,n7), *e7 = new Edge<T>(n7,n8), *e8 = new Edge<T>(n5, n8),
            *e9 = new Edge<T>(n1,n5), *e10 = new Edge<T>(n2,n6), *e11 = new Edge<T>(n3, n7), *e12 = new Edge<T>(n4,n8),
            *e13 = new Edge<T>(n2,n5), *e14 = new Edge<T>(n2,n7), *e15 = new Edge<T>(n4,n7), *e16 = new Edge<T>(n4,n5),
            *e17 = new Edge<T>(n5,n7), *e18 = new Edge<T>(n2,n4);

    listOfEdges.push_back(e1);
    listOfEdges.push_back(e2);
    listOfEdges.push_back(e3);
    listOfEdges.push_back(e4);
    listOfEdges.push_back(e5);
    listOfEdges.push_back(e6);
    listOfEdges.push_back(e7);
    listOfEdges.push_back(e8);
    listOfEdges.push_back(e9);
    listOfEdges.push_back(e10);
    listOfEdges.push_back(e11);
    listOfEdges.push_back(e12);
    listOfEdges.push_back(e13);
    listOfEdges.push_back(e14);
    listOfEdges.push_back(e15);
    listOfEdges.push_back(e16);
    listOfEdges.push_back(e17);
    listOfEdges.push_back(e18);

    Element<T> *el1 = new Element<T>(dimension);
    el1->addEdge(e1);
    el1->addEdge(e4);
    el1->addEdge(e9);
    el1->addEdge(e13);
    el1->addEdge(e16);
    el1->addEdge(e18);

    Element<T> *el2 = new Element<T>(dimension);
    el2->addEdge(e2);
    el2->addEdge(e3);
    el2->addEdge(e18);
    el2->addEdge(e11);
    el2->addEdge(e14);
    el2->addEdge(e15);

    Element<T> *el3 = new Element<T>(dimension);
    el3->addEdge(e14);
    el3->addEdge(e13);
    el3->addEdge(e17);
    el3->addEdge(e10);
    el3->addEdge(e5);
    el3->addEdge(e6);
    
    Element<T> *el4 = new Element<T>(dimension);
    el4->addEdge(e16);
    el4->addEdge(e15);
    el4->addEdge(e17);
    el4->addEdge(e12);
    el4->addEdge(e7);
    el4->addEdge(e8);

    listOfElements.push_back(el1);
    listOfElements.push_back(el2);
    listOfElements.push_back(el3);
    listOfElements.push_back(el4);
}


template <typename T>
void Mesh<T>::toGnuplot()
{
    std::string fileName = "constrainedN.txt";
    std::string fileName1 = "freeN.txt";
    std::string fileName2 = "plot.txt";

    std::string fileName3 = "constrainedE.txt";
    std::string fileName4 = "freeE.txt";
    
    std::ofstream outFile( fileName.c_str() );
    std::ofstream outFile1( fileName1.c_str() );
    std::ofstream outFile2( fileName2.c_str() );
    
    
    for( typename std::list<Node<T>*>::iterator it = listOfNodes.begin(); it != listOfNodes.end(); it++ )
    {
        for( int i = 0; i < (*it)->getDimension(); i++ )
        {
            if( (*it)->isFree() ) outFile1 << (*it)->getCoordinate(i).mid().leftBound() << "\t";
            else outFile << (*it)->getCoordinate(i).mid().leftBound() << "\t";
        }
        if( (*it)->isFree() ) outFile1 << "\n";
        else outFile << "\n";
    }
    
    if(dimension < 3 )
    outFile2 << "plot '" + fileName1 + "' with points pt 12, '" + fileName + "' with points pt 12\npause -1";
    else
    outFile2 << "splot '" + fileName1 + "' with points pt 12, '" + fileName + "' with points pt 12\npause -1";
    
    outFile.close();
    outFile1.close();
    outFile2.close();
}

template <typename T>
void Mesh<T>::printFreeNodesCoordinates()
{
    std::cout << "Mesh: ";
    for( typename std::list<Node<T>*>::iterator nodeIterator = listOfNodes.begin(); nodeIterator != listOfNodes.end(); nodeIterator++ )
    {
        if( (*nodeIterator)->isFree() )
        {
            if( dimension == 1 ) std::cout << (*nodeIterator)->getCoordinate(0) << ",";
            if( dimension == 2 ) std::cout << "p" << (*nodeIterator)->getCoordinate(0) << "," << (*nodeIterator)->getCoordinate(1) << ",";
            if( dimension > 2 ) std::cout << "Mesh.hpp - print Free Nodes Coordinates -> dimension > 2" << std::endl;
        }
    }
    std::cout << std::endl;
}

template <typename T>
void Mesh<T>::clearMemory()
{
    clearEdges();
    clearElements();
    clearNodes();
}

template <typename T>
void Mesh<T>::clearNodes()
{
    for( typename std::list<Node<T>*>::iterator it = listOfNodes.begin(); it != listOfNodes.end(); it++ ) delete (*it);
}

template <typename T>
void Mesh<T>::clearEdges()
{
    for( typename std::list<Edge<T>*>::iterator it = listOfEdges.begin(); it != listOfEdges.end(); it++ ) delete (*it);
}

template <typename T>
void Mesh<T>::clearElements()
{
    for( typename std::list<Element<T>*>::iterator it = listOfElements.begin(); it != listOfElements.end(); it++ ) delete (*it);
}

#endif // _MESH_HPP_
