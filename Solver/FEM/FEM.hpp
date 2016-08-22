#ifndef _FEM_HPP_
#define _FEM_HPP_

#include <algorithm>
#include <fstream>
#include <cmath>

#include "capd/capdlib.h"
#include "capd/vectalg/Matrix.hpp"
#include "capd/vectalg/Vector.hpp"

#include "../Mesh/Mesh.hpp"
#include "../Equation/Equation.hpp"
#include "../EigenTools/EigenTools.h"

//this class has forward-declaration in few files, so watch it when template argruments are altered
template <typename T>
class FEM
{
    public:
        FEM();
        ~FEM();

        void assembleStiffnessMatrix();
        void assembleMassMatrix();
        void assembleMassMatrix1D();
        void assembleMassMatrix2D();
        
        void invertMassMatrix();
        void invertsMassMatrixTimesStiffnessMatrix();

        void diagonalize();
        void diagonalize_simple();
        void diagonalize_for_candidate();

        void putGradientsIntoMatrix( capd::vectalg::Matrix<T,0,0> *matrixG, capd::vectalg::Matrix<T,0,0> *inverseOfMatrixM );
        //T dotProductOfVectors(int column1, int column2, capd::vectalg::Matrix<T,0,0> *inverseOfMatrixM);

        void solveLinearSystem();
 
        void computeCutCoordinateVector();


        void computeLoadVector(double time = 0 );
        void LoadVectorForInhomogeneousDirichletConditions( double time = 0 );
        void LoadVectorForHomogenousBoundaryCondition( double time = 0 );
        
        T timeDerivativeFactorOfLoadVector( Node<T> *centroidOfSimplex, T &areaOfSimplex, T timeDerivativeOfFunctionG );

        void EigenVectorsValues();

        bool courantCondition( double deltaT );
        void CAPD_EigenvectorsEigenvalues( capd::vectalg::Matrix<T,0,0> *matrix );

        void printAllData();


        Mesh<T> *mesh;
        Equation<T> *equation;

        int dimension;
        int size;

        capd::vectalg::Matrix<T,0,0> *P, *inverseOfP, *diagonalMatrix, *remainsFromDiagonal,*matrixForNonlinearPart;
        capd::vectalg::Vector<T,0> *diagonalVector;// vector made from diagonalMatrix so computation will be faster vector times vector instead matrix times vector

        //Au' + Ku + nlv = F
        capd::vectalg::Matrix<T,0,0> *massMatrix;//A
        capd::vectalg::Matrix<T,0,0> *stiffnessMatrix;//K
        capd::vectalg::Matrix<T,0,0> *inverseOfMassMatrix;//A^-1
        capd::vectalg::Matrix<T,0,0> *inverseOfMassMatrixTimesStiffnessMatrix;//A^{-1} * K

        capd::vectalg::Vector<T,0> *loadVector;//F
        capd::vectalg::Vector<T,0> *solutionVector;//U
        capd::vectalg::Vector<T,0> *nonLinearVector;//nlv

        capd::vectalg::Vector<int,0> *cutCoordinateVector;//helps to cut constrained nodes i.e. from stiffnessBeforeCutMatrix, massBeforeCutMatrix, etc...
};

template <typename T>
FEM<T>::FEM()
{
    stiffnessMatrix = NULL;
    loadVector = NULL;
    solutionVector = NULL;
    massMatrix = NULL;
    inverseOfMassMatrix = NULL;
    mesh = NULL;
    nonLinearVector = NULL;
    cutCoordinateVector = NULL;
    inverseOfMassMatrixTimesStiffnessMatrix = NULL;
    P = NULL; 
    inverseOfP = NULL;
    diagonalMatrix = NULL;
    matrixForNonlinearPart = NULL;
    diagonalVector = NULL;
    remainsFromDiagonal = NULL;
}

template <typename T>
FEM<T>::~FEM()
{
    delete stiffnessMatrix;
    delete loadVector;
    delete solutionVector;
    delete massMatrix;
    delete inverseOfMassMatrix;
    delete nonLinearVector;
    delete cutCoordinateVector;
    delete inverseOfMassMatrixTimesStiffnessMatrix;

    delete P;
    delete inverseOfP;
    delete diagonalMatrix;
    delete matrixForNonlinearPart;
    delete diagonalVector;
    delete remainsFromDiagonal;
    
    stiffnessMatrix = NULL;
    loadVector = solutionVector = NULL;
    massMatrix = inverseOfMassMatrix = NULL;
    mesh = NULL;
    nonLinearVector = NULL;
    cutCoordinateVector = NULL;
    inverseOfMassMatrixTimesStiffnessMatrix = NULL;
    P = inverseOfP = diagonalMatrix = NULL;
    diagonalVector = NULL;
    remainsFromDiagonal = NULL;
}

template <typename T>
void FEM<T>::computeCutCoordinateVector()
{
    std::vector<Node<T>*> tmpNodeVector;
    for( typename std::list<Node<T>*>::iterator nodeIterator = mesh->listOfNodes.begin(); nodeIterator != mesh->listOfNodes.end(); nodeIterator++ )
    {
        tmpNodeVector.push_back( (*nodeIterator) );
    }

    std::sort( tmpNodeVector.begin(), tmpNodeVector.end(), compareNodeByNumber<T>() );
    
    cutCoordinateVector = new capd::vectalg::Vector<int,0>( mesh->listOfNodes.size()+1 );//Nodes are numbered from 1 not 0
    
    int finalId;
    finalId = 0;

    for(int id = 0; id < (int)mesh->listOfNodes.size(); id++ )
    {
        if( tmpNodeVector[id]->isFree() )
        {
            (*cutCoordinateVector)[ tmpNodeVector[id]->getNumber() ] = finalId;
            finalId++;
        }
        else
        {
            (*cutCoordinateVector)[ tmpNodeVector[id]->getNumber() ] = -1;
        }
    }
}

/*
template <typename T>
T FEM<T>::dotProductOfVectors(int column1, int column2, capd::vectalg::Matrix<T,0,0> *inverseOfMatrixM)
{
    T tmp = 0;
    for( int row = 1; row < dimension + 1; row++ )
    {
        tmp += (*inverseOfMatrixM)[row][column1] * (*inverseOfMatrixM)[row][column2];
    }
    return tmp;
}
*/
template <typename T>
void FEM<T>::putGradientsIntoMatrix( capd::vectalg::Matrix<T,0,0> *matrixG, capd::vectalg::Matrix<T,0,0> *inverseOfMatrixM )
{
    for( int r = 0; r < dimension + 1; r++ )
    {
        for( int s = r; s < dimension + 1; s++ )
        {
            //(*matrixG)[r][s] = dotProductOfVectors( r, s, inverseOfMatrixM );
            for( int row = 1; row < dimension + 1; row++ )
            {
                (*matrixG)[r][s] += (*inverseOfMatrixM)[row][r] * (*inverseOfMatrixM)[row][s];
            }
        }
    }
}

template <typename T>
void FEM<T>::assembleStiffnessMatrix()
{
    stiffnessMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    for( typename std::list< Element<T>* >::iterator elementsIterator = mesh->listOfElements.begin(); elementsIterator != mesh->listOfElements.end(); elementsIterator++ )
    {
        std::list<Node<T>*> elementNodes;
        (*elementsIterator)->getElementNodes( elementNodes );

        int matrixRowIndex = 0;
        for( typename std::list<Node<T>*>::iterator it = elementNodes.begin(); it != elementNodes.end(); it++ )
        {
            (*(*elementsIterator)->matrixM)[matrixRowIndex][0] = 1;
            for( int i = 0; i < (*it)->getDimension(); i++ )
            {
                (*(*elementsIterator)->matrixM)[matrixRowIndex][i+1] = (*it)->getCoordinate(i);
            }
            ++matrixRowIndex;
        }
        (*(*elementsIterator)->inverseOfMatrixM) = capd::matrixAlgorithms::gaussInverseMatrix( (*(*elementsIterator)->matrixM) );

        capd::vectalg::Matrix<T,0,0> *matrixG = new capd::vectalg::Matrix<T,0,0>( elementNodes.size(), elementNodes.size() );

        putGradientsIntoMatrix( matrixG, (*elementsIterator)->inverseOfMatrixM );

        //step 1 page 163
        T areaOfSimplex = (*elementsIterator)->areaOfSimplex( elementNodes );

        //dot product of gradients
        int r,s;
        r = s = 0;
        for( typename std::list<Node<T>*>::iterator it1 = elementNodes.begin(); it1 != elementNodes.end(); it1++, r++ )
        {
            s = r;
            if( (*it1)->isFree() )
            {
                for( typename std::list<Node<T>*>::iterator it2 = it1; it2 != elementNodes.end(); it2++, s++ )
                {
                    if( (*it2)->isFree() )
                    {
                        (*stiffnessMatrix)[ (*cutCoordinateVector)[std::min( (*it1)->getNumber(), (*it2)->getNumber() )] ][ (*cutCoordinateVector)[std::max( (*it1)->getNumber(), (*it2)->getNumber() )] ] += (*matrixG)[r][s] * areaOfSimplex * equation->ni;
                        if( std::min( (*it1)->getNumber(), (*it2)->getNumber() ) !=  std::max( (*it1)->getNumber(), (*it2)->getNumber() ) )
                        {//thanks to this 'if' we will not sum twice same coefficient
                            (*stiffnessMatrix)[ (*cutCoordinateVector)[std::max( (*it1)->getNumber(), (*it2)->getNumber() )] ][ (*cutCoordinateVector)[std::min( (*it1)->getNumber(), (*it2)->getNumber() )] ] += (*matrixG)[r][s] * areaOfSimplex * equation->ni;
                        }
                    }
                }
            }

        }

        delete matrixG;
    }
}

template <typename T>
void FEM<T>::computeLoadVector( double time )
{
    LoadVectorForInhomogeneousDirichletConditions( time );
}

template <typename T>
void FEM<T>::LoadVectorForHomogenousBoundaryCondition( double time )
{
    if( loadVector != NULL ) delete loadVector;
    loadVector = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );
    capd::vectalg::Vector<T,0> *loadVectorBeforeCut = new capd::vectalg::Vector<T,0>( mesh->listOfNodes.size()+1 );//Nodes are numbered from 1 not 0

    T tmp;

    for( typename std::list< Element<T>* >::iterator elementsIterator = mesh->listOfElements.begin(); elementsIterator != mesh->listOfElements.end(); elementsIterator++ )
    {
        std::list<Node<T>*> elementNodes;
        (*elementsIterator)->getElementNodes( elementNodes );

        tmp = ( (*elementsIterator)->areaOfSimplex( elementNodes ) * equation->functionF( (*elementsIterator)->centroidOfSimplex( elementNodes ), time ) ) / ( dimension + 1 );

        for( typename std::list<Node<T>*>::iterator it = elementNodes.begin(); it != elementNodes.end(); it++ )
        {
            (*loadVectorBeforeCut)[(*it)->getNumber()] += tmp;
        }
    }

    int i = 0;
    for( typename std::list<Node<T>*>::iterator it = mesh->listOfNodes.begin(); it != mesh->listOfNodes.end(); it++ )
    {
        if( (*it)->isFree() ) (*loadVector)[i++] = (*loadVectorBeforeCut)[ (*it)->getNumber() ];
    }
    delete loadVectorBeforeCut;
}

template <typename T>
T FEM<T>::timeDerivativeFactorOfLoadVector( Node<T> *centroidOfSimplex, T &areaOfSimplex, T timeDerivativeOfFunctionG )
//centroid of simplex will be used in higher dimensional case (dim > 1)
{
    if( dimension != 1) return 0;
    return 1/4 * areaOfSimplex * timeDerivativeOfFunctionG;
}

template <typename T>
void FEM<T>::LoadVectorForInhomogeneousDirichletConditions( double time )
{
    if( loadVector != NULL ) delete loadVector;
    loadVector = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );
    capd::vectalg::Vector<T,0> *loadVectorBeforeCut = new capd::vectalg::Vector<T,0>( mesh->listOfNodes.size()+1 );//Nodes are numbered from 1 not 0

    for( typename std::list< Element<T>* >::iterator elementsIterator = mesh->listOfElements.begin(); elementsIterator != mesh->listOfElements.end(); elementsIterator++ )
    {
        std::list<Node<T>*> elementNodes;
        (*elementsIterator)->getElementNodes( elementNodes );
        Node<T> *centroidOfSimplex = (*elementsIterator)->centroidOfSimplex( elementNodes );

        //if f is nonzero
        if( true )
        {
            T tmp = ( (*elementsIterator)->areaOfSimplex( elementNodes ) * equation->functionF( centroidOfSimplex, time ) ) / ( dimension + 1 );
            for( typename std::list<Node<T>*>::iterator it = elementNodes.begin(); it != elementNodes.end(); it++ )
            {
                (*loadVectorBeforeCut)[(*it)->getNumber()] += tmp;
            }
        }

        //if g is nonzero and at least one node is constrained
	//attention, element should not have all nodes constrained
        bool condition = false;
        for( typename std::list<Node<T>*>::iterator it = elementNodes.begin(); it != elementNodes.end(); it++ )
        {
            if( (*it)->isConstrained() )
            {
                condition = true;
            }
        }
        if( condition )
        {
            //compute the gradient G on simplex p.170
            capd::vectalg::Vector<T,0> *gradientOfG = equation->gradientOfFunctionG( (*elementsIterator), elementNodes,dimension ,time );

            T areaOfSimplex = (*elementsIterator)->areaOfSimplex( elementNodes );
            int column = 0;

            for( typename std::list<Node<T>*>::iterator it = elementNodes.begin(); it != elementNodes.end(); it++, column++ )
            {
                if( (*it)->isFree() )
                {
                    T gradientDotProduct = 0;
                    for( int i = 0; i < dimension; i++ )
                    {
                        gradientDotProduct += (*gradientOfG)[i] * (*(*elementsIterator)->inverseOfMatrixM)[i+1][column];
                    }

                    (*loadVectorBeforeCut)[(*it)->getNumber()] -= gradientDotProduct * areaOfSimplex * equation->ni - timeDerivativeFactorOfLoadVector( centroidOfSimplex, areaOfSimplex, equation->timeDerivativeOfFunctionG( (*it), time ) );
                }
            }

            delete gradientOfG;
        }
        delete centroidOfSimplex;
    }

    int i = 0;
    for( typename std::list<Node<T>*>::iterator it = mesh->listOfNodes.begin(); it != mesh->listOfNodes.end(); it++ )
    {
        if( (*it)->isFree() ) (*loadVector)[i++] = (*loadVectorBeforeCut)[ (*it)->getNumber() ];
    }
    
    delete loadVectorBeforeCut;
}

template <typename T>
void FEM<T>::solveLinearSystem()
{
    solutionVector = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );
    (*solutionVector) = capd::matrixAlgorithms::gauss( *stiffnessMatrix, *loadVector);
    std::cout << "solutionVector \n" << *solutionVector <<std::endl;
}


template <typename T>
void FEM<T>::assembleMassMatrix()
{
    if( dimension == 1 ) assembleMassMatrix1D();
    if( dimension == 2 ) assembleMassMatrix2D();
    if( dimension > 2 ) std::cout << "FEM.hpp - Assemble Mass Matrix -> dimension > 2" << std::endl;
}

template <typename T>
void FEM<T>::assembleMassMatrix1D()
{
    massMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    T h = 0;
    for( typename std::list< Element<T>* >::iterator elementsIterator = mesh->listOfElements.begin(); elementsIterator != mesh->listOfElements.end(); elementsIterator++ )
    {
        std::list<Node<T>*> elementNodes;
        (*elementsIterator)->getElementNodes( elementNodes );

        if( elementNodes.front()->isFree() && elementNodes.back()->isFree() )
        {

            h = capd::abs( elementNodes.front()->getCoordinate(0) - elementNodes.back()->getCoordinate(0) );
            (*massMatrix)[ (*cutCoordinateVector)[ elementNodes.front()->getNumber() ] ][ (*cutCoordinateVector)[elementNodes.back()->getNumber()] ] = h / 6;
            (*massMatrix)[ (*cutCoordinateVector)[ elementNodes.back()->getNumber() ] ][ (*cutCoordinateVector)[elementNodes.front()->getNumber()] ] = h / 6;

            (*massMatrix)[ (*cutCoordinateVector)[ elementNodes.front()->getNumber() ] ][ (*cutCoordinateVector)[elementNodes.front()->getNumber()] ] += h / 3;
            (*massMatrix)[ (*cutCoordinateVector)[ elementNodes.back()->getNumber() ] ][ (*cutCoordinateVector)[elementNodes.back()->getNumber()] ] += h / 3;
        }
        else
        {
            h = capd::abs( elementNodes.front()->getCoordinate(0) - elementNodes.back()->getCoordinate(0) );
            if( elementNodes.front()->isFree() )
            {
                (*massMatrix)[ (*cutCoordinateVector)[ elementNodes.front()->getNumber() ] ][ (*cutCoordinateVector)[elementNodes.front()->getNumber()] ] += h / 3;
            }
            else
            {
                (*massMatrix)[ (*cutCoordinateVector)[ elementNodes.back()->getNumber() ] ][ (*cutCoordinateVector)[elementNodes.back()->getNumber()] ] += h / 3;
            }
        }
    }
}

template <typename T>
void FEM<T>::assembleMassMatrix2D()
{
    massMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );


    for( typename std::list< Element<T>* >::iterator elementsIterator = mesh->listOfElements.begin(); elementsIterator != mesh->listOfElements.end(); elementsIterator++ )
    {
        std::list<Node<T>*> elementNodes;
        (*elementsIterator)->getElementNodes( elementNodes );
        
        T areaOfSimplex = (*elementsIterator)->areaOfSimplex( elementNodes );
        
        for( typename std::list<Node<T>*>::iterator itI = elementNodes.begin(); itI != elementNodes.end(); itI++ )
        {
            for( typename std::list<Node<T>*>::iterator itJ = elementNodes.begin(); itJ != elementNodes.end(); itJ++ )
            {
                if( (*itI)->isFree() && (*itJ)->isFree() )
                {
                    int multiplier = itI == itJ ? 2 : 1;//3.66/3.67
                    (*massMatrix)[ (*cutCoordinateVector)[ (*itI)->getNumber() ] ][ (*cutCoordinateVector)[ (*itJ)->getNumber()] ] += ( multiplier * areaOfSimplex ) / 12;
                }
            }
        }
    }
}


template <typename T>
void FEM<T>::invertMassMatrix()
{
    inverseOfMassMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    (*inverseOfMassMatrix) = capd::matrixAlgorithms::gaussInverseMatrix( *massMatrix );
}

template <typename T>
void FEM<T>::invertsMassMatrixTimesStiffnessMatrix()
{
    inverseOfMassMatrixTimesStiffnessMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    
    if( inverseOfMassMatrix == NULL ) invertMassMatrix();
    
    (*inverseOfMassMatrixTimesStiffnessMatrix) = (*inverseOfMassMatrix) * (*stiffnessMatrix);
}


template <typename T>
void FEM<T>::diagonalize()
//TODO: I'm not proud of the following code - it works but is ugly as hell - needs to be rewritten in a spare time
{
    capd::vectalg::Vector<double,0> *eigenValueRealPart = new capd::vectalg::Vector<double,0>( mesh->numberOfFreeNodes() );
    capd::vectalg::Vector<double,0> *eigenValueImPart = new capd::vectalg::Vector<double,0>( mesh->numberOfFreeNodes() );

    capd::vectalg::Matrix<double,0,0> *eigenVectorRealPart = new capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<double,0,0> *eigenVectorImPart = new capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );


    capd::vectalg::Matrix<double,0,0> tmp = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            //make double matrix from interval matrix
            tmp[i][j] = (*inverseOfMassMatrixTimesStiffnessMatrix)[i][j].mid().leftBound();
        }

    capd::jaco::EigenTools<double,0> et;

    capd::vectalg::Matrix<double,0,0> schurT = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<double,0,0> schurS = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<double,0,0> inverseOfSchurS = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    //this procedure will fill matrices T(upper triangle) and S,  T = S^{-1}MS, where M = tmp
    et.schurDecompose( tmp, schurT, schurS );
    inverseOfSchurS = capd::matrixAlgorithms::gaussInverseMatrix( schurS );
    
    capd::vectalg::Matrix<T,0,0> intervalInverseOfSchurS = capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            intervalInverseOfSchurS[i][j] = (T)inverseOfSchurS[i][j] + T(-1,1);
        }
    
    et.krawczykRefineMatrix( schurS, inverseOfSchurS, intervalInverseOfSchurS );

    capd::alglib::computeEigenvaluesAndEigenvectors( schurT, (*eigenValueRealPart), (*eigenValueImPart), (*eigenVectorRealPart), (*eigenVectorImPart) );

    capd::vectalg::Matrix<double,0,0> E = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            E[i][j] = (*eigenVectorRealPart)[i][j];
        }


    capd::vectalg::Matrix<double,0,0> inverseOfE = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    inverseOfE = capd::matrixAlgorithms::gaussInverseMatrix( E );

    capd::vectalg::Matrix<T,0,0> intervalInverseOfE = capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            intervalInverseOfE[i][j] = (T)inverseOfE[i][j] + T(-1,1);
        }

    et.krawczykRefineMatrix( E, inverseOfE, intervalInverseOfE );

    diagonalMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    P = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    inverseOfP = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    (*inverseOfP) = intervalInverseOfE * intervalInverseOfSchurS;


    capd::vectalg::Matrix<double,0,0> doubleInverseOfP = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            doubleInverseOfP[i][j] = (*inverseOfP)[i][j].mid().leftBound();
        }

    tmp = capd::matrixAlgorithms::gaussInverseMatrix( doubleInverseOfP );

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            (*P)[i][j] = (T)tmp[i][j] + T(-1,1);
        }

    et.krawczykRefineMatrix( doubleInverseOfP, tmp, (*P) );

    (*diagonalMatrix) = (*inverseOfP) * (*inverseOfMassMatrixTimesStiffnessMatrix) * (*P);


    diagonalVector = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );

    std::vector<T> in;
    std::vector<T> out;

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        in.push_back( (*diagonalMatrix)[i][i] );

    out = in;
    //std::sort( out.begin(), out.end(), std::greater<T>() );//descending
    std::sort( out.begin(), out.end() );//ascending


    capd::vectalg::Matrix<T,0,0> transposedPermutationMatrix = capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<T,0,0> permutationMatrix = capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    //TODO: important!!! we assume that there are no two elements with the same value and all of them are greater then 0
    for( int i = 0; i < in.size(); i++ )
    {
        if( in[i] == out[i] )
        {
            transposedPermutationMatrix[i][i] = 1;
        }
        else
        {
            int j = 0;
            for( j = 0; j < in.size(); j++ )
            {
                if( in[i] == out[j] ) break;
            }
            transposedPermutationMatrix[i][j] = 1;
        }
    }

    permutationMatrix = capd::vectalg::transpose( transposedPermutationMatrix );

    (*diagonalMatrix) = permutationMatrix * (*diagonalMatrix) * transposedPermutationMatrix;
    std::cout << "diagonalMatrix: " << (*diagonalMatrix) << std::endl;

    (*inverseOfP) = permutationMatrix * (*inverseOfP);
    (*P) = (*P) * transposedPermutationMatrix;

    matrixForNonlinearPart = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    (*matrixForNonlinearPart) = (*inverseOfP) * (*inverseOfMassMatrix);

    remainsFromDiagonal = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            if( i != j )
            {
                (*remainsFromDiagonal)[i][j] += (*diagonalMatrix)[i][j];
                (*diagonalMatrix)[i][j] = 0;
            }
            else
            {
                (*diagonalVector)[i] = (*diagonalMatrix)[i][j];
            }
        }

    std::cout << "diag" << std::endl;
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        std::cout << (*diagonalVector)[i] << std::endl;

    delete eigenValueRealPart;
    delete eigenValueImPart;
    delete eigenVectorRealPart;
    delete eigenVectorImPart;
}


template <typename T>
void FEM<T>::diagonalize_simple()
{
    capd::vectalg::Vector<double,0> *eigenValueRealPart = new capd::vectalg::Vector<double,0>( mesh->numberOfFreeNodes() );
    capd::vectalg::Vector<double,0> *eigenValueImPart = new capd::vectalg::Vector<double,0>( mesh->numberOfFreeNodes() );

    capd::vectalg::Matrix<double,0,0> *eigenVectorRealPart = new capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<double,0,0> *eigenVectorImPart = new capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );


    capd::vectalg::Matrix<double,0,0> *tmp = new capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            (*tmp)[i][j] = (*inverseOfMassMatrixTimesStiffnessMatrix)[i][j].mid().leftBound();
        }

    capd::alglib::computeEigenvaluesAndEigenvectors( (*tmp), (*eigenValueRealPart), (*eigenValueImPart), (*eigenVectorRealPart), (*eigenVectorImPart) );

    delete tmp;

    diagonalMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    P = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    inverseOfP = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            (*P)[i][j] = (*eigenVectorRealPart)[i][j];
        }

    (*inverseOfP) = capd::matrixAlgorithms::gaussInverseMatrix( *P );

    (*diagonalMatrix) = (*inverseOfP) * (*inverseOfMassMatrixTimesStiffnessMatrix) * (*P);

    delete eigenValueRealPart;
    delete eigenValueImPart;
    delete eigenVectorRealPart;
    delete eigenVectorImPart;
}


template <typename T>
void FEM<T>::printAllData()
{
    if( stiffnessMatrix != NULL)
    {
        std::cout << "Stiffnes Matrix:" << std::endl;
        std::cout << *stiffnessMatrix << std::endl;
    }

    if( massMatrix != NULL )
    {
        std::cout << "Mass matrix:" << std::endl;
        std::cout << *massMatrix << std::endl;
    }

    if( loadVector != NULL )
    {
        std::cout << "Load Vector:" << std::endl;
        std::cout << *loadVector << std::endl;
    }

    if( nonLinearVector != NULL )
    {
        std::cout << "Non Linear Vector:" << std::endl;
        std::cout << *nonLinearVector << std::endl;
    }

    if( stiffnessMatrix != NULL && massMatrix != NULL )
    {
        std::cout << "Mass^-1 * Stiffness" << std::endl;
        std::cout << capd::matrixAlgorithms::gaussInverseMatrix( *massMatrix ) * (*stiffnessMatrix) << std::endl;
    }

    std::cout << "FEM: number of nodes (free + constrained): " << mesh->listOfNodes.size() << std::endl;
    std::cout << "FEM: number of nodes (free): " << mesh->numberOfFreeNodes() << std::endl;
    std::cout << "FEM: number of elements: " << mesh->listOfElements.size() << std::endl;
    //std::cout <<  << std::endl;


//    capd::vectalg::Matrix<T,0,0> *stiffnessBeforeCutMatrix;
//    capd::vectalg::Vector<T,0> *solutionVector;//U
//    capd::vectalg::Matrix<T,0,0> *massBeforeCutMatrix;

}


template <typename T>
void FEM<T>::EigenVectorsValues()
{
    //std::cout << capd::matrixAlgorithms::symMatrixDiagonalize(*stiffnessMatrix,*tmpMatrix, 1e-30) << std::endl;
    if( stiffnessMatrix == NULL || massMatrix == NULL ) return;

    std::cout << "Stiffness Matrix" << std::endl;
    CAPD_EigenvectorsEigenvalues( stiffnessMatrix );

    std::cout << "Mass Matrix" << std::endl;
    CAPD_EigenvectorsEigenvalues( massMatrix );

    if( inverseOfMassMatrixTimesStiffnessMatrix == NULL )
    {
        inverseOfMassMatrixTimesStiffnessMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
        if( inverseOfMassMatrix == NULL ) invertMassMatrix();
        (*inverseOfMassMatrixTimesStiffnessMatrix) = (*inverseOfMassMatrix) * (*stiffnessMatrix);
    }

    std::cout << "Mass^-1 * Stiffness Matrix" << std::endl;
    CAPD_EigenvectorsEigenvalues( inverseOfMassMatrixTimesStiffnessMatrix );
}

template <typename T>
void FEM<T>::CAPD_EigenvectorsEigenvalues( capd::vectalg::Matrix<T,0,0> *matrix )
{
    std::cout << "+CAPD eigen" << std::endl;

    capd::vectalg::Vector<T,0> *eigenRealPart = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );
    capd::vectalg::Vector<T,0> *eigenImPart = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );

    capd::vectalg::Matrix<T,0,0> *eigenVectorRealPart = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<T,0,0> *eigenVectorImPart = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    capd::alglib::computeEigenvaluesAndEigenvectors( (*matrix), (*eigenRealPart), (*eigenImPart), (*eigenVectorRealPart), (*eigenVectorImPart) );

    //tmp:
    std::cout << "Eigenvalues (Real part)" << std::endl;
    std::cout << (*eigenRealPart) << std::endl;

    std::cout << "Eigenvalues (Img part)" << std::endl;
    std::cout << (*eigenImPart) << std::endl;

    std::cout << "Eigenvectors" << std::endl;
    std::cout <<  capd::alglib::eigenvectorsToString( (*eigenVectorRealPart), (*eigenVectorImPart) ) << std::endl;

    std::cout << "+CAPD eigen" << std::endl;

    delete eigenRealPart;
    delete eigenImPart;
    delete eigenVectorRealPart;
    delete eigenVectorImPart;
}

template <typename T>
bool FEM<T>::courantCondition( double deltaT )
{
    if( stiffnessMatrix == NULL && massMatrix == NULL ) return false;

    capd::vectalg::Matrix<T,0,0> *massStiffness = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    (*massStiffness) = capd::matrixAlgorithms::gaussInverseMatrix( *massMatrix ) * (*stiffnessMatrix);

    capd::vectalg::Vector<T,0> *eigenRealPart = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );
    capd::vectalg::Vector<T,0> *eigenImPart = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );

    capd::alglib::computeEigenvalues( (*massStiffness), (*eigenRealPart), (*eigenImPart) );

    //TODO: attention - we are checking only real part of the vector
    T maxValue = capd::abs( (*eigenRealPart)[0] );
    for( int i = 0; i < (*eigenRealPart).dimension(); i++ )
    {
        maxValue = capd::max( capd::abs( (*eigenRealPart)[i] ), maxValue );
    }

    delete eigenRealPart;
    delete eigenImPart;
    delete massStiffness;

    if( deltaT < ( 0.5 * (1 / maxValue) ) ) return true;
    else return false;
}


///////////////////////////////////////////////////////////
template <typename T>
void FEM<T>::diagonalize_for_candidate()
{
    capd::vectalg::Vector<double,0> *eigenValueRealPart = new capd::vectalg::Vector<double,0>( mesh->numberOfFreeNodes() );
    capd::vectalg::Vector<double,0> *eigenValueImPart = new capd::vectalg::Vector<double,0>( mesh->numberOfFreeNodes() );

    capd::vectalg::Matrix<double,0,0> *eigenVectorRealPart = new capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<double,0,0> *eigenVectorImPart = new capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );


    capd::vectalg::Matrix<double,0,0> tmp = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            //make double matrix from interval matrix
            tmp[i][j] = mid( (*inverseOfMassMatrixTimesStiffnessMatrix)[i][j] );
        }

    capd::jaco::EigenTools<double,0> et;

    capd::vectalg::Matrix<double,0,0> schurT = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<double,0,0> schurS = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<double,0,0> inverseOfSchurS = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    //this procedure will fill matrices T(upper triangle) and S,  T = S^{-1}MS, where M = tmp
    et.schurDecompose( tmp, schurT, schurS );
    
    std::cout << "T" << schurT << std::endl;
    std::cout << "S "<< schurS << std::endl;
    
    inverseOfSchurS = capd::matrixAlgorithms::gaussInverseMatrix( schurS );
    
    capd::vectalg::Matrix<T,0,0> intervalInverseOfSchurS = capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            intervalInverseOfSchurS[i][j] = (T)inverseOfSchurS[i][j];
        }

    capd::alglib::computeEigenvaluesAndEigenvectors( schurT, (*eigenValueRealPart), (*eigenValueImPart), (*eigenVectorRealPart), (*eigenVectorImPart) );

    capd::vectalg::Matrix<double,0,0> E = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            E[i][j] = (*eigenVectorRealPart)[i][j];
        }


    capd::vectalg::Matrix<double,0,0> inverseOfE = capd::vectalg::Matrix<double,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    inverseOfE = capd::matrixAlgorithms::gaussInverseMatrix( E );

    capd::vectalg::Matrix<T,0,0> intervalInverseOfE = capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            intervalInverseOfE[i][j] = (T)inverseOfE[i][j];
        }

    diagonalMatrix = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    P = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    inverseOfP = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    (*inverseOfP) = intervalInverseOfE * intervalInverseOfSchurS;

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            tmp[i][j] = mid( (*inverseOfP)[i][j] );
        }

    tmp = capd::matrixAlgorithms::gaussInverseMatrix( tmp );

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            (*P)[i][j] = (T)tmp[i][j];
        }

    (*diagonalMatrix) = (*inverseOfP) * (*inverseOfMassMatrixTimesStiffnessMatrix) * (*P);

    remainsFromDiagonal = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    diagonalVector = new capd::vectalg::Vector<T,0>( mesh->numberOfFreeNodes() );

    std::vector<T> in;
    std::vector<T> out;

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        in.push_back( (*diagonalMatrix)[i][i] );

    out = in;
    //std::sort( out.begin(), out.end(), std::greater<T>() );//descending
    std::sort( out.begin(), out.end() );//ascending


    capd::vectalg::Matrix<T,0,0> transposedPermutationMatrix = capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    capd::vectalg::Matrix<T,0,0> permutationMatrix = capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );

    //TODO: important!!! we assume that there are no two elements with the same value
    for( int i = 0; i < in.size(); i++ )
    {
        if( in[i] == out[i] )
        {
            transposedPermutationMatrix[i][i] = 1;
        }
        else
        {
            int j = 0;
            for( j = 0; j < in.size(); j++ )
            {
                if( in[i] == out[j] ) break;
            }
            transposedPermutationMatrix[i][j] = 1;
        }
    }

    permutationMatrix = capd::vectalg::transpose( transposedPermutationMatrix );

    (*diagonalMatrix) = permutationMatrix * (*diagonalMatrix) * transposedPermutationMatrix;

    (*inverseOfP) = permutationMatrix * (*inverseOfP);
    (*P) = (*P) * transposedPermutationMatrix;

    matrixForNonlinearPart = new capd::vectalg::Matrix<T,0,0>( mesh->numberOfFreeNodes(), mesh->numberOfFreeNodes() );
    (*matrixForNonlinearPart) = (*inverseOfP) * (*inverseOfMassMatrix);

    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        for( int j = 0; j < mesh->numberOfFreeNodes(); j++ )
        {
            if( i != j )
            {
                (*remainsFromDiagonal)[i][j] += (*diagonalMatrix)[i][j];
                (*diagonalMatrix)[i][j] = 0;
            }
            else
            {
                (*diagonalVector)[i] = (*diagonalMatrix)[i][j];
            }
        }

    std::cout << "diag" << std::endl;
    for( int i = 0; i < mesh->numberOfFreeNodes(); i++ )
        std::cout << (*diagonalVector)[i] << std::endl;

    delete eigenValueRealPart;
    delete eigenValueImPart;
    delete eigenVectorRealPart;
    delete eigenVectorImPart;
}

#endif // _FEM_HPP_
