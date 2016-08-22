#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <list>
#include <algorithm>
#include <time.h>


#include "../Mesh/Mesh.hpp"

#include "../FEM/FEM.hpp"

#include "../Inclusion/Det.hpp"

#include "../FadMaps/FEMTaylor.h"
#include "../FadMaps/Inclusion/Data.hpp"
#include "../FadMaps/Inclusion/Head.hpp"
#include "../FadMaps/Inclusion/Perturbation.hpp"

#include "../Inclusion/DPDEInclRect2Set.h"
#include "../Inclusion/DPDEMultiMap.h"
#include "../Inclusion/DPDEInclusionCW.h"

#include "../Equation/BurgersNWave.hpp"

#include "../Inclusion/Algorithms2TailManager.h"
#include "../Inclusion/D1Real.h"
#include "../Inclusion/InclusionStepControl.hpp"

#include "capd/intervals/Interval.h"
#include "capd/poincare/TimeMap.hpp"


bool zawieranie( const capd::IVector &before, const capd::IVector &after )
{
    bool rtn = true;

    for( int i = 0; i < before.dimension(); i++ )
        if( ( before[i].leftBound() > after[i].leftBound() ) || ( before[i].rightBound() < after[i].rightBound() ) ) return false;

    return rtn;
}

template <typename Scalar>
capd::IVector intersection( const capd::IVector &lastValue, const capd::IVector &iv )
{
    capd::IVector rtn( iv.dimension() );

    for( int i = 0; i < lastValue.dimension(); i++ )
        capd::intervals::intersection( lastValue[i], iv[i], rtn[i] );

    return rtn;
}

int main( int argc, char *argv[] )
{
    std::cout.precision(15);
    if( argc == 10 )
    {
        typedef capd::intervals::Interval< double >  Scalar;
        typedef BurgersNWave<Scalar>  Equation;
        typedef Det<Scalar> PDET;

        typedef capd::jaco::PolyBd< Scalar, capd::jaco::D1Real<Scalar>, 0> PolyBd1Real;
        typedef capd::jaco::Algorithms2TailManager< PolyBd1Real, PDET> AlgorithmT;


        typedef Head<Scalar, 0, AlgorithmT>  FadMap;
        typedef Perturbation<Scalar,0, AlgorithmT> Perturbation;


        //typedef FEMTaylor< FadMap, capd::dynsys::NoStepControl >  FTaylor;// for fixed time step
        typedef FEMTaylor< FadMap, capd::dynsys::InclusionStepControl >  FTaylor;// for automatic time step


        typedef capd::jaco::DPDEMultiMap< FadMap, Perturbation >  MultiMap;
        typedef capd::jaco::DPDEInclusionCW< MultiMap, MultiMap::TailType, FTaylor> FEMDiffIncl;

        typedef capd::poincare::TimeMap<FEMDiffIncl> TMap;

        time_t curr = time(0);
        std::cout << "Ile free nodow: " << argv[1] << " ile nie dyssypatywnych: " << argv[2] << " rzad: " << argv[3] << " ni: " << argv[4] << " non autonomous parameter: " << argv[5] << " krok: " << argv[6] << " poczatkowy punkt czasowy: " << argv[7] << " hull: " << argv[8] << " initVector: " << argv[9] << " czas startu: " << ctime(&curr) << std::endl;

        int numberOfNodes = atoi(argv[1]);
        int nonDissipativeMods = atoi(argv[2]);
        int order = atoi(argv[3]);
        double ni = atof(argv[4]);
        double nonAutonomousParameter = atof(argv[5]);
        double timeStep = atof(argv[6]);
        double startTime = atof(argv[7]);
        double hull = atof(argv[8]);
        std::string str = argv[9];

        //Mesh
        Mesh<Scalar> mesh;
        mesh.make1DMesh( numberOfNodes, -10, 10 ); //dimension == 1
        mesh.renumerateNodes();

        //FEM

        FEM<Scalar> fem;
        fem.mesh = &mesh;
        fem.dimension = mesh.getDimension();

        Equation equation(hull, ni, nonAutonomousParameter);
        fem.equation = &equation;
        equation.fem = &fem;

        fem.computeCutCoordinateVector();
        fem.assembleStiffnessMatrix();
        fem.assembleMassMatrix();
        fem.invertMassMatrix();
        fem.invertsMassMatrixTimesStiffnessMatrix();
        fem.computeLoadVector();
        fem.diagonalize();


        double lambdaThreshold;
        lambdaThreshold = leftBound( (*fem.diagonalVector)[nonDissipativeMods-1] + ( (*fem.diagonalVector)[nonDissipativeMods] - (*fem.diagonalVector)[nonDissipativeMods-1] ) / 2 );


        std::cout << "\nStep: " << timeStep << " Initial condition interval: " << hull << std::endl;
        std::cout << "Lambda threshold: " << lambdaThreshold << " Ni: " << ni << std::endl;
        std::cout << "Non autonomous parameter: " << nonAutonomousParameter << std::endl;
        std::cout << std::endl;


        equation.assembleNonautonomousVector();//must be called after fem computeCutCoordinateVector
        //std::cout << "NAV: " << equation.NonautonomousVector << std::endl;

        //Numerics

        Data<Scalar> data( fem, lambdaThreshold );
        data.equation = &equation;

        FadMap map( data.numberOfNonDissipativeModes, data.numberOfModes,-equation.ni,0,4, data);
        Perturbation inclusion( data.numberOfNonDissipativeModes, data.numberOfModes,-equation.ni,0,4, data);

        capd::vectalg::Vector<Scalar,0> initialVector( data.numberOfModes );
        equation.candidate( &initialVector, str );

        capd::vectalg::Vector<Scalar,0> initTime;
        initTime = initialVector;

        std::cout << "Candidate in Diag base with hull: " << initialVector << std::endl;

        PolyBd1Real poly( 0, 4, data.numberOfNonDissipativeModes, data.numberOfModes );
        for( int i = 0; i < data.numberOfModes; i++ )
        {
            poly.set( i, initialVector[i] );
        }
        inclusion.setT0( poly );

        capd::vectalg::Vector<Scalar,0> initialVector1( data.numberOfNonDissipativeModes );
        initialVector1 = data.downsizeVector( initialVector );

        capd::jaco::DPDEInclRect2Set< capd::vectalg::Matrix<Scalar,0,0> > C0Set( initialVector1 );


        capd::IEuclNorm norm;
        MultiMap multiMap(map, inclusion);
        FEMDiffIncl solver( multiMap, order, timeStep, norm);

        fem.mesh->printFreeNodesCoordinates();
        std::cout << "Dbase: " << startTime << " " << initialVector << std::endl;
        std::cout << "Time: " << startTime << " " << initTime << std::endl;

        int j = 1;
        int periodMultiplier = 3;
        int periodCounter = 1;
        double i = 0;
        int timeStepCounter = 0;

        C0Set.setCurrentTime( startTime );


    capd::IVector iv;
        while( C0Set.getCurrentTime() < 2 * Scalar::pi() + startTime )
        {
            C0Set.move(solver, i, std::cout);

            //if( (j++ % 1000) == 0 )//uncomment while fixed time step is chosen, so only every 1000th time step will be printed
            {
                iv = (capd::IVector)C0Set;
                poly = inclusion.getT0();
                std::cout << "Step: " << solver.getStep() << std::endl;
                std::cout << "Diam: " << capd::vectalg::size( data.P * data.joinHeadWithTail(iv,poly) ) << std::endl;
                std::cout << "Dbase: " << C0Set.getCurrentTime() << " " << data.joinHeadWithTail(iv,poly) << std::endl;
                std::cout << "Time: " << C0Set.getCurrentTime() << " " << data.P * data.joinHeadWithTail(iv,poly) << std::endl;
                j = 1;
            }
            i++;
            timeStepCounter++;
            Scalar nextStep = solver.getStepControl().computeNextTimeStep(solver, C0Set.getCurrentTime(), C0Set, 1e-1);
            nextStep = capd::min( nextStep, 2 * Scalar::pi() + startTime - C0Set.getCurrentTime() );
            solver.setStep( nextStep );
        }
        std::cout << "No. of time steps: " << timeStepCounter << std::endl;
        std::cout << "Diag max diam: " << capd::vectalg::size( data.joinHeadWithTail(iv,poly) ).leftBound() << std::endl;
        std::cout << "Subset: " << zawieranie( initialVector, data.joinHeadWithTail(iv,poly) ) << std::endl;
        if( !zawieranie( initialVector, data.joinHeadWithTail(iv,poly) ) )
        {
            std::cout << "Intsec:" << intersection(initialVector, data.joinHeadWithTail(iv,poly) ) << std::endl;
        }
        mesh.clearMemory();
    }else{std::cout << "zla liczba paramatrow" << std::endl;}

    return 0;
}
