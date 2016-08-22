#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <list>
#include <algorithm>
#include <time.h>

#include "../Mesh/Mesh.hpp"

#include "../FEM/FEM.hpp"

#include "../FadMaps/CandidateFadMap.hpp"

#include "../Equation/BurgersNWave_candidate.hpp"

#include "../Inclusion/InclusionStepControl.hpp"

#include "capd/intervals/Interval.h"
#include "capd/poincare/TimeMap.hpp"
#include "capd/dynsys/BasicFadTaylor.hpp"
#include "capd/dynsys/TaylorHOE.hpp"


int main( int argc, char *argv[] )
{
std::cout.precision(15);
    if( argc == 7 )
    {
        typedef double  Scalar;
        typedef BurgersNWaveCandidate<Scalar>  Equation;

        time_t curr = time(0);
        std::cout << "Ile free nodow: " << argv[1] << " rzad: " << argv[2] << " ni: " << argv[3] << " non autonomous parameter: " << argv[4] << " krok: " << argv[5] << " dokad w czasie: " << argv[6] << " czas startu: " << ctime(&curr) << std::endl;


        int numberOfNodes = atoi(argv[1]);
        int order = atoi(argv[2]);
        double ni = atof(argv[3]);
        double nonAutonomousParameter = atof(argv[4]);
        double timeStep = atof(argv[5]);
        double howFarInTime = atof(argv[6]);

        bool automaticStep;
        if( timeStep == 0 )
        {
            automaticStep = true;
            timeStep = 1e-5;
        }
        else automaticStep = false;


        if( automaticStep ) std::cout << "\nStep: auto" << std::endl;
        else std::cout << "\nStep: " << timeStep << std::endl;

        std::cout << "Ni: " << ni << std::endl;
        std::cout << "Non autonomous parameter: " << nonAutonomousParameter << std::endl;
        std::cout << std::endl;

        //Mesh
        Mesh<Scalar> mesh;
        mesh.make1DMesh( numberOfNodes, -10, 10 ); //dimension == 1
        mesh.renumerateNodes();

        //FEM

        FEM<Scalar> fem;
        fem.mesh = &mesh;
        fem.dimension = mesh.getDimension();

        Equation equation( ni, nonAutonomousParameter );
        fem.equation = &equation;
        equation.fem = &fem;

        fem.computeCutCoordinateVector();
        fem.assembleStiffnessMatrix();
        fem.assembleMassMatrix();
        fem.invertMassMatrix();
        fem.invertsMassMatrixTimesStiffnessMatrix();
        fem.computeLoadVector();

        equation.assembleNonautonomousVector();//must be called after fem computeCutCoordinateVector
        //std::cout << "NAV: " << equation.NonautonomousVector << std::endl;

        //std::cout << "Mass matrix - A:" << *fem.massMatrix << std::endl;
        //std::cout << "Stiffness matrix - K:" << *fem.stiffnessMatrix << std::endl;
        //std::cout << "A^-1:" << *fem.inverseOfMassMatrix << std::endl;
        //std::cout << "A^-1 * K:" << *fem.inverseOfMassMatrixTimesStiffnessMatrix << std::endl;


        //Numerics
        CandidateFadMap<Scalar, 0, Equation> map;
        map.fem = &fem;
        map.equation = &equation;
        //capd::dynsys::BasicFadTaylor< CandidateFadMap<Scalar, 0, Equation>, capd::dynsys::NoStepControl > solver( map, 10, 0.0001 );
        capd::dynsys::BasicFadTaylor< CandidateFadMap<Scalar, 0, Equation>, capd::dynsys::DLastTermsStepControl > solver( map, order, timeStep );
        //capd::dynsys::TaylorHOE< CandidateFadMap<Scalar, 0, Equation> > solver( map, order, timeStep );

        capd::vectalg::Vector<Scalar,0> initialVector( fem.mesh->numberOfFreeNodes() );
        equation.setInitialVector( &initialVector, &(mesh.listOfNodes) );

        std::cout << "InitialVector: " << initialVector << std::endl;

        fem.mesh->printFreeNodesCoordinates();

        capd::poincare::TimeMap< capd::dynsys::BasicFadTaylor< CandidateFadMap<Scalar, 0, Equation>, capd::dynsys::DLastTermsStepControl> > timeMap( solver );
        timeMap.stopAfterStep(false);
        if( automaticStep ) timeMap.turnOnStepControl();
        else timeMap.turnOffStepControl();

    int numberOfPeriods = 0;
    double pi = 3.14159265359;
    double period =  2.0 * pi;

    double eps = 1e-9;

    capd::vectalg::Vector<Scalar,0> iv, iv1;
    iv = initialVector;
    
    bool run = true;
    while( run )
    {
        do
        {
            iv1 = timeMap( period, iv );
            std::cout << "Step: " << timeMap.getStep() << std::endl;
            std::cout << "Time: " << timeMap.getCurrentTime() << " " << iv1 << std::endl;
        }while( !timeMap.completed() );
        numberOfPeriods++;

        run = false;
        for( int i = 0; i < numberOfNodes; i++ )
        {
            if( capd::abs( iv[i] - iv1[i] ) > eps ) run = true;
        }
        iv = iv1;
    }

    std::cout << "Periods: " << numberOfPeriods << std::endl;

    mesh.clearMemory();
    }else{std::cout << "zla liczba paramatrow" << std::endl;}

    return 0;
}
