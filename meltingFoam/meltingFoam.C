/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantBoussinesqPimpleFoam

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "noRadiation.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "CSVReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	//Before entering the time loop decalre the hashTable
    	HashTable<scalar, scalar, Hash<scalar>> tlfValues;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

		//update the latent heat content after each loop        
		volScalarField DHnew = DH;


        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

		//Inside the time loop
	    scalar t = runTime.value();
		scalar tlf = lf.weightedAverage(mesh.V()).value();
		tlfValues.insert(t, tlf);

        runTime.write();

 		Info<< nl << endl;
		Info<< "Liquid fraction [-] = " << tlf << endl;
       	Info<< nl << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s";
		Info<< nl << endl;
        Info<< "  ClockTime = " << runTime.elapsedClockTime() << " s";
		Info<< nl << endl;
    }
    if (Pstream::master())
    {
	    // Outside the time loop, save the hash table to a file
	    OFstream myfile("./tlfvalues");
	    myfile() << "time\t" << "tlf" << endl;

	    //write the tlfValues to the file sorted by time
	    for (scalar ti : tlfValues.sortedToc())
	    {
	    	myfile() << ti << "\t" << tlfValues[ti] << endl;
	    }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
