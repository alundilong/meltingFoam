// Reading fields
Info<< "Reading fields\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Reading transport properties
Info<< "Reading thermophysical properties\n" << endl;
#include "readTransportProperties.H"

//- setup mask
// mask field for phase change zone (1 for cell meltable, 0 non meltable)
volScalarField mask
(
    IOobject
    (
        "mask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless,scalar(initialMaskValue))
);

CSVReader csvReader(fileName(transportProperties.lookup("PorousFile")));
scalar coordinateScaler(readScalar(transportProperties.lookup("coordinateScaler")));
bool loaded = csvReader.loadCSV();
if (loaded)
{
    const std::vector<CSVReader::DataRow> &dataRows = csvReader.getData();
    for (const auto& row : dataRows)
    {
        const scalar &x = row.point0;
        const scalar &y = row.point1;
        const scalar &z = row.point2;
        //Info << x << tab << y << tab << z << nl;
        label celli = mesh.findCell(point(x*coordinateScaler,y*coordinateScaler,z*coordinateScaler));
        if (celli != -1)
        {
            mask.primitiveFieldRef()[celli] = scalar(1.0);
        }
    }
}
mask.correctBoundaryConditions();
mask.write();

#include "createPhi.H"

// Calculatinf fit-parameters for phase change function
Info<< "Calculating phase change properties\n" << endl;

dimensionedScalar Tmelt
(
    "Tmelt", 
    (Tl+Ts)/2.0
);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    (1.0 - Foam::max(beta*(T - Tl),0.0))*mask + rhoPore/rho*(1-mask)
);
rhok.correctBoundaryConditions();

// D'arcy-type source term field
volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh
    ),
    DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs)
);
DC.correctBoundaryConditions();

// Thermal conductivity field
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh
    ),
    (alpha*lambdaL+(1.0-alpha)*lambdaS)*mask + lambdaPore*(1-mask)
);
lambda.correctBoundaryConditions();

// Heat capacity field
volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh
    ),
    (alpha*cpL+(1.0-alpha)*cpS)*mask + cpPore*(1-mask)
);
cp.correctBoundaryConditions();

// Kinematic viscosity field
volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh
    ),
    (alpha*nuL+(1.0-alpha)*nuS)*mask + nuPore*(1-mask)
);
nu.correctBoundaryConditions();

Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("ghf", g & mesh.Cf());

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    mesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

Info << "Writing mesh info" << endl;
// writeout mesh info
// C, magSf, weights, deltaCoeffs
// S, Cf, delta
mesh.C().write();
mesh.magSf().write();
mesh.weights().write();
mesh.deltaCoeffs().write();
mesh.Sf().write();
mesh.Cf().write();
mesh.delta().ref().write();

volVectorField U0("U0",U);
U0.write();
