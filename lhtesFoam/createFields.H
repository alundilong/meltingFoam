Info<< "Reading porous structure\n" << endl;

volScalarField lf
(
    IOobject
    (
        "lf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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
    dimensionedScalar("zero",dimless,scalar(0.0))
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

#include "createPhi.H"

#include "readTransportProperties.H"

CSVReader csvReader(fileName(materialProperties.lookup("PorousFile")));
scalar coordinateScaler(readScalar(materialProperties.lookup("coordinateScaler")));
csvReader.loadCSV();
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
mask.correctBoundaryConditions();
mask.write();

Info<< "Reading thermophysical properties\n" << endl;

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

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mask*(lf*rhoL+(1-lf)*rhoS)+(1-mask)*rhoPore
);

volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mask*(lf*cpL+(1-lf)*cpS)+(1-mask)*cpPore
);

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mask*(lf*kL+(1-lf)*kS)+(1-mask)*kPore
);

surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
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
    (1.0 - beta*(T - TRef))*rhoL
);

volScalarField h = rho*cp*T;

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField DH 
(
    IOobject
    (
        "DH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh/*,
    dimensionedScalar("DH", dimLength*dimLength/(dimTime*dimTime), 0)*/
);
DH.oldTime();


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


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
    pimple.dict(),
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

mesh.setFluxRequired(p_rgh.name());

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
