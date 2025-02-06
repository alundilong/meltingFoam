   IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // solid -> phase 1
    // liquid -> phase 2

    // Reading density rho
    dimensionedScalar rho(transportProperties.lookup("rho"));

    // Reading thermal conductivity lambda
    dimensionedScalar lambdaS(transportProperties.lookup("lambdaS"));
    dimensionedScalar lambdaL(transportProperties.lookup("lambdaL"));

    // Reading heat capacity cp
    dimensionedScalar cpS(transportProperties.lookup("cpS"));
    dimensionedScalar cpL(transportProperties.lookup("cpL"));

    // Reading kinematic viscosity
    dimensionedScalar nuS(transportProperties.lookup("nuS"));
    dimensionedScalar nuL(transportProperties.lookup("nuL"));

    // Reading latent heat of fusion hs
    dimensionedScalar hs(transportProperties.lookup("hs"));

    // Reading solid bound of melting temperature Ts
    dimensionedScalar Ts(transportProperties.lookup("Ts"));

    // Reading liquid bound of melting temperature Tl
    dimensionedScalar Tl(transportProperties.lookup("Tl"));

    // Reading phase front parameter Tdim
    dimensionedScalar Tdim(transportProperties.lookup("Tdim"));

    // Reading volume expansion factor beta
    dimensionedScalar beta(transportProperties.lookup("beta"));

    // Reading large D'arcy-type source term constant DCl
    dimensionedScalar DCl(transportProperties.lookup("DCl"));

    // Reading small D'arcy-type source term constant DCs
    dimensionedScalar DCs(transportProperties.lookup("DCs"));
