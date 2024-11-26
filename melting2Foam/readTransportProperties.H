singlePhaseTransportModel laminarTransport(U, phi);


// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);


IOdictionary materialProperties
(
    IOobject
    (
         "materialProperties",
         runTime.constant(),
         mesh,
         IOobject::MUST_READ_IF_MODIFIED,
         IOobject::NO_WRITE
    )
);

dimensionedScalar rhoRef
(
    materialProperties.lookup("rhoRef")
);

dimensionedScalar k
(
    materialProperties.lookup("k")
);

dimensionedScalar mu
(
    materialProperties.lookup("mu")
);

dimensionedScalar Cp
(
    materialProperties.lookup("Cp")
);

dimensionedScalar q
(
    materialProperties.lookup("q")
);

dimensionedScalar Tmelt
(
    materialProperties.lookup("Tmelt")
);

dimensionedScalar L
(
    materialProperties.lookup("L")
);

dimensionedScalar C
(
    materialProperties.lookup("C")
);

dimensionedScalar Tl
(
    materialProperties.lookup("Tl")
);

dimensionedScalar Ts
(
    materialProperties.lookup("Ts")
);

dimensionedScalar glRelax
(
    materialProperties.lookup("glRelax")
);

