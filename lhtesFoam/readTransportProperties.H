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

dimensionedScalar rhoL
(
    materialProperties.lookup("rhoL")
);

dimensionedScalar rhoS
(
    materialProperties.lookup("rhoS")
);

dimensionedScalar rhoPore
(
    materialProperties.lookup("rhoPore")
);

dimensionedScalar kL
(
    materialProperties.lookup("kL")
);

dimensionedScalar kS
(
    materialProperties.lookup("kS")
);

dimensionedScalar kPore
(
    materialProperties.lookup("kPore")
);

dimensionedScalar mu
(
    materialProperties.lookup("mu")
);

dimensionedScalar cpL
(
    materialProperties.lookup("cpL")
);

dimensionedScalar cpS
(
    materialProperties.lookup("cpS")
);

dimensionedScalar cpPore
(
    materialProperties.lookup("cpPore")
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

label TIter(materialProperties.lookupOrDefault("TIter",10));
scalar Ttolerance(materialProperties.lookupOrDefault("tolerance",1e-2));
