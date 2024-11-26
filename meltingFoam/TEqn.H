for (int i=0; i<10; i++)
{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    #include "sources.H"

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      - fvc::ddt(DH)/Cp
      - fvc::div(phi, DH)/Cp
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.setValues(cellSet, cellSetTemp);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);

    volVectorField gradT = fvc::grad(T);
        
    volScalarField dTdx = gradT.component(vector::X);

    volScalarField h = Cp*T;

    volScalarField Tstar = (Tl - Ts)*lf + Ts;

    DHnew = DH + glRelax*(h - (Cp*Tstar));
}
