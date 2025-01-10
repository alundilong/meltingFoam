for (int i=0; i<10; i++)
{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", pos(mask)*(turbulence->nu()/Pr + alphat) + pos(1.0-mask)*alphapore);

    #include "sources.H"

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      - fvc::ddt(DH)/Cp*pos(mask)
      - fvc::div(phi, DH)/Cp*pos(mask)
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);

    volScalarField h = Cp*T;

    volScalarField Tstar = (Tl - Ts)*lf + Ts;

    DHnew = DH + glRelax*(h - (Cp*Tstar));
}
