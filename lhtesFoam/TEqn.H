for (int i=0; i<10; i++)
{
    #include "sources.H"

    cp = mask*(lf*cpL+(1-lf)*cpS) + (1-mask)*cpPore;
    rho = mask*(lf*rhoL+(1-lf)*rhoS)+(1-mask)*rhoPore;
    volScalarField rhoCp = rho*cp;

    k = mask*(lf*kL+(1-lf)*kS)+(1-mask)*kPore;
    surfaceScalarField rhoCpPhi = fvc::interpolate(rhoCp)*phi;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp,T)
      + fvm::div(rhoCpPhi, T)
      - fvm::laplacian(k, T)
     ==
      - fvc::ddt(rho,DH)*mask
      - fvc::div(rhoPhi,DH)*mask
    );

    TEqn.relax();

    TEqn.solve();

    rhok = (1.0 - beta*(T - TRef))*rhoL;

    volScalarField h = cp*T;

    volScalarField Tstar = (Tl - Ts)*lf + Ts;

    DHnew = DH + glRelax*(h - (cp*Tstar));
}
