scalar error = GREAT;

label i = 0;
while (error > Ttolerance && i < TIter)
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

    rhok = mask*(1.0 - beta*(T - TRef))*rhoL + (1-mask)*rhoS;

    volScalarField h = cp*T;

    volScalarField Tstar = (Tl - Ts)*lf + Ts;

    // when (h - (cp*Tstar) < epsilon, DHnew == DH
    // the temperature converged

    volScalarField diff = h - (cp*Tstar);
    volScalarField absDiff = mag(diff)*mask;
    error = gMax(absDiff);
    DHnew = mask*(DH + glRelax*diff);
    i++;
}
Info << "TIter = " << i <<" max T = " << gMax(T) << tab << "min T = " << gMin(T) << tab <<" err: "<< error<< nl;
