{
    fvScalarMatrix hEqn
    (
        fvm::ddt(cp, T)
      + fvm::div(phi*fvc::interpolate(cp), T)
      + hs*exp(-pow((T-Tmelt)/Tdim,2))/Foam::sqrt(mathematical::pi)/Tdim*fvm::ddt(T)
      + hs*exp(-pow((T-Tmelt)/Tdim,2))/Foam::sqrt(mathematical::pi)/Tdim*(U & fvc::grad(T))
      - fvm::laplacian(lambda/rho, T)
    );

    hEqn.relax();
    hEqn.solve();

    alpha = 0.5*Foam::erf((T-Tmelt)/Tdim)+scalar(0.5);

    cp = alpha*cpL+(1.0-alpha)*cpS;
    lambda = alpha*lambdaL+(1.0-alpha)*lambdaS;
    nu = alpha*nuL+(1.0-alpha)*nuS;
    rhok = 1.0 - Foam::max(beta*(T - Tl),0.0);
    DC = DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs);
}
