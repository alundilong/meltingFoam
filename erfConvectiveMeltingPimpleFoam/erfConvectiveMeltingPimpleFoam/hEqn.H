{
    fvScalarMatrix hEqn
    (
        fvm::ddt(cp, T)
      + fvm::div(phi*fvc::interpolate(cp), T)
      + mask*hs*exp(-pow((T-Tmelt)/Tdim,2))/Foam::sqrt(mathematical::pi)/Tdim*fvm::ddt(T)
      + mask*hs*exp(-pow((T-Tmelt)/Tdim,2))/Foam::sqrt(mathematical::pi)/Tdim*(U & fvc::grad(T))
      - fvm::laplacian(lambda/rho, T)
    );

    hEqn.relax();
    hEqn.solve();

    alpha = (0.5*Foam::erf((T-Tmelt)/Tdim)+scalar(0.5))*mask + alpha*(1-mask);
    alpha.correctBoundaryConditions();

    cp = (alpha*cpL+(1.0-alpha)*cpS)*mask + cpPore*(1-mask);
    cp.correctBoundaryConditions();

    lambda = (alpha*lambdaL+(1.0-alpha)*lambdaS)*mask + lambdaPore*(1-mask);
    lambda.correctBoundaryConditions();

    nu = (alpha*nuL+(1.0-alpha)*nuS)*mask + nuPore*(1-mask);
    nu.correctBoundaryConditions();

    rhok = (1.0 - Foam::max(beta*(T - Tl),0.0))*mask + rhoPore/rho*(1-mask);
    rhok.correctBoundaryConditions();

    DC = DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs);
    DC.correctBoundaryConditions();
}
