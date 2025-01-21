    // Solve the momentum equation

    //- for stablization purpose
    volScalarField divRhoPhi
    (
        "divRhoPhi",
        fvc::div(phi)*rho
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho,U) + fvm::div(rhoPhi, U)
      + fvm::Sp(rho*C*pow(1-lf, 2)/(rhoRef*(pow(lf, 3)+q)), U) //added
      + turbulence->divDevRhoReff(rho,U)
      + fvm::SuSp(divRhoPhi, U) + divRhoPhi*U
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

    }
