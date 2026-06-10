if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        mesh.newIOobject("rhoPhiSum"),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );
    surfaceScalarField rhoCpPhiSum
    (
        mesh.newIOobject("rhoCpPhiSum"),
        mesh,
        dimensionedScalar(rhoCpPhi.dimensions(), Zero)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        rhoCpPhiSum += (runTime.deltaT()/totalDeltaT)*rhoCpPhi;
    }

    rhoPhi = rhoPhiSum;
    rhoCpPhi = rhoCpPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;
rhoCp ==
    alpha1*rho1*cp1
  + (scalar(1) - alpha1)*rho2*cp2;