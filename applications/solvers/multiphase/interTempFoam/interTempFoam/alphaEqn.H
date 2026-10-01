{
    // Temporarily making U relative to mesh motion
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    // Source terms handed to the geometric advection
    volScalarField::Internal Su
    (
        IOobject("Su", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("Su", dimless/dimTime, Zero)
    );
    volScalarField::Internal Sp
    (
        IOobject("Sp", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("Sp", dimless/dimTime, Zero)
    );

    if (pcVelocity)
    {
        // Interface moved by u_t = u + u_PC (alphaPCVelocity.H)
        #include "alphaPCVelocity.H"
        phiT == phi + phiPC;
        UT == U + UPC;
    }
    else
    {
        // Phase-change volume source (alphaSuSp.H)
        #include "alphaSuSp.H"
        phiPC == dimensionedScalar(phiPC.dimensions(), Zero);
        phiT == phi;
        UT == U;
    }

    advector.advect(Sp, Su);

    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    #include "rhofs.H"
    // getRhoPhi uses phiT = phi + phiPC; remove the phase-change part of the
    // carrier flux so momentum and enthalpy see the bulk flux
    rhoPhi = advector.getRhoPhi(rho1f, rho2f) - rho2f*phiPC;

    // Enthalpy flux consistent with the geometric alpha flux
    rhoCpPhi = advector.getRhoPhi(rho1*cp1, rho2*cp2) - rho2*cp2*phiPC;

    // Requested vs applied phase-change source (isoAdvector bounding can
    // drop part of it)
    {
        const scalarField& V = mesh.V();
        const scalar rDt = 1/runTime.deltaTValue();
        const scalar requested =
            gSum(phaseChangePtr->alpha1Gen()().primitiveField()*V);
        const scalar implicitForm =
            gSum((Su.field() + Sp.field()*alpha1.primitiveField())*V);
        const scalar applied = gSum
        (
            (
                (alpha1.primitiveField() - alpha1.oldTime().primitiveField())
               *rDt
              + fvc::surfaceIntegrate(advector.alphaPhi())().primitiveField()
            )*V
        );
        Info<< "alpha source [m3/s]: requested " << requested
            << ", implicit form " << implicitForm
            << ", applied " << applied << endl;
    }

    alpha2 = 1.0 - alpha1;
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
