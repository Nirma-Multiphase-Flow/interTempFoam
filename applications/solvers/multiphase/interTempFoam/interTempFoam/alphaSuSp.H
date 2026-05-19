// Phase-change sources for MULES (volScalarField::Internal, [1/s])
// Su  = alpha1Gen: positive for condensation, negative for evaporation
// Sp  = 0 (no implicit correction for interface-equilibrium model)
// divU = PCV: volume dilatation source, positive for evaporation
tmp<volScalarField> tDivU(phaseChangePtr->PCV());

// Phase-change sources for MULES
tmp<volScalarField> tSu
(
    volScalarField::New
    (
        "zeroSu",
        mesh,
        dimensionedScalar("zero", dimless/dimTime, Zero)
    )
);


volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimless/dimTime,
    tSu().primitiveField()
);

volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Sp", dimless/dimTime, Zero)
);

volScalarField::Internal divU
(
    IOobject
    (
        "divU",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimless/dimTime,
    tDivU().primitiveField()
);