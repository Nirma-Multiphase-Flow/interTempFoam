#include "createRDeltaT.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());


// -----------------------------------------------------------------
// Geometric VoF: one PLIC reconstruction per step, owned by isoAdvection
// (settings in the alpha solver dictionary of system/fvSolution).
// -----------------------------------------------------------------

// Must be constructed BEFORE the advector: plicRDF owns an internal
// reconstructedDistanceFunction also named "RDF"; the solver's RDF must win
// the registry name because the legacy phase-change model looks it up.
reconstructedDistanceFunction RDF(mesh);

// Total advecting flux/velocity of the interface, u_t = u + u_PC (Gada &
// Sharma 2009 Eq. 13; Shaikh et al. 2016 Eq. 8). isoAdvection keeps
// references to these two fields; alphaEqn.H refreshes them before advect().
// phiPC is the phase-change part of the flux (kept for the interface Courant
// number and the rhoPhi correction).
surfaceScalarField phiT
(
    IOobject("phiT", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    phi
);

volVectorField UT
(
    IOobject("UT", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    U
);

surfaceScalarField phiPC
(
    IOobject("phiPC", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

isoAdvection advector(alpha1, phiT, UT);

// Interface geometry (phase-signed distance psiRDF, curvature, CSF force);
// curvatureModel gradAlpha | RDF in transportProperties
interfaceCurvatureITF interface(mixture, alpha1);

zoneDistribute& exchangeFields =
    zoneDistribute::New(mesh);

Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


const word pcDictName("phaseChangeProperties");
autoPtr<thermalPhaseChangeModel> phaseChangePtr;

// Kept at outer scope: the model holds a reference to this dictionary.
autoPtr<IOdictionary> pcDictPtr;

if (mesh.time().db().found(pcDictName) || isFile(mesh.time().constant()/pcDictName))
{
    pcDictPtr.reset
    (
        new IOdictionary
        (
            IOobject
            (
                pcDictName,
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );

    phaseChangePtr = thermalPhaseChangeModel::New
    (
        pcDictName,
        *pcDictPtr,
        mixture,
        T,
        alpha1
    );
}
else
{
    Info<< "No phaseChangeProperties found -- using noPhaseChange\n" << endl;

    dictionary noDict;
    noDict.add("model", word("noPhaseChange"));
    noDict.add("T_sat",
        dimensionedScalar("T_sat", dimTemperature, scalar(373.15)).value());
    noDict.add("h_lv",
        dimensionedScalar("h_lv", dimEnergy/dimMass, scalar(2.26e6)).value());

    phaseChangePtr = thermalPhaseChangeModel::New
    (
        pcDictName,
        noDict,
        mixture,
        T,
        alpha1
    );
}

const dimensionedScalar& cp1 = phaseChangePtr->cp1();
const dimensionedScalar& cp2 = phaseChangePtr->cp2();

// Phase-change treatment in the alpha equation (alpha1 dict in fvSolution):
//   phaseChangeAdvection velocity;  interface moved by u_t = u + u_PC,
//                                   u_PC = mdot'' n_RDF/rho1 (default), or
//                       source;     explicit volumetric source alpha1Gen
const word pcAdvection
(
    mesh.solverDict(alpha1.name()).getOrDefault<word>
    (
        "phaseChangeAdvection",
        "velocity"
    )
);
const bool pcVelocity(pcAdvection == "velocity");

if (!pcVelocity && pcAdvection != "source")
{
    FatalErrorInFunction
        << "phaseChangeAdvection '" << pcAdvection
        << "' unknown; use velocity or source" << exit(FatalError);
}

Info<< "Alpha phase change: " << pcAdvection << nl << endl;

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)
    *(
        phi

    )
);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),

    alpha1*rho1*cp1
  + (scalar(1)-alpha1)*rho2*cp2
);
rhoCp.oldTime();

Info<< "Creating rhoCpPhi\n" << endl;

surfaceScalarField rhoCpPhi
(
    IOobject
    (
        "rhoCpPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),

    fvc::interpolate(rhoCp)*phi
);

Info<< "rhoCp min/max = "
    << gMin(rhoCp)
    << " "
    << gMax(rhoCp)
    << endl;

typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}



mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// Thermal phase-change model -- read from constant/phaseChangeProperties.
// Falls back to noPhaseChange if the dictionary is absent.


#include "createMRF.H"
#include "createFvOptions.H"
