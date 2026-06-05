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
// geometricVoF reconstruction infrastructure
// -----------------------------------------------------------------

IOdictionary reconDict
(
    IOobject
    (
        "reconstructionDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<reconstructionSchemes> surfPtr
(
    reconstructionSchemes::New
    (
        alpha1,
        phi,
        U,
        reconDict
    )
);

reconstructionSchemes& surf = surfPtr();

reconstructedDistanceFunction RDF(mesh);

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

const dimensionedScalar& cp1 = mixture.cp1();
const dimensionedScalar& cp2 = mixture.cp2();

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


const word pcDictName("phaseChangeProperties");
autoPtr<thermalPhaseChangeModel> phaseChangePtr;

// Kept at outer scope so TEqn.H can read flags (e.g. coldPhaseIsHighAlpha1)
// without a registry lookup that would fail after this block exits.
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

surfaceScalarField phiTotal
(
    IOobject
    (
        "phiTotal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phi
);

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


// Read cold-phase flag for TEqn hot-floor gate.
// true  → melting:      solid = alpha1=1 = cold phase
// false → evaporation:  vapor = alpha1=0 = cold phase
// REPLACE the entire coldPhaseIsHighAlpha1 block with:
const bool coldPhaseIsHighAlpha1 =
(
    pcDictPtr.valid()
    ? pcDictPtr->lookupOrDefault<bool>("coldPhaseIsHighAlpha1", false)
    : false
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// Thermal phase-change model -- read from constant/phaseChangeProperties.
// Falls back to noPhaseChange if the dictionary is absent.


#include "createMRF.H"
#include "createFvOptions.H"
