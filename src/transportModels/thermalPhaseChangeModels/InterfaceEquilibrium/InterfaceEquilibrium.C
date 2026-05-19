/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "InterfaceEquilibrium.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(InterfaceEquilibrium, 0);
    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel, 
        InterfaceEquilibrium, 
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::InterfaceEquilibrium::InterfaceEquilibrium
(
        const word& name,
        const dictionary& dict,
        const immiscibleIncompressibleTwoPhaseMixture& mixture,
        const volScalarField& T,
        const volScalarField& alpha1
)
:
    thermalPhaseChangeModel
    (
        name, 
        dict,
        mixture, 
        T, alpha1
    ),
    mesh_(T.mesh()),
    Q_pc_
    (
        IOobject
        (
            "PhaseChangeHeat",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "Q_pc", dimensionSet(1,-1,-3,0,0,0,0), Zero)   
    ),
    
    InterfaceField_
    (
        IOobject
        (
            "InterfaceField",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("I_gamma", dimless, Zero)
    ),
    WallField_
    (
        IOobject
        (
            "WallField",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("wall", dimless, Zero)
    ),

    CondThresh_(dict.lookupOrDefault<scalar>("CondThresh", 0.95)),
    EvapThresh_(dict.lookupOrDefault<scalar>("EvapThresh", 0.05)),
    RelaxFac_(dict.lookupOrDefault<scalar>("RelaxFac", 1.0))
{
    // Build wall mask once -- mesh topology does not change
    for (const polyPatch& pp : mesh_.boundaryMesh())
    {
        if (isA<wallPolyPatch>(pp))
        {
            for (const label ci : pp.faceCells())
            {
                WallField_.primitiveFieldRef()[ci] = scalar(1);
            }
        }
    }

    correct();
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::InterfaceEquilibrium::calcQ_pc()
{
    // ------------------------------------------------------------------
    // 1. Interface band indicator I_gamma (dimensionless {0,1})
    // ------------------------------------------------------------------
    {
        const scalarField& a1 = alpha1_.primitiveField();
        scalarField& Ig       = InterfaceField_.primitiveFieldRef();

        forAll(a1, ci)
        {
            Ig[ci] =
                (a1[ci] > EvapThresh_ && a1[ci] < CondThresh_)
                ? scalar(1)
                : scalar(0);
        }

        // Wall cells also participate (for condensation detection near hot walls)
        const scalarField& wf = WallField_.primitiveField();
        forAll(wf, ci)
        {
            if (wf[ci] > scalar(0.5))
            {
                Ig[ci] = scalar(1);
            }
        }
    }

    // ------------------------------------------------------------------
    // 2. Time step and phase constants
    // ------------------------------------------------------------------
    const dimensionedScalar dT
    (
        "dT",
        dimTime,
        max(mesh_.time().deltaTValue(), SMALL)
    );

    const dimensionedScalar& rho1 = mixture_.rho1();
    const dimensionedScalar& rho2 = mixture_.rho2();
    const dimensionedScalar& cp1  = mixture_.cp1();
    const dimensionedScalar& cp2  = mixture_.cp2();

    // ------------------------------------------------------------------
    // 3. Unlimited equilibrium heat Q_pc* = I_gamma * rhoCp * (T-Tsat)/dt
    //    Units: [kg/(m^3)] * [J/(kg*K)] * [K/s] = [W/m^3]
    // ------------------------------------------------------------------
    const volScalarField rhoCp
    (
        alpha1_*rho1*cp1 + (scalar(1) - alpha1_)*rho2*cp2
    );

    Q_pc_ = InterfaceField_ * rhoCp * (T_ - T_sat_) / dT;

    // ------------------------------------------------------------------
    // 4. Fluid-availability limiters
    //    LimCond [W/m^3]: max condensation heat = all vapour condenses in dT
    //    LimEvap [W/m^3]: max evaporation heat  = all liquid evaporates in dT
    // ------------------------------------------------------------------
    const volScalarField LimCond( (scalar(1) - alpha1_)*rho2*h_lv_/dT );
    const volScalarField LimEvap( alpha1_*rho1*h_lv_/dT );

    const volScalarField Q_fluid =
        neg(Q_pc_)*max(Q_pc_, -LimCond)
      + pos(Q_pc_)*min(Q_pc_,  LimEvap);

    // ------------------------------------------------------------------
    // 5. Volume-generation (CFL-like) limiter
    //    PCV_fac = dt * (Q_pc/h_lv) * (1/rho2 - 1/rho1)  [dimless]
    //    Require |PCV_fac| <= 1: scale Q_pc by min(1/|PCV_fac|, 1).
    // ------------------------------------------------------------------
    const volScalarField PCV_fac =
        dT*(Q_pc_/h_lv_)*(scalar(1)/rho2 - scalar(1)/rho1);

    const dimensionedScalar eps("eps", dimless, SMALL);

    const volScalarField Q_vol =
        Q_pc_ * min
        (
            dimensionedScalar("one", dimless, scalar(1)),
            dimensionedScalar("one", dimless, scalar(1)) / (mag(PCV_fac) + eps)
        );

    // ------------------------------------------------------------------
    // 6. Composite limit: most restrictive of the three
    //    Condensation (Q_pc < 0): max picks least-negative = smallest |rate|
    //    Evaporation  (Q_pc > 0): min picks smallest value
    // ------------------------------------------------------------------
    Q_pc_ =
        neg(Q_pc_)*max(max(Q_pc_, Q_fluid), Q_vol)
      + pos(Q_pc_)*min(min(Q_pc_, Q_fluid), Q_vol);

    // ------------------------------------------------------------------
    // 7. Under-relaxation
    // ------------------------------------------------------------------
    Q_pc_ *= RelaxFac_;

    Info<< "Phase-change energy integral: "
        << gSum(Q_pc_.primitiveField() * mesh_.V()) << " W" << nl
        << "Volume-change integral: "
        << gSum((Q_pc_/h_lv_*(scalar(1)/rho2 - scalar(1)/rho1))().primitiveField()*mesh_.V())
        << " m^3/s" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::InterfaceEquilibrium::Q_pc() const
{
    return Q_pc_;
}


void Foam::thermalPhaseChangeModels::InterfaceEquilibrium::correct()
{
    calcQ_pc();
}


bool Foam::thermalPhaseChangeModels::InterfaceEquilibrium::read
(
    const dictionary& dict
)
{
    thermalPhaseChangeModel::read(dict);
    CondThresh_ = dict.lookupOrDefault<scalar>("CondThresh", 0.95);
    EvapThresh_ = dict.lookupOrDefault<scalar>("EvapThresh", 0.05);
    RelaxFac_   = dict.lookupOrDefault<scalar>("RelaxFac",   1.0);
    return true;
}


// ************************************************************************* //