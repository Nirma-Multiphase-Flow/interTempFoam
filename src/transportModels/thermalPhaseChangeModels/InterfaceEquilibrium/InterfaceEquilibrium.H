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

Class
    Foam::thermalPhaseChangeModel::InterfaceEquilibrium

Description
    This model forces wall and interface cells to the equilibrium state,
    causing evaporation or condensation as appropriate. The phase change
    heat term is limited to not produce/sink more than the volume of a 
    cell in a single time step(like the courantAlpha limit).

    This model uses different interface thresholds for identifying
    evaporating and condensing cell pairs.

    This is similar to the previous HiLoFacePair code, but the phase change
    rate is numerically relaxed by some user-defined value (should be <= 1,
    where 0 is no phase change, and 1 is no relaxation). This seems to help
    with overshooting where the temperature in a cell can go from superheated
    to subcooled in one time step.


SourceFiles
    InterfaceEquilibrium.C

\*--------------------------------------------------------------------*/

#ifndef InterfaceEquilibrium_H
#define InterfaceEquilibrium_H

#include "thermalPhaseChangeModel.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*---------------------------------------------------------------------------*\
                   Class InterfaceEquilibrium Declaration
\*---------------------------------------------------------------------------*/

class InterfaceEquilibrium
:
    public thermalPhaseChangeModel
{
    const fvMesh& mesh_;

    volScalarField Q_pc_;          // phase-change heat [W/m^3], AUTO_WRITE
    volScalarField InterfaceField_; // interface indicator {0,1}, AUTO_WRITE
    volScalarField WallField_;      // wall-adjacent mask {0,1}

    scalar CondThresh_;   // alpha1 threshold for condensation-side detection
    scalar EvapThresh_;   // alpha1 threshold for evaporation-side detection
    scalar RelaxFac_;     // under-relaxation factor in (0,1]

    void calcQ_pc();

public:

    TypeName("InterfaceEquilibrium");

    InterfaceEquilibrium
    (
        const word& name,
        const dictionary& dict,
        const immiscibleIncompressibleTwoPhaseMixture& mixture,
        const volScalarField& T,
        const volScalarField& alpha1
    );

    virtual ~InterfaceEquilibrium() = default;

    virtual tmp<volScalarField> Q_pc() const override;
    virtual void correct()            override;
    virtual bool read(const dictionary& dict) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // thermalPhaseChangeModels
} // Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //