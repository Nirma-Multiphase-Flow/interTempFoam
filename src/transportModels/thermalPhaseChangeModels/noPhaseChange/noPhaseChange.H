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
    Foam::thermalPhaseChangeModel::noPhaseChange

Description
    Dummy model for no phase change. Just set none as the model type.

SourceFiles
    noPhaseChange.C

\*--------------------------------------------------------------------*/

#ifndef noPhaseChange_H
#define noPhaseChange_H

#include "thermalPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                         Class noPhaseChange
\*--------------------------------------------------------------------*/

class noPhaseChange
:
    public thermalPhaseChangeModel
{

// Private Member Functions


    //- Phase change heat
    volScalarField Q_pc_;

public:

    //- Runtime type information
    TypeName("noPhaseChange");


    // Constructors

    //- construct from components
    noPhaseChange
    (
        const word& name,
        const dictionary& dict,
        const immiscibleIncompressibleTwoPhaseMixture& mixture,
        const volScalarField& T,
        const volScalarField& alpha1
    );



    virtual ~noPhaseChange() = default;

    virtual tmp<volScalarField> Q_pc()       const override;
    virtual tmp<volScalarField> PCV()        const override;
    virtual tmp<volScalarField> alpha1Gen()  const override;
    virtual void correct()                        override;
    virtual bool read(const dictionary& dict)     override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
