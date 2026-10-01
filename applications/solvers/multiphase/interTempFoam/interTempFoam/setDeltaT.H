/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    scalar maxDeltaTFact =
        min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));

    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    // Capillary (Brackbill-Kothe-Zemach) limit: dt <= maxCapillaryNum*dt_c
    const scalar maxCapillaryNum =
        runTime.controlDict().getOrDefault<scalar>("maxCapillaryNum", 1);

    // Phase-change limit: dt <= maxPhaseChangeCo*h/max(|mdot''|/rho1)
    const scalar maxPhaseChangeCo =
        runTime.controlDict().getOrDefault<scalar>("maxPhaseChangeCo", 0.25);
    const scalar uPC = phaseChangePtr->maxInterfaceSpeed();
    const scalar dtPC =
        uPC > SMALL ? maxPhaseChangeCo*gMin(interface.cellSize())/uPC : GREAT;

    runTime.setDeltaT
    (
        min
        (
            min(deltaTFact*runTime.deltaTValue(), maxDeltaT),
            min(maxCapillaryNum*interface.capillaryDeltaT(), dtPC)
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
