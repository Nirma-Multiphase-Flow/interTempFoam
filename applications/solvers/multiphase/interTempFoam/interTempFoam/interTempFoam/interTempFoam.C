/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Application
    interTempFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with energy equation for heat transfer. Optional mesh motion and mesh
    topology changes including adaptive re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "thermalPhaseChangeModel.H"
#include "reconstructionSchemes.H"
#include "reconstructedDistanceFunction.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, non-isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing with energy equation.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            surf.reconstruct();

            // Diagnostic — kept for monitoring, NOT passed to constructRDF
            {
                volScalarField magGradAlpha(mag(fvc::grad(alpha1)));
                scalar gMaxAlpha = gMax(magGradAlpha);
                Info<< "gMax(|grad(alpha1)|) = " << gMaxAlpha << endl;
                label nInterface = 0;
                forAll(magGradAlpha, celli)
                    if (magGradAlpha[celli] > 1e-3*gMaxAlpha) ++nInterface;
                Info<< "interfaceMask cells = "
                    << returnReduce(nInterface, sumOp<label>()) << endl;
            }

            
            const scalar alphaTolRDF = 1e-4;
            {
                boolList filteredInterfaceCell(mesh.nCells(), false);
                forAll(alpha1, celli)
                {
                    if (   alpha1[celli] > alphaTolRDF
                        && alpha1[celli] < scalar(1) - alphaTolRDF)
                    {
                        filteredInterfaceCell[celli] = true;
                    }
                }
                RDF.markCellsNearSurf(filteredInterfaceCell, 2);
            }
            RDF.constructRDF
            (
                RDF.nextToInterface(),
                surf.centre(),
                surf.normal(),
                exchangeFields,
                true
            );

            // Bootstrap fallback: seed RDF from (alpha-0.5)/|grad(alpha)| when
            // no PLIC surface cells exist (step-function IC or first iteration).
            // Threshold matches interfaceMask (1e-3 * gMax) to restrict seeding
            // to interface-adjacent cells only and avoid near-zero division.
            if (gMax(mag(RDF.primitiveField())) < SMALL)
            {
                const volVectorField gAlpha(fvc::grad(alpha1));
                const scalar gMaxAlphaBS =
                    max(gMax(mag(gAlpha.primitiveField())), SMALL);
                const scalar gradThresh = 1e-3 * gMaxAlphaBS;

                forAll(RDF, celli)
                {
                    const scalar magG = mag(gAlpha[celli]);
                    if (magG > gradThresh)
                        RDF[celli] = (alpha1[celli] - scalar(0.5)) / magG;
                }
                RDF.correctBoundaryConditions();
                Info<< "RDF bootstrap: seeded from grad(alpha), "
                    << "max(|RDF|) = "
                    << gMax(mag(RDF.primitiveField())) << " m" << endl;
            }

            {
                label nBand = 0;
                forAll(RDF.nextToInterface(), celli)
                    if (RDF.nextToInterface()[celli]) ++nBand;
                Info<< "RDF: nextToInterface = "
                    << returnReduce(nBand, sumOp<label>())
                    << " cells, max(|RDF|) = "
                    << gMax(mag(RDF.primitiveField())) << " m" << endl;
            }

            // Phase-change evaluated once per timestep (firstIter only).
            // Calling correct() on every outer iteration causes surf.reconstruct()
            // to run on a partially-moved alpha each time, producing inconsistent
            // PLIC normals between iterations and a standing-wave instability in
            // phiStefan.  The single-per-timestep evaluation is consistent with
            // the explicit-source treatment of latent heat in interFoam-family
            // solvers and introduces at most a one-timestep lag in T coupling.
            if (pimple.firstIter())
            {
                phaseChangePtr->correct();
            }

            phiTotal =
                phi
              + (coldPhaseIsHighAlpha1 ? scalar(-1) : scalar(1))
              * phaseChangePtr->phiStefan();

            
            mixture.correct();


            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            // Update mixture thermal properties with current alpha after
            // alphaEqn has moved the interface and pEqn has updated phi.
            rho =
                alpha1*rho1
            + (scalar(1) - alpha1)*rho2;

            rhoCp =
                alpha1*rho1*cp1
            + (scalar(1) - alpha1)*rho2*cp2;


            #include "TEqn.H"

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
