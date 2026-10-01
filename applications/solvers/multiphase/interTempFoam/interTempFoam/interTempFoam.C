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
#include "localEulerDdtScheme.H"
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
#include "isoAdvection.H"
#include "interfaceCurvatureITF.H"

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

    if (LTS)
    {
        FatalErrorInFunction
            << "Local time stepping (localEuler) is not supported:"
            << " isoAdvector and the phase-change source need a uniform"
            << " time step. Use ddtSchemes { default Euler; }."
            << exit(FatalError);
    }

    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        // Explicit surface tension is unstable beyond the capillary limit
        // (V2 at dt = 36 s: 15000x over it, round-off grew into a 2D mode)
        if
        (
            !adjustTimeStep
         && runTime.deltaTValue() > interface.capillaryDeltaT()
         && runTime.timeIndex() < 1
        )
        {
            WarningInFunction
                << "fixed deltaT " << runTime.deltaTValue()
                << " exceeds the capillary time-step limit "
                << interface.capillaryDeltaT()
                << ": explicit surface tension is unstable." << endl;
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

            // PLIC of the advected (new) alpha; advect() reconstructed the
            // start-of-step alpha. Then mixture properties (and the stock
            // curvature used by curvatureModel gradAlpha), then the interface
            // geometry, which also marks the RDF band for the legacy RDF.
            advector.surf().reconstruct();
            mixture.correct();
            interface.correct(advector.surf(), RDF);

            RDF.constructRDF
            (
                RDF.nextToInterface(),
                advector.surf().centre(),
                advector.surf().normal(),
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

            // Sharp models: interfacial mass transfer from the converged
            // temperature, used by the next time step
            if (pimple.finalIter())
            {
                phaseChangePtr->correctAfterT();
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        Info<< "Total mass sum(rho V) = "
            << gSum(rho.primitiveField()*mesh.V()) << " kg" << endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
