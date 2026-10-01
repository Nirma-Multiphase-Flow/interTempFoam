/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::interfaceCurvatureITF

Description
    Interface geometry for interTempFoam, built once per outer iteration from
    the isoAdvector PLIC reconstruction of the advected alpha:

    - psiRDF: signed distance to the PLIC interface, psi > 0 in phase 1,
      valid on the interface cells (ring 0) and 3 point-neighbour rings.
      Each ring takes the facets carried by the previous ring; the distance
      is the weighted plane distance of reconstructedDistanceFunction.
      Only "genuine" interface cells (a face neighbour on the other side of
      alpha = 0.5) contribute facets; isolated wisps are counted and ignored.
    - aInterface: PLIC area density |surf.normal()|/V of genuine cells [1/m].
    - kappaI: interface curvature, constant across the interface:
        kappa_j = -div(grad(psi)/|grad(psi)|)       (rings 0-1)
        kappa_i = weighted mean over the stencil of interface cell i of
                  kappa_j shifted to the interface,
                  kappa_j/(1 + kappa_j psi_j/(D - 1)),
                  weights 1/(|psi_j| + 0.1 h)
        ring-1 cells take kappa of the nearest interface cell.
    - surfaceTensionForce(): sigma kappa_f snGrad(alpha1), with
      kappa_f = (w_P k_P + w_N k_N)/(w_P + w_N), w = alpha(1 - alpha) + 1e-6
      ("RDF" model), or the stock interfaceProperties force ("gradAlpha").

    - heightFunction: Popinet (2009) height functions on interface cells
      whose 3-column stencil (up to 7 layers, full end layers of opposite
      phase) lies on this processor in a structured part of a 2D planar or
      axisymmetric wedge mesh; the RDF value elsewhere. On wedges the radial
      heights are volume weighted and the hoop term -n_r/r is added.

    Selected by "curvatureModel gradAlpha | RDF | heightFunction | constant;"
    in transportProperties (default gradAlpha). "constant" applies the prescribed
    "constantCurvature" on rings 0-1 (balanced-force verification only).

SourceFiles
    interfaceCurvatureITF.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCurvatureITF_H
#define interfaceCurvatureITF_H

#include "volFields.H"
#include "surfaceFields.H"
#include "reconstructionSchemes.H"
#include "reconstructedDistanceFunction.H"
#include "surfaceTensionModel.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interfaceCurvatureITF Declaration
\*---------------------------------------------------------------------------*/

class interfaceCurvatureITF
{
    // Private data

        const fvMesh& mesh_;

        const immiscibleIncompressibleTwoPhaseMixture& mixture_;

        const volScalarField& alpha1_;

        //- gradAlpha | RDF | heightFunction | constant
        word model_;

        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- In-plane cell size V/max(non-wedge/empty face area) [m]
        scalarField h_;

        //- Number of curvature directions: 1 (2D planar) or 2 (3D, wedge)
        scalar nCurvDir_;

        //- Axisymmetric wedge mesh
        bool wedge_;

        //- In-plane unit directions (2D planar: the two solution
        //  directions; wedge: e1_ = radial, e2_ = axis). Unused in 3D.
        vector e1_;
        vector e2_;

        //- Smallest in-plane cell size for the capillary time-step limit [m]
        scalar hCapMin_;

        //- Point-neighbour ring level of the interface band (-1 outside)
        volScalarField level_;

        //- Signed distance, > 0 in phase 1
        volScalarField psi_;

        //- PLIC interface area density of genuine interface cells [1/m]
        volScalarField aInterface_;

        //- Extended interface curvature [1/m]
        volScalarField kappa_;

        //- Genuine interface cells
        boolList interfaceCell_;

        //- Diagnostics of the last correct()
        label nInterface_;
        label nWisp_;

        //- Prescribed curvature of curvatureModel constant (verification)
        scalar constantKappa_;


    // Private Member Functions

        void calcCellSize();

        //- Mark genuine interface cells, count wisps
        void markInterfaceCells(const reconstructionSchemes& surf);

        //- Height-function curvature on interface cells with a valid
        //  structured column stencil (2D planar or axisymmetric wedge);
        //  overwrites kappa_ there, other cells keep the RDF value.
        //  Returns the number of cells with a valid height function.
        label heightFunctionKappa(const reconstructionSchemes& surf);

        //- Face neighbour of celli across the internal face whose outward
        //  normal is aligned with e (cos > 0.99); -1 if none
        label neighbourAlong(const label celli, const vector& e) const;

        //- Build psi on rings 0-3 from the PLIC facets
        void calcPsi
        (
            const reconstructionSchemes& surf,
            reconstructedDistanceFunction& RDF
        );

        //- Curvature from psi, interface averaging and extension
        void calcKappa
        (
            const reconstructionSchemes& surf,
            reconstructedDistanceFunction& RDF
        );


public:

    // Constructors

        interfaceCurvatureITF
        (
            const immiscibleIncompressibleTwoPhaseMixture& mixture,
            const volScalarField& alpha1
        );

        interfaceCurvatureITF(const interfaceCurvatureITF&) = delete;
        void operator=(const interfaceCurvatureITF&) = delete;


    // Member Functions

        //- Update from the current PLIC reconstruction. Must be called after
        //  surf.reconstruct() on the advected alpha and before
        //  surfaceTensionForce(). Also marks the RDF band (3 rings) so the
        //  legacy RDF can be built on it.
        void correct
        (
            const reconstructionSchemes& surf,
            reconstructedDistanceFunction& RDF
        );

        const word& model() const noexcept { return model_; }

        const volScalarField& psi() const noexcept { return psi_; }

        const volScalarField& aInterface() const noexcept
        {
            return aInterface_;
        }

        const volScalarField& kappa() const noexcept { return kappa_; }

        const boolList& interfaceCell() const noexcept
        {
            return interfaceCell_;
        }

        const scalarField& cellSize() const noexcept { return h_; }

        //- Surface tension face force [N/m^3 * m^2 / m^2] as used in UEqn/pEqn
        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Brackbill-Kothe-Zemach capillary time-step limit
        //  sqrt(rho_avg h_min^3/(2 pi sigma_max)); GREAT when sigma = 0
        scalar capillaryDeltaT() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
