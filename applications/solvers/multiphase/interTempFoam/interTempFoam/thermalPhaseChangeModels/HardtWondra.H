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
    Foam::thermalPhaseChangeModels::HardtWondra

Description
    Sharp-interface thermal phase change for interTempFoam.

    Conventions: psi (field "psiRDF", from interfaceCurvatureITF) > 0 in
    phase 1; mdot > 0 converts phase 1 into phase 2; L = h_lv (1 -> 2) > 0.

    Temperature (ghost fluid, Gibou et al. 2002): on every face f whose two
    cells lie on different sides of the interface,
        theta_P = psi_P/(psi_P - psi_N)   clipped to [thetaMin, 1]
        a_P     = k_P |S_f|/(theta_P d)   (k_P of the phase of cell P)
    the conduction across f is removed (kappaf = 0) and cell P gets the
    implicit Dirichlet term a_P (T_P - Tsat) (same for N).

    Mass transfer, evaluated with the converged T of the step (one-step lag):
        Q_f = a_P (T_P - Tsat) + a_N (T_N - Tsat)       [W]
        mdot_f = Q_f/L                                   [kg/s]
    assigned to the interface cell(s) of the pair (split by theta when both
    are), else to the phase-1 cell; mdot''' = sum mdot_f / V.
    The heat removed from the T equation equals L sum mdot_f by construction.

    Sources:
        alpha:  S_alpha = -mdot'''/rho1 in the receiving cells (alpha1Gen)
        volume: Helmholtz-smoothed mdot''' (length lambdaCells*h), moved into
                the non-interface cells of the lighter phase, rescaled to
                the total, times (1/rho2 - 1/rho1) (Hardt & Wondra 2008
                shift); 0 for equal densities.

Usage
    constant/phaseChangeProperties:
    \verbatim
    model       HardtWondra;
    T_sat       373.15;     // [K]
    h_lv        2.26e6;     // [J/kg], phase 1 -> 2
    k1          0.68;       // phase-1 conductivity [W/m/K]
    k2          0.025;      // phase-2 conductivity [W/m/K]
    thetaMin    0.05;       // optional
    lambdaCells 1;          // optional, dilatation smoothing length / h
    fluxMethod  probe;      // optional, probe | ghostFluid
    \endverbatim

SourceFiles
    HardtWondra.C

\*---------------------------------------------------------------------------*/

#ifndef HardtWondra_H
#define HardtWondra_H

#include "thermalPhaseChangeModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*---------------------------------------------------------------------------*\
                         Class HardtWondra Declaration
\*---------------------------------------------------------------------------*/

class HardtWondra
:
    public thermalPhaseChangeModel
{
    // Private data

        const fvMesh& mesh_;

        dimensionedScalar k1_;
        dimensionedScalar k2_;
        scalar thetaMin_;
        scalar lambdaCells_;

        //- probe (default): mass flux from pure-cell (T - Tsat)/d samples
        //  next to each interface cell (Malan et al. 2021); ghostFluid: from
        //  the ghost-fluid face heat (exactly energy consistent, but uses
        //  mixed-cell temperatures with 1/theta)
        word fluxMethod_;

        //- Interface region count of the last warning (report changes only)
        label nRegionsLast_;

        //- probe method: heat [W] drained by the ghost-fluid pin but not
        //  converted into mass, returned to the probe cells next step
        scalarField heatReturn_;

        //- Volumetric mass-transfer rate, > 0 for phase 1 -> 2 [kg/m3/s]
        volScalarField mdot_;

        //- Dilatation source for pEqn [1/s]
        volScalarField PCV_;

        //- In-plane cell size [m]
        scalarField h_;

        //- |psi| < psiTol_ h is treated as "on the interface"; the phase
        //  side is then taken from alpha1 >= 0.5 - 1e-8 (a centred facet
        //  gives alpha = 0.5 +- round-off: the tie goes to phase 1)
        static constexpr scalar psiTol_ = 1e-6;

        //- Diagnostics of the last TSource() call
        mutable label nCrossing_;
        mutable label nClipped_;


    // Private Member Functions

        const volScalarField& psi() const;

        //- Phase side of every cell: true = phase 1
        boolList cellSide() const;

        //- theta of cell P for a crossing face (unclipped, in (0,1))
        static scalar theta
        (
            const scalar psiP,
            const scalar psiN,
            const scalar aP,
            const scalar aN
        );

        void calcCellSize();

        //- Shifted, normalised dilatation from mdot_
        void calcPCV();

        HardtWondra(const HardtWondra&) = delete;
        void operator=(const HardtWondra&) = delete;


public:

    TypeName("HardtWondra");


    // Constructor

        HardtWondra
        (
            const word& name,
            const dictionary& dict,
            const immiscibleIncompressibleTwoPhaseMixture& mixture,
            const volScalarField& T,
            const volScalarField& alpha1
        );


    //- Destructor
    virtual ~HardtWondra() = default;


    // Member functions

        //- Latent heat rate L mdot''' [W/m3] (diagnostic)
        virtual tmp<volScalarField> Q_pc() const override;

        //- Sharp face conductivity, zero on interface-crossing faces
        virtual tmp<surfaceScalarField> kappaf() const override;

        //- Ghost-fluid Tsat terms
        virtual tmp<fvScalarMatrix> TSource
        (
            const volScalarField& T
        ) const override;

        //- -mdot'''/rho1 [1/s]
        virtual tmp<volScalarField> alpha1Gen() const override;

        //- Shifted dilatation [1/s]
        virtual tmp<volScalarField> PCV() const override;

        //- Nothing at the start of the step (see correctAfterT)
        virtual void correct() override {}

        //- mdot''' and PCV from the converged temperature
        virtual void correctAfterT() override;

        virtual scalar maxInterfaceSpeed() const override;

        virtual bool read(const dictionary& dict) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
