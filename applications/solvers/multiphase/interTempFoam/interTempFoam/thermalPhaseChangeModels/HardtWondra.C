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

\*---------------------------------------------------------------------------*/

#include "HardtWondra.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "regionSplit.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "zoneDistribute.H"
#include "reconstructionSchemes.H"
#include "reconstructedDistanceFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(HardtWondra, 0);

    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel,
        HardtWondra,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::volScalarField&
Foam::thermalPhaseChangeModels::HardtWondra::psi() const
{
    return mesh_.lookupObject<volScalarField>("psiRDF");
}


Foam::boolList
Foam::thermalPhaseChangeModels::HardtWondra::cellSide() const
{
    const scalarField& p = psi().primitiveField();
    const scalarField& a = alpha1_.primitiveField();
    boolList side(mesh_.nCells());
    forAll(side, celli)
    {
        // |psi| below round-off of the facet distance (interface through
        // the cell centre, degenerate class a): decide by alpha, otherwise
        // the side follows the sign of round-off noise and differs from
        // row to row of a planar front
        side[celli] =
            (mag(p[celli]) > psiTol_*h_[celli])
          ? (p[celli] > 0)
          : (a[celli] >= 0.5 - 1e-8);
    }
    return side;
}


Foam::scalar Foam::thermalPhaseChangeModels::HardtWondra::theta
(
    const scalar psiP,
    const scalar psiN,
    const scalar aP,
    const scalar aN
)
{
    if (psiP*psiN < 0)
    {
        return psiP/(psiP - psiN);
    }
    // Interface not resolved by psi (e.g. a step-function initial field
    // with the interface on the face): locate it from alpha
    if (mag(aP - aN) > SMALL)
    {
        return min(max((aP - 0.5)/(aP - aN), scalar(0)), scalar(1));
    }
    return 0.5;
}


void Foam::thermalPhaseChangeModels::HardtWondra::calcCellSize()
{
    h_.setSize(mesh_.nCells());
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(h_, celli)
    {
        scalar maxA = SMALL;
        for (const label facei : mesh_.cells()[celli])
        {
            if (mesh_.isInternalFace(facei))
            {
                maxA = max(maxA, mesh_.magSf()[facei]);
            }
            else
            {
                const label patchi = pbm.whichPatch(facei);
                const polyPatch& pp = pbm[patchi];
                if (isA<wedgePolyPatch>(pp) || isA<emptyPolyPatch>(pp))
                {
                    continue;
                }
                maxA = max
                (
                    maxA,
                    mesh_.magSf().boundaryField()[patchi][facei - pp.start()]
                );
            }
        }
        h_[celli] = mesh_.V()[celli]/maxA;
    }
}


void Foam::thermalPhaseChangeModels::HardtWondra::calcPCV()
{
    PCV_ = dimensionedScalar(PCV_.dimensions(), Zero);

    const scalar rho1 = mixture_.rho1().value();
    const scalar rho2 = mixture_.rho2().value();
    const scalar dRhoInv = 1/rho2 - 1/rho1;
    const scalarField& V = mesh_.V();

    const scalar mTot = gSum(mdot_.primitiveField()*V);
    if (mag(rho1 - rho2) < SMALL*rho1 || gMax(mag(mdot_.primitiveField())) < VSMALL)
    {
        Info<< "HardtWondra: PCV = 0 (equal densities or no mass transfer)"
            << endl;
        return;
    }

    const scalarField& aInt =
        mesh_.lookupObject<volScalarField>("aInterface").primitiveField();

    // Helmholtz smoothing of mdot''' with lambda = lambdaCells*h (mean h of
    // the source cells; constant coefficient for parallel consistency)
    scalar hSum = 0;
    label nSrc = 0;
    forAll(mdot_, celli)
    {
        if (mdot_[celli] != 0)
        {
            hSum += h_[celli];
            ++nSrc;
        }
    }
    reduce(hSum, sumOp<scalar>());
    reduce(nSrc, sumOp<label>());
    const dimensionedScalar lambdaSqr
    (
        dimArea, sqr(lambdaCells_*hSum/max(nSrc, label(1)))
    );

    volScalarField s
    (
        IOobject("mdotSmoothHW", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_, dimensionedScalar(mdot_.dimensions(), Zero), "zeroGradient"
    );
    s.primitiveFieldRef() = mdot_.primitiveField();
    s.correctBoundaryConditions();
    {
        fvScalarMatrix sEqn
        (
            fvm::Sp(scalar(1), s) - fvm::laplacian(lambdaSqr, s) == mdot_
        );
        // fvSolution entry "mdotSmoothHW" if present, else PCG/DIC to 1e-14
        dictionary solverDict;
        if
        (
            mesh_.solutionDict().subDict("solvers").found(s.name())
        )
        {
            solverDict = mesh_.solverDict(s.name());
        }
        else
        {
            solverDict.add("solver", "PCG");
            solverDict.add("preconditioner", "DIC");
            solverDict.add("tolerance", 1e-14);
            solverDict.add("relTol", 0);
            solverDict.add("maxIter", 2000);
        }
        sEqn.solve(solverDict);
    }

    // Lighter-phase, non-interface cells receive the dilatation
    const bool lighterIsPhase1 = rho1 < rho2;
    const boolList side(cellSide());
    scalar den = 0;
    forAll(side, celli)
    {
        if (side[celli] == lighterIsPhase1 && aInt[celli] <= 0)
        {
            den += s[celli]*V[celli];
        }
    }
    reduce(den, sumOp<scalar>());
    if (mag(den) < VSMALL)
    {
        WarningInFunction
            << "no smoothed source in the lighter phase; PCV set to zero"
            << endl;
        return;
    }
    const scalar N = mTot/den;
    forAll(side, celli)
    {
        if (side[celli] == lighterIsPhase1 && aInt[celli] <= 0)
        {
            PCV_[celli] = N*s[celli]*dRhoInv;
        }
    }
    PCV_.correctBoundaryConditions();

    // Single connected interface assumed by the global normalisation
    boolList blocked(mesh_.nFaces(), true);
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(nei, facei)
    {
        blocked[facei] = !(aInt[own[facei]] > 0 && aInt[nei[facei]] > 0);
    }
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        blocked[facei] = false;
    }
    const regionSplit regions(mesh_, blocked);
    labelHashSet ifRegions;
    forAll(aInt, celli)
    {
        if (aInt[celli] > 0) ifRegions.insert(regions[celli]);
    }
    List<labelList> procRegions(Pstream::nProcs());
    procRegions[Pstream::myProcNo()] = ifRegions.toc();
    Pstream::gatherList(procRegions);
    Pstream::scatterList(procRegions);
    labelHashSet allRegions;
    for (const labelList& lst : procRegions)
    {
        allRegions.insert(lst);
    }
    // Face-connected interface cells: a sloped front touching only at cell
    // corners also counts as several regions. Report changes only.
    if (allRegions.size() > 1 && allRegions.size() != nRegionsLast_)
    {
        WarningInFunction
            << allRegions.size() << " disconnected interface regions: the"
            << " global PCV normalisation moves volume between them"
            << " (TODO: region-wise normalisation)" << endl;
    }
    nRegionsLast_ = allRegions.size();

    Info<< "HardtWondra: sum(PCV V) = " << gSum(PCV_.primitiveField()*V)
        << " m3/s, target sum(mdot V)(1/rho2 - 1/rho1) = " << mTot*dRhoInv
        << " m3/s" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::HardtWondra::HardtWondra
(
    const word& name,
    const dictionary& dict,
    const immiscibleIncompressibleTwoPhaseMixture& mixture,
    const volScalarField& T,
    const volScalarField& alpha1
)
:
    thermalPhaseChangeModel(name, dict, mixture, T, alpha1),
    mesh_(T.mesh()),
    k1_("k1", dimPower/dimLength/dimTemperature, Zero),
    k2_("k2", dimPower/dimLength/dimTemperature, Zero),
    thetaMin_(0.05),
    lambdaCells_(1),
    fluxMethod_("probe"),
    nRegionsLast_(1),
    heatReturn_(),
    mdot_
    (
        IOobject("mdot", mesh_.time().timeName(), mesh_,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh_, dimensionedScalar(dimDensity/dimTime, Zero), "zeroGradient"
    ),
    PCV_
    (
        IOobject("PCV", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh_, dimensionedScalar(dimless/dimTime, Zero), "zeroGradient"
    ),
    h_(),
    nCrossing_(0),
    nClipped_(0)
{
    read(dict);
    calcCellSize();
    heatReturn_.setSize(mesh_.nCells(), Zero);

    // Restart: rebuild the dilatation from the stored mdot at the first
    // correctAfterT (the interface geometry does not exist yet here)
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::Q_pc() const
{
    return mdot_*h_lv_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::kappaf() const
{
    const boolList side(cellSide());
    const scalar k1 = k1_.value();
    const scalar k2 = k2_.value();

    auto tkf = surfaceScalarField::New("kappafHW", mesh_, k1_);
    surfaceScalarField& kf = tkf.ref();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(nei, facei)
    {
        const bool sP = side[own[facei]];
        kf[facei] = (sP != side[nei[facei]]) ? 0 : (sP ? k1 : k2);
    }

    const volScalarField& p = psi();
    forAll(kf.boundaryField(), patchi)
    {
        fvsPatchScalarField& pkf = kf.boundaryFieldRef()[patchi];
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        if (p.boundaryField()[patchi].coupled())
        {
            const scalarField pN(p.boundaryField()[patchi].patchNeighbourField());
            const scalarField aN
            (
                alpha1_.boundaryField()[patchi].patchNeighbourField()
            );
            forAll(pkf, i)
            {
                const bool sP = side[fc[i]];
                const bool sN = (mag(pN[i]) > psiTol_*h_[fc[i]]) ? (pN[i] > 0) : (aN[i] >= 0.5 - 1e-8);
                pkf[i] = (sP != sN) ? 0 : (sP ? k1 : k2);
            }
        }
        else
        {
            forAll(pkf, i)
            {
                pkf[i] = side[fc[i]] ? k1 : k2;
            }
        }
    }

    return tkf;
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::thermalPhaseChangeModels::HardtWondra::TSource
(
    const volScalarField& T
) const
{
    auto tM = tmp<fvScalarMatrix>::New(T, dimPower);
    fvScalarMatrix& M = tM.ref();
    scalarField& diag = M.diag();
    scalarField& src = M.source();

    const scalar Tsat = T_sat_.value();
    const scalar k1 = k1_.value();
    const scalar k2 = k2_.value();
    const boolList side(cellSide());
    const scalarField& p = psi().primitiveField();
    const scalarField& a = alpha1_.primitiveField();
    const vectorField& C = mesh_.C();
    const scalarField& magSf = mesh_.magSf();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    nCrossing_ = 0;
    nClipped_ = 0;
    scalar thetaRawMin = 1;

    forAll(nei, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        if (side[P] == side[N]) continue;

        const scalar thP0 = theta(p[P], p[N], a[P], a[N]);
        const scalar thP = max(thP0, thetaMin_);
        const scalar thN = max(1 - thP0, thetaMin_);
        thetaRawMin = min(thetaRawMin, min(thP0, 1 - thP0));
        if (thP0 < thetaMin_ || 1 - thP0 < thetaMin_) ++nClipped_;
        ++nCrossing_;

        const scalar d = mag(C[N] - C[P]);
        const scalar aP = (side[P] ? k1 : k2)*magSf[facei]/(thP*d);
        const scalar aN = (side[N] ? k1 : k2)*magSf[facei]/(thN*d);
        diag[P] += aP;
        src[P] += aP*Tsat;
        diag[N] += aN;
        src[N] += aN*Tsat;
    }

    const volScalarField& pf = psi();
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fp = mesh_.boundary()[patchi];
        if (!fp.coupled()) continue;

        const scalarField pN(pf.boundaryField()[patchi].patchNeighbourField());
        const scalarField aNb
        (
            alpha1_.boundaryField()[patchi].patchNeighbourField()
        );
        const vectorField delta(fp.delta());
        const scalarField& mS = fp.magSf();
        const labelUList& fc = fp.faceCells();
        forAll(fc, i)
        {
            const label P = fc[i];
            const bool sN = (mag(pN[i]) > psiTol_*h_[fc[i]]) ? (pN[i] > 0) : (aNb[i] >= 0.5 - 1e-8);
            if (side[P] == sN) continue;
            const scalar thP0 = theta(p[P], pN[i], a[P], aNb[i]);
            const scalar thP = max(thP0, thetaMin_);
            thetaRawMin = min(thetaRawMin, min(thP0, 1 - thP0));
            if (thP0 < thetaMin_) ++nClipped_;
            ++nCrossing_;
            const scalar aP = (side[P] ? k1 : k2)*mS[i]/(thP*mag(delta[i]));
            diag[P] += aP;
            src[P] += aP*Tsat;
        }
    }

    // explicit return of the heat drained but not converted last step [W]
    src += heatReturn_;

    Info<< "HardtWondra: ghost-fluid faces = "
        << returnReduce(nCrossing_, sumOp<label>())
        << ", clipped (theta < " << thetaMin_ << ") = "
        << returnReduce(nClipped_, sumOp<label>())
        << ", min raw theta = " << returnReduce(thetaRawMin, minOp<scalar>())
        << endl;

    return tM;
}


void Foam::thermalPhaseChangeModels::HardtWondra::correctAfterT()
{
    const scalar Tsat = T_sat_.value();
    const scalar L = h_lv_.value();
    const scalar k1 = k1_.value();
    const scalar k2 = k2_.value();
    const boolList side(cellSide());
    const volScalarField& pf = psi();
    const scalarField& p = pf.primitiveField();
    const scalarField& a = alpha1_.primitiveField();
    const scalarField& T = T_.primitiveField();
    const volScalarField& aIntF =
        mesh_.lookupObject<volScalarField>("aInterface");
    const scalarField& aInt = aIntF.primitiveField();
    const vectorField& C = mesh_.C();
    const scalarField& magSf = mesh_.magSf();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    scalarField mCell(mesh_.nCells(), Zero);   // [kg/s]
    scalar Qsum = 0;

    // Share of Q_f received by cell P
    auto shareP = [](bool ifP, bool ifN, bool sideP, scalar thP0)
    {
        if (ifP && ifN) return 1 - thP0;     // closer cell gets more
        if (ifP) return scalar(1);
        if (ifN) return scalar(0);
        return sideP ? scalar(1) : scalar(0); // neither: the phase-1 cell
    };

    forAll(nei, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        if (side[P] == side[N]) continue;

        const scalar thP0 = theta(p[P], p[N], a[P], a[N]);
        const scalar thP = max(thP0, thetaMin_);
        const scalar thN = max(1 - thP0, thetaMin_);
        const scalar d = mag(C[N] - C[P]);
        const scalar aP = (side[P] ? k1 : k2)*magSf[facei]/(thP*d);
        const scalar aN = (side[N] ? k1 : k2)*magSf[facei]/(thN*d);
        const scalar Q = aP*(T[P] - Tsat) + aN*(T[N] - Tsat);
        const scalar sP = shareP(aInt[P] > 0, aInt[N] > 0, side[P], thP0);
        mCell[P] += sP*Q/L;
        mCell[N] += (1 - sP)*Q/L;
        Qsum += Q;
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fp = mesh_.boundary()[patchi];
        if (!fp.coupled()) continue;

        const scalarField pN(pf.boundaryField()[patchi].patchNeighbourField());
        const scalarField aNb
        (
            alpha1_.boundaryField()[patchi].patchNeighbourField()
        );
        const scalarField TN(T_.boundaryField()[patchi].patchNeighbourField());
        const scalarField ifN
        (
            aIntF.boundaryField()[patchi].patchNeighbourField()
        );
        const vectorField delta(fp.delta());
        const scalarField& mS = fp.magSf();
        const labelUList& fc = fp.faceCells();
        forAll(fc, i)
        {
            const label P = fc[i];
            const bool sN = (mag(pN[i]) > psiTol_*h_[fc[i]]) ? (pN[i] > 0) : (aNb[i] >= 0.5 - 1e-8);
            if (side[P] == sN) continue;
            const scalar thP0 = theta(p[P], pN[i], a[P], aNb[i]);
            const scalar thP = max(thP0, thetaMin_);
            const scalar thN = max(1 - thP0, thetaMin_);
            const scalar dd = mag(delta[i]);
            const scalar aP = (side[P] ? k1 : k2)*mS[i]/(thP*dd);
            const scalar aN = (sN ? k1 : k2)*mS[i]/(thN*dd);
            const scalar Q = aP*(T[P] - Tsat) + aN*(TN[i] - Tsat);
            mCell[P] += shareP(aInt[P] > 0, ifN[i] > 0, side[P], thP0)*Q/L;
            // each processor adds only its own cell's share; Qsum counts
            // this side's Dirichlet heat only (the other side counts its own)
            Qsum += aP*(T[P] - Tsat);
        }
    }

    if (fluxMethod_ == "probe")
    {
        // Malan et al. (2021) JCP 426:109920, Sec. 3.2: heat flux from pure
        // cells (T - Tsat)/d_Gamma next to each interface cell; mixed-cell
        // temperatures (volume averages) are not used. The ghost-fluid
        // estimate above uses them with 1/theta and drains the sensible heat
        // of a liquid cell whose centre the interface approaches (V4 +20-40 %).
        const scalarField& Vc = mesh_.V();
        scalarField mProbe(mesh_.nCells(), Zero);
        label nFallback = 0;

        const reconstructionSchemes& surf =
            mesh_.lookupObject<reconstructionSchemes>("reconstructionScheme");
        const volVectorField& nS = surf.normal();

        // pure-cell flux samples: G = (T - Tsat)/max(|psi|, 0.05 h), valid
        // only in non-interface cells with alpha away from (0, 1)
        volScalarField G
        (
            IOobject("probeGHW", mesh_.time().timeName(), mesh_,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh_, dimensionedScalar(dimless, Zero), "zeroGradient"
        );
        volScalarField pure
        (
            IOobject("probePureHW", mesh_.time().timeName(), mesh_,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh_, dimensionedScalar(dimless, Zero), "zeroGradient"
        );
        auto isPure = [&](const label c)
        {
            return
                aInt[c] <= 0 && p[c] != 0
             && (a[c] < 1e-3 || a[c] > 1 - 1e-3);
        };
        label nSecond = 0, nFirst = 0;
        forAll(G, celli)
        {
            if (!isPure(celli)) continue;
            const scalar d1 = max(mag(p[celli]), 0.05*h_[celli]);
            // second probe: face neighbour of the same phase farthest from
            // the interface -> quadratic through (0, Tsat), (d1, T1), (d2, T2)
            label qq = -1;
            scalar d2 = d1;
            for (const label facei : mesh_.cells()[celli])
            {
                if (!mesh_.isInternalFace(facei)) continue;
                const label nb = (own[facei] == celli ? nei[facei] : own[facei]);
                if (!isPure(nb) || side[nb] != side[celli]) continue;
                if (mag(p[nb]) > d2)
                {
                    d2 = mag(p[nb]);
                    qq = nb;
                }
            }
            if (qq >= 0 && d2 - d1 > 0.3*h_[celli])
            {
                G[celli] =
                    ((T[celli] - Tsat)*sqr(d2) - (T[qq] - Tsat)*sqr(d1))
                   /(d1*d2*(d2 - d1));
                ++nSecond;
            }
            else
            {
                G[celli] = (T[celli] - Tsat)/d1;
                ++nFirst;
            }
            pure[celli] = side[celli] ? 1 : -1;
        }
        G.correctBoundaryConditions();
        pure.correctBoundaryConditions();

        const reconstructedDistanceFunction& RDF =
            mesh_.lookupObject<reconstructedDistanceFunction>("RDF");
        zoneDistribute& dist = zoneDistribute::New(mesh_);
        const boolList& zone = RDF.nextToInterface();
        dist.setUpCommforZone(zone, false);
        const labelListList& stencil = dist.getStencil();
        const Map<scalar> mapG = dist.getDatafromOtherProc(zone, G);
        const Map<scalar> mapPure = dist.getDatafromOtherProc(zone, pure);
        const Map<vector> mapC = dist.getDatafromOtherProc(zone, mesh_.C());

        const globalIndex& gi = dist.globalNumbering();
        heatReturn_ = 0;
        forAll(aInt, celli)
        {
            if (aInt[celli] <= 0) continue;
            const vector n = nS[celli]/max(mag(nS[celli]), VSMALL);
            scalar sum[2] = {0, 0}, wsum[2] = {0, 0};
            DynamicList<label> localProbe;
            DynamicList<scalar> localW;
            for (const label g : stencil[celli])
            {
                const scalar pr = dist.getValue(pure, mapPure, g);
                if (pr == 0) continue;
                const vector dx = dist.getValue(mesh_.C(), mapC, g) - mesh_.C()[celli];
                const scalar mdx = mag(dx);
                if (mdx < VSMALL) continue;
                const scalar w = mag(n & dx)/mdx/sqr(mdx);
                const label s1 = pr > 0 ? 0 : 1;
                sum[s1] += w*dist.getValue(G, mapG, g);
                wsum[s1] += w;
                // zoneDistribute numbers cells and boundary faces: only
                // local cell entries can receive heat
                if (gi.isLocal(g) && gi.toLocal(g) < mesh_.nCells())
                {
                    localProbe.append(gi.toLocal(g));
                    localW.append(w*(s1 == 0 ? k1 : k2));
                }
            }
            // phase side without a pure neighbour: keep the ghost-fluid
            // share of that cell (logged)
            if (wsum[0] <= 0 || wsum[1] <= 0)
            {
                ++nFallback;
                mProbe[celli] = mCell[celli]/Vc[celli];
                continue;
            }
            const scalar q = k1*sum[0]/wsum[0] + k2*sum[1]/wsum[1];   // W/m2
            mProbe[celli] = q/L*aInt[celli];

            // Heat drained by the ghost-fluid pin but not converted into
            // mass is returned to the probe cells in the next step, so the
            // temperature equation loses exactly L mdot (energy consistency)
            const scalar D = L*(mCell[celli] - mProbe[celli]*Vc[celli]);   // W
            scalar wl = 0;
            forAll(localW, k) wl += localW[k];
            if (wl > 0)
            {
                forAll(localProbe, k)
                {
                    heatReturn_[localProbe[k]] += D*localW[k]/wl;
                }
            }
            else
            {
                heatReturn_[celli] += D;
            }
        }
        mdot_.primitiveFieldRef() = mProbe;
        const scalar mGF = gSum(mCell);
        const scalar mPr = gSum(mProbe*Vc);
        Info<< "HardtWondra: probe flux: mdotTotal = " << mPr
            << " kg/s (ghost-fluid estimate " << mGF << "), interface cells"
            << " without pure neighbours on one side = "
            << returnReduce(nFallback, sumOp<label>())
            << ", probes 2nd/1st order = " << returnReduce(nSecond, sumOp<label>())
            << "/" << returnReduce(nFirst, sumOp<label>())
            << ", heat returned (pin heat - L mdot) = " << L*(mGF - mPr)
            << " W" << endl;
    }
    else
    {
        mdot_.primitiveFieldRef() = mCell/mesh_.V();
        heatReturn_ = 0;
    }
    mdot_.correctBoundaryConditions();

    reduce(Qsum, sumOp<scalar>());
    const scalar mTot = gSum(mdot_.primitiveField()*mesh_.V());
    Info<< "HardtWondra: interface heat sum(Q_f) = " << Qsum
        << " W, L sum(mdot) = " << L*mTot << " W, mdotTotal = " << mTot
        << " kg/s" << endl;

    calcPCV();
}


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::alpha1Gen() const
{
    if (!sw_alpha1Gen_)
    {
        return volScalarField::New
        (
            "alpha1Gen", mesh_, dimensionedScalar(dimless/dimTime, Zero)
        );
    }
    return volScalarField::New("alpha1Gen", -mdot_/mixture_.rho1());
}


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::PCV() const
{
    if (!sw_PCV_)
    {
        return volScalarField::New
        (
            "PCV0", mesh_, dimensionedScalar(dimless/dimTime, Zero)
        );
    }
    return PCV_;
}


Foam::scalar
Foam::thermalPhaseChangeModels::HardtWondra::maxInterfaceSpeed() const
{
    const scalarField& aInt =
        mesh_.lookupObject<volScalarField>("aInterface").primitiveField();
    const scalar rho1 = mixture_.rho1().value();
    scalar u = 0;
    forAll(mdot_, celli)
    {
        const scalar m = mag(mdot_[celli]);
        if (m > 0)
        {
            // mdot'' = mdot'''/a_i; cells without a facet: a = 1/h
            const scalar ai = aInt[celli] > 0 ? aInt[celli] : 1/h_[celli];
            u = max(u, m/(rho1*ai));
        }
    }
    return returnReduce(u, maxOp<scalar>());
}


bool Foam::thermalPhaseChangeModels::HardtWondra::read
(
    const dictionary& dict
)
{
    thermalPhaseChangeModel::read(dict);

    k1_.value() = dict.get<scalar>("k1");
    k2_.value() = dict.get<scalar>("k2");
    thetaMin_ = dict.getOrDefault<scalar>("thetaMin", 0.05);
    lambdaCells_ = dict.getOrDefault<scalar>("lambdaCells", 1);
    fluxMethod_ = dict.getOrDefault<word>("fluxMethod", "probe");
    if (fluxMethod_ != "probe" && fluxMethod_ != "ghostFluid")
    {
        FatalIOErrorInFunction(dict)
            << "fluxMethod " << fluxMethod_ << " unknown; valid: probe ghostFluid"
            << exit(FatalIOError);
    }

    for
    (
        const word key
      : {
            "k_liq", "k_vap", "coldPhaseIsHighAlpha1", "interfaceWidthCells",
            "betaThermal", "kinFloorCells", "liquidBiasCoeff", "harmonicBias",
            "RelaxFac", "useEnthalpyCorrection", "nExtrapIter", "extrapCFL",
            "Qcorr", "nSmoothIter", "alphaSmoothWidth", "lambdaSmearCells",
            "lambdaAlphaCells", "mdotMax", "AiFloorAbs", "maxGradT",
            "alphaPureLiquid", "alphaPureVapor", "useTangentialSmoothing"
        }
    )
    {
        if (dict.found(key))
        {
            WarningInFunction
                << "phaseChangeProperties entry '" << key << "' is not used"
                << " by model HardtWondra (legacy model: HardtWondraLegacy)"
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //
