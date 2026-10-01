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

#include "interfaceCurvatureITF.H"
#include "zoneDistribute.H"
#include "leastSquaresGrad.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interfaceCurvatureITF::calcCellSize()
{
    // V / (largest own face area), skipping wedge/empty faces: the in-plane
    // spacing on 2D planar and axisymmetric wedge meshes, h on a cube
    h_.setSize(mesh_.nCells());
    const cellList& cells = mesh_.cells();
    const scalarField& V = mesh_.V();
    const scalarField& magSf = mesh_.magSf();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(cells, celli)
    {
        scalar maxA = SMALL;
        for (const label facei : cells[celli])
        {
            scalar a = 0;
            if (mesh_.isInternalFace(facei))
            {
                a = magSf[facei];
            }
            else
            {
                const label patchi = pbm.whichPatch(facei);
                const polyPatch& pp = pbm[patchi];
                if (isA<wedgePolyPatch>(pp) || isA<emptyPolyPatch>(pp))
                {
                    continue;
                }
                a = mesh_.magSf().boundaryField()[patchi][facei - pp.start()];
            }
            maxA = max(maxA, a);
        }
        h_[celli] = V[celli]/maxA;
    }
}


Foam::label Foam::interfaceCurvatureITF::neighbourAlong
(
    const label celli,
    const vector& e
) const
{
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    for (const label facei : mesh_.cells()[celli])
    {
        if (!mesh_.isInternalFace(facei)) continue;
        vector n = mesh_.faceAreas()[facei];
        const scalar mn = mag(n);
        if (mn < VSMALL) continue;
        n /= mn;
        if (own[facei] != celli) n = -n;
        if ((n & e) > 0.99)
        {
            return (own[facei] == celli ? nei[facei] : own[facei]);
        }
    }
    return -1;
}


Foam::label Foam::interfaceCurvatureITF::heightFunctionKappa
(
    const reconstructionSchemes& surf
)
{
    if (!wedge_ && mesh_.nSolutionD() != 2)
    {
        return 0;   // 3D: RDF everywhere
    }

    const scalarField& a = alpha1_.primitiveField();
    const vectorField& C = mesh_.C();
    const pointField& pts = mesh_.points();
    const scalar tol = 1e-6;

    // Extent of a cell along direction e (from its points)
    auto extent = [&](const label celli, const vector& e, scalar& lo, scalar& hi)
    {
        lo = GREAT; hi = -GREAT;
        for (const label pi : mesh_.cellPoints()[celli])
        {
            scalar x = pts[pi] & e;
            if (wedge_ && e == e1_)
            {
                // radial coordinate, independent of the wedge angle
                x = mag(pts[pi] - (pts[pi] & e2_)*e2_);
            }
            lo = min(lo, x);
            hi = max(hi, x);
        }
    };

    // Height function along d, tangential t, both positive unit axes.
    auto hf = [&](const label celli, const vector& d, const vector& t, scalar& k)
    {
        const bool radialCols = wedge_ && d == e1_;

        // Base cells of the three columns (j = 0, 1, 2 <-> -t, 0, +t)
        label base[3] = {neighbourAlong(celli, -t), celli, neighbourAlong(celli, t)};
        if (base[2] < 0) return false;
        bool mirror = false;
        if (base[0] < 0)
        {
            // Axis: the -r neighbour of an axial column is its mirror image
            if (wedge_ && t == e1_ && (C[celli] & e1_) < 2*h_[celli])
            {
                mirror = true;
                base[0] = celli;
            }
            else
            {
                return false;
            }
        }

        // cells[j][3 + m], m = -3..3
        label cells[3][7];
        for (label j = 0; j < 3; ++j)
        {
            for (label m = 0; m < 7; ++m) cells[j][m] = -1;
            cells[j][3] = base[j];
        }

        auto layerAvg = [&](const label m)
        {
            return (a[cells[0][m]] + a[cells[1][m]] + a[cells[2][m]])/3;
        };
        auto full = [&](const scalar x) { return x < tol || x > 1 - tol; };

        label top = -1, bot = -1;
        for (label sgnDir : {1, -1})
        {
            label end = -1;
            for (label k = 1; k <= 3; ++k)
            {
                const label m = 3 + sgnDir*k;
                for (label j = 0; j < 3; ++j)
                {
                    cells[j][m] =
                        neighbourAlong(cells[j][m - sgnDir], scalar(sgnDir)*d);
                    if (cells[j][m] < 0) return false;
                }
                if (full(layerAvg(m)))
                {
                    end = m;
                    break;
                }
            }
            if (end < 0) return false;
            (sgnDir > 0 ? top : bot) = end;
        }

        const scalar aTop = layerAvg(top);
        const scalar aBot = layerAvg(bot);
        if ((aTop > 0.5) == (aBot > 0.5)) return false;
        const bool phase1Up = aTop > 0.5;
        const scalar sgn = phase1Up ? 1 : -1;

        // Interface position of each column
        scalar f[3], x[3];
        for (label j = 0; j < 3; ++j)
        {
            scalar lo, hi;
            extent(cells[j][bot], d, lo, hi);
            const scalar start = lo;
            scalar sum = 0;
            for (label m = bot; m <= top; ++m)
            {
                const label c = cells[j][m];
                extent(c, d, lo, hi);
                const scalar below = phase1Up ? 1 - a[c] : a[c];
                // radial columns: volume-weighted (r^2) heights
                sum += below*(radialCols ? (sqr(hi) - sqr(lo)) : (hi - lo));
            }
            f[j] = radialCols ? Foam::sqrt(sqr(start) + sum) : start + sum;
            if (wedge_ && t == e1_)
            {
                // distance from the wedge axis (the mesh is a 3D wedge)
                const vector& c = C[base[j]];
                x[j] = mag(c - (c & e2_)*e2_);
            }
            else
            {
                x[j] = C[base[j]] & t;
            }
        }
        if (mirror)
        {
            x[0] = -x[1];
            f[0] = f[1];
        }

        const scalar dm = x[1] - x[0];
        const scalar dp = x[2] - x[1];
        if (dm < SMALL || dp < SMALL) return false;
        const scalar f1 = ((f[2] - f[1])/dp*dm + (f[1] - f[0])/dm*dp)/(dm + dp);
        const scalar f2 = 2*((f[2] - f[1])/dp - (f[1] - f[0])/dm)/(dm + dp);
        const scalar q = Foam::sqrt(1 + sqr(f1));

        k = sgn*f2/pow3(q);
        if (wedge_)
        {
            if (radialCols)
            {
                // g(y) = r of the interface; hoop = -n_r/r
                k += -sgn/(f[1]*q);
            }
            else
            {
                // f(r) = y of the interface at the column radius x[1]
                k += sgn*f1/(x[1]*q);
            }
        }
        return true;
    };

    label nValid = 0;
    forAll(interfaceCell_, celli)
    {
        if (!interfaceCell_[celli]) continue;
        const vector n = surf.normal()[celli];
        const bool first1 = mag(n & e1_) >= mag(n & e2_);
        const vector dA = first1 ? e1_ : e2_;
        const vector dB = first1 ? e2_ : e1_;
        scalar k = 0;
        if (hf(celli, dA, dB, k) || hf(celli, dB, dA, k))
        {
            kappa_[celli] = k;
            ++nValid;
        }
    }
    return nValid;
}



void Foam::interfaceCurvatureITF::markInterfaceCells
(
    const reconstructionSchemes& surf
)
{
    const volVectorField& normal = surf.normal();
    const scalarField& a = alpha1_.primitiveField();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    // Cells with a face neighbour on the other side of alpha = 0.5
    boolList touchesOther(mesh_.nCells(), false);
    forAll(nei, facei)
    {
        if ((a[own[facei]] >= 0.5) != (a[nei[facei]] >= 0.5))
        {
            touchesOther[own[facei]] = true;
            touchesOther[nei[facei]] = true;
        }
    }
    forAll(alpha1_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pa = alpha1_.boundaryField()[patchi];
        if (pa.coupled())
        {
            const scalarField an(pa.patchNeighbourField());
            const labelUList& fc = pa.patch().faceCells();
            forAll(fc, i)
            {
                if ((a[fc[i]] >= 0.5) != (an[i] >= 0.5))
                {
                    touchesOther[fc[i]] = true;
                }
            }
        }
    }

    nInterface_ = 0;
    nWisp_ = 0;
    aInterface_ = dimensionedScalar(aInterface_.dimensions(), Zero);
    forAll(interfaceCell_, celli)
    {
        const bool hasFacet = mag(normal[celli]) > 0;
        interfaceCell_[celli] = hasFacet && touchesOther[celli];
        if (interfaceCell_[celli])
        {
            ++nInterface_;
            aInterface_[celli] = mag(normal[celli])/mesh_.V()[celli];
        }
        else if (hasFacet)
        {
            ++nWisp_;
        }
    }
    aInterface_.correctBoundaryConditions();
    reduce(nInterface_, sumOp<label>());
    reduce(nWisp_, sumOp<label>());
}


void Foam::interfaceCurvatureITF::calcPsi
(
    const reconstructionSchemes& surf,
    reconstructedDistanceFunction& RDF
)
{
    const label nRings = 3;
    RDF.markCellsNearSurf(interfaceCell_, nRings);
    const boolList& zone = RDF.nextToInterface();
    level_.primitiveFieldRef() = RDF.cellDistLevel().primitiveField();
    level_.correctBoundaryConditions();

    zoneDistribute& dist = zoneDistribute::New(mesh_);
    dist.setUpCommforZone(zone, true);
    const labelListList& stencil = dist.getStencil();
    const Map<scalar> mapLevel = dist.getDatafromOtherProc(zone, level_);

    // Facet carried by each band cell (centre, unit normal)
    volVectorField fc
    (
        IOobject("facetCentreITF", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_, dimensionedVector(dimLength, Zero), "zeroGradient"
    );
    volVectorField fn
    (
        IOobject("facetNormalITF", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_, dimensionedVector(dimless, Zero), "zeroGradient"
    );

    const volVectorField& centre = surf.centre();
    const volVectorField& normal = surf.normal();
    const vectorField& C = mesh_.C();
    scalarField dRaw(mesh_.nCells(), Zero);

    forAll(interfaceCell_, celli)
    {
        if (interfaceCell_[celli])
        {
            fc[celli] = centre[celli];
            fn[celli] = normal[celli]/mag(normal[celli]);
            dRaw[celli] = (C[celli] - fc[celli]) & fn[celli];
        }
    }

    for (label ring = 1; ring <= nRings; ++ring)
    {
        fc.correctBoundaryConditions();
        fn.correctBoundaryConditions();
        const Map<vector> mapC = dist.getDatafromOtherProc(zone, fc);
        const Map<vector> mapN = dist.getDatafromOtherProc(zone, fn);

        // Collect first, assign after: the ring-(L) cells must only see
        // facets of rings < L
        DynamicList<label> cellsL;
        DynamicList<vector> newC, newN;

        forAll(zone, celli)
        {
            if (!zone[celli] || label(level_[celli]) != ring) continue;

            const point& p = C[celli];
            scalar dmin = GREAT;
            forAll(stencil[celli], k)
            {
                const label g = stencil[celli][k];
                const label lg = label(dist.getValue(level_, mapLevel, g));
                if (lg < 0 || lg >= ring) continue;
                const vector n = dist.getValue(fn, mapN, g);
                if (mag(n) < SMALL) continue;
                dmin = min(dmin, mag(dist.getValue(fc, mapC, g) - p));
            }
            if (dmin >= GREAT) continue;

            scalar sumD = 0, sumW = 0, best = GREAT;
            vector bestC(Zero), bestN(Zero);
            forAll(stencil[celli], k)
            {
                const label g = stencil[celli][k];
                const label lg = label(dist.getValue(level_, mapLevel, g));
                if (lg < 0 || lg >= ring) continue;
                const vector n = dist.getValue(fn, mapN, g);
                if (mag(n) < SMALL) continue;
                const vector c = dist.getValue(fc, mapC, g);
                vector dc(p - c);
                const scalar dist2c = mag(dc);
                if (dist2c > dmin + h_[celli]) continue;
                // reconstructedDistanceFunction weighting
                const scalar w =
                    dist2c > SMALL ? sqr((dc/dist2c) & n) : scalar(1);
                sumD += w*(dc & n);
                sumW += w;
                if (dist2c < best)
                {
                    best = dist2c;
                    bestC = c;
                    bestN = n;
                }
            }
            if (sumW > SMALL)
            {
                dRaw[celli] = sumD/sumW;
                cellsL.append(celli);
                newC.append(bestC);
                newN.append(bestN);
            }
        }
        forAll(cellsL, i)
        {
            fc[cellsL[i]] = newC[i];
            fn[cellsL[i]] = newN[i];
        }
    }

    // Orientation: psi > 0 in phase 1. Vote on the band cells outside the
    // interface (pure cells), log the agreement.
    scalar vote = 0;
    label nAgree = 0, nVote = 0;
    forAll(dRaw, celli)
    {
        if (level_[celli] >= 1 && dRaw[celli] != 0)
        {
            const scalar s = alpha1_[celli] - 0.5;
            if (mag(s) > 0.49)
            {
                vote += dRaw[celli]*s;
                ++nVote;
                if (dRaw[celli]*s > 0) ++nAgree;
            }
        }
    }
    reduce(vote, sumOp<scalar>());
    reduce(nVote, sumOp<label>());
    reduce(nAgree, sumOp<label>());
    const scalar sgn = (vote < 0 ? -1 : 1);

    psi_.primitiveFieldRef() = sgn*dRaw;
    psi_.correctBoundaryConditions();

    Info<< "interfaceCurvatureITF: interface cells = " << nInterface_
        << ", wisps ignored = " << nWisp_
        << ", psi orientation " << sgn << " (agreement "
        << (nVote ? scalar(sgn > 0 ? nAgree : nVote - nAgree)/nVote : 1)
        << " of " << nVote << " pure band cells)" << endl;
}


void Foam::interfaceCurvatureITF::calcKappa
(
    const reconstructionSchemes& surf,
    reconstructedDistanceFunction& RDF
)
{
    const boolList& zone = RDF.nextToInterface();
    const scalarField& lev = level_.primitiveField();

    const volVectorField gradPsi
    (
        fv::leastSquaresGrad<scalar>(mesh_).calcGrad(psi_, "grad(psiRDF)")
    );

    volVectorField nHat
    (
        IOobject("nHatITF", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_, dimensionedVector(dimless, Zero), "zeroGradient"
    );
    forAll(nHat, celli)
    {
        if (lev[celli] >= 0 && lev[celli] <= 2)
        {
            const vector& g = gradPsi[celli];
            const scalar mg = mag(g);
            if (mg > SMALL) nHat[celli] = g/mg;
        }
    }
    nHat.correctBoundaryConditions();

    volScalarField kCell
    (
        IOobject("kappaCellITF", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_, dimensionedScalar(dimless/dimLength, Zero), "zeroGradient"
    );
    {
        const volScalarField divN
        (
            fvc::div(fvc::interpolate(nHat) & mesh_.Sf())
        );
        forAll(kCell, celli)
        {
            if (lev[celli] >= 0 && lev[celli] <= 1)
            {
                kCell[celli] = -divN[celli];
            }
        }
    }
    kCell.correctBoundaryConditions();

    zoneDistribute& dist = zoneDistribute::New(mesh_);
    const labelListList& stencil = dist.getStencil();
    const Map<scalar> mapK = dist.getDatafromOtherProc(zone, kCell);
    const Map<scalar> mapPsi = dist.getDatafromOtherProc(zone, psi_);
    const Map<scalar> mapLev = dist.getDatafromOtherProc(zone, level_);
    const globalIndex& gi = dist.globalNumbering();

    kappa_ = dimensionedScalar(kappa_.dimensions(), Zero);

    // Interface cells: shifted, distance-weighted stencil mean
    forAll(interfaceCell_, celli)
    {
        if (!interfaceCell_[celli]) continue;

        const scalar hc = h_[celli];
        scalar sumK = 0, sumW = 0;
        bool selfSeen = false;
        const label gSelf = gi.toGlobal(celli);

        auto add = [&](const scalar kj, const scalar pj)
        {
            const scalar denom = 1 + kj*pj/nCurvDir_;
            if (denom < 0.5) return;
            const scalar w = 1/(mag(pj) + 0.1*hc);
            sumK += w*kj/denom;
            sumW += w;
        };

        for (const label g : stencil[celli])
        {
            if (g == gSelf) selfSeen = true;
            const scalar lj = dist.getValue(level_, mapLev, g);
            if (lj < 0 || lj > 1) continue;
            add
            (
                dist.getValue(kCell, mapK, g),
                dist.getValue(psi_, mapPsi, g)
            );
        }
        if (!selfSeen)
        {
            add(kCell[celli], psi_[celli]);
        }
        if (sumW > SMALL)
        {
            kappa_[celli] = sumK/sumW;
        }
    }
    if (model_ == "heightFunction")
    {
        const label nHF = returnReduce
        (
            heightFunctionKappa(surf),
            sumOp<label>()
        );
        Info<< "interfaceCurvatureITF: height function on " << nHF
            << " of " << nInterface_ << " interface cells (RDF elsewhere)"
            << endl;
    }
    kappa_.correctBoundaryConditions();

    // Ring-1 cells: curvature of the nearest interface cell (facet centroid)
    {
        volScalarField isIf
        (
            IOobject("isInterfaceITF", mesh_.time().timeName(), mesh_,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh_, dimensionedScalar(dimless, Zero), "zeroGradient"
        );
        forAll(isIf, celli)
        {
            isIf[celli] = interfaceCell_[celli] ? 1 : 0;
        }
        isIf.correctBoundaryConditions();
        const Map<scalar> mapIf = dist.getDatafromOtherProc(zone, isIf);
        const Map<scalar> mapKI = dist.getDatafromOtherProc(zone, kappa_);
        const Map<vector> mapCC =
            dist.getDatafromOtherProc(zone, mesh_.C());

        scalarField kNew(kappa_.primitiveField());
        forAll(kNew, celli)
        {
            if (interfaceCell_[celli] || lev[celli] != 1) continue;
            scalar best = GREAT;
            for (const label g : stencil[celli])
            {
                if (dist.getValue(isIf, mapIf, g) < 0.5) continue;
                const scalar d =
                    mag(dist.getValue(mesh_.C(), mapCC, g) - mesh_.C()[celli]);
                if (d < best)
                {
                    best = d;
                    kNew[celli] = dist.getValue(kappa_, mapKI, g);
                }
            }
        }
        kappa_.primitiveFieldRef() = kNew;
        kappa_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCurvatureITF::interfaceCurvatureITF
(
    const immiscibleIncompressibleTwoPhaseMixture& mixture,
    const volScalarField& alpha1
)
:
    mesh_(alpha1.mesh()),
    mixture_(mixture),
    alpha1_(alpha1),
    model_(mixture.getOrDefault<word>("curvatureModel", "gradAlpha")),
    sigmaPtr_(surfaceTensionModel::New(mixture, alpha1.mesh())),
    h_(),
    nCurvDir_(2),
    wedge_(false),
    e1_(Zero),
    e2_(Zero),
    hCapMin_(GREAT),
    level_
    (
        IOobject("bandLevelITF", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_, dimensionedScalar(dimless, -1), "zeroGradient"
    ),
    psi_
    (
        IOobject("psiRDF", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh_, dimensionedScalar(dimLength, Zero), "zeroGradient"
    ),
    aInterface_
    (
        IOobject("aInterface", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_, dimensionedScalar(dimless/dimLength, Zero), "zeroGradient"
    ),
    kappa_
    (
        IOobject("kappaI", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh_, dimensionedScalar(dimless/dimLength, Zero), "zeroGradient"
    ),
    interfaceCell_(mesh_.nCells(), false),
    nInterface_(0),
    nWisp_(0),
    constantKappa_(0)
{
    if
    (
        model_ != "gradAlpha" && model_ != "RDF"
     && model_ != "heightFunction" && model_ != "constant"
    )
    {
        FatalIOErrorInFunction(mixture)
            << "curvatureModel " << model_
            << " unknown; valid: gradAlpha RDF heightFunction constant"
            << exit(FatalIOError);
    }
    if (model_ == "constant")
    {
        // Balanced-force verification only: prescribed exact curvature
        constantKappa_ = mixture.get<scalar>("constantCurvature");
    }

    calcCellSize();

    // An OpenFOAM axisymmetric case is a 3D wedge (nSolutionD = 3):
    // take the in-plane frame from the wedge geometry, not solutionD.
    for (const polyPatch& pp : mesh_.boundaryMesh())
    {
        if (isA<wedgePolyPatch>(pp))
        {
            const wedgePolyPatch& wp = refCast<const wedgePolyPatch>(pp);
            wedge_ = true;
            e2_ = normalised(wp.axis());                    // axis
            e1_ = normalised(e2_ ^ wp.centreNormal());      // radial
            if ((e1_ & (mesh_.C()[0] - (mesh_.C()[0] & e2_)*e2_)) < 0)
            {
                e1_ = -e1_;
            }
            break;
        }
    }
    nCurvDir_ = (!wedge_ && mesh_.nSolutionD() == 2) ? 1 : 2;

    if (!wedge_ && mesh_.nSolutionD() == 2)
    {
        DynamicList<vector> axes;
        const Vector<label>& sd = mesh_.solutionD();
        for (direction c = 0; c < vector::nComponents; ++c)
        {
            if (sd[c] == 1)
            {
                vector e(Zero);
                e[c] = 1;
                axes.append(e);
            }
        }
        e1_ = axes[0];
        e2_ = axes[1];
    }

    // Capillary length scale. On a wedge the triangular cells at the axis
    // have V/max(face area) = h/2 although their radial and axial extents
    // are h: use the point extents along the in-plane axes there.
    hCapMin_ = gMin(h_);
    if (wedge_)
    {
        const pointField& pts = mesh_.points();
        scalar hMin = GREAT;
        forAll(h_, celli)
        {
            scalar lo1 = GREAT, hi1 = -GREAT, lo2 = GREAT, hi2 = -GREAT;
            for (const label pi : mesh_.cellPoints()[celli])
            {
                const scalar d1 = pts[pi] & e1_;
                const scalar d2 = pts[pi] & e2_;
                lo1 = min(lo1, d1);
                hi1 = max(hi1, d1);
                lo2 = min(lo2, d2);
                hi2 = max(hi2, d2);
            }
            hMin = min(hMin, min(hi1 - lo1, hi2 - lo2));
        }
        hCapMin_ = returnReduce(hMin, minOp<scalar>());
    }

    Info<< "interfaceCurvatureITF: curvatureModel " << model_
        << ", curvature directions " << nCurvDir_
        << ", nSolutionD " << mesh_.nSolutionD()
        << ", wedge " << wedge_ << ", in-plane axes " << e1_ << " " << e2_
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interfaceCurvatureITF::correct
(
    const reconstructionSchemes& surf,
    reconstructedDistanceFunction& RDF
)
{
    markInterfaceCells(surf);
    calcPsi(surf, RDF);

    if (model_ == "RDF" || model_ == "heightFunction")
    {
        calcKappa(surf, RDF);
    }
    else if (model_ == "constant")
    {
        forAll(kappa_, celli)
        {
            kappa_[celli] =
                (level_[celli] >= 0 && level_[celli] <= 1) ? constantKappa_ : 0;
        }
        kappa_.correctBoundaryConditions();
    }
    else
    {
        // Stock interfaceProperties curvature (calculateK): -div(nHatf)
        kappa_ = -fvc::div(mixture_.nHatf());
    }

    // Diagnostics over interface cells
    scalar kMin = GREAT, kMax = -GREAT, kSum = 0;
    label n = 0;
    forAll(interfaceCell_, celli)
    {
        if (interfaceCell_[celli])
        {
            kMin = min(kMin, kappa_[celli]);
            kMax = max(kMax, kappa_[celli]);
            kSum += kappa_[celli];
            ++n;
        }
    }
    reduce(kMin, minOp<scalar>());
    reduce(kMax, maxOp<scalar>());
    reduce(kSum, sumOp<scalar>());
    reduce(n, sumOp<label>());
    Info<< "interfaceCurvatureITF: kappa over interface cells"
        << " min/mean/max = " << kMin << " / " << kSum/max(n, 1)
        << " / " << kMax << endl;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceCurvatureITF::surfaceTensionForce() const
{
    if (model_ == "gradAlpha")
    {
        return mixture_.surfaceTensionForce();
    }

    surfaceScalarField kf
    (
        IOobject("kappafITF", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_, dimensionedScalar(dimless/dimLength, Zero)
    );

    const scalarField& a = alpha1_.primitiveField();
    const scalarField& k = kappa_.primitiveField();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    auto wt = [](const scalar x) { return x*(1 - x) + 1e-6; };

    forAll(nei, facei)
    {
        const scalar wP = wt(a[own[facei]]);
        const scalar wN = wt(a[nei[facei]]);
        kf[facei] = (wP*k[own[facei]] + wN*k[nei[facei]])/(wP + wN);
    }
    forAll(kf.boundaryField(), patchi)
    {
        fvsPatchScalarField& pkf = kf.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pk = kappa_.boundaryField()[patchi];
        const fvPatchScalarField& pa = alpha1_.boundaryField()[patchi];
        if (pk.coupled())
        {
            const scalarField kI(pk.patchInternalField());
            const scalarField kN(pk.patchNeighbourField());
            const scalarField aI(pa.patchInternalField());
            const scalarField aN(pa.patchNeighbourField());
            forAll(pkf, i)
            {
                const scalar wP = wt(aI[i]);
                const scalar wN = wt(aN[i]);
                pkf[i] = (wP*kI[i] + wN*kN[i])/(wP + wN);
            }
        }
        else
        {
            pkf = pk;
        }
    }

    return fvc::interpolate(sigmaPtr_->sigma())*kf*fvc::snGrad(alpha1_);
}


Foam::scalar Foam::interfaceCurvatureITF::capillaryDeltaT() const
{
    const scalar sigmaMax = gMax(sigmaPtr_->sigma()().primitiveField());
    if (sigmaMax < SMALL)
    {
        return GREAT;
    }
    const scalar rhoAvg =
        0.5*(mixture_.rho1().value() + mixture_.rho2().value());
    return Foam::sqrt
    (
        rhoAvg*pow3(hCapMin_)/(constant::mathematical::twoPi*sigmaMax)
    );
}


// ************************************************************************* //
