// Phase-change advection velocity of the interface (velocity mode).
//
// Gada & Sharma (2009), Eq. 13; Shaikh et al. (2016), Eq. 8: the interface is
// advected with the total velocity u_t = u + u_PC,
//
//     d(alpha1)/dt + u_t . grad(alpha1) = 0,
//     u_PC = C mdot'' n,     n = grad(psi)/|grad(psi)|  (psi = psiRDF, > 0 in
//                                                        phase 1, n into 1)
//
// mdot'' [kg/m2/s] > 0 converts phase 1 into phase 2. In the flux form of
// isoAdvection this reads
//
//     d(alpha1)/dt + div(alpha1 phiT) = alpha1 div(phiPC),   phiT = phi + phiPC
//
// so Su = alpha1 div(phiPC) [1/s] is passed to advect() (Sp = 0). In a pure
// cell alpha_f = alpha1 and the two cancel exactly; the volume change appears
// only where alpha_f differs from the cell value, i.e. at the interface. The
// bulk dilatation alpha1 div(u) stays out (it sits in pEqn through PCV).
//
// C = 1/rho1: the HardtWondra dilatation is shifted off the interface, so the
// bulk velocity at the interface is the liquid velocity and the interface
// moves at mdot''/rho1 relative to it (Gada & Sharma's 0.5(1/rho1 + 1/rho2)
// assumes the dilatation sits at the interface).
//
// mdot''' = -rho1 alpha1Gen is converted to a flux per interface area with the
// PLIC area density (cells without a facet: 1/h). Written for sharp models
// (mdot''' confined to interface cells). On a face the area flux is the mean of
// the neighbouring cells that carry one.

const dimensionedScalar CPC("CPC", dimVolume/dimMass, 1/rho1.value());

tmp<volScalarField> tAlpha1Gen(phaseChangePtr->alpha1Gen());

// mdot'' on the cells that carry the source, and a 0/1 marker of those cells
volScalarField mdotArea
(
    IOobject("mdotArea", mesh.time().timeName(), mesh),
    mesh,
    dimensionedScalar(dimMass/dimArea/dimTime, Zero)
);
volScalarField mdotMask
(
    IOobject("mdotMask", mesh.time().timeName(), mesh),
    mesh,
    dimensionedScalar(dimless, Zero)
);

{
    const scalarField& S = tAlpha1Gen().primitiveField();
    const scalarField& aInt = interface.aInterface().primitiveField();
    const scalarField& h = interface.cellSize();
    scalarField& mA = mdotArea.primitiveFieldRef();
    scalarField& mM = mdotMask.primitiveFieldRef();

    forAll(S, celli)
    {
        if (S[celli] == 0) continue;
        const scalar ai = (aInt[celli] > 0 ? aInt[celli] : 1/h[celli]);
        mA[celli] = -rho1.value()*S[celli]/ai;
        mM[celli] = 1;
    }
}
mdotArea.correctBoundaryConditions();
mdotMask.correctBoundaryConditions();

// RDF-based interface normal (into phase 1)
volVectorField nRDF(fvc::grad(interface.psi()));
nRDF /= max(mag(nRDF), dimensionedScalar(dimless, SMALL));

// Face flux of the phase-change velocity [m3/s]
{
    const surfaceScalarField mdotAreaf
    (
        fvc::interpolate(mdotArea)
       /max(fvc::interpolate(mdotMask), dimensionedScalar(dimless, SMALL))
    );

    surfaceVectorField nf(fvc::interpolate(nRDF));
    nf /= max(mag(nf), dimensionedScalar(dimless, SMALL));

    phiPC = CPC*mdotAreaf*(nf & mesh.Sf());

    surfaceScalarField::Boundary& phiPCb = phiPC.boundaryFieldRef();
    forAll(phiPCb, patchi)
    {
        if (!phiPCb[patchi].coupled())
        {
            phiPCb[patchi] = 0;
        }
    }
}

// Cell velocity (isoAdvection times the face crossings with it)
const volVectorField UPC(CPC*mdotArea*nRDF);

Sp.field() = Zero;
{
    const scalarField divPhiPC(fvc::surfaceIntegrate(phiPC)().primitiveField());
    const scalarField& a = alpha1.primitiveField();
    scalarField& su = Su.field();
    forAll(su, celli)
    {
        su[celli] = min(max(a[celli], scalar(0)), scalar(1))*divPhiPC[celli];
    }

    Info<< "alphaPCVelocity: max|u_PC| = " << gMax(mag(UPC.primitiveField()))
        << " m/s, requested " << gSum(tAlpha1Gen().primitiveField()*mesh.V())
        << " m3/s" << endl;
}
