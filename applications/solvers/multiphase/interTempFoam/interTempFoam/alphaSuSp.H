// Phase-change source of the alpha equation (phaseChangeAdvection source),
// passed to isoAdvector as an explicit Su [1/s] (Sp = 0). Su and Sp are
// declared in alphaEqn.H.
//
// S = phaseChangePtr->alpha1Gen() [1/s]: < 0 consumes phase 1, > 0 produces it.
//
// Budget: in one step a cell can give at most the phase-1 volume it holds
// (alpha V) and take at most the room it has ((1 - alpha) V). The implicit
// sink Sp = S/max(alpha, floor) used before delivered only 45-60% of the
// requested sink on V2 (dt = 36 s: ~0.6 cell volumes per step), so the
// latent heat removed by the temperature equation melted less phase 1 than
// it paid for. Here the excess of a cell is passed, over up to nPass sweeps,
// to the face neighbour best aligned with the interface normal (grad(alpha1)
// for sinks, -grad(alpha1) for sources; pure cells inherit the direction of
// their donor), so a planar front stays planar; what cannot be placed
// (processor boundaries, no neighbour) is logged.
//
// No divU term: isoAdvector is conservative flux-form advection and the
// dilatation (PCV) sits outside the receiving cells (see git history of
// cbf786e for the spurious growth an alpha*divU term caused).

tmp<volScalarField> tAlpha1Gen(phaseChangePtr->alpha1Gen());

{
    const scalarField& S = tAlpha1Gen().primitiveField();
    const scalarField& a = alpha1.primitiveField();
    const scalarField& V = mesh.V();
    const scalar dt = runTime.deltaTValue();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const cellList& cells = mesh.cells();
    const label nPass = 5;

    // Requested phase-1 volume change per cell over the step [m3]
    scalarField dVol(S*V*dt);

    // Push direction: unit grad(alpha1) (towards phase 1)
    vectorField dir(fvc::grad(alpha1)().primitiveField());
    forAll(dir, celli)
    {
        const scalar m = mag(dir[celli]);
        dir[celli] = (m > SMALL ? dir[celli]/m : vector::zero);
    }
    scalar unplaced = 0;

    for (label pass = 0; pass <= nPass; ++pass)
    {
        scalarField excess(dVol.size(), Zero);
        forAll(dVol, celli)
        {
            const scalar lo = -a[celli]*V[celli];          // max removal
            const scalar hi = (1 - a[celli])*V[celli];     // max addition
            if (dVol[celli] < lo)
            {
                excess[celli] = dVol[celli] - lo;
                dVol[celli] = lo;
            }
            else if (dVol[celli] > hi)
            {
                excess[celli] = dVol[celli] - hi;
                dVol[celli] = hi;
            }
        }

        if (pass == nPass)
        {
            unplaced += sum(mag(excess));
            break;
        }

        // Pass each excess to the local face neighbour best able to take it
        forAll(excess, celli)
        {
            const scalar e = excess[celli];
            if (e == 0) continue;
            // sinks move towards phase 1 (+dir), sources towards phase 2
            const vector d = (e < 0 ? 1 : -1)*dir[celli];
            label best = -1;
            scalar bestScore = -GREAT;
            for (const label facei : cells[celli])
            {
                if (!mesh.isInternalFace(facei)) continue;
                const label nb = (own[facei] == celli ? nei[facei] : own[facei]);
                // remaining capacity of the neighbour after its own request
                const scalar cap =
                    e < 0
                  ? a[nb]*V[nb] + dVol[nb]
                  : (1 - a[nb])*V[nb] - dVol[nb];
                if (cap <= 0) continue;
                const vector cn(mesh.C()[nb] - mesh.C()[celli]);
                // alignment first, capacity as tie-break
                const scalar score = (cn/mag(cn) & d) + 1e-6*cap/V[nb];
                if (score > bestScore)
                {
                    bestScore = score;
                    best = nb;
                }
            }
            if (best >= 0)
            {
                dVol[best] += e;
                if (mag(dir[best]) < 0.5)
                {
                    dir[best] = dir[celli];
                }
            }
            else
            {
                unplaced += mag(e);
            }
        }
    }

    Su.field() = dVol/(V*dt);

    const scalar requested = gSum(S*V);
    Info<< "alphaSuSp: requested " << requested << " m3/s, placed "
        << gSum(Su.field()*V) << " m3/s, unplaced "
        << returnReduce(unplaced, sumOp<scalar>())/dt << " m3/s" << endl;
}
