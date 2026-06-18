{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh,
                mesh.ddtScheme("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            if (nAlphaSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                alphaRestart
             || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Add the optional shear compression contribution
    if (scAlpha > 0)
    {
        phic +=
            scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }


    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi);

    // For coldPhaseIsHighAlpha1 (melting), nHat points from liquid into solid
    // so the raw Stefan flux is negative (−x).  Negating it makes MULES advect
    // alpha.solid in the +x direction, correctly eroding the solid.
    const surfaceScalarField phiStefan
    (
        (coldPhaseIsHighAlpha1 ? scalar(-1) : scalar(1))
      * phaseChangePtr->phiStefan()
    );

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN =
            cnCoeff*(phi)
        + (1.0 - cnCoeff)
        *(phi.oldTime());
    }
    else
    {
        phiCN = phi;
    }

    if (MULESCorr)
    {
        #include "alphaSuSp.H"

        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
         ==
            Su + fvm::Sp(Sp + divU, alpha1) - divU*alpha1.oldTime()
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        alphaPhi10 = talphaPhi1UD();

        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                alpha1,
                alphaPhi10,
                talphaPhi1Corr0.ref(),
                oneField(),
                zeroField()
            );

            alphaPhi10 += talphaPhi1Corr0();
        }

        // Cache the upwind-flux
        talphaPhi1Corr0 = talphaPhi1UD;

        alpha1.max(0.0);
        alpha1.min(1.0);
        alpha1.correctBoundaryConditions();

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }


    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        #include "alphaSuSp.H"

        surfaceScalarField phir(phic*mixture.nHatf());

        // Pure conservative convection: fvmDiv(phi, alpha1) is exact flux-form
        // (Gauss's theorem), so only the boundary flux and the explicit
        // Su/Sp = alpha1Gen (RDF-Stefan mdot/rho_l) source can change
        // total liquid volume -- no dilatation-cancellation term is used.
        //
        // HISTORY: an earlier attempt added a -fvm::Sp(divU, alpha1) term
        // (divU = fvc::div(phi) = PCV) to the predictor in the style of
        // interPhaseChangeFoam's "Sp + divU".  That pattern is only safe
        // there because its alpha-source coefficient (vDotAlphal, see
        // phaseChangeTwoPhaseMixture.C) is itself built from the SAME
        // (1/rho_l - 1/rho_v) dilatation-scale factor as divU, so adding
        // them is a self-consistent rewrite.  Our alpha1Gen = mdot/rho_l
        // (pure liquid mass-loss rate) is NOT that scale -- it is
        // ~1626x smaller than PCV = mdot*(1/rho_l-1/rho_v) for water/steam.
        // Adding divU to Sp therefore injected a real, uncompensated
        // +alpha1*divU source (large and positive exactly where the
        // droplet sits), which is what caused the 35-38x volume growth.
        // Converting fvmDiv(phi,a)-Sp(divU,a) to advective form and back
        // out is only a no-op if the RHS also carries the matching
        // -Sp(divU,a); doing that algebra exactly cancels divU on both
        // sides, leaving plain fvmDiv(phi,a) == Su + Sp(Sp,a) -- i.e. no
        // divU term at all is the exact, mass-conservative form here.
        alphaPhiUn =
                fvc::flux
                (
                    phi,
                    alpha1,
                    alphaScheme
                )
            + fvc::flux
                (
                -fvc::flux(-phir, alpha2, alpharScheme),
                    alpha1,
                    alpharScheme
                );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhi1Corr(alphaPhiUn - alphaPhi10);

            volScalarField alpha10("alpha10", alpha1);

            // Pass Sp and the matching Su = -Sp*alpha1 so the MULES limiter
            // is source-aware.  Canonical interFoam VoF pattern:
            //   VoF/alphaEqn.H MULESCorr branch, OpenFOAM v2412.
            // Without these args the limiter did not account for the
            // consumed-phase sink and the antidiffusive flux injected
            // spurious liquid (α>0.5 volume grew 38× in ~0.1 ms).
            MULES::correct
            (
                geometricOneField(),
                alpha1,
                alphaPhiUn,
                talphaPhi1Corr.ref(),
                Sp,
                (-Sp*alpha1)(),
                oneField(),
                zeroField()
            );

            alpha1.max(0.0);
            alpha1.min(1.0);
            alpha1.correctBoundaryConditions();

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi10 += talphaPhi1Corr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi10 += 0.5*talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi10 = alphaPhiUn;

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                alphaPhi10,
                Sp,
                Su,
                oneField(),
                zeroField()
            );

            
        }
        alpha1.max(0.0);
        alpha1.min(1.0);
        alpha1.correctBoundaryConditions();

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhi1Corr0 = alphaPhi10 - talphaPhi1Corr0;
        talphaPhi1Corr0.ref().rename("alphaPhi1Corr0");
    }
    else
    {
        talphaPhi1Corr0.clear();
    }

    #include "rhofs.H"

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
     || word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::localEulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
        // rhoCpPhi uses fluid flux only: Stefan velocity carries latent heat
        // (handled by Q_pc/Qcorr in TEqn), not sensible enthalpy.
        // Using alphaPhi10 (which includes phiStefan) would inject
        // ρCp_liq*T_liq into low-ρCp vapor cells → catastrophic T spike.
        rhoCpPhi = phi * fvc::interpolate(rhoCp);
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi10 =
                (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
        rhoCpPhi = phi * fvc::interpolate(rhoCp);
    }

    

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
