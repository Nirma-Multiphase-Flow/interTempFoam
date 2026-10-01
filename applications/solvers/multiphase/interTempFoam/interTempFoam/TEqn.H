{
    // Face conductivity of the phase-change model (sharp models: zero on
    // interface-crossing faces, whose conduction is replaced by TSource)
    const surfaceScalarField kappaf(phaseChangePtr->kappaf());

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T)
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
      - fvm::laplacian(kappaf, T)
      + phaseChangePtr->TSource(T)
    );

    TEqn.relax();
    TEqn.solve();
    T.correctBoundaryConditions();

    Info<< "min/max(T) = " << min(T).value() << " " << max(T).value()
        << endl;
}
