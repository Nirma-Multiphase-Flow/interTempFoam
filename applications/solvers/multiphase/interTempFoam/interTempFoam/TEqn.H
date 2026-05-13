surfaceScalarField kappaf = mixture.kappaf();

fvScalarMatrix TEqn
(
    fvm::ddt(rhoCp, T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
);

TEqn.relax();

TEqn.solve();

T.correctBoundaryConditions();

Info<< "min/max(T) = "
    << min(T).value()
    << " "
    << max(T).value()
    << endl;