{
        // update the boundary absolute velocity
        mrfZones.correctBoundaryVelocity(U);

        // set the pseudo face velocity
        Godunov.dotX() = mrfZones.faceU();

        // activate sub-cycle loop
        TimeState subCycleTimeState = runTime.subCycle(numberSubCycles);

        // get access to physical deltaT in ddt scheme!
        physDeltaT[1] = subCycleTimeState.deltaT().value();
        physDeltaT[2] = subCycleTimeState.deltaT0().value();

        // Begin sub-cycling - PseudoTime Integration
        // adjustTimeStep and numberSubCycles > 1 does not make any sense
        for (label j=0; j < numberSubCycles; j++)
        {
            // Update subCycle Time
            runTime++;

            // update local time step sizes just once each iteration for
            // all ddt schemes
            if (!adjustTimeStep)
            {  
                localTimeStep_Real.update(maxCo,adjustTimeStep);
            }

            // Plot Pseudo Time here, so that one can analyse the residuls with
            // pyFoam*
            Info<< "\n Time = " << runTime.value() << nl << endl;

            // Begin outer Loop for Runge - Kutta
            forAll (beta,i)
            {
                // Store the previous iteration data for p to allow utility of waveTransmissive BC
                p.storePrevIter();

                // Update primitive boundaries
                p.correctBoundaryConditions();
                U.correctBoundaryConditions();
                h.correctBoundaryConditions();
                e.correctBoundaryConditions();
 

                Godunov.update(secondOrder);

                // get access to multi-stage coefficient in ddt scheme
                physDeltaT[0] = beta[i];

                // \f$ \mu \left( \nabla \vec{U} + \left( \nabla \vec{U}
                // \right)^T - \frac{2}{nD} \left( \nabla \bullet \vec{U}
                // \right) I \right) \f$
                // nD = 3 in three dimensional space
                // for two- and one-dimensional computations this
                // implementation is wrong
                // is equal to -turbulence->devRhoReff(), but we do not need to
                // evaluate gradU twice
                const volSymmTensorField tau
                (
                    "tau",
                    -turbulence->devRhoReff()
                    -((2.0/3.0)*I)*rho*turbulence->k()
                );

                volScalarField rhoFlux = -fvc::div(Godunov.rhoFlux());

                volVectorField rhoUFlux = -fvc::div(Godunov.rhoUFlux())
                    + fvc::div(tau);

                volScalarField rhoEFlux = -fvc::div(Godunov.rhoEFlux())
                    // Viscous heating with
                    + fvc::div( tau & U )
                    // Fourier law with static enthalpy
                    // with alphaEff() - the effective turbulent
                    // thermal diffusivity.
                    + fvc::laplacian(turbulence->alphaEff(), h)
                    // Molecular Diffusion and Turbulent Transport closure
                    // Wilcox (2006): Chapter 5.4.3
                    // should be better used DkEff(F1) instead of muEff(), but
                    // this function is not virtual, now it is assumed that
                    // \sigma_k = 5/3 is hard coded
                    + fvc::laplacian
                      (
                          (turbulence->mu()+0.6*turbulence->mut()),
                          turbulence->k()
                      )
                    ;

                // Till here, the rho internal field is updated.

                // Add source terms to the numerical fluxes due to the rotation
                // Is not needed if mesh is really rotating, as this effect is
                // captured then by the time stepping of the velocity field. As
                // the vector is not adapted, when the mesh is rotating.
                mesh.type() == "staticFvMesh" ?
                    mrfZones.addCoriolis(rho,U,rhoUFlux):void();

                // Pseudo time integration
                solve
                (
                    fvm::ddt(rho) == rhoFlux
                );

                // time integration
                solve
                (
                    fvm::ddt(rhoU) == rhoUFlux
                );

                // time integration
                solve
                (
                    fvm::ddt(rhoE) == rhoEFlux
                );
                
                // bound density
                boundMinMax(rho,rhoMin,rhoMax);

                // bound rhoE
                boundMinMax(rhoE,rhoEMin,rhoEMax);

                // Compute internal field of U
                U.dimensionedInternalField() =
                    rhoU.dimensionedInternalField()
                   /rho.dimensionedInternalField();


                // Update boundary conditions of U
                U.correctBoundaryConditions();

                // Bound the velocity
                volScalarField magU = mag(U);

                if (max(magU) > UMax)
                {

                    volScalarField Ulimiter = pos(magU - UMax)*UMax
                        /(magU + smallU) + neg(magU - UMax);
                    Ulimiter.max(scalar(0));
                    Ulimiter.min(scalar(1));

                    U *= Ulimiter;
                    U.correctBoundaryConditions();
                }

                // Update static energy:
                // The turbulent kinetic energy k is part of the total energy E
                // Therefore it needs to be subtracted from E in order to get
                // the static energy e

                // Updating the static energy field with updated rho internal field.
                e = rhoE/rho - 0.5*magSqr(U) - turbulence->k(); 
                // Due to the e has boundary conditions, thus we can update static
                // energy boundary values with 
                e.correctBoundaryConditions();
                // bound static energy
                boundMinMax(e,eMin,eMax);


                // Inorder to use same tables, the pressure is updated
                // with secant method.
                // The table used here is rho = L(e,p) table
                forAll(p,celli)
                 {   
                    //target density
                    scalar rhof = rho[celli]; 
                    // x_n value
                    scalar p2 = p[celli];
                    // xn-1 value, with small increment of p2. 
                    scalar p1 = p2 + 0.00001*p2; 
                    // initialize x_final here, which is pf.
                    scalar pf = p2;
 
                    scalar rho2  = rhoepTable(e[celli],p2);
                    scalar rho1  = rhoepTable(e[celli],p1);


                    scalar drho = mag(rho2 - rho1);
                    scalar e_celli = e[celli];

                    // start secant method and convergence check here
                    while(drho > 1e-07)
                    {

                        pf     = p2 - (p2-p1)/(rho2-rho1)*(rho2-rhof);

                        p1     = p2; 
                        p2     = pf; 

                        rho2   = rhoepTable(e_celli,p2);
                        rho1   = rhoepTable(e_celli,p1);

                        drho   = mag(rho2 - rho1);

                    }
                    // find the final value and update p field.
                    p[celli]   = pf; 
                }

                // correct p boundary conditions
                p.correctBoundaryConditions();
                // bound pressure
                boundMinMax(p,pMin,pMax);


                // Update static enthalpy field
                h = e + p/rho;
                // Correct static enthalpy boudnary values
                h.correctBoundaryConditions();
                // bound static enthalpy
                boundMinMax(h,hMin,hMax);

                // correct thermo physical properties
                // Temperature fields and boundary values are corrected here.
                thermo.correct();


                // At the final stage, update the acoustic speed
                forAll(A,celli)                                 
                {                                               
                    A[celli] = AepTable(e[celli],p[celli]);     
                }                                              
                // update the boundary value for acoustic speed with
                // look-up table method.   A = L(e,p).
                forAll(A.boundaryField(), patchI)
                {
                    fvPatchScalarField& Ab = A.boundaryField()[patchI];
                    const fvPatchScalarField& eb = e.boundaryField()[patchI];
                    const fvPatchScalarField& pb = p.boundaryField()[patchI];
                    forAll(Ab, faceI)
                    {
                        Ab[faceI] = AepTable(eb[faceI], pb[faceI]);
                    }
                }



                // rho's boundary conditions are given to calculated.
                // Thus force updating of boundary values is required.
                // Here the look-up table  rho = L(e,p) is required.
                // Update boundary field of rho
                forAll(rho.boundaryField(), patchI)
                {

                        fvPatchScalarField& rhob = rho.boundaryField()[patchI];
                   const fvPatchScalarField& pb  = p.boundaryField()[patchI];
                   const fvPatchScalarField& eb  = e.boundaryField()[patchI];
                   forAll(rhob, faceI)
                   {
                       rhob[faceI] = rhoepTable(eb[faceI], pb[faceI]);
                   }

                }
                // The correcBoundaryCondition expression is valied for processor boundary condition
                // to avoid discontinuity when carrying on parallization. 
                rho.correctBoundaryConditions(); 

                // Update boundary field of rhoU
                rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
                

                // Update boundary field of rhoE
                rhoE.boundaryField() =
                    rho.boundaryField()*
                    (
                        e.boundaryField() +
                       0.5*magSqr(U.boundaryField())
                      + turbulence->k()().boundaryField()
                    ) ;

                // needed for turbulence and CoEuler ddt scheme
                // and maybe elsewhere;
                // phi is already realtive to the mesh motion, if we are using
                // the ALE approach
                phi = Godunov.rhoFlux();
//                 mrfZones.relativeFlux(fvc::interpolate(rho),phi);

                // Convergence check
                // Averaged L2-Norm of fluxes
                scalar L2NormRho = max(Foam::sqrt(sum(sqr(rhoFlux.internalField()))
                    /mesh.nCells()),SMALL);
                scalar LInfNormRho = max(mag(rhoFlux.internalField()));
//                 scalar L2NormRhoU = max(Foam::sqrt(sum(magSqr(rhoUFlux.internalField()))
//                     /mesh.nCells()),SMALL);
//                 scalar L2NormRhoE = max(Foam::sqrt(sum(sqr(rhoEFlux.internalField()))
//                     /mesh.nCells()),SMALL);

  

                // Averaged L2-Norm of fluxes
                Info<< "rho L2 Residual: "
                    << Foam::log10(L2NormRho)  << endl
                    << "rho LInf Residual: "
                    << Foam::log10(LInfNormRho) << endl
//                     << "rhoU residual: "
//                     << Foam::log10(L2NormRhoU) << endl
//                     << "rhoE residual: "
//                     << Foam::log10(L2NormRhoE) << endl
                    << endl;

            // End outer Loop for Runge - Kutta
            }

            // Update relative velocity
            mrfZones.relativeVelocity(U,URel);

            //Update turbulence after the multi-stage time integration
            //Info <<"before turbulence" << endl;
            turbulence->correct();
            //Info << "afeer turbuelnce" << endl;
        }
        // End dual-time stepping
        runTime.endSubCycle();
#       include "updateDualTimeSteppingFields.H"
}
