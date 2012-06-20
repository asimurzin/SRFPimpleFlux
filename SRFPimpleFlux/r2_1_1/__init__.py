#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV, Andrey SIMURZIN
##


#---------------------------------------------------------------------------
from Foam import ref, man


#---------------------------------------------------------------------------
def _createFields( runTime, mesh ):
    
    ref.ext_Info() << "Reading field p\n" << ref.nl
    p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )
    
    ref.ext_Info() << "Reading field Urel\n" << ref.nl

    Urel = man.volVectorField( man.IOobject( ref.word( "Urel" ),
                                             ref.fileName( runTime.timeName() ),
                                             mesh,
                                             ref.IOobject.MUST_READ,
                                             ref.IOobject.AUTO_WRITE ),
                               mesh )
  
    ref.ext_Info() << "Reading/calculating face flux field phi\n" << ref.nl
    phi = man.surfaceScalarField( man.IOobject( ref.word( "phi" ),
                                                ref.fileName( runTime.timeName() ),
                                                mesh,
                                                ref.IOobject.READ_IF_PRESENT,
                                                ref.IOobject.AUTO_WRITE ), 
                                  man.surfaceScalarField( ref.linearInterpolate( Urel ) & mesh.Sf(), man.Deps( mesh, Urel ) ) )
    
    pRefCell = 0
    pRefValue = 0.0
    
    pRefCell, pRefValue = ref.setRefCell( p, mesh.solutionDict().subDict( ref.word( "PIMPLE" ) ), pRefCell, pRefValue )
    
    laminarTransport = man.singlePhaseTransportModel( Urel, phi )
    
    turbulence = man.incompressible.turbulenceModel.New( Urel, phi, laminarTransport )
    
    ref.ext_Info() << "Creating SRF model\n" << ref.nl
    SRF = man.SRF.SRFModel.New( Urel ) 
    
    sources = man.IObasicSourceList( mesh )
    
    # Create the absolute velocity
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.NO_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            man.volVectorField( Urel() + SRF.U(), man.Deps( Urel, SRF ) ) ) # mixed  calculations

    return p, U, Urel, SRF, phi, turbulence, pRefCell, pRefValue, laminarTransport, sources


#--------------------------------------------------------------------------------------
def _UrelEqn( mesh, pimple, phi, Urel, p, turbulence, SRF, sources ):

    # Relative momentum predictor
    # The initial C++ expression does not work properly, because of
    #  1. turbulence.divDevRhoReff( U ) - changes values for the U boundaries
    #  2. the order of expression arguments computation differs with C++
    # UrelEqn = man.fvVectorMatrix( ref.fvm.ddt( Urel ) + ref.fvm.div( phi, Urel ) + turbulence.divDevReff( Urel )  + SRF.Su(), \
    #                              man.Deps( phi, Urel, turbulence, SRF ) );
    UrelEqn = man.fvVectorMatrix( turbulence.divDevReff( Urel )  + ref.fvm.div( phi, Urel ) + ref.fvm.ddt( Urel ) + SRF.Su(), \
                                  man.Deps( phi, Urel, turbulence, SRF ) );

    UrelEqn.relax()

    sources.constrain( UrelEqn )

    ref.solve( UrelEqn == -ref.fvc.grad( p ) + sources( Urel ) )
    
    return UrelEqn


#--------------------------------------------------------------------------------------
def pEqn( runTime, mesh, pimple, Urel, UrelEqn, phi, p, pRefCell, pRefValue, cumulativeContErr, sources ): 
    rAUrel = 1.0 / ( UrelEqn == sources( Urel ) ).A()
    Urel << rAUrel * UrelEqn.H()
    if ( pimple.nCorrPISO() <= 1 ):
       # UrelEqn.clear()
       pass
       
    phi << ( ref.fvc.interpolate( Urel ) & mesh.Sf() ) + ref.fvc.ddtPhiCorr( rAUrel, Urel, phi )
 
    ref.adjustPhi( phi, Urel, p )

    # Non-orthogonal pressure corrector loop
    while pimple.correctNonOrthogonal():
        #Pressure corrector
        pEqn = ref.fvm.laplacian( rAUrel, p ) == ref.fvc.div( phi )
        pEqn.setReference( pRefCell, pRefValue )
        
        pEqn.solve( mesh.solver( p.select( pimple.finalInnerIter() ) ) )
           
        if pimple.finalNonOrthogonalIter():
           phi -= pEqn.flux()
           pass
        pass
    cumulativeContErr = ref.ContinuityErrs( phi(), runTime, mesh, cumulativeContErr )
    
    # Explicitly relax pressure for momentum corrector
    p.relax()

    Urel -= rAUrel * ref.fvc.grad( p )
    Urel.correctBoundaryConditions()
    
    sources.correct( Urel )

    return cumulativeContErr


#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )

    p, U, Urel, SRF, phi, turbulence, pRefCell, pRefValue, laminarTransport, sources = _createFields( runTime, mesh )
    
    cumulativeContErr = ref.initContinuityErrs()
    
    pimple = man.pimpleControl( mesh )
    
    ref.ext_Info() << "\nStarting time loop\n" <<ref.nl
    
    while runTime.run() :
        adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )

        CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )
      
        runTime = ref.setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )
        
        runTime.increment()
                
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        # --- Pressure-velocity PIMPLE corrector loop
        while pimple.loop():

            UrelEqn = _UrelEqn( mesh, pimple, phi, Urel, p, turbulence, SRF, sources )
            
            # --- Pressure corrector loop
            while pimple.correct():
               cumulativeContErr = pEqn( runTime, mesh, pimple, Urel, UrelEqn, phi, p, pRefCell, pRefValue, cumulativeContErr, sources )
               pass
            
            # Update the absolute velocity
            U << Urel() + SRF.U() # mixed calculations
            
            if pimple.turbCorr():                             
                turbulence.correct()
                pass
            pass

        runTime.write();
        
        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
        pass

    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020101" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   ref.ext_Info() << "\n\n To use this solver it is necessary to SWIG OpenFOAM-2.1.1 or higher\n"
   pass


#--------------------------------------------------------------------------------------
