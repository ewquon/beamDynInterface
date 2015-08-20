# beamDynInterface
Interface boundary condition to be used with pimpleBeamDynFoam

add the following to Make/files
  $(derivedPoint)/beamDynInterface/beamDynInterfacePointPatchVectorField.C

add the following to Make/options
    -I$(FOAM_SOLVERS)/incompressible/pimpleBeamDynFoam
