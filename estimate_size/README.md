These scripts use the Arepo .out and cpu.txt files respectively to give an estimate of how large the corresponding tracer HDF5 file will be.
(Useful if the simulation doesn't take too long, but you are worried about space constraints)

This requires running the simulation with:
 - COSMIC_RAYS_ELECTRONS in place, but having commented out lines where CREST writes to file.
 - TRACER_PARTICLE in place (COSMIC_RAYS_ELECTRONS not required)
(respectively)

The latter option does not account for data stored from COSMIC_RAYS_SHOCK_ACCELERATION.

The files will need to be edited if you will be storing extra data (e.g. using config flags COSMIC_RAYS_SN_INJECTION, COSMIC_RAYS_MAGNETIC_OBLIQUITY, or simply editing the existing fields)
