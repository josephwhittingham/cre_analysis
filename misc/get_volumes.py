# To be run in /lustre/joseph/Sims/ShocktubeSims/
# script saves volumes for later reading

import os
import freud
import numpy as np
import astropy.units as units
from mpi4py import MPI
from crest import *

# ==============================================================================
# Start up with parameter loading
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
freud.parallel.set_num_threads(size)

# ==============================================================================
# Setting size of high-res box used in simulation (in kpc)
# (will be used to set size of Freud box, within which volumes are found)

#number of gas cells (each direction)
nx = 90
ny = 60
nz = 60

#box size (original size)
long_x = nx/ny * 4
long_y = 1.
long_z = 1.

box_size = 300

sx = long_x * box_size
sy = long_y * box_size
sz = long_z * box_size


# ==============================================================================

# Crest snapshot prefix
cre_filebase = "Snap_CRe_"

# Simulation variations and directories
mach_variation   = ["output_crest", "output_crest--mach-no-3"]

sim_dir = \
[
#"Flat-noCRs-ramped_buffer-cartesian-equal_density",
"Flat-noCRs-ramped_buffer-cartesian-equal_density-b_rms_1.6e-07",
#"Turb-noCRs-ramped_buffer-cartesian-equal_density-inj_150.0-sigma_3.9e-07-powerlaw_-3.666667",
#"Turb-noCRs-ramped_buffer-cartesian-equal_density-inj_150.0-sigma_3.9e-07-powerlaw_-3.666667-b_rms_1.6e-07",
#"Turb-noCRs-ramped_buffer-cartesian-equal_density-inj_150.0-sigma_1.9e-07-powerlaw_-3.666667-b_rms_1.6e-07",
#"Turb-noCRs-ramped_buffer-cartesian-equal_density-inj_150.0-sigma_7.7e-07-powerlaw_-3.666667-b_rms_1.6e-07",
#"Turb-noCRs-ramped_buffer-cartesian-equal_density-inj_150.0-sigma_3.9e-07-powerlaw_-2.833333-b_rms_1.6e-07",
#"Turb-noCRs-ramped_buffer-cartesian-equal_density-inj_150.0-sigma_3.9e-07-powerlaw_-5.333333-b_rms_1.6e-07",
#"Turb-noCRs-ramped_buffer-cartesian-equal_density-inj_75.0-sigma_3.9e-07-powerlaw_-3.666667-b_rms_1.6e-07",
#"Turb-noCRs-ramped_buffer-cartesian-equal_density-inj_200.0-sigma_3.9e-07-powerlaw_-3.666667-b_rms_1.6e-07"]
]

def tracer_volumes(snap_no, crest_file):
  """
  Use Freud functionality to return Voronoi volumes
  Return: volume associated with each (injected) tracer
  """

  # Convert tracer positions from "cm" to "kpc"
  tracer_pos = (crest_file.pos * units.cm).to(units.kpc).value

  # Centre the box on [0,0,0]
  tracer_pos = np.subtract(tracer_pos, [sx/2, sy/2, sz/2])

  # Create "Freud" box
  box = freud.box.Box(Lx=sx, Ly=sy, Lz=sz)

  # Compute Voronoi volumes
  voro = freud.locality.Voronoi()
  voro.compute((box, tracer_pos))
  tracer_vol = np.copy(voro.volumes)

  MPI.Finalize()

  return tracer_vol


def load_tracer_volumes(snap_no, directory, mach_variation, crest_file):
    """
    Load volume data or create and save if not available
    """

    # Path where volume data is saved
    volume_dir = directory + "/" + mach_variation + "/volumes/"
    volume_filename = "volume_{:03d}.npy".format(snap_no)

    # Make path if it doesn't exist
    if not(os.path.exists(volume_dir)):
        os.makedirs(volume_dir)

    # If the file exists, load it. Otherwise, make the data and save it.
    if not(os.path.exists(volume_dir + volume_filename)):

        tracer_vol = tracer_volumes(snap_no, crest_file)
        np.save(volume_dir + volume_filename, tracer_vol)
    else:
        tracer_vol = np.load(volume_dir + volume_filename)

    return tracer_vol


if __name__ == "__main__":
  """
  Loop over directories and save tracer data for all snapshots
  """
  for i in range(len(sim_dir)):
      print("Reading from: {}".format(sim_dir[i]))

      for j in range(2):
          print("  Current run variation: {}".format(mach_variation[j]))

          for snap_no in range(25+1):
              print("    Snap no: {}".format(snap_no))
              crest_file   = CrestSnapshot("./{}/{}/{}{:03d}.dat".format(sim_dir[i], mach_variation[j], cre_filebase, snap_no), high_tracer_number=True)

              tracer_vol = load_tracer_volumes(snap_no, sim_dir[i], mach_variation[j], crest_file)
