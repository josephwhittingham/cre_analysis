import numpy as np
import h5py

for i in range(306):
  hf = h5py.File("tracer_file_{:03d}.hdf5".format(i),"r")
  idx = np.where(hf['TracerData/ParticleIDs'][()] == 2000055001)
  try:
    time_idx = np.where(hf['TracerData/TimestepLastIndex'][()] >= idx[0][0])[0][0]
  except:
    time_idx = np.where(hf['TracerData/TimestepLastIndex'][()] >= idx[0])[0][0]
  print(i, hf['TracerData/InjectionEnergy'][idx], np.linalg.norm(np.vstack([hf['TracerData/MagneticField/X'][idx], hf['TracerData/MagneticField/Y'][idx], hf['TracerData/MagneticField/Z'][idx]]), axis=0), hf['TracerData/Time'][time_idx])
