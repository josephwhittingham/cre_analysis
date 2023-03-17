import numpy as np
import h5py

for i in range(306):
  hf3 = h5py.File("tracer_file_{:03d}.hdf5".format(i),"r")
  idx2 = np.where(hf3['TracerData/ParticleIDs'][()] == 2000055001)
  try:
    time_idx = np.where(hf3['TracerData/TimestepLastIndex'][()] > idx2[0][0])[0][0]
  except:
    time_idx = np.where(hf3['TracerData/TimestepLastIndex'][()] > idx2[0])[0][0]
  print(i, hf3['TracerData/InjectionEnergy'][idx2], np.linalg.norm(np.vstack([hf3['TracerData/MagneticField/X'][idx2], hf3['TracerData/MagneticField/Y'][idx2], hf3['TracerData/MagneticField/Z'][idx2]]), axis=0), hf3['TracerData/Time'][time_idx])
