# To be run in the directory where the tracer HDF5 files are saved
# creates a copy of the old file with TimestepLastIndex renamed to NextTimestepStartIndex
# copy has format tracer_file_prefix_new_XXX.hdf5 during write
# copy is renamed on success to tracer_file_prefix_XXX.hdf5, whilst the source file is renamed to tracer_file_prefix_old_XXX.hdf5

import h5py
import os

last_tracer_no = 1e9                      # Last tracer snap to read in (alternatively choose a very large number at the script will crash when it runs out of files to convert)
tracer_file_prefix =  "./tracer_file"     # called InputTracerDataFileBase in the CREST parameter file

# Loop over tracer files
for i in range(last_tracer_no+1):

  # Data we will read from / write to
  read_filepath = "{}_{:03d}.hdf5".format(tracer_file_prefix, i)
  write_filepath = "{}_new_{:03d}.hdf5".format(tracer_file_prefix, i)

  f_src  = h5py.File(read_filepath,'r')
  f_dest = h5py.File(write_filepath,'w')

  # If this if the first file we also need the Header data
  if(i == 0):
    f_src.copy(f_src["Header"], f_dest, "Header")

  f_src.copy(f_src["TracerData"], f_dest, "TracerData")

  # Rename group
  f_dest["TracerData/NextTimestepStartIndex"] = f_dest["TracerData/TimestepLastIndex"]
  del(f_dest["TracerData/TimestepLastIndex"])

  f_src.close()
  f_dest.close()

  # Change the names
  old_data_filepath = "{}_old_{:03d}.hdf5".format(tracer_file_prefix, i)
  os.rename(read_filepath, old_data_filepath)
  os.rename(write_filepath, read_filepath)
