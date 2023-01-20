with open("21240.out") as file:
  lines = [line.rstrip() for line in file]

tot_num_tracers = 0
tot_num_shocked_tracers = 0
tot_num_active_timesteps = 0
max_num_tracers = 0

for i in range(len(lines)):
  if(len(lines[i]) == 0):   # ignore blank lines
    continue

  elements = lines[i].split()           # look at elements in line (separated by blank space)

  if(elements[0] in ['COSMIC_RAYS_ELECTRONS:']):   # output from COSMIC_RAYS_ELECTRONS
    if(elements[2] in ['tracers']):
      if(elements[3] in ['active']):
        no_tracers = int(elements[1])
        tot_num_tracers += int(elements[1])
        tot_num_active_timesteps += 1
      else:
        tot_num_shocked_tracers += int(elements[1])

  if(no_tracers > max_num_tracers):           # for an estimate if we had recorded every tracer for each timestep
    max_num_tracers = no_tracers

chunk_size = 44             # data size recorded for one tracer (4 bytes * 11 fields ) - standard data
shocked_chunk_size = 32     # data size recorded for one tracer (4 bytes * 8 fields ) - shock acceleration data (only stored for non-zero shock flags)

file_size = tot_num_tracers * chunk_size + tot_num_shocked_tracers * shocked_chunk_size

print("{0:.2E}".format(tot_num_tracers))    # total number of active tracers over course of simulation
print("{0:2f}".format(file_size / 1E9))     # probable tracer data size in GB
print("{0:2f}".format(file_size /(tot_num_active_timesteps * max_num_tracers * (chunk_size + shocked_chunk_size))))   # percentage relative to if we saved every timestep (incl. shock acceleration module data saved)
