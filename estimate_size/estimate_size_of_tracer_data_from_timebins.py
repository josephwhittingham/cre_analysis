with open("timebins.txt") as file:
  lines = [line.rstrip() for line in file]

tot_num_tracers = 0
tot_num_active_timesteps = 0
max_num_tracers = 0

for i in range(len(lines)):
  if(len(lines[i]) == 0):   # ignore blank lines
    continue

  elements = lines[i].split()           # look at elements in line (separated by blank space)

  no_tracers = 0

  if(elements[0] in ['PM-Step.','Total']):   # this is the summary for the current timebin
    no_tracers = int(elements[4])
    tot_num_active_timesteps += 1

  tot_num_tracers += no_tracers

  if(no_tracers > max_num_tracers):           # for an estimate if we had recorded every tracer for each timestep
    max_num_tracers = no_tracers

chunk_size = 40         # data size recorded for one tracer (4 bytes * 10 fields )

print("{0:.2E}".format(tot_num_tracers))    # total number of active tracers over course of simulation
print("{0:2f}".format(tot_num_tracers * chunk_size / 1E9))    # probable tracer data size in GB
print("{0:2f}".format(tot_num_tracers/(tot_num_active_timesteps * max_num_tracers)))   # percentage relative to if we saved every timestep (excl. shock acceleration module data saved)
