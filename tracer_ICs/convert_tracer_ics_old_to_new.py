# Converts the old style tracer ICs into the new style (i.e. removes masses)

import struct
file = open("arepo_tracer_initial_positions_0100_xParticles.bin","rb")
file_write = open("arepo_tracer_initial_positions_0100_xParticles-new.bin","wb")

tot_num_tracers = struct.unpack('i', file.read(4))[0]
file_write.write(struct.pack('i', tot_num_tracers))

for i in range(tot_num_tracers):
  struct.unpack('d', file.read(8))[0]   # mass (ignore this)
  file_write.write(struct.pack('d', struct.unpack('d', file.read(8))[0]))  # coordinates (write these to new file)
  file_write.write(struct.pack('d', struct.unpack('d', file.read(8))[0]))
  file_write.write(struct.pack('d', struct.unpack('d', file.read(8))[0]))

file.close()
file_write.close()
