import numpy as np
import struct
import argparse
import time
#import IPython
from os.path import isdir

# command line arguments
parser = argparse.ArgumentParser(description="Create binary file containing tracer initial conditions (IDs and positions)")
parser.add_argument('--output_dir', default="./input", type=str, help="Output folder for the tracer file")
parser.add_argument('--num_tracer_x', default = 10, help="Number of tracers in x-direction", type = int)
parser.add_argument('--num_tracer_y', default = 10, help="Number of tracers in y-direction", type = int)
parser.add_argument('--num_tracer_z', default = 10, help="Number of tracers in z-direction", type = int)
parser.add_argument('--box_size', help="Size of sim box (use --long_x, --long_y, --long_z for cuboid boxes). Must match parameter file of sim.", type = float)
parser.add_argument('--long_x', default = 1, help="Elongation in x-direction (multiple of box size). Must match parameter file of sim.", type = float)
parser.add_argument('--long_y', default = 1, help="Elongation in y-direction (multiple of box size). Must match parameter file of sim.", type = float)
parser.add_argument('--long_z', default = 1, help="Elongation in z-direction (multiple of box size). Must match parameter file of sim.", type = float)
args = parser.parse_args()

t_start = time.perf_counter()

if not isdir(args.output_dir):
    print("Given output directory:{}, does not exist!\n".format(output_dir))
    exit()

if(args.box_size is None):
    print("Must give box size as input parameter! (e.g., python create_tracer_ics --box_size=x)")
    exit()

tot_num_tracers = args.num_tracer_x * args.num_tracer_y * args.num_tracer_z

# actual box size
side_x = args.long_x * args.box_size
side_y = args.long_y * args.box_size
side_z = args.long_z * args.box_size

# calculate tracer spacing
dx = side_x / args.num_tracer_x
dy = side_y / args.num_tracer_y
dz = side_z / args.num_tracer_z

# tracer positions along each axis (including shift)
tracer_x_pos = np.arange(0, side_x, dx) + 0.5*dx
tracer_y_pos = np.arange(0, side_y, dy) + 0.5*dy
tracer_z_pos = np.arange(0, side_z, dz) + 0.5*dz

t_1 = time.perf_counter()

# make list of positions in the way Arepo requires
meshed_grid = np.meshgrid(tracer_x_pos, tracer_y_pos, tracer_z_pos)

t_2 = time.perf_counter()

positions = np.swapaxes(np.vstack((meshed_grid[0].flatten(), meshed_grid[1].flatten(), meshed_grid[2].flatten())), 0, 1)

t_3 = time.perf_counter()

# give IDs to tracers
IDOffset = 2000000000       #this is from Arepo -  makes sure tracer IDs are unique from other particle IDs
IDs =  np.arange(IDOffset, IDOffset + tot_num_tracers)

t_4 = time.perf_counter()

# write total number of tracers, IDs, and positions to binary file
file_name = '{}/tracer_ics--box_size-{}-{}-{}--num_tracers-{}-{}-{}.bin'.format(args.output_dir, int(side_x), int(side_y), int(side_z), args.num_tracer_x, args.num_tracer_y, args.num_tracer_z)

file = open(file_name, "wb")
file.write(struct.pack('i', tot_num_tracers))
for i in range(tot_num_tracers):
    file.write(struct.pack('d', positions[i,0]))
    file.write(struct.pack('d', positions[i,1]))
    file.write(struct.pack('d', positions[i,2]))
for i in range(tot_num_tracers):
    file.write(struct.pack('i', int(IDs[i])))           # IDs are given last so they can be ignored if need be
file.close()

t_end = time.perf_counter()

#IPython.embed()

print("\nWrote initial positions for {:d} particles to file: {:s}\n".format(tot_num_tracers, file_name))
print("------ Timings (in secs) ------\nArrange axes:\t\t{:.3g}\nMesh grid:\t\t{:.3g}\nSwap axes:\t\t{:.3g}\nMake list of IDs:\t{:.3g}\nWrite to file:\t\t{:.3g}\nTotal:\t\t\t{:.3g}\n".format(t_1-t_start, t_2-t_1, t_3-t_2, t_4-t_3, t_end-t_4, t_end-t_start))
