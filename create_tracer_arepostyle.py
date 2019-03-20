import numpy as np
import struct
import sys
import os
sys.path.append('../../CRelectron/Analysis/')
from ParamClass import *
from Physics import *

fname = 'tracer.dat'

size_i = struct.calcsize('i')
if size_i != 4:
    print("Type 'int' is not 32 bit on this platform")

size_f = struct.calcsize('f')
if size_f !=4:
    print("Type 'float' is not 32 bit on this platform")
    
size_d = struct.calcsize('d')
if size_d !=8:
    print("Type 'double' is not 32 bit on this platform")


# Number of Snapshots and tracer particles
nSnap = 9 + 1
nPart = 1

# Create necessary arrays
time	= np.zeros(nSnap)
pos_x = np.ndarray((nPart,nSnap),dtype=float)
pos_y = np.ndarray((nPart,nSnap),dtype=float)
pos_z = np.ndarray((nPart,nSnap),dtype=float)
rho = np.ndarray((nPart,nSnap),dtype=float)
temp = np.ndarray((nPart,nSnap),dtype=float)
Utherm = np.ndarray((nPart,nSnap),dtype=float)
b_field = np.ndarray((nPart,nSnap),dtype=float)
dist_to_shock = np.ndarray((nPart,nSnap),dtype=float)
compression_ratio = np.ndarray((nPart,nSnap),dtype=float)
b_field_ratio = np.ndarray((nPart,nSnap),dtype=float)
eInjection = np.ndarray((nPart,nSnap),dtype=int)
eEnergyPerMass = np.ndarray((nPart,nSnap),dtype=float)
V_pre_Shock = np.ndarray((nPart,nSnap),dtype=int)

# Fill the arrays with something (everything in cgs)
time = np.array([1.e7 * t * YEARS_IN_SECONDS for t in [0,10.,100.,400.,700.,700.8,701.6,703.2,706.4,712.8]])

for p in np.arange(nPart):
    eInjection[p]		= np.array([1,1,1,1,0,0,0,0,0,0])

pos_x.fill(0.5)
pos_y.fill(0.5)
pos_z.fill(0.5)
rho.fill(1.e-3)
temp.fill(1.e5)
Utherm.fill(1.)
b_field.fill(1.e-6)
dist_to_shock.fill(0.)
compression_ratio.fill(4.)
b_field_ratio.fill(1.)
eEnergyPerMass.fill(1.e10)
V_pre_Shock.fill(CLIGHT/6000)

float_buffer = np.ndarray(nPart,dtype=float)
int_buffer = np.ndarray(nPart,dtype=int)

BaseDir = '../'  # main folder from where Input, Output, etc can accessed

with open(fname,'wb') as f:
    dummy = nPart * (12 * size_f + size_i) + 8
    f.write(struct.pack('i',dummy))

    for s in np.arange(nSnap):
        f.write(struct.pack('d',time[s]))
        float_buffer[:] = pos_x[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = pos_y[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = pos_z[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = rho[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = temp[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = Utherm[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = b_field[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        int_buffer[:] = eInjection[:,s]
        f.write(struct.pack('{:d}i'.format(nPart),*int_buffer))
        float_buffer[:] = eEnergyPerMass[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = dist_to_shock[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = compression_ratio[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = b_field_ratio[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))
        float_buffer[:] = V_pre_Shock[:,s]
        f.write(struct.pack('{:d}f'.format(nPart),*float_buffer))

        f.write(struct.pack('i',dummy))
        if s < nSnap - 1:
            f.write(struct.pack('i',dummy))            
        
f.close()

