#
# Script to isolate the data from a single tracer ID
#

import h5py
import numpy as np
import fnmatch
import os
import datetime

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# TODO: JW: Add config options here? (rather than using the try: except: approach)
make_single_ID_file = True

# ID to find
ID =  # <-- this tracer!

# Loop over snapshots until this one
snap_no_last = len(fnmatch.filter(os.listdir("."), 'tracer_file*'))     # All snapshots in file
arepo_no_last = len(fnmatch.filter(os.listdir("."), 'snapshot*'))     # All snapshots in file
#snap_no_last = 293 #0 +1
start_time = datetime.datetime.now()
print("Start:", start_time)
print("Last tracer file found: {:d}".format(snap_no_last))

# Arrays to fill
times = np.array([])
next_timestep_start_index = np.array([])
particle_ids = np.array([])
pos_x = np.array([])
pos_y = np.array([])
pos_z = np.array([])
mag_x = np.array([])
mag_y = np.array([])
mag_z = np.array([])
density = np.array([])
internal_energy = np.array([])
mag_field_strength = np.array([])
shock_flag = np.array([])
CRp_injected_energy = np.array([])
mag_obliquity = np.array([])
post_shock_density = np.array([])
post_shock_speed = np.array([])
pre_shock_density = np.array([])
shock_crossing_time = np.array([])
shock_direction_x = np.array([])
shock_direction_y = np.array([])
shock_direction_z = np.array([])
mach_number = np.array([])

arepo_snap_no = 0

def get_arepo_time(arepo_snap_no):
    hf_arepo = h5py.File("snapshot_{:03d}.hdf5".format(arepo_snap_no), "r")
    arepo_time = hf_arepo["Header"].attrs["Time"]
    return arepo_time

arepo_time = get_arepo_time(arepo_snap_no)
old_arepo_time = arepo_time
eps_time = 0

print("\nFinding data for ID: {}\n".format(ID))

for snap_no in range(snap_no_last):
    print("Collecting data from snap: {:d}".format(snap_no))

    hf = h5py.File("tracer_file_{:03d}.hdf5".format(snap_no), "r")

    particle_ids_all = hf["TracerData/ParticleIDs"][()]
    idx = np.where(particle_ids_all == ID)[0]
    particle_ids = np.append(particle_ids, particle_ids_all[idx])

    next_timestep_start_index = np.append(0, hf["TracerData/NextTimestepStartIndex"][()])

    time = np.array([])

    for i in range(len(idx)):
        time = np.append(time, hf["TracerData/Time"][np.where(idx[i] < next_timestep_start_index)[0][0] -1])

    times = np.append(times, time)

    # Get variables that are always written out
    mag_field_strength = np.append(mag_field_strength, hf["TracerData/MagneticFieldStrength"][idx])
    density = np.append(density, hf["TracerData/Density"][idx])
    internal_energy = np.append(internal_energy, hf["TracerData/InternalEnergy"][idx])

    # Only captured if we're running with COSMIC_RAYS_SHOCK_ACCELERATION
    try:
        # Get shock data
        shock_flag_all = hf["TracerData/ShockFlag"][()]
        shock_idx = np.where(shock_flag_all > 1)
        shock_flag = np.append(shock_flag, shock_flag_all[idx])

        _, _, shock_idx_particle = np.intersect1d(idx, shock_idx, return_indices=True)  # Shock data is ordered differently

        if(len(shock_idx_particle) != 0):
            CRp_injected_energy = np.append(CRp_injected_energy, hf["TracerData/CRpInjectedEnergy"][shock_idx_particle])
            pre_shock_density = np.append(pre_shock_density, hf["TracerData/PreShockDensity"][shock_idx_particle])
            post_shock_density = np.append(post_shock_density, hf["TracerData/PostShockDensity"][shock_idx_particle])
            post_shock_speed = np.append(post_shock_speed, hf["TracerData/PostShockSpeed"][shock_idx_particle])
            shock_crossing_time = np.append(post_shock_density, hf["TracerData/ShockCrossingTime"][shock_idx_particle])

            # Only read in sometimes
            try:
                mag_obliquity = np.append(mag_obliquity, hf["TracerData/MagneticObliquity"][shock_idx_particle])
            except:
                pass
            try:
                mach_number = np.append(hf["TracerData/MachNumber"][shock_idx_particle], mach_number)
            except:
                pass
    except:
        print("No shock data")

    # Get parts we only store on snapshot times
    for i in range(len(time)):
        if((time[i] + eps_time >= arepo_time) or ((i == len(time)-1) and (snap_no == snap_no_last-1))):
            idx = np.where(particle_ids_all[next_timestep_start_index[i]: next_timestep_start_index[i+1]] == ID)[0]

            if(len(idx) == 0):
                print("No data for this ID")
                break

            pos_x = np.append(pos_x, hf["TracerData/Coordinates/X"][idx])
            pos_y = np.append(pos_y, hf["TracerData/Coordinates/Y"][idx])
            pos_z = np.append(pos_z, hf["TracerData/Coordinates/Z"][idx])

            mag_x = np.append(mag_x, hf["TracerData/MagneticField/X"][idx])
            mag_y = np.append(mag_y, hf["TracerData/MagneticField/Y"][idx])
            mag_z = np.append(mag_z, hf["TracerData/MagneticField/Z"][idx])

            # Only captured if we're running with COSMIC_RAYS_SHOCK_ACCELERATION
            try:
                shock_direction_x = np.append(shock_direction_x, hf["TracerData/ShockDirection/X"][shock_idx_particle])
                shock_direction_y = np.append(shock_direction_y, hf["TracerData/ShockDirection/Y"][shock_idx_particle])
                shock_direction_z = np.append(shock_direction_z, hf["TracerData/ShockDirection/Z"][shock_idx_particle])
            except:
                pass

            eps_time = 1e-5 * (arepo_time - old_arepo_time)
            old_arepo_time = arepo_time
            arepo_snap_no += 1

            if(arepo_snap_no == arepo_no_last):
                print("No more Arepo snapshots!")
                break

            arepo_time = get_arepo_time(arepo_snap_no)

    hf.close()


if make_single_ID_file == True:
    print("Saving data:", datetime.datetime.now())

    hf = h5py.File("tracer_file_{:03d}.hdf5".format(0), "r")

    # Duplicate into new file
    with h5py.File("ID_{:d}--tracer_file_{:03d}.hdf5".format(ID, 0), "w") as h5w:
        for obj in hf.keys():
            hf.copy(obj, h5w)
    hf.close()

    hf2 = h5py.File("ID_{:d}--tracer_file_{:03d}.hdf5".format(ID, 0), "r+")

    del hf2["Header/AllTracerParticleIDs"]
    hf2.create_dataset("Header/AllTracerParticleIDs", data=[ID], dtype=np.uint32)

    del hf2["TracerData/ParticleIDs"]
    hf2.create_dataset("TracerData/ParticleIDs", data=particle_ids, dtype=np.uint32)

    del hf2["TracerData/Time"]
    hf2.create_dataset("TracerData/Time", data=times, dtype=np.float64)

    del hf2["TracerData/NextTimestepStartIndex"]
    hf2.create_dataset("TracerData/NextTimestepStartIndex", data=np.arange(1, len(times)+1), dtype=np.int32)

    del hf2["TracerData/Coordinates/X"]
    del hf2["TracerData/Coordinates/Y"]
    del hf2["TracerData/Coordinates/Z"]
    hf2.create_dataset("TracerData/Coordinates/X", data=pos_x, dtype=np.float32)
    hf2.create_dataset("TracerData/Coordinates/Y", data=pos_y, dtype=np.float32)
    hf2.create_dataset("TracerData/Coordinates/Z", data=pos_z, dtype=np.float32)

    del hf2["TracerData/MagneticField/X"]
    del hf2["TracerData/MagneticField/Y"]
    del hf2["TracerData/MagneticField/Z"]
    hf2.create_dataset("TracerData/MagneticField/X", data=mag_x, dtype=np.float32)
    hf2.create_dataset("TracerData/MagneticField/Y", data=mag_y, dtype=np.float32)
    hf2.create_dataset("TracerData/MagneticField/Z", data=mag_z, dtype=np.float32)

    del hf2["TracerData/MagneticFieldStrength"]
    del hf2["TracerData/Density"]
    del hf2["TracerData/InternalEnergy"]
    hf2.create_dataset("TracerData/MagneticFieldStrength", data=mag_field_strength, dtype=np.float32)
    hf2.create_dataset("TracerData/Density", data=density, dtype=np.float32)
    hf2.create_dataset("TracerData/InternalEnergy", data=internal_energy, dtype=np.float32)

    if(len(CRp_injected_energy) != 0):
        del hf2["TracerData/ShockFlag"]
        del hf2["TracerData/CRpInjectedEnergy"]
        del hf2["TracerData/PreShockDensity"]
        del hf2["TracerData/PostShockDensity"]
        del hf2["TracerData/PostShockSpeed"]
        del hf2["TracerData/ShockCrossingTime"]

        hf2.create_dataset("TracerData/ShockFlag", data=shock_flag, dtype=np.int32)
        hf2.create_dataset("TracerData/CRpInjectedEnergy", data=CRp_injected_energy, dtype=np.float32)
        hf2.create_dataset("TracerData/PreShockDensity", data=pre_shock_density, dtype=np.float32)
        hf2.create_dataset("TracerData/PostShockDensity", data=post_shock_density, dtype=np.float32)
        hf2.create_dataset("TracerData/PostShockSpeed", data=post_shock_speed, dtype=np.float32)
        hf2.create_dataset("TracerData/ShockCrossingTime", data=shock_crossing_time, dtype=np.float32)

        del hf2["TracerData/ShockDirection/X"]
        del hf2["TracerData/ShockDirection/Y"]
        del hf2["TracerData/ShockDirection/Z"]
        hf2.create_dataset("TracerData/ShockDirection/X", data=shock_direction_x, dtype=np.float32)
        hf2.create_dataset("TracerData/ShockDirection/Y", data=shock_direction_y, dtype=np.float32)
        hf2.create_dataset("TracerData/ShockDirection/Z", data=shock_direction_z, dtype=np.float32)

        if(len(mag_obliquity) != 0):
            del hf2["TracerData/MagneticObliquity"]
            hf2.create_dataset("TracerData/MagneticObliquity", data=mag_obliquity, dtype=np.float32)
        if(len(mach_number) != 0):
            del hf2["TracerData/MachNumber"]
            hf2.create_dataset("TracerData/MachNumber", data=mach_number, dtype=np.float32)

    hf2.close()

end_time = datetime.datetime.now()
print("End:", end_time)
print("Duration:", end_time - start_time)
