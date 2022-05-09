import numpy as np
import struct
import sys
from Physics import PROTONMASS
from os.path import isfile
from matplotlib.cbook import flatten
import h5py

####################################################################################################
# class which handles the parameters needed for calcualting the steady state solutions
# needed because these factors were phased out of Crest output
class cre_inj:
	def __init__(self):
		self.injRate = 0
		self.alphaInj = 0
		self.pInj = 0

####################################################################################################
# class which handles the parameters given via files to the C program
# parameters are needed for calculating the plots
class CrestParameters:
	"""Read the parameter file for a CREST simulation"""

	def __init__(self, file_name = None, verbose = False):
		"""
		Initialize all possible variables which can be defined by the parameter file.
		Note that it contains some legacy variables.
		If a path to the parameter file is given the file will be read in
		"""

		# General Settings for I/O
		self.OutputDir                   = ''
		self.InputDataFile               = ''
		self.InputDataFileBase           = ''
		self.InputTracerDataFileBase	 = ''
		self.InputArepoDataFileBase		 = ''
		self.InputFileFirstNum			 = 0
		self.NumInputFiles				 = 0
		self.MaxStepsReadIn				 = 0
		self.SnapshotFileBase            = ''
		self.OutputFileBase            	 = ''
		self.OutputVerbose				 = 0

		# Settings for Discretization
		self.NumberOfMomentumBins        = 0
		self.CourantFac                  = 0.0
		self.AlphaCoefficientMaximum     = 0.0
		self.MaximumSubcycles_in_log2    = 0
		self.MinimumMomentum             = 0.0
		self.MaximumMomentum             = 0.0
		self.TimeStepsUpdate             = 1
		self.IncludeMaximumMomentum      = 0

		# CR density distribution function
		self.AlphaSpectralIndex          = 0.0
		self.MomentumLowCutoff           = 0.0
		self.MomentumHighCutoff          = 0.0
		self.NormalizationFactor         = 0.0
		self.InitialSpectrumFile         = ''
		self.UseInitialSpectrumFile      = 0

		# System of units (does not apply for electron energy/momentum)
		self.UnitLength_in_cm            = 1.
		self.UnitMass_in_g               = 1.
		self.UnitVelocity_in_cm_per_s    = 1.

		# Output Settings
		self.OutputEverySnapshotOn       = 0
		self.TimeOfFirstSnapshot         = 0.0
		self.TimeBetSnapshot             = 0.01
		self.ComovingIntegrationOn   	 = 0

		# Flags
		self.FlagAllowSubcycles          = 1
		self.FlagCooling                 = 1
		self.Flag_Fermi_I_Reacceleration = 1
		self.Flag_Fermi_I_Acceleration   = 1
		self.Flag_Fermi_II_Reacceleration= 1
		self.FlagExternalInjection       = 0

		# Cooling & Diffusion
		self.n_elec                      = 1.157
		self.HydrogenMassFrac            = 0.76
		self.DiffusionTimeInGyr          = 0.
		self.Lambda                      = 0.
		self.Radiation_Field_in_eps_CMB  = -1. # if this value is -1 then it was not set
		self.Magnetic_Field_Amplification = -1. # if this value is -1 then it was not set
		self.Amplification_Flag			  = -1. # if this value is -1 then it was not set
		self.MeanMolecularWeight          = 1.

		# parameters for shock injection
		self.ShockParamA                 = 0.
		self.ShockParamB                 = 0.
		self.ShockParamC                 = 0.
		self.Acceleration_Max_Momentum   = 1.e20
		self.zeta_pe                     = 0.
		self.obliquity_critAngle         = -1.
		self.obliquity_width             = -1.
		self.obliquity_minEff            = -1.
		self.obliquity_maxEff            = -1.

		# new parameters with tracer data
		self.FunctionValueChop           = 1.e-30

		# Semi analytic treatement
		self.UseSemiAnalyticSolutionLow  = 1
		self.UseSemiAnalyticSolutionHigh = 1
		self.FlagFixedNumericsBoundaries = 1
		self.MaximumMomentumNumerics 	 = 0.
		self.MinimumMomentumNumerics 	 = 0.

		if file_name is not None:
			self.read_data(file_name, verbose)

	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def show(self):
		for var in vars(self):
			if type(getattr(self, var)) is int:
				print("{:25} {:d}".format(var,getattr(self,var)))
			if type(getattr(self, var)) is float:
				print("{:25} {:.5e}".format(var,getattr(self,var)))
			if type(getattr(self, var)) is str:
				print("{:25} {:}".format(var,getattr(self,var)))

	# read in the parameter file and set the private class variables accordingly
	def read_data(self, file_name, verbose = False):
		fParam = open(file_name, 'r')
		if verbose:
			print("Reading parameters from file '{:}'\n".format(file_name))
		for line in fParam:
			lineParam = (line.strip()).lstrip()

			# ignore lines beginning with a '%' sign
			if(lineParam != ''):
				if(lineParam[0] != '%'):
					columnParam = line.split()
					# only take lines which are of the following format: ParameterTag Space Value
					if(len(columnParam) == 2):
						# loop over all variables to find the one corresponding to the one read from the parameter fiel
						for var in vars(self):
							if(var == columnParam[0]):
								if type(getattr(self, var)) is int:
									setattr(self,var,int(columnParam[1]))
									if verbose:
										print("\t{:25} {:}".format(columnParam[0],columnParam[1]))
									continue
								elif type(getattr(self, var)) is float:
									setattr(self,var,float(columnParam[1]))
									if verbose:
										print("\t{:25} {:}".format(columnParam[0],columnParam[1]))
									continue
								elif type(getattr(self, var)) is str:
									setattr(self,var,columnParam[1])
									if verbose:
										print("\t{:25} {:}".format(columnParam[0],columnParam[1]))
									continue
		if self.OutputDir[-1] != '/':
			self.OutputDir += '/'
		if verbose:
			print("\n")
		line = None
		lineParam = None
		columnParam = None
		fParam.close()

	def BinsPerDec(self):
		return (self.NumberOfMomentumBins - self.IncludeMaximumMomentum) / int(np.log10(self.MaximumMomentum) - np.log10(self.MinimumMomentum))



####################################################################################################
def check_encoding():
	""" Check the size of integers and floats and give back their size in bytes. """

	error = [0, 0, 0, 0]

	# check integers
	size_i = struct.calcsize('i')
	if size_i != 4:
		error[0] = 1

	# check integers
	size_I = struct.calcsize('I')
	if size_i != 4:
		error[1] = 1

	# check single precision floats
	size_f = struct.calcsize('f')
	if size_f !=4:
		error[2] = 1

	# check double precision floats
	size_d = struct.calcsize('d')
	if size_d !=8:
		error[3] = 1

	if sum(error) > 0 :
		sys.exit("Data types ({}{}{}{}{}{}{}) not correctly encoded on this machine!".format(
			["", "int"][error[0]==1],
			["", ", "][error[0]==1 and (error[1]==1 or (error[2]==1 or error[3]==1))],
			["", "unsigned int"][error[1]==1],
			["", ", "][error[1]==1 and (error[2]==1 or error[3]==1)],
			["", "float"][error[2]==1],
			["", ", "][error[2] == 1 and error[3] == 1],
			["", "double"][error[3]==1]))

	else:
		return size_i, size_I, size_f, size_d


####################################################################################################
class CrestSnapshot:
	""" class for spectral snapshots of CREST """

	def __init__(self, file_name = None, verbose = False, get_only_header = False, specific_fields=None, use_HDF5=True):
		"""
		Initialize an instance of CREST snapshot.

		If no parameters are given an empty instance is initialized.
		If a path to a file is provided the file will be read in

		Args:
		   file_name (str): Path of file to be read

		   verbose (bool): Print more information on disply

		   get_only_header (bool): Read only header

		   specific_fields (list): List of strings of the variable which should be stored,
		                           e.g., ['ID', 'time']. Of no list is given, all variables are read.
		"""
		self._var_name = None
		self._var_dtype = None
		self._var_store = None

		self._use_hdf5 = use_HDF5		# By default use new HDF5 format; set = 0 to use original binary Arepo output instead

		if file_name is not None:
			self.read_data(file_name, verbose=verbose, get_only_header=get_only_header, specific_fields=specific_fields)

	def read_data(self, file_name, verbose = False, get_only_header = False, specific_fields=None):
		size_i, size_I, size_f, size_d = check_encoding()
		with open(file_name,'rb') as f:
			if(verbose):
				print("Reading snapshot data from file '{:}'".format(file_name))


			# Header information
			blocksize = int(struct.unpack('I', f.read(size_i))[0])
			headersize = 3 * size_i + size_d

			self.version = int(  struct.unpack('i', f.read(size_i))[0])
			if(verbose):
				print("Version {:d}-{:02d}".format(self.version//100, self.version%100))

			if self.version <= 201902:
				headersize = 3 * size_i + size_d
			else: # version >= 201903
				headersize = 4 * size_i + size_d

			if headersize != blocksize:
				sys.exit("Size of header is {:d} bytes, but expexted {:d}".format(blocksize, headersize))

			self.time    = float(struct.unpack('d', f.read(size_d))[0])
			self.nPart   = int(  struct.unpack('i', f.read(size_i))[0])
			self.nBins   = int(  struct.unpack('i', f.read(size_i))[0])

			if self.version >= 201903:
				self.flag_shock_acceleration = int(  struct.unpack('i', f.read(size_i))[0])
			else:
				self.flag_shock_acceleration = 1 # Default for older versions

			if int(struct.unpack('I',f.read(size_i))[0]) != blocksize:
				sys.exit("header data block not correctly enclosed")

			# Momentum Bins
			blocksize = int(struct.unpack('I', f.read(size_i))[0])
			momentumsize = self.nBins * size_d
			if momentumsize != blocksize:
				sys.exit("Block size is {:d} bytes, but expexted {:d}".format(blocksize, momentumsize))
			self.p = np.ndarray(self.nBins, dtype=float)
			self.p[:] = struct.unpack('{:d}d'.format(self.nBins), f.read(size_d * self.nBins))[:]

			if  int(struct.unpack('I',f.read(size_i))[0]) != blocksize:
				sys.exit("2nd data block not correctly enclosed")

			# Data
			if not get_only_header:

				# # Version specific setups
				# if self.version == 201903:
				# 	self._var_name = ['id','mass', 'n_gas', 'u_therm',\
				# 					  'eps_photon', 'B', 'pos', 'f']

				# 	# types of the variable
				# 	# we assume the types 'np.array' to be of lenght 3 and dtype np.float32
				# 	self._var_dtype = [np.uint32,  np.float64, np.float64, np.float64,\
				# 					   np.float64, np.float64, np.float64, np.float64]

				# 	# For scalars length is 1, otherwise it is the length of the array
				# 	self._var_length = [1, 1, 1, 1,\
				# 						1, 3, 3, self.nBins]


				# elif self.version == 201902 or self.version == 201901:
				# 	self._var_name = ['id', 'parent_cell_id', 'mass', 'n_gas',\
				# 					  'u_therm', 'eps_photon', 'B', 'pos', 'f']

				# 	# types of the variable
				# 	# we assume the types 'np.array' to be of lenght 3 and dtype np.float32
				# 	self._var_dtype = [np.uint32,  np.unit32, np.float64,  np.float64, \
				# 					   np.float64, np.float64, np.float64, np.float64, np.float64]

				# 	# For scalars length is 1, otherwise it is the length of the array
				# 	self._var_length = [1, 1, 1, 1,\
				# 						1, 1, 3, 3, self.nBins]

				# else:
				# 	sys.exit("Version {:d} not supported".format(self.version))

				# if len(self._var_name) != len(self._var_dtype) or len(self._var_name) != len(self._var_dtype):
				# 	sys.exit("Arrays of variable names and types need to have the same length")

				# # will the variable be stored or not
				# self._var_store = np.ones(len(self._var_name), dtype=bool) # array filled with True

				# if specific_fields is not None:
				# # make no storage default, invert var_store array
				# 	np.logical_not(self._var_store, out=self._var_store)

				# 	for sf in specific_fields:
				# 		if sf in self._var_name:
				# 			self._var_store[ self._var_name.index(sf) ] = True
				# 		else:
				# 			sys.exit("Variable '{:}' does not exist in tracer output!".format(sf))

				# for i in np.arange(len(self._var_name)):
				# 	# Create nPart x nSnap size arrays for all variables if var_store element is true
				# 	if self._var_store[i] and self._var_length[i] > 1:
				# 		setattr(self, self._var_name[i], np.ndarray((self.nPart, self._var_length), dtype=self._var_dtype[i]))
				# 	elif self._var_store[i] and self._var_length[i] == 1:
				# 		setattr(self, self._var_name[i], np.ndarray(self.nPart, dtype=self._var_dtype[i]))

				# 	# loop over all possible variables
				# 	for i in np.arange(len(self._var_name)):

				# 		# if variable should be stored, check the type and read from file in the correct format
				# 		if self._var_store[i]:

				# 			# Array like variables with dimensions stored separately
				# 			if self._var_length[i] > 1 and self._var_name[i] != 'f':
				# 				for j in np.arange(self._var_length[i]): # iterate over dimensions
				# 					if self._var_dtype[i] == np.uint32:
				# 						getattr(self, self._var_name[i])[:, j] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				# 					elif self._var_dtype[i] == np.int32:
				# 						getattr(self, self._var_name[i])[:, j] = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
				# 					elif self._var_dtype[i] == np.float32:
				# 						getattr(self, self._var_name[i])[:, j] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				# 					elif self._var_dtype[i] == np.float64:
				# 						getattr(self, self._var_name[i])[:, j] = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				# 					else:
				# 						sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))

				# 			# Spectrum where entire array for one particle is stored at once
				# 			elif self._var_length[i] > 1 and self._var_name[i] == 'f':
				# 				for k in np.arange(self.nPart): #iterate over particles
				# 					if self._var_dtype[i] == np.uint32:
				# 						getattr(self, self._var_name[i])[k, :] = struct.unpack('{:d}I'.format(self.nBins), f.read(size_I * self.nBins))
				# 					elif self._var_dtype[i] == np.int32:
				# 						getattr(self, self._var_name[i])[k, :] = struct.unpack('{:d}i'.format(self.nBins), f.read(size_i * self.nBins))
				# 					elif self._var_dtype[i] == np.float32:
				# 						getattr(self, self._var_name[i])[k, :] = struct.unpack('{:d}f'.format(self.nBins), f.read(size_f * self.nBins))
				# 					elif self._var_dtype[i] == np.float64:
				# 						getattr(self, self._var_name[i])[k, :] = struct.unpack('{:d}d'.format(self.nBins), f.read(size_d * self.nBins))
				# 					else:
				# 						sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))

				# 			elif self._var_length[i] == 1:
				# 				if self._var_dtype[i] == np.uint32:
				# 					getattr(self, self._var_name[i])[:] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				# 				elif self._var_dtype[i] == np.int32:
				# 					getattr(self, self._var_name[i])[:] = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
				# 				elif self._var_dtype[i] == np.float32:
				# 					getattr(self, self._var_name[i])[:] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				# 				elif self._var_dtype[i] == np.float64:
				# 					getattr(self, self._var_name[i])[:] = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				# 				else:
				# 					sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))
				# 			else:
				# 				sys.exit("Something is not correctly defined for variable '{:}' with type '{:}' and length '{:d}'".format(self._var_name[i], self._var_dtype, self._var_length))


				# 			if self._var_dtype[i] == np.uint32:
				# 				getattr(self, self._var_name[i])[n, :] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart)) # equivalent to e.g. self.ID[n, :] = struct.unpack( ... )
				# 			elif self._var_dtype[i] == np.int32:
				# 				getattr(self, self._var_name[i])[n, :] = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
				# 			elif self._var_dtype[i] == np.float32:
				# 				getattr(self, self._var_name[i])[n, :] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				# 			elif self._var_dtype[i] == np.float64:
				# 				getattr(self, self._var_name[i])[n, :] = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				# 			elif self._var_dtype[i] == np.ndarray:
				# 				for j in np.arange(3):
				# 					getattr(self, self._var_name[i])[n, :, j] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				# 			else:
				# 				sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))
				# 		else:
				# 			#if variable should not be stored, skip right number of bytes in file
				# 			if self._var_dtype[i] == np.uint32:
				# 				f.seek(size_I * self.nPart * self._var_length, 1)
				# 			elif self._var_dtype[i] == np.int32:
				# 				f.seek(size_i * self.nPart * self._var_length, 1)
				# 			elif self._var_dtype[i] == np.float32:
				# 				f.seek(size_f * self.nPart * self._var_length, 1)
				# 			elif self._var_dtype[i] == np.float64:
				# 				f.seek(size_d * self.nPart * self._var_length, 1)
				# 			else:
				# 				sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))

				# Spectrum Data
				blocksize = int(struct.unpack('I', f.read(size_i))[0])
				if self.version >= 201903:
					datasize = self.nPart * ( self.nBins * size_d + 1 * size_I + 10 * size_d)
				elif self.version == 201902:
					datasize = self.nPart * ( self.nBins * size_d + 2 * size_I + 10 * size_d)
				elif self.version == 201901:
					datasize = self.nPart * ( self.nBins * size_d + 2 * size_I + 8 * size_d)
				else:
					sys.exit("Version {:d} not supported".format(self.version))

				if blocksize != datasize:
					sys.exit("Block size is {:d} bytes, but expexted {:d}".format(blocksize, datasize))

				self.f = np.ndarray((self.nPart, self.nBins), dtype = float)
				self.id = np.ndarray(self.nPart, dtype=np.uint32)
				if self.version <= 201902:
					self.parent_cell_id = np.ndarray(self.nPart, dtype=np.uint32)

				if self._use_hdf5 == False:
					self.mass = np.ndarray(self.nPart, dtype=float)

				self.n_gas = np.ndarray(self.nPart, dtype=float)
				self.u_therm = np.ndarray(self.nPart, dtype=float)
				self.eps_photon = np.ndarray(self.nPart, dtype=float)
				self.pos = np.ndarray((self.nPart, 3), dtype=float)

				if self.version>=201902:
					self.B = np.ndarray((self.nPart, 3), dtype=float)

				self.id[:]             = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				if self.version <= 201902:
					self.parent_cell_id[:] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))

				if self._use_hdf5 == False:
					self.mass[:]           = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))

				self.n_gas[:]          = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.u_therm[:]        = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.eps_photon[:]     = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))

				if self.version>=201902:
					for j in np.arange(3):
						self.B[:, j]   = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				else:
					self.B[:]          = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))

				for j in np.arange(3):
					self.pos[:, j]     = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))

				for i in np.arange(self.nPart):
					self.f[i, :]      = struct.unpack('{:d}d'.format(self.nBins), f.read(size_d * self.nBins))

				blocksize_end = int(struct.unpack('I',f.read(size_i))[0])
				if blocksize_end != blocksize:
					sys.exit("3rd data block not correctly enclosed")



####################################################################################################

class ArepoTracerOutput:
	"""
	Class to read in and handle the output for the tracer particles from an AREPO simulation.
	The header file will be read in to provide basic information about the used unit system, as
	all quantities will be in code units of the AREPO simulation.
	If numbers of data files are given, data is read in.
	All variables will be stored as (nPart, nSnap) arrays unless only one specific particle
	is read in then the arrays have the shape (nSnap, ).

	Alternatively, a new empty tracer instance can be created, variables can be assigned, and
	data can be stored in file.

	Example:
	   Load a tracer output with path and name to file by, e.g.,

	      $ data = ArepoTracerOutput("path_to_file/fileBase")

	   and access all x positions

	      $ data.x

	   The data structure can be sliced or single particles/snapshots be picked by

	      $ data[1, 2] # particle 1 and snapshot 2
	      $ data[:, 1] # all particles and snapshot 1

	   Create new empty tracer data for 10 snapshots and 20 particles

	     $ data2 = ArepoTracerOutput()
	     $ data2.

	"""


	# instance variables
	def __init__(self, file_base = None, file_numbers=None, version=None, cgs_units = False, verbose = False, read_only_ic= False, specific_particles=None, first_snap=None, last_snap=None, specific_fields=None, splitted_files=True, use_HDF5=True, reshape_output=True):
		"""
		Initialize an instance of ArepoTracerOutput.

		If no parameters are given an empty instance is initialized.
		If a path to a file is provided the file will be read in

		Args:
		   file_base (str): Base file name (version >= 2020-01) or full file name (version <= 2019-03)

		   file_numbers (int or list): File numbers of data files to be read (version >= 2020-01)
		       Default: None (No data files to be read)

		   cgs_units (bool): Flag if the values should be converted to cgs units immediately

		   read_only_ic (bool): Read only header and 0th snapshot/initial conditions

		   first_snap (int): First snapshot to be read in (relative to first file number)

		   last_snap (int): Last snapshot to be read in (exclusive, relative to first file number)

		   specific_fields (list): List of strings of the variable which should be stored,
		      e.g., ['ID', 'time']. Of no list is given, all variables are read.
		      Possible variables names are:
		       - standard: ['time', 'pos', 'B',  'n_gas', 'u_therm', 'eps_photon']
		       - shock acceleration: ['ShockFlag', 'eps_CRp_acc',
		                              'n_gasPreShock', 'n_gasPostShock',
		                              'VShock', 'timeShockCross', 'ShockDir']
		       - magnetic obliquity: ['theta']
		       - SN injection: ['eps_CRp_inj']
		      Please note that some variable blocks are only available if the code was compiled
		      and run with these configurations.

		    splitted_files (bool): Separate header and data files (default, version >= 2020-01).
		       Chose 'False' if file of version <= 2019-03. Default: True

		"""

		# with_cr_electrons is set to 1 if arepo was compiled with #COSMIC_RAYS_ELECTRONS
		# need to set dummy values as these determine the types

		self.nSnap = 0
		self.nPart = 0
		self.All_Units_in_cgs = False
		self.UnitLength_in_cm = 1.
		self.UnitMass_in_g = 1.
		self.UnitVelocity_in_cm_per_s = 1.
		self._var_name = None
		self._var_dtype = None
		self._var_store = None
		self._var_cgs_factors = None
		self._version = None # current default version
		self._traceroutput_tracersize = None
		self._traceroutput_headersize = None
		self._use_hdf5 = use_HDF5			# By default use new HDF5 format; set = 0 to use original binary Arepo output instead

		if self._use_hdf5 and file_base is not None:
			self.read_header_hdf5(file_base, verbose=verbose)

			self.read_data_hdf5(file_base, reshape_output=reshape_output, file_numbers=file_numbers, cgs_units=cgs_units, verbose=verbose)

		elif file_base is not None:
			self.read_header(file_base, verbose=verbose, splitted_files=splitted_files)

			self.read_data(file_base, file_numbers=file_numbers, cgs_units=cgs_units, verbose=verbose, read_only_ic=read_only_ic,
						   specific_particles=specific_particles, first_snap=first_snap, last_snap=last_snap, specific_fields=specific_fields)

	@property
	def var_name(self):
		""" List of names of the variables which are in the file. Variables are stored if corresponding value of var_store is True. """
		return self._var_name

	@property
	def version(self):
		""" List of names of the variables which are in the file. Variables are stored if corresponding value of var_store is True. """
		return self._version

	def set_version(self, version):
		""" Manually set the version number"""
		self._version = version

	@property
	def var_dtype(self):
		""" List of data types of the variables """
		return self._var_dtype

	@property
	def var_store(self):
		""" List of flags to store of the variables """
		return self._var_store

	@property
	def var_cgs_factor(self):
		""" List of cgs conversion factors """
		return self._var_cgs_factors

	def define_variables(self, new=False):
		""" Define names, types, and unit conversion factor

		Args:

		new (bool) : Flag whether we create new tracer data instead of reading from a file
		    default False.

        """

		size_i, size_I, size_f, size_d = check_encoding()

		def variable_size(type_arr):
			def typesize(type_str):
				if type_str == np.int32:
					return size_i
				elif type_str == np.uint32:
					return size_I
				elif type_str == np.float32:
					return size_f
				elif type_str == np.float64:
					return size_d
				elif type_str == np.ndarray:
					return 3*size_f
				else:
					sys.exit("error")
			return np.sum(np.array([typesize(type_str) for type_str in type_arr]))

		L = self.UnitLength_in_cm
		M = self.UnitMass_in_g
		V = self.UnitVelocity_in_cm_per_s

		# cgs scaling of variables
		B = np.sqrt(M * V**2 / L**3) # B scaling
		N = M / (PROTONMASS * L**3) # Gas Density Scaling without mean molecular weight!
		E = M * V**2 / L**3 # Energy Density Scaling


		# names of all possible variables
		# for better readability: new line after 5 elements
		# the order of the variables has to be same as in the file


		if self._version == 201901:
			self._var_name = ['ID', 'time', 'ParentCellID', 'TracerMass', 'pos',\
							  'n_gas', 'u_therm', 'B', 'eps_photon', 'ShockFlag',\
							  'eps_CRp_acc', 'n_gasPreShock', 'n_gasPostShock', 'VShock', 'timeShockCross',\
							  'theta', 'CReInjection', 'injRate', 'alphaInj', 'pInj']\


			# types of the variable
			# we assume the types 'np.array' to be of lenght 3 and dtype np.float32
			self._var_dtype = [np.uint32,  np.float64, np.uint32,  np.float32, np.ndarray,\
							   np.float32, np.float32, np.float32, np.float32, np.int32,\
							   np.float32, np.float32, np.float32, np.float32, np.float32,\
							   np.float32, np.int32,   np.float32, np.float32, np.float32]


			if not new and variable_size(self._var_dtype) != self._traceroutput_tracersize:
				sys.exit("Size of tracer data block given in file is not the same as in code!")
			else:
				self._traceroutput_tracersize = variable_size(self._var_dtype)

			# cgs scaling of the variable
			self._var_cgs_factor = [1, L/V, 1, M, L,\
									N, V**2, B, E, 1,\
									E, N, N, V, L/V,\
									1, 1, V/L, 1, 1]


		elif self._version == 201902:
			self._var_name = ['ID', 'time', 'ParentCellID', 'TracerMass', 'pos',\
							  'n_gas', 'u_therm', 'B', 'eps_photon', 'ShockFlag',\
							  'eps_CRp_acc', 'n_gasPreShock', 'n_gasPostShock', 'VShock', 'timeShockCross',\
							  'ShockDir',  'theta', 'CReInjection', 'injRate', 'alphaInj',\
							  'pInj']

			# types of the variable
			self._var_dtype = [np.uint32,  np.float64, np.uint32,  np.float32, np.ndarray,\
							   np.float32, np.float32, np.ndarray, np.float32, np.int32,\
							   np.float32, np.float32, np.float32, np.float32, np.float32,\
							   np.ndarray, np.float32, np.int32,   np.float32, np.float32,\
							   np.float32]


			if not new and variable_size(self._var_dtype) != self._traceroutput_tracersize:
				sys.exit("Size of tracer data block given in file is not the same as in code!")
			else:
				self._traceroutput_tracersize = variable_size(self._var_dtype)

			# cgs scaling of the variable
			self._var_cgs_factor = [1, L/V, 1, M, L,\
									N, V**2, B, E, 1,\
									E, N, N, V, L/V,\
									1, 1, 1, V/L, 1,\
									1]

		elif self._version >= 201903:
			self._var_name = ['time', 'pos', 'B', 'n_gas', 'u_therm', 'eps_photon']
			self._var_dtype = [np.float64, np.ndarray, np.ndarray, np.float32, np.float32, np.float32]
			self._var_cgs_factor = [L/V, L, B, N, V**2, E]

			if self._use_hdf5:
				if self.flag_comoving_integration_on:
					self._var_name.append('dtValues')
					self._var_dtype.append(np.float64)
					self._var_cgs_factor.append(L/V)

					self._var_cgs_factor[0] = 1		# This value is now the scale factor

			if self.flag_cosmic_ray_shock_acceleration:
				self._var_name.append(['ShockFlag', 'eps_CRp_acc', 'n_gasPreShock',
										   'n_gasPostShock', 'VShock', 'timeShockCross', 'ShockDir'])

				self._var_dtype.append([np.int32, np.float32, np.float32,
										   np.float32, np.float32, np.float32, np.ndarray])

				self._var_cgs_factor.append([1, E, N,
												 N, V, L/V,	1])

				if self.flag_cosmic_ray_magnetic_obliquity:
					  self._var_name.append('theta')
					  self._var_dtype.append(np.float32)
					  self._var_cgs_factor.append(1)

			if self.flag_cosmic_ray_sn_injection:
				self._var_name.append('eps_CRp_inj')
				self._var_dtype.append(np.float32)
				self._var_cgs_factor.append(E)

			self._var_name = list(flatten(self._var_name))
			self._var_dtype = list(flatten(self._var_dtype))
			self._var_cgs_factor = list(flatten(self._var_cgs_factor))

		else:
			sys.exit("Version '{:d}-{:02d}' not supported or implemented!".format(self._version // 100, self._version % 100))

		if new:
			self._traceroutput_tracersize = (variable_size(self._var_dtype) - size_d) * self.nPart + size_d


	def initialize_variables(self, specific_fields=None, specific_particles=None):
		""" Initialize are variable arrays.

		The function define_variables() has to be called before!

		Args:
		   specific_fields (list of strings): Variable names if only specific fields should be stored

		   specific_particles (int or list of ints): Particles to be read/stored

		"""


		# will the variable be stored or not
		self._var_store = np.ones(len(self._var_name), dtype=bool) # array filled with True

		if specific_fields is not None:
			# make no storage default, invert var_store array
			np.logical_not(self._var_store, out=self._var_store)

			if type(specific_fields) is str:
				specific_fields = [specific_fields]

			for sf in specific_fields:
				if sf in self._var_name:
					self._var_store[ self._var_name.index(sf) ] = True
				else:
					sys.exit("Variable '{:}' does not exist in tracer output!".format(sf))


		if type(specific_particles) is not int:
			for i in np.arange(len(self._var_name)):
				# Create nPart x nSnap size arrays for all variables if var_store element is true
				if self._var_store[i]:
					if self._var_name[i] == 'time' and self._version >= 201903:
						setattr(self, self._var_name[i], np.ndarray(self.nSnap, dtype=self._var_dtype[i]))
					elif self._var_dtype[i] is not np.ndarray: # scalar variable
						setattr(self, self._var_name[i], np.ndarray((self.nSnap, self.nPart), dtype=self._var_dtype[i]))
					else: # 3D Vector
						setattr(self, self._var_name[i], np.ndarray((self.nSnap, self.nPart, 3), dtype=np.float32))

		else:
			for i in np.arange(len(self._var_name)):
				# Create nPart x nSnap size arrays for all variables if var_store element is true
				if self._var_store[i]:
					if self._var_name[i] == 'time' and self._version >= 201903:
						setattr(self, self._var_name[i], np.ndarray(self.nSnap, dtype=self._var_dtype[i]))
					elif self._var_dtype[i] is not np.ndarray: # scalar
						setattr(self, self._var_name[i], np.ndarray((self.nSnap, 1), dtype=self._var_dtype[i]))
					else: # 3D Vector
						setattr(self, self._var_name[i], np.ndarray((self.nSnap, 1, 3), dtype=np.float32))


	def read_header_hdf5(self, file_base, verbose = False):

		file_name = "{:}_000.hdf5".format(file_base)

		hf = h5py.File(file_name, 'r')

		self._version = hf['Header'].attrs['TracerOutputVersion']
		self.flag_cosmic_ray_shock_acceleration = hf['Header'].attrs['CosmicRaysShockAccelerationFlag']
		self.flag_cosmic_ray_magnetic_obliquity = hf['Header'].attrs['CosmicRaysMagneticObliquityFlag']
		self.flag_cosmic_ray_sn_injection = hf['Header'].attrs['CosmicRaysSNInjectionFlag']
		self.flag_comoving_integration_on = hf['Header'].attrs['ComovingIntegrationOnFlag']

		self.AllIDs = hf['Header/AllTracerParticleIDs'][()]

		self.nPart = len(self.AllIDs)

		self.UnitLength_in_cm = hf['Header'].attrs['UnitLength_in_cm']
		self.UnitMass_in_g = hf['Header'].attrs['UnitMass_in_g']
		self.UnitVelocity_in_cm_per_s = hf['Header'].attrs['UnitVelocity_in_cm_per_s']

		if self.flag_comoving_integration_on:
			self.hubble_param = hf['Header'].attrs['HubbleParam']

		hf.close()

		self.define_variables()

		if verbose:
			print("Header was read successfully")


	def read_data_hdf5(self, file_base, reshape_output, file_numbers=None, cgs_units = False, verbose = False):

		self.initialize_variables()

		if file_numbers is not None:
			if type(file_numbers) is int:
			  file_numbers = [file_numbers]
			  file_names = ["{:}_{:03d}.hdf5".format(file_base, num) for num in file_numbers]
		else:
			file_names = ["{:}_{:03d}.hdf5".format(file_base, 0)]

		if verbose:
			print("Read Arepo's tracer output from file '{}'".format(file_name))

		snapRead = 0 # number of already read in snapshots
		for file_name, file_num in zip(file_names, np.arange(len(file_names))):
			if verbose:
				print("Opening data file #{:d}".format(file_num))

			hf = h5py.File(file_name, 'r')

			self.time = hf['TracerData/Time'][()]

			self.nSnap = self.time.shape[0]

			self.ID = hf['TracerData/ParticleIDs'][()]
			pos_x = hf['TracerData/Coordinates/X'][()]
			pos_y = hf['TracerData/Coordinates/Y'][()]
			pos_z = hf['TracerData/Coordinates/Z'][()]
			self.pos = np.array([pos_x.T, pos_y.T, pos_z.T]).T
			mag_x = hf['TracerData/MagneticField/X'][()]
			mag_y = hf['TracerData/MagneticField/Y'][()]
			mag_z = hf['TracerData/MagneticField/Z'][()]
			self.B = np.array([mag_x.T, mag_y.T, mag_z.T]).T
			self.n_gas = hf['TracerData/Density'][()]
			self.u_therm = hf['TracerData/InternalEnergy'][()]
			self.eps_photon = hf['TracerData/PhotonEnergyDensity'][()]

			if self.flag_cosmic_ray_shock_acceleration:
				self.ShockFlag = hf['TracerData/ShockFlag'][()]
				shock_x = hf['TracerData/ShockDirection/X'][()]
				shock_y = hf['TracerData/ShockDirection/Y'][()]
				shock_z = hf['TracerData/ShockDirection/Z'][()]
				self.ShockDir = np.array([shock_x.T, shock_y.T, shock_z.T]).T
				self.eps_CRp_acc = hf['TracerData/ShockDissipatedThermalEnergy'][()]
				self.n_gasPreShock = hf['TracerData/PreShockDensity'][()]
				self.n_gasPostShock = hf['TracerData/PostShockDensity'][()]
				self.VShock = hf['TracerData/ShockVelocity'][()]
				self.timeShockCross = hf['TracerData/ShockCrossingTime'][()]

				if self.flag_cosmic_ray_sn_injection:
					self.theta = hf['TracerData/MagneticObliquity'][()]

			if self.flag_cosmic_ray_sn_injection:
				self.eps_CRp_inj = hf['TracerData/InjectionEnergy'][()]

			if self.flag_comoving_integration_on:
				self.dtValues = hf['TracerData/dtValues'][()]

			hf.close()

			if(reshape_output):
				print("CRE_ANALYSIS: reshape_output=True")
				print("CRE_ANALYSIS: reshaping Arepo tracer output into shapes of (nSnap, nPart)\n")

				self.ID = self.ID.reshape(self.nSnap, self.nPart)
				self.pos = self.pos.reshape(self.nSnap, self.nPart, 3)
				self.B = self.B.reshape(self.nSnap, self.nPart, 3)
				self.n_gas = self.n_gas.reshape(self.nSnap, self.nPart)
				self.u_therm = self.u_therm.reshape(self.nSnap, self.nPart)
				self.eps_photon = self.eps_photon.reshape(self.nSnap, self.nPart)

				if self.flag_cosmic_ray_shock_acceleration:
					self.ShockFlag = self.ShockFlag.reshape(self.nSnap, self.nPart)
					self.ShockDir = self.ShockDir.reshape(self.nSnap, self.nPart, 3)
					self.eps_CRp_acc = self.eps_CRp_acc.reshape(self.nSnap, self.nPart)
					self.n_gasPreShock = self.n_gasPreShock.reshape(self.nSnap, self.nPart)
					self.n_gasPostShock = self.n_gasPostShock.reshape(self.nSnap, self.nPart)
					self.VShock = self.VShock.reshape(self.nSnap, self.nPart)
					self.timeShockCross = self.timeShockCross.reshape(self.nSnap, self.nPart)

					if self.flag_cosmic_ray_sn_injection:
						self.theta = self.theta.reshape(self.nSnap, self.nPart)

				if self.flag_cosmic_ray_sn_injection:
					self.eps_CRp_inj = self.eps_CRp_inj.reshape(self.nSnap, self.nPart)

				if self.flag_comoving_integration_on:
					self.dtValues = self.dtValues.reshape(self.nSnap)

		if verbose:
			print("Data was read successfully")

		if cgs_units:
			self.scale_to_cgs_units(verbose)


	def read_header(self, file_base, verbose = False, splitted_files=True):
		""" Read in the header data from file. This function is automatically called if class constructor is called with the file name"""

		if splitted_files: # Standard for version >= 2020-01
			file_name = "{:}_header.dat".format(file_base)
		else:
			file_name = file_base

		with open(file_name,'rb') as f:
			if verbose:
				print("Read only initial conditions: {:}".format(read_only_ic))
				print("Read Arepo's tracer output from file '{}'".format(file_name))
			size_i, size_I, size_f, size_d = check_encoding()

			def variable_size(type_arr):
				def typesize(type_str):
					if type_str == np.int32:
						return size_i
					elif type_str == np.uint32:
						return size_I
					elif type_str == np.float32:
						return size_f
					elif type_str == np.float64:
						return size_d
					elif type_str == np.ndarray:
						return 3*size_f
					else:
						sys.exit("error")
				return np.sum(np.array([typesize(type_str) for type_str in type_arr]))

			# Version dependend configuration
			# we need to make sure that old simulations can be read in
			# python tool for tracer output conversion??

			# definitions from arepo's Config.sh file
			# only from version 2019-03 stored in tracer particle file
			self.flag_cosmic_ray_shock_acceleration = True
			self.flag_cosmic_ray_magnetic_obliquity = True
			self.flag_cosmic_ray_sn_injection = True

			# Reading first block with unit system
			self._traceroutput_headersize = int(struct.unpack('i', f.read(size_i))[0])
			if self._traceroutput_headersize == 3 * size_d:

				print("Warning: Old tracer output, data equivalent to version '2019-01'")
				self.UnitLength_in_cm         = struct.unpack('d', f.read(size_d))[0]
				self.UnitMass_in_g            = struct.unpack('d', f.read(size_d))[0]
				self.UnitVelocity_in_cm_per_s = struct.unpack('d', f.read(size_d))[0]

				# Hard coded to guarantee readability of old tracer output files
				self._version = 201901
				f.seek(size_i, 1) # jump to beginning of next block to get extra information
				self._traceroutput_tracersize = 2 * size_i + 2 * size_I + 17 * size_f + 1 * size_d
				self.nPart = int(struct.unpack('i',f.read(size_i))[0]) // self._traceroutput_tracersize
				f.seek( -2 * size_i, 1) # jump back to previous position after the unit block

			elif self._traceroutput_headersize == 3 * size_d + 3 * size_i:
				# version 2019-02
				self._version                 = struct.unpack('i', f.read(size_i))[0]
				self.UnitLength_in_cm         = struct.unpack('d', f.read(size_d))[0]
				self.UnitMass_in_g            = struct.unpack('d', f.read(size_d))[0]
				self.UnitVelocity_in_cm_per_s = struct.unpack('d', f.read(size_d))[0]
				self.nPart                    = struct.unpack('i', f.read(size_i))[0]
				self._traceroutput_tracersize = struct.unpack('i', f.read(size_i))[0]

			else:
				# version 2019-03 or newer
				self._version                 = struct.unpack('i', f.read(size_i))[0]
				self.UnitLength_in_cm         = struct.unpack('d', f.read(size_d))[0]
				self.UnitMass_in_g            = struct.unpack('d', f.read(size_d))[0]
				self.UnitVelocity_in_cm_per_s = struct.unpack('d', f.read(size_d))[0]
				self.nPart                    = struct.unpack('i', f.read(size_i))[0]
				self.flag_cosmic_ray_shock_acceleration = bool(struct.unpack('i', f.read(size_i))[0])
				self.flag_cosmic_ray_magnetic_obliquity = bool(struct.unpack('i', f.read(size_i))[0])
				self.flag_cosmic_ray_sn_injection = bool(struct.unpack('i', f.read(size_i))[0])
				self.TracerMass = np.ndarray(self.nPart, dtype=np.float32)
				self.TracerMass[:] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))[:]
				self.ID = np.ndarray(self.nPart, dtype=np.uint32)
				self.ID[:] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))[:]
				self._traceroutput_tracersize = struct.unpack('i', f.read(size_i))[0]

			end_integer = int(struct.unpack('i', f.read(size_i))[0])
			if end_integer != self._traceroutput_headersize:
				sys.exit("Header block not correctly enclosed")

			self.define_variables()

			f.close()
			if verbose:
				print("Header was successfully read")

	def read_data(self, file_base, file_numbers=None, cgs_units = False, verbose = False, read_only_ic = False, specific_particles = None, first_snap = None, last_snap = None, specific_fields=None):
		""" Read in the header data from file. This function is automatically called if class constructor is called with the file name"""

		if self.version <= 201903:
			file_names = [file_base]
		else: # version >= 202001
			if file_numbers is not None:
				if type(file_numbers) is int:
					file_numbers = [file_numbers]

				file_names = ["{:}_file_{:03d}.dat".format(file_base, num) for num in file_numbers]
			else:
				file_names = []


		if verbose:
			print("Read Arepo's tracer output from file '{}'".format(file_names))

		size_i, size_I, size_f, size_d = check_encoding()


		# Version dependend configuration
		# we need to make sure that old simulations can be read in

		# Reading block with data values
		nPartInFile = self.nPart
		self.nSnap	= 0

		for file_name in file_names:
			with open(file_name,'rb') as f:
				if self._version <= 201903:
					# Jump over headerblock
					f.seek(2*size_i + self._traceroutput_headersize, 0)

				buf = f.read(size_i) # Read starting integer of first data block
				if buf: #if buf > 0 or True next block will be read
					blocksize = int(struct.unpack('i', buf)[0])


					# Jump over data blocks and check if starting and trailing integer are the same
					while(buf):
						# move pointer forward
						if self._version <= 201902:
							f.seek(self.nPart * self._traceroutput_tracersize, 1)
						else: # version >= 201903
							f.seek(self._traceroutput_tracersize, 1)

						if  int(struct.unpack('i',f.read(size_i))[0]) != blocksize:
							sys.exit("data not in block #{:d} not correctly enclosed".format(self.nSnap))

						self.nSnap += 1
						if read_only_ic or last_snap == self.nSnap:
							buf = False
							break # stop reading
						else:
							buf = f.read(size_i)
				else:
					print("No data block was found")

				f.close()

		if specific_particles is not None:
			if type(specific_particles) is int:
				self.nPart = 1
			else:
				self.nPart = len(np.arange(self.nPart)[specific_particles])

		if first_snap is not None:
			self.nSnap -= first_snap

		if verbose:
			print("Number of particles: {:d}".format(self.nPart))
			print("Number of snapshots: {:d}".format(self.nSnap))

		if len(self._var_name) != len(self._var_dtype) or len(self._var_name) != len(self._var_cgs_factor):
			sys.exit("Arrays of variable names and types need to have the same length. Check read_header function")


		self.initialize_variables(specific_fields=specific_fields, specific_particles=specific_particles)

		snapRead = 0 # number of already read in snapshots
		for file_name, file_num in zip(file_names, np.arange(len(file_names))):
			with open(file_name,'rb') as f:
				if verbose:
					print("Opening data file #{:d}".format(file_num))

				if self._version <= 201903:
					# jump over header block
					f.seek(2*size_i + self._traceroutput_headersize, 0)

				if file_num == 0 and first_snap is not None:
					# skip some lines
					blocksize = int(struct.unpack('i', f.read(size_i))[0])
					f.seek(first_snap * (blocksize + 2*size_i) - size_i, 1)

				buf = f.read(size_i) # Read starting integer of first data block
				if buf: #if buf > 0 or True next block will be read
					blocksize = int(struct.unpack('i', buf)[0])

					if specific_particles is None:
						# read all the data spec
						for n in np.arange(snapRead, self.nSnap):

							# loop over all possible variables
							for i in np.arange(len(self._var_name)):

								if self._var_store[i]:
									# if variable should be stored, check the type and read from file in the correct format
									if self._var_name[i] == 'time' and self._version >= 201903:
										self.time[n] = struct.unpack('d', f.read(size_d))[0] #We assume time to be of dtype=np.float64

									elif self._var_dtype[i] == np.uint32:
										getattr(self, self._var_name[i])[n, :] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart)) # equivalent to e.g. self.ID[n, :] = struct.unpack( ... )
									elif self._var_dtype[i] == np.int32:
										getattr(self, self._var_name[i])[n, :] = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
									elif self._var_dtype[i] == np.float32:
										getattr(self, self._var_name[i])[n, :] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
									elif self._var_dtype[i] == np.float64:
										getattr(self, self._var_name[i])[n, :] = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
									elif self._var_dtype[i] == np.ndarray:
										for j in np.arange(3):
											getattr(self, self._var_name[i])[n, :, j] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
									else:
										sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))
								else:
									#if variable should not be stored, skip right number of bytes in file
									if self._var_name[i] == 'time' and self._version >= 201903:
										f.seek(size_d, 1) # We assume time to be of dtype=np.float64
									elif self._var_dtype[i] == np.uint32:
										f.seek(size_I * self.nPart, 1)
									elif self._var_dtype[i] == np.int32:
										f.seek(size_i * self.nPart, 1)
									elif self._var_dtype[i] == np.float32:
										f.seek(size_f * self.nPart, 1)
									elif self._var_dtype[i] == np.float64:
										f.seek(size_d * self.nPart, 1)
									elif self._var_dtype[i] == np.ndarray:
										f.seek(3 * size_f * self.nPart, 1)
									else:
										sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))

							if  int(struct.unpack('i', f.read(size_i))[0]) != blocksize:
								sys.exit("Data in block #{:d} not correctly enclosed.".format(n))

							# Read integer of next data block
							buf = f.read(size_i)
							if buf:
								blocksize_next = int(struct.unpack('i', buf)[0])
								if blocksize_next != blocksize:
									sys.exit("Starting integer of block #{:d} differs from block before.".format(n+1))
							else:
								break

					else: # specific particles are chosen to be read in
						if type(specific_particles) is np.ndarray or type(specific_particles) is list:
							if len(specific_particles) == 1:
								specific_particles = specific_particles[0]

						if type(specific_particles) is int:
							pos = specific_particles
							# read single data
							for n in np.arange(snapRead, self.nSnap):
								# loop over all possible variables
								for i in np.arange(len(self._var_name)):

									if self._var_store[i]:
										# if variable should be stored, check the type and read from file in the correct format
										# jump to the position of the desired particle in the variable block
										# read in the right amount of bytes
										# jump to end of the variable block
										if self._var_name[i] == 'time' and self._version >= 201903:
											self.time[n] = struct.unpack('d', f.read(size_d))[0] #We assume time to be of dtype=np.float64

										elif self._var_dtype[i] == np.uint32:
											f.seek(pos * size_I, 1)
											getattr(self, self._var_name[i])[n, :] = struct.unpack('I', f.read(size_I))[0]
											f.seek(size_I * (nPartInFile - pos - 1), 1)

										elif self._var_dtype[i] == np.int32:
											f.seek(pos * size_i, 1)
											getattr(self, self._var_name[i])[n, :] = struct.unpack('i', f.read(size_i))[0]
											f.seek(size_i * (nPartInFile - pos - 1), 1)

										elif self._var_dtype[i] == np.float32:
											f.seek(pos * size_f, 1)
											getattr(self, self._var_name[i])[n, :] = struct.unpack('f', f.read(size_f))[0]
											f.seek(size_f * (nPartInFile - pos - 1), 1)

										elif self._var_dtype[i] == np.float64:
											f.seek(pos * size_d, 1)
											getattr(self, self._var_name[i])[n, :] = struct.unpack('d', f.read(size_d))[0]
											f.seek(size_d * (nPartInFile - pos - 1), 1)

										elif self._var_dtype[i] == np.ndarray:
											for j in np.arange(3):
												f.seek(pos * size_f, 1)
												getattr(self, self._var_name[i])[n, :, j] = struct.unpack('f', f.read(size_f))[0]
												f.seek(size_f * (nPartInFile - pos - 1), 1)

										else:
											sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))
									else:
										#if variable should not be stored, skip right number of bytes in file
										if self._var_name[i] == 'time' and self._version >= 201903:
											f.seek(size_d, 1) #We assume time to be of dtype=np.float64
										elif self._var_dtype[i] == np.uint32:
											f.seek(size_I * nPartInFile, 1)
										elif self._var_dtype[i] == np.int32:
											f.seek(size_i * nPartInFile, 1)
										elif self._var_dtype[i] == np.float32:
											f.seek(size_f * nPartInFile, 1)
										elif self._var_dtype[i] == np.float64:
											f.seek(size_d * nPartInFile, 1)
										elif self._var_dtype[i] == np.ndarray:
											f.seek(3 * size_f * nPartInFile, 1)
										else:
											sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))

								if  int(struct.unpack('i', f.read(size_i))[0]) != blocksize:
									sys.exit("Data in block #{:d} not correctly enclosed.".format(n))

								f.seek(size_i, 1)


						else:
							# read all the data but just take a slice out of it
							for n in np.arange(snapRead, self.nSnap):

								# loop over all possible variables
								for i in np.arange(len(self._var_name)):

									if self._var_store[i]:
										# if variable should be stored, check the entire block of the variable
										# convert to numpy array which supports indexing with arrays, slices, etc.
										# (return type of struct.unpack is 'tuple')
										if self._var_dtype[i] == np.float64:
											self.time[n] = np.array(struct.unpack('d', f.read(size_d)), dtype=np.float64)[0] #We assume time to be of dtype=np.float64
										elif self._var_dtype[i] == np.uint32:
											getattr(self, self._var_name[i])[n, :] = np.array(struct.unpack('{:d}I'.format(nPartInFile), f.read(size_I * nPartInFile)), dtype=np.uint32)[specific_particles]
										elif self._var_dtype[i] == np.int32:
											getattr(self, self._var_name[i])[n, :] = np.array(struct.unpack('{:d}i'.format(nPartInFile), f.read(size_i * nPartInFile)), dtype=np.int32)[specific_particles]
										elif self._var_dtype[i] == np.float32:
											getattr(self, self._var_name[i])[n, :] = np.array(struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile)), dtype=np.float32)[specific_particles]
										elif self._var_dtype[i] == np.float64:
											getattr(self, self._var_name[i])[n, :] = np.array(struct.unpack('{:d}d'.format(nPartInFile), f.read(size_d * nPartInFile)), dtype=np.float64)[specific_particles]
										elif self._var_dtype[i] == np.ndarray:
											for j in np.arange(3):
												getattr(self, self._var_name[i])[n, :, j] = np.array(struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile)), dtype=np.float32)[specific_particles]
										else:
											sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))

									else:
										#if variable should not be stored, skip right number of bytes in file
										if self._var_name[i] == 'time' and self._version >= 201903:
											f.seek(size_d, 1) #We assume time to be of dtype=np.float64
										elif self._var_dtype[i] == np.uint32:
											f.seek(size_I * nPartInFile, 1)
										elif self._var_dtype[i] == np.int32:
											f.seek(size_i * nPartInFile, 1)
										elif self._var_dtype[i] == np.float32:
											f.seek(size_f * nPartInFile, 1)
										elif self._var_dtype[i] == np.float64:
											f.seek(size_d * nPartInFile, 1)
										elif self._var_dtype[i] == np.ndarray:
											f.seek(3 * size_f * nPartInFile, 1)
										else:
											sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))

								if  int(struct.unpack('i', f.read(size_i))[0]) != blocksize:
									sys.exit("Data in block #{:d} not correctly enclosed.".format(n))

								# Read integer of next data block
								buf = f.read(size_i)
								if buf:
									blocksize_next = int(struct.unpack('i', buf)[0])
									if blocksize_next != blocksize:
										sys.exit("Starting integer of block #{:d} differs from block before.".format(n+1))
								else:
									break

					snapRead += n + 1

				else:
					print("No data block was found. If first_snap larger than number of Snapshots in first file, chose smaller first_snap and larger first file number.")

				f.close()

		if verbose:
			print("Data was successfully read")

		if cgs_units:
			self.scale_to_cgs_units(verbose)


	def scale_to_cgs_units(self, verbose=False):
		if not self.All_Units_in_cgs:
			if verbose:
				print("Scale to cgs with UnitLength_in_cm = {:.3e}, UnitMass_in_g = {:.3e}, UnitVeloctiy_in_cm_per_s = {:.3e}".format(self.UnitLength_in_cm, self.UnitMass_in_g, self.UnitVelocity_in_cm_per_s))
			for i in np.arange(len(self._var_name)):
				if self._var_store[i]:
					setattr(self, self._var_name[i], np.multiply(self._var_cgs_factor[i], getattr(self,self._var_name[i])).astype(self._var_dtype[i]))

			self.All_Units_in_cgs = True

		else:
			print("Variables are already stored in cgs units")

	def scale_to_code_units(self, verbose=False):
		if self.All_Units_in_cgs:
			if verbose:
				print("Scale to code units with UnitLength_in_cm = {:.3e}, UnitMass_in_g = {:.3e}, UnitVeloctiy_in_cm_per_s = {:.3e}".format(self.UnitLength_in_cm, self.UnitMass_in_g, self.UnitVelocity_in_cm_per_s))
			for i in np.arange(len(self._var_name)):
				if self._var_store[i]:
					setattr(self, self._var_name[i], np.divide(self._var_cgs_factor[i], getattr(self,self._var_name[i])).astype(self._var_dtype[i]))

			self.All_Units_in_cgs = False

		else:
			print("Variables are already stored in code units")


	def create_new_data(self, nSnap, nPart, version=202201, UnitLength_in_cm = 1., UnitMass_in_g = 1., UnitVelocity_in_cm_per_s = 1.,
						flag_cosmic_ray_shock_acceleration = False, flag_cosmic_ray_magnetic_obliquity = False,
						flag_cosmic_ray_sn_injection = False, flag_comoving_integration_on = False, hubble_param = 1):
		""" Create new empty tracer data

		Args:
           nSnap (int) : Number of snapshots

		   nPart (int) : Number of tracer particles

		   UnitLength_in_cm (float) : default 1

		   UnitMass_in_g (float) : default 1

		   UnitVelocity_in_cm_per_s (float) : default 1

		   flag_cosmic_ray_sn_injection (bool) : default False

		   flag_cosmic_ray_magnetic_obliquity (bool) : default False

		   flag_cosmic_ray_sn_injection (bool) : default False

		   flag_comoving_integration_on (bool) : default False

		"""

		size_i, size_I, size_f, size_d = check_encoding()

		self._version = version
		self.nSnap = nSnap
		self.nPart = nPart

		self.flag_cosmic_ray_shock_acceleration = flag_cosmic_ray_shock_acceleration
		self.flag_cosmic_ray_magnetic_obliquity = flag_cosmic_ray_magnetic_obliquity
		self.flag_cosmic_ray_sn_injection = flag_cosmic_ray_sn_injection
		self.flag_comoving_integration_on = flag_comoving_integration_on

		self.UnitLength_in_cm = UnitLength_in_cm
		self.UnitMass_in_g = UnitMass_in_g
		self.UnitVelocity_in_cm_per_s = UnitVelocity_in_cm_per_s

		if self.flag_comoving_integration_on:
			self.hubble_param = hubble_param

		self.define_variables(new=True)
		self.initialize_variables()

		self.TracerMass = np.ndarray(self.nPart, dtype=np.float32)

		if self._use_hdf5:
			self.AllIDs = np.ndarray(self.nPart, dtype=np.uint32)
			self.ID = np.ndarray([self.nSnap,self.nPart], dtype=np.uint32)
		else:
			self.ID = np.ndarray(self.nPart, dtype=np.uint32)

		if self._version >= 201903:
			self._traceroutput_headersize = 6 * size_i + 3 * size_d + self.nPart * (size_f + size_I)


		else:
			sys.exit("Version {:}-{:02d} not implemented yet for writing")



	def __getitem__(self, key):
		# check dimensions of return
		ret = ArepoTracerOutput()
		ret.All_Units_in_cgs = self.All_Units_in_cgs
		ret.UnitLength_in_cm = self.UnitLength_in_cm
		ret.UnitMass_in_g = self.UnitMass_in_g
		ret.UnitVelocity_in_cm_per_s = self.UnitVelocity_in_cm_per_s
		ret._version = self._version
		ret.flag_cosmic_ray_shock_acceleration = self.flag_cosmic_ray_shock_acceleration
		ret.flag_cosmic_ray_magnetic_obliquity = self.flag_cosmic_ray_magnetic_obliquity
		ret.flag_cosmic_ray_sn_injection = self.flag_cosmic_ray_sn_injection


		ret._var_name = self._var_name
		ret._var_dtype = self._var_dtype
		ret._var_cgs_factor = self._var_cgs_factor
		ret._var_store = self._var_store

		if isinstance(key, int) or isinstance(key, np.int)\
				or isinstance(key, np.int32) or isinstance(key, np.int64):
			ret.nSnap = 1
			ret.nPart = self.nPart
		elif isinstance(key, slice):
			start, stop, step = key.indices(self.nSnap)
			ret.nSnap = np.arange(start, stop, step).size
			ret.nPart = self.nPart
		elif isinstance(key, tuple):
			if len(key) == 2:
				if isinstance(key[0], int) or isinstance(key[0], np.int)\
						or isinstance(key[0], np.int32) or isinstance(key[0], np.int64):
					ret.nSnap = 1
				elif isinstance(key[0], slice):
					start, stop, step = key[0].indices(self.nSnap)
					ret.nSnap = np.arange(start, stop, step).size
				else:
					raise TypeError('Index must be int or slice, not {}'.format(type(key[0]).__name__))

				if isinstance(key[1], int) or isinstance(key[1], np.int)\
						or isinstance(key[1], np.int32) or isinstance(key[1], np.int64):
					ret.nPart = 1
				elif isinstance(key[1], slice):
					start, stop, step = key[1].indices(self.nPart)
					ret.nPart = np.arange(start, stop, step).size
				else:
					raise TypeError('Index must be int or slice, not {}'.format(type(key[1]).__name__))
		else:
			raise TypeError('Tuple Index must be of length 2, not {}'.format(len(key)))


		for i in np.arange(len(self._var_name)):
			if self._var_store[i]:
				if self._var_name[i] == 'time' and self._version >= 201903:
					if len(key) == 2:
						ret.time = self.time.__getitem__(key[1])
				else:
					setattr(ret, ret._var_name[i], getattr(self, self._var_name[i]).__getitem__(key))

		return ret

	def save(self, file_name):
		"""
		Write a new output file with the subset of currently stored particles.

		Args:

		   file_name (str) : Base for file name (version >= 2020-01) or full file name (version <= 2019-03)

		"""

		if not self._var_store.all():
			print("Warning: No output was created! Not every possible field is defined, therefore an output is unreadable by CREST.")
			return None

		if self._use_hdf5:

			file_name_data = "{:}_000.hdf5".format(file_name)

			if isfile(file_name_data):
				print("Warning: File '{:}' already exists and will be overwritten".format(file_name_data))			# Should we ask for user input here? - will this ever need to run on a cluster?

			hf = h5py.File(file_name_data, 'w')

			# Store header block
			header = hf.create_group('Header')
			group_dat = hf.create_group('TracerData')
			group_pos = group_dat.create_group('Coordinates')
			group_mag = group_dat.create_group('MagneticField')

			if self.flag_cosmic_ray_shock_acceleration:
				group_shock = group_dat.create_group('ShockDirection')

			# Datatypes
			int = np.int32
			float = np.float32
			double =np.float64

			len_tracer_data = self.nSnap * self.nPart

			# Header info
			header.attrs['TracerOutputVersion'] = self._version
			header.attrs['CosmicRaysShockAccelerationFlag'] = int(self.flag_cosmic_ray_shock_acceleration)
			header.attrs['CosmicRaysMagneticObliquityFlag'] = int(self.flag_cosmic_ray_magnetic_obliquity)
			header.attrs['CosmicRaysSNInjectionFlag'] = int(self.flag_cosmic_ray_sn_injection)
			header.attrs['ComovingIntegrationOnFlag'] = int(self.flag_comoving_integration_on)
			header.create_dataset('AllTracerParticleIDs', data = np.unique(self.ID), dtype=int)

			header.attrs['UnitLength_in_cm'] = self.UnitLength_in_cm
			header.attrs['UnitMass_in_g'] = self.UnitMass_in_g
			header.attrs['UnitVelocity_in_cm_per_s'] = self.UnitVelocity_in_cm_per_s

			if self.flag_comoving_integration_on:
			    header.attrs['HubbleParam'] = self.hubble_param

			# Tracer data
			d1 = group_dat.create_dataset('ParticleIDs', (len_tracer_data) , dtype=int)

			d2 = group_pos.create_dataset('X', (len_tracer_data) , dtype=float)
			d3 = group_pos.create_dataset('Y', (len_tracer_data) , dtype=float)
			d4 = group_pos.create_dataset('Z', (len_tracer_data) , dtype=float)

			d5 = group_mag.create_dataset('X', (len_tracer_data,) , dtype=float)
			d6 = group_mag.create_dataset('Y', (len_tracer_data) , dtype=float)
			d7 = group_mag.create_dataset('Z', (len_tracer_data) , dtype=float)

			d8 = group_dat.create_dataset('Density', (len_tracer_data) , dtype=float)
			d9 = group_dat.create_dataset('InternalEnergy', (len_tracer_data) , dtype=float)
			d10 = group_dat.create_dataset('PhotonEnergyDensity', (len_tracer_data) , dtype=float)

			if self.flag_cosmic_ray_shock_acceleration:
			    d11 = group_dat.create_dataset('ShockFlag', (len_tracer_data) , dtype=int)

			    d12 = group_shock.create_dataset('X', (len_tracer_data) , dtype=float)
			    d13 = group_shock.create_dataset('Y', (len_tracer_data) , dtype=float)
			    d14 = group_shock.create_dataset('Z', (len_tracer_data) , dtype=float)

			    d15 = group_dat.create_dataset('ShockDissipatedThermalEnergy', (len_tracer_data) , dtype=float)
			    d16 = group_dat.create_dataset('PreShockDensity', (len_tracer_data) , dtype=float)
			    d17 = group_dat.create_dataset('PostShockDensity', (len_tracer_data) , dtype=float)
			    d18 = group_dat.create_dataset('ShockVelocity', (len_tracer_data) , dtype=float)
			    d19 = group_dat.create_dataset('ShockCrossingTime', (len_tracer_data) , dtype=float)

			    if self.flag_cosmic_ray_magnetic_obliquity:
			        d20 = group_dat.create_dataset('MagneticObliquity', (len_tracer_data) , dtype=float)

			if self.flag_cosmic_ray_sn_injection:
			    d21 = group_dat.create_dataset('InjectionEnergy', (len_tracer_data) , dtype=float)

			d22 = group_dat.create_dataset('Time', (self.nSnap) , dtype=double)

			if self.flag_comoving_integration_on:
			    d23 = group_dat.create_dataset('dtValues', (self.nSnap) , dtype=double)

			d24 = group_dat.create_dataset('TimestepLastIndex', (self.nSnap), dtype=int)

			last_index = 0

			for i in range(self.nSnap):
				d1[i] = self.ID[i].flatten()
				d2[i] = self.pos[i,:,0].flatten()
				d3[i] = self.pos[i,:,1].flatten()
				d4[i] = self.pos[i,:,2].flatten()
				d5[i] = self.B[i,:,0].flatten()
				d6[i] = self.B[i,:,1].flatten()
				d7[i] = self.B[i,:,2].flatten()
				d8[i] = self.n_gas[i].flatten()
				d9[i] = self.u_therm[i].flatten()
				d10[i] = self.eps_photon[i].flatten()

				if self.flag_cosmic_ray_shock_acceleration:
					d11[i] = self.ShockFlag[i].flatten()
					d12[i] = self.ShockDir[i,:,0].flatten()
					d13[i] = self.ShockDir[i,:,1].flatten()
					d14[i] = self.ShockDir[i,:,2].flatten()
					d15[i] = self.eps_CRp_acc[i].flatten()
					d16[i] = self.n_gasPreShock[i].flatten()
					d17[i] = self.n_gasPostShock[i].flatten()
					d18[i] = self.VShock[i].flatten()
					d19[i] = self.timeShockCross[i].flatten()
					if self.flag_cosmic_ray_magnetic_obliquity:
						d20[i] = self.theta[i].flatten()

				if self.flag_cosmic_ray_sn_injection:
					d21[i] = self.eps_CRp_inj[i].flatten()

				d22[i] = self.time[i].flatten()

				if self.flag_comoving_integration_on:
					d23[i] = self.dtValues[i].flatten()

				last_index += len(self.ID[i])
				d24[i] = last_index

			hf.close()

			print("Arepo tracer output data written to '{}'".format(file_name_data))


		elif self._version >= 202001:
			file_name_header = "{:}_header.dat".format(file_name)
			file_name_data = "{:}_file_000.dat".format(file_name)

			print("I create the files '{:}' and '{:}'".format(file_name_header, file_name_data))

			if isfile(file_name_header):
				print("Warning: File '{:}' existed and will be overwritten".format(file_name_header))
				f_op = 'wb'
			else:
				f_op = 'xb'

			with open(file_name_header, f_op) as f:
				# Store header block
				f.write(struct.pack('i', self._traceroutput_headersize))
				f.write(struct.pack('i', self._version))
				f.write(struct.pack('d', self.UnitLength_in_cm))
				f.write(struct.pack('d', self.UnitMass_in_g))
				f.write(struct.pack('d', self.UnitVelocity_in_cm_per_s))
				f.write(struct.pack('i', self.nPart))
				f.write(struct.pack('i', int(self.flag_cosmic_ray_shock_acceleration)))
				f.write(struct.pack('i', int(self.flag_cosmic_ray_magnetic_obliquity)))
				f.write(struct.pack('i', int(self.flag_cosmic_ray_sn_injection)))
				f.write(struct.pack('{:d}f'.format(self.nPart), *self.TracerMass.astype(np.float)[:]))
				f.write(struct.pack('{:d}I'.format(self.nPart), *self.ID.astype(np.uint32)[:]))
				f.write(struct.pack('i', self._traceroutput_tracersize))
				f.write(struct.pack('i', self._traceroutput_headersize))

				f.close()

			if isfile(file_name_data):
				print("Warning: File '{:}' existed and will be overwritten".format(file_name_data))
				f_op = 'wb'
			else:
				f_op = 'xb'

			with open(file_name_data, f_op) as f:
				for n in np.arange(self.nSnap):
					f.write(struct.pack('i', self._traceroutput_tracersize))
					for i in np.arange(len(self._var_name)):
						if self._var_name[i] == 'time' and self._version >= 201903:
							f.write(struct.pack('d', self.time[n].astype(self._var_dtype[i])))
						elif self._var_dtype[i] == np.uint32:
							f.write(struct.pack('{:d}I'.format(self.nPart), *getattr(self, self._var_name[i])[n, :].astype(self._var_dtype[i]))) # equivalent to e.g. struct.pack(... , self.ID[n, :] )
						elif self._var_dtype[i] == np.int32:
							f.write(struct.pack('{:d}i'.format(self.nPart), *getattr(self, self._var_name[i])[n, :].astype(self._var_dtype[i])))
						elif self._var_dtype[i] == np.float32:
							f.write(struct.pack('{:d}f'.format(self.nPart), *getattr(self, self._var_name[i])[n, :].astype(self._var_dtype[i])))
						elif self._var_dtype[i] == np.float64:
							f.write(struct.pack('{:d}d'.format(self.nPart), *getattr(self, self._var_name[i])[n, :].astype(self._var_dtype[i])))
						elif self._var_dtype[i] == np.ndarray:
							for j in np.arange(3):
								f.write(struct.pack('{:d}f'.format(self.nPart), *getattr(self, self._var_name[i])[n, :, j].astype(np.float32)))
						else:
							sys.exit("Data of type '{:}' not supported".format(self._var_dtype[i]))

					f.write(struct.pack('i', self._traceroutput_tracersize))

				f.close()

			print("Files were written successfully.")


		else:
			sys.exit("Version {:d}-{:02d} not yet implemented".format(self._version // 100, self._version % 100))



def SplitTracerOutput(file_name_base, file_name_in, MaxStepsPerFile = 100, start=0, stop=-1, truncate_from_start=False):
	file_in = open(file_name_in,'rb+')
	header_fname = "{:}_header.dat".format(file_name_base)

	file_out = open(header_fname, 'wb')
	size_i, size_I, size_f, size_d = check_encoding()

	print("Splitting tracer file")
	print("Input: {:}".format(file_name_in))
	print("Header Output: {:}".format(header_fname))

	# Write header block
	traceroutput_headersize = int(struct.unpack('i', file_in.read(size_i))[0])
	version = struct.unpack('i', file_in.read(size_i))[0]
	file_out.write(struct.pack('i', traceroutput_headersize))
	version_new = 202001
	file_out.write(struct.pack('i', version_new))
	block = file_in.read(traceroutput_headersize - size_i)
	file_out.write(block)
	file_out.write(struct.pack('i', traceroutput_headersize))
	del block
	trailing_integer = struct.unpack('i', file_in.read(size_i))[0]
	if trailing_integer != traceroutput_headersize:
		sys.exit("Input Header area not correctly enclosed")

	file_out.close()

	### jump over sum blocks
	if start > 0:
		blocksize = int(struct.unpack('i', file_in.read(size_i))[0])
		file_in.seek(- size_i, 1)
		jumps = start // MaxStepsPerFile
		for i in np.arange(jumps):
			print("Jumping over block {:d}".format(i*MaxStepsPerFile))
			file_in.seek(MaxStepsPerFile*(blocksize + 2*size_i), 1)

	current_file_num = start // MaxStepsPerFile
	current_step = start
	buf = file_in.read(size_i) # Read starting integer of first data block
	if buf:
		data_fname = "{:}_file_{:03d}.dat".format(file_name_base, current_file_num)
		file_out = open(data_fname, 'wb')

		print("Starting data block transfer")
		print("--------------------------------------------------")
		print("File {:}".format(data_fname))


		while(buf):
			blocksize = int(struct.unpack('i', buf)[0])
			print("Block {:d}".format(current_step))
			datablock = file_in.read(blocksize)
			trailing_integer = int(struct.unpack('i', file_in.read(size_i))[0])

			file_out.write(struct.pack('i', blocksize))
			file_out.write(datablock)
			file_out.write(struct.pack('i', blocksize))

			del datablock


			if trailing_integer != blocksize:
				sys.exit("block not correctly enclosed")

			current_step += 1
			if stop >= start and current_step == stop:
				break

			buf = file_in.read(size_i)
			if buf:
				blocksize = int(struct.unpack('i', buf)[0])
				if current_step % MaxStepsPerFile == 0:
					file_out.close()
					current_file_num += 1
					data_fname = "{:}_file_{:03d}.dat".format(file_name_base, current_file_num)
					file_out = open(data_fname, 'wb')
					print("--------------------------------------------------")
					print("File {:}".format(current_file_num, data_fname))

		if truncate_from_start and start > 0:
			file_in.seek(traceroutput_headersize + 2*size_i, 0)
			jumps = start // MaxStepsPerFile
			for i in np.arange(jumps):
				print("Jumping over block {:d}".format(i*MaxStepsPerFile))
				file_in.seek(MaxStepsPerFile*(blocksize + 2*size_i), 1)

			print("Truncate to size: {:}".format(file_in.truncate()))


		file_out.close()

	else:
		data_fname = None
		print("Nothing to transfer")

	file_in.close()








####################################################################################################
def ConvertTracerOutput(file_name, out_version=201901):
	""" Convert pre 2019-01  legacy Arepo tracer output to version 2019-01 and 2019-02"""

	with open(file_name, 'rb') as file_in:
		size_i, size_I, size_f, size_d = check_encoding()

		# Reading first block with unit system
		traceroutput_headersize = int(struct.unpack('i',file_in.read(size_i))[0])
		print("Converting tracer output to version '2019-01'")
		if traceroutput_headersize == 3 * size_d:

			UnitLength_in_cm         = struct.unpack('d', file_in.read(size_d))[0]
			UnitMass_in_g            = struct.unpack('d', file_in.read(size_d))[0]
			UnitVelocity_in_cm_per_s = struct.unpack('d', file_in.read(size_d))[0]


			version = out_version
			if struct.unpack('i', file_in.read(size_i))[0] != traceroutput_headersize:
				sys.exit("Expected header block with size of 3 doubles")

			# now change to new traceroutput_headersize
			traceroutput_headersize = 3 * size_i + 3 * size_d

			traceroutput_tracersize = 2 * size_i + 2 * size_I + 17 * size_f + 1 * size_d
			tracersize_before_temp  =              2 * size_I +  4 * size_f + 1 * size_d
			tracersize_after_temp   = 2 * size_i              + 13 * size_f
			blocksize = int(struct.unpack('i',file_in.read(size_i))[0])

			nPart = blocksize // (traceroutput_tracersize + 1 * size_f) # the old version contains one additional float
			if out_version == 201902:
				traceroutput_tracersize += 2 * size_f # 2019-02 contains 3D magnetic field

		elif traceroutput_headersize == 3 * size_d + 3 * size_i:
			# first 2019-01 version with the variable temp included
			version = struct.unpack('i', f.read(size_i))[0]
			UnitLength_in_cm         = struct.unpack('d', file_in.read(size_d))[0]
			UnitMass_in_g            = struct.unpack('d', file_in.read(size_d))[0]
			UnitVelocity_in_cm_per_s = struct.unpack('d', file_in.read(size_d))[0]
			nPart                    = struct.unpack('i', f.read(size_i))[0]
			traceroutput_tracersize = struct.unpack('i', f.read(size_i))[0]
			tracersize_before_temp  =              2 * size_I +  4 * size_f + 1 * size_d
			tracersize_after_temp   = 2 * size_i              + 13 * size_f
			if traceroutput_tracersize ==  2 * size_i + 2 * size_I + 18 * size_f + 1 * size_d:
				traceroutput_tracersize -= size_f
			else:
				sys.exit("This file already is in version 2019-01 without temperature or 2019-02")

			if out_version == 201902:
				traceroutput_tracersize += 5 * size_f # 2019-02 contains 3D magnetic field and ShockDir
				version = 201902

		else:
			blocksize = 0
			sys.exit("Cannot convert this version")

		# write the new file
		file_out = open(file_name[:-4] + "_new.dat", "xb")

		# header block
		file_out.write(struct.pack('i', 3 * size_i + 3 * size_d))
		file_out.write(struct.pack('i', version))
		file_out.write(struct.pack('d', UnitLength_in_cm))
		file_out.write(struct.pack('d', UnitMass_in_g))
		file_out.write(struct.pack('d', UnitVelocity_in_cm_per_s))
		file_out.write(struct.pack('i', nPart))
		file_out.write(struct.pack('i', traceroutput_tracersize))
		file_out.write(struct.pack('i', traceroutput_headersize))

		# tracer particle blocks
		buf = 1 #if buf > 0 or True next block will be read

		while(buf):

			file_out.write(struct.pack('i', traceroutput_tracersize * nPart)) # starting integer
			file_out.write(file_in.read(tracersize_before_temp * nPart)) # data block before temp
			file_in.seek(nPart * size_f, 1) # jump over the temperature block
			file_out.write(file_in.read(tracersize_after_temp * nPart)) # data block after temp
			file_out.write(struct.pack('i', traceroutput_tracersize * nPart)) # closing integer

			if  int(struct.unpack('i',file_in.read(size_i))[0]) != blocksize:
					sys.exit("data not correctly enclosed")

			buf = file_in.read(size_i) # closing integer of one big block

		file_out.close()


####################################################################################################
# class which handles the parameters given via files to the C program
# parameters are needed for calculating the plots
class ArepoParameters:
	"""Read the parameter file for AREPO simulation"""

	def __init__(self, ParameterFileName = None, verbose = False):
		# Unit system
		self.UnitLength_in_cm = 3.085678e18
		self.UnitMass_in_g = 1.989e33
		self.UnitVelocity_in_cm_per_s = 1.e5
		self.AccelerationEfficiency = 0.1

		self.OutputDir = ''
		self.SnapshotFileBase = ''

		if ParameterFileName is not None:
			self.read_data(ParameterFileName, verbose)

	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def show(self):
		for var in vars(self):
			if type(getattr(self, var)) is int:
				print("{:25} {:d}".format(var,getattr(self,var)))
			if type(getattr(self, var)) is float:
				print("{:25} {:.5e}".format(var,getattr(self,var)))
			if type(getattr(self, var)) is str:
				print("{:25} {:}".format(var,getattr(self,var)))

	# read in the parameter file and set the private class variables accordingly
	def read_data(self,ParameterFileName, verbose = False):
		fParam = open(ParameterFileName,'r')
		if verbose:
			print("Reading parameters from file '{:}'\n".format(ParameterFileName))
		for line in fParam:
			lineParam = (line.strip()).lstrip()

			# ignore lines beginning with a '%' sign
			if(lineParam != ''):
				if(lineParam[0] != '%'):
					columnParam = line.split()
					# only take lines which are of the following format: ParameterTag Space Value
					if(len(columnParam) == 2):
						# loop over all variables to find the one corresponding to the one read from the parameter fiel
						for var in vars(self):
							if(var == columnParam[0]):
								if type(getattr(self, var)) is int:
									setattr(self,var,int(columnParam[1]))
									continue
								elif type(getattr(self, var)) is float:
									setattr(self,var,float(columnParam[1]))
									continue
								elif type(getattr(self, var)) is str:
									setattr(self,var,columnParam[1])
									continue

		if self.OutputDir[-1] != '/':
			self.OutputDir += '/'
		if verbose:
			print("")
			self.show()
			print("")

		line = None
		lineParam = None
		columnParam = None
		fParam.close()
