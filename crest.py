import numpy as np
import struct
import sys
from Physics import PROTONMASS


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
		self.SnapshotFileBase            = ''

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
		self.CosmologicalIntegrationOn   = 0

		# Flags
		self.FlagAllowSubcycles          = 1
		self.FlagCooling                 = 1
		self.Flag_Fermi_I_Reacceleration = 1
		self.Flag_Fermi_I_injection      = 1
		self.Flag_Fermi_II_Reacceleration= 1
		self.FlagExternalInjection       = 0

		# Cooling & Diffusion
		self.n_elec                      = 1.157
		self.HydrogenMassFrac            = 0.76
		self.DiffusionTimeInGyr          = 0.
		self.Lambda                      = 0.
		self.Radiation_Field_in_eps_CMB  = -1. # if this value is -1 then it was not set
		self.Magnetic_Field_Amplification = -1. # if this value is -1 then it was not set

		# parameters for shock injection
		self.ShockParamA                 = 0.
		self.ShockParamB                 = 0.
		self.ShockParamC                 = 0.
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
		self.MaximumMomentumNumerics = 0.
		self.MinimumMomentumNumerics = 0.

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
	
	def __init__(self, file_name = None, verbose = False, get_only_header = False):
		if file_name is not None:
			self.read_data(file_name, verbose, get_only_header)

	def read_data(self, file_name, verbose = False, get_only_header = False):
		size_i, size_I, size_f, size_d = check_encoding()
		with open(file_name,'rb') as f:
			if(verbose):
				print("Reading snapshot data from file '{:}'".format(file_name))

			
			# Header information
			blocksize = int(struct.unpack('I', f.read(size_i))[0])
			headersize = 3 * size_i + size_d
			if headersize != blocksize:
				sys.exit("Block size is {:d} bytes, but expexted {:d}".format(blocksize, headersize))

			self.version = int(  struct.unpack('i', f.read(size_i))[0])
			if(verbose):
				print("Version {:d}-{:02d}".format(self.version//100, self.version%100))
			
			self.time    = float(struct.unpack('d', f.read(size_d))[0])
			self.nPart   = int(  struct.unpack('i', f.read(size_i))[0])
			self.nBins   = int(  struct.unpack('i', f.read(size_i))[0])

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
				# Spectrum Data
				blocksize = int(struct.unpack('I', f.read(size_i))[0])
				if self.version == 201902:
					datasize = self.nPart * ( self.nBins * size_d + 2 * size_I + 10 * size_d)
				else:
					datasize = self.nPart * ( self.nBins * size_d + 2 * size_I + 8 * size_d)
				if blocksize != datasize:
					sys.exit("Block size is {:d} bytes, but expexted {:d}".format(blocksize, datasize))

				self.f = np.ndarray((self.nPart, self.nBins), dtype = float)
				self.id = np.ndarray(self.nPart, dtype=np.uint32)
				self.parent_cell_id = np.ndarray(self.nPart, dtype=np.uint32)
				self.mass = np.ndarray(self.nPart, dtype=float)
				self.n_gas = np.ndarray(self.nPart, dtype=float)
				self.u_therm = np.ndarray(self.nPart, dtype=float)
				self.eps_photon = np.ndarray(self.nPart, dtype=float)
				self.pos = np.ndarray((self.nPart, 3), dtype=float)				

				if self.version==201902:
					self.B = np.ndarray((self.nPart, 3), dtype=float)
				else:
					self.B = np.ndarray(self.nPart, dtype=float)
				

				self.id[:]             = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				self.parent_cell_id[:] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				
				self.mass[:]           = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.n_gas[:]          = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.u_therm[:]        = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.eps_photon[:]     = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				
				if self.version==201902:
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
	The header will be read in to provide basic information about the used unit system, as
	all quantities will be in code units of the AREPO simulation.
	All variables will be stored as (nPart, nSnap) arrays unless only one specific particle
	is read in then the arrays have the shape (nSnap, ).

	Example:
	   Load a tracer output with path and name to file by, e.g.,

	      $ data = ArepoTracerOutput("path_to_file/file_name.dat")

	   and access all x positions

	      $ data.x
	
	   The data structure can be sliced or single particles/snapshots be picked by

	      $ data[1, 2] # particle 1 and snapshot 2
	      $ data[:, 1] # all particles and snapshot 1

	"""

	
	# instance variables
	def __init__(self, file_name = None, version=None, cgs_units = False, verbose = False, read_only_ic= False, specific_particles=None, first_snap=None, last_snap=None, specific_fields=None):
		"""
		Initialize an instance of ArepoTracerOutput.

		If no parameters are given an empty instance is initialized.
		If a path to a file is provided the file will be read in 
		
		Args:
		   file_name (str): Path of file to be read
		   
		   cgs_units (bool): Flag if the values should be converted to cgs units immediately

		   read_only_ic (bool): Read only header and 0th snapshot/initial conditions

		   first_snap (int): First snapshot to be read in

		   last_snap (int): Last snapshot to be read in (exclusive)

		   specific_fields (list): List of strings of the variable which should be stored, 
		                           e.g., ['ID', 'time']
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

		if file_name is not None:
			self.read_data(file_name, version=version, cgs_units=cgs_units, verbose=verbose, read_only_ic=read_only_ic, specific_particles=specific_particles, first_snap=first_snap, last_snap=last_snap, specific_fields=specific_fields)

	@property
	def var_name(self):
		""" List of names of the variables which are in the file. Variables are stored if corresponding value of var_store is True. """
		return self._var_name

	@property
	def version(self):
		""" List of names of the variables which are in the file. Variables are stored if corresponding value of var_store is True. """
		return self._version

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

	def read_data(self, file_name, version=None, cgs_units = False, verbose = False, read_only_ic = False, specific_particles = None, first_snap = None, last_snap = None, specific_fields=None, UnitLength_in_cm = 1., UnitMass_in_g = 1., UnitVelocity_in_cm_per_s = 1.):
		""" Read in the data from file. This function is automatically called if class constructor is called with the file name"""
		
		with open(file_name,'rb') as f:
			if verbose:
				print("Read only initial conditions: {:}".format(read_only_ic))
				print("Read Arepo's tracer output from file '{}'".format(file_name))
			size_i, size_I, size_f, size_d = check_encoding()

			def blocksize(type_arr): 
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

			# Version dependend configurati
			# we need to make sure that old simulations can be read in
			# python tool for tracer output conversion??
			

			
			# definitions 

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

				self._version                 = struct.unpack('i', f.read(size_i))[0]
				self.UnitLength_in_cm         = struct.unpack('d', f.read(size_d))[0]
				self.UnitMass_in_g            = struct.unpack('d', f.read(size_d))[0]
				self.UnitVelocity_in_cm_per_s = struct.unpack('d', f.read(size_d))[0]
				self.nPart                    = struct.unpack('i', f.read(size_i))[0]
				self._traceroutput_tracersize = struct.unpack('i', f.read(size_i))[0]

			else:
							
				sys.exit("Expected header block with size of 3 doubles (old style) or 3 doubles plus 3 integers.")

			if  int(struct.unpack('i', f.read(size_i))[0]) != self._traceroutput_headersize:
				sys.exit("Expected header block with size of 3 doubles (old style) or 3 doubles plus 3 integers.")

			L = self.UnitLength_in_cm
			M = self.UnitMass_in_g
			V = self.UnitVelocity_in_cm_per_s


			# names of all possible variables
			# for better readability: new line after 5 elements
			# the order of the variables has to be same as in the file


			if self._version == 201901:
				self._var_name = ['ID', 'time', 'ParentCellID', 'TracerMass', 'pos',\
								  'n_gas', 'u_therm', 'B', 'eps_photon', 'ShockFlag',\
								  'eps_CRpShockInj', 'n_gasPreShock', 'n_gasPostShock', 'VShock', 'timeShockCross',\
								  'theta', 'CReInjection', 'injRate', 'alphaInj', 'pInj']\
							 

				# types of the variable
				# we assume the types 'np.array' to be of lenght 3 and dtype np.float32
				self._var_dtype = [np.uint32,  np.float64, np.uint32,  np.float32, np.ndarray,\
								   np.float32, np.float32, np.float32, np.float32, np.int32,\
								   np.float32, np.float32, np.float32, np.float32, np.float32,\
								   np.float32, np.int32,   np.float32, np.float32, np.float32]
								   

				if blocksize(self._var_dtype) != self._traceroutput_tracersize:
					sys.exit("Size of tracer data block given in file is not the same as in code!")

				# cgs scaling of the variable
				B = np.sqrt(M * V**2 / L**3) # B scaling
				N = M / (PROTONMASS * L**3) # Gas Density Scaling
				E = M * V**2 / L**3 # Energy Density Scaling
				self._var_cgs_factor = [1, L/V, 1, M, L,\
										N, V**2, B, E, 1,\
										E, N, N, V, L/V,\
										1, 1, V/L, 1, 1]
										

			elif self._version == 201902:
				self._var_name = ['ID', 'time', 'ParentCellID', 'TracerMass', 'pos',\
								  'n_gas', 'u_therm', 'B', 'eps_photon', 'ShockFlag',\
								  'eps_CRpShockInj', 'n_gasPreShock', 'n_gasPostShock', 'VShock', 'timeShockCross',\
								  'ShockDir',  'theta', 'CReInjection', 'injRate', 'alphaInj',\
								  'pInj']

				# types of the variable
				self._var_dtype = [np.uint32,  np.float64, np.uint32,  np.float32, np.ndarray,\
								   np.float32, np.float32, np.ndarray, np.float32, np.int32,\
								   np.float32, np.float32, np.float32, np.float32, np.float32,\
								   np.ndarray, np.float32, np.int32,   np.float32, np.float32,\
								   np.float32]

				if blocksize(self._var_dtype) != self._traceroutput_tracersize:
					sys.exit("Size of tracer data block given in file is not the same as in code!")

				# cgs scaling of the variable
				B = np.sqrt(M * V**2 / L**3) # B scaling
				N = M / (PROTONMASS * L**3) # Gas Density Scaling
				E = M * V**2 / L**3 # Energy Density Scaling
				self._var_cgs_factor = [1, L/V, 1, M, L,\
										N, V**2, B, E, 1,\
										E, N, N, V, L/V,\
										1, 1, 1, V/L, 1,\
										1]

			else:
				sys.exit("Version '{:d}-{:02d}' not supported or implemented!".format(self._version // 100, self._version % 100))
			
			
			# Reading block with data values
			blocksize = int(struct.unpack('i',f.read(size_i))[0])
			nPartInFile = self.nPart
			self.nSnap	= 0
			buf 	= 1 #if buf > 0 or True next block will be read		   					
 
			while(buf):
				# move pointer forward
				f.seek(self.nPart * self._traceroutput_tracersize, 1) 
				if  int(struct.unpack('i',f.read(size_i))[0]) != blocksize:
					sys.exit("data not in block #{:i} not correctly enclosed".format(self.nSnap))

				self.nSnap += 1
				if read_only_ic or last_snap == self.nSnap:
					buf = False
				else:
					buf = f.read(size_i)

			if specific_particles is not None:
				if type(specific_particles) is int:
					self.nPart = 1
				else:
					self.nPart = len(np.arange(self.nPart)[specific_particles])

			# go back to the beginning of the file
			f.seek(3*size_i + self._traceroutput_headersize, 0)

			if first_snap is not None:
				self.nSnap -= first_snap
		
			buf = 0
			if verbose:
				print("Number of particles: {:d}".format(self.nPart))
				print("Number of snapshots: {:d}".format(self.nSnap))

			if len(self._var_name) != len(self._var_dtype) or len(self._var_name) != len(self._var_cgs_factor):
				sys.exit("Arrays of variable names and types need to have the same length")


			# will the variable be stored or not
			self._var_store = np.ones(len(self._var_name), dtype=bool) # array filled with True


			if specific_fields is not None:
				# make no storage default, invert var_store array
				np.logical_not(self._var_store, out=self._var_store)

				for sf in specific_fields:
					if sf in self._var_name:
						self._var_store[ self._var_name.index(sf) ] = True
					else:
						sys.exit("Variable '{:}' does not exist in tracer output!".format(sf))
					
			
			if type(specific_particles) is not int:
				for i in np.arange(len(self._var_name)):
					# Create nPart x nSnap size arrays for all variables if var_store element is true
					if self._var_store[i]:
						if self._var_dtype[i] is not np.ndarray: # scalar variable
							setattr(self, self._var_name[i], np.ndarray((self.nSnap, self.nPart), dtype=self._var_dtype[i]))
						else: # 3D Vector
							setattr(self, self._var_name[i], np.ndarray((self.nSnap, self.nPart, 3), dtype=np.float32))

			else:
				for i in np.arange(len(self._var_name)):
					# Create nPart x nSnap size arrays for all variables if var_store element is true
					if self._var_store[i]:
						if self._var_dtype[i] is not np.ndarray: # scalar
							setattr(self, self._var_name[i], np.ndarray((self.nSnap, 1), dtype=self._var_dtype[i]))
						else: # 3D Vector
							setattr(self, self._var_name[i], np.ndarray((self.nSnap, 1, 3), dtype=np.float32))



			if first_snap is not None:
				# skip some lines
				f.seek(first_snap * (blocksize + 2*size_i), 1) 

			if specific_particles is None:
				# read all the data spec
				for n in np.arange(self.nSnap):

					# loop over all possible variables
					for i in np.arange(len(self._var_name)):
						
						if self._var_store[i]:
							# if variable should be stored, check the type and read from file in the correct format
							if self._var_dtype[i] == np.uint32:
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
							if self._var_dtype[i] == np.uint32:
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

					f.seek(size_i, 1)

			elif type(specific_particles) is int:
				pos = specific_particles
				# read single data
				for n in np.arange(self.nSnap):
					# loop over all possible variables
					for i in np.arange(len(self._var_name)):
						
						if self._var_store[i]:
							# if variable should be stored, check the type and read from file in the correct format
							# jump to the position of the desired particle in the variable block
							# read in the right amount of bytes
							# jump to end of the variable block
							if self._var_dtype[i] == np.uint32:
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
							if self._var_dtype[i] == np.uint32:
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
				for n in np.arange(self.nSnap):

					# loop over all possible variables
					for i in np.arange(len(self._var_name)):
						
						if self._var_store[i]:
							# if variable should be stored, check the entire block of the variable
							# convert to numpy array which supports indexing with arrays, slices, etc.
							# (return type of struct.unpack is 'tuple')
							if self._var_dtype[i] == np.uint32:
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
							if self._var_dtype[i] == np.uint32:
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

			f.close()
			if verbose:
				print("Data was successfully read")


			if cgs_units:
				self.scale_to_cgs_units(verbose)

	def scale_to_cgs_units(self, verbose=False):
		if not self.All_Units_in_cgs:
			if verbose:
				print("Scale to cgs with UnitLenght_in_cm = {:.3e}, UnitMass_in_g = {:.3e}, UnitVeloctiy_in_cm_per_s = {:.3e}".format(self.UnitLength_in_cm, self.UnitMass_in_g, self.UnitVelocity_in_cm_per_s))
			for i in np.arange(len(self._var_name)):
				if self._var_store[i]:
					setattr(self, self._var_name[i], np.multiply(self._var_cgs_factor[i], getattr(self,self._var_name[i])).astype(self._var_dtype[i]))

			self.All_Units_in_cgs = True

		else:
			print("Variables are already stored in cgs units")

	def scale_to_code_units(self, verbose=False):
		if self.All_Units_in_cgs:
			if verbose:
				print("Scale to code units with UnitLenght_in_cm = {:.3e}, UnitMass_in_g = {:.3e}, UnitVeloctiy_in_cm_per_s = {:.3e}".format(self.UnitLength_in_cm, self.UnitMass_in_g, self.UnitVelocity_in_cm_per_s))
			for i in np.arange(len(self._var_name)):
				if self._var_store[i]:
					setattr(self, self._var_name[i], np.divide(self._var_cgs_factor[i], getattr(self,self._var_name[i])).astype(self._var_dtype[i]))

			self.All_Units_in_cgs = False

		else:
			print("Variables are already stored in code units")


				

	def __getitem__(self, key):
		# check dimensions of return
		ret = ArepoTracerOutput()
		ret.All_Units_in_cgs = self.All_Units_in_cgs
		ret.UnitLength_in_cm = self.UnitLength_in_cm
		ret.UnitMass_in_g = self.UnitMass_in_g
		ret.UnitVelocity_in_cm_per_s = self.UnitVelocity_in_cm_per_s

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
				setattr(ret, ret._var_name[i], getattr(self, self._var_name[i]).__getitem__(key))

		return ret

	
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

