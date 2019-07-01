import numpy as np
import struct
import sys


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
			blocksize = int(struct.unpack('i', f.read(size_i))[0])
			headersize = 3 * size_i + size_d
			if headersize != blocksize:
				sys.exit("Block size is {:d} bytes, but expexted {:d}".format(blocksize, headersize))

			self.version = int(  struct.unpack('i', f.read(size_i))[0])
			if(verbose):
				print("Version {:d}-{:02d}".format(self.version//100, self.version%100))
			
			self.time    = float(struct.unpack('d', f.read(size_d))[0])
			self.nPart   = int(  struct.unpack('i', f.read(size_i))[0])
			self.nBins   = int(  struct.unpack('i', f.read(size_i))[0])

			if int(struct.unpack('i',f.read(size_i))[0]) != blocksize:
				sys.exit("header data block not correctly enclosed")

			# Momentum Bins
			blocksize = int(struct.unpack('i', f.read(size_i))[0])
			momentumsize = self.nBins * size_d
			if momentumsize != blocksize:
				sys.exit("Block size is {:d} bytes, but expexted {:d}".format(blocksize, momentumsize))
			self.p = np.ndarray(self.nBins, dtype=float)
			self.p[:] = struct.unpack('{:d}d'.format(self.nBins), f.read(size_d * self.nBins))[:]
			
			if  int(struct.unpack('i',f.read(size_i))[0]) != blocksize:
				sys.exit("2nd data block not correctly enclosed")

			# Data
			if not get_only_header:
				# Spectrum Data
				blocksize = int(struct.unpack('i', f.read(size_i))[0])
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
				self.B = np.ndarray(self.nPart, dtype=float)
				self.pos = np.ndarray((self.nPart, 3), dtype=float)

				self.id[:]             = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				self.parent_cell_id[:] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				
				self.mass[:]           = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.n_gas[:]          = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.u_therm[:]        = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.eps_photon[:]     = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.B[:]             = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))

				self.pos[:, 0]        = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.pos[:, 1]        = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.pos[:, 2]        = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				
				for i in np.arange(self.nPart):
					self.f[i, :]      = struct.unpack('{:d}d'.format(self.nBins), f.read(size_d * self.nBins))

				if int(struct.unpack('i',f.read(size_i))[0]) != blocksize:
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

	      $ data = TracerOutput("path_to_file/file_name.dat")

	   and access all x positions

	      $ data.x
	
	   The data structure can be sliced or single particles/snapshots be picked by

	      $ data[1, 2] # particle 1 and snapshot 2
	      $ data[:, 1] # all particles and snapshot 1

	"""

	
	# instance variables
	def __init__(self, file_name = None, version=None, cgs_units = False, verbose = False, read_only_ic= False, specific_particles=None, first_snap=None, last_snap=None, specific_fields=None):
		"""
		Initialize an instance of the TracerOutput.

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
				self._traceroutput_tracersize = 2 * size_i + 2 * size_I + 18 * size_f + 1 * size_d
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
				self._var_name = ['ID', 'time', 'ParentCellID', 'TracerMass', 'x',\
							 'y', 'z', 'n_gas', 'temp', 'u_therm',\
							 'B', 'eps_photon', 'ShockFlag', 'eps_CRpShockInj', 'n_gasPreShock',\
							 'n_gasPostShock', 'VShock', 'timeShockCross', 'theta', 'CReInjection',\
							 'injRate', 'alphaInj', 'pInj']

				# types of the variable
				self._var_dtype = [np.uint32,  np.float64, np.uint32,  np.float32, np.float32,\
							 np.float32, np.float32, np.float32, np.float32, np.float32,\
							 np.float32, np.float32, np.int32,   np.float32, np.float32,\
							 np.float32, np.float32, np.float32, np.float32, np.int32,\
							 np.float32, np.float32, np.float32]

				# cgs scaling of the variable
				self._var_cgs_factor = [1, L/V, 1, M, L,\
										L, L, M / (PROTONMASS * L**3), V**2, V**2,\
										np.sqrt(M * V**2 / L**3), M * V**2 / L**3, 1, M * V**2 / L**3,  M / (PROTONMASS * L**3),\
										M / (PROTONMASS * L**3), V, L/V, 1., 1,\
										V/L, 1., 1.]
			else:
				sys.exit("Version '{:d}-{:d}' not supported or implemented!".format(self._version // 100, self._version % 100))

			
			
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
						setattr(self, self._var_name[i], np.ndarray((self.nPart, self.nSnap), dtype=self._var_dtype[i]))

			else:
				for i in np.arange(len(self._var_name)):
					# Create nPart x nSnap size arrays for all variables if var_store element is true
					if self._var_store[i]:
						setattr(self, self._var_name[i], np.ndarray((self.nSnap), dtype=self._var_dtype[i]))

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
								getattr(self, self._var_name[i])[:, n] = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart)) # equivalent to e.g. self.ID[:, n] = struct.unpack( ... )
							elif self._var_dtype[i] == np.int32:
								getattr(self, self._var_name[i])[:, n] = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
							elif self._var_dtype[i] == np.float32:
								getattr(self, self._var_name[i])[:, n] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
							elif self._var_dtype[i] == np.float64:
								getattr(self, self._var_name[i])[:, n] = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
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
								getattr(self, self._var_name[i])[n] = struct.unpack('I', f.read(size_I))[0]
								f.seek(size_I * (nPartInFile - pos - 1), 1)
								
							elif self._var_dtype[i] == np.int32:
								f.seek(pos * size_i, 1)
								getattr(self, self._var_name[i])[n] = struct.unpack('i', f.read(size_i))[0]
								f.seek(size_i * (nPartInFile - pos - 1), 1)
								
							elif self._var_dtype[i] == np.float32:
								f.seek(pos * size_f, 1)
								getattr(self, self._var_name[i])[n] = struct.unpack('f', f.read(size_f))[0]
								f.seek(size_f * (nPartInFile - pos - 1), 1)
								
							elif self._var_dtype[i] == np.float64:
								f.seek(pos * size_d, 1)
								getattr(self, self._var_name[i])[n] = struct.unpack('d', f.read(size_d))[0]
								f.seek(size_d * (nPartInFile - pos - 1), 1)
								
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
								getattr(self, self._var_name[i])[:, n] = np.array(struct.unpack('{:d}I'.format(nPartInFile), f.read(size_I * nPartInFile)), dtype=np.uint32)[specific_particles]
							elif self._var_dtype[i] == np.int32:
								getattr(self, self._var_name[i])[:, n] = np.array(struct.unpack('{:d}i'.format(nPartInFile), f.read(size_i * nPartInFile)), dtype=np.int32)[specific_particles]
							elif self._var_dtype[i] == np.float32:
								getattr(self, self._var_name[i])[:, n] = np.array(struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile)), dtype=np.float32)[specific_particles]
							elif self._var_dtype[i] == np.float64:
								getattr(self, self._var_name[i])[:, n] = np.array(struct.unpack('{:d}d'.format(nPartInFile), f.read(size_d * nPartInFile)), dtype=np.float64)[specific_particles]
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
		ret = TracerOutput()
		ret.All_Units_in_cgs = self.All_Units_in_cgs
		ret.UnitLength_in_cm = self.UnitLength_in_cm
		ret.UnitMass_in_g = self.UnitMass_in_g
		ret.UnitVelocity_in_cm_per_s = self.UnitVelocity_in_cm_per_s

		ret._var_name = self._var_name
		ret._var_dtype = self._var_dtype
		ret._var_cgs_factor = self._var_cgs_factor
		ret._var_store = self._var_store

		if isinstance(key, int):
			ret.nPart = 1
			ret.nSnap = self.nSnap
		elif isinstance(key, slice):
			start, stop, step = key.indices(self.nPart)
			ret.nPart = np.arange(start, stop, step).size
			ret.nSnap = self.nSnap
		elif isinstance(key, tuple):
			if len(key) == 2:
				if isinstance(key[0], int):
					ret.nPart = 1
				elif isinstance(key[0], slice):
					start, stop, step = key[0].indices(self.nPart)
					ret.nPart = np.arange(start, stop, step).size
				else:
					raise TypeError('Index must be int or slice, not {}'.format(type(key[0]).__name__))

				if isinstance(key[1], int):
					ret.nSnap = 1
				elif isinstance(key[1], slice):
					start, stop, step = key[1].indices(self.nSnap)
					ret.nSnap = np.arange(start, stop, step).size
				else:
					raise TypeError('Index must be int or slice, not {}'.format(type(key[1]).__name__))
		else:
			raise TypeError('Tuple Index must be of length 2, not {}'.format(len(key)))

		
		for i in np.arange(len(self._var_name)):
			if self._var_store[i]:
				setattr(ret, ret._var_name[i], getattr(self, self._var_name[i]).__getitem__(key))

		return ret


	
####################################################################################################
def ConvertTracerOutput(file_name):
	""" Convert pre 2019-01  legacy Arepo tracer output to version 2019-01"""
	
	with open(file_name, 'rb') as file_in:
		size_i, size_I, size_f, size_d = check_encoding()

		# Reading first block with unit system
		traceroutput_headersize = int(struct.unpack('i',file_in.read(size_i))[0])
		if traceroutput_headersize == 3 * size_d:
			print("Converting tracer output to equivalent version '2019-01'")
			UnitLength_in_cm         = struct.unpack('d', file_in.read(size_d))[0]
			UnitMass_in_g            = struct.unpack('d', file_in.read(size_d))[0]
			UnitVelocity_in_cm_per_s = struct.unpack('d', file_in.read(size_d))[0]

	
			version = 201901
			if struct.unpack('i', file_in.read(size_i))[0] != traceroutput_headersize:
				sys.exit("Expected header block with size of 3 doubles")

			# now change to new traceroutput_headersize
			traceroutput_headersize = 3 * size_i + 3 * size_d
				
			traceroutput_tracersize = 2 * size_i + 2 * size_I + 18 * size_f + 1 * size_d
			blocksize = int(struct.unpack('i',file_in.read(size_i))[0]) 
			
			nPart = blocksize // traceroutput_tracersize
			traceroutput_tracerblocksize = traceroutput_tracersize * nPart

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
			file_out.write(file_in.read(traceroutput_tracerblocksize)) # data block
			file_out.write(struct.pack('i', traceroutput_tracersize * nPart)) # closing integer

			if  int(struct.unpack('i',file_in.read(size_i))[0]) != blocksize:
					sys.exit("data not correctly enclosed")

			buf = file_in.read(size_i) # closing integer of one big block

		file_out.close()
