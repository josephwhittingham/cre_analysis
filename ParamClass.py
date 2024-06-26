import matplotlib.pyplot as plt
import numpy as np
import struct
import sys
from Physics import *


####################################################################################################
# class which handles the parameters given via files to the C program
# parameters are needed for calculating the plots
class CRelectronParameters:
	"""Read the parameter file for a CREST simulation"""

	def __init__(self,ParameterFileName = None, RunID = None, verbose = False):
		"""
		Initialize all possible variables which can be defined by the parameter file.
		Note that it contains some legacy variables.
		If a path to the parameter file is given the file will be read in
		"""

		self.Run_ID                      = 0

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
		self.x_e                         = 1.157
		self.HydrogenMassFrac            = 0.76
		self.DiffusionTimeInGyr          = 0.
		self.Lambda                      = 0.
		self.Radiation_Field_in_eps_CMB  = -1. # if this value is -1 then it was not set

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

		if ParameterFileName is not None:
			self.read_data(ParameterFileName, RunID, verbose)

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
	def read_data(self,ParameterFileName, RunID = None, verbose = False):
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
		if RunID != None:
			self.Run_ID = RunID
		line = None
		lineParam = None
		columnParam = None
		fParam.close()

	def BinsPerDec(self):
		return (self.NumberOfMomentumBins - self.IncludeMaximumMomentum) / int(np.log10(self.MaximumMomentum) - np.log10(self.MinimumMomentum))




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



####################################################################################################
# class which handles all the Snapshot which was also provided to the C program
# by setting the variable 'InputDataFile'.
# Parameters can be read out by giving the name of this file
# it is useful to take the value provided by the CRelectronParameters class
class SnapshotData:
	NumberOfSnapshotsMax = 50
	
	# instance variables
	def __init__(self,DataFileName = None):
		self.NumberOfSnapshots = 0
		self.time    = np.zeros(self.NumberOfSnapshotsMax)
		#self.pos_x    = np.zeros(self.NumberOfSnapshotsMax)
		#self.pos_y    = np.zeros(self.NumberOfSnapshotsMax)
		#self.pos_z    = np.zeros(self.NumberOfSnapshotsMax)
		self.n_gas    = np.zeros(self.NumberOfSnapshotsMax)
		#self.u        = np.zeros(self.NumberOfSnapshotsMax)
		#self.vel_x    = np.zeros(self.NumberOfSnapshotsMax)
		#self.vel_y    = np.zeros(self.NumberOfSnapshotsMax)
		#self.vel_z    = np.zeros(self.NumberOfSnapshotsMax)
		#self.dedt    = np.zeros(self.NumberOfSnapshotsMax)
		self.B        = np.zeros(self.NumberOfSnapshotsMax)
		if DataFileName is not None:
			self.read_data(DataFileName)
	
	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def read_data(self,DataFileName):
		fData = open(DataFileName)
		print("Reading data from file '{:}'\n".format(DataFileName))
		print("\t time (s) \t n_gas (cm^-3) \t B (G) \n")        
		#print("\t time \t\t pos_x \t\t pos_y \t\t pos_z \t\t n_gas \t\t u \t\t vel_x \t\t vel_y \t\t vel_z \t\t dedt \n")
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.time[self.NumberOfSnapshots]     = float(columnData[0])
					#self.pos_x[self.NumberOfSnapshots]    = float(columnData[1])
					#self.pos_y[self.NumberOfSnapshots]     = float(columnData[2])
					#self.pos_z[self.NumberOfSnapshots]     = float(columnData[3])
					self.n_gas[self.NumberOfSnapshots]     = float(columnData[4])
					#self.u[self.NumberOfSnapshots]         = float(columnData[5])
					#self.vel_x[self.NumberOfSnapshots]     = float(columnData[6])
					#self.vel_y[self.NumberOfSnapshots]     = float(columnData[7])
					#self.vel_z[self.NumberOfSnapshots]     = float(columnData[8])
					#self.dedt[self.NumberOfSnapshots]     = float(columnData[9])
					self.B[self.NumberOfSnapshots]     = float(columnData[10])
					#print("\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[self.NumberOfSnapshots], self.pos_x[self.NumberOfSnapshots], self.pos_y[self.NumberOfSnapshots], self.pos_z[self.NumberOfSnapshots], self.n_gas[self.NumberOfSnapshots], self.u[self.NumberOfSnapshots], self.vel_x[self.NumberOfSnapshots], self.vel_y[self.NumberOfSnapshots], self.vel_z[self.NumberOfSnapshots], self.dedt[self.NumberOfSnapshots]))
					print("\t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[self.NumberOfSnapshots], self.n_gas[self.NumberOfSnapshots],self.B[self.NumberOfSnapshots]))
					self.NumberOfSnapshots += 1
		self.time.resize(self.NumberOfSnapshots)
		self.n_gas.resize(self.NumberOfSnapshots)
		#self.pos_x.resize(self.NumberOfSnapshots)
		#self.pos_y.resize(self.NumberOfSnapshots)
		#self.pos_z.resize(self.NumberOfSnapshots)
		#self.u.resize(self.NumberOfSnapshots)
		#self.vel_x.resize(self.NumberOfSnapshots)
		#self.vel_y.resize(self.NumberOfSnapshots)
		#self.vel_z.resize(self.NumberOfSnapshots)
		#self.dedt.resize(self.NumberOfSnapshots)


		print("\n")
		line = None
		lineData = None
		columnData = None
		fData.close()

	def show(self):
		print("Number of Snapshots: {:d}".format(self.NumberOfSnapshots))
		print("\t time (s) \t n_gas (cm^-3) \t B (G) \n")
		for i in np.arange(NumberOfSnapshots):
			print("\t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[i], self.n_gas[i],self.B[i]))
		print("")

####################################################################################################
# class which handles all the Snapshot which was also provided to the C program
# by setting the variable 'InputDataFile'.
# Parameters can be read out by giving the name of this file
# it is useful to take the value provided by the CRelectronParameters class
# 
# 2nd Version: Parameter file for CRelectronTracer Programm
class SnapshotData2:
	NumberOfSnapshotsMax = 50
	
	# instance variables
	def __init__(self,DataFileName = None):
		self.NumberOfSnapshots = 0
		self.time    = np.ndarray(self.NumberOfSnapshotsMax,dtype=float)
		self.n_gas    = np.ndarray(self.NumberOfSnapshotsMax,dtype=float)
		self.T        = np.ndarray(self.NumberOfSnapshotsMax,dtype=float)
		self.B        = np.ndarray(self.NumberOfSnapshotsMax,dtype=float)
		self.inj    = np.ndarray(self.NumberOfSnapshotsMax,dtype=int)
		if DataFileName is not None:
			self.read_data(DataFileName)
	
	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def read_data(self,DataFileName):
		fData = open(DataFileName)
		print("Reading data from file '{:}'\n".format(DataFileName))
		print("\t time (s) \t n_gas (cm^-3) \t T (K) \t B (G) \t inj \n")
		#print("\t time \t\t pos_x \t\t pos_y \t\t pos_z \t\t n_gas \t\t u \t\t vel_x \t\t vel_y \t\t vel_z \t\t dedt \n")
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.time[self.NumberOfSnapshots]     = float(columnData[0])
					self.n_gas[self.NumberOfSnapshots]     = float(columnData[1])
					self.T[self.NumberOfSnapshots]         = float(columnData[2])
					self.B[self.NumberOfSnapshots]         = float(columnData[3])
					self.inj[self.NumberOfSnapshots]     = int(columnData[4])
					print("\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:d}".format(self.time[self.NumberOfSnapshots], self.n_gas[self.NumberOfSnapshots],self.T[self.NumberOfSnapshots], self.B[self.NumberOfSnapshots], self.inj[self.NumberOfSnapshots]))
					self.NumberOfSnapshots += 1
		self.time.resize(self.NumberOfSnapshots)
		self.n_gas.resize(self.NumberOfSnapshots)
		self.T.resize(self.NumberOfSnapshots)
		self.B.resize(self.NumberOfSnapshots)
		self.inj.resize(self.NumberOfSnapshots)

		print("\n")
		line = None
		lineData = None
		columnData = None
		fData.close()

	def show(self):
		print("Number of Snapshots: {:d}".format(self.NumberOfSnapshots))
		print("\t time (s) \t n_gas (cm^-3) \t T (K) \t B (G) \t inj \n")       
		for i in np.arange(NumberOfSnapshots):
			print("\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:d}".format(self.time[self.NumberOfSnapshots], self.n_gas[self.NumberOfSnapshots],self.T[self.NumberOfSnapshots], self.B[self.NumberOfSnapshots], self.inj[self.NumberOfSnapshots]))
		print("")

	
if __name__ == "__main__":
	run_ID = 99

	filename = '../Input/param_{:02d}.txt'.format(run_ID)
	
	#testParams = CRelectronParameters()
	#testParams.read_parameters(filename)
	testParams = CRelectronParameters(filename)

	#testData = SnapshotData()
	#testData.read_data('../{:}'.format(testParams.InputDataFile))
	testData = SnapshotData('../{:}'.format(testParams.InputDataFile))



####################################################################################################
# Simple Class to Read Error File
class runErrorRelative:
	nPointsMax = 30
	def __init__(self,runErrorFile = None):
		self.nPoints = 0
		self.time    = np.zeros(self.nPointsMax)
		self.n_gas    = np.zeros(self.nPointsMax)
		self.err    = np.zeros(self.nPointsMax)
		
		if runErrorFile is not None:
			self.read_data(runErrorFile)
	
	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def read_data(self, runErrorFile):
		print("Reading data from file '{:}'".format(runErrorFile))
		fData = open(runErrorFile,'r')
		self.nPoints = 0
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.time[self.nPoints] = columnData[0]
					self.n_gas[self.nPoints]    = columnData[1]
					self.err[self.nPoints]    = columnData[2]
					print("\t {:E} \t {:E} \t {:E}".format(self.time[self.nPoints], self.n_gas[self.nPoints], self.err[self.nPoints]))
					self.nPoints+=1
		self.time.resize(self.nPoints)
		self.n_gas.resize(self.nPoints)
		self.err.resize(self.nPoints)

		print("\n")
		line = None
		lineData = None
		columnData = None
		fData.close()

####################################################################################################
# Simple Class to Read Error for CR properties
class PropertyError:
	nPointsMax = 30
	def __init__(self,runErrorFile = None):
		self.nPoints = 0
		self.time    = np.zeros(self.nPointsMax)
		self.del_energy    = np.zeros(self.nPointsMax)
		self.del_pressure    = np.zeros(self.nPointsMax)
		self.del_number    = np.zeros(self.nPointsMax)
		
		if runErrorFile is not None:
			self.read_data(runErrorFile)
	
	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def read_data(self, runErrorFile):
		print("Reading data from file '{:}'".format(runErrorFile))
		fData = open(runErrorFile,'r')
		self.nPoints = 0
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.time[self.nPoints] = columnData[0]
					self.del_energy[self.nPoints]    = columnData[1]
					self.del_pressure[self.nPoints]    = columnData[2]
					self.del_number[self.nPoints]    = columnData[2]
					print("\t {:E} \t {:E} \t {:E}".format(self.time[self.nPoints], self.del_energy[self.nPoints], self.del_pressure[self.nPoints], self.del_number[self.nPoints]))
					self.nPoints+=1
		self.time.resize(self.nPoints)
		self.del_energy.resize(self.nPoints)
		self.del_pressure.resize(self.nPoints)
		self.del_number.resize(self.nPoints)

		print("\n")
		line = None
		lineData = None
		columnData = None
		fData.close()

####################################################################################################
# Simple Class to Read Approximate Solution for Protons
class ApproxData:
	nSnapsMax = 50
	def __init__(self,approxFile = None):
		self.time = np.zeros(self.nSnapsMax)
		self.Norm = np.zeros(self.nSnapsMax)
		self.Cut  = np.zeros(self.nSnapsMax)
		self.nSnaps = 0

		if approxFile is not None:
			self.read_data(approxFile)

	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def read_data(self,approxFile):
		print("Reading approximate solution from file '{:}'".format(approxFile))
		fData = open(approxFile,'r')
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.time[self.nSnaps] = float(columnData[0])
					self.Norm[self.nSnaps] = float(columnData[1])
					self.Cut[self.nSnaps]  = float(columnData[2])
					print("\t {:E} \t {:E} \t {:E}".format(self.time[self.nSnaps], self.Norm[self.nSnaps], self.Cut[self.nSnaps]))
					self.nSnaps += 1
		self.time.resize(self.nSnaps)
		self.Norm.resize(self.nSnaps)
		self.Cut.resize(self.nSnaps)

		print("\n")
		line = None
		lineData = None
		columnData = None
		fData.close()

####################################################################################################
# Simple Class for the distribution function
# takes the snapshot file to fill the data
class DistributionFunction:
	def __init__(self,SnapshotDataFile,NumberOfMomentumBins):
		self.p = np.zeros(NumberOfMomentumBins)
		self.f = np.zeros(NumberOfMomentumBins)
		print("Reading snapshot data from file '{:}'".format(SnapshotDataFile))
		fSnap = open(SnapshotDataFile,'r')
		i = 0
		for line in fSnap:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.p[i] = columnData[0]
					self.f[i] = columnData[1]
					i+=1

		self.p.resize(i)
		self.f.resize(i)
		i = None
		fSnap.close()


####################################################################################################
# Read the distribution function
def SpectrumSnapshot(fname,  nBinsIn=None, NoIdTimeHeader=False):

	size_i, size_I, size_f, size_d = checkNumberEncoding()
	with open(fname,'rb') as file:
		print("Reading snapshot data from file '{:}'".format(fname))

		# Read the first block with id and time
		if not(NoIdTimeHeader):
			dummy = int(struct.unpack('i', file.read(size_i))[0])
			if dummy != size_i + size_d:
				sys.exit("First block size is {:d} bytes, but expexted {:d}".format(dummy, size_i + size_d))

			id = int(struct.unpack('i', file.read(size_i))[0])
			time = float(struct.unpack('d', file.read(size_d))[0])

			file.seek(size_i, 1)
		else:
			id = -1
			time = -1.
				


		# Read first information block with number of particles and momentum bins
		dummy = int(struct.unpack('i', file.read(size_i))[0])
		if nBinsIn != None:
			if dummy != nBinsIn * size_d:
				sys.exit("Block size is {:d} bytes, but expected {:d}".format(dummy, nBins * size_d))

		nBins = dummy // size_d
	
		f = np.array(struct.unpack('{:d}d'.format(nBins), file.read(size_d * nBins)))

		dummy = int(struct.unpack('i',file.read(size_i))[0])
		if dummy != nBins * size_d:
			sys.exit("1st data block not correctly enclosed")

		file.close()

	return f, id, time


####################################################################################################
#
class cre_snapshot:
	def __init__(self, fname = None, verbose = False, get_only_header = False):
		if fname is not None:
			self.read_data(fname, verbose, get_only_header)

	def read_data(self, fname, verbose = False, get_only_header = False):
		size_i, size_I, size_f, size_d = checkNumberEncoding()
		with open(fname,'rb') as file:
			if(verbose):
				print("Reading snapshot data from file '{:}'".format(fname))

			
			# Header information
			dummy = int(struct.unpack('i', file.read(size_i))[0])
			blocksize = 2 * size_i + size_d
			if dummy != blocksize:
				sys.exit("Block size is {:d} bytes, but expexted {:d}".format(dummy, blocksize))

			self.time = float(struct.unpack('d', file.read(size_d))[0])
			self.nPart = int(struct.unpack('i', file.read(size_i))[0])
			self.nBins = int(struct.unpack('i', file.read(size_i))[0])

			dummy = int(struct.unpack('i',file.read(size_i))[0])
			if dummy != blocksize:
				sys.exit("1st data block not correctly enclosed")

			# Momentum Bins
			dummy = int(struct.unpack('i', file.read(size_i))[0])
			blocksize = self.nBins * size_d
			if dummy != blocksize:
				sys.exit("Block size is {:d} bytes, but expexted {:d}".format(dummy, blocksize))			

			self.p = np.ndarray(self.nBins, dtype=float)
			self.p[:] = struct.unpack('{:d}d'.format(self.nBins), file.read(size_d * self.nBins))[:]
			
			dummy = int(struct.unpack('i',file.read(size_i))[0])
			if dummy != blocksize:
				sys.exit("2nd data block not correctly enclosed")

			if not get_only_header:
				# Spectrum Data
				dummy = int(struct.unpack('i', file.read(size_i))[0])
				blocksize = self.nPart * ( self.nBins * size_d + size_I + 6 * size_d)
				if dummy != blocksize:
					sys.exit("Block size is {:d} bytes, but expexted {:d}".format(dummy, blocksize))

				self.f = np.ndarray((self.nPart, self.nBins), dtype = float)
				self.id = np.ndarray(self.nPart, dtype=np.uint32)
				self.mass = np.ndarray(self.nPart, dtype=float)
				self.n_gas = np.ndarray(self.nPart, dtype=float)
				self.u_therm = np.ndarray(self.nPart, dtype=float)
				self.pos = np.ndarray((self.nPart, 3), dtype=float)
				for i in np.arange(self.nPart):
					self.id[i] = struct.unpack('I', file.read(size_I))[0]
					self.mass[i] = struct.unpack('d', file.read(size_d))[0]
					self.n_gas[i] = struct.unpack('d', file.read(size_d))[0]
					self.u_therm[i] = struct.unpack('d', file.read(size_d))[0]
					self.pos[i, :] = struct.unpack('3d', file.read(3 * size_d))[:]
					self.f[i, :] = struct.unpack('{:d}d'.format(self.nBins), file.read(size_d * self.nBins))[:]

				dummy = int(struct.unpack('i',file.read(size_i))[0])
				if dummy != blocksize:
					sys.exit("3rd data block not correctly enclosed")
			

			

		
			



####################################################################################################
## gives two numbers for representing a float in scientific notation
## 23. will return [2.3,1]
def exp_rep(f):
	if f==0:
		f_exp = int(0)
		f_coeff = 0.0
	else:
		f_exp = int(np.floor(np.log10(abs(f)) ))
		f_coeff = f / np.power(10.,f_exp)
	return [f_coeff, f_exp]


####################################################################################################
# Simple Class for the distribution function
# takes the snapshot file to fill the data
class OtherSolution:
	def __init__(self,DataFile,nMax=400):
		self.p = np.zeros(nMax)
		self.f = np.zeros(nMax)
		print("Reading other solution data from file '{:}'".format(DataFile))
		fData = open(DataFile,'r')
		i = 0
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.p[i] = columnData[0]
					self.f[i] = columnData[1]
					i+=1

		self.p.resize(i)
		self.f.resize(i)
		i = None
		fData.close()


####################################################################################################
# New Tracer Output with more options

class TracerOutput:
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
	def __init__(self, fname = None, version=None, cgs_units = False, verbose = False, read_only_ic= False, specific_particles=None, firstSnap=None, lastSnap=None, specific_fields=None):
		"""
		Initialize an instance of the TracerOutput.

		If no parameters are given an empty instance is initialized.
		If a path to a file is provided the file will be read in 
		
		Args:
		   fname (str): Path of file to be read
		   
		   cgs_units (bool): Flag if the values should be converted to cgs units immediately

		   read_only_ic (bool): Read only header and 0th snapshot/initial conditions

		   firstSnap (int): First snapshot to be read in

		   lastSnap (int): Last snapshot to be read in (exclusive)

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

		if fname is not None:
			self.read_data(fname, version=version, cgs_units=cgs_units, verbose=verbose, read_only_ic=read_only_ic, specific_particles=specific_particles, firstSnap=firstSnap, lastSnap=lastSnap, specific_fields=specific_fields)

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

	def read_data(self, fname, version=None, cgs_units = False, verbose = False, read_only_ic = False, specific_particles = None, firstSnap = None, lastSnap = None, specific_fields=None, UnitLength_in_cm = 1., UnitMass_in_g = 1., UnitVelocity_in_cm_per_s = 1.):
		with open(fname,'rb') as f:
			if verbose:
				print("Read only initial conditions: {:}".format(read_only_ic))
				print("Read Arepo's tracer output from file '{}'".format(fname))
			size_i, size_I, size_f, size_d = checkNumberEncoding()


			# Version dependend configurati
			# we need to make sure that old simulations can be read in
			# python tool for tracer output conversion??
			

			
			# definitions 

			# Reading first block with unit system
			self._traceroutput_headersize = int(struct.unpack('i',f.read(size_i))[0])
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
			dummy = int(struct.unpack('i',f.read(size_i))[0])
			nPartInFile = self.nPart
			self.nSnap	= 0
			buf 	= 1 #if buf > 0 or True next block will be read		   					
 
			while(buf):
				# move pointer forward
				f.seek(self.nPart * self._traceroutput_tracersize, 1) 
				if  int(struct.unpack('i',f.read(size_i))[0]) != dummy:
					sys.exit("data not in block #{:i} not correctly enclosed".format(self.nSnap))

				self.nSnap += 1
				if read_only_ic or lastSnap == self.nSnap:
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

			if firstSnap is not None:
				self.nSnap -= firstSnap
		
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

			if firstSnap is not None:
				# skip some lines
				f.seek(firstSnap * (dummy + 2*size_i), 1) 

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

					if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
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

					if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
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

					if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
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

#### Function to convert legacy tracer output
def ConvertTracerOutput(fname):
	with open(fname, 'rb') as file_in:
		size_i, size_I, size_f, size_d = checkNumberEncoding()

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
			dummy = int(struct.unpack('i',file_in.read(size_i))[0]) 
			
			nPart = dummy // traceroutput_tracersize
			traceroutput_tracerblocksize = traceroutput_tracersize * nPart

		else:
			dummy = 0
			sys.exit("Cannot convert this version")

		# write the new file
		file_out = open(fname[:-4] + "_new.dat", "xb")

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

			if  int(struct.unpack('i',file_in.read(size_i))[0]) != dummy:
					sys.exit("data not correctly enclosed")

			buf = file_in.read(size_i) # closing integer of one big block

		file_out.close()
		
		
			

# ####################################################################################################
# # Class to read the tracer output data
# # takes the binary output file of arepo to extract data
# # gets directly number of snapshots
# # needs packages: struct
		
# class TracerOutput:
# 	# instance variables
# 	def __init__(self, fname = None, cgs_units = False, verbose = False, read_only_ic= False, specific_particles=None, firstSnap=None, lastSnap=None, specific_fields=None):
# 		""" e.g. specific fields = ['time', 'ShockFlag']""" 
		
# 		# with_cr_electrons is set to 1 if arepo was compiled with #COSMIC_RAYS_ELECTRONS
# 		# need to set dummy values as these determine the types
# 		self.nSnap = 0
# 		self.nPart = 0
# 		self.All_Units_in_cgs = False
# 		self.UnitLength_in_cm = 1.
# 		self.UnitMass_in_g = 1.
# 		self.UnitVelocity_in_cm_per_s = 1.

# 		if fname is not None:
# 			self.read_data(fname, cgs_units=cgs_units, verbose=verbose, read_only_ic=read_only_ic, specific_particles=specific_particles, firstSnap=firstSnap, lastSnap=lastSnap, specific_fields=specific_fields)


# 	def __del__(self):
# 		for var in vars(self):
# 			setattr(self,var,None)

# 	def read_data(self, fname, cgs_units = False, verbose = False, read_only_ic = False, specific_particles = None, firstSnap = None, lastSnap = None, specific_fields=None, UnitLength_in_cm = 1., UnitMass_in_g = 1., UnitVelocity_in_cm_per_s = 1.):
# 		with open(fname,'rb') as f:
# 			if verbose:
# 				print("Read only initial conditions: {:}".format(read_only_ic))
# 				print("Read Arepo's tracer output from file '{}'".format(fname))
# 			size_i, size_I, size_f, size_d = checkNumberEncoding()

# 			# Reading first block with unit system
# 			dummy = int(struct.unpack('i',f.read(size_i))[0])
# 			if dummy != 3 * size_d:
# 				sys.exit("Expected 3 double values at beginning, ")

# 			self.UnitLength_in_cm         = struct.unpack('d', f.read(size_d))[0]
# 			self.UnitMass_in_g            = struct.unpack('d', f.read(size_d))[0]
# 			self.UnitVelocity_in_cm_per_s = struct.unpack('d', f.read(size_d))[0]

# 			if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
# 				sys.exit("Expected 3 double values at beginning, ")
			
			
# 			# Reading block with data values
# 			dummy = int(struct.unpack('i',f.read(size_i))[0])
# 			self.nPart = dummy // (2 * size_i + 2 * size_I + 18 * size_f + 1 * size_d)
# 			nPartInFile = self.nPart
# 			self.nSnap	= 0
# 			buf 	= 1		   					

# 			while(buf):
# 				# move pointer forward
# 				f.seek(self.nPart * (2 * size_i + 2 * size_I + 18 * size_f + 1 * size_d), 1) 
# 				if  int(struct.unpack('i',f.read(size_i))[0]) != dummy:
# 					sys.exit("data not correctly enclosed 1, ")

# 				self.nSnap += 1
# 				if read_only_ic or lastSnap == self.nSnap:
# 					buf = False
# 				else:
# 					buf = f.read(size_i)

# 			if specific_particles is not None:
# 				# set attributes for reading a single particle
# 				if type(specific_particles) is int:
# 					if specific_particles <= self.nPart:
# 						self.nPart = 1
# 					else:
# 						sys.exit("Cannot read particle {:d} as there are only {:d} particles".format(specific_particles, self.nPart))
# 				elif type(specific_particles) is slice:
# 					# set attributes for reading a slice
# 					if specific_particles.stop <= self.nPart:
# 						if type(specific_particles) is int:
# 							self.nPart = (specific_particles.stop - specific_particles.start) // specific_particles.step
# 							if (specific_particles.stop - specific_particles.start) % specific_particles.step > 0:
# 								self.nPart += 1
# 						else:
# 							self.nPart = (specific_particles.stop - specific_particles.start)



# 					else:
# 						sys.exit("Cannot read particles until {:d} as there are only {:d} particles".format(specific_particles.stop, self.nPart))
# 				else:
# 					sys.exit("specific_particles = {:} not supported".format(specific_particles))

# 			# go back to the beginning of the file
# 			f.seek(3*size_i + 3*size_d, 0)

# 			if firstSnap is not None:
# 				self.nSnap -= firstSnap
		
# 			buf = 0
# 			if verbose:
# 				print("Number of particles: {:d}".format(self.nPart))
# 				print("Number of snapshots: {:d}".format(self.nSnap))

# 			if specific_fields is None:
# 				self.store_ID             = True
# 				self.store_time           = True
# 				self.store_ParentCellID   = True
# 				self.store_TracerMass     = True

# 				self.store_x	            = True
# 				self.store_y	            = True
# 				self.store_z	            = True
# 				self.store_n_gas	        = True
# 				self.store_temp	        = True
# 				self.store_u_therm        = True

# 				# parameters for cooling
# 				self.store_B			    = True
# 				self.store_eps_photon       = True

# 				# parameters for diffusive shock acceleration
# 				self.store_ShockFlag      = True
# 				self.store_eps_CRpShockInj  = True
# 				self.store_n_gasPreShock  = True
# 				self.store_n_gasPostShock = True
# 				self.store_VShock         = True
# 				self.store_timeShockCross = True
# 				self.store_theta       = True

# 				# parameters for electron injection (apart from DSA given above)
# 				self.store_CReInjection   = True
# 				self.store_injRate        = True
# 				self.store_alphaInj       = True
# 				self.store_pInj           = True
				
# 			else:
# 				self.store_ID             = False
# 				self.store_time           = False
# 				self.store_ParentCellID   = False
# 				self.store_TracerMass     = False

# 				self.store_x	            = False
# 				self.store_y	            = False
# 				self.store_z	            = False
# 				self.store_n_gas	        = False
# 				self.store_temp	        = False
# 				self.store_u_therm        = False

# 				# parameters for cooling
# 				self.store_B			    = False
# 				self.store_eps_photon       = False

# 				# parameters for diffusive shock acceleration
# 				self.store_ShockFlag      = False
# 				self.store_eps_CRpShockInj  = False
# 				self.store_n_gasPreShock  = False
# 				self.store_n_gasPostShock = False
# 				self.store_VShock         = False
# 				self.store_timeShockCross = False
# 				self.store_theta       = False

# 				# parameters for electron injection (apart from DSA given above)
# 				self.store_CReInjection   = False
# 				self.store_injRate        = False
# 				self.store_alphaInj       = False
# 				self.store_pInj           = False

# 				for sf in specific_fields:
# 					if 'store_'+sf in vars(self):
# 						setattr(self, 'store_' + sf, True)
# 					else:
# 						sys.exit("Variable '{:}' does not exist in tracer output!".format(sf))
					
			

# 			if type(specific_particles) is not int:
# 				# create the arrays
# 				if self.store_ID :
# 					self.ID             = np.ndarray((self.nPart, self.nSnap), dtype=np.uint32)
# 				if self.store_time :
# 					self.time           = np.ndarray((self.nPart, self.nSnap), dtype=float) # time or scale parameter
# 				if self.store_ParentCellID :
# 					self.ParentCellID   = np.ndarray((self.nPart, self.nSnap), dtype=np.uint32)
# 				if self.store_TracerMass :
# 					self.TracerMass     = np.ndarray((self.nPart, self.nSnap), dtype=float)

# 				if self.store_x :
# 					self.x	            = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_y :
# 					self.y	            = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_z :
# 					self.z	            = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_n_gas :
# 					self.n_gas	        = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_temp :
# 					self.temp	        = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_u_therm :
# 					self.u_therm        = np.ndarray((self.nPart, self.nSnap), dtype=float)

# 				# parameters for cooling
# 				if self.store_B :
# 					self.B			    = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_eps_photon :
# 					self.eps_photon       = np.ndarray((self.nPart, self.nSnap), dtype=float)

# 				# parameters for diffusive shock acceleration
# 				if self.store_ShockFlag :
# 					self.ShockFlag      = np.ndarray((self.nPart, self.nSnap), dtype=int)
# 				if self.store_eps_CRpShockInj :
# 					self.eps_CRpShockInj  = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_n_gasPreShock :
# 					self.n_gasPreShock  = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_n_gasPostShock :
# 					self.n_gasPostShock = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_VShock :
# 					self.VShock         = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_timeShockCross :
# 					self.timeShockCross = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_theta :
# 					self.theta       = np.ndarray((self.nPart, self.nSnap), dtype=float)

# 				# parameters for electron injection (apart from DSA given above)
# 				if self.store_CReInjection :
# 					self.CReInjection   = np.ndarray((self.nPart, self.nSnap), dtype=int)
# 				if self.store_injRate :
# 					self.injRate        = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_alphaInj :
# 					self.alphaInj       = np.ndarray((self.nPart, self.nSnap), dtype=float)
# 				if self.store_pInj :
# 					self.pInj           = np.ndarray((self.nPart, self.nSnap), dtype=float)

# 			else:
# 				# create the arrays
# 				if self.store_ID :
# 					self.ID             = np.ndarray(self.nSnap, dtype=np.uint32)
# 				if self.store_time :
# 					self.time           = np.ndarray(self.nSnap, dtype=float) # time or scale parameter
# 				if self.store_ParentCellID :
# 					self.ParentCellID   = np.ndarray(self.nSnap, dtype=np.uint32)
# 				if self.store_TracerMass :
# 					self.TracerMass     = np.ndarray(self.nSnap, dtype=float)

# 				if self.store_x :
# 					self.x	            = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_y :
# 					self.y	            = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_z :
# 					self.z	            = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_n_gas :
# 					self.n_gas	        = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_temp :
# 					self.temp	        = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_u_therm :
# 					self.u_therm        = np.ndarray(self.nSnap, dtype=float)

# 					# parameters for cooling
# 				if self.store_B :
# 					self.B			    = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_eps_photon :
# 					self.eps_photon       = np.ndarray(self.nSnap, dtype=float)

# 					# parameters for diffusive shock acceleration
# 				if self.store_ShockFlag :
# 					self.ShockFlag      = np.ndarray(self.nSnap, dtype=int)
# 				if self.store_eps_CRpShockInj :
# 					self.eps_CRpShockInj  = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_n_gasPreShock :
# 					self.n_gasPreShock  = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_n_gasPostShock :
# 					self.n_gasPostShock = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_VShock :
# 					self.VShock         = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_timeShockCross :
# 					self.timeShockCross = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_theta :
# 					self.theta       = np.ndarray(self.nSnap, dtype=float)

# 					# parameters for electron injection (apart from DSA given above)
# 				if self.store_CReInjection :
# 					self.CReInjection   = np.ndarray(self.nSnap, dtype=int)
# 				if self.store_injRate :
# 					self.injRate        = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_alphaInj :
# 					self.alphaInj       = np.ndarray(self.nSnap, dtype=float)
# 				if self.store_pInj :
# 					self.pInj           = np.ndarray(self.nSnap, dtype=float)

# 			if firstSnap is not None:
# 				# skip some lines
# 				f.seek(firstSnap * (dummy + 2*size_i), 1) 




# 			if specific_particles is None:
# 				# read all the data
# 				for n in np.arange(self.nSnap):

# 					# if store_VARNAME is set the block for VARNAME will be read in and stored into a numpy array of the name VARNAME
# 					# if not the block will be skipped
# 					if self.store_ID :
# 						self.ID[:, n]             = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
# 					else:
# 						f.seek(size_I * self.nPart, 1)
						
# 					if self.store_time :
# 						self.time[:, n]		      = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
# 					else:
# 						f.seek(size_d * self.nPart, 1)
						
# 					if self.store_ParentCellID :
# 						self.ParentCellID[:, n]   = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
# 					else:
# 						f.seek(size_I * self.nPart, 1)
						
# 					if self.store_TracerMass :
# 						self.TracerMass[:, n]     = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_x :
# 						self.x[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_y :
# 						self.y[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_z :
# 						self.z[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_n_gas :
# 						self.n_gas[:, n]		  = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_temp :
# 						self.temp[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_u_therm :
# 						self.u_therm[:, n]	      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))

# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_B :
# 						self.B[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_eps_photon :
# 						self.eps_photon[:, n]       = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_ShockFlag :
# 						self.ShockFlag[:, n]      = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
# 					else:
# 						f.seek(size_i * self.nPart, 1)
						
# 					if self.store_eps_CRpShockInj :
# 						self.eps_CRpShockInj[:, n]  = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_n_gasPreShock :
# 						self.n_gasPreShock[:, n]  = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_n_gasPostShock :
# 						self.n_gasPostShock[:, n] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_VShock :
# 						self.VShock[:, n]         = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_timeShockCross :
# 						self.timeShockCross[:, n] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_theta :
# 						self.theta[:, n]          = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))

# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_CReInjection :
# 						self.CReInjection[:, n]   = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
# 					else:
# 						f.seek(size_i * self.nPart, 1)
						
# 					if self.store_injRate :
# 						self.injRate[:, n]        = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_alphaInj :
# 						self.alphaInj[:, n]       = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)
						
# 					if self.store_pInj :
# 						self.pInj[:, n]           = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
# 					else:
# 						f.seek(size_f * self.nPart, 1)


# 					if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
# 						sys.exit("data not correctly enclosed, ")

# 					f.seek(size_i, 1)

# 			elif type(specific_particles) is int:
# 				pos = specific_particles
# 				# read single data
# 				for n in np.arange(self.nSnap):

# 					if self.store_ID :
# 						f.seek(pos * size_I, 1)
# 						self.ID[n]             = struct.unpack('I', f.read(size_I))[0]
# 						f.seek((nPartInFile - pos - 1) * size_I, 1)
# 					else:
# 						f.seek(size_I * nPartInFile, 1)
						
# 					if self.store_time :
# 						f.seek(pos * size_d, 1)
# 						self.time[n]		      = struct.unpack('d', f.read(size_d))[0]
# 						f.seek((nPartInFile - pos - 1) * size_d, 1)
# 					else:
# 						f.seek(size_d * nPartInFile, 1)
						
# 					if self.store_ParentCellID :
# 						f.seek(pos * size_I, 1)
# 						self.ParentCellID[n]   = struct.unpack('I', f.read(size_I))[0]
# 						f.seek((nPartInFile - pos - 1) * size_I, 1)
# 					else:
# 						f.seek(size_I * nPartInFile, 1)
						
# 					if self.store_TracerMass :
# 						f.seek(pos * size_f, 1)
# 						self.TracerMass[n]     = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_x :
# 						f.seek(pos * size_f, 1)
# 						self.x[n]		      = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_y :
# 						f.seek(pos * size_f, 1)
# 						self.y[n]		      = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_z :
# 						f.seek(pos * size_f, 1)
# 						self.z[n]		      = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_n_gas :
# 						f.seek(pos * size_f, 1)
# 						self.n_gas[n]		  = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_temp :
# 						f.seek(pos * size_f, 1)
# 						self.temp[n]		      = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_u_therm :
# 						f.seek(pos * size_f, 1)
# 						self.u_therm[n]	      = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_B :
# 						f.seek(pos * size_f, 1)
# 						self.B[n]		      = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_eps_photon :
# 						f.seek(pos * size_f, 1)
# 						self.eps_photon[n]       = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_ShockFlag :
# 						f.seek(pos * size_i, 1)
# 						self.ShockFlag[n]      = struct.unpack('i', f.read(size_i))[0]
# 						f.seek((nPartInFile - pos - 1) * size_i, 1)
# 					else:
# 						f.seek(size_i * nPartInFile, 1)
						
# 					if self.store_eps_CRpShockInj :
# 						f.seek(pos * size_f, 1)
# 						self.eps_CRpShockInj[n]  = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_n_gasPreShock :
# 						f.seek(pos * size_f, 1)
# 						self.n_gasPreShock[n]  = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_n_gasPostShock :
# 						f.seek(pos * size_f, 1)
# 						self.n_gasPostShock[n] = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_VShock :
# 						f.seek(pos * size_f, 1)
# 						self.VShock[n]         = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_timeShockCross :
# 						f.seek(pos * size_f, 1)
# 						self.timeShockCross[n] = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_theta :
# 						f.seek(pos * size_f, 1)
# 						self.theta[n]       = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_CReInjection :
# 						f.seek(pos * size_i, 1)
# 						self.CReInjection[n]   = struct.unpack('i', f.read(size_i))[0]
# 						f.seek((nPartInFile - pos - 1) * size_i, 1)
# 					else:
# 						f.seek(size_i * nPartInFile, 1)
						
# 					if self.store_injRate :
# 						f.seek(pos * size_f, 1)
# 						self.injRate[n]        = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_alphaInj :
# 						f.seek(pos * size_f, 1)
# 						self.alphaInj[n]       = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_pInj :
# 						f.seek(pos * size_f, 1)
# 						self.pInj[n]           = struct.unpack('f', f.read(size_f))[0]
# 						f.seek((nPartInFile - pos - 1) * size_f, 1)
# 					else:
# 						f.seek(size_f * nPartInFile, 1)


# 					if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
# 						sys.exit("data not correctly enclosed, ")

# 					f.seek(size_i, 1)

# 			else:
# 				# read all the data but just take a slice out of it
# 				for n in np.arange(self.nSnap):
# 					if self.store_ID :
# 						self.ID[:, n]             = struct.unpack('{:d}I'.format(nPartInFile), f.read(size_I * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_I * nPartInFile, 1)
						
# 					if self.store_time :
# 						self.time[:, n]		      = struct.unpack('{:d}d'.format(nPartInFile), f.read(size_d * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_d * nPartInFile, 1)
						
# 					if self.store_ParentCellID :
# 						self.ParentCellID[:, n]   = struct.unpack('{:d}I'.format(nPartInFile), f.read(size_I * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_I * nPartInFile, 1)
						
# 					if self.store_TracerMass :
# 						self.TracerMass[:, n]     = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]

# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_x :
# 						self.x[:, n]		      = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_y :
# 						self.y[:, n]		      = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_z :
# 						self.z[:, n]		      = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_n_gas :
# 						self.n_gas[:, n]		  = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_temp :
# 						self.temp[:, n]		      = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_u_therm :
# 						self.u_therm[:, n]	      = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]

# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_B :
# 						self.B[:, n]		      = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_eps_photon :
# 						self.eps_photon[:, n]       = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_ShockFlag :
# 						self.ShockFlag[:, n]      = struct.unpack('{:d}i'.format(nPartInFile), f.read(size_i * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_i * nPartInFile, 1)
						
# 					if self.store_eps_CRpShockInj :
# 						self.eps_CRpShockInj[:, n]  = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_n_gasPreShock :
# 						self.n_gasPreShock[:, n]  = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_n_gasPostShock :
# 						self.n_gasPostShock[:, n] = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]

# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_VShock :
# 						self.VShock[:, n]         = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_timeShockCross :
# 						self.timeShockCross[:, n] = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_theta :
# 						self.theta[:, n]       = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]

# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_CReInjection :
# 						self.CReInjection[:, n]   = struct.unpack('{:d}i'.format(nPartInFile), f.read(size_i * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_i * nPartInFile, 1)
						
# 					if self.store_injRate :
# 						self.injRate[:, n]        = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_alphaInj :
# 						self.alphaInj[:, n]       = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)
						
# 					if self.store_pInj :
# 						self.pInj[:, n]           = struct.unpack('{:d}f'.format(nPartInFile), f.read(size_f * nPartInFile))[specific_particles]
# 					else:
# 						f.seek(size_f * nPartInFile, 1)


# 					if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
# 						sys.exit("data not correctly enclosed, ")

# 					f.seek(size_i, 1)

				

# 			f.close()
# 			if verbose:
# 				print("Data was successfully read")
			
# 			if cgs_units:
# 				self.scale_to_cgs(verbose)

# 	def scale_to_cgs(self, verbose=False):	
# 		if verbose:
# 			print("Scale to cgs with UnitLenght_in_cm = {:.3e}, UnitMass_in_g = {:.3e}, UnitVeloctiy_in_cm_per_s = {:.3e}".format(self.UnitLength_in_cm, self.UnitMass_in_g, self.UnitVelocity_in_cm_per_s))
# 		self.time           = np.multiply(self.time, self.UnitLength_in_cm / self.UnitVelocity_in_cm_per_s)
# 		self.TracerMass     = np.multiply(self.TracerMass, self.UnitMass_in_g)
# 		self.x              = np.multiply(self.x, self.UnitLength_in_cm)
# 		self.y              = np.multiply(self.y, self.UnitLength_in_cm)
# 		self.z              = np.multiply(self.z, self.UnitLength_in_cm)
# 		self.n_gas          = np.multiply(self.n_gas, self.UnitMass_in_g / (PROTONMASS * np.power(self.UnitLength_in_cm, 3)))
# 		self.temp           = np.multiply(self.temp, np.square(self.UnitVelocity_in_cm_per_s))
# 		self.u_therm        = np.multiply(self.u_therm, np.square(self.UnitVelocity_in_cm_per_s))

# 		self.B              = np.multiply(self.B, np.sqrt(self.UnitMass_in_g / self.UnitLength_in_cm) * self.UnitVelocity_in_cm_per_s / self.UnitLength_in_cm)
# 		self.eps_photon       = np.multiply(self.eps_photon, self.UnitMass_in_g * np.square(self.UnitVelocity_in_cm_per_s) / np.power(self.UnitLength_in_cm, 3))

# 		self.eps_CRpShockInj  = np.multiply(self.eps_CRpShockInj, self.UnitMass_in_g * np.square(self.UnitVelocity_in_cm_per_s) / np.power(self.UnitLength_in_cm, 3))

# 		self.n_gasPreShock  = np.multiply(self.n_gasPreShock,  self.UnitMass_in_g / (PROTONMASS * np.power(self.UnitLength_in_cm, 3)))	
# 		self.n_gasPostShock = np.multiply(self.n_gasPostShock, self.UnitMass_in_g / (PROTONMASS * np.power(self.UnitLength_in_cm, 3)))

# 		self.VShock         = np.multiply(self.VShock, self.UnitVelocity_in_cm_per_s)		
# 		self.timeShockCross = np.multiply(self.timeShockCross, self.UnitLength_in_cm / self.UnitVelocity_in_cm_per_s)

# 		self.injRate        = np.multiply(self.injRate, self.UnitVelocity_in_cm_per_s / self.UnitLength_in_cm)

# 	def __getitem__(self, key):
# 		# check dimensions of return
# 		ret = TracerOutput()
# 		ret.All_Units_in_cgs = self.All_Units_in_cgs
# 		ret.UnitLength_in_cm = self.UnitLength_in_cm
# 		ret.UnitMass_in_g = self.UnitMass_in_g
# 		ret.UnitVelocity_in_cm_per_s = self.UnitVelocity_in_cm_per_s

# 		if isinstance(key, int):
# 			ret.nPart = 1
# 			ret.nSnap = self.nSnap
# 		elif isinstance(key, slice):
# 			start, stop, step = key.indices(self.nPart)
# 			ret.nPart = np.arange(start, stop, step).size
# 			ret.nSnap = self.nSnap
# 		elif isinstance(key, tuple):
# 			if len(key) == 2:
# 				if isinstance(key[0], int):
# 					ret.nPart = 1
# 				elif isinstance(key[0], slice):
# 					start, stop, step = key[0].indices(self.nPart)
# 					ret.nPart = np.arange(start, stop, step).size
# 				else:
# 					raise TypeError('Index must be int or slice, not {}'.format(type(key[0]).__name__))

# 				if isinstance(key[1], int):
# 					ret.nSnap = 1
# 				elif isinstance(key[1], slice):
# 					start, stop, step = key[1].indices(self.nSnap)
# 					ret.nSnap = (stop - start + 1)//step
# 				else:
# 					raise TypeError('Index must be int or slice, not {}'.format(type(key[1]).__name__))
# 		else:
# 			raise TypeError('Tuple Index must be of length 2, not {}'.format(len(key)))	

# 		# create the arrays
# 		ret.ID             = self.ID.__getitem__(key)
# 		ret.time           = self.time.__getitem__(key)
# 		ret.ParentCellID   = self.ParentCellID.__getitem__(key)
# 		ret.TracerMass     = self.TracerMass.__getitem__(key)

# 		ret.x	            = self.x.__getitem__(key)
# 		ret.y	            = self.y.__getitem__(key)
# 		ret.z	            = self.z.__getitem__(key)
# 		ret.n_gas	        = self.n_gas.__getitem__(key)
# 		ret.temp	        = self.temp.__getitem__(key)
# 		ret.u_therm        = self.u_therm.__getitem__(key)

# 		# parameters for cooling
# 		ret.B			    = self.B.__getitem__(key)
# 		ret.eps_photon       = self.eps_photon.__getitem__(key)

# 		# parameters for diffusive shock acceleration
# 		ret.ShockFlag      = self.ShockFlag.__getitem__(key)
# 		ret.eps_CRpShockInj  = self.eps_CRpShockInj.__getitem__(key)
# 		ret.n_gasPreShock  = self.n_gasPreShock.__getitem__(key)
# 		ret.n_gasPostShock = self.n_gasPostShock.__getitem__(key)
# 		ret.VShock         = self.VShock.__getitem__(key)
# 		ret.timeShockCross = self.timeShockCross.__getitem__(key)
# 		ret.theta       = self.theta.__getitem__(key)

# 		# parameters for electron injection (apart from DSA given above)
# 		ret.CReInjection   = self.CReInjection.__getitem__(key)
# 		ret.injRate        = self.injRate.__getitem__(key)
# 		ret.alphaInj       = self.alphaInj.__getitem__(key)
# 		ret.pInj           = self.pInj.__getitem__(key)		

# 		return ret




####################################################################################################
# Class to read the tracer output data
# takes the binary output file of arepo to extract data
# gets directly number of snapshots
# needs packages: struct
		
class TracerOutputOld:
	# instance variables
	def __init__(self, fname = None, with_cr_electrons=1):
		# with_cr_electrons is set to 1 if arepo was compiled with #COSMIC_RAYS_ELECTRONS
		# need to set dummy values as these determine the types
		checkNumberEncoding()
		self.nSnap = 0
		self.nPart = 0
		if fname is not None:
			self.read_data(fname,with_cr_electrons)

	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def read_data(self, fname,with_cr_electrons=1):
		with open(fname,'rb') as f:
			print("Read Arepo's tracer output from file '{}'".format(fname))
			size_i, size_I, size_f, size_d = checkNumberEncoding()
		
			dummy = int(struct.unpack('i',f.read(size_i))[0])
			if with_cr_electrons == 1:
				self.nPart = (dummy - size_d) // (12 * size_f + size_i)
			else:
				self.nPart = (dummy - size_d) // (6 * size_f)

			self.nSnap	= 0
			buf 	= 1

			while(buf):
				# move pointer forward
				f.seek(size_d + 6 * self.nPart * size_f, 1) 
				if with_cr_electrons:
					# move pointer forward
					f.seek(6 * self.nPart * size_f + self.nPart * size_i, 1)

				if  int(struct.unpack('i',f.read(size_i))[0]) != dummy:
					sys.exit("data not correctly enclosed 1, ")

				self.nSnap += 1
				buf = f.read(size_i)

			# go back to the beginning of the file
			f.seek(size_i, 0)
			buf = 0
			print("Number of particles: {:d}".format(self.nPart))
			print("Number of snapshots: {:d}".format(self.nSnap))

			# create the arrays    
			self.time	= np.ndarray(self.nSnap, dtype = float)
			self.redshift = np.ndarray(self.nSnap, dtype = float) # TODO: Read real data or calculate them
			self.x	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.y	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.z	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.n_gas	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.temp	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.u_therm = np.ndarray((self.nPart, self.nSnap), dtype=float)

			if with_cr_electrons:
				self.eps_photon   = np.ndarray((self.nPart, self.nSnap), dtype=float) # TODO: Read real data or calculate them
				self.B			= np.ndarray((self.nPart, self.nSnap), dtype=float)
				self.Dist		= np.ndarray((self.nPart, self.nSnap), dtype=float)
				self.Comp		= np.ndarray((self.nPart, self.nSnap), dtype=float)
				self.BRat		= np.ndarray((self.nPart, self.nSnap), dtype=float)
				self.eInjection		= np.ndarray((self.nPart, self.nSnap), dtype=int)
				self.eEnergyPerMass	= np.ndarray((self.nPart, self.nSnap), dtype=float)
				self.V			= np.ndarray((self.nPart, self.nSnap), dtype=float)

			# read the data
			for n in np.arange(self.nSnap):
				self.time[n]		= float(struct.unpack('d', f.read(size_d))[0])
				self.x[:, n]		= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.y[:, n]		= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.z[:, n]		= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.n_gas[:, n]		= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.temp[:, n]		= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.u_therm[:, n]	= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				if with_cr_electrons:
					self.B[:, n]			= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
					self.eInjection[:, n]		= struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
					self.eEnergyPerMass[:, n]		= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))	
					self.Dist[:, n]			= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
					self.Comp[:, n]			= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
					self.BRat[:, n]			= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
					self.V[:, n]			= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))

				if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
					sys.exit("data not correctly enclosed, ")

				f.seek(size_i, 1)

			f.close()

			print("Data was successfully read")


	


####################################################################################################
# Function to check whether system encoding of int, float and double is correct
# returns the sizes of int, float and double if everything is alright
def checkNumberEncoding():
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
# Function for writing out the tracer data in arepostyle
# returns 0 if everything is alright
def writeTracerArepo(fileName, nSnap, nPart, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s, ID, time, ParentCellID, TracerMass, x, y, z, n_gas, temp, u_therm, B, eps_photon, ShockFlag, eps_CRpShockInj, n_gasPreShock, n_gasPostShock, VShock, timeShockCross, theta, CReInjection, injRate, alphaInj, pInj):
	
	size_i, size_I, size_f, size_d = checkNumberEncoding()

	# do some consistency checks
	if ID.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('ID', ID.shape, nPart, nSnap))

	if time.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('time', time.shape, nPart, nSnap))

	if ParentCellID.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('ParentCellID', ParentCellID.shape, nPart, nSnap))

	if TracerMass.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('TracerMass', TracerMass.shape, nPart, nSnap))
		
	if x.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('x', pos_x.shape, nPart, nSnap))

	if y.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('y', y.shape, nPart, nSnap))

	if z.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('z', z.shape, nPart, nSnap))

	if n_gas.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('n_gas', n_gas.shape, nPart, nSnap))

	if temp.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('temp', temp.shape, nPart, nSnap))

	if u_therm.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('u_therm', u_therm.shape, nPart, nSnap))

	if B.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('B', B.shape, nPart, nSnap))

	if eps_photon.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('eps_photon', eps_photon.shape, nPart, nSnap))

	if ShockFlag.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('ShockFlag', ShockFlag.shape, nPart, nSnap))

	if eps_CRpShockInj.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('eps_CRpShockInj', eps_CRpShockInj.shape, nPart, nSnap))
		
	if n_gasPreShock.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('n_gasPreShock', n_gasPreShock.shape, nPart, nSnap))

	if n_gasPostShock.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('n_gasPostShock', n_gasPostShock.shape, nPart, nSnap))

	if VShock.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('VShock', VShock.shape, nPart, nSnap))

	if timeShockCross.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('timeShockCross', timeShockCross.shape, nPart, nSnap))

	if theta.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('theta', theta.shape, nPart, nSnap))

	if CReInjection.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('CReInjection', CReInjection.shape, nPart, nSnap))

	if injRate.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('injRate', injRate.shape, nPart, nSnap))

	if alphaInj.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('alphaInj', alphaInj.shape, nPart, nSnap))

	if pInj.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('pInj', pInj.shape, nPart, nSnap))

	float_buffer = np.ndarray(nPart,dtype=float)
	int_buffer   = np.ndarray(nPart,dtype=int)
	uint_buffer  = np.ndarray(nPart,dtype=np.uint32)

	with open(fileName,'wb') as f:
		dummy = 3 * size_d
		f.write(struct.pack('i',dummy))
		f.write(struct.pack('d',UnitLength_in_cm))
		f.write(struct.pack('d',UnitMass_in_g))
		f.write(struct.pack('d',UnitVelocity_in_cm_per_s))
		f.write(struct.pack('i',dummy))


		dummy = nPart * (2 * size_i + 2 * size_I + 18 * size_f + 1 * size_d)
		f.write(struct.pack('i',dummy))

		for s in np.arange(nSnap):
			uint_buffer[:] = ID[:, s]
			f.write(struct.pack('{:d}I'.format(nPart), *uint_buffer))			

			float_buffer[:] = time[:, s]
			f.write(struct.pack('{:d}d'.format(nPart), *float_buffer))

			uint_buffer[:] = ParentCellID[:, s]
			f.write(struct.pack('{:d}I'.format(nPart), *uint_buffer))			

			float_buffer[:] = TracerMass[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))			

			float_buffer[:] = x[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = y[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = z[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = n_gas[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = temp[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = u_therm[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			# parameters for cooling
			float_buffer[:] = B[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = eps_photon[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			# parameters for diffusive shock acceleration
			int_buffer[:] = ShockFlag[:, s]
			f.write(struct.pack('{:d}i'.format(nPart), *int_buffer))

			float_buffer[:] = eps_CRpShockInj[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))


			float_buffer[:] = n_gasPreShock[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = n_gasPostShock[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = VShock[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = timeShockCross[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = theta[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			# parameters for electron injection (apart from DSA given above)
			int_buffer[:] = CReInjection[:, s]
			f.write(struct.pack('{:d}i'.format(nPart), *int_buffer))

			float_buffer[:] = injRate[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = alphaInj[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = pInj[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))


			f.write(struct.pack('i',dummy))
			if s < nSnap - 1:
				f.write(struct.pack('i',dummy))            

		f.close()

	return 0

####################################################################################################
# Writes a binary file with the initial spectrum data

def writeInitialSpectrumFile(fname, nPart, nBins, f, sameSpectra=False):
	# if sameSpectra = True, we assume a 1D momentum array for f which is the same for all particles
	size_i, size_I, size_f, size_d = checkNumberEncoding()

	if sameSpectra:
		if f.shape == (nPart, nBins):
			sys.exit("Set sameSpectra=False in function call for different spectra for every particle")
		if f.shape != (nBins, ):
			sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, )".format('f', f.shape, nBins))
	else:
		if f.shape != (nPart, nBins):
			sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('f', f.shape, nPart, nBins))	

	float_buffer = np.ndarray(nBins,dtype=float)

	with open(fname,'wb') as file:

		# first block with basic information
		file.write(struct.pack('i',2 * size_i))
		file.write(struct.pack('i',nPart))
		file.write(struct.pack('i',nBins))
		file.write(struct.pack('i',2 * size_i))

		# second block with actual data
		file.write(struct.pack('i', nPart * nBins * size_d))
		if sameSpectra:
			for i in np.arange(nPart):
				file.write(struct.pack('{:d}d'.format(nBins), *f))
		else:
			for i in np.arange(nPart):
				file.write(struct.pack('{:d}d'.format(nBins), *f[i]))

		file.write(struct.pack('i',nPart * nBins * size_d))
		file.close()
		
		print("Wrote initial spectrum for {:d} particles with {:d} momentum bins file to '{:s}'".format(nPart, nBins, fname))

	return 0

####################################################################################################
# Reads a binary file with the initial spectrum data

def readInitialSpectrumFile(fname, nPartIn=None, nBinsIn=None, sameSpectra=False):

	size_i, size_I, size_f, size_d = checkNumberEncoding()
	with open(fname,'rb') as file:
		print("Read initial spectrum for tracer particles")

		# Read first information block with number of particles and momentum bins
		dummy = int(struct.unpack('i',file.read(size_i))[0])
		if dummy != 2 * size_i:
			sys.exit("Block size is {:d}, but expected {:d}".format(dummy, 2 * size_i))

		nPart = int(struct.unpack('i',file.read(size_i))[0])
		if nPartIn != None :
			if nPart != nPartIn:
				sys.exit("nPart is {:d}, but expected {:d}".format(nPart, nPartIn))

		nBins = int(struct.unpack('i',file.read(size_i))[0])
		if nBinsIn != None:
			if nBins != nBinsIn:
				sys.exit("nBins is {:d}, but expected {:d}".format(nBins, nBinsIn))

		dummy = int(struct.unpack('i',file.read(size_i))[0])
		if dummy != 2 * size_i:
			sys.exit("1st data block not correctly enclosed")

		# if we assume that all particles contain the same spectrum we just read the spectrum of the first particle and we're done
		if sameSpectra:
			f = np.ndarray(nBins, dtype=float)
		else:
			f = np.ndarray((nPart, nBins), dtype=float)

		# Read now the actual spectral data
		dummy = int(struct.unpack('i',file.read(size_i))[0])
		if dummy != nPart * nBins * size_d:
			sys.exit("Block size is {:d} dummy, but expected {:d}".format(dummy, nPart * nBins * size_d))

		if sameSpectra:
			f = np.copy(struct.unpack('{:d}d'.format(nBins), file.read(size_d * nBins)))
			file.seek((nPart - 1) * size_d * nBins, 1) 
		else:
			for i in np.arange(nPart):
				f[i] = np.copy(struct.unpack('{:d}d'.format(nBins), file.read(size_d * nBins)))

		dummy = int(struct.unpack('i',file.read(size_i))[0])
		if dummy != nPart * nBins * size_d:
			sys.exit("1st data block not correctly enclosed")

		file.close()

	return f
