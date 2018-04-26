import matplotlib.pyplot as plt
import numpy as np
import struct
import sys
from Physics import *


####################################################################################################
# class which handles the parameters given via files to the C program
# parameters are needed for calculating the plots
class CRelectronParameters:
	# instance variables

	def __init__(self,ParameterFileName = None, RunID = None, verbose = False):
		# need to set dummy values as these determine the types
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
		self.n_elec                      = 1.157
		self.HydrogenMassFrac            = 0.76
		self.DiffusionTimeInGyr          = 0.

		# parameters for shock injection
		self.ShockParamA                 = 0.
		self.ShockParamB                 = 0.
		self.ShockParamC                 = 0.

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
				print '{:25} {:d}'.format(var,getattr(self,var))
			if type(getattr(self, var)) is float:
				print '{:25} {:.5e}'.format(var,getattr(self,var))
			if type(getattr(self, var)) is str:
				print '{:25} {:}'.format(var,getattr(self,var))
	
	# read in the parameter file and set the private class variables accordingly
	def read_data(self,ParameterFileName, RunID = None, verbose = False):
		fParam = open(ParameterFileName,'r')
		if verbose:
			print "Reading parameters from file '{:}'\n".format(ParameterFileName)
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
										print '\t{:25} {:}'.format(columnParam[0],columnParam[1])
									continue
								elif type(getattr(self, var)) is float:
									setattr(self,var,float(columnParam[1]))
									if verbose:
										print '\t{:25} {:}'.format(columnParam[0],columnParam[1])
									continue
								elif type(getattr(self, var)) is str:
									setattr(self,var,columnParam[1])
									if verbose:
										print '\t{:25} {:}'.format(columnParam[0],columnParam[1])
									continue
		if self.OutputDir[-1] != '/':
			self.OutputDir += '/'
		if verbose:
			print '\n'
		if RunID != None:
			self.Run_ID = RunID
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
		print "Reading data from file '{:}'\n".format(DataFileName)
		print "\t time (s) \t n_gas (cm^-3) \t B (G) \n"        
		#print "\t time \t\t pos_x \t\t pos_y \t\t pos_z \t\t n_gas \t\t u \t\t vel_x \t\t vel_y \t\t vel_z \t\t dedt \n"
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
					#print "\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[self.NumberOfSnapshots], self.pos_x[self.NumberOfSnapshots], self.pos_y[self.NumberOfSnapshots], self.pos_z[self.NumberOfSnapshots], self.n_gas[self.NumberOfSnapshots], self.u[self.NumberOfSnapshots], self.vel_x[self.NumberOfSnapshots], self.vel_y[self.NumberOfSnapshots], self.vel_z[self.NumberOfSnapshots], self.dedt[self.NumberOfSnapshots])
					print "\t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[self.NumberOfSnapshots], self.n_gas[self.NumberOfSnapshots],self.B[self.NumberOfSnapshots])
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


		print "\n"
		line = None
		lineData = None
		columnData = None
		fData.close()

	def show(self):
		print "Number of Snapshots: {:d}".format(self.NumberOfSnapshots)
		print "\t time (s) \t n_gas (cm^-3) \t B (G) \n"
		for i in np.arange(NumberOfSnapshots):
			print "\t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[i], self.n_gas[i],self.B[i])
		print ''

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
		print "Reading data from file '{:}'\n".format(DataFileName)
		print "\t time (s) \t n_gas (cm^-3) \t T (K) \t B (G) \t inj \n"        
		#print "\t time \t\t pos_x \t\t pos_y \t\t pos_z \t\t n_gas \t\t u \t\t vel_x \t\t vel_y \t\t vel_z \t\t dedt \n"
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
					print "\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:d}".format(self.time[self.NumberOfSnapshots], self.n_gas[self.NumberOfSnapshots],self.T[self.NumberOfSnapshots], self.B[self.NumberOfSnapshots], self.inj[self.NumberOfSnapshots])
					self.NumberOfSnapshots += 1
		self.time.resize(self.NumberOfSnapshots)
		self.n_gas.resize(self.NumberOfSnapshots)
		self.T.resize(self.NumberOfSnapshots)
		self.B.resize(self.NumberOfSnapshots)
		self.inj.resize(self.NumberOfSnapshots)

		print "\n"
		line = None
		lineData = None
		columnData = None
		fData.close()

	def show(self):
		print "Number of Snapshots: {:d}".format(self.NumberOfSnapshots)
		print "\t time (s) \t n_gas (cm^-3) \t T (K) \t B (G) \t inj \n"        
		for i in np.arange(NumberOfSnapshots):
			print "\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:d}".format(self.time[self.NumberOfSnapshots], self.n_gas[self.NumberOfSnapshots],self.T[self.NumberOfSnapshots], self.B[self.NumberOfSnapshots], self.inj[self.NumberOfSnapshots])
		print ''

	
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
		print "Reading data from file '{:}'".format(runErrorFile)
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
					print '\t {:E} \t {:E} \t {:E}'.format(self.time[self.nPoints], self.n_gas[self.nPoints], self.err[self.nPoints])
					self.nPoints+=1
		self.time.resize(self.nPoints)
		self.n_gas.resize(self.nPoints)
		self.err.resize(self.nPoints)

		print '\n'
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
		print "Reading data from file '{:}'".format(runErrorFile)
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
					print '\t {:E} \t {:E} \t {:E}'.format(self.time[self.nPoints], self.del_energy[self.nPoints], self.del_pressure[self.nPoints], self.del_number[self.nPoints])
					self.nPoints+=1
		self.time.resize(self.nPoints)
		self.del_energy.resize(self.nPoints)
		self.del_pressure.resize(self.nPoints)
		self.del_number.resize(self.nPoints)

		print '\n'
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
		print "Reading approximate solution from file '{:}'".format(approxFile)
		fData = open(approxFile,'r')
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.time[self.nSnaps] = float(columnData[0])
					self.Norm[self.nSnaps] = float(columnData[1])
					self.Cut[self.nSnaps]  = float(columnData[2])
					print '\t {:E} \t {:E} \t {:E}'.format(self.time[self.nSnaps], self.Norm[self.nSnaps], self.Cut[self.nSnaps])
					self.nSnaps += 1
		self.time.resize(self.nSnaps)
		self.Norm.resize(self.nSnaps)
		self.Cut.resize(self.nSnaps)

		print '\n'
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
		print "Reading snapshot data from file '{:}'".format(SnapshotDataFile)
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
		print "Reading snapshot data from file '{:}'".format(fname)

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

		nBins = dummy / size_d
	
		f = np.array(struct.unpack('{:d}d'.format(nBins), file.read(size_d * nBins)))

		dummy = int(struct.unpack('i',file.read(size_i))[0])
		if dummy != nBins * size_d:
			sys.exit("1st data block not correctly enclosed")

		file.close()

	return f, id, time


####################################################################################################
#
class cre_snapshot:
	def __init__(self, fname = None, verbose = False):
		if fname is not None:
			self.read_data(fname, verbose)

	def read_data(self, fname, verbose = False):
		size_i, size_I, size_f, size_d = checkNumberEncoding()
		with open(fname,'rb') as file:
			if(verbose):
				print "Reading snapshot data from file '{:}'".format(fname)

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

			# Spectrum Data
			dummy = int(struct.unpack('i', file.read(size_i))[0])
			blocksize = self.nPart * ( self.nBins * size_d + size_I)
			if dummy != blocksize:
				sys.exit("Block size is {:d} bytes, but expexted {:d}".format(dummy, blocksize))

			self.f = np.ndarray((self.nPart, self.nBins), dtype = float)
			self.id = np.ndarray(self.nPart, dtype=np.uint32)
			for i in np.arange(self.nPart):
				self.id[i] = struct.unpack('I', file.read(size_I))[0]
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
		print "Reading other solution data from file '{:}'".format(DataFile)
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
# Class to read the tracer output data
# takes the binary output file of arepo to extract data
# gets directly number of snapshots
# needs packages: struct
		
class TracerOutput:
	# instance variables
	def __init__(self, fname = None, cgs_units = False, verbose = False):
		# with_cr_electrons is set to 1 if arepo was compiled with #COSMIC_RAYS_ELECTRONS
		# need to set dummy values as these determine the types
		self.nSnap = 0
		self.nPart = 0
		self.All_Units_in_cgs = False
		self.UnitLength_in_cm = 1.
		self.UnitMass_in_g = 1.
		self.UnitVelocity_in_cm_per_s = 1.

		if fname is not None:
			self.read_data(fname, cgs_units, verbose)


	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def read_data(self, fname, cgs_units = False, param = None, verbose = False, UnitLength_in_cm = 1., UnitMass_in_g = 1., UnitVelocity_in_cm_per_s = 1.):
		with open(fname,'rb') as f:
			if verbose:
				print "Read Arepo's tracer output from file '{}'".format(fname)
			size_i, size_I, size_f, size_d = checkNumberEncoding()

			# Reading first block with unit system
			dummy = int(struct.unpack('i',f.read(size_i))[0])
			if dummy != 3 * size_d:
				sys.exit("Expected 3 double values at beginning, ")

			self.UnitLength_in_cm         = struct.unpack('d', f.read(size_d))[0]
			self.UnitMass_in_g            = struct.unpack('d', f.read(size_d))[0]
			self.UnitVelocity_in_cm_per_s = struct.unpack('d', f.read(size_d))[0]

			if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
				sys.exit("Expected 3 double values at beginning, ")
			
			
			# Reading block with data values
			dummy = int(struct.unpack('i',f.read(size_i))[0])
			self.nPart = dummy / (2 * size_i + 2 * size_I + 18 * size_f + 1 * size_d)
			self.nSnap	= 0
			buf 	= 1

			while(buf):
				# move pointer forward
				f.seek(self.nPart * (2 * size_i + 2 * size_I + 19 * size_f + 1 * size_d), 1) 

				if  int(struct.unpack('i',f.read(size_i))[0]) != dummy:
					sys.exit("data not correctly enclosed 1, ")

				self.nSnap += 1
				buf = f.read(size_i)

			# go back to the beginning of the file
			f.seek(3*size_i + 3*size_d, 0)
			buf = 0
			if verbose:
				print 'Number of particles: {:d}'.format(self.nPart)
				print 'Number of snapshots: {:d}'.format(self.nSnap)

			# create the arrays
			self.ID             = np.ndarray((self.nPart, self.nSnap), dtype=np.uint32)
			self.time           = np.ndarray((self.nPart, self.nSnap), dtype=float) # time or scale parameter
			self.ParentCellID   = np.ndarray((self.nPart, self.nSnap), dtype=np.uint32)
			self.TracerDensity  = np.ndarray((self.nPart, self.TracerDensity), dtype=float)

			self.x	            = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.y	            = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.z	            = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.n_gas	        = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.temp	        = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.u_therm        = np.ndarray((self.nPart, self.nSnap), dtype=float)

			# parameters for cooling
			self.B			    = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.u_photon       = np.ndarray((self.nPart, self.nSnap), dtype=float)

			# parameters for diffusive shock acceleration
			self.ShockFlag      = np.ndarray((self.nPart, self.nSnap), dtype=int)
			self.n_gasPreShock  = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.n_gasPostShock = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.BpreShock      = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.BpostShock     = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.VShock         = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.timeShockCross = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.cosTheta       = np.ndarray((self.nPart, self.nSnap), dtype=float)
			
			# parameters for electron injection (apart from DSA given above)
			self.CReInjection   = np.ndarray((self.nPart, self.nSnap), dtype=int)
			self.injRate        = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.alphaInj       = np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.pInj           = np.ndarray((self.nPart, self.nSnap), dtype=float)

			# read the data
			for n in np.arange(self.nSnap):
				self.ID[:, n]             = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				self.time[:, n]		      = struct.unpack('{:d}d'.format(self.nPart), f.read(size_d * self.nPart))
				self.ParentCellID[:, n]   = struct.unpack('{:d}I'.format(self.nPart), f.read(size_I * self.nPart))
				self.TracerDensity[:, n]  = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))

				self.x[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.y[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.z[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.n_gas[:, n]		  = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.temp[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.u_therm[:, n]	      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))

				self.B[:, n]		      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.u_photon[:, n]       = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.ShockFlag[:, n]      = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
				self.n_gasPreShock[:, n]  = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.n_gasPostShock[:, n] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))

				self.BpreShock[:, n]      = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.BpostShock[:, n]     = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.VShock[:, n]         = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.timeShockCross[:, n] = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.cosTheta[:, n]       = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))

				self.CReInjection[:, n]   = struct.unpack('{:d}i'.format(self.nPart), f.read(size_i * self.nPart))
				self.injRate[:, n]        = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.alphaInj[:, n]       = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.pInj[:, n]           = struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))


				if  int(struct.unpack('i', f.read(size_i))[0]) != dummy:
					sys.exit("data not correctly enclosed, ")

				f.seek(size_i, 1)

			f.close()
			if verbose:
				print "Data was successfully read"
			
			if cgs_units:
				self.scale_to_cgs(verbose)

	def scale_to_cgs(self, verbose=False):	
		if verbose:
			print "Scale to cgs with UnitLenght_in_cm = {:.3e}, UnitMass_in_g = {:.3e}, UnitVeloctiy_in_cm_per_s = {:.3e}".format(self.UnitLength_in_cm, self.UnitMass_in_g, self.UnitVelocity_in_cm_per_s)
		self.time           = np.multiply(self.time, self.UnitLength_in_cm / self.UnitVelocity_in_cm_per_s)
		self.TracerDensity  = np.multiply(self.TracerDensity, np.power(self.UnitLength_in_cm, -3.))
		self.x              = np.multiply(self.x, self.UnitLength_in_cm)
		self.y              = np.multiply(self.y, self.UnitLength_in_cm)
		self.z              = np.multiply(self.z, self.UnitLength_in_cm)
		self.n_gas          = np.multiply(self.n_gas, self.UnitMass_in_g / (PROTONMASS * np.power(self.UnitLength_in_cm, 3)))
		self.temp           = np.multiply(self.temp, np.square(self.UnitVelocity_in_cm_per_s))
		self.u_therm        = np.multiply(self.u_therm, self.UnitMass_in_g * np.square(self.UnitVelocity_in_cm_per_s) / np.power(self.UnitLength_in_cm, 3))

		self.B              = np.multiply(self.B, np.sqrt(self.UnitMass_in_g / self.UnitLength_in_cm) * self.UnitVelocity_in_cm_per_s / self.UnitLength_in_cm)
		self.u_photon       = np.multiply(self.u_photon, self.UnitMass_in_g * np.square(self.UnitVelocity_in_cm_per_s) / np.power(self.UnitLength_in_cm, 3))

		self.n_gasPreShock  = np.multiply(self.n_gasPreShock,  self.UnitMass_in_g / (PROTONMASS * np.power(self.UnitLength_in_cm, 3)))	
		self.n_gasPostShock = np.multiply(self.n_gasPostShock, self.UnitMass_in_g / (PROTONMASS * np.power(self.UnitLength_in_cm, 3)))

		self.BpreShock      = np.multiply(self.B, np.sqrt(self.UnitMass_in_g / self.UnitLength_in_cm) * self.UnitVelocity_in_cm_per_s / self.UnitLength_in_cm)
		self.BpostShock     = np.multiply(self.B, np.sqrt(self.UnitMass_in_g / self.UnitLength_in_cm) * self.UnitVelocity_in_cm_per_s / self.UnitLength_in_cm)
		self.VShock         = np.multiply(self.VShock, self.UnitVelocity_in_cm_per_s)		
		self.timeShockCross = np.multiply(self.timeShockCross, self.UnitLength_in_cm / self.UnitVelocity_in_cm_per_s)

		self.injRate        = np.multiply(self.injRate, self.UnitVelocity_in_cm_per_s / self.UnitLength_in_cm)

	def __getitem__(self, key):
		# check dimensions of return
		ret = TracerOutput()
		ret.All_Units_in_cgs = self.All_Units_in_cgs
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
					ret.nSnap = (stop - start + 1)/step
				else:
					raise TypeError('Index must be int or slice, not {}'.format(type(key[1]).__name__))
		else:
			raise TypeError('Tuple Index must be of length 2, not {}'.format(len(key)))	

		# create the arrays
		ret.ID             = self.ID.__getitem__(key)
		ret.time           = self.time.__getitem__(key)
		ret.ParentCellID   = self.ParentCellID.__getitem__(key)
		ret.TracerDensity  = self.TracerDensity.__getitem__(key)

		ret.x	            = self.x.__getitem__(key)
		ret.y	            = self.y.__getitem__(key)
		ret.z	            = self.z.__getitem__(key)
		ret.n_gas	        = self.n_gas.__getitem__(key)
		ret.temp	        = self.temp.__getitem__(key)
		ret.u_therm        = self.u_therm.__getitem__(key)

		# parameters for cooling
		ret.B			    = self.B.__getitem__(key)
		ret.u_photon       = self.u_photon.__getitem__(key)

		# parameters for diffusive shock acceleration
		ret.ShockFlag      = self.ShockFlag.__getitem__(key)
		ret.n_gasPreShock  = self.n_gasPreShock.__getitem__(key)
		ret.n_gasPostShock = self.n_gasPostShock.__getitem__(key)
		ret.BpreShock      = self.BpreShock.__getitem__(key)
		ret.BpostShock     = self.BpostShock.__getitem__(key)
		ret.VShock         = self.VShock.__getitem__(key)
		ret.timeShockCross = self.timeShockCross.__getitem__(key)
		ret.cosTheta       = self.cosTheta.__getitem__(key)

		# parameters for electron injection (apart from DSA given above)
		ret.CReInjection   = self.CReInjection.__getitem__(key)
		ret.injRate        = self.injRate.__getitem__(key)
		ret.alphaInj       = self.alphaInj.__getitem__(key)
		ret.pInj           = self.pInj.__getitem__(key)		

		return ret

			
		
		
			
	# 	a = self.ShockFlag.__getitem__(key)
	# 	#return a

	# 	if isinstance(key, int):
	# 		return self.getitem(key)
	# 	elif isinstance(key, slice):
	# 		return self.getslice(key)
	# 	elif isinstance(key, tuple):
	# 		if len(key) > 2:
	# 			raise TypeError('Tuple Index must be of length 2, not {}'.format(len(key)))
	# 		else:
	# 			ret = self.__getitem__(key[0])
	# 			if isinstance(key[1], int):
	# 				return ret.getsnapitem(key[1])
	# 			elif isinstance(key[1], slice):
	# 				return ret.getsnapslice(key[1])
	# 			else:
	# 				raise TypeError('Index must be int, not {}'.format(type(key[1]).__name__))
	# 	else:
	# 		raise TypeError('Index must be int, not {}'.format(type(key).__name__))

	# def getitem(self, key):
	# 	if key < self.nPart:
	# 		ret = TracerOutput()
	# 		ret.nPart = 1
	# 		ret.nSnap = self.nSnap

	# 		# create the arrays
	# 		ret.ID             = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 		ret.time           = np.ndarray((ret.nPart, ret.nSnap), dtype = float) # time or scale parameter

	# 		ret.x	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.y	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.z	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.n_gas	        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.temp	        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.u_therm        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 		# parameters for cooling
	# 		ret.B			    = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.u_photon       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 		# parameters for diffusive shock acceleration
	# 		ret.ShockFlag      = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 		ret.n_gasPreShock    = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.n_gasPostShock   = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.BpreShock      = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.BpostShock     = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.VShock      = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.timeShockCross = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.cosTheta       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
			
	# 		# parameters for electron injection (apart from DSA given above)
	# 		ret.CReInjection   = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 		ret.injRate        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.alphaInj       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.pInj           = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 		ret.ID[0, :]             = self.ID[key, :]
	# 		ret.time[0, :]		     = self.time[key, :]

	# 		ret.x[0, :]		         = self.x[key, :] 
	# 		ret.y[0, :]		         = self.y[key, :] 
	# 		ret.z[0, :]		         = self.z[key, :] 
	# 		ret.n_gas[0, :]		     = self.n_gas[key, :] 
	# 		ret.temp[0, :]		     = self.temp[key, :] 
	# 		ret.u_therm[0, :]	     = self.u_therm[key, :] 

	# 		ret.B[0, :]		         = self.B[key, :] 
	# 		ret.u_photon[0, :]       = self.u_photon[key, :] 

	# 		ret.ShockFlag[0, :]      = self.ShockFlag[key, :] 
	# 		ret.n_gasPreShock[0, :]    = self.n_gasPreShock[key, :] 
	# 		ret.n_gasPostShock[0, :]   = self.n_gasPostShock[key, :] 
	# 		ret.BpreShock[0, :]      = self.BpreShock[key, :] 
	# 		ret.BpostShock[0, :]     = self.BpostShock[key, :] 
	# 		ret.VShock[0, :]      = self.VShock[key, :] 
	# 		ret.timeShockCross[0, :] = self.timeShockCross[key, :] 
	# 		ret.cosTheta[0, :]       = self.cosTheta[key, :] 

	# 		ret.CReInjection[0, :]   = self.CReInjection[key, :] 
	# 		ret.injRate[0, :]        = self.injRate[key, :] 
	# 		ret.alphaInj[0, :]       = self.alphaInj[key, :] 
	# 		ret.pInj[0, :]           = self.pInj[key, :] 

	# 		return ret

	# def getslice(self, key):
	# 	start, stop, step = key.indices(self.nPart)
	# 	ret = TracerOutput()
	# 	ret.nPart = (stop - start + 1)/step
	# 	ret.nSnap = self.nSnap

	# 	# create the arrays
	# 	ret.ID             = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 	ret.time           = np.ndarray((ret.nPart, ret.nSnap), dtype = float) # time or scale parameter

	# 	ret.x	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.y	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.z	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.n_gas	        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.temp	        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.u_therm        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 	# parameters for cooling
	# 	ret.B			    = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.u_photon       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 	# parameters for diffusive shock acceleration
	# 	ret.ShockFlag      = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 	ret.n_gasPreShock    = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.n_gasPostShock   = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.BpreShock      = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.BpostShock     = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.VShock      = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.timeShockCross = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.cosTheta       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 	# parameters for electron injection (apart from DSA given above)
	# 	ret.CReInjection   = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 	ret.injRate        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.alphaInj       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.pInj           = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 	# Copy now the data
	# 	j = 0
	# 	for i in np.arange(start, stop, step):
	# 		ret.ID[j, :]             = self.ID[i, :]
	# 		ret.time[j, :]		     = self.time[i, :]

	# 		ret.x[j, :]		         = self.x[i, :] 
	# 		ret.y[j, :]		         = self.y[i, :] 
	# 		ret.z[j, :]		         = self.z[i, :] 
	# 		ret.n_gas[j, :]		     = self.n_gas[i, :] 
	# 		ret.temp[j, :]		     = self.temp[i, :] 
	# 		ret.u_therm[j, :]	     = self.u_therm[i, :] 

	# 		ret.B[j, :]		         = self.B[i, :] 
	# 		ret.u_photon[j, :]       = self.u_photon[i, :] 

	# 		ret.ShockFlag[j, :]      = self.ShockFlag[i, :] 
	# 		ret.n_gasPreShock[j, :]    = self.n_gasPreShock[i, :] 
	# 		ret.n_gasPostShock[j, :]   = self.n_gasPostShock[i, :] 
	# 		ret.BpreShock[j, :]      = self.BpreShock[i, :] 
	# 		ret.BpostShock[j, :]     = self.BpostShock[i, :] 
	# 		ret.VShock[j, :]      = self.VShock[i, :] 
	# 		ret.timeShockCross[j, :] = self.timeShockCross[i, :] 
	# 		ret.cosTheta[j, :]       = self.cosTheta[i, :] 

	# 		ret.CReInjection[j, :]   = self.CReInjection[i, :] 
	# 		ret.injRate[j, :]        = self.injRate[i, :] 
	# 		ret.alphaInj[j, :]       = self.alphaInj[i, :] 
	# 		ret.pInj[j, :]           = self.pInj[i, :]

	# 		j += 1

	# 	j = None
	# 	start = None
	# 	stop = None
	# 	step = None

	# 	return ret

	# def getsnapitem(self, key):
	# 	if key < self.nSnap:
	# 		ret = TracerOutput()
	# 		ret.nPart = self.nPart
	# 		ret.nSnap = 1

	# 		# create the arrays
	# 		ret.ID             = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 		ret.time           = np.ndarray((ret.nPart, ret.nSnap), dtype = float) # time or scale parameter

	# 		ret.x	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.y	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.z	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.n_gas	        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.temp	        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.u_therm        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 		# parameters for cooling
	# 		ret.B			    = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.u_photon       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 		# parameters for diffusive shock acceleration
	# 		ret.ShockFlag      = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 		ret.n_gasPreShock    = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.n_gasPostShock   = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.BpreShock      = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.BpostShock     = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.VShock      = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.timeShockCross = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.cosTheta       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
			
	# 		# parameters for electron injection (apart from DSA given above)
	# 		ret.CReInjection   = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 		ret.injRate        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.alphaInj       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 		ret.pInj           = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 		for p in np.arange(ret.nPart):
	# 			ret.ID[p, :]             = self.ID[p, key]
	# 			ret.time[p, :]		     = self.time[p, key]

	# 			ret.x[p, :]		         = self.x[p, key] 
	# 			ret.y[p, :]		         = self.y[p, key] 
	# 			ret.z[p, :]		         = self.z[p, key] 
	# 			ret.n_gas[p, :]		     = self.n_gas[p, key] 
	# 			ret.temp[p, :]		     = self.temp[p, key] 
	# 			ret.u_therm[p, :]	     = self.u_therm[p, key] 

	# 			ret.B[p, :]		         = self.B[p, key] 
	# 			ret.u_photon[p, :]       = self.u_photon[p, key] 

	# 			ret.ShockFlag[p, :]      = self.ShockFlag[p, key] 
	# 			ret.n_gasPreShock[p, :]    = self.n_gasPreShock[p, key] 
	# 			ret.n_gasPostShock[p, :]   = self.n_gasPostShock[p, key] 
	# 			ret.BpreShock[p, :]      = self.BpreShock[p, key] 
	# 			ret.BpostShock[p, :]     = self.BpostShock[p, key] 
	# 			ret.VShock[p, :]      = self.VShock[p, key] 
	# 			ret.timeShockCross[p, :] = self.timeShockCross[p, key] 
	# 			ret.cosTheta[p, :]       = self.cosTheta[p, key] 

	# 			ret.CReInjection[p, :]   = self.CReInjection[p, key] 
	# 			ret.injRate[p, :]        = self.injRate[p, key] 
	# 			ret.alphaInj[p, :]       = self.alphaInj[p, key] 
	# 			ret.pInj[p, :]           = self.pInj[p, key] 

	# 		return ret

	# def getsnapslice(self, key):
	# 	start, stop, step = key.indices(self.nSnap)
	# 	ret = TracerOutput()
	# 	ret.nPart = self.nPart
	# 	ret.nSnap = (stop - start + 1)/step

	# 	# create the arrays
	# 	ret.ID             = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 	ret.time           = np.ndarray((ret.nPart, ret.nSnap), dtype = float) # time or scale parameter

	# 	ret.x	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.y	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.z	            = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.n_gas	        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.temp	        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.u_therm        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 	# parameters for cooling
	# 	ret.B			    = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.u_photon       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 	# parameters for diffusive shock acceleration
	# 	ret.ShockFlag      = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 	ret.n_gasPreShock    = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.n_gasPostShock   = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.BpreShock      = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.BpostShock     = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.VShock      = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.timeShockCross = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.cosTheta       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 	# parameters for electron injection (apart from DSA given above)
	# 	ret.CReInjection   = np.ndarray((ret.nPart, ret.nSnap), dtype=int)
	# 	ret.injRate        = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.alphaInj       = np.ndarray((ret.nPart, ret.nSnap), dtype=float)
	# 	ret.pInj           = np.ndarray((ret.nPart, ret.nSnap), dtype=float)

	# 	for p in np.arange(ret.nPart):
	# 		# Copy now the data
	# 		j = 0
	# 		for i in np.arange(start, stop, step):
	# 			ret.ID[p, j]             = self.ID[p, i]
	# 			ret.time[p, j]		     = self.time[p, i]

	# 			ret.x[p, j]		         = self.x[p, i] 
	# 			ret.y[p, j]		         = self.y[p, i] 
	# 			ret.z[p, j]		         = self.z[p, i] 
	# 			ret.n_gas[p, j]		     = self.n_gas[p, i] 
	# 			ret.temp[p, j]		     = self.temp[p, i] 
	# 			ret.u_therm[p, j]	     = self.u_therm[p, i] 

	# 			ret.B[p, j]		         = self.B[p, i] 
	# 			ret.u_photon[p, j]       = self.u_photon[p, i] 

	# 			ret.ShockFlag[p, j]      = self.ShockFlag[p, i] 
	# 			ret.n_gasPreShock[p, j]    = self.n_gasPreShock[p, i] 
	# 			ret.n_gasPostShock[p, j]   = self.n_gasPostShock[p, i] 
	# 			ret.BpreShock[p, j]      = self.BpreShock[p, i] 
	# 			ret.BpostShock[p, j]     = self.BpostShock[p, i] 
	# 			ret.VShock[p, j]      = self.VShock[p, i] 
	# 			ret.timeShockCross[p, j] = self.timeShockCross[p, i] 
	# 			ret.cosTheta[p, j]       = self.cosTheta[p, i] 

	# 			ret.CReInjection[p, j]   = self.CReInjection[p, i] 
	# 			ret.injRate[p, j]        = self.injRate[p, i] 
	# 			ret.alphaInj[p, j]       = self.alphaInj[p, i] 
	# 			ret.pInj[p, j]           = self.pInj[p, i] 
	# 			j += 1

	# 	j = None
	# 	start = None
	# 	stop = None
	# 	step = None

	# 	return ret
		



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
			print "Read Arepo's tracer output from file '{}'".format(fname)
			size_i, size_I, size_f, size_d = checkNumberEncoding()
		
			dummy = int(struct.unpack('i',f.read(size_i))[0])
			if with_cr_electrons == 1:
				self.nPart = (dummy - size_d) / (12 * size_f + size_i)
			else:
				self.nPart = (dummy - size_d) / (6 * size_f)

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
			print 'Number of particles: {:d}'.format(self.nPart)
			print 'Number of snapshots: {:d}'.format(self.nSnap)

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
				self.u_photon   = np.ndarray((self.nPart, self.nSnap), dtype=float) # TODO: Read real data or calculate them
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

			print "Data was successfully read"


	


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
def writeTracerArepo(fileName, nSnap, nPart, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s, ID, time, ParentCellID, x, y, z, n_gas, temp, u_therm, B, u_photon, ShockFlag, n_gasPreShock, n_gasPostShock, BpreShock, BpostShock, VShock, timeShockCross, cosTheta, CReInjection, injRate, alphaInj, pInj):
	
	size_i, size_I, size_f, size_d = checkNumberEncoding()

	# do some consistency checks
	if ID.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('ID', ID.shape, nPart, nSnap))

	if time.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('time', time.shape, nPart, nSnap))

	if ParentCellID.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('ParentCellID', ParentCellID.shape, nPart, nSnap))

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

	if u_photon.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('u_photon', u_photon.shape, nPart, nSnap))

	if ShockFlag.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('ShockFlag', ShockFlag.shape, nPart, nSnap))

	if n_gasPreShock.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('n_gasPreShock', n_gasPreShock.shape, nPart, nSnap))

	if n_gasPostShock.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('n_gasPostShock', n_gasPostShock.shape, nPart, nSnap))

	if BpreShock.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('BpreShock', BpreShock.shape, nPart, nSnap))

	if BpostShock.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('BpostShock', BpostShock.shape, nPart, nSnap))

	if VShock.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('VShock', VShock.shape, nPart, nSnap))

	if timeShockCross.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('timeShockCross', timeShockCross.shape, nPart, nSnap))

	if cosTheta.shape != (nPart, nSnap):
		sys.exit("Dimensions error: shape of '{:s}' is {:s}, expected is ({:d}, {:d})".format('cosTheta', cosTheta.shape, nPart, nSnap))

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
	uint_buffer  = np.ndarray(nPart,dtype=uint32)

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

			float_buffer[:] = u_photon[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			# parameters for diffusive shock acceleration
			int_buffer[:] = ShockFlag[:, s]
			f.write(struct.pack('{:d}i'.format(nPart), *int_buffer))

			float_buffer[:] = n_gasPreShock[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = n_gasPostShock[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = BpreShock[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = BpostShock[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = VShock[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = timeShockCross[:, s]
			f.write(struct.pack('{:d}f'.format(nPart), *float_buffer))

			float_buffer[:] = cosTheta[:, s]
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

def writeInitialSpectrumFile(fname, nPart, nBins, f):
	size_i, size_I, size_f, size_d = checkNumberEncoding()

	float_buffer = np.ndarray(nBins,dtype=float)

	with open(fname,'wb') as file:

		# first block with basic information
		file.write(struct.pack('i',2 * size_i))
		file.write(struct.pack('i',nPart))
		file.write(struct.pack('i',nBins))
		file.write(struct.pack('i',2 * size_i))

		# second block with actual data
		file.write(struct.pack('i', nPart * nBins * size_d))
		for i in np.arange(nPart):
			file.write(struct.pack('{:d}d'.format(nBins), *f[i]))

		file.write(struct.pack('i',nPart * nBins * size_d))
		file.close()

	return 0

####################################################################################################
# Reads a binary file with the initial spectrum data

def readInitialSpectrumFile(fname, nPartIn=None, nBinsIn=None):

	size_i, size_I, size_f, size_d = checkNumberEncoding()
	with open(fname,'rb') as file:
		print "Read initial spectrum for tracer particles"

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


		f = np.ndarray((nPart, nBins), dtype=float)

		# Read now the actual spectral data
		dummy = int(struct.unpack('i',file.read(size_i))[0])
		if dummy != nPart * nBins * size_d:
			sys.exit("Block size is {:d} dummy, but expected {:d}".format(dummy, nPart * nBins * size_d))

		for i in np.arange(nPart):
			f[i] = struct.unpack('{:d}d'.format(nBins), file.read(size_d * nBins))

		dummy = int(struct.unpack('i',file.read(size_i))[0])
		if dummy != nPart * nBins * size_d:
			sys.exit("1st data block not correctly enclosed")

		file.close()

	return f
