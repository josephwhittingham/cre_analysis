import matplotlib.pyplot as plt
import numpy as np
import struct


######################################################################################
# class which handles the parameters given via files to the C program
# parameters are needed for calculating the plots
class CRelectronParameters:
	# instance variables

	modeList=['CN: A(n+1)F(n+1) = B(n) F(n)','CN: A(n+0.5)F(n+1) = B(n+0.5) F(n)','CN: A(n+1)F(n+1) = B(n+1) F(n)','Explicit: F(n+1) = B(n) F(n)', 'Donor-Cell','Lax-Wendroff','Beam-Warming','Fromm','SL Minmod','SL Superbee','SL MC','SL van Leer']
	def __init__(self,ParameterFileName = None):
		# need to set dummy values as these determine the types

		# General Settings for I/O
		self.InputDataFile               = ''  
		self.Run_ID                      = 0  # obsolete
		self.OutputDir                   = ''
		self.SnapshotFileBase            = ''

		# Settings for Discretization
		self.NumberOfMomentumBins        = 0
		self.CourantFac                  = 0.0
		self.AlphaCoefficientMaximum     = 0.0
		self.MinimumMomentum             = 0.0
		self.MaximumMomentum             = 0.0
		self.TimeStepsUpdate             = 1

		# CR density distribution function
		self.AlphaSpectralIndex          = 0.0
		self.MomentumLowCutoff           = 0.0
		self.MomentumHighCutoff          = 0.0
		self.NormalizationFactor         = 0.0

		# Flags
		self.FlagLogGrid                 = 1 # obsolete
		self.FlagAllowSubcycles          = 1
		self.FlagSolverForAdvection      = 12 # obsolete
		self.mode                        = '' # this is directly given by the parameter file, but depends on FlagSolverForAdvection
		self.FlagProtons                 = 0 # obsolete
		self.Flag3DFunction              = 0 # obsolete

		# Cooling
		self.FlagCooling                 = 1
		self.n_elec                      = 1.157
		self.HydrogenMassFrac            = 0.76

		# Source Parameters (ALL obsolete)
		self.FlagSourceFunction          = 0
		self.SourceLowCutoff             = 0.0
		self.SourceHighCutoff            = 0.0
		self.SourceSpectralIndex         = 0.0
		self.SourceNormalization         = 0.0

		# parameters for shock injection
		self.ShockParamA                 = 0.
		self.ShockParamB                 = 0.
		self.ShockParamC                 = 0.

		# new parameters
		self.NumberOfTracerParticles     = 1
		self.UseSemiAnalyticSolutionLow  = 0
		self.UseSemiAnalyticSolutionHigh = 0
		self.FunctionValueChop           = 1.e-30

		if ParameterFileName is not None:
			self.read_data(ParameterFileName)

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
	def read_data(self,ParameterFileName):
		fParam = open(ParameterFileName,'r')
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
									print '\t{:25} {:}'.format(columnParam[0],columnParam[1])
									continue
								elif type(getattr(self, var)) is float:
									setattr(self,var,float(columnParam[1]))
									print '\t{:25} {:}'.format(columnParam[0],columnParam[1])
									continue
								elif type(getattr(self, var)) is str:
									setattr(self,var,columnParam[1])
									print '\t{:25} {:}'.format(columnParam[0],columnParam[1])
									continue
		if self.OutputDir[-1] != '/':
			self.OutputDir += '/'
		self.mode = self.modeList[self.FlagSolverForAdvection - 1]
		print '\n'
		line = None
		lineParam = None
		columnParam = None
		fParam.close()

######################################################################################
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
		self.rho    = np.zeros(self.NumberOfSnapshotsMax)
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
		print "\t time (s) \t rho (cm^-3) \t B (G) \n"        
		#print "\t time \t\t pos_x \t\t pos_y \t\t pos_z \t\t rho \t\t u \t\t vel_x \t\t vel_y \t\t vel_z \t\t dedt \n"
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.time[self.NumberOfSnapshots]     = float(columnData[0])
					#self.pos_x[self.NumberOfSnapshots]    = float(columnData[1])
					#self.pos_y[self.NumberOfSnapshots]     = float(columnData[2])
					#self.pos_z[self.NumberOfSnapshots]     = float(columnData[3])
					self.rho[self.NumberOfSnapshots]     = float(columnData[4])
					#self.u[self.NumberOfSnapshots]         = float(columnData[5])
					#self.vel_x[self.NumberOfSnapshots]     = float(columnData[6])
					#self.vel_y[self.NumberOfSnapshots]     = float(columnData[7])
					#self.vel_z[self.NumberOfSnapshots]     = float(columnData[8])
					#self.dedt[self.NumberOfSnapshots]     = float(columnData[9])
					self.B[self.NumberOfSnapshots]     = float(columnData[10])
					#print "\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[self.NumberOfSnapshots], self.pos_x[self.NumberOfSnapshots], self.pos_y[self.NumberOfSnapshots], self.pos_z[self.NumberOfSnapshots], self.rho[self.NumberOfSnapshots], self.u[self.NumberOfSnapshots], self.vel_x[self.NumberOfSnapshots], self.vel_y[self.NumberOfSnapshots], self.vel_z[self.NumberOfSnapshots], self.dedt[self.NumberOfSnapshots])
					print "\t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[self.NumberOfSnapshots], self.rho[self.NumberOfSnapshots],self.B[self.NumberOfSnapshots])
					self.NumberOfSnapshots += 1
		self.time.resize(self.NumberOfSnapshots)
		self.rho.resize(self.NumberOfSnapshots)
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
		print "\t time (s) \t rho (cm^-3) \t B (G) \n"
		for i in np.arange(NumberOfSnapshots):
			print "\t {:1.2E} \t {:1.2E} \t {:1.2E}".format(self.time[i], self.rho[i],self.B[i])
		print ''

######################################################################################
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
		self.rho    = np.ndarray(self.NumberOfSnapshotsMax,dtype=float)
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
		print "\t time (s) \t rho (cm^-3) \t T (K) \t B (G) \t inj \n"        
		#print "\t time \t\t pos_x \t\t pos_y \t\t pos_z \t\t rho \t\t u \t\t vel_x \t\t vel_y \t\t vel_z \t\t dedt \n"
		for line in fData:
			lineData = (line.strip()).lstrip()
			if(lineData != ''): #ignore empty lines
				if(lineData != '%' and lineData != '#'): # ignore the symbols %,#
					columnData = line.split()
					self.time[self.NumberOfSnapshots]     = float(columnData[0])
					self.rho[self.NumberOfSnapshots]     = float(columnData[1])
					self.T[self.NumberOfSnapshots]         = float(columnData[2])
					self.B[self.NumberOfSnapshots]         = float(columnData[3])
					self.inj[self.NumberOfSnapshots]     = int(columnData[4])
					print "\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:d}".format(self.time[self.NumberOfSnapshots], self.rho[self.NumberOfSnapshots],self.T[self.NumberOfSnapshots], self.B[self.NumberOfSnapshots], self.inj[self.NumberOfSnapshots])
					self.NumberOfSnapshots += 1
		self.time.resize(self.NumberOfSnapshots)
		self.rho.resize(self.NumberOfSnapshots)
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
		print "\t time (s) \t rho (cm^-3) \t T (K) \t B (G) \t inj \n"        
		for i in np.arange(NumberOfSnapshots):
			print "\t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:1.2E} \t {:d}".format(self.time[self.NumberOfSnapshots], self.rho[self.NumberOfSnapshots],self.T[self.NumberOfSnapshots], self.B[self.NumberOfSnapshots], self.inj[self.NumberOfSnapshots])
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



######################################################################################
# Simple Class to Read Error File
class runErrorRelative:
	nPointsMax = 30
	def __init__(self,runErrorFile = None):
		self.nPoints = 0
		self.time    = np.zeros(self.nPointsMax)
		self.rho    = np.zeros(self.nPointsMax)
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
					self.rho[self.nPoints]    = columnData[1]
					self.err[self.nPoints]    = columnData[2]
					print '\t {:E} \t {:E} \t {:E}'.format(self.time[self.nPoints], self.rho[self.nPoints], self.err[self.nPoints])
					self.nPoints+=1
		self.time.resize(self.nPoints)
		self.rho.resize(self.nPoints)
		self.err.resize(self.nPoints)

		print '\n'
		line = None
		lineData = None
		columnData = None
		fData.close()

######################################################################################
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

######################################################################################
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

######################################################################################
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


######################################################################################
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

######################################################################################
# Class to read the tracer output data
# takes the binary output file of arepo to extract data
# gets directly number of snapshots
# needs packages: struct
		
class TracerOutput:
	# instance variables
	def __init__(self, fname = None, with_cr_electrons=1):
		# with_cr_electrons is set to 1 if arepo was compiled with #COSMIC_RAYS_ELECTRONS
		# need to set dummy values as these determine the types
		checkNumberEncoding()
		self.nSnap = 99
		self.nPart = 0
		if fname is not None:
			self.read_data(fname,with_cr_electrons)

	def __del__(self):
		for var in vars(self):
			setattr(self,var,None)

	def read_data(self, fname,with_cr_electrons=1):
		with open(fname,'rb') as f:
			print "Read Arepo's tracer output from file '{}'".format(fname)
			size_i = struct.calcsize('i')
			size_f = struct.calcsize('f')
			size_d = struct.calcsize('d')

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
			self.x	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.y	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.z	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.rho   	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.temp	= np.ndarray((self.nPart, self.nSnap), dtype=float)
			self.Utherm	= np.ndarray((self.nPart, self.nSnap), dtype=float)

			if with_cr_electrons:
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
				self.rho[:, n]		= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.temp[:, n]		= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
				self.Utherm[:, n]	= struct.unpack('{:d}f'.format(self.nPart), f.read(size_f * self.nPart))
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


	


####################################################################################
# Function to check whether system encoding of int, float and double is correct
# returns 0 if everything is alright
def checkNumberEncoding():
	error = [0,0,0]

	# check integers
	size = struct.calcsize('i')
	if size != 4:
		error[0] = 1

	# check single precision floats
	size = struct.calcsize('f')
	if size !=4:
		error[1] = 1

	# check double precision floats
	size = struct.calcsize('d')
	if size !=8:
		error[2] = 1

	if sum(error) > 0 :
		sys.exit("Data types ({}{}{}{}{}) not correctly encoded on this machine!".format(
			["", "int"][error[0]==1],
			["", ", "][error[0]==1 and error[1]==1],
			["", "float"][error[1]==1],
			["", ", "][(error[0] == 1 and error[2]==1) or (error[1] == 1 and error[2] == 1)],
			["", "double"][error[2]==1]))

	else:
		return 0
