import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
import argparse
import sys
import os
sys.path.append('../../CRelectron/Analysis/')
from ParamClass import *
from Physics import *
from matplotlib import rc

parser = argparse.ArgumentParser(description="Makes the plots of the provided run id")
parser.add_argument('Run_ID', type=int, help="pass a Run_ID into the script")
parser.add_argument('-np','--noplot', help="optional no plotting", action="store_true")
 
######################################################################################
# Function to make the plots
def makeplot(snapIni, n,param,SnapData, pGrid, noplots=False):
	
	#read current Distribution
	snapCurr = DistributionFunction('{:}{:}{:}tp000_{num:03d}.dat'.format(BaseDir,param.OutputDir,param.SnapshotFileBase,num=n),param.NumberOfMomentumBins)

	if not(param.Flag3DFunction == 0):
		# calculate to the effective 1D function
		snapCurr.f = np.multiply(4. * pi * np.power(snapCurr.p,2),snapCurr.f)

	fA = SolutionAdiabaticChanges(snap0.p, SnapData.rho[0], SnapData.rho[n], param)
	pL,pH = CutoffsAdiabaticChanges(snapCurr.p, SnapData.rho[0], SnapData.rho[n], param)
	

	#find the lowest bin which is above the cutoff
	imin = np.where(snapCurr.p	>=	pL * 50.)[0][0]
	imax = np.where(snapCurr.p	<=	pH / 50.)[0][-1]

	err1 = simps(np.absolute(np.ones(param.NumberOfMomentumBins)[imin:imax+1] - np.divide(fA[imin:imax+1],snapCurr.f[imin:imax+1])),snapCurr.p[imin:imax+1])/(snapCurr.p[imax+1] - snapCurr.p[imin])

	imin = 2
	imax = -2
	err2 = simps(np.absolute(np.ones(param.NumberOfMomentumBins)[imin:imax+1] - np.divide(fA[imin:imax+1],snapCurr.f[imin:imax+1])),snapCurr.p[imin:imax+1])/(snapCurr.p[imax] - snapCurr.p[imin])
		

	if not(noplots):
		#create a new figure
		fig = plt.figure(figsize=[12,8], facecolor='white')
		axes = fig.gca()

		#set x and y axis
		if param.FlagLogGrid == 1:
			axes.set_xscale('log')
		else:
			axes.set_xscale('linear')
		axes.set_yscale('log')

		#set the plot range in x and y direction
		axes.set_ylim([1E-6,1E1])
		axes.set_xlim([param.MinimumMomentum,param.MaximumMomentum])

		
		#simulated solution
		#axes.plot(snapCurr.p,np.multiply(snapCurr.f,np.power(snapCurr.p,param.AlphaSpectralIndex)),linestyle='-',color='b',linewidth = 2, label='simulated spectrum')
		axes.plot(snapCurr.p,SpectralEnergyDensity(snapCurr.p,snapCurr.f),linestyle='-',color='b',linewidth = 2, label='simulated spectrum')
		#analytical solution
		#axes.plot(snapIni.p,np.multiply(fA,np.power(snapIni.p,param.AlphaSpectralIndex)),linestyle='--', color='r', linewidth = 2, label='analytical solution')
		axes.plot(snapIni.p,SpectralEnergyDensity(snapCurr.p,fA),linestyle='--', color='r', linewidth = 2, label='analytical solution')


		plt.xlabel(r'momentum $p = \frac{{P}}{{m_\mathrm{{e}}c}}$',fontsize=16)
		if param.Flag3DFunction == 0:
			plt.ylabel(r'$f^{{\mathrm{{(1D)}}}} p^{{\alpha}}$',fontsize=16)
		else:
			plt.ylabel(r'$4 \pi p^2 f^{{\mathrm{{(3D)}}}} p^{{\alpha}}$'%param.AlphaSpectralIndex)
		#add a time and error label
		modesub=['no','yes']

		axes.text(0.4,0.88, 
			r"""Mode: {:}
$T = {:.1f}$
$\rho = {:.1f}$
$\delta= {:.2E}$""".format(param.mode,SnapData.time[n] ,SnapData.rho[n],err1), transform = fig.transFigure, fontsize=16, horizontalalignment='left', verticalalignment='top')

		#show important paramters
		axes.text(0.15,0.88, 
			r"""run {:2d}
$\alpha = {:.2f}$
$C_{{cfl}} = {:.2f}$
$N_{{p}} = {:d}$""".format(param.Run_ID, param.AlphaSpectralIndex,param.CourantFac,param.NumberOfMomentumBins), transform = fig.transFigure, fontsize=16, horizontalalignment='left', verticalalignment='top')


		#legend
		legend = plt.legend(loc=1)
		
		#save the figure
		fig.savefig('{:}{:}{:}{num:03d}.png'.format(BaseDir,param.OutputDir,'plt_func_',num=n),dpi=100)

		axes = None
		fig = None
		plt.close()
		



	snapCurr = None
	fA		= None

	return err1, err2

	
######################################################################################
# Main Function
if __name__ == "__main__":

	
	args = parser.parse_args()
	Run_ID = int(args.Run_ID)	
	
	BaseDir = '../' # main folder from where Input, Output, etc can accessed

	# Read now the parameter file
	param = CRelectronParameters('{:}Input/param_{:03d}.txt'.format(BaseDir,Run_ID))

	if(Run_ID != param.Run_ID):
		#print "Provided Run_ID is not the same as that given in the parameter file!"
		a=1
	else:
		# Read now the tracer data file
		snapData = SnapshotData2('{:}{:}'.format(BaseDir,param.InputDataFile))

		#read in inital distribution
		snap0 = DistributionFunction('{:}{:}{:}tp000_{num:03d}.dat'.format(BaseDir,param.OutputDir,param.SnapshotFileBase,num=0),param.NumberOfMomentumBins)

		if param.Flag3DFunction != 0:
			# calculate to the effective 1D function
			snap0.f = np.multiply(4. * np.pi * np.power(snap0.p,2),snap0.f)

		#Calculate the Momentum grid for simple integration
		pGrid = np.zeros(param.NumberOfMomentumBins +1)
		if param.FlagLogGrid == 1:
			del_p = (np.log(param.MaximumMomentum) - np.log(param.MinimumMomentum)) / param.NumberOfMomentumBins
			for i in range(param.NumberOfMomentumBins + 1):
				pGrid[i] = param.MinimumMomentum * np.exp(del_p * i)
		else:
			del_p = (param.MaximumMomentum - param.MinimumMomentum)/ param.NumberOfMomentumBins
			for i in range(param.NumberOfMomentumBins + 1):
				pGrid[i] = param.MinimumMomentum + del_p * i

		err1 = np.zeros(snapData.NumberOfSnapshots)
		err2 = np.zeros(snapData.NumberOfSnapshots)
		for n in np.arange(snapData.NumberOfSnapshots):
			e = makeplot(snap0, n,param,snapData, pGrid,args.noplot)
			err1[n] 	= e[0]
			err2[n] 	= e[1]
	
		print "\nResulting Errors for run {:d}:".format(param.Run_ID)
		print err1
		print ""
		print err2
		print ""
		fout_name = '{:}{:}Error_Relative.dat'.format(BaseDir,param.OutputDir,param.Run_ID,param.Run_ID)
		fout = open(fout_name,'w')
	
		for n in np.arange(snapData.NumberOfSnapshots):
			fout.write('{:f}\t{:f}\t{:f}\n'.format(snapData.time[n],snapData.rho[n],err1[n]))

		fout.close()

		fout_name = '{:}{:}Error_Relative_FullRange.dat'.format(BaseDir,param.OutputDir,param.Run_ID,param.Run_ID)
		fout = open(fout_name,'w')
	
		for n in np.arange(snapData.NumberOfSnapshots):
			fout.write('{:f}\t{:f}\t{:f}\n'.format(snapData.time[n],snapData.rho[n],err2[n]))

		fout.close()
		snapData = None
		snap0 = None
		pGrid = None
		err = None
		




