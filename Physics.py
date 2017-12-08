import numpy as np
from scipy.integrate import simps, quad

THOMPSON                = 6.65245873e-25
ELECTRONMASS            = 9.1093829e-28
CLIGHT                    = 2.99792458e10
PLANCK                    = 6.6260695e-27
ELECTRONCHARGE            = 4.8032042e-10
PROTONMASS                = 1.67262178e-24
YEARS_IN_SECONDS         = 3.1536e7
ELECTRONVOLT_IN_ERGS    = 1.60217656e-12
BARN_IN_CM2                = 1.E-24
BOLTZMANN               = 1.38065e-16
CMB_MAGNETIC_FIELD        = 3.24e-6         # (Gauss)
CMB_ENERGY_DENSITY        = 4.165659e-13        # (ergs / cm^3)
ALPHA_FINESTRUCTURE       = 0.0072973526    # CODATA 2014

# physical functions
def plasma_frequency(n_e):
	# calculate the plasma frequency $ \omega = \sqrt{ 4 \pi n_e e^2 / m_e }$
	return np.sqrt(np.pi * n_e / ELECTRONMASS ) * 2 * ELECTRONCHARGE


def coulomb_loss_rate(p, n_e ):
	# Exact form for Coulomb loss rate after: Gould, R.J., 1972, Physica, 50, 145
	# Exact form can be reduced to that given in Schlickeiser if p >= 10
	# remind that PLANCK = h, but $\hbar = h / (2 \pi)$
	return (3. * THOMPSON * n_e * CLIGHT ) / (2. * (beta_factor(p))**2 ) * (
				+ np.log(ELECTRONMASS * CLIGHT**2 * beta_factor(p) * np.sqrt(gamma_factor(p) - 1. ) * 2. * np.pi  / ( PLANCK * plasma_frequency(n_e) ))
				- np.log(2.) * ((beta_factor(p))**2/2. + 1./gamma_factor(p))  + 0.5 + ( (gamma_factor(p)-1.) / (4 * gamma_factor(p)) )**2
				)

def ic_sync_loss_rate(p, u_photon, magnetic_field, redshift_z):
	cmb_energy_density     = 0.26# [eV cm^-3]
	magnetic_field_CMB     = 3.24E-6    # [Gauss]
	return 4./3. * THOMPSON / (ELECTRONMASS * CLIGHT) * p**2 / beta_factor(p) * ( (np.power(1+redshift_z,4) + (magnetic_field / CMB_MAGNETIC_FIELD)**2)  * CMB_ENERGY_DENSITY + u_photon) 

def ic_loss_rate(p, u_photon, redshift_z):
	cmb_energy_density     = 0.26# [eV cm^-3]
	magnetic_field_CMB     = 3.24E-6    # [Gauss]
	return 4./3. * THOMPSON / (ELECTRONMASS * CLIGHT) * p**2 / beta_factor(p) * ( np.power(1+redshift_z,4) * CMB_ENERGY_DENSITY + u_photon) 

def sync_loss_rate(p, magnetic_field):
	cmb_energy_density     = 0.26# [eV cm^-3]
	magnetic_field_CMB     = 3.24E-6    # [Gauss]
	return 4./3. * THOMPSON / (ELECTRONMASS * CLIGHT) * p**2 / beta_factor(p) * (magnetic_field / CMB_MAGNETIC_FIELD )**2  * CMB_ENERGY_DENSITY


def gamma_factor(p):
	# calculates the relativistic gamma factor
	if type(p)==np.ndarray or type(p)==list:
		return np.array([np.sqrt(1. + pp**2) for pp in p ])
	else:
		return np.sqrt(1. + p**2)

def beta_factor(p):
	if type(p)==np.ndarray or type(p)==list:
		return np.array([(pp / np.sqrt( 1. + pp**2 )) for pp in p ])
	else:
		return (p / np.sqrt( 1. + p**2 ))

	return 

def electron_density(rho, x_e, HydrogenMassFrac):
	# calculates the electron density $n_e = \rho / m_p  X_H x_e$
	# X_H Hydrogen mass fraction, x_e free electron mass fraction
	return rho * HydrogenMassFrac * x_e;
	# rho is given in cm^(-3), therefore no devision by PROTONMASS necessary ;

def coulomb_loss_rate_proton(p, n_e):
	return 4. * np.pi * ELECTRONCHARGE**4 * n_e / ( ELECTRONMASS * PROTONMASS * beta_factor(p)**2 * CLIGHT**3 ) * (
				    np.log(4. * np.pi * ELECTRONMASS * CLIGHT**2 * beta_factor(p) * p  / (PLANCK * plasma_frequency(n_e))  ) -0* beta_factor(p)**2/2.
				) 
#

def ionisation_loss_rate_proton(p, rho, HydrogenMassFrac):
	return  4. * np.pi * ELECTRONCHARGE**4 / ( ELECTRONMASS * PROTONMASS * beta_factor(p)**2 * CLIGHT**3 ) * rho * (
				   HydrogenMassFrac                * np.log(2. * ELECTRONMASS * CLIGHT**2 * p**2 / ( 13.6 * ELECTRONVOLT_IN_ERGS)  )
				+ (1. - HydrogenMassFrac)/4.     * np.log(2. * ELECTRONMASS * CLIGHT**2 * p**2 / ( 24.6 * ELECTRONVOLT_IN_ERGS)  )
				- (1. + HydrogenMassFrac)/2. *  beta_factor(p)**2
				)

def hadronic_loss_rate(p, rho, AlphaSpectralIndex, HydrogenMassFrac):
	if p >= 0.78E9* ELECTRONVOLT_IN_ERGS / (PROTONMASS * CLIGHT**2):
		return (CLIGHT * 2. * rho / (1+ HydrogenMassFrac) * 32. * (0.96 + np.exp(4.4 - 2.4 * AlphaSpectralIndex) ) * 1.E-3 * BARN_IN_CM2  * 0.5 * (np.sqrt(1 + p**2)-1)) / beta_factor(p)
	else:
		return 0.0
 
def CoulombCoefficient(n_e):
	return 4.0 * np.pi * ELECTRONCHARGE**4 * n_e / (ELECTRONMASS * PROTONMASS * CLIGHT**3) * np.log(4. * np.pi *  ELECTRONMASS * CLIGHT**2 / (PLANCK * plasma_frequency(n_e)))

def HadronicCoefficient(rho,alpha, HydrogenMassFrac):
	return 32. * (0.96 + np.exp(4.4 - 2.4 * alpha)) * 1.E-3 * BARN_IN_CM2 * 0.5 * CLIGHT * 2 *  rho / (1+ HydrogenMassFrac)


# Bremsstrahlung and associated functions
def F_int(x):
	return np.log(1 + x) / x

def F(x):
	return quad(F_int,0,x, limit=100,maxp1=100,limlst=100)[0]

def F_err(x):
	return quad(F_int,0,x)[1]


def chi(p):
	#if type(p) is np.ndarray:
	#    return np.array([3./4. * (np.log(2 * (np.square(pp) + 1.)) - 1./3.) for pp in p])
	#else:
	#    return 3./4. * (np.log(2 * (np.square(p) + 1.)) - 1./3.)
	
	if type(p) is np.ndarray:
		return np.array([( + (12. * np.square(pp) + 16. )/(3. * np.sqrt(pp**2 + 1.)  * pp ) * np.log( np.sqrt(pp**2 + 1.)  + pp )
				      - (8 * np.sqrt(pp**2 + 1.) + 6 * pp)/( 3. * np.sqrt(pp**2 + 1.) * pp**2) * np.square(np.log(np.sqrt(pp**2 + 1.) + pp ))
				      - 4./3.
				           + 2. /  (np.sqrt(pp**2 + 1.) * pp) * F( 2. * pp * (np.sqrt(pp**2 + 1.) + pp))) for pp in p])
	else:
		return ( + (12. * np.square(p) + 16. )/(3. * np.sqrt(p**2 + 1.)  * p ) * np.log( np.sqrt(p**2 + 1.)  + p )
				 - (8 * np.sqrt(p**2 + 1.) + 6 * p)/( 3. * np.sqrt(p**2 + 1.) * p**2) * np.square(np.log(np.sqrt(p**2 + 1.) + p ))
				 - 4./3.
				 + 2. /  (np.sqrt(p**2 + 1.) * p) * F( 2. * p * (np.sqrt(p**2 + 1.) + p)))

	
def bremsstrahlung_loss_rate(p, n_gas):
	return n_gas * gamma_factor(p) * ELECTRONCHARGE**4 / (ELECTRONMASS**2 * CLIGHT**3) * ALPHA_FINESTRUCTURE * chi(p)



######################################################################################
def SolutionSteadyStateElectrons(p, param, ne, u_photon, B, z):
	# p         Momenta
	# param     instance of CRelectronParameters class
	# ne        electron density        
	# B        magnetic field (G)
	# z        redshift
	fA = np.zeros(p.size)

	for i in np.arange(p.size):
		if p[i] >= param.SourceLowCutoff:
			fA[i] = param.SourceNormalization * pow(p[i],1-param.SourceSpectralIndex) / abs(coulomb_loss_rate(p[i],ne) + ic_sync_loss_rate(p[i], u_photon, B, z)) / (param.SourceSpectralIndex - 1)
		else:
			fA[i] = param.SourceNormalization * pow(param.SourceLowCutoff,1-param.SourceSpectralIndex) / abs(coulomb_loss_rate(p[i],ne) + ic_sync_loss_rate(p[i], u_photon, B, z)) / (param.SourceSpectralIndex - 1)
	
	return fA

######################################################################################
def SolutionSteadyStateElectronsWithCutoff(p, param, cInj, comp, V, B, BRat, ne, n_gas, u_photon, z):
	# p         Momenta
	# param     instance of CRelectronParameters class
	# c_inj     injection rate
	# comp      shock compression ratio
	# B         magnetic field
	# ne        electron number density
	# n_gas     gas number density
	# u_star_cmb_ratio    ratio of star light to CMB photon energies
	# z         redshift
	# V         upstream gas velocity

	fA = np.zeros(p.size)

	# first calculate the parameters of the steady state spectrum
	pCut = 0.5 * V / CLIGHT *np.sqrt( (comp - 1.) * 3. * ELECTRONCHARGE * B / (THOMPSON * comp * 
				  (CMB_ENERGY_DENSITY * ( ( 1. + comp * BRat)*np.power(1+z,4) + np.square(B / CMB_MAGNETIC_FIELD) * (1. + comp / BRat) )
				   + u_photon * ( 1 + + comp * BRat) ))) 
	alphaInj = (comp + 2.) / ( comp - 1.) # correct 1D slope

	fA = cInj / ((alphaInj - 1.) * (coulomb_loss_rate(p, ne) + bremsstrahlung_loss_rate(p, n_gas) + ic_sync_loss_rate(p, u_photon, B, z))) * np.power(p, - alphaInj + 1) * np.power((1. + param.ShockParamA * np.power(p/pCut, param.ShockParamB)), param.ShockParamC) * np.exp( - np.square( p / pCut))
	
	return fA


######################################################################################
def SolutionSteadyStateProtons(p, param, ne, alpha, rho):
	# p         Momenta
	# param     instance of CRelectronParameters class
	# ne        electron density
	# alpha        Spectral Index of CRp population
	# rho        (number) density of plasma

	fA = np.zeros(p.size)

	for i in np.arange(p.size):
		if p[i] >= param.SourceLowCutoff:
			fA[i] = param.SourceNormalization * pow(p[i],(1-param.SourceSpectralIndex)) / (abs(coulomb_loss_rate_proton(p[i],ne) + ionisation_loss_rate_proton(p[i], rho, param.HydrogenMassFrac) + hadronic_loss_rate(p[i], rho, param.SourceSpectralIndex,param.HydrogenMassFrac)) * (param.SourceSpectralIndex - 1))
		else:
			fA[i] = param.SourceNormalization * pow(param.SourceLowCutoff,1-param.SourceSpectralIndex) / abs(coulomb_loss_rate_proton(p[i],ne) + ionisation_loss_rate_proton(p[i], rho, param.HydrogenMassFrac) + hadronic_loss_rate(p[i], rho, param.SourceSpectralIndex,param.HydrogenMassFrac) ) / (param.SourceSpectralIndex - 1)
	
	return fA

######################################################################################
def SolutionSteadyStateProtonsAsymptotic(p, param, A_C, pStar):
	# p         Momenta
	# param     instance of CRelectronParameters class
	# A_C        Coulomb loss rate coefficient
	# pStar    cross-over-momentum

	return np.array([  param.SourceNormalization * np.power(pp,-param.SourceSpectralIndex) / ( (param.SourceSpectralIndex - 1.) * A_C * (np.power(pStar,-3) + np.power(pp,-3)))       for pp in p])


############
# calculate the spectral energy distribution function based on the cr density distribution function
def SpectralEnergyDensity(p,f):
	E_Dens = np.zeros(p.size)
	for i in np.arange(p.size):
		E_Dens[i] = f[i] * p[i] * (np.sqrt(1 + p[i]**2)-1)

	return E_Dens

############
# calculate the spectral energy distribution function based on the cr density distribution function
def EnergyDensity(p,f):
	E_Dens = np.zeros(p.size)
	for i in np.arange(p.size):
		E_Dens[i] = f[i] * (np.sqrt(1 + p[i]**2)-1)

	return E_Dens

##########
# Analytic Solution for Adiabatic Only Changes
def SolutionAdiabaticChanges(p, rho_0, rho_curr, param):
	# rho_0, C_0, q_0 initial values for density, normalization, cutoff
	fA = np.zeros(p.size)
	# higher and lower cutoffs
	if param.MomentumLowCutoff < p[2]:
		pL = 0.
	else:
		pL = ( param.MomentumLowCutoff * pow(rho_curr/rho_0,1./3.))

	if param.MomentumHighCutoff > p[-2]:
		pH = 1e99
	else:
		pH = ( param.MomentumHighCutoff * pow(rho_curr/rho_0,1./3.))
	for i in np.arange(p.size):
		if p[i] >= pL and p[i] <= pH:
			fA[i] = pow((rho_curr/rho_0),(param.AlphaSpectralIndex+2.)/3.) * param.NormalizationFactor * np.power(p[i],-param.AlphaSpectralIndex)
		else:
			fA[i] = 1.E-100
	
	return fA

##########
# Analytic Solution for Adiabatic Only Changes
def SolutionApproxProtons(p, C, q, alpha):
	fA = np.zeros(p.size)
	for i in np.arange(p.size):
		if p[i] >= q:
			fA[i] = C * np.power(p[i],-alpha)
		else:
			fA[i] = 1.E-100
	
	return fA

##########
# Adiabatic Cutoff
def CutoffsAdiabaticChanges(p, rho_0, rho_curr, param):
	pL = (max(param.MomentumLowCutoff, p[2]        )     * pow(rho_curr/rho_0,1./3.))
	pH = (min(param.MomentumHighCutoff, p[-2]    )     * pow(rho_curr/rho_0,1./3.))
	return pL, pH
	
#############
# Kinetic Energy of Protons
def kinetic_energy_proton(p):
	if type(p)==np.ndarray or type(p)==list:
		return np.array([(np.sqrt(1 + pp**2) - 1)*PROTONMASS*CLIGHT**2  for pp in p])
	else:
		return (np.sqrt(1 + p**2) - 1)*PROTONMASS*CLIGHT**2

def kinetic_energy_electron(p):
	if type(p)==np.ndarray or type(p)==list:
		return np.array([(np.sqrt(1 + pp**2) - 1)*ELECTRONMASS*CLIGHT**2  for pp in p])
	else:
		return (np.sqrt(1 + p**2) - 1)*ELECTRONMASS*CLIGHT**2


#######################################################################################
# Class for Calculating characteristic values of CR distribution
# energy, pressure and number
# Analytic, Simulatex, Approximated
# Calcs the error relative to the analytic solution

class CharacteristicValues:
	def __init__(self,FlagProtons,imin,imax,p,fAna,fSim,fApp=None):
		self.FlagProtons = FlagProtons
		self.imin        = imin
		self.imax         = imax

		# Calculate the characteristic properties for electrons        
		if self.FlagProtons == 0:        
			self.energy_ana   = self.calc_energy_electron(p,fAna)
			self.pressure_ana = self.calc_pressure_electron(p,fAna)
			
			self.energy_sim   = self.calc_energy_electron(p,fSim)
			self.pressure_sim = self.calc_pressure_electron(p,fSim)
			
			if fApp is not None:
				self.energy_app   = self.calc_energy_electron(p,fApp)
				self.pressure_app = self.calc_pressure_electron(p,fApp)


		# Calculate the characteristic properties for protons        
		else:        
			self.energy_ana   = self.calc_energy_proton(p,fAna)
			self.pressure_ana = self.calc_pressure_proton(p,fAna)
			
			self.energy_sim   = self.calc_energy_proton(p,fSim)
			self.pressure_sim = self.calc_pressure_proton(p,fSim)
			
			if fApp is not None:
				self.energy_app   = self.calc_energy_proton(p,fApp)
				self.pressure_app = self.calc_pressure_proton(p,fApp)

		
		# Calculate the number and the relative errors
		self.number_ana = self.calc_number(p,fAna)
		self.number_sim    = self.calc_number(p,fSim)

		if self.energy_ana == 0.0:
			self.del_energy_sim     = None
		else:
			self.del_energy_sim     = (self.energy_sim   - self.energy_ana)   / self.energy_ana

		if self.pressure_ana == 0.0:
			self.del_pressure_sim   = None
		else:
			self.del_pressure_sim   = (self.pressure_sim - self.pressure_ana) / self.pressure_ana
		
		if self.number_ana == 0.0:
			self.del_number_sim     = None
		else:
			self.del_number_sim     = (self.number_sim   - self.number_ana)   / self.number_ana

		if fApp is None:
			self.energy_app         = None
			self.pressure_app       = None
			self.number_app         = None

			self.del_energy_app     = None
			self.del_pressure_app   = None
			self.del_number_app     = None
		else:
			self.number_app = self.calc_number(p,fApp)
			
			self.del_energy_app     = (self.energy_app   - self.energy_ana)   / self.energy_ana
			self.del_pressure_app   = (self.pressure_app - self.pressure_ana) / elf.pressure_ana        
			self.del_number_app     = (self.number_app   - self.number_ana)   / self.number_ana


	def calc_energy_proton(self, p, f):
		return simps(np.multiply(f[self.imin:self.imax], kinetic_energy_proton(p[self.imin:self.imax])) , p[self.imin:self.imax])

	def calc_energy_electron(self, p, f):
		return simps(np.multiply(f[self.imin:self.imax], kinetic_energy_electron(p[self.imin:self.imax])) , p[self.imin:self.imax])
	
	def calc_pressure_proton(self, p, f):
		return PROTONMASS * CLIGHT**2 / 3. * simps(np.multiply(f[self.imin:self.imax],np.multiply(p[self.imin:self.imax], beta_factor(p)[self.imin:self.imax])), p[self.imin:self.imax])

	def calc_pressure_electron(self, p, f):
		return ELECTRONMASS * CLIGHT**2 / 3. * simps(np.multiply(f[self.imin:self.imax],np.multiply(p[self.imin:self.imax], beta_factor(p)[self.imin:self.imax])), p[self.imin:self.imax])

	def calc_number(self, p, f):
		return simps(f[self.imin:self.imax], p[self.imin:self.imax])

def CalculateMomentumGrid(MinimumMomentum, MaximumMomentum, NumberOfMomentumBins, IncludeMaximumMomentum=False):
	if IncludeMaximumMomentum:
		del_p = (np.log(MaximumMomentum) - np.log(MinimumMomentum)) / (NumberOfMomentumBins - 1)
	else:
		del_p = (np.log(MaximumMomentum) - np.log(MinimumMomentum)) / (NumberOfMomentumBins)

	p = np.array([MinimumMomentum * np.exp( del_p * i) for i in np.arange(NumberOfMomentumBins)])

	return p
