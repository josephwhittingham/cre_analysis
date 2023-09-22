import sys
cre_analysis_path = '/home/lena/code/crest/cre_analysis/'
sys.path.append(cre_analysis_path)
import numpy as np
from matplotlib.colors import LogNorm
import astropy
import matplotlib.pyplot as plt
from crest import *
from Physics import *
from logutils import *

output_path = '/home/lena/analysis/figures/figures_crest/theory/'

p = np.logspace(-4, 8, num=121)
gamma = gamma_factor(p)

E_erg = p * ELECTRONMASS * CLIGHT**2 # erg 
# E_erg = gamma * ELECTRONMASS * CLIGHT**2
E_eV = E_erg / ELECTRONVOLT_IN_ERGS  # eV

v = CLIGHT * beta_factor(p)
# v = CLIGHT * np.sqrt(1 - (1./gamma**2)) # cm/s

seconds_in_yr = 60.*60.*24.*365.25
seconds_in_Myr = 60.*60.*24.*365.25*1.e6
seconds_in_Gyr = 60.*60.*24.*365.25*1.e9


# kT = 1.e4 # keV
# n_e = 1.e-3 # cm-3
n_e = 15.e-2 # cm-3
# n_gas = n_e / (HYDROGEN_MASS_FRAC * N_ELEC)
B_cmb = CMB_MAGNETIC_FIELD # Gauss
# B = 1.e-6 # Gauss
B = 1 # Gauss
B_microG = B * 1.e6 # microGauss

# classical electron radius
r_0 = 1 * ELECTRONCHARGE**2 / (ELECTRONMASS * CLIGHT**2)


def b_coul(gamma, n_e):
    return 1.2e-12 * n_e * (1. + (np.log(gamma/n_e)/75.))

def b_icsync(gamma, B_sync):
    return (THOMPSON * gamma**2 / (6. * np.pi * ELECTRONMASS * CLIGHT)) * (B_sync**2 + B_cmb**2)

b_total_gamma = b_coul(gamma, n_e) + b_icsync(gamma, B)
b_total_energy = b_total_gamma * ELECTRONMASS * CLIGHT**2
b_total_mom = b_total_energy / (beta_factor(p) * ELECTRONMASS * CLIGHT**2)

t_gamma = gamma / b_total_gamma
t_gamma /= seconds_in_Myr

t_E = E_erg / b_total_energy
t_E /= seconds_in_Myr

t_p = p / b_total_mom
t_p /= seconds_in_Myr

E_GeV = E_eV / 1.e9

fig = plt.figure()

plt.title("n = {:.1e} cm-3     B = {:.1e} microG".format(n_e, B_microG))
plt.plot(gamma, t_gamma, color='red')
plt.xlabel(r'$ \gamma $')
plt.ylabel(r'electron loss timescales $\tau = \frac{\gamma}{b(\gamma)}$ [Myr]')
plt.xlim(1.e0, 1.e5)
# plt.ylim(1.e-3, 1.e4)
plt.xscale('log')
plt.yscale('log')
fig.savefig(output_path + f't_cool_gamma.png', dpi=500)


fig = plt.figure()

plt.title("n = {:.1e} cm-3     B = {:.1e} microG".format(n_e, B_microG))
plt.plot(E_GeV, t_E, color='green')
plt.xlabel(r'kinetic electron energy $ T = pmc^2$ [GeV]')
plt.ylabel(r'electron loss timescales $\tau = \frac{E}{b(E)}$ [Myr]')
plt.xlim(1.e-5, 1.e2)
# plt.ylim(1.e-3, 1.e4)
plt.xscale('log')
plt.yscale('log')
fig.savefig(output_path + f't_cool_E.png', dpi=500)


fig = plt.figure()

plt.title("n = {:.1e} cm-3     B = {:.1e} microG".format(n_e, B_microG))
plt.plot(p, t_p, color='black')
plt.xlabel(r'normalised electron momentum $ p = P / mc$')
plt.ylabel(r'electron loss timescales $\tau = \frac{p}{b(p)}$ [Myr]')
plt.xlim(p.min(), p.max())
# plt.ylim(1.e-3, 1.e4)
plt.xscale('log')
plt.yscale('log')
fig.savefig(output_path + f't_cool_p.png', dpi=500)


p = np.logspace(-3, 8, num=121)
p_low = p[np.where(p < 1.e3)]
p_high = p[np.where(p >= 1.e3)]
alpha = 2.2

f_low = np.zeros((t.shape[0], p_low.shape[0]))
f_high = np.zeros((t.shape[0], p_high.shape[0]))


# b_icsync = sync_loss_rate(p_high, B) + ic_loss_rate(p_low, 0, 0)
b_icsync = ic_sync_loss_rate(p_high, 0, B, 0)
b_coul = coulomb_loss_rate(p_low, n_e)

fig = plt.figure()

for i in range(1, t.shape[0]):
    # f_low[i, :] = E_to_inject*np.exp(-t[i] / injection_timescale) / ((1+p_low**(-2))) * (p_low**(-alpha+1) / (alpha - 1)) 
    f_low[i, :] = E_to_inject*np.exp(-t[i] / injection_timescale) / (b_coul * (1+p_low**(-2))) * (p_low**(-alpha+1) / (alpha - 1)) 
    # f_high[i, :] = E_to_inject*np.exp(-t[i] / injection_timescale) / (p_high**2) * (p_high**(-alpha+1) / (alpha - 1)) 
    f_high[i, :] = E_to_inject*np.exp(-t[i] / injection_timescale) / (b_icsync) * (p_high**(-alpha-1) / (alpha - 1)) 
    plt.plot(p_low, p_low**alpha * f_low[i, :])
    plt.plot(p_high, p_high**alpha * f_high[i, :])

# plt.plot(p_low, p_low**alpha * f_low[-1, :])
# plt.plot(p_high, p_high**alpha * f_high[-1, :])
plt.yscale('log')
plt.xscale('log')