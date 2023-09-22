"""
This script plots the cooling timescales for electrons
against momentum, Lorentz factor and kinetic energy.
One can change physical quantities to see the effect on the cooling times.

It uses the constants and functions defined in Physics.py

September 2023, Léna Jlassi
"""
import sys
import os
cre_analysis_path = '../'
sys.path.append(cre_analysis_path)
import numpy as np
import matplotlib.pyplot as plt
from crest import *
from Physics import *
from logutils import *

output_path = './figures/'

if not os.path.isdir(output_path):
    os.mkdir(output_path)

# Define physical quantities that affect the cooling
n_e = 1.e-4 # cm-3
B_cmb = CMB_MAGNETIC_FIELD # Gauss
B = 1.e-6 # Gauss
B_microG = B * 1.e6 # microGauss


# Define momentum bins and convert to useful quantities
p = np.logspace(-4, 8, num=121)

gamma = gamma_factor(p)

E_erg = p * ELECTRONMASS * CLIGHT**2 # erg 

E_eV = E_erg / ELECTRONVOLT_IN_ERGS  # eV

v = CLIGHT * beta_factor(p)

# classical electron radius
r_0 = 1 * ELECTRONCHARGE**2 / (ELECTRONMASS * CLIGHT**2)

# I used the following relations which use the Lorentz factor
# but there are predefined functions in Physics.py which use the momentum

def b_coul(gamma, n_e):
    # Sarazin 1999 eq 9 (from Rephaeli 1979 eq 2)
    return 1.2e-12 * n_e * (1. + (np.log(gamma/n_e)/75.))

def b_icsync(gamma, B_sync):
    # Sarazin 1999 eq 7
    # converted from energy density field to equivalent magnetic field
    return (THOMPSON * gamma**2 / (6. * np.pi * ELECTRONMASS * CLIGHT)) * (B_sync**2 + B_cmb**2)

# calculate total loss rate
b_total_gamma = b_coul(gamma, n_e) + b_icsync(gamma, B)

# convert total loss rate from lorentz factor to energy and momentum
b_total_energy = b_total_gamma * ELECTRONMASS * CLIGHT**2
b_total_mom = b_total_energy / (beta_factor(p) * ELECTRONMASS * CLIGHT**2)

# can easily change to yr, Myr, Gyr
seconds_in_yr = 60.*60.*24.*365.25
seconds_in_Myr = 60.*60.*24.*365.25*1.e6
seconds_in_Gyr = 60.*60.*24.*365.25*1.e9

# calculate cooling timescales
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
