#!/usr/bin/env python3

# written by Grey Christoforo <first name [at] last name [not] net>
# Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances)
# from Generalized Detailed Balance Theory of Solar Cells By Thomas Kirchartz
# Chapter 2.2, the Shockley-Queisser limit

from numpy import *
import scipy
import scipy.integrate as integrate
import functools
import csv
import io
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances, perfect absorption above the band gap)')

parser.add_argument("--t_cell", default=30, type=float_, help="Temperature of the solar cell [deg C]")
parser.add_argument("--band_gap", default=1.34, type=float_, help="Band gap of the solar cell [eV]")

args = parser.parse_args()

lambd = array([]) #[nm] wavelength
etr = array([]) #[W/(m^2*nm)]  extraterrestrial
am15 = array([]) #[W/(m^2*nm)] AM1.5
tot = array([]) #[W/(m^2*nm)] AM1.5
with open('ASTMG173.csv') as csvfile:
  stream = csv.reader(csvfile)
  for row in stream:
    lambd = append(lambd,float(row[0])) #[nm] wavelength
    etr = append(etr,float(row[1])) #[W/m^2/nm] extra terrestrial radiation 
    am15 = append(am15,float(row[2])) #[W/m^2/nm] "Global Tilt" = spectral radiation from solar disk plus sky diffuse and diffuse reflected from ground on south facing surface tilted 37 deg from horizontal
    tot = append(tot,float(row[3])) #[W/m^2/nm] Direct + Circumsolar
    # where Direct = Direct Normal Irradiance Nearly parallel (0.5 deg divergent cone) radiation on surface with surface normal tracking (pointing to) the sun, excluding scattered sky and reflected ground radiation
    # and Circumsolar = Spectral irradiance within +/- 2.5 degree (5 degree diameter) field of view centered on the 0.5 deg diameter solar disk, but excluding the radiation from the disk
csvfile.close()

T_cell = args.t_cell #[degC]
#T_cell = 26.85 #[degC]
K_offset = float_(273.15)
T_cell = K_offset + T_cell # [K]

k = 1.3806488e-23 #[J/K] boltzman constant
q = 1.60217657e-19 #[C] or [A*s] elementary charge
h = 6.62607004081e-34 #[m^2*kg/s]planck constant
c = 299792458 #[m/s]speed of light
hc = h*c #[J*m]
hc_evnm = hc/q*1e9 #[eV*nm]
E_solar = hc/(lambd*1e-9) #[J] x axis as joules
eV_solar = hc_evnm/lambd #[eV] x axis as eV
photonDensity = am15/E_solar # [photons/s/m^2/nm]

E_min = 0
E_max = max(E_solar) #highest energy for absorption (=5 micron wavelength electromagnetic radiation) [J]

# takes energy in joules and returns absorption
# must return values for inputs on [E_min E_max] joules
def a(E, E_BG=1.14*q):
  #TODO: insert a real absorption spectrum here
  
  # let's assume absorption of 1 for photons the bandgap
  a_above = 1
  return (E > E_BG) * a_above
  
# takes energy cutoff in joules and returns current in amps per square meter
def current(E_cutoff):
  cutoffi = argmax(E_cutoff>=E_solar)
  absorptions = a(E_solar[0:cutoffi], E_BG=E_cutoff)
  photonFlux = integrate.trapezoid(photonDensity[0:cutoffi]*absorptions, lambd[0:cutoffi])# [photons/s/m^2]
  return photonFlux * q # [A/m^2]

# takes energy in joules and returns emitted photon flux of a black body of temperature T
def psi_bb(E,T):
  return 2*pi*E**2/(h**3*c**2)*exp(-E/(k*T))
  
#def psi_approx(a,V,E):
#  return a(E)*psi_bb(E)*exp(qV/(k*T))

#def psi(a,V,E):
#  return 2*pi*E**2/(h**3*c**2)*a(E)/(exp((E-q*V)/(k*T)-1))

# inputs:
# the solar cell's absorption (a function which takes energy in joules)
# temperature of the solar cell [K]
# outputs:
# the cell's radiative saturation current [A]
def radiativeSaturationCurrent(a, T):
  return q*integrate.quad(lambda x: a(x)*psi_bb(x,T), E_min, E_max)[0]

# inputs:
# the solar cell's absorption (a function which takes energy in joules)
# voltage applied to the solar cell [V]
# temperature of the solar cell [K]
# output:
# current through the solar cell when it's dark [A]
def darkCurrent(T,I0,V):
  return I0*(exp(q*V/(k*T))-1)

# inputs:
# radiativeSaturationCurrent [A]
# temperature of the cell [K]
# photogenerated current [A]
# output:
# open circuit voltage [V]
def openCircuitVoltage (I0, T, Iph):
  return(k*T/q)*log(Iph/I0+1)

# voltage at max power point
def V_mpp (I0,T,Iph):
  return (k*T/q)*(scipy.special.lambertw(((I0+Iph)*e)/I0)-1)

E_BG = args.band_gap*q; #[J] band gap energy
print("We've assumed our perfect solar cell is at", T_cell, "degrees kelvin and has a band gap")
print("of", E_BG/q, "electron volts.")

print("")
print("That means")
print("")

# device absorption
aDevice = functools.partial(a, E_BG=E_BG)
J0 = radiativeSaturationCurrent(aDevice, T_cell)
print("its radiative saturation current density")
print("is", J0/10, "mA/cm^2")

print("")
print("and if we shine AM1.5 illumination (as defined by ASTM G173) at it,")
print("")

J_ph = current(E_BG)
print("its photocurrent density")
print("is", J_ph/10, "mA/cm^2,")

print("")
print("which makes:")
print("")

Voc = openCircuitVoltage(J0, T_cell, J_ph)
print("its open circuit voltage")
print(Voc, "volts.")

#print("")

J_dark = functools.partial(darkCurrent, T_cell, J0)
Jsc = J_dark(0) - J_ph
#print("its short circuit current density")
#print(J_sc/10, "mA/cm^2.")

print("")

Vmpp = real_if_close(V_mpp(J0,T_cell,J_ph))
print("the voltage at its maximum power point")
print(Vmpp, "volts.")

print("")

Jmpp = J_dark(Vmpp)-J_ph
print("the current density at its maximum power point")
print(Jmpp/10*-1, "mA/cm^2.")

print("")

FF = Jmpp*Vmpp/(Jsc*Voc)
print("its fill factor")
print(FF*100, "percent.")

print("")
print("and")
print("")

Pmax = Jmpp*Vmpp*-1
print("its power conversion efficency")
print(f"{Pmax/10} percent ({J_ph*Voc/10} if FF was 1.0).")

nPoints = 1000
vMin = -0.2 #[V]
vMax = Voc*1.2 #[V]
v = linspace(vMin,vMax,nPoints)


plt.plot(v, -1*J_dark(v)/10, v, -1*J_dark(v)/10+J_ph/10)
plt.xlabel('Terminal Voltage [V]')
plt.ylabel('Current Density [mA/cm^2]')
buffer = io.StringIO()
print("The Current Through A Perfect,", E_BG/q, "eV Solar Cell at",T_cell-K_offset, "deg C",file=buffer,end='')
plt.title(buffer.getvalue())
plt.legend(['In the Dark','Under AM1.5'], loc='best')
plt.ylim([-J_ph*1.1/10,J_ph*1.1/10])
#plt.xlim([vMin,vMax])
plt.grid('on')
plt.show()
