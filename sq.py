#!/usr/bin/env python3

# written by Grey Christoforo <first name [at] last name [not] net>
# Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances)
# from Generalized Detailed Balance Theory of Solar Cells By Thomas Kirchartz
# Chapter 2.2, the Shockley-Queisser limit

import numpy as np
import scipy
import scipy.integrate as integrate
import functools
import csv
import io
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--t-cell", default=25, type=np.double, help="Temperature of the solar cell [deg C]")
parser.add_argument("--band-gap", default=1.35, type=np.double, help="Band gap of the solar cell [eV] (ignored if --device-absorption-file is given)")
parser.add_argument("--no-plot", default=False, action='store_true', help="Disable plot")
parser.add_argument("--solar-spectra-file", default="ASTMG173.csv", help="File to read the solar spectra from")
parser.add_argument("--device-absorption-file", default="", help="File to read the absorption spectrum from")
parser.add_argument("--sqcm", type=np.double, help="Set the device's area [cm^2]")

args = parser.parse_args()

lambd = np.array([]) #[nm] wavelength
etr = np.array([]) #[W/(m^2*nm)]  extraterrestrial
am15 = np.array([]) #[W/(m^2*nm)] AM1.5
tot = np.array([]) #[W/(m^2*nm)] AM1.5
with open(args.solar_spectra_file) as sunfile:
  stream = csv.reader(sunfile)
  for row in stream:
    lambd = np.append(lambd,np.double(row[0])) #[nm] wavelength
    etr = np.append(etr,np.double(row[1])) #[W/m^2/nm] extra terrestrial radiation 
    am15 = np.append(am15,np.double(row[2])) #[W/m^2/nm] "Global Tilt" = spectral radiation from solar disk plus sky diffuse and diffuse reflected from ground on south facing surface tilted 37 deg from horizontal
    tot = np.append(tot,np.double(row[3])) #[W/m^2/nm] Direct + Circumsolar
    # where Direct = Direct Normal Irradiance Nearly parallel (0.5 deg divergent cone) radiation on surface with surface normal tracking (pointing to) the sun, excluding scattered sky and reflected ground radiation
    # and Circumsolar = Spectral irradiance within +/- 2.5 degree (5 degree diameter) field of view centered on the 0.5 deg diameter solar disk, but excluding the radiation from the disk

T_cell = args.t_cell #[degC]
#T_cell = 26.85 #[degC]
K_offset = np.double(273.15)
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
def a_gap(E, E_BG=1.14*q):
  # let's assume absorption of 1 for photons above the bandgap
  a_above = 1
  return (E > E_BG) * a_above

def a_file(E, E_BG=0, E_vctr=np.array([]), a_vctr=np.array([])):
  return np.interp(E, E_vctr, a_vctr, left=0, right=0)

daf = Path(args.device_absorption_file)
if daf.is_file():
  cell_type = "semi-perfect"
  nm_abs = np.array([]) # [nm] wavelength
  a_abs = np.array([]) # absorption/emission
  with open(daf) as absfile:
    stream = csv.reader(absfile)
    for row in stream:
      nm_abs = np.append(nm_abs, np.double(row[0]))  # [nm] wavelength
      a_abs = np.append(a_abs, np.double(row[1]))  # [W/m^2/nm] extra terrestrial radiation
  E_abs = hc/(nm_abs*1e-9)  # [J] x axis as joules
  a = functools.partial(a_file, E_vctr=np.flip(E_abs), a_vctr=np.flip(a_abs))
else:
  cell_type = "perfect"
  a = a_gap
  
# takes energy cutoff in joules and returns current in amps per square meter
def current(E_cutoff):
  cutoffi = np.argmax(E_cutoff>=E_solar)
  absorptions = a(E_solar[0:cutoffi], E_BG=E_cutoff)
  photonFlux = integrate.trapezoid(photonDensity[0:cutoffi]*absorptions, lambd[0:cutoffi])  # [photons/s/m^2]
  return photonFlux * q # [A/m^2]

# takes energy in joules and returns emitted photon flux of a black body of temperature T
def psi_bb(E,T):
  return 2*np.pi*E**2/(h**3*c**2)*np.exp(-E/(k*T))
  
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
  return I0*(np.exp(q*V/(k*T))-1)

# inputs:
# radiativeSaturationCurrent [A]
# temperature of the cell [K]
# photogenerated current [A]
# output:
# open circuit voltage [V]
def openCircuitVoltage (I0, T, Iph):
  return(k*T/q)*np.log(Iph/I0+1)

# voltage at max power point
def V_mpp (I0,T,Iph):
  return (k*T/q)*(scipy.special.lambertw(((I0+Iph)*np.e)/I0)-1)

E_BG = args.band_gap*q; #[J] band gap energy
print(f"We've assumed that our {cell_type} solar cell is at", T_cell, "degrees kelvin")
if daf.is_file():
  print(f"and has an absorption/emission spectra given in {daf}")
else:
  print("and has a band gap of", E_BG/q, "electron volts.")

if args.sqcm:
  print(f"and has an area of {args.sqcm} cm^2.")
  density_str = ""
  per_str = ""
  area_str = f"{args.sqcm} cm^2 "
  area = args.sqcm
else:
  density_str = " density"
  per_str = "/cm^2"
  area_str = f""
  area = 1

print("")
print("That means")
print("")

# device absorption
aDevice = functools.partial(a, E_BG=E_BG)
J0 = radiativeSaturationCurrent(aDevice, T_cell)
print(f"its radiative saturation current{density_str}")
print("is", J0/10*area, f"mA{per_str}")

print("")
print("and if we shine AM1.5 illumination (as defined by ASTM G173) at it,")
print("")

J_ph = current(E_BG)
print(f"its photocurrent{density_str}")
print("is", J_ph/10*area, f"mA{per_str},")

print("")
print("which makes:")
print("")

Voc = openCircuitVoltage(J0, T_cell, J_ph)
print("its open circuit voltage")
print(Voc, "volts.")

print("")

J_dark = functools.partial(darkCurrent, T_cell, J0)
Jsc = J_dark(0) - J_ph
print(f"its short circuit current{density_str}")
print(Jsc/10*-1*area, f"mA{per_str}.")

print("")

Vmpp = np.real_if_close(V_mpp(J0, T_cell, J_ph))
print("the voltage at its maximum power point")
print(Vmpp, "volts.")

print("")

Jmpp = J_dark(Vmpp)-J_ph
print(f"the current{density_str} at its maximum power point")
print(Jmpp/10*-1*area, f"mA{per_str}.")

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

if args.no_plot == False:
  nPoints = 1000
  vMin = -0.2 #[V]
  vMax = Voc*1.2 #[V]
  v = np.linspace(vMin,vMax,nPoints)

  plt.plot(v, -1*J_dark(v)/10*area, v, -1*J_dark(v)/10+J_ph/10*area)
  plt.xlabel('Terminal Voltage [V]')
  plt.ylabel(f'Current{density_str} [mA{per_str}]')
  buffer = io.StringIO()
  if daf.is_file():
    print(f"The Current Through A {area_str}Semi-Perfect Solar Cell at",T_cell-K_offset, "deg C",file=buffer,end='')
  else:
    print(f"The Current Through A {area_str}Perfect,", E_BG/q, "eV Solar Cell at",T_cell-K_offset, "deg C",file=buffer,end='')
  plt.title(buffer.getvalue())
  plt.legend(['In the Dark','Under AM1.5'], loc='best')
  plt.ylim([-J_ph*1.1/10*area,J_ph*1.1/10*area])
  #plt.xlim([vMin,vMax])
  plt.grid('on')
  plt.show()
