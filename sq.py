#!/usr/bin/env python3

# written by Grey Christoforo <first name [at] last name [not] net>
# Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances)
# from Generalized Detailed Balance Theory of Solar Cells By Thomas Kirchartz
# Chapter 2.2, the Shockley-Queisser limit

import numpy as np
from scipy import special
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
etr = np.array([])  # [W/(m^2*nm)]  extraterrestrial
am15 = np.array([])  # [W/(m^2*nm)] AM1.5
dc = np.array([])  # [W/(m^2*nm)] DC

with open(args.solar_spectra_file) as sunfile:
  if '\t' in sunfile.readline():
    dlm = '\t'
  else:
    dlm = ','
  sunfile.seek(0)
  stream = csv.reader(sunfile, delimiter=dlm)
  for row in stream:
    if len(row) == 2:
      lambd = np.append(lambd, np.double(row[0]))
      am15 = np.append(am15, np.double(row[1]))
    else:
      lambd = np.append(lambd,np.double(row[0]))  # [nm] wavelength
      etr = np.append(etr,np.double(row[1]))  # [W/m^2/nm] extra terrestrial radiation
      am15 = np.append(am15,np.double(row[2]))  # [W/m^2/nm] "Global Tilt" = spectral radiation from solar disk plus sky diffuse and diffuse reflected from ground on south facing surface tilted 37 deg from horizontal
      dc = np.append(dc,np.double(row[3]))  # [W/m^2/nm] Direct + Circumsolar
      # where Direct = Direct Normal Irradiance Nearly parallel (0.5 deg divergent cone) radiation on surface with surface normal tracking (pointing to) the sun, excluding scattered sky and reflected ground radiation
      # and Circumsolar = Spectral irradiance within +/- 2.5 degree (5 degree diameter) field of view centered on the 0.5 deg diameter solar disk, but excluding the radiation from the disk

T_cellc = args.t_cell  # [degC]
K_offset = np.double(273.15)
T_cell = K_offset + T_cellc  # [K]

k = 1.3806488e-23  # [J/K] boltzman constant
q = 1.60217657e-19  # [C] or [A*s] elementary charge
h = 6.62607004081e-34  # [m^2*kg/s]planck constant
c = 299792458  # [m/s]speed of light
hc = h*c  # [J*m]
hc_evnm = hc/q*1e9  # [eV*nm]
E_solar = hc/(lambd*1e-9)  # [J] x axis as joules
eV_solar = hc_evnm/lambd  # [eV] x axis as eV
photonDensity = am15/E_solar  # [photons/s/m^2/nm]
Pin = np.trapezoid(am15, lambd)  # input power [W/m^2]
E_max = max(E_solar)  # highest photon energy in the input data [J]

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
  E_BG = 0;  # [J] band gap energy
  cell_type = "semi-perfect"
  nm_abs = np.array([])  # [nm] wavelength
  a_abs = np.array([])  # absorption/emission
  with open(daf) as absfile:
    stream = csv.reader(absfile)
    for row in stream:
      nm_abs = np.append(nm_abs, np.double(row[0]))  # [nm] wavelength
      a_abs = np.append(a_abs, np.double(row[1]))  # [W/m^2/nm] extra terrestrial radiation
  E_abs = 1e9*hc/nm_abs  # [J] x axis as joules
  a = functools.partial(a_file, E_vctr=np.flip(E_abs), a_vctr=np.flip(a_abs))
else:
  E_BG = args.band_gap*q;  # [J] band gap energy
  cell_type = "perfect"
  a = a_gap
  
# takes energy cutoff in joules and returns current in amps per square meter
def current(E_cutoff):
  absorptions = a(E_solar)
  fluxes = photonDensity*absorptions
  bool_vec = E_cutoff<E_solar
  these_fluxes = fluxes[bool_vec]
  these_lambds = lambd[bool_vec]
  if not all(bool_vec):  # we're some place in the middle
    these_fluxes = np.append(these_fluxes, np.interp(E_cutoff, E_solar, fluxes, left=0, right=0))
    these_lambds = np.append(these_lambds, 1e9*hc/E_cutoff)
  photonFlux = np.trapezoid(these_fluxes, these_lambds)
  return photonFlux * q  # [A/m^2]

# takes energy in joules and returns emitted photon flux of a black body of temperature T
def psi_bb(E, T):
  return 2*np.pi*E**2/(h**3*c**2)*np.exp(-E/(k*T))

# inputs:
# the solar cell's absorption (a function which takes energy in joules)
# temperature of the solar cell [K]
# outputs:
# the cell's radiative saturation current [A]
def radiativeSaturationCurrent(a, T):
  absorptions = a(E_solar)
  psis = psi_bb(E_solar, T)
  rscs = absorptions*psis
  r_sat_part = np.trapezoid(np.flip(rscs), np.flip(E_solar))
  return r_sat_part * q  # [A/m^2]

# inputs:
# the solar cell's absorption (a function which takes energy in joules)
# voltage applied to the solar cell [V]
# temperature of the solar cell [K]
# output:
# current through the solar cell when it's dark [A]
def darkCurrent(T, I0, V):
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
def V_mpp (I0, T, Iph):
  return (k*T/q)*(special.lambertw(((I0+Iph)*np.e)/I0)-1)

print(f"We've assumed that our {cell_type} solar cell is at", T_cell, "degrees kelvin")
if daf.is_file():
  print(f"and has an absorption/emission spectra given in the file {daf.name}")
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
print(f"its radiative saturation current{density_str} is")
print(J0/10*area, f"mA{per_str}")

print("")
print(f"and if we shine light as defined by the file {Path(args.solar_spectra_file).name} at it, then")
print("")

J_ph = current(E_BG)
print(f"its photocurrent{density_str} is")
print(J_ph/10*area, f"mA{per_str},")

print("")
print("which makes")
print("")

Voc = openCircuitVoltage(J0, T_cell, J_ph)
print("its open circuit voltage")
print(Voc*1000, "mV")

print("")

J_dark = functools.partial(darkCurrent, T_cell, J0)
Jsc = J_dark(0) - J_ph
print(f"its short circuit current{density_str}")
print(Jsc/10*-1*area, f"mA{per_str}")

print("")

Vmpp = np.real_if_close(V_mpp(J0, T_cell, J_ph))
print("the voltage at its maximum power point")
print(Vmpp*1000, "mV")

print("")

Jmpp = J_dark(Vmpp)-J_ph
print(f"the current{density_str} at its maximum power point")
print(Jmpp/10*-1*area, f"mA{per_str}")

print("")

FF = Jmpp*Vmpp/(Jsc*Voc)
print("its fill factor")
print(FF*100, "percent")

print("")
print("and")
print("")

Pmax = Jmpp*Vmpp*-1
print("its power production is")
print(f"{Pmax/10*area} mW{per_str} ({J_ph*Voc/10*area} if FF was 1.0)")
print("or")
print(f"{Pmax/Pin*100} percent ({J_ph*Voc/Pin*100} if FF was 1.0).")

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
