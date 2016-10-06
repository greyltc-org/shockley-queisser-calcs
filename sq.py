#!/usr/bin/env python3

# written by Grey Christoforo <first name [at] last name [not] net>
# Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances)
# from Generalized Detailed Balance Theory of Solar Cells By Thomas Kirchartz
# Chapter 2.2, the Shockley-Queisser limit

from numpy import *
import scipy.integrate as integrate
import functools

T_cell = 29 #[degC]
T_cell = 273.15 + T_cell # [K]

k = 1.3806488e-23 #boltzman constant
q = 1.60217657e-19 #electron charge
h = 6.62607004081e-34 #planck constant
c = 299792458 #speed of light


#E_min = 3.9733979e-20 #lower bound for energy calculation (=0.1 micron wavelength electromagnetic radiation) [J]
E_min = 0
E_max = 1.98637857e-18 #highest energy for absorption (=5 micron wavelength electromagnetic radiation) [J]

# takes energy in joules and returns absorption
# must return values for inputs on [E_min E_max] joules
def a(E):
  #TODO: insert a real absorption spectrum here
  SiBG = 1.14*q # let's assume absorption of 1 below Si's bandgap for now
  value = SiBG - E_min
  return value

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
def radiativeSaturationCurrent(a,T):
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
def openCircuitVoltage (I0,T,Iph):
  return(k*T/q)*log(Iph/I0+1)

nPoints = 1000
vMin = -0.2 #[V]
vMax = 1 #[V]
v = linspace(vMin,vMax,nPoints)

I0 = radiativeSaturationCurrent(a,T_cell)

print("The radiative saturation current of our solar cell")
print("is", I0, "amperes.")

print("")
print("and")
print("")

I_ph = 50e-3 #photocurrent [A]
Voc = openCircuitVoltage(I0, T_cell, I_ph)
print("If the photocurrent is", I_ph, "amperes, then")
print("the maximum open circuit voltage is", Voc, "volts.")

i_dark = fromiter(map(functools.partial(darkCurrent, T_cell,I0), v),float)
i_light = i_dark - I_ph

import matplotlib.pyplot as plt
plt.plot(v, i_dark, v, i_light)
plt.xlabel('Volts')
plt.ylabel('Amperes')
plt.title('The Current Through A Solar Cell')
plt.legend(['In the Dark','When Illuminated'], loc='best')
plt.ylim([-I_ph*1.1,I_ph*3])
plt.grid('on')
plt.show()
