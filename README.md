# shockley-queisser-calcs
Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances)

## Usage
```bash
We've assumed our perfect solar cell is at 303.15 degrees kelvin and has a band gap
of 1.34 electron volts.

That means

its radiative saturation current density
is 4.0801766261360677e-17 mA/cm^2

and if we shine AM1.5 illumination (as defined by ASTM G173) at it,

its photocurrent density
is 36.3717514429 mA/cm^2,

which makes:

its open circuit voltage
1.07972366497 volts.

the voltage at its maximum power point
0.9842365031715747 volts.

the current density at its maximum power point
35.4313386912 mA/cm^2.

its fill factor
88.7994334123 percent.

and

its power conversion efficency
34.8728168961 percent.
```
![Some Graph](/figure_1.png)
