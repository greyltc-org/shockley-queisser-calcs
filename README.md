# shockley-queisser-calcs
Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances, perfect absorption above the band gap)

## Usage
```bash
$ git clone https://github.com/AFMD/shockley-queisser-calcs.git
$ cd shockley-queisser-calcs
$ sq.py -h
usage: sq.py [-h] [--t_cell T_CELL] [--band_gap BAND_GAP]

Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic
resistances, perfect absorption above the band gap)

optional arguments:
  -h, --help           show this help message and exit
  --t_cell T_CELL      Temperature of the solar cell [deg C]
  --band_gap BAND_GAP  Band gap of the solar cell [eV]
$ ./sq.py --t_cell 50 --band_gap 1.5
We've assumed our perfect solar cell is at 323.15 degrees kelvin and has a band gap
of 1.5 electron volts.

That means

its radiative saturation current density
is 4.159217767932097e-18 mA/cm^2

and if we shine AM1.5 illumination (as defined by ASTM G173) at it,

its photocurrent density
is 30.2901034907 mA/cm^2,

which makes:

its open circuit voltage
1.2094476896 volts.

the voltage at its maximum power point
1.1062234690451316 volts.

the current density at its maximum power point
29.5463350376 mA/cm^2.

its fill factor
89.2192647139 percent.

and

its power conversion efficency
32.6848492429 percent.
```
![Some Graph](/figure_1.png)
