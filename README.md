# shockley-queisser-calcs
Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances, perfect absorption above the band gap)

## Usage
```bash
$ git clone https://github.com/AFMD/shockley-queisser-calcs.git
$ cd shockley-queisser-calcs
$ ./sq.py --help
usage: sq.py [-h] [--t-cell T_CELL] [--band-gap BAND_GAP] [--no-plot]

Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances, perfect absorption above the band gap)

optional arguments:
  -h, --help           show this help message and exit
  --t-cell T_CELL      Temperature of the solar cell [deg C]
  --band-gap BAND_GAP  Band gap of the solar cell [eV]
  --no-plot            Disable plot
$ ./sq.py --t-cell 50 --band-gap 1.5 --no-plot
We've assumed our perfect solar cell is at 323.15 degrees kelvin and has a band gap
of 1.5 electron volts.

That means

its radiative saturation current density
is 4.159217767932097e-18 mA/cm^2

and if we shine AM1.5 illumination (as defined by ASTM G173) at it,

its photocurrent density
is 28.937181142314238 mA/cm^2,

which makes:

its open circuit voltage
1.2081752604715137 volts.

its short circuit current density
28.937181142314238 mA/cm^2.

the voltage at its maximum power point
1.1049815516988595 volts.

the current density at its maximum power point
28.22585450065439 mA/cm^2.

its fill factor
89.21049812068259 percent.

and

its power conversion efficency
31.189048504159324 percent (34.96118636392688 if FF was 1.0).
```
![Some Graph](/figure_1.svg)
