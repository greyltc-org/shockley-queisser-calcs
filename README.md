# shockley-queisser-calcs
Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances, optionally perfect absorption above the band gap)

## Usage
```
$ git clone https://github.com/AFMD/shockley-queisser-calcs.git
$ cd shockley-queisser-calcs
$ ./sq.py --help
usage: sq.py [-h] [--t-cell T_CELL] [--band-gap BAND_GAP] [--no-plot] [--solar-spectra-file SOLAR_SPECTRA_FILE] [--device-absorption-file DEVICE_ABSORPTION_FILE] [--sqcm SQCM]

Shockley-Queisser calcs for an ideal solar cell (n=1, no parasitic resistances)

options:
  -h, --help            show this help message and exit
  --t-cell T_CELL       Temperature of the solar cell [deg C] (default: 25)
  --band-gap BAND_GAP   Band gap of the solar cell [eV] (ignored if --device-absorption-file is given) (default: 1.35)
  --no-plot             Disable plot (default: False)
  --solar-spectra-file SOLAR_SPECTRA_FILE
                        File to read the solar spectra from (default: ASTMG173.csv)
  --device-absorption-file DEVICE_ABSORPTION_FILE
                        File to read the absorption spectrum from (default: )
  --sqcm SQCM           Set the device's area [cm^2] (default: None)
$ ./sq.py --t-cell 50 --band-gap 1.5 --no-plot
We've assumed that our perfect solar cell is at 323.15 degrees kelvin
and has a band gap of 1.5 electron volts.

That means

its radiative saturation current density is
4.142578546268278e-18 mA/cm^2

and if we shine light as defined by the file ASTMG173.csv at it, then

its photocurrent density is
28.954642177247273 mA/cm^2,

which makes

its open circuit voltage
1208.303685366321 mV

its short circuit current density
28.954642177247273 mA/cm^2

the voltage at its maximum power point
1105.1068955914961 mV

the current density at its maximum power point
28.24296505756455 mA/cm^2

its fill factor
89.21138363250819 percent

and

its power production is
31.211495437064265 mW/cm^2 (34.986000851231 if FF was 1.0)
or
31.19993100874486 percent (34.97303789981322 if FF was 1.0).
```
![Some Graph](/figure_1.svg)
