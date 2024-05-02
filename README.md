# lcurvemaker

Python code for working with light curves of variable stars

## Dependencies

* [Matplotlib](https://matplotlib.org)
* [pandas](https://pandas.pydata.org)
* [Astropy](https://www.astropy.org)

## Installation

`pip install matplotlib pandas astropy requests`  
or  
`pip install -r requirements.txt`

## Usage

`python plot_merged_phase_data.py [-v] [-l] [-s] [-c RA DEC] [-p PERIOD] [-e EPOCH] [-z MIN MAX] [-o] nickname`

### Positional arguments

`nickname` is an alias of the object, optionally with the directory name. It is
used to search for files and assign names to data processing products.

The script will search for a settings file named  `nickname.json` in the
`objects` directory. It should not contain spaces.

### Options

* `-h, --help` show help message and exit
* `-v, --verbose` be more verbose
* `-l, --lines` draw lines on the light curve and phase plot to mark the epoch, maximum value, and max/min phase
* `-s, --show` show interactive plots instead of saving figures
* `-c RA DEC, --coord RA DEC` set the coordinates of the object in degrees
* `-p PERIOD, --period PERIOD` set the period for phase plot in days
* `-e EPOCH, --epoch EPOCH` set the epoch for phase plot in [HJD](https://en.wikipedia.org/wiki/Heliocentric_Julian_Day)
* `-z MIN MAX, --ztfran MIN MAX` delete all ZTF data out of range
* `-o, --localps` use local PS1 data instead of requesting it via the API

### Examples

```bash
python plot_merged_phase_data.py gusev4
python plot_merged_phase_data.py minkovskiy24 -l
```

Result:

* [Gusev 4](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227045) [settings](objects/gusev4.json), [phase plot](lc/gusev4-ps1-ztf-atlas-poss-ph.png), [light curve](lc/gusev4-ps1-ztf-atlas.png)
* [Minkovskiy 24](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2387050) [settings](objects/minkovskiy24.json), [phase plot](lc/minkovskiy24-ps1-ztf-ph.png), [light curve](lc/minkovskiy24-ps1-ztf.png)

## Object settings files

The settings for the object are located in a [JSON](https://en.wikipedia.org/wiki/JSON) file.
[sample.json](objects/sample.json) contains most of the possible settings.

Basic settings for an object with [ZTF](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?projshort=ZTF),
[PS1](https://ps1images.stsci.edu/ps1_dr2_api.html), [ASAS-SN](https://asas-sn.osu.edu/),
[ATLAS](https://fallingstar.com/), [CRTS](http://nunuku.caltech.edu/cgi-bin/getcssconedb_priv.cgi),
[Gaia DR3](https://ui.adsabs.harvard.edu/abs/2022yCat.1355....0G),
[SDSS](https://skyserver.sdss.org/dr18/en/tools/search/radial.aspx),
[OGLE](https://ogledb.astrouw.edu.pl/~ogle/OCVS/catalog_query.php) and
[GDS](https://ui.adsabs.harvard.edu/abs/2015AN....336..590H) photometry data available are:

```json
{
"sample-nickname": {
  "name": "ZTF19acdncga",
  "other": "USNO-B1.0 1462-0437198, 2MASS J22470763+5617523, GSC2.3 N1CQ180004",
  "coord": "22 47 07.629 +56 17 52.26",
  "coordeg": [341.7817888, 56.2978492],
  "max": 14.35,
  "min": 14.93,
  "period": 24.8435,
  "epoch": 2460161.091,
  "2ndmin": 0.399,
  "atlasfnam": "sample-nickname-atlas.txt", "atlaslim": {"o": 0.013, "c": 0.017},
  "asasfnam": "sample-nickname-asas.csv", "asaslim": {"V": 0.05, "g": 0.06},
  "pslim": 0.1, "pslocal": false,
  "gaiadr3fnam": "sample-nickname-gaiadr3.dat",
  "crtsfnam": "sample-nickname-crts.csv", "crtslim": 0.06,
  "oglefnam": "sample-nickname-ogle.dat",
  "sdssfnam": "sample-nickname-sdss.dat",
  "curveshift": true,
  "clrshift": {
    "g": 0.3, "psg": 0.32, "r": 0, "psr": 0.01, "i": -0.1, "psi": -0.1, "I": 0.1,
    "psz": -0.2, "psy": -0.3, "o": -0.05, "c": 0.1, "V": 0.3, "asasg": 0.7,
    "G": 0.05, "CV": -0.1, "gr": 0.02, "gi": -0.1
  },
  "plot": {
    "ztffilt": "gri", "psfilt": "grizy", "atlasfilt": "oc", "asasfilt": "gV", "gdsfilt": "r",
    "atlasms": {"o": 2.8, "c": 2.8}, "atlaselw": 0.26, "atledgclr": "darkred",
    "asasms": {"V": 3.5, "g": 3.5}, "asaselw": 0.26,
    "zms": 4, "gms": 4, "psms": 6, "crtsms": 5,
    "xmal": 500, "xmil": 170, "ymal": 0.1, "ymil": 0.01, "leg": "lower left",
    "xlim": [-0.5, 1.0], "xedges": 90, "xlima": [56963, 60331],
    "ylim": [14.7, 13.94], "ylima": [14.88, 13.8]
  }
}}
```

## References for files in the `data` directory

* [Masci, F. J.; et al., 2019, The Zwicky Transient Facility: Data Processing, Products, and Archive](https://ui.adsabs.harvard.edu/abs/2019PASP..131a8003M), [ZTF DR20 query form](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd?catalog=ztf_objects_dr20)
* [Kochanek, C. S.; et al., 2017, The All-Sky Automated Survey for Supernovae (ASAS-SN) Light Curve Server v1.0](https://ui.adsabs.harvard.edu/abs/2017PASP..129j4502K)
* [Tonry, J. L.; et al., 2018, ATLAS: A High-cadence All-sky Survey System](https://ui.adsabs.harvard.edu/abs/2018PASP..130f4505T)
* [Chambers, K. C.; et al., 2016, The Pan-STARRS1 Surveys](https://ui.adsabs.harvard.edu/abs/2016arXiv161205560C)
* [Gaia collaboration; et al., 2022, Gaia Data Release 3 (Gaia DR3) Part 1 Main source](https://ui.adsabs.harvard.edu/abs/2022yCat.1355....0G)
* [Gaia Photometric Science Alerts](http://gsaweb.ast.cam.ac.uk/alerts)

## Examples of light curves (presented in the [VSX](https://aavso.org/vsx/) detail sheets)

### Eclipsing binaries

* [Gusev 4](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227045) [settings](objects/gusev4.json), [phase plot](https://www.aavso.org/vsx_docs/2227045/5780/gusev4-phase-ps1-ztf-atlas.png)
* [Minkovskiy 16](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2225441) [settings](objects/minkovskiy16.json), [phase plot](https://www.aavso.org/vsx_docs/2225441/5780/minkovskiy16-phase-ps1-ztf-atlas.png)
* [GSC 02150-01562](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359207) [phase plot](https://www.aavso.org/vsx_docs/359207/5780/a293-29-asas-gaia-ph_3.422652.png)
* [UCAC3 170-058819](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359210) [phase plot](https://www.aavso.org/vsx_docs/359210/5780/u170-ztf-gaia-ph_3.20341.png)
* [UCAC3 240-187355](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359214) [phase plot](https://www.aavso.org/vsx_docs/359214/5780/u240-ztf-ph_0.699445.png)
* [UCAC3 248-199991](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359217) [phase plot](https://www.aavso.org/vsx_docs/359217/5780/u248-199-ztf-gaia-ps1-ph_1.976275.png)
* [UCAC3 248-205306](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359219) [phase plot](https://www.aavso.org/vsx_docs/359219/5780/u248-ztf-ps1-ph_1.328946.png)
* [USNO-B1.0 1534-0126575](https://www.aavso.org/vsx/index.php?view=detail.top&oid=256288) [phase plot](https://www.aavso.org/vsx_docs/256288/5780/u1534-0126575-ztf-ps1-gaia-ph_1.853704.png)
* [USNO-B1.0 1534-0125222](https://www.aavso.org/vsx/index.php?view=detail.top&oid=256287) [phase plot](https://www.aavso.org/vsx_docs/256287/5780/u1534-222-ps1-ztf-atlas-gaia-ph_0_1.664525), [light curve](https://www.aavso.org/vsx_docs/256287/5780/u1534-222-ps1-ztf-atlas-gaia.png)
* [ZTF19acdncga](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2344600) [phase plot](https://www.aavso.org/vsx_docs/2344600/5780/19acdncga-ztf-ph_24.8435.png)

### Cataclysmic variables

* [NSV 14686](https://www.aavso.org/vsx/index.php?view=detail.top&oid=53310) [phase plot](https://www.aavso.org/vsx_docs/53310/5780/nsv14686-ztf-atlas-asas-gaia-ph_0_1.2101554), [light curve](https://www.aavso.org/vsx_docs/53310/5780/nsv14686-ztf-atlas-asas-gaia_1.png)

### Pulsating variables

* [DAV V8](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227371) [phase plot](https://www.aavso.org/vsx_docs/2227371/5780/davv8-ps1-ztf-atlas-ogle-ph.png), [light curve](https://www.aavso.org/vsx_docs/2227371/5780/davv8-ps1-ztf-atlas-ogle.png)
* [DAV V10](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227460) [phase plot](https://www.aavso.org/vsx_docs/2227438/5780/davv10-ztf-atlas-ph_378_1_2.png), [light curve](https://www.aavso.org/vsx_docs/2227438/5780/davv10-ps1-ztf-atlas_1_2_3.png)
* [DAV V11](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227460) [phase plot](https://www.aavso.org/vsx_docs/2227460/5780/davv11-ps1-ztf-atlas-gaia-ph_360.png), [light curve](https://www.aavso.org/vsx_docs/2227460/5780/davv11-ps1-ztf-atlas-gaia.png)

## Image optimization applied

* [TinyPNG: WebP, PNG, JPEG optimization](https://tinypng.com/)
* [OptiPNG](https://optipng.sourceforge.net/), see [guide to PNG optimization](https://optipng.sourceforge.net/pngtech/optipng.html)
* [Jpegoptim](https://www.kokkonen.net/tjko/projects.html), [for Windows](https://github.com/XhmikosR/jpegoptim-windows)
* [JPEGoptim + OptiPNG + TinyPNG - image optimization (in russian)](https://open-networks.ru/d/14-jpegoptim-optipng-tinypng-optimizaciya-izobrazenii)
