# lcurvemaker object settings files

## Names and designations

```jsonc
"gdr3": 3120034358879078016,
"tic": 43319633,
"hip": 25316,
"ascc": 1203076,
"hd": 74588,
"kic": 9153621,
"epic": 247612547,
"corot": 315189250,
```

## Gaia DR3 related settings

```jsonc
"epphot": true,  // Presence of Gaia DR3 epoch photometry
"veb": true,  // Presence in Gaia DR3 eclipsing binaries data files
"vst": false,  // Presence in Gaia DR3 short-timescale sources data files
"G": 7.7493,  // G magnitude from Gaia DR3
"BP-RP": 0.2086,  // BP-RP color from Gaia DR3
"E(BP-RP)": 0.0287,  // Reddening, from I/355/paramp table
"r": 115.63,  // Median of the geometric distance posterior, pc; Bailer-Jones+, 2021
"rp": 116.65,  // Median of the photogeometric distance posterior, pc; Bailer-Jones+, 2021
```

## Data

```jsonc
"errlim": {"as3": 0.039, "asn": {"g": 0.015, "V": 0.025}, "ps1": 0.1, "crts": 0.1},  // Limitations for filtering data based on accuracy
"data": {"as3": "V", "asn": "gV"},  // Data to be saved and merged
```

## Plotting

```jsonc
"plot": {
  "ztf": "gri", "ps1": "grizy", "atl": "oc", "as3": "V", "asn": "gV", "gds": "ri", "tess": "QLP",  // What data to plot
  "ms": {"as3": 4.5, "tess": 5, "ztf": {"g": 4, "r": 4, "i": 4}, "g3": 8, "ps1": 6, "atl":{"o": 3.5, "c": 3.5}, "asn": {"g": 3, "V": 4}},  // Matplotlib marker sizes on plots, markersize setting
  "leg": "lower left",  // The location of a legend within a plot (Matplotlib loc parameter)
  "xmal": 500,  // Sets the locator of the major ticker on the x-axis
  "xmil": 100,  // Sets the locator of the minor ticker on the x-axis
  "ymal": 0.002,  // Sets the locator of the major ticker on the y-axis
  "ymil": 0.0005,  // Sets the locator of the minor ticker on the y-axis
  "xedges": 50,  // The value of trimming range of x-axis limits in days, controls the x-axis view limits of a JD plot
  "xlim": [-0.3, 1],  // The range for the x-axis of the phase plot
  "xlima": [60663.9, 60687.8],  // The range for the x-axis of the plot
  "ylim": [7.911, 7.733],  // The range of magnitudes displayed on the phase plot
  "ylima": [7.924, 7.755],  // The range of magnitudes displayed on the JD plot
}

"zoom": {
  "xlim": [0.518, 1.035],  // The phase range displayed on the zoomed phase plot
  "ylim": [12.723, 12.575],  // The magnitude range displayed on the zoomed phase plot
  "xlim1": [-0.045, 0.045],  // The phase range displayed in the first subplot of zoomed phase plot with two subplots
  "xlim2": [0.5134, 0.6034],  // The phase range displayed in the second subplot of zoomed phase plot with two subplots
  "xmalp": 0.02,  // Sets the locator of the major ticker on the x-axis of the subplots
  "xmilp": 0.002,  // Sets the locator of the minor ticker on the x-axis of the subplots
}
```

## TESS specific settings

```jsonc
"tess": {
  "sect": [59, 73, 19],  // Sector numbers to be used. The data will be plotted on the graph in the specified order.
  "mags": {"52": [8.538, 0], "87": [8.538, 10]},  // Magnitude intervals for the specified sectors to be cut. Zero means that this boundary is not used.
  "cuts": {"52": [59721.9, 59723.64], "87": [60688.97, 0]},  // Time intervals in MJD for the specified sectors to be cut. Zero means that this boundary is not used.
}
"plot": {
  "tessclrs": ["r", "g", "b", "#777", "gold", "darkorange", "lawngreen", "#04D8B2", "darkmagenta"],  // The colors that are used to plot TESS data
}
```

## An example of some notable settings

```jsonc
{"hd54896": {
  "name": "HD 54896",  // Primary designation. Note that this is not necessarily the variable name in GCVS or in the NSV catalog.
  "other": "HIP 34806, TIC 38215399, AGASC 322837264, BD+36 1582, ASCC 483473",  // Other object designations, used for the second header
  "desig": "UCAC4 631-040638, 2MASS J07121855+3607377, GSC 02463-00866, SAO 59882, TYC 2463-866-1",  //  Designations not contained in the other field
  "gdr3": 897789765642169984,  // Gaia DR3 source identifier
  "tic": 38215399,  // TESS Input Catalog identifier
  "ascc": 483473,  // ASCC-2.5 number
  "epphot": true,  // Presence of Gaia DR3 epoch photometry
  "veb": true,  // Presence in Gaia DR3 eclipsing binaries data files
  "vst": false,  // Presence in Gaia DR3 short-timescale sources data files
  "coord": "07 12 18.535 +36 07 37.56",  // Coordinates (J2000)
  "coordeg": [108.077229, 36.127100],  // RA, DEC coordinates in degrees
  "type": "EA",  // The object's variability type
  "max": 7.76,  // Magnitude in maximum brightness
  "min": 7.91,  // Magnitude in minimum brightness
  "amp": 0.15,  // Amplitude of variations (eclipse depth for an eclipsing binary)
  "mint": 0.055,  // Totality of the primary eclipse
  "min2": 7.786,  // Magnitude of the secondary minimum (for an eclipsing binary)
  "amp2": 0.008,  // Secondary eclipse depth for an eclipsing binary
  "d": 0.067,  // Primary eclipse duration as a fraction of the period (for an eclipsing binary)
  "d2": 0.059,  // Secondary eclipse duration as a fraction of the period (for an eclipsing binary)
  "min2t": 0.025,  // Totality of the secondary eclipse
  "min2ph": 0.5305,  // Phase of the secondary eclipse
  "system": "V",  // The photometric passband of the magnitudes
  "period": 2.433558,  // The variability period in days
  "epoch": 2459580.242,  // Epoch, HJD
  "ell": 0.004,  // Ellipsoidal variations amplitude, mag
  "ellsys": "TESS",  // Ellipsoidal variations passband
  "G": 7.7493,  // G magnitude from Gaia DR3
  "BP-RP": 0.2086,  // BP-RP color from Gaia DR3
  "E(BP-RP)": 0.0287,  // Reddening, from I/355/paramp table
  "V": 7.78,  // V magnitude
  "sp": "A2",  // Spectral type
  "spref": "1993yCat.3135....0C",  // Reference for spectral type
  "teff": 8165,  // Effective temperature estimation
  "teffref": "2019AJ....158..138S",  // Reference for the effective temperature value
  "logg": 4.27,  // Surface gravity estimation
  "loggref": "2019AJ....158..138S",  // Reference for the surface gravity value
  "comm": "Eccentric system. Min II at phase 0.5305, amplitude 0.008 TESS, duration 5.9%. Min II is Total. Ellipsoidal variations with amplitude 0.004 TESS.",  // Comments
  "r": 115.63,  // Median of the geometric distance posterior, pc; Bailer-Jones+, 2021
  "rp": 116.65,  // Median of the photogeometric distance posterior, pc; Bailer-Jones+, 2021
  "vsxoid": 249919,  // AAVSO VSX object identifier
  "lclnk": "https://vsx.aavso.org/vsx_docs/249919/5780/hd54896-hip-g3-tess-ph.png",  // Link to light curve image
  "lczoom": "https://vsx.aavso.org/vsx_docs/249919/5780/hd54896-g3-hip-ts-phz.png",  // Link to zoomed light curve image
  "lcs": ["https://vsx.aavso.org/vsx_docs/249919/5780/hd54896-g3-hip-ts-phzmax.png"],  // An array of links to light curves
  "curveshift": true,  // Is the shift of the light curves in different bandwidths applied
  "clrshift": {"g": 0, "r": 0, "i": 0, "Kp": -0.06, "TESS": 0.035, "V": 0, "G": -0.02, "Hp": 0.053, "MSC": 0.05},  // The shift of the light curves in different bandwidths
  "plot": {
    "xmal": 500,  // Sets the locator of the major ticker of x-axis
    "xmil": 100,  // Sets the locator of the minor ticker of x-axis
    "ymal": 0.002,  // Sets the locator of the major ticker of y-axis
    "ymil": 0.0005,  // Sets the locator of the minor ticker of y-axis
    "leg": "lower left",  // The location of a legend within a plot (Matplotlib loc parameter)
    "xedges": 50,  // The value of trimming range of x-axis limits in days, controls the x-axis view limits of a JD plot
    "xlim": [-0.3, 1],  // The range for the x-axis of the phase plot
    "ylim": [7.911, 7.733],  // The range of magnitudes displayed on the phase plot
    "ylima": [7.924, 7.755]  // The range of magnitudes displayed on the JD plot
  }
}}
```

Note the trailing zeros of numeric fields, they are read using the [simplejson](https://pypi.org/project/simplejson/) library.
