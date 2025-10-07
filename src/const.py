DATA_PATH = "../data/"
DATA_RAW = f"{DATA_PATH}raw/"
EXT = "png"
HTML_TABLE_URL = "../index.html"
TYPES = {
    "ecl": ("EA", "EA/WD", "EA/HW", "EA+HB"),
    "cv": (
        "AM", "CV", "NL", "NL+E", "NL/VY", "SN", "SN|UG", "UG", "UG+VY",
        "UGSU", "UGSU+E", "UGZ", "UGZ/IW", "UGZ/IW+VY", "UG|YSO", "VAR",
    ),
    "pls": (
        "L", "M", "SR", "SRA", "SRB",
    ),
    "er": (
        "UV", "UV+BY",
    ),
}
AT_URL = "https://www.wis-tns.org/object/"
VSX_URL = "https://www.aavso.org/vsx/index.php?view=detail.top&amp;oid="

ABBR = {
    "ZTF": "ztf",
    "Pan-STARRS1": "ps1",
    "ASAS-SN": "asn",
    "ASAS-3": "as3",
    "ATLAS": "atl",
    "CRTS": "css",
    "Gaia DR3": "g3",
    "Gaia Alerts": "gaia",
    "GDS": "gds",
    "SuperWASP": "wsp",
    "OGLE": "ogl",
    "MASCARA": "msc",
    "SDSS": "sdss",
    "Hipparcos": "hp",
    "TESS": "tess",
    "CoRoT": "cr",
    "Kepler": "kp",
    "K2": "k2",
    "DSS": "dss",
}

DEFAULTS = {
    "zcatf": "gr",
    "pslocal": False,
    "curveshift": True,
    "clrshift": {
        "u": 0, "g": 0, "psg": 0, "r": 0, "psr": 0, "i": 0, "I": 0, "psi": 0,
        "z": 0, "psz": 0, "psy": 0, "o": 0, "c": 0, "V": 0, "AS3V": 0, "G": 0,
        "asng": 0, "gr": 0, "gi": 0, "W": 0, "CV": 0, "MSC": 0, "TESS": 0,
        "Kp": 0, "K2": 0, "CRT": 0, "Hp": 0,
    },
    "errlim": {  # Constants for filtering data by magnitude error
        "ps1": 0.5, "asn": {"V": 0.5, "g": 0.5}, "atl": {"o": 0.6, "c": 0.6},
        "crts": 0.5, "as3": 0.07, "msc": 0.1, "wsp": 0.1, "hip": 0.2,
    },
    "plot": {  # What data to plot
        "ztf": "gri", "ps1": "grizy", "asn": "gV", "atl": "oc", "ogl": "", "gds": "ri",
        "crts": True, "as3": True, "msc": True, "wsp": True, "g3": True, "tess": "QLP",
        "hip": True, "kp": True, "k2": True, "crtruns": [],
        # Plot settings
        "xmal": 500, "xmil": 100,
        "ymil": 0.01, "ymal": 0.02,
        "xedges": 20,
        "xlim": [-0.5, 1],
        "atledgclr": "darkred", "atlelw": 0.26,
        "ms": {  # Marker sizes
            "ztf": {"g": 6, "r": 6, "i": 6}, "ps1": 7, "asn": {"g": 4, "V": 5}, "atl": {"o": 4, "c": 4},
            "ogl": {"I": 6, "V": 6}, "gds": {"r": 8, "i": 8}, "crts": 5, "as3": 7,
            "msc": 4.5, "wsp": 3.5, "dss": 9, "gaia": 7, "g3": 7, "tess": 3.5, "hip": 8.5, "kp": 5, "crt": 3.5,
        },
    },
    "zoom": {
        "xlim": [-0.15, 0.15],
        "xmalp": 0.1,
        "xmilp": 0.05
    },
}

ATLAS_RAWDATA = True
ATLAS_SAVEDATA = True
ZTF_RAWDATA = True

GAIADR3_JD_SHIFT = 2455197.5

PS_CONE_RADIUS = 0.53

PSFILT = "grizy"
ZTFFILT = "gri"
ATLASFILT = "oc"
ASASFILT = "Vg"
GDSFILT = "ri"

DATA_TO_MERGE = {
    "ZTF": "gri",
    "ATLAS": "oc",
    "ALeRCE": "ri",
    "OGLE": "I",
    "ASAS": "Vg",
    "PS1": "grizy",
    "Gaia": "G"
}

# Plotting parameters
ZEROFILTSHIFT = {"u": 0, "g": 0, "psg": 0, "r": 0, "psr": 0, "i": 0, "I": 0, "psi": 0,
                 "z": 0, "psz": 0, "psy": 0, "o": 0, "c": 0, "V": 0, "G": 0, "asasg": 0,
                 "gr": 0, "gi": 0,
                 }

# Marker colors for given filters
FILTER_COLORS = {
    "u": "#00f",
    "g": "#080",
    "r": "#f33",
    "i": "#d61fff",
    "zg": "#070",
    "zr": "#f00",
    "zi": "#d61fff",
    "asasg": "#ff8300",
    "V": "#00cc68",
    "I": "#1f71ff", "z": "k", "crts": "#666",
    "ps_g": "#009F72",
    "ps_r": "#f88", "ps_i": "#33f", "ps_z": "k", "ps_y": "darkgrey",
    "y": "darkgrey",
    "gr": "#f00",
    "gi": "#730bdb",
    "o": "darkorange",  # "gold",
    "c": "#11f",
    "G": "#2f2",
    "ALeRCE r": "r",
    }
FILT_CLRS = {
    "g": "g", "asasg": "lightgreen", "r": "r", "i": "m", "I": "indigo", "z": "k",
    "y": "darkgrey", "crts": "#666", "o": "darkorange", "c": "darkblue",
    "G": "darkgrey", "V": "lime", "ALeRCE r": "maroon",
    }

MARKERS = "DpshX"
# Marker sizes
PSMS = 5
ZMS = 3
ASASMS = 3
ATLASMS = {"o": 3, "c": 2}
CRTSMS = 4
OGLEMS = 5
GAIAMS = 7.5
GDSMS = 7
DSSMS = 9

ELINWDTH = 0.5
PSELINWDTH = 1
ZELINWDTH = 0.3
ATELINWDTH = {"o": 0.5, "c": 0.3}
ATELINWDTHT = 0.3
ASELINWDTH = 0.3
CELINWDTH = 0.3
OELINWDTH = 0.8
GDSELINWDTH = 1

LW_MIN = 0.6

XMAL = 500
XMIL = 100
YMAL = 0.5
YMIL = 0.02
XEDGES = 100
XLIMP = [-0.5, 1]
XMALP = 0.5
XMILP = 0.05

MEW = 0.4
CMEW = 0.5
GMEW = 0.8
DMS = 0.8

TITLEFNTSZ = 18
LEGFNTSZ = 14

ATLAS_MAGLIM_UP = 10

PLT_EDGES = [0.06, 0.09, 0.985, 0.95]
FIGSZ = [16, 9]
DPI = 120
DPI_SCREEN = 125
