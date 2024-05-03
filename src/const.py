DATA_PATH = "../data/"
DATA_RAW = f"{DATA_PATH}raw/"
EXT = "png"
HTML_TABLE_URL = "../index.html"
TYPES = {
    "ecl": ("EA", "EA/WD"),
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

ATLAS_RAWDATA = True
ATLAS_SAVEDATA = True
ZTF_RAWDATA = True

GAIADR3_JD_SHIFT = 2455197.5

PS_CONE_RADIUS = 0.53

PSFILT = "griz"
ZTFFILT = "gri"
ATLASFILT = "oc"
ASASFILT = "Vg"

# Constants for filtering data by magnitude error
PS_MAG_LIM = 0.11
ASAS_MAG_LIM = 0.75
ATLAS_MAG_LIM = 0.19
ATLAS_MAG_LIMA = {"o": 0.19, "c": 0.15}

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
ZEROFILTSHIFT = {"g": 0, "psg": 0, "r": 0, "psr": 0, "i": 0, "I": 0, "psi": 0,
                 "z": 0, "psz": 0, "psy": 0, "o": 0, "c": 0, "V": 0, "G": 0, "asasg": 0,
                 }

# Marker colors for given filters
FILTER_COLORS = {
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
    "gi": "brown",
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

ELINWDTH = 0.5
PSELINWDTH = 1
ZELINWDTH = 1
ATELINWDTH = {"o": 0.5, "c": 0.3}
ATELINWDTHT = 0.3
ASELINWDTH = 0.3
CELINWDTH = 0.3
OELINWDTH = 0.8
GDSELINWDTH = 1

MEW = 0.4

TITLEFNTSZ = 18
LEGFNTSZ = 14

ATLAS_MAGLIM_UP = 10
