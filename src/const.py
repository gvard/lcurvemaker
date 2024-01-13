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

PSFILT = "griz"
ZTFFILT = "gri"
ATLASFILT = "oc"
ASASFILT = "Vg"

# Constants for filtering data by magnitude error
PS_MAG_LIM = 0.11
ASAS_MAG_LIM = 0.75
ATLAS_MAG_LIM = 0.19


# Plotting parameters
ZEROFILTSHIFT = {"g": 0, "psg": 0, "r": 0, "psr": 0, "i": 0, "I": 0, "psi": 0,
             "z": 0, "psz": 0, "o": 0, "c": 0, "V": 0, "G": 0, "asasg": 0,
             }

# Marker colors for given filters
FILT_CLRS = {"g": "g", "asasg": "lightgreen", "r": "r", "i": "m", "I": "indigo", "z": "k",
             "y": "darkgrey", "o": "darkorange", "c": "darkblue", "G": "darkgrey",
             }

MARKERS = "DpshX"
# Marker sizes
PSMS = 5
ZMS = 3
ASASMS = 3
ATLASMS = 3
CRTSMS = 3
OGLEMS = 3

ELINWDTH = 0.5
PSELINWDTH = 1
ZELINWDTH = 0.75
ATELINWDTHT = 0.3
ATELINWDTH = 0.5
ASELINWDTH = 0.75
CELINWDTH = 0.8
OELINWDTH = 0.8
