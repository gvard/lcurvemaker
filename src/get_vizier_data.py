from os import path
from argparse import ArgumentParser
from decimal import Decimal, ROUND_HALF_UP

import colorama
import simplejson as json
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from astroquery.vizier import Vizier
import pandas as pd

from serializer import CompactJSONEncoder


parser = ArgumentParser(description="Python script for getting data using astroquery.vizier")
parser.add_argument("alias", type=str, nargs='?', default="",
                    help="Object alias")
parser.add_argument("-g", "--gdr3", type=int,
                    help="set the Gaia DR3 source identifier")
parser.add_argument("-t", "--tic", type=int, default=0,
                    help="set the TIC identifier")
parser.add_argument("-v", "--vsx", type=int, default=0,
                    help="set the VSXOID")
parser.add_argument("-n", "--name", nargs='*', type=str,
                    help="set the primary source designation")
parser.add_argument("-V", "--verbose", action="store_true",
                    help="be more verbose")

args = parser.parse_args()
alias = args.alias
og3 = args.gdr3
name = " ".join(args.name)
odir = "../objects"

fnam = f"{odir}/{alias}"
add = ""
if path.isfile(f"{fnam}.json"):
    add = "-vi"
    print(f"file {fnam}.json exists, write to {fnam}{add}.json")
    fnam = f"{fnam}{add}"

print(f"{alias}: {name}, Gaia DR3 {og3}, TIC {args.tic}, VSX OID {args.vsx}")


def custom_round(value):
    if value > 900:
        return Decimal(value).quantize(Decimal('1'), rounding=ROUND_HALF_UP)
    elif 200 <= value <= 900:
        return Decimal(value).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP)
    else:
        return Decimal(value).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP)


def query_vizier(oid, cat='I/358/varisum', cols=["Source", "VST", "VEB"], rl=2, verbose=False):
    if verbose:
        print(f"{cat} with ID {oid}, cols {cols}, row_limit {rl}")
    v = Vizier(catalog=cat, columns=cols, row_limit=rl)
    if cat == "B/vsx/vsx":
        result = v.query_constraints(OID=oid)
    elif cat == "IV/39/tic82":
        result = v.query_constraints(TIC=oid)
    else:
        result = v.query_constraints(Source=oid)
    if len(result) == 0 or len(result[0]) == 0:
        print("Object", oid, "not found in", cat)
        return None
    return result[0]

if args.verbose:
    print(f'"name": {name},')
    print(f'"gdr3": {args.gdr3},')
    if args.tic:
        print(f'"tic": {args.tic},')
    if args.vsx:
        print(f'"vsxoid": {args.vsx},')


cat = "I/355/gaiadr3"
cols = ["Source", "RAJ2000", "DEJ2000", "Gmag", "BP-RP"]
rl = 2
res = query_vizier(og3, cat=cat, cols=cols, rl=rl)[0]
if len(res):
    gmagv = Decimal(res["Gmag"]).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP)
    bprpv = Decimal(res["BP-RP"]).quantize(Decimal('0.0001'), rounding=ROUND_HALF_UP)
    rar = Decimal(res["RAJ2000"]).quantize(Decimal('0.000001'), rounding=ROUND_HALF_UP)
    der = Decimal(res["DEJ2000"]).quantize(Decimal('0.000001'), rounding=ROUND_HALF_UP)
    ra = Angle(res["RAJ2000"] * u.degree)
    dec = Angle(res["DEJ2000"] * u.degree)
    coord_str = f"{ra.to_string(unit=u.hourangle, sep=' ', pad=True, precision=3)} {dec.to_string(unit=u.deg, sep=' ', pad=True, precision=2)}"
    if args.verbose:
        print(f'"coord": "{coord_str}",')
        print(f'"coordeg": [{rar}, {der}],')
        print(f'"G": {gmagv},')
        print(f'"BP-RP": {bprpv},')


cat = "I/355/epphot"
cols = ["Source", "TimeG", "Gmag"]  # "TimeBP", "BPmag", "TimeRP", "RPmag", "noisyFlag"
rl = 299
res = query_vizier(og3, cat=cat, cols=cols, rl=rl)
epphot = False
if res:
    epphot = True
    g3phfnam = f"../data/{alias}-g3.dat"
    if not path.isfile(g3phfnam):
        print(f"write Gaia DR3 epoch photometry to {g3phfnam}")
        hjds = pd.array(res["TimeG"])
        mags = res["Gmag"]
        g3p = pd.DataFrame({
            "hjd": pd.array(hjds),
            "mag": pd.array(mags),
        }).dropna().to_csv(g3phfnam, sep=" ", float_format='%.6f', index=False, header=False)
        with open(g3phfnam, 'r') as f:
            lines = f.readlines()
        lines[-1] = lines[-1].rstrip("\r\n")
        with open(g3phfnam, 'w') as f:
            f.writelines(lines)
if args.verbose:
    print(f'"epphot": {str(epphot).lower()},')


cat = "I/358/veb"
cols = ["Source", "TimeRef", "Freq"]
rl = 2
veb = False
vebper = None
res = query_vizier(og3, cat=cat, cols=cols, rl=rl)
if res:
    res = res[0]
    veb = True
    vebper = round(1/float(res["Freq"]), 8)
if args.verbose:
    print(f'"veb": {str(veb).lower()},')
    if veb:
        print(f'"vebper": {vebper},')


cat = "I/358/vst"
cols = ["Source", "Freq"]
rl = 2
vst = False
vstper = None
res = query_vizier(og3, cat=cat, cols=cols, rl=rl)
if res:
    res = res[0]
    vst = True
    vstper = round(1/float(res["Freq"]), 8)
if args.verbose:
    print(f'"vst": {str(vst).lower()},')
    if vst:
        print(f'"vstper": {vstper},')


# Distances to 1.47 billion stars in Gaia EDR3 (Bailer-Jones+, 2021)
cat = "I/352/gedr3dis"
cols = ["Source", "rgeo", "rpgeo"]
rl = 2
res = query_vizier(og3, cat=cat, cols=cols, rl=rl)[0]
if len(res):
    try:
        rgeov = custom_round(res["rgeo"])
    except TypeError:
        rgeov = 0
    try:
        rpgeov = custom_round(res["rpgeo"])
    except TypeError:
        rpgeov = 0
if args.verbose:
    print(f'"r": {rgeov},')
    print(f'"rp": {rpgeov},')


cat = 'I/355/paramp'
cols = ["Source", "Teff", "logg", "[Fe/H]", "Dist", "E(BP-RP)", "GMAG", "SpType-ELS"]  # "Lum-Flame", "Mass-Flame"
rl = 1
res = query_vizier(og3, cat=cat, cols=cols, rl=rl)[0]
if len(res):
    teffv = Decimal(res["Teff"]).quantize(Decimal('1'), rounding=ROUND_HALF_UP)
    teffref = "2023A&A...674A..26C"
    loggv = Decimal(res["logg"]).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP)
    loggref = "2023A&A...674A..26C"
    ebprp = Decimal(float(res["E(BP-RP)"])).quantize(Decimal('0.0001'), rounding=ROUND_HALF_UP)
    absg = Decimal(res["GMAG"]).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP)
    dist = Decimal(res["Dist"]).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP)
    spv = str(res["SpType-ELS"])
else:
    teffv = 0
    teffref = ""
    loggv = 0
    loggref = ""
    ebprp = 0
    absg = 0
    dist = 0
    spv = ""
gmg = float(gmagv)
if args.verbose and spv:
    print(f'"sp": "{spv}",')
    print('"spref": "2023A&A...674A..26C",')
if args.verbose and teffv:
    print(f'"teff": "{teffv}",')
    print('"teffref": "2023A&A...674A..26C",')
if args.verbose and loggv:
    print(f'"logg": "{loggv}",')
    print('"loggref": "2023A&A...674A..26C",')


if args.tic:
    cat = "IV/39/tic82"
    cols = ["TIC", "RAJ2000", "DEJ2000", "Teff", "logg"]
    rl = 2
    res = query_vizier(args.tic, cat=cat, cols=cols, rl=rl)[0]
    if res:
        try:
            tteffv = Decimal(float(res["Teff"])).quantize(Decimal('1'), rounding=ROUND_HALF_UP)
            tteffref = "2019AJ....158..138S"
        except TypeError:
            tteffv = 0
            tteffref = ""
        tloggv = Decimal(float(res["logg"])).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP)
        tloggref = "2019AJ....158..138S"
        trar = Decimal(res["RAJ2000"]).quantize(Decimal('0.000001'), rounding=ROUND_HALF_UP)
        tder = Decimal(res["DEJ2000"]).quantize(Decimal('0.000001'), rounding=ROUND_HALF_UP)
        tra = Angle(res["RAJ2000"] * u.degree)
        tdec = Angle(res["DEJ2000"] * u.degree)
        tcoord_str = f"{tra.to_string(unit=u.hourangle, sep=' ', pad=True, precision=3)} {tdec.to_string(unit=u.deg, sep=' ', pad=True, precision=2)}"
    else:
        tteffv = 0
        tloggv = 0
        trar = 0
        tder = 0
        tra = 0
        tdec = 0
        tcoord_str = ""
    if args.verbose and tteffv:
        print(f'"tteff": "{tteffv}",')
        print(f'"tteffref": {tteffref},')
    if args.verbose and tloggv:
        print(f'"tlogg": "{tloggv}",')
        print(f'"tloggref": {tloggref},')


data_json = {alias:
{
    "name": name,
    "other": "",
    "desig": "",
    "gdr3": og3,
    "tic": args.tic,
}}
data_json[alias]["epphot"] = epphot
data_json[alias]["veb"] = veb
if vebper:
    data_json[alias]["vebper"] = vebper
data_json[alias]["vst"] = vst
if vstper:
    data_json[alias]["vstper"] = vstper
data_json[alias]["coord"] = coord_str
data_json[alias]["coordeg"] = [rar, der]
data_json[alias]["type"] = "EA"
data_json[alias]["max"] = 0
data_json[alias]["min"] = 0
data_json[alias]["system"] = ""
data_json[alias]["period"] = 0
data_json[alias]["epoch"] = 0
if args.tic and coord_str != tcoord_str:
    data_json[alias]["tcoord"] = tcoord_str
if args.tic and (abs(rar - trar) > 0.000002 or abs(der - tder) > 0.000002):
    data_json[alias]["tcoordeg"] = [trar, tder]
data_json[alias]["G"] = gmagv
data_json[alias]["BP-RP"] = bprpv
if ebprp:
    data_json[alias]["E(BP-RP)"] = ebprp
data_json[alias]["etype"] = "min"
data_json[alias]["sp"] = spv
if spv:
    data_json[alias]["spref"] = "2023A&A...674A..26C"
if teffv:
    data_json[alias]["teff"] = teffv
    data_json[alias]["teffref"] = teffref
if loggv:
    data_json[alias]["logg"] = loggv
    data_json[alias]["loggref"] = loggref
if args.tic and tteffv:
    data_json[alias]["tteff"] = tteffv
    data_json[alias]["tteffref"] = tteffref
if args.tic and tloggv:
    data_json[alias]["tlogg"] = tloggv
    data_json[alias]["tloggref"] = tloggref
data_json[alias]["r"] = rgeov
data_json[alias]["rp"] = rpgeov
data_json[alias]["comm"] = ""


if args.vsx:
    data_json[alias]["vsxoid"] = args.vsx
    cat = "B/vsx/vsx"
    cols = ["OID", "RAJ2000", "DEJ2000", "Name", "V", "Type", "max", "n_max", "min", "n_min", "Period", "Epoch", "Sp"]
    rl = 2
    res = query_vizier(args.vsx, cat=cat, cols=cols, rl=rl)[0]
    if res:  # and args.verbose:
        print(colorama.Fore.RED + f'"vsxnam": "{res["Name"]}",' + colorama.Style.RESET_ALL)
        data_json[alias]["vsxnam"] = res["Name"]
        print(colorama.Fore.RED + f'"vsxfl": {res["V"]},' + colorama.Style.RESET_ALL)
        data_json[alias]["vsxfl"] = int(res["V"])
        print(colorama.Fore.RED + f'"vsxcrd": [{Decimal(str(res["RAJ2000"]))}, {Decimal(str(res["DEJ2000"]))}],' + colorama.Style.RESET_ALL)
        vra = Angle(res["RAJ2000"] * u.degree)
        vdec = Angle(res["DEJ2000"] * u.degree)
        vcoord_str = f"{vra.to_string(unit=u.hourangle, sep=' ', pad=True, precision=2)} {vdec.to_string(unit=u.deg, sep=' ', pad=True, precision=1)}"
        print(colorama.Fore.RED + f'"vsxcr": "{vcoord_str}",' + colorama.Style.RESET_ALL)
        data_json[alias]["vsxcr"] = vcoord_str
        data_json[alias]["vsxcrd"] = [Decimal(str(res["RAJ2000"])), Decimal(str(res["DEJ2000"]))]
        if res["Type"]:
            print(colorama.Fore.RED + f'"vsxtyp": "{res["Type"]}",' + colorama.Style.RESET_ALL)
            data_json[alias]["vsxtyp"] = res["Type"]
        if not (np.isnan(res["max"]) or str(res["max"]) in ("--", "?")):
            print(colorama.Fore.RED + f'"vsxmax": {Decimal(str(res["max"]))},' + colorama.Style.RESET_ALL)
            data_json[alias]["vsxmax"] = Decimal(str(res["max"]))
        if res["n_max"] and str(res["n_max"]) not in ("--", "?"):
            print(colorama.Fore.RED + f'"vsxmaxpb": {str(res["n_max"])},' + colorama.Style.RESET_ALL)
            data_json[alias]["vsvsxmaxpb"] = str(res["n_max"])
        if not (np.isnan(res["min"]) or str(res["min"]) in ("--", "?")):
            print(colorama.Fore.RED + f'"vsxmin": {Decimal(str(res["min"]))},' + colorama.Style.RESET_ALL)
            data_json[alias]["vsxmin"] = Decimal(str(res["min"]))
        if not (np.isnan(res["Period"]) or str(res["Period"]) == "--"):
            print(colorama.Fore.RED + f'"vsxper": {Decimal(str(res["Period"]))},' + colorama.Style.RESET_ALL)
            data_json[alias]["vsxper"] = Decimal(str(res["Period"]))
        if not (np.isnan(res["Epoch"]) or str(res["Epoch"]) == "--"):
            print(colorama.Fore.RED + f'"vsxep": {Decimal(str(res["Epoch"]))},' + colorama.Style.RESET_ALL)
            data_json[alias]["vsxep"] = Decimal(str(res["Epoch"]))
        if res["Sp"]:
            print(colorama.Fore.RED + f'"vsxsp": {res["Sp"]},' + colorama.Style.RESET_ALL)
            data_json[alias]["vsxsp"] = res["Sp"]


data_json[alias]["lc"] = ""
data_json[alias]["lczoom"] = ""
data_json[alias]["lcs"] = []
data_json[alias]["curveshift"] = True
data_json[alias]["clrshift"] = {}
data_json[alias]["data"] = {}
data_json[alias]["errlim"] = {}
data_json[alias]["ms"] = {}
data_json[alias]["plot"] = {}

with open(f"{odir}/{alias}{add}.json", "w") as f:
    json.dump(data_json, f, cls=CompactJSONEncoder, ensure_ascii=False, indent=0)
