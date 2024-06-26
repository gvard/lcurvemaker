"""Python script for working with light curves of variable stars
"""

import warnings
from sys import exit
import json
from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
from astropy.utils.exceptions import AstropyWarning

import const
from curves import (
    Ztf,
    Zobj,
    # Alerce,
    Ps,
    Atlas,
    Ogle,
    Crts,
    read_gds_data,
    read_ps_data,
    mk_phased,
    mk_fsh,
    save_datafile,
    save_merged,
    save_gaia_datafile,
    mk_hjd_corr,
)


warnings.simplefilter("ignore", category=AstropyWarning)
pd.options.mode.chained_assignment = None  # default="warn"

parser = ArgumentParser(
    description="Python script for working with light curves of variable stars"
)
parser.add_argument("nickname", type=str, default="",
                    help="alias of the object, optionally with the directory name")
parser.add_argument("savedir", type=str, nargs="?", default="",
                    help="set default directory for saving plots")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="be more verbose")
parser.add_argument("-l", "--lines", action="store_true",
    help="draw lines on the light curve and phase plot to mark the epoch, and max/min values and phases")
parser.add_argument("-s", "--show", action="store_true",
                    help="show interactive plots instead of saving figures")
parser.add_argument("-c", "--coord", nargs="+", type=float,
                    help="set the coordinates of the object in degrees")
parser.add_argument("-p", "--period", type=float,
                    help="set the period for phase plot in days")
parser.add_argument("-e", "--epoch", type=float,
                    help="set the epoch for phase plot")
parser.add_argument("-r", "--ztfran", nargs="+", type=float,
                    help="delete all ZTF data out of range")
parser.add_argument("-o", "--localps", action="store_true",
                    help="use local PS1 data instead of requesting it via the API")
parser.add_argument("-z", "--zoom", action="store_true",
                    help="use settings for zoomed plot")
parser.add_argument("-m", "--model", action="store_true",
                    help="draw a simple light curve model")
parser.add_argument("-t", "--plot", nargs="+", type=str,
    help="what data to plot. Possible values are: zt ps as at cs ga og gd")
args = parser.parse_args()

try:
    # Try to read settings file, otherwise use args
    with open(f"../objects/{args.nickname}.json", "r", encoding="utf8") as loc_file:
        settings = json.load(loc_file)
    name = next(iter(settings))
    # Use data filenames prefix as object name from settings
    obj = settings[name]
    if "plot" not in obj:
        obj["plot"] = {}
    # Set title for plots
    titlename = obj["name"] if "name" in obj else name
    name_add = obj["other"] if "other" in obj else ""
except FileNotFoundError:
    if args.nickname and args.coord:
        print(
            f"Settings file not found, continue with {args.nickname} name and"
            f"{args.coord} coordinates"
        )
        name = args.nickname
        obj = {"plot": {}}
        titlename = name
        name_add = ""
    else:
        exit("Specify the name of the settings file or the name of the object and its coordinates")

objtitle = f"{titlename} = {name_add}" if name_add else titlename
fnam = name.lower().replace(" ", "")
objpsfnam = f"{fnam}-ps"

if args.verbose:
    print("args:", args)

if args.ztfran:
    obj["ztfran"] = args.ztfran

if args.plot and "zt" in args.plot:
    const.ZTFFILT = "gri"
elif obj["plot"].get("ztffilt"):
    const.ZTFFILT = obj["plot"].get("ztffilt")
else:
    const.ZTFFILT = ""
if args.plot and "ps" in args.plot:
    const.PSFILT = "grizy"
else:
    const.PSFILT = ""
if args.plot and "as" in args.plot:
    const.ASASFILT = "gV"
else:
    const.ASASFILT = ""
if args.plot and "at" in args.plot:
    const.ATLASFILT = "oc"
elif obj["plot"].get("atlasfilt"):
    const.ATLASFILT = obj["plot"].get("atlasfilt")
else:
    const.ATLASFILT = ""
if args.plot and "gd" in args.plot:
    const.GDSFILT = "ri"
elif obj["plot"].get("gdsfilt"):
    const.GDSFILT = obj["plot"].get("gdsfilt")
else:
    const.GDSFILT = ""

if args.coord:
    ra, dec = args.coord
else:
    ra, dec = obj["coordeg"]
if args.period:
    period = args.period
else:
    period = obj.get("period")
if args.epoch:
    epoch = args.epoch
else:
    epoch = obj.get("epoch")

JD_SHIFT = 2400000.5

filtshift = const.ZEROFILTSHIFT
curveshift = filtshift
if obj.get("clrshift"):
    filtshift.update(obj.get("clrshift"))
if obj.get("curveshift") and obj.get("clrshift"):
    curveshift.update(obj.get("clrshift"))

if obj["plot"].get("crtsms"):
    const.CRTSMS = obj["plot"]["crtsms"]

print(f"-== {name} ==-")

if "psfilt" in obj["plot"]:
    const.PSFILT = obj["plot"]["psfilt"]
    if obj["plot"].get("psms"):
        const.PSMS = obj["plot"]["psms"]
    MSS = [const.PSMS - 0.5, const.PSMS + 1, const.PSMS, const.PSMS + 1, const.PSMS]

# Get PS1 data
if (
    not args.localps
    and obj.get("plot")
    and const.PSFILT
    and not obj.get("pslocal")
):
    Psdat = Ps(ra, dec, radius=const.PS_CONE_RADIUS)  # , radius=0.09
    Psdat.cone_search()
    if obj["plot"].get("psfilt"):
        Psdat.PSFILTRS = const.PSFILT
    if Psdat.results:
        Psdat.get_table()
    else:
        print(f"No PS1 data for {name}")

# Light curve figure
fig, ax = plt.subplots(figsize=const.FIGSZ)
fig.subplots_adjust(*const.PLT_EDGES)

# ASAS Stuff
if "asasfilt" in obj["plot"]:
    const.ASASFILT = obj["plot"].get("asasfilt")
if obj.get("asasfnam"):
    asasd = pd.read_csv(f"{const.DATA_RAW}{obj['asasfnam']}")
    asasd = asasd.rename(
        columns={"HJD": "hjd", "mag_err": "magerr", "Filter": "filter"}
    )
    # For plotting upper limits
    if obj.get("asuplim"):
        asasd_bad = asasd[asasd["magerr"] == 99.99]
        try:
            asasd_bad = asasd_bad[asasd_bad["mag"] < 21]
        except TypeError:
            asasd_bad["mag"] = asasd_bad["mag"].str.lstrip(">").astype(float)
            asasd_bad = asasd_bad[asasd_bad["mag"] < 21]
    # Strip upper limits
    asasd["mag"] = asasd["mag"].astype("str").str.lstrip(">").astype(float)
    ASASERRLIM = obj.get("asaslim") if obj.get("asaslim") else const.ASAS_MAG_LIM
    asasd = asasd[asasd["magerr"] != 99.99]
    # asasd = asasd[asasd["magerr"] < ASASERRLIM]
    asdata = {}
    for filtr in const.ASASFILT:
        asdata[filtr] = asasd[asasd["filter"] == filtr]
        asdata[filtr]["mag"] = asdata[filtr]["mag"].astype(float)
        asdata[filtr] = asdata[filtr][asdata[filtr]["magerr"] < ASASERRLIM[filtr]]
        if obj.get("asasuplim"):
            asdata[filtr] = asdata[filtr][asdata[filtr]["mag"] >= obj.get("asasuplim") - curveshift[filtr]]
        if period:
            asdata[filtr] = mk_phased(asdata[filtr], epoch, period, jdnam="hjd")

    try:
        print(
            f"ASAS range for {name}: {min(asdata['V']['mag'])}-{max(asdata['V']['mag'])} V,",
            f"{min(asdata['g']['mag'])}-{max(asdata['g']['mag'])} g",
        )
    except (ValueError, KeyError):
        pass
    save_datafile(asasd, f"{const.DATA_PATH}{fnam}-asassn.dat", hjdprec=6, magprec=3)

    if "asasms" in obj["plot"]:
        const.ASASMS = obj["plot"].get("asasms")
    pd.DataFrame({
        "hjd": asasd["hjd"].astype("str").str.ljust(13, "0"),
        "mag": asasd["mag"].round(6).astype("str").str.ljust(6, "0"),
        "magerr": asasd["magerr"].astype("str").str.ljust(5, "0"),
        "filter": asasd["filter"],
    }).to_csv(f"{const.DATA_PATH}{fnam}-asassn.dat", sep=" ", index=False)
    if obj["plot"].get("asasupperlim") and len(asasd_bad["mag"]):
        plt.plot(asasd_bad["hjd"] - JD_SHIFT, asasd_bad["mag"], marker="v", ls="none",
                 c="#bbb", ms=2, zorder=-2)  # , label="ASAS-SN upper limits"
    if "g" in const.ASASFILT and len(asdata["g"]["mag"]):
        plt.errorbar(
            asdata["g"]["hjd"] - JD_SHIFT, (asdata["g"]["mag"] - curveshift["asasg"]), asdata["g"]["magerr"],
            marker="d", ls="none", c=const.FILTER_COLORS["asasg"], markeredgecolor="k", zorder=-1,
            mew=const.MEW, label=f"ASAS-SN g{mk_fsh(curveshift['asasg'])}", ms=const.ASASMS,
            elinewidth=const.ASELINWDTH,
        )
        print(
            f'ASAS-SN range for {name}, g: {min(asdata["g"]["mag"])}-{max(asdata["g"]["mag"])}'
        )
    if "V" in const.ASASFILT and len(asdata["V"]["mag"]):
        plt.errorbar(
            asdata["V"]["hjd"] - JD_SHIFT, (asdata["V"]["mag"] - curveshift["V"]), asdata["V"]["magerr"],
            marker="s", ls="none", c=const.FILTER_COLORS["V"], markeredgecolor="k", zorder=-1,
            mew=const.MEW, label=f"ASAS-SN V{mk_fsh(curveshift['V'])}", ms=const.ASASMS-1,
            elinewidth=const.ASELINWDTH,
        )
        print(
            f'ASAS-SN range for {name}, V: {min(asdata["V"]["mag"])}-{max(asdata["V"]["mag"])}'
        )

if obj.get("gaiaobj"):
    gaiadata = pd.read_csv(f"{const.DATA_RAW}{obj['gaiafnam']}", skiprows=1)
    gaiadata = gaiadata[gaiadata["averagemag"].notna()]
    gaiadata = gaiadata[gaiadata["averagemag"] != "untrusted"]
    gaiadata["averagemag"] = gaiadata["averagemag"].astype(float)
    # Time of observation is in barycentric coordinate time (TCB)
    save_gaia_datafile(gaiadata, f"../data/{obj['gaiaobj'].removesuffix('.csv')}.dat")
    plt.plot(
        gaiadata["JD(TCB)"] - JD_SHIFT, gaiadata.get("averagemag") - filtshift["G"], "*",
        markeredgecolor="k", mew=const.MEW, ms=const.GAIAMS, c=const.FILTER_COLORS["G"],
        label=f"{obj['gaiaobj']}{mk_fsh(curveshift['G'])}",
        zorder=10,
    )
    if period:
        gaiadata = mk_phased(gaiadata, epoch, period, jdnam="JD(TCB)")

if obj.get("gaiadr3fnam"):
    gdata = pd.read_csv(
        f"../data/{obj.get('gaiadr3fnam')}",
        sep=r"\s+",
        names=["TimeG", "Gmag"]
        # names=["TimeG", "Gmag", "TimeBP", "BPmag", "TimeRP", "RPmag"]
    )
    gdata["TimeG"] += const.GAIADR3_JD_SHIFT

if obj.get("gaiadr3fnam"):
    plt.plot(
        gdata["TimeG"] - JD_SHIFT, gdata["Gmag"] - filtshift["G"],
        marker="d", ls="none",
        c=const.FILTER_COLORS["G"], label=f"Gaia G{mk_fsh(filtshift['G'])}",
        markeredgecolor="k", mew=const.GMEW, zorder=10,
        ms=7)

# Atlas Stuff
const.ATLAS_SAVEDATA = True
if obj.get("atlaslim"):
    const.ATLAS_MAG_LIM = obj.get("atlaslim")
Atl = Atlas(fnam, atlaslim=const.ATLAS_MAG_LIM)
if not const.ATLAS_RAWDATA and Atl.is_data_exist("o", lim=const.ATLAS_MAG_LIM["o"]):
    Atl.read_prepared_data(filtlims={"o": const.ATLAS_MAG_LIM["o"], "c": const.ATLAS_MAG_LIM["c"]})
    if args.verbose:
        print("read atlas data for", fnam, "with lim", const.ATLAS_MAG_LIM)
elif const.ATLAS_RAWDATA and Atl.is_raw_data_exist():
    Atl.read_raw_data()
    if obj.get("atlasuplim"):
        const.ATLAS_MAGLIM_UP = obj.get("atlasuplim")
    Atl.prepare_data(ra, dec, maglim_up=const.ATLAS_MAGLIM_UP, obs="Atlas-CHL")
    if args.verbose:
        print("read and prepare atlas raw data for", fnam, "with lim", const.ATLAS_MAG_LIM)
    if const.ATLAS_SAVEDATA:
        if args.verbose:
            print(f"save atlas data for {const.ATLASFILT} filters")
        for filtr in const.ATLASFILT:
            Atl.save_datafile(filtr)
else:
    if args.verbose:
        print("No Atlas data", hasattr(Atl, "data"))

if obj["plot"].get("atlasms"):
    const.ATLASMS = obj["plot"].get("atlasms")
if hasattr(Atl, "data") and len(Atl.data.index):
    for afiltr in const.ATLASFILT:
        alclr = const.FILTER_COLORS[afiltr]
        amaglim = None
        if obj.get("period"):
            Atl.mk_phased(epoch, period)
        alcurve = Atl.get_data(afiltr, maglim=amaglim)
        if obj.get("curveshift") and obj.get("clrshift") and const.ATLASFILT:
            plt.errorbar(
                alcurve["hjd"] - JD_SHIFT, alcurve["mag"] - curveshift[afiltr], alcurve["magerr"],
                marker="o", ls="none", elinewidth=const.ATELINWDTH[afiltr], c=alclr,
                ms=const.ATLASMS[afiltr], label=f"ATLAS {afiltr}{mk_fsh(curveshift[afiltr])}",
                zorder=0, markeredgecolor="k", mew=const.MEW,
            )
        elif const.ATLASFILT:
            plt.errorbar(
                alcurve["hjd"] - JD_SHIFT, alcurve["mag"], alcurve["magerr"], marker="o",
                ls="none", elinewidth=const.ATELINWDTH[afiltr], c=alclr, ms=const.ATLASMS[afiltr],
                label=f"ATLAS {afiltr}", zorder=0, markeredgecolor="k", mew=const.MEW,
            )

# -== ZTF stuff ==-
if obj["plot"].get("zms"):
    const.ZMS = obj["plot"].get("zms")
zmss = {"g": const.ZMS, "r": const.ZMS, "i": const.ZMS}

Zt = Ztf(fnam, ztflim=obj.get("ztflim"))
zdata = {}
for filtr in const.ZTFFILT:
    if Zt.is_raw_data_exist(add=f"{filtr}1", ext="csv"):
        data = Zt.read_raw_data_filtr(filtr=filtr)
        if obj.get("zcatf") and filtr in obj.get("zcatf"):
            if args.verbose:
                print(f"Apply zcat filter for ZTF {filtr} data")
            data = data[data["catflags"] != 32768]
        if obj.get("ztfran"):
            data = data[data["mag"] >= obj["ztfran"][1] + filtshift[filtr]]
            data = data[data["mag"] <= obj["ztfran"][0] + filtshift[filtr]]
        try:
            save_datafile(data, f"{const.DATA_PATH}{fnam}-ztf{filtr}.dat", hjdprec=7, magprec=6)
            try:
                print(
                    f'ZTF range for {name} {filtr}: {round(min(data["mag"]), 2)}-{round(max(data["mag"]), 2)}'
                )
            except ValueError:
                pass
            zdata[filtr] = data
            if period:
                zdata[filtr] = mk_phased(data, epoch, period)
            if obj.get("curveshift") and obj.get("clrshift"):
                plt.errorbar(
                    zdata[filtr]["hjd"] - JD_SHIFT, zdata[filtr]["mag"] - curveshift[filtr], zdata[filtr]["magerr"],
                    marker="o", ls="none", c=const.FILTER_COLORS[f"z{filtr}"],
                    elinewidth=const.ZELINWDTH-0.3, label=f"ZTF {filtr}{mk_fsh(curveshift[filtr])}",
                    ms=zmss[filtr], markeredgecolor="k", mew=const.MEW,
                )
            else:
                plt.errorbar(
                    zdata[filtr]["hjd"] - JD_SHIFT, zdata[filtr]["mag"], zdata[filtr]["magerr"],
                    marker="o", ls="none", c=const.FILTER_COLORS[f"z{filtr}"],
                    elinewidth=const.ZELINWDTH, label=f"ZTF {filtr}",
                    ms=zmss[filtr], markeredgecolor="k", mew=const.MEW,
                )
        except TypeError:
            continue
        except NameError:
            print(f"Err {filtr}")

if obj.get("ztfobj") and obj["plot"].get("zobj"):
    Z = Zobj(obj["ztfobj"])
    Z.read_raw_data(add="")
    Z.prepare_data()
    if period:
        Z.mk_phased(epoch, period)
    for filt in "gri":
        Zdata_filt = Z.data[Z.data["filter"] == filt]
        if len(Zdata_filt.index):
            plt.errorbar(
                Zdata_filt["hjd"] - JD_SHIFT, Zdata_filt["mag"] - curveshift[f"{filt}"], Zdata_filt["magerr"],
                marker="o", ls="none", c=const.FILTER_COLORS[f'z{filt}'],
                elinewidth=const.ZELINWDTH-0.3, ms=zmss[filt], markeredgecolor="k", mew=const.MEW,
            )

# if obj.get("ztfobj") and obj["plot"].get("zobj"):
#     Al = Alerce(obj["ztfobj"])
#     Al.read_raw_data(add="")
#     Al.prepare_data()
#     if period:
#         Al.mk_phased(epoch, period)
#     for filt in "gri":
#         Aldata_filt = Al.data[Al.data["filter"] == filt]
#         if len(Aldata_filt.index):
#             plt.errorbar(
#                 Aldata_filt["hjd"] - JD_SHIFT, Aldata_filt["mag"] - curveshift[f"ALeRCE {filt}"], Aldata_filt["magerr"],
#                 marker="o", ls="none", c=const.FILTER_COLORS[f'z{filt}'],
#                 elinewidth=const.ZELINWDTH-0.3, ms=zmss[filt], markeredgecolor="k", mew=const.MEW,
#             )

if const.GDSFILT:
    for filt in const.GDSFILT:
        gdsdata = read_gds_data(fnam, filt, ra, dec)
        plt.errorbar(
            gdsdata["hjd"] - JD_SHIFT, gdsdata["mag"] - filtshift[f"g{filt}"], gdsdata["magerr"],
            marker="o", ls="none", c=const.FILTER_COLORS[f"g{filt}"],
            label=f"GDS {filt}{mk_fsh(filtshift[f'g{filt}'])}", ms=const.GDSMS,
            markeredgecolor="k", mew=const.MEW, elinewidth=const.GDSELINWDTH,
        )


if obj.get("oglefnam"):
    Ogl = Ogle(obj.get("oglefnam"))
    Ogl.read_raw_data()
    plt.errorbar(
        Ogl.data["hjd"] - JD_SHIFT, Ogl.data["mag"] - filtshift["I"], Ogl.data["magerr"],
        marker="o", ls="none", c=const.FILTER_COLORS["I"],
        label=f"OGLE I{mk_fsh(filtshift['I'])}", ms=const.OGLEMS,
        markeredgecolor="k", mew=const.MEW, elinewidth=const.OELINWDTH,
    )

    if period:
        Ogl.mk_phased(epoch, period)

    imerg = pd.concat((zdata["i"], Ogl.data)).sort_values(by=["hjd"])
    pd.cfDataFrame({
        "hjd": imerg["hjd"].astype("str").str.ljust(18, "0"),
        "mag": imerg["mag"].round(6).astype("str").str.ljust(9, "0"),
        "magerr": imerg["magerr"].round(6).astype("str").str.ljust(8, "0"),
    }).to_csv(f"{const.DATA_PATH}{fnam}-ogle-ztf-imerg.dat", sep=" ", index=False)
    print(f"OGLE range for {name}: {min(Ogl.data['mag'])}-{max(Ogl.data['mag'])}")

# Plot PS1 data
psnam = ""
psdata = {}
if obj.get("pslim"):
    const.PS_MAG_LIM = obj.get("pslim")
if not args.localps and obj["plot"].get("psfilt") and Psdat.results:
    psnam = "-ps1"
    for i, filter in enumerate(const.PSFILT):
        w = np.where(Psdat.dtab["filter"] == filter)
        try:
            # if not obj.get("period"):
            #     mk_lc_ps1(t_corr, mag[w], magerr[w], f"{DATA_PATH}{objpsfnam}{filter}.dat", jds="hjd")
            print(
                f"PS1 for {name}, filter {filter}: {str(round(min(Psdat.mag[w]), 2)).ljust(5, '0')}-{str(round(max(Psdat.mag[w]), 2)).ljust(5, '0')}"
            )
            # plt.plot(t[w], mag[w], f'{line}{MARKERS[i]}{colors[i]}', markeredgewidth=0.5, markeredgecolor='k',
            #          label=f"PS1 {filter}", lw=0.75)
            # psms = obj.get("plot").get("psms") if obj.get("plot") and obj.get("plot").get("psms") else 5
            psdata[filter] = Psdat.mk_data(w)
            save_datafile(
                psdata[filter], f"{const.DATA_PATH}{objpsfnam}{filter}.dat",
                hjdprec=7, magprec=4,
            )
            psdata[filter] = psdata[filter][psdata[filter]["magerr"] < const.PS_MAG_LIM]

            if obj.get("curveshift") and obj.get("clrshift"):
                plt.errorbar(
                    psdata[filter]["hjd"] - JD_SHIFT, psdata[filter]["mag"] - curveshift[f"ps{filter}"],
                    psdata[filter]["magerr"], marker=const.MARKERS[i], ls="none",
                    c=const.FILTER_COLORS[f"ps_{filter}"], markeredgewidth=0.5,
                    markeredgecolor="k", ms=MSS[i], lw=0.75,
                    label=f"PS1 {filter}{mk_fsh(curveshift[f'ps{filter}'])}",
                )
            else:
                plt.errorbar(
                    psdata[filter]["hjd"] - JD_SHIFT, psdata[filter]["mag"], psdata[filter]["magerr"],
                    marker=const.MARKERS[i], ls="none", c=const.FILTER_COLORS[f"ps_{filter}"],
                    markeredgewidth=0.5, markeredgecolor="k", label=f"PS1 {filter}", lw=0.75,
                    ms=MSS[i],
                )
            if period:
                # mk_lc_ps1(t_corr, mag[w], magerr[w], f"{DATA_PATH}{objpsfnam}{filter}.dat", jds="hjd")
                psdata[filter] = mk_phased(
                    psdata[filter], epoch, period
                )
        except (IndexError, ValueError):
            print(f"no PS1 plot for {filter}")
            # psdata[filter] = None
if args.localps and obj["plot"].get("psfilt"):
    psnam = "-ps1"
    for i, filter in enumerate(const.PSFILT):
        try:
            ps_fnam = f"{const.DATA_PATH}{fnam}-ps{filter}.dat"
            pslc = read_ps_data(ps_fnam)
            print(
                f"PS1 for {name}, filter {filter}:",
                f"{str(round(min(pslc['mag']), 2)).ljust(5, '0')}-{str(round(max(pslc['mag']), 2)).ljust(5, '0')}"
            )
            psfiltshift = filtshift[filter]
            if obj.get("clrshift").get("ps" + filter):
                psfiltshift = obj.get("clrshift").get("ps" + filter)
            plt.errorbar(
                pslc["hjd"] - JD_SHIFT, pslc["mag"] - psfiltshift, pslc["magerr"],
                marker=const.MARKERS[i], ls="none", c=const.FILTER_COLORS[f"ps_{filter}"],
                mew=const.MEW, markeredgecolor="k",
                label=f"PS1 {filter}{mk_fsh(psfiltshift)}", lw=0.75,
                ms=MSS[i],
            )
            if period:
                psdata[filter] = pd.DataFrame({
                    "hjd": pslc["hjd"],
                    "mag": pslc["mag"],
                    "magerr": pslc["magerr"],
                })
                psdata[filter] = psdata[filter][psdata[filter]["magerr"] < const.PS_MAG_LIM]
                psdata[filter] = mk_phased(
                    psdata[filter], epoch, period
                )
        except (IndexError, ValueError, FileNotFoundError):
            print(f"no PS1 plot for {filter}")

if obj.get("crtsfnam"):
    Cr = Crts(fnam, maglim=obj.get("crtslim"))
    Cr.read_raw_data()
    Cr.prepare_data(ra, dec)
    if period:
        Cr.mk_phased(epoch, period)
    plt.errorbar(
        Cr.data["hjd"] - JD_SHIFT, Cr.data["mag"] - curveshift["CV"], Cr.data["magerr"],
        marker="o", ls="none", c=const.FILTER_COLORS["crts"],
        elinewidth=const.CELINWDTH, mew=const.CMEW, ms=const.CRTSMS,
        label=f"CRTS{mk_fsh(curveshift['CV'])}", markeredgecolor="k",
    )

plt.title(objtitle, fontsize=const.TITLEFNTSZ)
if obj["plot"].get("xmal"):
    const.XMAL = obj["plot"]["xmal"]
    ax.xaxis.set_major_locator(MultipleLocator(const.XMAL))
if obj["plot"].get("xmil"):
    const.XMIL = obj["plot"]["xmil"]
    ax.xaxis.set_minor_locator(MultipleLocator(const.XMIL))
if obj["plot"].get("ymal"):
    const.YMAL = obj["plot"].get("ymal")
    ax.yaxis.set_major_locator(MultipleLocator(const.YMAL))
if obj["plot"].get("ymil"):
    const.YMIL = obj["plot"].get("ymil")
    ax.yaxis.set_minor_locator(MultipleLocator(const.YMIL))

if obj["plot"].get("xlima"):
    xmin, xmax = obj["plot"].get("xlima")
    plt.xlim(xmin, xmax)
    if args.lines and obj.get("max"):
        plt.plot([xmin, xmax], [obj["max"], obj["max"]],
                 "--k", lw=0.95, zorder=-20)
    if args.lines and epoch and obj["plot"].get("ylima"):
        plt.plot([epoch-JD_SHIFT, epoch-JD_SHIFT], obj["plot"]["ylima"],
                 "--k", lw=0.95, zorder=-20)
elif obj["plot"].get("xedges"):
    xmin, xmax = ax.get_xlim()
    xmin, xmax = xmin + obj["plot"]["xedges"], xmax - obj["plot"]["xedges"]
    plt.xlim(xmin, xmax)
    if args.lines and obj.get("max"):
        plt.plot([xmin, xmax], [obj["max"], obj["max"]],
                 "--k", lw=0.95, zorder=-20)
    if args.lines and epoch and obj["plot"].get("ylima"):
        plt.plot([epoch-JD_SHIFT, epoch-JD_SHIFT], obj["plot"]["ylima"],
                 "--k", lw=0.95, zorder=-20)

plt.ylabel("Mag", fontsize=16)
plt.xlabel(f"HJD - {JD_SHIFT}", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
if obj["plot"].get("leg"):
    plt.legend(fontsize=const.LEGFNTSZ, loc=obj["plot"].get("leg"))
else:
    plt.legend(fontsize=const.LEGFNTSZ)
if obj["plot"].get("ylima"):
    ax.set_ylim(obj["plot"].get("ylima"))
else:
    ax.invert_yaxis()

ztfnam = "-ztf" if const.ZTFFILT else ""
crtsnam = "-css" if obj.get("crtsfnam") else ""
atlasfnam = "-atlas" if hasattr(Atl, "data") and const.ATLASFILT else ""
asasfnam = "-asas" if obj.get("asasfnam") and const.ASASFILT else ""
oglenam = "-ogle" if obj.get("oglefnam") else ""
gdsnam = "-gds" if const.GDSFILT else ""
gaianam = "-gaia" if (obj.get("gaiaobj") or obj.get("gaiadr3fnam")) else ""
fname = f"{fnam}{psnam}{ztfnam}{crtsnam}{atlasfnam}{asasfnam}{oglenam}{gdsnam}{gaianam}"
if args.show:
    fig.set_dpi(const.DPI_SCREEN)
    plt.get_current_fig_manager().window.state('zoomed')
    plt.show()
else:
    plt.savefig(
        f"../lc/{args.savedir}{fname}.{const.EXT}",
        dpi=const.DPI,
    )
# End of lightcurve plot
plt.close()

if period:
    # Phased plot figure
    fig2, ax2 = plt.subplots(figsize=const.FIGSZ)
    fig2.subplots_adjust(*const.PLT_EDGES)
    data_to_merge = []  # data in HJD to be merged

    if "g" in const.ASASFILT:
        plt.errorbar(
            asdata["g"]["phased"], asdata["g"]["mag"] - curveshift["asasg"], asdata["g"]["magerr"],
            marker="d", ls="none", c=const.FILTER_COLORS["asasg"], markeredgecolor="k",
            mew=const.MEW, label=f"ASAS-SN g{mk_fsh(curveshift['asasg'])}",
            ms=const.ASASMS, elinewidth=const.ASELINWDTH,
        )
        if const.DATA_TO_MERGE.get("ASAS") and "g" in const.DATA_TO_MERGE["ASAS"]:
            data_to_merge.append(
                pd.DataFrame({
                    "hjd": asdata["g"]["hjd"],
                    "mag": asdata["g"]["mag"] - curveshift["asasg"],
                    "magerr": asdata["g"]["magerr"],
                })
            )
    if "V" in const.ASASFILT:
        plt.errorbar(
            asdata["V"]["phased"], asdata["V"]["mag"] - curveshift["V"], asdata["V"]["magerr"],
            marker="s", ls="none", c=const.FILTER_COLORS["V"], markeredgecolor="k",
            mew=const.MEW, ms=const.ASASMS, elinewidth=const.ASELINWDTH,
            label=f"ASAS-SN V{mk_fsh(curveshift['V'])}",
        )

    if obj["plot"].get("atlaselw"):
        const.ATELINWDTH = obj["plot"].get("atlaselw")
    if hasattr(Atl, "data") and len(Atl.data.index) and const.ATLASFILT:
        Atl.mk_phased(epoch, period)
        for afiltr in const.ATLASFILT:
            alclr = const.FILTER_COLORS[afiltr]
            alcurve = Atl.get_data(afiltr, maglim=amaglim)
            plt.errorbar(
                alcurve["phased"], alcurve["mag"] - filtshift[afiltr], alcurve["magerr"],
                marker="o", ls="none", elinewidth=const.ATELINWDTHT, c=alclr,
                label=f"ATLAS {afiltr}{mk_fsh(filtshift[afiltr])}", zorder=0,
                ms=const.ATLASMS[afiltr], mew=const.MEW,
                markeredgecolor=obj["plot"].get("atledgclr"),
            )
            if const.DATA_TO_MERGE.get("ATLAS") and afiltr in const.DATA_TO_MERGE["ATLAS"]:
                data_to_merge.append(
                    pd.DataFrame({
                        "hjd": alcurve["hjd"],
                        "mag": alcurve["mag"] - filtshift[afiltr],
                        "magerr": alcurve["magerr"],
                    })
                )

    for filtr in const.ZTFFILT:
        if Zt.is_raw_data_exist(add=f"{filtr}1", ext="csv"):
            zsorted = zdata[filtr].sort_values(by="phased")
            plt.errorbar(
                zdata[filtr]["phased"], zdata[filtr]["mag"] - filtshift[filtr],
                zdata[filtr]["magerr"], mew=const.MEW, markeredgecolor="k",
                marker="o", ls="none", elinewidth=const.ZELINWDTH,
                c=const.FILTER_COLORS[f"z{filtr}"], ms=zmss[filtr],
                label=f"ZTF {filtr}{mk_fsh(filtshift[filtr])}",
            )
            if const.DATA_TO_MERGE.get("ZTF") and filtr in const.DATA_TO_MERGE["ZTF"]:
                data_to_merge.append(
                    pd.DataFrame({
                        "hjd": zdata[filtr]["hjd"],
                        "mag": zdata[filtr]["mag"] - filtshift[filtr],
                        "magerr": zdata[filtr]["magerr"],
                    })
                )
    if obj.get("ztfobj") and obj["plot"].get("zobj"):
        for filt in "gri":
            Zdata_filt = Z.data[Z.data["filter"] == filt]
            if len(Zdata_filt.index):
                plt.errorbar(
                    Zdata_filt["phased"], Zdata_filt["mag"] - curveshift[f"{filt}"],
                    Zdata_filt["magerr"], marker="o", ls="none", ms=zmss[filt],
                    c=const.FILTER_COLORS[f"z{filt}"], markeredgecolor="k",
                    label=f"ZTF {filt}{mk_fsh(curveshift[f'{filt}'])}",
                    mew=const.MEW, elinewidth=const.ZELINWDTH-0.3,
                )
            # if const.DATA_TO_MERGE.get("ALeRCE") and filt in const.DATA_TO_MERGE["ALeRCE"]:
            #     data_to_merge.append(
            #         pd.DataFrame({
            #             "hjd": Aldata_filt["hjd"],
            #             "mag": Aldata_filt["mag"] - filtshift[f"ALeRCE {filt}"],
            #             "magerr": Aldata_filt["magerr"],
            #         })
            #     )

    if const.GDSFILT:
        for filt in const.GDSFILT:
            gdsdata = read_gds_data(fnam, filt, ra, dec)
            gdsdata = mk_phased(gdsdata, epoch, period)
            plt.errorbar(
                gdsdata["phased"], gdsdata["mag"] - filtshift[f"g{filt}"], gdsdata["magerr"],
                marker="o", ls="none", c=const.FILTER_COLORS[f"g{filt}"],
                label=f"GDS {filt}{mk_fsh(filtshift[f'g{filt}'])}", ms=const.GDSMS,
                markeredgecolor="k", mew=const.MEW, elinewidth=const.GDSELINWDTH,
            )

    if obj.get("oglefnam"):
        plt.errorbar(
            Ogl.data["phased"], Ogl.data["mag"] - filtshift["I"], Ogl.data["magerr"],
            marker="o", ls="none", c=const.FILTER_COLORS["I"],
            label=f"OGLE I{mk_fsh(filtshift['I'])}", ms=const.OGLEMS,
            markeredgecolor="k", mew=const.MEW, elinewidth=const.OELINWDTH,
        )
        if const.DATA_TO_MERGE.get("OGLE"):
            data_to_merge.append(
                pd.DataFrame({
                    "hjd": Ogl.data["hjd"],
                    "mag": Ogl.data["mag"] - filtshift["I"],
                    "magerr": Ogl.data["magerr"],
                })
            )

    if obj.get("gaiaobj"):
        plt.plot(
            gaiadata["phased"], gaiadata.get("averagemag") - filtshift["G"], "*",
            markeredgecolor="k", mew=const.MEW, ms=const.GAIAMS, zorder=10,
            c=const.FILTER_COLORS["G"], label=f"{obj['gaiaobj']}{mk_fsh(filtshift['G'])}",
                )
        if const.DATA_TO_MERGE.get("Gaia"):
            data_to_merge.append(
                pd.DataFrame({
                    "hjd": gaiadata["JD(TCB)"],
                    "mag": gaiadata["averagemag"] - filtshift["G"],
                    "magerr": 0.1,
                })
            )
    if obj.get("gaiadr3fnam"):
        gdata = mk_phased(gdata, epoch, period, jdnam="TimeG")
        plt.plot(gdata["phased"], gdata["Gmag"] - filtshift["G"], marker="d", ls="none",
                 c=const.FILTER_COLORS["G"], label=f"Gaia G{mk_fsh(filtshift['G'])}",
                 markeredgecolor="k", mew=const.GMEW, zorder=30, ms=const.GAIAMS)

        data_to_merge.append(
            pd.DataFrame({
                "hjd": gdata["TimeG"],
                "mag": gdata["Gmag"] - filtshift["G"],
                "magerr": 0,
            })
        )

    if obj.get("crtsfnam"):
        data_to_merge.append(
            pd.DataFrame({
                "hjd": Cr.data["hjd"],
                "mag": Cr.data["mag"] - filtshift["CV"],
                "magerr": Cr.data["magerr"],
            })
        )
        plt.errorbar(
            Cr.data["phased"], Cr.data["mag"] - filtshift["CV"], Cr.data["magerr"],
            marker="o", ls="none", c=const.FILTER_COLORS["crts"],
            label=f"CRTS{mk_fsh(filtshift['CV'])}", ms=const.CRTSMS,
            markeredgecolor="k", mew=const.MEW, elinewidth=const.CELINWDTH,
        )

    if obj.get("dssdata"):
        poss = pd.DataFrame.from_dict({"hjd": [obj["dssdata"]["R"][0]], "mag": [obj["dssdata"]["R"][1]]})
        poss["hjd"] = mk_hjd_corr(poss["hjd"]-JD_SHIFT, ra, dec)
        poss = mk_phased(poss, epoch, period, jdnam="hjd")
        plt.plot(poss["phased"], poss["mag"], "v", c="darkred", markeredgecolor="k",
                 markeredgewidth=const.DMS, ms=const.DSSMS, label="POSS II $R_c$")

    if const.PSFILT:
        for i, filter in enumerate(const.PSFILT):
            if filter in psdata:
                # try:
                psfiltshift = filtshift["ps" + filter]
                if filtshift.get("ps" + filter):
                    psfiltshift = filtshift.get("ps" + filter)
                plt.errorbar(
                    psdata[filter]["phased"], psdata[filter]["mag"] - psfiltshift,
                    psdata[filter]["magerr"], marker=const.MARKERS[i], ls="none",
                    c=const.FILTER_COLORS[f"ps_{filter}"], mew=const.MEW,
                    elinewidth=const.PSELINWDTH, markeredgecolor="k", ms=MSS[i],
                    label=f"PS1 {filter}{mk_fsh(psfiltshift)}",
                )
                if const.DATA_TO_MERGE.get("PS1") and filter in const.DATA_TO_MERGE["PS1"]:
                    data_to_merge.append(
                        pd.DataFrame({
                            "hjd": psdata[filter]["hjd"],
                            "mag": psdata[filter]["mag"] - psfiltshift,
                            "magerr": psdata[filter]["magerr"],
                        })
                    )
                    if args.verbose:
                        print("save the PS1!", filter)
    try:
        merged = pd.concat(data_to_merge)  # .sort_values(by=["hjd"])
        save_merged(merged, fnam)
    except ValueError:
        pass

    if args.zoom and obj["plot"].get("zoom") and obj["plot"]["zoom"].get("xlim"):
        XLIMM = obj["plot"]["zoom"]["xlim"]
    elif obj["plot"].get("xlim"):
        XLIMM = obj["plot"]["xlim"]
    else:
        XLIMM = const.XLIMP

    if obj.get("min") and isinstance(obj["min"], str):
        obj["min"] = float(obj["min"].strip("<").strip(":"))
    if args.model and obj.get("minecl"):
        obj["min"] = obj["minecl"]
    if obj.get("min2") and isinstance(obj["min2"], str):
        obj["min2"] = float(obj["min2"].strip("<").strip(":"))
    if args.model:
        if all([obj.get(x) for x in ("max", "min", "d")]):
            datax = [
                XLIMM[0],
                -obj["d"]/2,
                0,
                obj["d"]/2,
            ]
            datay = [
                obj["max"],
                obj["max"],
                obj["min"],
                obj["max"],
            ]
        if all([obj.get(x) for x in ("max", "2ndmin", "min2", "d2")]):
            datax.extend((
                obj["2ndmin"]-obj["d2"]/2,
                obj["2ndmin"],
                obj["2ndmin"]+obj["d2"]/2,
                1-obj["d"]/2,
                1,
                1+obj["d"]/2,
                max(XLIMM[1], 1+obj["d"]/2),
            ))
            datay.extend((
                obj["max"],
                obj["min2"],
                obj["max"],
                obj["max"],
                obj["min"],
                obj["max"],
                obj["max"],
            ))
        elif all([obj.get(x) for x in ("max", "min", "d")]):
            datax.extend((
                1-obj["d"]/2,
                1,
                1+obj["d"]/2,
                max(XLIMM[1], 1+obj["d"]/2),
            ))
            datay.extend((
                obj["max"],
                obj["min"],
                obj["max"],
                obj["max"],
            ))
        try:
            plt.plot(datax, datay, "-ok", lw=3)
        except NameError:
            pass

    # Phased plot parameters
    if obj["plot"].get("leg"):
        leg = plt.legend(fontsize=const.LEGFNTSZ, loc=obj["plot"].get("leg"))
    else:
        leg = plt.legend(fontsize=const.LEGFNTSZ)
    leg.set_zorder(99)
    plt.title(objtitle, fontsize=const.TITLEFNTSZ)
    plt.ylabel("Mag", fontsize=16)
    plt.xlabel("Phase", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim(XLIMM)

    YLIM = obj["plot"].get("ylim")
    if args.zoom and obj["plot"].get("zoom") and obj["plot"]["zoom"].get("ylim"):
        YLIM = obj["plot"]["zoom"]["ylim"]
    if YLIM:
        ax2.set_ylim(YLIM)
        if args.lines:
            plt.plot([0, 0], YLIM, "--k", lw=const.LW_MIN, zorder=-20)
            if XLIMM[1] != 1:
                plt.plot([1, 1], YLIM, "--k", lw=const.LW_MIN, zorder=-20)
            if obj.get("2ndmin"):
                plt.plot([obj["2ndmin"], obj["2ndmin"]], YLIM, "--k", lw=0.95, zorder=-20)
            if obj.get("max"):
                plt.plot([XLIMM[0], XLIMM[1]], [obj["max"], obj["max"]], "--k", lw=0.95, zorder=-20)
            if obj.get("min"):
                plt.plot([XLIMM[0], XLIMM[1]], [obj["min"], obj["min"]], "--",
                         c="grey", lw=0.95, zorder=-20)
            if obj.get("min2"):
                plt.plot([XLIMM[0], XLIMM[1]], [obj["min2"], obj["min2"]], "--",
                         c="grey", lw=0.95, zorder=-20)
            if obj.get("d"):
                plt.plot([obj["d"]/2, obj["d"]/2], YLIM, "--b", lw=0.95, zorder=-10)
                plt.plot([-obj["d"]/2, -obj["d"]/2], YLIM, "--b", lw=0.95, zorder=-10)
                plt.plot([1+obj["d"]/2, 1+obj["d"]/2], YLIM, "--b", lw=0.95, zorder=-10)
                plt.plot([1-obj["d"]/2, 1-obj["d"]/2], YLIM, "--b", lw=0.95, zorder=-10)
            if obj.get("d2") and obj.get("2ndmin"):
                plt.plot([obj["2ndmin"]+obj["d2"]/2, obj["2ndmin"]+obj["d2"]/2], YLIM, "--",
                         c="violet", lw=0.95, zorder=-10)
                plt.plot([obj["2ndmin"]-obj["d2"]/2, obj["2ndmin"]-obj["d2"]/2], YLIM, "--",
                         c="violet", lw=0.95, zorder=-10)
    else:
        plt.gca().invert_yaxis()
    if obj["plot"].get("ymal"):
        ax2.yaxis.set_major_locator(MultipleLocator(const.YMAL))
    if args.zoom and obj["plot"].get("zoom") and obj["plot"]["zoom"].get("ymil"):
        ax2.yaxis.set_minor_locator(MultipleLocator(obj["plot"]["zoom"]["ymil"]))
    elif obj["plot"].get("ymil"):
        ax2.yaxis.set_minor_locator(MultipleLocator(const.YMIL))

    zmed = "z" if args.zoom else ""
    if args.zoom and obj["plot"].get("zoom") and obj["plot"]["zoom"].get("xmalp"):
        const.XMALP = obj["plot"]["zoom"]["xmalp"]
    elif obj["plot"].get("xmalp"):
        const.XMALP = obj["plot"].get("xmalp")
    if args.zoom and obj["plot"].get("zoom") and obj["plot"]["zoom"].get("xmilp"):
        const.XMILP = obj["plot"]["zoom"]["xmilp"]
    elif obj["plot"].get("xmilp"):
        const.XMILP = obj["plot"].get("xmilp")
    ax2.xaxis.set_major_locator(MultipleLocator(const.XMALP))
    ax2.xaxis.set_minor_locator(MultipleLocator(const.XMILP))
    ax2.tick_params(which="major", length=6)
    ax2.tick_params(which="minor", length=4)
    if args.show:
        fig2.set_dpi(const.DPI_SCREEN)
        plt.get_current_fig_manager().window.state('zoomed')
        plt.show()
    else:
        plt.savefig(
            f"../lc/{args.savedir}{fname}-ph_{round(period, 8)}{zmed}.{const.EXT}",
            dpi=const.DPI,
        )
