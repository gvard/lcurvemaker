import json
import os
import warnings
from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
from astropy.utils.exceptions import AstropyWarning

from curves import (
    Ps,
    Atlas,
    Ogle,
    read_ps_data,
    read_crts_data,
    mk_phased,
    mk_hjd_corr,
    mk_fsh,
    save_datafile,
    save_merged,
    save_gaia_datafile,
)


warnings.simplefilter("ignore", category=AstropyWarning)
pd.options.mode.chained_assignment = None  # default='warn'

parser = ArgumentParser(
    description="Graph plotter script with number of people in space"
)
parser.add_argument("-l", "--localps", action="store_true", help="local PS1 files")
parser.add_argument("-v", "--verln", action="store_true",
                    help="draw vertical line on phased plot")
parser.add_argument("-s", "--show", action="store_true",
                    help="Show interactive plot instead of saving figure")
parser.add_argument("-f", "--filename", type=str, default="",
                    help="file with object parameters")
args = parser.parse_args()

print("localps", args.localps, "filename", args.filename)
FN = "minkovskiy16.json"
if args.filename:
    FN = f"{args.filename}.json"
with open(f"../objects/{FN}", "r", encoding="utf8") as loc_file:
    objs = json.load(loc_file)

DATA_PATH = "../data/"
DATA_RAW = f"{DATA_PATH}raw/"
EXT = "png"
JD_SHIFT = 2400000.5

# Data filtering constants, change it in objects JSON parameters file
PS_MAG_LIM = 0.11
ASAS_MAG_LIM = 0.75
ATLAS_MAG_LIM = 0.19

# Plotting parameters
ELINEWDTH = 0.5
MARKERS = "DpshX"
# colors = "grmkb"
FILT_CLRS = {"g": "g", "r": "r", "i": "m", "z": "k", "y": "darkgrey"}

for name, obj in objs.items():
    fnam = name.lower().replace(" ", "")
    objpsfnam = f"{fnam}-ps"
    name_add = obj["other"]
    objname = f"{name} = {name_add}"
    ra, dec = obj["coordeg"]
    filtshift = {"g": 0, "r": 0, "psr": 0, "i": 0, "I": 0, "psi": 0, "z": 0, "o": 0, "c": 0, "V": 0}
    if obj.get("clrshift"):
        filtshift = obj.get("clrshift")
    print(f"-=={name}==-")
    psms = (
        obj.get("plot").get("psms")
        if obj.get("plot") and obj.get("plot").get("psms")
        else 5
    )
    MSS = [psms - 0.5, psms + 1, psms, psms + 1, psms]
    ZMS = 3  # 5
    if obj.get("curveshift") and obj.get("clrshift"):
        curveshift = obj.get("clrshift")
    # Get PS1 data
    if (
        not args.localps
        and obj.get("plot")
        and obj["plot"].get("psfilt")
        and not obj.get("pslocal")
    ):
        Psdat = Ps(ra, dec)
        Psdat.cone_search()
        if obj.get("plot") and obj.get("plot").get("psfilt"):
            Psdat.PSFILTRS = obj["plot"]["psfilt"]
        if Psdat.results:
            Psdat.get_table()
        else:
            print(Psdat.results, f"No PS1 data for {name}")

    # Light curve figure
    fig, ax = plt.subplots(figsize=(16, 9))
    fig.subplots_adjust(0.06, 0.09, 0.985, 0.95)

    # ASAS Stuff
    if obj.get("asasfnam"):
        asasd = pd.read_csv(f"{DATA_RAW}{obj['asasfnam']}")
        asasd = asasd.rename(
            columns={"HJD": "hjd", "mag_err": "magerr", "Filter": "filter"}
        )
        # For plotting upper limits
        asasd_bad = asasd[asasd["magerr"] == 99.99]
        try:
            asasd_bad = asasd_bad[asasd_bad["mag"] < 21]
        except TypeError:
            asasd_bad["mag"] = asasd_bad["mag"].str.lstrip(">").astype(float)
            asasd_bad = asasd_bad[asasd_bad["mag"] < 21]
        ASASERRLIM = obj.get("asaslim") if obj.get("asaslim") else ASAS_MAG_LIM
        asasd = asasd[asasd["magerr"] != 99.99]
        asasd = asasd[asasd["magerr"] < ASASERRLIM]
        asasd_V = asasd[asasd["filter"] == "V"]
        asasd_V["mag"] = asasd_V["mag"].astype(float)
        if obj.get("period"):
            mk_phased(asasd_V, obj["epoch"], obj.get("period"), jdnam="hjd")

        asasd_g = asasd[asasd["filter"] == "g"]
        asasd_g["mag"] = asasd_g["mag"].astype(float)
        if obj.get("period"):
            mk_phased(asasd_g, obj["epoch"], obj.get("period"), jdnam="hjd")
        # Strip upper limits
        asasd["mag"] = asasd["mag"].astype("str").str.lstrip(">").astype(float)
        try:
            print(
                f"ASAS range for {name}: {min(asasd_V['mag'])}-{max(asasd_V['mag'])} V,",
                f"{min(asasd_g['mag'])}-{max(asasd_g['mag'])} g",
            )
        except ValueError:
            pass
        save_datafile(asasd, f"{DATA_PATH}{fnam}-asassn.dat", hjdprec=6, magprec=3)
        pd.DataFrame({
            "hjd": asasd["hjd"].astype("str").str.ljust(13, "0"),
            "mag": asasd["mag"].round(6).astype("str").str.ljust(6, "0"),
            "magerr": asasd["magerr"].astype("str").str.ljust(5, "0"),
            "Filter": asasd["filter"],
        }).to_csv(f"{DATA_PATH}{fnam}-asassn.dat", sep=" ", index=False)
        if obj["plot"].get("asasupperlim") and len(asasd_bad["mag"]):
            plt.plot(asasd_bad["hjd"] - JD_SHIFT, asasd_bad["mag"], marker="v", ls="none", c="#bbb", ms=2, zorder=-2)  # , label="ASAS-SN"
        if (
            obj["plot"].get("asasfilt")
            and "V" in obj["plot"].get("asasfilt")
            and len(asasd_V["mag"])
        ):
            plt.errorbar(
                asasd_V["hjd"] - JD_SHIFT, asasd_V["mag"], asasd_V["magerr"],
                marker="s", ls="none", c="g", markeredgecolor="k", mew=0.5,
                label="ASAS-SN V", ms=4, elinewidth=ELINEWDTH,
            )
            print(
                f'ASASSN range for {name}, V: {min(asasd_V["mag"])}-{max(asasd_V["mag"])}'
            )
        if "g" in obj["plot"].get("asasfilt") and len(asasd_g["mag"]):
            plt.errorbar(
                asasd_g["hjd"] - JD_SHIFT, asasd_g["mag"], asasd_g["magerr"],
                marker="d", ls="none", c="lightgreen", markeredgecolor="k",  # lime
                mew=0.5, label="ASAS-SN g", ms=4, elinewidth=ELINEWDTH,
            )
            print(
                f'ASASSN range for {name}, g: {min(asasd_g["mag"])}-{max(asasd_g["mag"])}'
            )

    if obj.get("gaiafnam"):
        gaiadata = pd.read_csv(f"{DATA_RAW}{obj['gaiafnam']}", skiprows=1)
        gaiadata = gaiadata[gaiadata["averagemag"].notna()]
        gaiadata = gaiadata[gaiadata["averagemag"] != "untrusted"]
        gaiadata["averagemag"] = gaiadata["averagemag"].astype(float)
        # Time of observation is in barycentric coordinate time (TCB)
        gaiadata["JD(TCB)"] = gaiadata["JD(TCB)"]  # - JD_SHIFT
        gaiadata["hmjd"] = gaiadata["JD(TCB)"] - JD_SHIFT
        save_gaia_datafile(gaiadata, f"../data/{obj['gaiafnam'].removesuffix('.csv')}.dat")
        plt.plot(
            gaiadata["hmjd"], gaiadata.get("averagemag"), "*", markeredgecolor="k",
            mew=0.5, ms=6, c="darkgrey", label=obj["gaiaobj"],
        )
        if obj.get("period"):
            gaiadata = mk_phased(gaiadata, obj["epoch"], obj.get("period"), jdnam="JD(TCB)")

    # Atlas Stuff
    ATLAS_RAWDATA = True
    ATLAS_SAVEDATA = True
    ATLAS_CLRS = {"o": "darkorange", "c": "darkblue"}
    if obj.get("atlaslim"):
        ATLAS_MAG_LIM = obj.get("atlaslim")
    atlasfilt = "oc"
    if obj["plot"].get("atlasfilt"):
        atlasfilt = obj["plot"].get("atlasfilt")
    Atl = Atlas(fnam, atlaslim=ATLAS_MAG_LIM)
    if not ATLAS_RAWDATA and Atl.is_data_exist("o", lim=ATLAS_MAG_LIM):
        Atl.read_prepared_data(filtlims={"o": ATLAS_MAG_LIM, "c": ATLAS_MAG_LIM})
        print("read atlas data for", fnam, "with lim", ATLAS_MAG_LIM)
    elif ATLAS_RAWDATA and Atl.is_raw_data_exist():
        Atl.read_raw_data()
        Atl.prepare_data(ra, dec)
        print("read and prepare atlas raw data for", fnam, "with lim", ATLAS_MAG_LIM)
        if ATLAS_SAVEDATA:
            print(f"save atlas data for {atlasfilt} filters")
            for filtr in atlasfilt:
                Atl.save_datafile(filtr)
    else:
        print("No Atlas data", hasattr(Atl, "data"))

    atlasms = 3
    if obj["plot"].get("atlasms"):
        atlasms = obj["plot"].get("atlasms")
    if hasattr(Atl, "data") and len(Atl.data.index):
        for afiltr in atlasfilt:
            alclr = ATLAS_CLRS[afiltr]
            alms = {"o": atlasms, "c": atlasms - 1}[afiltr]
            alcurve = Atl.get_data(afiltr)
            if obj.get("curveshift") and obj.get("clrshift"):
                plt.errorbar(
                    alcurve["hjd"] - JD_SHIFT, alcurve["mag"] - curveshift[afiltr], alcurve["magerr"],
                    marker="o", ls="none", elinewidth=ELINEWDTH, c=alclr,
                    ms=alms, label=f"ATLAS {afiltr}{mk_fsh(curveshift[afiltr])}", zorder=0, markeredgecolor="k", mew=0.5,
                )
            else:
                plt.errorbar(
                    alcurve["hjd"] - JD_SHIFT, alcurve["mag"], alcurve["magerr"],
                    marker="o", ls="none", elinewidth=ELINEWDTH, c=alclr,
                    ms=alms, label=f"ATLAS {afiltr}", zorder=0, markeredgecolor="k", mew=0.5,
                )

    if obj.get("plot") and obj["plot"].get("ylim"):
        plt.ylim(obj["plot"].get("ylim"))
    if obj.get("period"):
        obj["period"] = float(obj["period"])
    if obj.get("epoch"):
        obj["epoch"] = float(obj["epoch"])

    # -== ZTF stuff ==-

    all_datar = []
    for i in range(1, 4):
        ztf_fnam = f"{DATA_RAW}{fnam}-zr{i}.csv"
        if os.path.isfile(ztf_fnam):
            all_datar.append(pd.read_csv(ztf_fnam))
    if all_datar:
        datar = pd.concat(all_datar)
        datar = datar.sort_values(by=["hjd"])
        datar = datar.reset_index(drop=True)
        save_datafile(datar, f"{DATA_PATH}{fnam}-ztfr.dat", hjdprec=7, magprec=6)
        if obj.get("period"):
            datar = mk_phased(datar, obj["epoch"], obj.get("period"))
    all_datag = []
    for i in range(1, 4):
        ztf_fnam = f"{DATA_RAW}{fnam}-zg{i}.csv"
        if os.path.isfile(ztf_fnam):
            all_datag.append(pd.read_csv(ztf_fnam))
    if all_datag:
        datag = pd.concat(all_datag)
        datag = datag.sort_values(by=["hjd"])
        datag = datag.reset_index(drop=True)
        save_datafile(datag, f"{DATA_PATH}{fnam}-ztfg.dat", hjdprec=7, magprec=6)
        if obj.get("period"):
            datag = mk_phased(datag, obj["epoch"], obj.get("period"))
    if (
        obj.get("plot")
        and obj["plot"].get("ztffilt")
        and "i" in obj["plot"].get("ztffilt")
    ):
        all_datai = []
        for i in range(1, 4):
            ztf_fnam = f"{DATA_RAW}{fnam}-zi{i}.csv"
            if os.path.isfile(ztf_fnam):
                all_datai.append(pd.read_csv(ztf_fnam))
        if all_datai:
            datai = pd.concat(all_datai)
            datai = datai.sort_values(by=["hjd"])
            datai = datai.reset_index(drop=True)
            save_datafile(
                datai, f"{DATA_PATH}{fnam}-ztfi.dat", hjdprec=7, magprec=6
            )
            print(
                f'ZTF range for {name}: {round(min(datai["mag"]), 2)}-{round(max(datai["mag"]), 2)} i'
            )
            if obj.get("period"):
                datai = mk_phased(datai, obj["epoch"], obj.get("period"))
    if all_datag and all_datar:
        xes = [
            min(min(datar["hjd"]), min(datag["hjd"])),
            max(max(datar["hjd"]), max(datag["hjd"])),
        ]
        ykas = (
            (round(min(datag["mag"]), 2), round(max(datag["mag"]), 2)),
            (round(min(datar["mag"]), 2), round(max(datar["mag"]), 2)),
        )
        print(
            f"ZTF range for {name}: {ykas[0][0]}-{ykas[0][1]} g, {ykas[1][0]}-{ykas[1][1]} r"
        )

    if obj.get("plot") and obj["plot"].get("zms"):
        ZMS = obj["plot"].get("zms")
    # if obj.get("plot") and obj["plot"].get("ztffilt") and "g" in obj["plot"].get("ztffilt"):
    ELINEWDTH = 0.75
    try:
        if obj.get("curveshift") and obj.get("clrshift"):
            plt.errorbar(
                datag["hjd"] - JD_SHIFT, datag["mag"] - curveshift['g'], datag["magerr"],
                marker="o", ls="none", c="g", elinewidth=ELINEWDTH, label=f"ZTF g{mk_fsh(curveshift['g'])}",
                ms=ZMS, markeredgecolor="k", mew=0.5,
            )
        else:
            plt.errorbar(
                datag["hjd"] - JD_SHIFT, datag["mag"], datag["magerr"],
                marker="o", ls="none", c="g", elinewidth=ELINEWDTH, label="ZTF g",
                ms=ZMS, markeredgecolor="k", mew=0.5,
            )
    except NameError:
        pass
    # if obj.get("plot") and obj["plot"].get("ztffilt") and "r" in obj["plot"].get("ztffilt"):
    try:
        if obj.get("curveshift") and obj.get("clrshift"):
            plt.errorbar(
                datar["hjd"] - JD_SHIFT, datar["mag"] - curveshift['r'], datar["magerr"],
                marker="o", ls="none", c="r", elinewidth=ELINEWDTH, label=f"ZTF r{mk_fsh(curveshift['r'])}",
                ms=4, markeredgecolor="k", mew=0.5, zorder=1
            )
        else:
            plt.errorbar(
                datar["hjd"] - JD_SHIFT, datar["mag"], datar["magerr"],
                marker="o", ls="none", c="r", elinewidth=ELINEWDTH, label="ZTF r",
                ms=ZMS, markeredgecolor="k", mew=0.5,
            )
    except NameError:
        pass
    # if obj.get("plot") and obj["plot"].get("ztffilt") and "i" in obj["plot"].get("ztffilt"):
    try:
        plt.errorbar(
            datai["hjd"] - JD_SHIFT, datai["mag"], datai["magerr"], marker="o", ls="none",
            c="m", elinewidth=ELINEWDTH, label="ZTF i", markeredgecolor="k", mew=0.5,
            ms=3,
        )
    except NameError:
        pass
    ztfnam = "-ztf" if obj.get("plot") and obj["plot"].get("ztffilt") else ""

    if obj.get("oglefnam"):
        Ogl = Ogle(obj.get("oglefnam"))
        Ogl.read_raw_data()
        plt.errorbar(
            Ogl.data["hjd"] - JD_SHIFT, Ogl.data["mag"], Ogl.data["magerr"], marker="o",
            ls="none", c="indigo", label="OGLE I",
            ms=4, markeredgecolor="k", mew=0.5, elinewidth=0.8,
        )

        if obj.get("period"):
            Ogl.mk_phased(obj["epoch"], obj["period"])

        imerg = pd.concat((datai, Ogl.data)).sort_values(by=["hjd"])
        pd.DataFrame({
            "hjd": imerg["hjd"].astype("str").str.ljust(18, "0"),
            "mag": imerg["mag"].round(6).astype("str").str.ljust(9, "0"),
            "magerr": imerg["magerr"].round(6).astype("str").str.ljust(8, "0"),
        }).to_csv(f"{DATA_PATH}{fnam}-ogle-ztf-imerg.dat", sep=" ", index=False)
        print(f"OGLE range for {name}: {min(Ogl.data['mag'])}-{max(Ogl.data['mag'])}")

    # Plot PS1 data
    psnam = ""
    psdata = {}
    if obj.get("pslim"):
        PS_MAG_LIM = obj.get("pslim")
    if not args.localps and obj.get("plot") and obj.get("plot").get("psfilt") and Psdat.results:
        psnam = "-ps1"
        for i, filter in enumerate(obj["plot"]["psfilt"]):
            w = np.where(Psdat.dtab["filter"] == filter)
            try:
                print(
                    f"PS1 for {name}, filter {filter}: {str(round(min(Psdat.mag[w]), 2)).ljust(5, '0')}-{str(round(max(Psdat.mag[w]), 2)).ljust(5, '0')}"
                )
                plt.errorbar(
                    Psdat.t[w], Psdat.mag[w], Psdat.magerr[w],
                    marker=MARKERS[i], ls="none", c=FILT_CLRS[filter],
                    markeredgewidth=0.5, markeredgecolor="k", label=f"PS1 {filter}", lw=0.75,
                    ms=MSS[i],
                )
                psdata[filter] = Psdat.mk_data(w)
                save_datafile(
                    psdata[filter], f"{DATA_PATH}{objpsfnam}{filter}.dat",
                    hjdprec=7, magprec=4,
                )
                psdata[filter] = psdata[filter][psdata[filter]["magerr"] < PS_MAG_LIM]
                if obj.get("period"):
                    psdata[filter] = mk_phased(
                        psdata[filter], obj["epoch"], obj.get("period")
                    )
            except (IndexError, ValueError):
                print(f"no PS1 plot for {filter}")
    if args.localps and obj.get("plot") and obj.get("plot").get("psfilt"):
        psnam = "-ps1"
        for i, filter in enumerate(obj["plot"]["psfilt"]):
            try:
                ps_fnam = f"{DATA_PATH}{fnam}-ps{filter}.dat"
                pslc = read_ps_data(ps_fnam)
                print(
                    f"PS1 for {name}, filter {filter}:",
                    f"{str(round(min(pslc['mag']), 2)).ljust(5, '0')}-{str(round(max(pslc['mag']), 2)).ljust(5, '0')}"
                )
                psfiltshift = filtshift[filter]
                if obj.get("clrshift").get("ps" + filter):
                    psfiltshift = obj.get("clrshift").get("ps" + filter)
                plt.errorbar(
                    pslc["mjd"], pslc["mag"] - psfiltshift, pslc["magerr"],
                    marker=MARKERS[i], ls="none", c=FILT_CLRS[filter],
                    markeredgewidth=0.5, markeredgecolor="k",
                    label=f"PS1 {filter}{mk_fsh(psfiltshift)}", lw=0.75,
                    ms=MSS[i],
                )
                if obj.get("period"):
                    t_corr = mk_hjd_corr(pslc["mjd"], ra, dec, obs="Haleakala")
                    psdata[filter] = pd.DataFrame({
                        "hjd": t_corr,
                        "mjd": pslc["mjd"],
                        "mag": pslc["mag"],
                        "magerr": pslc["magerr"],
                    })
                    psdata[filter] = psdata[filter][psdata[filter]["magerr"] < PS_MAG_LIM]
                    psdata[filter] = mk_phased(
                        psdata[filter], obj["epoch"], obj.get("period")
                    )
            except (IndexError, ValueError, FileNotFoundError):
                print(f"no PS1 plot for {filter}")

    if obj.get("crtsfnam"):
        hjds, mags, magerrs = read_crts_data(f"{DATA_PATH}{obj['crtsfnam']}")
        hjds -= JD_SHIFT
        ELINEWDTH = 0.8
        plt.errorbar(
            hjds, mags, magerrs,
            marker="o", ls="none", c="k", elinewidth=ELINEWDTH, label="CRTS",
            ms=4,
        )
    # del hjds, mags

    plt.title(objname, fontsize=18)
    # ax.ticklabel_format(useOffset=False, style="plain")
    xmal = 1000
    if obj.get("plot"):
        if obj["plot"].get("xmal"):
            ax.xaxis.set_major_locator(MultipleLocator(obj["plot"]["xmal"]))
        if obj["plot"].get("xmil"):
            ax.xaxis.set_minor_locator(MultipleLocator(obj["plot"]["xmil"]))
        # ax.yaxis.set_major_locator(MultipleLocator(0.5))
        if obj["plot"].get("ymil"):
            ax.yaxis.set_minor_locator(MultipleLocator(obj["plot"]["ymil"]))
        if obj["plot"].get("xedges"):
            xmin, xmax = ax.get_xlim()
            plt.xlim(xmin + obj["plot"]["xedges"], xmax - obj["plot"]["xedges"])
    plt.ylabel("Mag", fontsize=16)
    plt.xlabel("HMJD", fontsize=16)
    if not obj.get("atlasfnam") and not obj.get("asasfnam"):
        plt.xlabel(f"HJD - {JD_SHIFT}", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.tight_layout()
    # plt.grid(axis="y", ls=":")
    plt.legend(fontsize=14)  # , loc="upper left"
    if obj["plot"].get("ylim"):
        plt.ylim(obj["plot"].get("ylim"))
    plt.gca().invert_yaxis()

    crtsnam = "-css" if obj.get("crtsfnam") else ""
    atlasfnam = "-atlas" if hasattr(Atl, "data") else ""
    asasfnam = "-asas" if obj.get("asasfnam") else ""
    oglenam = "-ogle" if obj.get("oglefnam") else ""
    gaianam = "-gaia" if obj.get("gaiafnam") else ""
    if args.show:
        plt.show()
    else:
        plt.savefig(
            f"../lc/{fnam}{psnam}{ztfnam}{crtsnam}{atlasfnam}{asasfnam}{oglenam}{gaianam}.{EXT}",
            dpi=120,
        )
    # End of lightcurve plot
    plt.close()

    if obj.get("period"):
        # Phased plot figure
        fig, ax = plt.subplots(figsize=(16, 9))
        fig.subplots_adjust(0.06, 0.09, 0.985, 0.95)
        data_to_merge = []  # data in HJD to be merged
        ELINEWDTH = 0.8
        if "g" in obj["plot"].get("ztffilt") and "datag" in globals():
            plt.errorbar(
                datag["phased"], datag["mag"] - filtshift["g"], datag["magerr"],
                markeredgewidth=0.4, markeredgecolor="k", marker="o", ls="none",
                elinewidth=ELINEWDTH, c="g", label=f"ZTF g{mk_fsh(filtshift['g'])}",
                ms=4,
            )
            data_to_merge.append(
                pd.DataFrame({
                    "hjd": datag["hjd"],
                    "mag": datag["mag"] - filtshift["g"],
                    "magerr": datag["magerr"],
                })
            )
        if "r" in obj["plot"].get("ztffilt") and "datar" in globals():
            plt.errorbar(
                datar["phased"], datar["mag"] - filtshift["r"], datar["magerr"],
                marker="o", markeredgewidth=0.4, markeredgecolor="k", ls="none",
                elinewidth=ELINEWDTH, c="r", label=f"ZTF r{mk_fsh(filtshift['r'])}",
                ms=4,
            )
            data_to_merge.append(
                pd.DataFrame({
                    "hjd": datar["hjd"],
                    "mag": datar["mag"] - filtshift["r"],
                    "magerr": datar["magerr"],
                })
            )
        if "i" in obj["plot"].get("ztffilt"):
            try:
                plt.errorbar(
                    datai["phased"], datai["mag"] - filtshift["i"], datai["magerr"],
                    marker="o", markeredgewidth=0.4, markeredgecolor="k",
                    elinewidth=ELINEWDTH, ls="none", c="m",
                    label=f"ZTF i{mk_fsh(filtshift['i'])}",
                    ms=4,
                )
                data_to_merge.append(
                    pd.DataFrame({
                        "hjd": datai["hjd"],
                        "mag": datai["mag"] - filtshift["i"],
                        "magerr": datai["magerr"],
                    })
                )
            except NameError:
                pass
        ATMS = 2
        if obj["plot"].get("atlasms"):
            ATMS = obj["plot"].get("atlasms")
        ELINEWDTH = 0.2  # 0.5
        if obj["plot"].get("atlaselw"):
            ELINEWDTH = obj["plot"].get("atlaselw")
        if hasattr(Atl, "data") and len(Atl.data.index) and atlasfilt:
            Atl.mk_phased(obj["epoch"], obj.get("period"))
            for afiltr in atlasfilt:
                alclr = ATLAS_CLRS[afiltr]
                alms = {"o": atlasms, "c": atlasms - 1}[afiltr]
                alcurve = Atl.get_data(afiltr)
                plt.errorbar(
                    alcurve["phased"], alcurve["mag"] - filtshift[afiltr], alcurve["magerr"],
                    marker="o", ls="none", elinewidth=ELINEWDTH, c=alclr,
                    label=f"ATLAS {afiltr}{mk_fsh(filtshift[afiltr])}", zorder=0,
                    ms=alms, markeredgewidth=0.4, markeredgecolor="k",
                )
                data_to_merge.append(
                    pd.DataFrame({
                        "hjd": alcurve["hjd"],
                        "mag": alcurve["mag"] - filtshift[afiltr],
                        "magerr": alcurve["magerr"],
                    })
                )
        if obj.get("asasfnam"):
            plt.errorbar(
                asasd_V["phased"], asasd_V["mag"], asasd_V["magerr"],
                marker="s", ls="none", c="g", markeredgecolor="k", mew=0.5,
                label="ASAS-SN V", ms=4, elinewidth=ELINEWDTH,
            )
        if obj["plot"].get("asasfilt") and "g" in obj["plot"].get("asasfilt") and len(asasd_g["mag"]):
            plt.errorbar(
                asasd_g["phased"], asasd_g["mag"], asasd_g["magerr"],
                marker="d", ls="none", c="lightgreen", markeredgecolor="k",  # lime
                mew=0.5, label="ASAS-SN g", ms=4, elinewidth=ELINEWDTH,
            )

        if obj.get("oglefnam"):
            plt.errorbar(
                Ogl.data["phased"], Ogl.data["mag"] - filtshift["I"], Ogl.data["magerr"],
                marker="o", ls="none", c="indigo", label=f"OGLE I{mk_fsh(filtshift['I'])}", ms=5,
                markeredgecolor="k", mew=0.5, elinewidth=0.8,
            )
        if obj.get("gaiafnam"):
            plt.plot(
                        gaiadata["phased"], gaiadata.get("averagemag") - filtshift["G"], "*", markeredgecolor="k",
                        mew=0.5, ms=6, c="darkgrey", label=f"{obj['gaiaobj']}{mk_fsh(filtshift['G'])}",
                    )
            # data_to_merge.append(
            #     pd.DataFrame({
            #         "mjd": gaiadata["mjd"],
            #         "mag": gaiadata["mag"] - psfiltshift,
            #         "magerr": gaiadata["magerr"],
            #     })
            # )

        if obj.get("plot") and obj["plot"].get("psfilt"):
            for i, filter in enumerate(obj["plot"]["psfilt"]):
                if filter in psdata:
                    try:
                        psfiltshift = filtshift[filter]
                        if filtshift.get("ps" + filter):
                            psfiltshift = filtshift.get("ps" + filter)
                        ELINEWDTH = 1
                        plt.errorbar(
                            psdata[filter]["phased"], psdata[filter]["mag"] - psfiltshift,
                            psdata[filter]["magerr"],
                            marker=MARKERS[i], ls="none", c=FILT_CLRS[filter],
                            markeredgewidth=0.5, elinewidth=ELINEWDTH,
                            markeredgecolor="k", label=f"PS1 {filter}{mk_fsh(psfiltshift)}",
                            ms=MSS[i],
                        )
                        # data_to_merge.append(
                        #     pd.DataFrame({
                        #         "mjd": psdata[filter]["mjd"],
                        #         "mag": psdata[filter]["mag"] - psfiltshift,
                        #         "magerr": psdata[filter]["magerr"],
                        #     })
                        # )
                    except KeyError:
                        print("pass ps filter", filter)
                        pass
        merged = pd.concat(data_to_merge)  # .sort_values(by=["mjd"])
        save_merged(merged, fnam)

        # Phased plot parameters
        plt.legend(fontsize=14, loc="lower left")  # loc="lower left"
        plt.title(objname, fontsize=18)
        plt.gca().invert_yaxis()
        if obj["plot"].get("ymil"):
            ax.yaxis.set_minor_locator(MultipleLocator(obj["plot"].get("ymil")))
        plt.ylabel("Mag", fontsize=16)
        plt.xlabel("Phase", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim(-0.5, 1)
        YLIM = obj["plot"].get("ylim")
        if YLIM:
            plt.ylim(YLIM)
            if args.verln:
                plt.plot([0, 0], YLIM, "--k")
        YML = 0.5
        if obj["plot"].get("yml"):
            YML = obj["plot"].get("yml")
        ax.yaxis.set_major_locator(MultipleLocator(YML))
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.tick_params(which="major", length=6)
        ax.tick_params(which="minor", length=4)
        plt.savefig(
            f"../lc/{fnam}{ztfnam}{atlasfnam}{asasfnam}{oglenam}{gaianam}-ph_{round(obj.get('period'), 6)}.{EXT}",
            dpi=120,
        )
    # if loop on different objects, don't forget to delete user-defined
    # variables to avoid plotting wrong data
    # del psnam, ztfnam, crtsnam, atlasfnam, asasfnam
    # del datag, datar, datai
    # del atlas_c, atlas_o
