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
    read_ps_data,
    read_crts_data,
    mk_phased,
    mk_hjd_corr,
    mk_fsh,
    save_datafile,
    save_merged,
)


warnings.simplefilter("ignore", category=AstropyWarning)
pd.options.mode.chained_assignment = None  # default='warn'

parser = ArgumentParser(
    description="Graph plotter script with number of people in space"
)
parser.add_argument("-l", "--localps", action="store_true", help="local PS1 files")
parser.add_argument("-v", "--verln", action="store_true",
                    help="draw vertical line on phased plot")
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
    filtshift = obj.get("clrshift")
    print(f"-=={name}==-")
    psms = (
        obj.get("plot").get("psms")
        if obj.get("plot") and obj.get("plot").get("psms")
        else 5
    )
    MSS = [psms - 0.5, psms + 1, psms, psms + 1, psms]
    ZMS = 3  # 5
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
        ASASERRLIM = obj.get("asaserrlim") if obj.get("asaserrlim") else ASAS_MAG_LIM
        asasd = asasd[asasd["magerr"] != 99.99]
        asasd = asasd[asasd["magerr"] < ASASERRLIM]
        asasd_V = asasd[asasd["filter"] == "V"]
        asasd_V["mag"] = asasd_V["mag"].astype(float)
        asasd_g = asasd[asasd["filter"] == "g"]
        asasd_g["mag"] = asasd_g["mag"].astype(float)
        # Strip upper limits
        asasd["mag"] = asasd["mag"].astype("str").str.lstrip(">").astype(float)
        print(
            f"ASAS range for {name}: {min(asasd_V['mag'])}-{max(asasd_V['mag'])} V,",
            f"{min(asasd_g['mag'])}-{max(asasd_g['mag'])} g",
        )
        save_datafile(asasd, f"{DATA_PATH}{fnam}-asassn.dat", hjdprec=6, magprec=3)
        pd.DataFrame({
            "hjd": asasd["hjd"].astype("str").str.ljust(13, "0"),
            "mag": asasd["mag"].round(6).astype("str").str.ljust(6, "0"),
            "magerr": asasd["magerr"].astype("str").str.ljust(5, "0"),
            "Filter": asasd["filter"],
        }).to_csv(f"{DATA_PATH}{fnam}-asassn.dat", sep=" ", index=False)
        # if len(asasd_bad["mag"]):
        #     plt.plot(asasd_bad["HJD"] - JD_SHIFT, asasd_bad["mag"], marker="v", ls="none", c="grey", ms=2)  # , label="ASAS-SN"
        if (
            obj["plot"].get("asasfilt")
            and "V" in obj["plot"].get("asasfilt")
            and len(asasd_V["mag"])
        ):
            plt.errorbar(
                asasd_V["HJD"] - JD_SHIFT, asasd_V["mag"], asasd_V["mag_err"],
                marker="s", ls="none", c="g", markeredgecolor="k", mew=0.5,
                label="ASAS-SN V", ms=4,
            )
            print(
                f'ASASSN range for {name}, V: {min(asasd_V["mag"])}-{max(asasd_V["mag"])}'
            )
        if "g" in obj["plot"].get("asasfilt") and len(asasd_g["mag"]):
            plt.errorbar(
                asasd_g["HJD"] - JD_SHIFT, asasd_g["mag"], asasd_g["mag_err"],
                marker="d", ls="none", c="orange", markeredgecolor="k",
                mew=0.5, label="ASAS-SN g", ms=4,
            )
            print(
                f'ASASSN range for {name}, g: {min(asasd_g["mag"])}-{max(asasd_g["mag"])}'
            )

    if obj.get("gaiafnam"):
        gaiadata = pd.read_csv(f"{DATA_RAW}{obj['gaiafnam']}", skiprows=1)
        gaiadata = gaiadata[gaiadata["averagemag"].notna()]
        gaiadata = gaiadata[gaiadata["averagemag"] != "untrusted"]
        gaiadata["averagemag"] = gaiadata["averagemag"].astype(float)
        gaiadata["JD(TCB)"] = gaiadata["JD(TCB)"] - JD_SHIFT
        plt.plot(
            gaiadata["JD(TCB)"], gaiadata.get("averagemag"),
            "*", ms=5, c="orange", label=obj["gaiaobj"],
        )

    # Atlas Stuff
    if obj.get("atlasfnam"):
        if obj.get("atlaslim"):
            ATLAS_MAG_LIM = obj.get("atlaslim")
        if obj["plot"].get("atlaslim"):
            ATLAS_MAG_LIM = obj["plot"].get("atlaslim")
        print("atlasfnam!", obj.get("atlasfnam"), "Mag lim:", ATLAS_MAG_LIM)
        alcurve = pd.read_csv(f"{DATA_RAW}{obj['atlasfnam']}", delim_whitespace=True)
        alcurve = alcurve.rename(
            columns={"###MJD": "mjd", "m": "mag", "dm": "magerr", "F": "filter"}
        )
        alcurve = alcurve[alcurve["magerr"] < ATLAS_MAG_LIM]
        alcurve = alcurve[alcurve["mag"] > 11]
        alcurve["hjd"] = mk_hjd_corr(
            alcurve["mjd"], ra, dec, obs="Haleakala"
        )
        afiltr = "o"
        alcurve_o = alcurve[alcurve["filter"] == afiltr]
        atlas_data_fnam = f"{DATA_PATH}{obj['atlasfnam'].removeprefix('.txt')}-{afiltr}-cleaned-{ATLAS_MAG_LIM}m.dat"
        save_datafile(alcurve_o, atlas_data_fnam, hjdprec=6, magprec=3)

        afiltr = "c"
        alcurve_c = alcurve[alcurve["filter"] == afiltr]
        atlas_data_fnam = f"{DATA_PATH}{obj['atlasfnam'].removeprefix('.txt')}-{afiltr}-cleaned-{ATLAS_MAG_LIM}m.dat"
        save_datafile(alcurve_c, atlas_data_fnam, hjdprec=6, magprec=3)
        # Plot filtered raw data
        if obj["plot"].get("atlasfilt") and "c" in obj["plot"].get("atlasfilt"):
            plt.errorbar(
                alcurve_c["mjd"], alcurve_c["mag"], alcurve_c["magerr"],
                marker="o", ls="none", elinewidth=ELINEWDTH, c="darkblue",
                ms=2, label="ATLAS c", zorder=0,
            )
        if obj["plot"].get("atlasfilt") and "o" in obj["plot"].get("atlasfilt"):
            plt.errorbar(
                alcurve_o["mjd"], alcurve_o["mag"], alcurve_o["magerr"],
                marker="o", ls="none", elinewidth=ELINEWDTH, c="darkorange",
                ms=3, label="ATLAS o", zorder=0,
            )

        if obj.get("plot") and obj["plot"].get("ylim"):
            plt.ylim(obj["plot"].get("ylim"))

    if obj.get("period"):
        obj["period"] = float(obj["period"])
    if obj.get("epoch"):
        obj["epoch"] = float(obj["epoch"])

    # -== ZTF stuff ==-

    if obj.get("rawdata"):
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
            obj.get("rawdata")
            and obj.get("plot")
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
        plt.errorbar(
            datag["hjd"] - JD_SHIFT, datag["mag"], datag["magerr"],
            marker="o", ls="none", c="g", elinewidth=ELINEWDTH, label="ZTF g",
            ms=ZMS,
        )
    except NameError:
        pass
    # if obj.get("plot") and obj["plot"].get("ztffilt") and "r" in obj["plot"].get("ztffilt"):
    try:
        plt.errorbar(
            datar["hjd"] - JD_SHIFT, datar["mag"], datar["magerr"],
            marker="o", ls="none", c="r", elinewidth=ELINEWDTH, label="ZTF r",
            ms=ZMS,
        )
    except NameError:
        pass
    # if obj.get("plot") and obj["plot"].get("ztffilt") and "i" in obj["plot"].get("ztffilt"):
    try:
        plt.errorbar(
            datai["hjd"] - JD_SHIFT, datai["mag"], datai["magerr"], marker="o", ls="none",
            c="darkred", elinewidth=ELINEWDTH, label="ZTF i",
            ms=3,
        )
    except NameError:
        pass
    ztfnam = "-ztf" if obj.get("plot") and obj["plot"].get("ztffilt") else ""

    if obj.get("oglefnam"):
        ogled = pd.read_csv(
            f"{DATA_RAW}{obj['oglefnam']}", delim_whitespace=True
        )
        ogled.iloc[:, 0] = ogled.iloc[:, 0] + 50000
        plt.errorbar(
            ogled.iloc[:, 0], ogled.iloc[:, 1], ogled.iloc[:, 2], marker="o",
            ls="none", c="k", label="OGLE I",
            ms=3,
        )

        if obj.get("period"):
            ogled["phased"] = (
                (ogled.iloc[:, 0] + 2400000 - obj["epoch"]) % obj.get("period")
            ) / obj.get("period")
            print(ogled.iloc[:, 0] - obj["epoch"])
            dsh = ogled[ogled["phased"] > 0.5]
            dsh["phased"] = dsh["phased"] - 1
            ogled = pd.concat((ogled, dsh))

        oglei = pd.DataFrame({
            "hjd": ogled.iloc[:, 0] + 2400000,
            "mag": ogled.iloc[:, 1],
            "magerr": ogled.iloc[:, 2],
        })
        imerg = pd.concat((datai, oglei)).sort_values(by=["hjd"])
        pd.DataFrame({
            "hjd": imerg["hjd"].astype("str").str.ljust(18, "0"),
            "mag": imerg["mag"].round(6).astype("str").str.ljust(9, "0"),
            "magerr": imerg["magerr"].round(6).astype("str").str.ljust(8, "0"),
        }).to_csv(f"{DATA_PATH}{fnam}-ogle-ztf-imerg.dat", sep=" ", index=False)
        print(f"OGLE range for {name}: {min(ogled.iloc[:, 1])}-{max(ogled.iloc[:, 1])}")

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
    plt.xlabel("HJD", fontsize=16)
    if not obj.get("atlasfnam") and not obj.get("asasfnam"):
        plt.xlabel(f"HJD - {JD_SHIFT}", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.tight_layout()
    # plt.grid(axis="y", ls=":")
    plt.legend(fontsize=14)  # , loc="upper left"
    if obj.get("ylim"):
        plt.ylim(obj.get("ylim"))
    plt.gca().invert_yaxis()

    crtsnam = "-css" if obj.get("crtsfnam") else ""
    atlasfnam = "-atlas" if obj.get("atlasfnam") else ""
    asasfnam = "-asas" if obj.get("asasfnam") else ""
    oglenam = "-ogle" if obj.get("oglefnam") else ""
    gaianam = "-gaia" if obj.get("gaiafnam") else ""
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
        data_to_merge = []  # data in JD to be merged!
        ELINEWDTH = 0.8
        if "g" in obj["plot"].get("ztffilt"):
            plt.errorbar(
                datag["phased"], datag["mag"] - filtshift["g"], datag["magerr"],
                markeredgewidth=0.4, markeredgecolor="k", marker="o", ls="none",
                elinewidth=ELINEWDTH, c="g", label=f"ZTF g{mk_fsh(filtshift['g'])}",
                ms=4,
            )
            data_to_merge.append(
                pd.DataFrame({
                    "mjd": datag["mjd"],
                    "mag": datag["mag"] - filtshift["g"],
                    "magerr": datag["magerr"],
                })
            )
        if "r" in obj["plot"].get("ztffilt"):
            plt.errorbar(
                datar["phased"], datar["mag"], datar["magerr"],
                marker="o", markeredgewidth=0.4, markeredgecolor="k", ls="none",
                elinewidth=ELINEWDTH, c="r", label="ZTF r",
                ms=4,
            )
            data_to_merge.append(
                pd.DataFrame({
                    "mjd": datar["mjd"],
                    "mag": datar["mag"],
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
                        "mjd": datai["mjd"],
                        "mag": datai["mag"] - filtshift["i"],
                        "magerr": datai["magerr"],
                    })
                )
            except NameError:
                pass
        ATMS = obj["plot"].get("atlasms")
        ELINEWDTH = 0.2  # 0.5
        if obj["plot"].get("atlaselw"):
            ELINEWDTH = obj["plot"].get("atlaselw")
        if (
            obj.get("atlasfnam")
            and obj["plot"].get("atlasfilt")
            and "o" in obj["plot"].get("atlasfilt")
        ):
            alcurve_o["phased"] = (
                (alcurve_o["hjd"] - obj["epoch"]) % obj.get("period")
            ) / obj.get("period")
            dsh = alcurve_o[alcurve_o["phased"] > 0.5]
            dsh["phased"] = dsh["phased"] - 1
            alcurve_o = pd.concat((alcurve_o, dsh))
            data_to_merge.append(
                pd.DataFrame({
                    "mjd": alcurve_o["mjd"],
                    "mag": alcurve_o["mag"] - filtshift["o"],
                    "magerr": alcurve_o["magerr"],
                })
            )
            plt.errorbar(
                alcurve_o["phased"], alcurve_o["mag"] - filtshift["o"], alcurve_o["magerr"],
                marker="o", ls="none", elinewidth=ELINEWDTH, c="darkorange",
                label=f"ATLAS o{mk_fsh(filtshift['o'])}", zorder=0,
                ms=ATMS,
            )
        if (
            obj.get("atlasfnam")
            and obj["plot"].get("atlasfilt")
            and "c" in obj["plot"].get("atlasfilt")
        ):
            alcurve_c["phased"] = (
                (alcurve_c["hjd"] - obj["epoch"]) % obj.get("period")
            ) / obj.get("period")
            data_to_merge.append(
                pd.DataFrame({
                    "mjd": alcurve_c["mjd"],
                    "mag": alcurve_c["mag"] - filtshift["c"],
                    "magerr": alcurve_c["magerr"],
                })
            )
            dsh = alcurve_c[alcurve_c["phased"] > 0.5]
            dsh["phased"] = dsh["phased"] - 1
            alcurve_c = pd.concat((alcurve_c, dsh))
            plt.errorbar(
                alcurve_c["phased"], alcurve_c["mag"] - filtshift["c"], alcurve_c["magerr"],
                marker="o", ls="none", elinewidth=ELINEWDTH, c="darkblue",
                label=f"ATLAS c{mk_fsh(filtshift['c'])}", zorder=0,
                ms=ATMS,
            )
        if obj.get("oglefnam"):
            plt.errorbar(
                ogled["phased"], ogled.iloc[:, 1], ogled.iloc[:, 2], marker="o",
                ls="none", c="k", label="OGLE I",
                ms=5,
            )
        if obj.get("plot") and obj["plot"].get("psfilt"):
            for i, filter in enumerate(obj["plot"]["psfilt"]):
                if filter in psdata:
                    try:
                        psfiltshift = filtshift[filter]
                        if obj.get("clrshift").get("ps" + filter):
                            psfiltshift = obj.get("clrshift").get("ps" + filter)
                        ELINEWDTH = 1
                        plt.errorbar(
                            psdata[filter]["phased"], psdata[filter]["mag"] - psfiltshift,
                            psdata[filter]["magerr"],
                            marker=MARKERS[i], ls="none", c=FILT_CLRS[filter],
                            markeredgewidth=0.5, elinewidth=ELINEWDTH,
                            markeredgecolor="k", label=f"PS1 {filter}{mk_fsh(psfiltshift)}",
                            ms=MSS[i],
                        )
                        data_to_merge.append(
                            pd.DataFrame({
                                "mjd": psdata[filter]["mjd"],
                                "mag": psdata[filter]["mag"] - psfiltshift,
                                "magerr": psdata[filter]["magerr"],
                            })
                        )
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
            f"../lc/{fnam}{ztfnam}{atlasfnam}{oglenam}-phased_{obj.get('period')}.{EXT}",
            dpi=120,
        )
    # if loop on different objects, don't forget to delete user-defined
    # variables to avoid plotting wrong data
    # del psnam, ztfnam, crtsnam, atlasfnam, asasfnam
    # del datag, datar, datai
    # del atlas_c, atlas_o
