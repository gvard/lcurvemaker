import json
import warnings
from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
from astropy.utils.exceptions import AstropyWarning

import const
from curves import (
    Ps,
    Ztf,
    Atlas,
    Ogle,
    read_ps_data,
    read_crts_data,
    mk_phased,
    mk_fsh,
    save_datafile,
    save_merged,
    save_gaia_datafile,
)


warnings.simplefilter("ignore", category=AstropyWarning)
pd.options.mode.chained_assignment = None  # default="warn"

parser = ArgumentParser(
    description="Python script for working with light curves of variable stars"
)
parser.add_argument("-l", "--localps", action="store_true", help="local PS1 files")
parser.add_argument("-v", "--verln", action="store_true",
                    help="draw vertical line on phased plot")
parser.add_argument("-s", "--show", action="store_true",
                    help="Show interactive plot instead of saving figure")
parser.add_argument("-f", "--filename", type=str, default="",
                    help="json filename with object parameters")
args = parser.parse_args()

print("localps", args.localps, "filename", args.filename)
FN = "minkovskiy16.json"
if args.filename:
    FN = f"{args.filename}.json"
with open(f"../objects/{FN}", "r", encoding="utf8") as loc_file:
    objs = json.load(loc_file)

JD_SHIFT = 2400000.5
filtshift = const.ZEROFILTSHIFT

for name, obj in objs.items():
    fnam = name.lower().replace(" ", "")
    objpsfnam = f"{fnam}-ps"
    name_add = obj["other"]
    objname = f"{name} = {name_add}"
    ra, dec = obj["coordeg"]
    curveshift = filtshift
    if obj.get("clrshift"):
        filtshift.update(obj.get("clrshift"))
    if obj.get("curveshift") and obj.get("clrshift"):
        curveshift.update(obj.get("clrshift"))
    print(f"-=={name}==-")
    if obj.get("plot") and obj["plot"].get("psms"):
        const.PSMS = obj["plot"].get("psms")
    MSS = [const.PSMS - 0.5, const.PSMS + 1, const.PSMS, const.PSMS + 1, const.PSMS]

    # Get PS1 data
    const.PSFILT
    if obj.get("plot") and obj["plot"].get("psfilt"):
        const.PSFILT = obj["plot"].get("psfilt")
    if (
        not args.localps
        and obj.get("plot")
        and const.PSFILT
        and not obj.get("pslocal")
    ):
        Psdat = Ps(ra, dec)
        Psdat.cone_search()
        if obj.get("plot") and obj.get("plot").get("psfilt"):
            Psdat.PSFILTRS = const.PSFILT
        if Psdat.results:
            Psdat.get_table()
        else:
            print(Psdat.results, f"No PS1 data for {name}")

    # Light curve figure
    fig, ax = plt.subplots(figsize=(16, 9))
    fig.subplots_adjust(0.06, 0.09, 0.985, 0.95)

    # ASAS Stuff
    if obj.get("plot") and obj["plot"].get("asasfilt"):
        const.ASASFILT = obj["plot"].get("asasfilt")
    if obj.get("asasfnam"):
        asasd = pd.read_csv(f"{const.DATA_RAW}{obj['asasfnam']}")
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
        ASASERRLIM = obj.get("asaslim") if obj.get("asaslim") else const.ASAS_MAG_LIM
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
        save_datafile(asasd, f"{const.DATA_PATH}{fnam}-asassn.dat", hjdprec=6, magprec=3)
        pd.DataFrame({
            "hjd": asasd["hjd"].astype("str").str.ljust(13, "0"),
            "mag": asasd["mag"].round(6).astype("str").str.ljust(6, "0"),
            "magerr": asasd["magerr"].astype("str").str.ljust(5, "0"),
            "filter": asasd["filter"],
        }).to_csv(f"{const.DATA_PATH}{fnam}-asassn.dat", sep=" ", index=False)
        if obj["plot"].get("asasupperlim") and len(asasd_bad["mag"]):
            plt.plot(asasd_bad["hjd"] - JD_SHIFT, asasd_bad["mag"], marker="v", ls="none",
                     c="#bbb", ms=2, zorder=-2)  # , label="ASAS-SN"
        if "V" in const.ASASFILT and len(asasd_V["mag"]):
            plt.errorbar(
                asasd_V["hjd"] - JD_SHIFT, asasd_V["mag"] - curveshift["V"], asasd_V["magerr"],
                marker="s", ls="none", c=const.FILT_CLRS["V"], markeredgecolor="k", mew=const.MEW, zorder=-1,
                label=f"ASAS-SN V{mk_fsh(curveshift['V'])}", ms=const.ASASMS-1, elinewidth=const.ASELINWDTH,
            )
            print(
                f'ASAS-SN range for {name}, V: {min(asasd_V["mag"])}-{max(asasd_V["mag"])}'
            )
        if "g" in const.ASASFILT and len(asasd_g["mag"]):
            plt.errorbar(
                asasd_g["hjd"] - JD_SHIFT, asasd_g["mag"] - curveshift["asasg"], asasd_g["magerr"],
                marker="d", ls="none", c=const.FILT_CLRS["asasg"], markeredgecolor="k", zorder=-1,
                mew=const.MEW, label=f"ASAS-SN g{mk_fsh(curveshift['asasg'])}", ms=const.ASASMS,
                elinewidth=const.ASELINWDTH,
            )
            print(
                f'ASAS-SN range for {name}, g: {min(asasd_g["mag"])}-{max(asasd_g["mag"])}'
            )

    if obj.get("gaiafnam"):
        gaiadata = pd.read_csv(f"{const.DATA_RAW}{obj['gaiafnam']}", skiprows=1)
        gaiadata = gaiadata[gaiadata["averagemag"].notna()]
        gaiadata = gaiadata[gaiadata["averagemag"] != "untrusted"]
        gaiadata["averagemag"] = gaiadata["averagemag"].astype(float)
        # Time of observation is in barycentric coordinate time (TCB)
        save_gaia_datafile(gaiadata, f"../data/{obj['gaiafnam'].removesuffix('.csv')}.dat")
        plt.plot(
            gaiadata["JD(TCB)"] - JD_SHIFT, gaiadata.get("averagemag"), "*", markeredgecolor="k",
            mew=const.MEW, ms=const.GAIAMS, c=const.FILT_CLRS["G"], label=obj["gaiaobj"],
        )
        if obj.get("period"):
            gaiadata = mk_phased(gaiadata, obj["epoch"], obj.get("period"), jdnam="JD(TCB)")

    # Atlas Stuff
    const.ATLAS_SAVEDATA = True
    if obj.get("atlaslim"):
        const.ATLAS_MAG_LIM = obj.get("atlaslim")
    if obj["plot"].get("atlasfilt"):
        const.ATLASFILT = obj["plot"].get("atlasfilt")
    Atl = Atlas(fnam, atlaslim=const.ATLAS_MAG_LIM)
    if not const.ATLAS_RAWDATA and Atl.is_data_exist("o", lim=const.ATLAS_MAG_LIM):
        Atl.read_prepared_data(filtlims={"o": const.ATLAS_MAG_LIM, "c": const.ATLAS_MAG_LIM})
        print("read atlas data for", fnam, "with lim", const.ATLAS_MAG_LIM)
    elif const.ATLAS_RAWDATA and Atl.is_raw_data_exist():
        Atl.read_raw_data()
        Atl.prepare_data(ra, dec)
        print("read and prepare atlas raw data for", fnam, "with lim", const.ATLAS_MAG_LIM)
        if const.ATLAS_SAVEDATA:
            print(f"save atlas data for {const.ATLASFILT} filters")
            for filtr in const.ATLASFILT:
                Atl.save_datafile(filtr)
    else:
        print("No Atlas data", hasattr(Atl, "data"))

    if obj["plot"].get("atlasms"):
        const.ATLASMS = obj["plot"].get("atlasms")
    if hasattr(Atl, "data") and len(Atl.data.index):
        for afiltr in const.ATLASFILT:
            alclr = const.FILT_CLRS[afiltr]
            alms = {"o": const.ATLASMS, "c": const.ATLASMS - 1}[afiltr]
            alcurve = Atl.get_data(afiltr)
            if obj.get("curveshift") and obj.get("clrshift") and const.ATLASFILT:
                plt.errorbar(
                    alcurve["hjd"] - JD_SHIFT, alcurve["mag"] - curveshift[afiltr], alcurve["magerr"],
                    marker="o", ls="none", elinewidth=const.ATELINWDTH, c=alclr,
                    ms=alms, label=f"ATLAS {afiltr}{mk_fsh(curveshift[afiltr])}", zorder=0,
                    markeredgecolor="k", mew=const.MEW,
                )
            elif const.ATLASFILT:
                plt.errorbar(
                    alcurve["hjd"] - JD_SHIFT, alcurve["mag"], alcurve["magerr"],
                    marker="o", ls="none", elinewidth=const.ATELINWDTH, c=alclr,
                    ms=alms, label=f"ATLAS {afiltr}", zorder=0, markeredgecolor="k", mew=const.MEW,
                )
    if obj.get("period"):
        obj["period"] = float(obj["period"])
    if obj.get("epoch"):
        obj["epoch"] = float(obj["epoch"])

    # -== ZTF stuff ==-
    if obj["plot"].get("ztffilt"):
        const.ZTFFILT = obj["plot"].get("ztffilt")
    if obj.get("plot") and obj["plot"].get("zms"):
        const.ZMS = obj["plot"].get("zms")
    zmss = {"g": const.ZMS, "r": const.ZMS, "i": const.ZMS}

    Zt = Ztf(fnam)
    zdata = {}
    for filtr in const.ZTFFILT:
        data = Zt.read_raw_data_filtr(filtr=filtr)
        try:
            save_datafile(data, f"{const.DATA_PATH}{fnam}-ztf{filtr}.dat", hjdprec=7, magprec=6)
            print(
                f'ZTF range for {name} {filtr}: {round(min(data["mag"]), 2)}-{round(max(data["mag"]), 2)}'
            )
            zdata[filtr] = data
            if obj.get("period"):
                zdata[filtr] = mk_phased(data, obj["epoch"], obj.get("period"))

            if obj.get("curveshift") and obj.get("clrshift"):
                plt.errorbar(
                    zdata[filtr]["hjd"] - JD_SHIFT, zdata[filtr]["mag"] - curveshift[filtr], zdata[filtr]["magerr"],
                    marker="o", ls="none", c=const.FILT_CLRS[filtr], elinewidth=const.ZELINWDTH-0.3,
                    label=f"ZTF {filtr}{mk_fsh(curveshift[filtr])}",
                    ms=zmss[filtr], markeredgecolor="k", mew=const.MEW,
                )
            else:
                plt.errorbar(
                    zdata[filtr]["hjd"] - JD_SHIFT, zdata[filtr]["mag"], zdata[filtr]["magerr"],
                    marker="o", ls="none", c=const.FILT_CLRS[filtr], elinewidth=const.ZELINWDTH,
                    label=f"ZTF {filtr}",
                    ms=zmss[filtr], markeredgecolor="k", mew=const.MEW,
                )
        except TypeError:
            continue
        except NameError:
            print(f"Err {filtr}")
    ztfnam = "-ztf" if obj.get("plot") and const.ZTFFILT else ""

    if obj.get("oglefnam"):
        Ogl = Ogle(obj.get("oglefnam"))
        Ogl.read_raw_data()
        plt.errorbar(
            Ogl.data["hjd"] - JD_SHIFT, Ogl.data["mag"] - filtshift["I"], Ogl.data["magerr"],
            marker="o", ls="none", c=const.FILT_CLRS["I"], label=f"OGLE I{mk_fsh(filtshift['I'])}",
            ms=const.OGLEMS, markeredgecolor="k", mew=const.MEW, elinewidth=const.OELINWDTH,
        )

        if obj.get("period"):
            Ogl.mk_phased(obj["epoch"], obj["period"])

        imerg = pd.concat((zdata["i"], Ogl.data)).sort_values(by=["hjd"])
        pd.DataFrame({
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
    if not args.localps and obj.get("plot") and obj.get("plot").get("psfilt") and Psdat.results:
        psnam = "-ps1"
        for i, filter in enumerate(const.PSFILT):
            w = np.where(Psdat.dtab["filter"] == filter)
            try:
                print(
                    f"PS1 for {name}, filter {filter}: {str(round(min(Psdat.mag[w]), 2)).ljust(5, '0')}-{str(round(max(Psdat.mag[w]), 2)).ljust(5, '0')}"
                )
                plt.errorbar(
                    Psdat.t[w], Psdat.mag[w], Psdat.magerr[w],
                    marker=const.MARKERS[i], ls="none", c=const.FILT_CLRS[filter],
                    mew=const.MEW, markeredgecolor="k", label=f"PS1 {filter}", lw=0.75,
                    ms=MSS[i],
                )
                psdata[filter] = Psdat.mk_data(w)
                save_datafile(
                    psdata[filter], f"{const.DATA_PATH}{objpsfnam}{filter}.dat",
                    hjdprec=7, magprec=4,
                )
                psdata[filter] = psdata[filter][psdata[filter]["magerr"] < const.PS_MAG_LIM]
                if obj.get("period"):
                    psdata[filter] = mk_phased(
                        psdata[filter], obj["epoch"], obj.get("period")
                    )
            except (IndexError, ValueError):
                print(f"no PS1 plot for {filter}")
    if args.localps and obj.get("plot") and obj.get("plot").get("psfilt"):
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
                    marker=const.MARKERS[i], ls="none", c=const.FILT_CLRS[filter],
                    mew=const.MEW, markeredgecolor="k",
                    label=f"PS1 {filter}{mk_fsh(psfiltshift)}", lw=0.75,
                    ms=MSS[i],
                )
                if obj.get("period"):
                    psdata[filter] = pd.DataFrame({
                        "hjd": pslc["hjd"],
                        "mag": pslc["mag"],
                        "magerr": pslc["magerr"],
                    })
                    psdata[filter] = psdata[filter][psdata[filter]["magerr"] < const.PS_MAG_LIM]
                    psdata[filter] = mk_phased(
                        psdata[filter], obj["epoch"], obj.get("period")
                    )
            except (IndexError, ValueError, FileNotFoundError):
                print(f"no PS1 plot for {filter}")

    if obj.get("crtsfnam"):
        crts_data = read_crts_data(f"{const.DATA_PATH}{obj['crtsfnam']}")
        plt.errorbar(
            crts_data["hjd"] - JD_SHIFT, crts_data["mag"], crts_data["magerr"],
            marker="o", ls="none", c="#888", elinewidth=const.CELINWDTH, label="CRTS",
            markeredgecolor="k", mew=const.MEW+0.1, ms=const.CRTSMS,
        )

    plt.title(objname, fontsize=18)
    # ax.ticklabel_format(useOffset=False, style="plain")
    if obj.get("plot"):
        if obj["plot"].get("xmal"):
            ax.xaxis.set_major_locator(MultipleLocator(obj["plot"]["xmal"]))
        if obj["plot"].get("xmil"):
            ax.xaxis.set_minor_locator(MultipleLocator(obj["plot"]["xmil"]))
        # ax.yaxis.set_major_locator(MultipleLocator(0.5))
        if obj["plot"].get("ymil"):
            YMIL = obj["plot"].get("ymil")
            ax.yaxis.set_minor_locator(MultipleLocator(YMIL))
        if obj["plot"].get("ymal"):
            YMAL = obj["plot"].get("ymal")
            ax.yaxis.set_major_locator(MultipleLocator(YMAL))
        if obj["plot"].get("xedges"):
            xmin, xmax = ax.get_xlim()
            plt.xlim(xmin + obj["plot"]["xedges"], xmax - obj["plot"]["xedges"])
        elif obj["plot"].get("xlima"):
            plt.xlim(obj["plot"].get("xlima"))
    plt.ylabel("Mag", fontsize=16)
    plt.xlabel(f"HJD - {JD_SHIFT}", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.tight_layout()
    # plt.grid(axis="y", ls=":")
    plt.legend(fontsize=14)  # , loc="upper left"
    if obj["plot"].get("ylima"):
        ax.set_ylim(obj["plot"].get("ylima"))
    else:
        ax.invert_yaxis()

    crtsnam = "-css" if obj.get("crtsfnam") else ""
    atlasfnam = "-atlas" if hasattr(Atl, "data") and const.ATLASFILT else ""
    asasfnam = "-asas" if obj.get("asasfnam") else ""
    oglenam = "-ogle" if obj.get("oglefnam") else ""
    gaianam = "-gaia" if obj.get("gaiafnam") else ""
    if args.show:
        plt.show()
    else:
        plt.savefig(
            f"../lc/{fnam}{psnam}{ztfnam}{crtsnam}{atlasfnam}{asasfnam}{oglenam}{gaianam}.{const.EXT}",
            dpi=120,
        )
    # End of lightcurve plot
    plt.close()

    if obj.get("period"):
        # Phased plot figure
        fig2, ax2 = plt.subplots(figsize=(16, 9))
        fig2.subplots_adjust(0.06, 0.09, 0.985, 0.95)
        data_to_merge = []  # data in HJD to be merged
        for filtr in const.ZTFFILT:
            plt.errorbar(
                zdata[filtr]["phased"], zdata[filtr]["mag"] - filtshift[filtr], zdata[filtr]["magerr"],
                mew=const.MEW, markeredgecolor="k", marker="o", ls="none",
                elinewidth=const.ZELINWDTH, c=const.FILT_CLRS[filtr],
                label=f"ZTF {filtr}{mk_fsh(filtshift[filtr])}", ms=zmss[filtr],
            )
            data_to_merge.append(
                pd.DataFrame({
                    "hjd": zdata[filtr]["hjd"],
                    "mag": zdata[filtr]["mag"] - filtshift[filtr],
                    "magerr": zdata[filtr]["magerr"],
                })
            )

        if obj["plot"].get("atlaselw"):
            const.ATELINWDTH = obj["plot"].get("atlaselw")
        if hasattr(Atl, "data") and len(Atl.data.index) and const.ATLASFILT:
            Atl.mk_phased(obj["epoch"], obj.get("period"))
            for afiltr in const.ATLASFILT:
                alclr = const.FILT_CLRS[afiltr]
                alms = {"o": const.ATLASMS, "c": const.ATLASMS - 1}[afiltr]
                alcurve = Atl.get_data(afiltr)
                plt.errorbar(
                    alcurve["phased"], alcurve["mag"] - filtshift[afiltr], alcurve["magerr"],
                    marker="o", ls="none", elinewidth=const.ATELINWDTHT, c=alclr,
                    label=f"ATLAS {afiltr}{mk_fsh(filtshift[afiltr])}", zorder=0,
                    ms=alms, mew=const.MEW, markeredgecolor="k",
                )
                data_to_merge.append(
                    pd.DataFrame({
                        "hjd": alcurve["hjd"],
                        "mag": alcurve["mag"] - filtshift[afiltr],
                        "magerr": alcurve["magerr"],
                    })
                )
        if "V" in const.ASASFILT and len(asasd_V["mag"]):
            plt.errorbar(
                asasd_V["phased"], asasd_V["mag"], asasd_V["magerr"],
                marker="s", ls="none", c=const.FILT_CLRS["V"], markeredgecolor="k", mew=const.MEW,
                label="ASAS-SN V", ms=const.ASASMS, elinewidth=const.ASELINWDTH,
            )
        if "g" in const.ASASFILT and len(asasd_g["mag"]):
            plt.errorbar(
                asasd_g["phased"], asasd_g["mag"] - curveshift["asasg"], asasd_g["magerr"],
                marker="d", ls="none", c=const.FILT_CLRS["asasg"], markeredgecolor="k",
                mew=const.MEW, label=f"ASAS-SN g {mk_fsh(curveshift['asasg'])}",
                ms=const.ASASMS, elinewidth=const.ASELINWDTH,
            )

        if obj.get("oglefnam"):
            plt.errorbar(
                Ogl.data["phased"], Ogl.data["mag"] - filtshift["I"], Ogl.data["magerr"],
                marker="o", ls="none", c=const.FILT_CLRS["I"], label=f"OGLE I{mk_fsh(filtshift['I'])}", ms=const.OGLEMS,
                markeredgecolor="k", mew=const.MEW, elinewidth=const.OELINWDTH,
            )
        if obj.get("gaiafnam"):
            plt.plot(
                        gaiadata["phased"], gaiadata.get("averagemag") - filtshift["G"], "*", markeredgecolor="k",
                        mew=const.MEW, ms=const.GAIAMS, c=const.FILT_CLRS["G"], label=f"{obj['gaiaobj']}{mk_fsh(filtshift['G'])}",
                    )
            # data_to_merge.append(
            #     pd.DataFrame({
            #         "hjd": gaiadata["JD(TCB)"],
            #         "mag": gaiadata["averagemag"] - filtshift["G"],
            #         "magerr": 0.1,
            #     })
            # )

        if const.PSFILT:
            for i, filter in enumerate(const.PSFILT):
                if filter in psdata:
                    try:
                        psfiltshift = filtshift[filter]
                        if filtshift.get("ps" + filter):
                            psfiltshift = filtshift.get("ps" + filter)
                        plt.errorbar(
                            psdata[filter]["phased"], psdata[filter]["mag"] - psfiltshift,
                            psdata[filter]["magerr"], marker=const.MARKERS[i], ls="none",
                            c=const.FILT_CLRS[filter], mew=const.MEW, elinewidth=const.PSELINWDTH,
                            markeredgecolor="k", label=f"PS1 {filter}{mk_fsh(psfiltshift)}",
                            ms=MSS[i],
                        )
                        data_to_merge.append(
                            pd.DataFrame({
                                "hjd": psdata[filter]["hjd"],
                                "mag": psdata[filter]["mag"] - psfiltshift,
                                "magerr": psdata[filter]["magerr"],
                            })
                        )
                    except KeyError:
                        print("pass ps filter", filter)
        merged = pd.concat(data_to_merge)  # .sort_values(by=["hjd"])
        save_merged(merged, fnam)

        # Phased plot parameters
        plt.legend(fontsize=14)  # loc="lower left"
        plt.title(objname, fontsize=18)
        plt.ylabel("Mag", fontsize=16)
        plt.xlabel("Phase", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim(-0.5, 1)
        YLIM = obj["plot"].get("ylim")
        if YLIM:
            ax2.set_ylim(YLIM)
            if args.verln:
                plt.plot([0, 0], YLIM, "--k")
        else:
            plt.gca().invert_yaxis()
        if obj["plot"].get("ymal"):
            ax2.yaxis.set_major_locator(MultipleLocator(YMAL))
        if obj["plot"].get("ymil"):
            ax2.yaxis.set_minor_locator(MultipleLocator(YMIL))
        ax2.xaxis.set_major_locator(MultipleLocator(0.5))
        ax2.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax2.tick_params(which="major", length=6)
        ax2.tick_params(which="minor", length=4)
        plt.savefig(
            f"../lc/{fnam}{ztfnam}{atlasfnam}{asasfnam}{oglenam}{gaianam}-ph_{round(obj.get('period'), 6)}.{const.EXT}",
            dpi=120,
        )
    # if loop on different objects, don't forget to delete user-defined
    # variables to avoid plotting wrong data
    # del psnam, ztfnam, crtsnam, atlasfnam, asasfnam
    # del datag, datar, datai
    # del atlas_c, atlas_o
