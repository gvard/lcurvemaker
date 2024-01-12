import os

import numpy as np
import pandas as pd

from astropy.io import ascii
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

from panstarrs import ps1cone, ps1search, addfilter


HALEAKALA = EarthLocation(lon=-156.257, lat=20.71, height=3048)
PALOMAR = EarthLocation(lon=-116.86194, lat=33.3581, height=1712)
JD_SHIFT = 2400000.5
OGLE_SHIFT = 2450000
ABZPMAG_JY = 8.9
LGE_25 = 2.5 / np.log(10.0)


class Ps:
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec
        self.RADIUS = 2.1 / 3600.0 / 4
        self.release = "dr2"

    def cone_search(self):
        # strip blanks and weed out blank and commented-out values
        columns = """objID,raMean,decMean,nDetections,ng,nr,ni,nz,ny,gMeanPSFMag,gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag""".split(
            ","
        )
        # gMeanPSFMag,gMeanPSFMagErr,rMeanPSFMag,rMeanPSFMagErr,iMeanPSFMag,iMeanPSFMagErr,zMeanPSFMag,zMeanPSFMagErr,yMeanPSFMag,yMeanPSFMagErr""".split(',')
        columns = [x.strip() for x in columns]
        columns = [x for x in columns if x and not x.startswith("#")]
        CONSTRAINTS = {"nDetections.gt": 1}
        self.PSFILTRS = "griz"  # possible value: "grizy"
        self.results = ps1cone(
            self.ra, self.dec, self.RADIUS, release=self.release, columns=columns, **CONSTRAINTS
        )

    def get_table(self):
        tab = ascii.read(self.results)
        for filter in self.PSFILTRS:
            col = filter + "MeanPSFMag"
            # print(f"{tab[col+'Err']}")
            tab[col].format = ".4f"
            tab[col][tab[col] == -999.0] = np.nan
        objid = tab["objID"][0]
        dconstraints = {"objID": objid}
        dcolumns = (
            """objID,detectID,filterID,obsTime,ra,dec,psfFlux,psfFluxErr,psfMajorFWHM,psfMinorFWHM,psfQfPerfect,apFlux,apFluxErr,infoFlag,infoFlag2,infoFlag3"""
        ).split(",")
        # strip blanks and weed out blank and commented-out values
        dcolumns = [x.strip() for x in dcolumns]
        dcolumns = [x for x in dcolumns if x and not x.startswith("#")]
        dresults = ps1search(
            table="detection", release=self.release, columns=dcolumns, **dconstraints
        )
        self.dtab = addfilter(ascii.read(dresults))
        self.dtab.sort("obsTime")
        self.t = self.dtab["obsTime"]
        self.dtab = self.dtab[self.dtab["psfFlux"] != -999.0]
        self.mag = -2.5 * np.log10(self.dtab["psfFlux"]) + 8.90
        self.magerr = LGE_25 * self.dtab["psfFluxErr"] / self.dtab["psfFlux"]

    def mk_data(self, w):
        t_corr = mk_hjd_corr(self.t[w], self.ra, self.dec, obs="Haleakala")
        return pd.DataFrame({
                "mjd": self.t[w],
                "hjd": t_corr,
                "mag": self.mag[w],
                "magerr": self.magerr[w],
            })


class Data:
    def __init__(self):
        self.data_dir = "../data"
        self.raw_dir = f"{self.data_dir}/raw"
        self.name = ""
        self.survey = ""
        self.maglim = None
        self.columns = None

    def is_raw_data_exist(self, add="", ext="csv"):
        return os.path.isfile(f"{self.raw_dir}/{self.name}-{self.survey}{add}.{ext}")

    def is_data_exist(self, filtr, lim="", ext="dat"):
        return os.path.isfile(f"{self.data_dir}/{self.name}-{self.survey}{filtr}{lim}.{ext}")

    def read_raw_data(self, add="", ext="csv"):
        self.raw_data = pd.read_csv(f"{self.raw_dir}/{self.name}-{self.survey}{add}.{ext}",
                                    delim_whitespace=True)
        if self.columns:
            self.raw_data = self.raw_data.rename(columns=self.columns)


class Atlas(Data):
    def __init__(self, name, atlaslim=0.5):
        super().__init__()
        self.name = name
        self.survey = "atlas"
        self.maglim = atlaslim
        self.columns = {"###MJD": "mjd", "m": "mag", "dm": "magerr", "F": "filter"}

    def is_raw_data_exist(self, add="", ext="txt"):
        return super().is_raw_data_exist(add=add, ext=ext)

    def is_data_exist(self, filtr, lim=0.5, ext="dat"):
        return super().is_data_exist(filtr, lim=f"-cleaned-{lim}m", ext=ext)

    def read_raw_data(self, add="", ext="txt"):
        super().read_raw_data(add=add, ext=ext)

    def prepare_data(self, ra, dec, maglim_up=10, maglim_low=None):
        data = self.raw_data[self.raw_data["magerr"] < self.maglim]
        data = data[data["mag"] > maglim_up]
        if maglim_low:
            data = data[data["mag"] < maglim_low]
        data["hjd"] = mk_hjd_corr(data["mjd"], ra, dec, obs="Haleakala")
        self.data = pd.DataFrame({
            "hjd": data["hjd"],
            # "mjd": data["mjd"],
            "filter": data["filter"],
            "mag": data["mag"],
            "magerr": data["magerr"]
        })

    def mk_phased(self, epoch, period):
        self.data = mk_phased(self.data, epoch, period, jdnam="hjd")

    def save_datafile(self, filtr="o", ext="dat", hjdprec=6, magprec=3):
        atlas_data_fnam = f"{self.data_dir}/{self.name}-{self.survey}{filtr}-cleaned-{round(self.maglim, 2)}m.{ext}"
        data = self.data[self.data["filter"] == filtr]
        if len(data.index):
            save_datafile(data, atlas_data_fnam, hjdprec=hjdprec, magprec=magprec)

    def read_prepared_data(self, filtlims={"o": 0.5, "c": 0.5}, ext="dat"):
        alldata = []
        for filtr, lim in filtlims.items():
            if self.is_data_exist(filtr, lim=lim):
                data = pd.read_csv(
                    f"{self.data_dir}/{self.name}-{self.survey}{filtr}-cleaned-{round(lim, 2)}m.{ext}",
                    delim_whitespace=True
                    )
                data["filter"] = filtr
                alldata.append(data)
        self.data = pd.concat(alldata)

    def get_data(self, filtr="o"):
        return self.data[self.data["filter"] == filtr]


class Ogle(Data):
    def __init__(self, name):
        super().__init__()
        self.name = name
        self.survey = "ogle"
        self.maglim = 1.6

    def read_raw_data(self, add="", ext="dat"):
        columns = ["hjd", "mag", "magerr"]
        self.data = pd.read_csv(
            f"{self.raw_dir}/{self.name}", delim_whitespace=True,
            names=columns)
        self.data["hjd"] += OGLE_SHIFT

    def mk_phased(self, epoch, period):
        self.data = mk_phased(self.data, epoch, period, jdnam="hjd")


class Ztf(Data):
    def __init__(self, name, ztflim=None):
        super().__init__()
        self.name = name
        self.survey = "ztf"
        self.maglim = ztflim

    def is_raw_data_exist(self, add="", ext="csv"):
        return os.path.isfile(f"{self.raw_dir}/{self.name}-z{add}.{ext}")

    def is_data_exist(self, filtr, lim="", ext="dat"):
        return super().is_data_exist(filtr, lim=lim, ext=ext)

    def read_raw_data_filtr(self, filtr="r", ext="csv"):
        all_data = []
        for i in range(1, 4):
            try:
                ztf_fnam = f"{self.raw_dir}/{self.name}-z{filtr}{i}.{ext}"
                all_data.append(pd.read_csv(ztf_fnam))
            except FileNotFoundError:
                pass
        if all_data:
            data = pd.concat(all_data)
            data = data.sort_values(by=["hjd"])
            data = data.reset_index(drop=True)
            return pd.DataFrame({
                "hjd": data["hjd"],
                "mag": data["mag"],
                "magerr": data["magerr"],
                "filter": filtr,
            })
        else:
            return False

    def read_raw_data(self, filtrs="gri", ext="csv"):
        for filtr in filtrs:
            data = self.read_raw_data_filtr(filtr=filtr)
            if data:
                try:
                    self.data = pd.concat(self.data, data)
                except ValueError:
                    self.data = data

    def prepare_data(self, maglim_up=None, maglim_low=None):
        self.data = self.data[self.data["magerr"] < self.maglim]
        if maglim_up:
            self.data = self.data[self.data["mag"] > maglim_up]
        if maglim_low:
            self.data = self.data[self.data["mag"] < maglim_low]

    def mk_phased(self, epoch, period):
        self.data = mk_phased(self.data, epoch, period, jdnam="hjd")

    def save_datafile(self, filtr="r", ext="dat", hjdprec=7, magprec=6):
        data = self.data[self.data["filter"] == filtr]
        if len(data.index):
            save_datafile(data, f"{self.data_dir}/{self.name}-ztf{filtr}.{ext}", hjdprec=hjdprec, magprec=magprec)

    def read_prepared_data_filtr(self, filtr="r", ext="dat"):
        if self.is_data_exist(filtr):
            data = pd.read_csv(
                f"{self.data_dir}/{self.name}-{self.survey}{filtr}.{ext}",
                delim_whitespace=True
                )
            data["filter"] = filtr
            return data
        else:
            return False

    def read_prepared_data(self, filtrs="gri", ext="dat"):
        for filtr in filtrs:
            data = self.read_prepared_data_filtr(filtr=filtr)
            if data:
                try:
                    self.data = pd.concat(data, self.data)
                except ValueError:
                    self.data = data

    def get_data(self, filtr="r"):
        return self.data[self.data["filter"] == filtr]


def read_ps_data(filename):
    """read data: mjd mag magerr"""
    return pd.read_csv(filename, delim_whitespace=True)


def read_crts_data(filename):
    """Read CRTS data"""
    columns = ["hjd", "mag", "magerr"]
    return pd.read_csv(filename, delim_whitespace=True, names=columns)


def save_gaia_datafile(curve, fnsav, hjdprec=5, magprec=2):
    """Save lightcurve data as ASCII file with space separated values"""
    pd.DataFrame({
        "hjd": curve["JD(TCB)"].round(hjdprec).astype("str").str.ljust(hjdprec + 8, "0"),
        "mag": curve["averagemag"].astype("str").str.ljust(magprec + 3, "0"),
    }).to_csv(fnsav, sep=" ", index=False)


def save_datafile(curve, fnsav, hjdprec=7, magprec=6, addfilter=False):
    """Save lightcurve data as ASCII file with space separated values"""
    data = pd.DataFrame({
        "hjd": curve["hjd"].round(hjdprec).astype("str").str.ljust(hjdprec + 8, "0"),
        "mag": curve["mag"].round(magprec).astype("str").str.ljust(magprec + 3, "0"),
        "magerr": curve["magerr"].round(magprec).astype("str").str.ljust(magprec + 2, "0"),
    })
    if addfilter and "filter" in curve:
        data["filter"] = curve["filter"]
    data.to_csv(fnsav, sep=" ", index=False)


def mk_phased(data, epoch, period, jdnam="hjd"):
    """Make phased lightcurve data"""
    data["phased"] = ((data[jdnam] - epoch) % period) / period
    dsh = data[data["phased"] > 0.5]
    dsh["phased"] = dsh["phased"] - 1
    return pd.concat((data, dsh))


def save_merged(data, fnam):
    "save merged data to .dat file"
    pd.DataFrame({
        "hjd": (data["hjd"]).round(6).astype("str").str.ljust(14, "0"),
        "mag": data["mag"].round(5).astype("str").str.ljust(8, "0"),
        "magerr": data["magerr"].round(5).astype("str").str.ljust(7, "0"),
    }).to_csv(f"../data/{fnam}-allmerged.dat", sep=" ", index=False)


def mk_hjd_corr(dates, ra, dec, obs="Palomar"):
    """Add heliocentric correction to MJD dates, return HJD dates"""
    coord = SkyCoord(ra=ra, dec=dec, unit="deg")
    time = Time(Time(dates, format="mjd"), scale="utc")  # scale="tai"
    loc = HALEAKALA
    if obs == "Palomar":
        loc = PALOMAR
    helio_time = time + time.light_travel_time(coord, "heliocentric", location=loc)
    # Add 15s shift derived from ZTF raw data
    return helio_time.mjd + JD_SHIFT + 0.00017361


def mk_fsh(filtshift):
    """Generate text for plot labels"""
    if filtshift == 0:
        return ""
    elif float(filtshift) < 0:
        return f" +{abs(filtshift)}"
    else:
        return f" −{filtshift}"  # :+1.1f
