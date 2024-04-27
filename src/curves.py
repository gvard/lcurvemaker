import os

import numpy as np
import pandas as pd

from astropy.io import ascii
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

from panstarrs import ps1cone, ps1search, addfilter


HALEAKALA = EarthLocation(lon=-156.257, lat=20.71, height=3040)
ATLAS_MLO = EarthLocation(lon=-155.5761, lat=19.5361, height=3397)
ATLAS_CHL = EarthLocation(lon=-70.76498, lat=-30.47103, height=1575)
ATLAS_STH = EarthLocation(lon=20.8105, lat=-32.3783, height=1766)
PALOMAR = EarthLocation(lon=-116.86194, lat=33.3581, height=1712)
CATALINA = EarthLocation(lon=-110.789, lat=32.442, height=2791)
SIDINGSPRINGS = EarthLocation(lon=149.1, lat=-31.3, height=1150)
JD_SHIFT = 2400000.5
OGLE_SHIFT = 2450000
ABZPMAG_JY = 8.9
LGE_25 = 2.5 / np.log(10.0)


class Ps:
    def __init__(self, ra, dec, radius=0.53):
        """radius in arcseconds"""
        self.ra = ra
        self.dec = dec
        self.RADIUS = radius / 3600.0
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
                                    sep="\s+")
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

    def prepare_data(self, ra, dec, maglim_up=10, maglim_low=None, obs="Haleakala"):
        data = self.raw_data[self.raw_data["magerr"] < max(self.maglim.values())]
        data = data.drop(data[(data["filter"] == "o") & (data["magerr"] >= self.maglim["o"])].index)
        data = data.drop(data[(data["filter"] == "c") & (data["magerr"] >= self.maglim["c"])].index)
        data = data[data["mag"] > maglim_up]
        if maglim_low:
            data = data[data["mag"] < maglim_low]
        data["hjd"] = mk_hjd_corr(data["mjd"], ra, dec, obs=obs)
        self.data = pd.DataFrame({
            "hjd": data["hjd"],
            "mag": data["mag"],
            "magerr": data["magerr"],
            "filter": data["filter"],
        })

    def clean_by_filter(self, filtr="c", maglim=0.035):
        data = self.data[self.data["filter"] == filtr]
        return data[data["magerr"] < maglim]

    def mk_phased(self, epoch, period):
        self.data = mk_phased(self.data, epoch, period, jdnam="hjd")

    def save_datafile(self, filtr="o", ext="dat", hjdprec=6, magprec=3):
        atlas_data_fnam = f"{self.data_dir}/{self.name}-{self.survey}{filtr}-clean{round(self.maglim[filtr], 2)}.{ext}"
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

    def get_data(self, filtr="o", maglim=None):
        data = self.data[self.data["filter"] == filtr]
        if maglim:
            data = data[data["magerr"] < maglim]
        return data


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
                "catflags": data["catflags"]
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

    def prepare_data(self, maglim_up=None, maglim_low=None, catfilt=True):
        if self.maglim:
            self.data = self.data[self.data["magerr"] < self.maglim]
        if maglim_up:
            self.data = self.data[self.data["mag"] > maglim_up]
        if maglim_low:
            self.data = self.data[self.data["mag"] < maglim_low]
        if catfilt:
            self.data = self.data[self.data["catflags"] != 32768]

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
            data = self.read_prepared_data_filtr(filtr=filtr, ext=ext)
            if data:
                try:
                    self.data = pd.concat(data, self.data)
                except ValueError:
                    self.data = data

    def get_data(self, filtr="r"):
        return self.data[self.data["filter"] == filtr]


class Alerce(Data):
    def __init__(self, objs):
        super().__init__()
        self.filters = {1: "g", 2: "r",  3: "i"}
        self.objects = objs

    def read_raw_data(self, add="detections", ext="csv"):
        for obj in self.objects:
            self.raw_data = pd.read_csv(f"{self.raw_dir}/{obj}{add}.{ext}")
            self.raw_data["filter"] = self.raw_data["fid"].replace(self.filters)

    def prepare_data(self):
        self.data = pd.DataFrame({
            "hjd": self.raw_data["mjd"] + JD_SHIFT,
            "mag": self.raw_data["magpsf"],
            "magerr": self.raw_data["sigmapsf"],
            "filter": self.raw_data["filter"],
        })

    def mk_phased(self, epoch, period):
        self.data = mk_phased(self.data, epoch, period, jdnam="hjd")


class Crts(Data):
    def __init__(self, name, maglim=None):
        super().__init__()
        self.name = name
        self.survey = "crts"
        self.maglim = maglim
        self.columns = {"Mag": "mag", "Magerr": "magerr"}

    def is_raw_data_exist(self, add="", ext="csv"):
        return os.path.isfile(f"{self.raw_dir}/{self.name}-{self.survey}.{ext}")

    def is_data_exist(self, filtr="", lim="", ext="dat"):
        return super().is_data_exist(filtr, lim=lim, ext=ext)

    def read_raw_data(self, add="", ext="csv"):
        self.raw_data = pd.read_csv(f"{self.raw_dir}/{self.name}-{self.survey}{add}.{ext}")

    def prepare_data(self, ra, dec, maglim_up=None, maglim_low=None):
        self.data = pd.DataFrame({
            "mag": self.raw_data["Mag"],
            "magerr": self.raw_data["Magerr"],
        })
        self.data["hjd"] = mk_hjd_corr(self.raw_data["MJD"], ra, dec, obs="Catalina")  # obs="Siding Springs" obs="Catalina"
        if self.maglim:
            self.data = self.data[self.data["magerr"] < self.maglim]
        if maglim_up:
            self.data = self.data[self.data["mag"] > maglim_up]
        if maglim_low:
            self.data = self.data[self.data["mag"] < maglim_low]

    def read_data(self, filename):
        """Read CRTS data"""
        columns = ["hjd", "mag", "magerr"]
        self.data = pd.read_csv(filename, delim_whitespace=True, names=columns)

    def save_datafile(self, ext="dat", hjdprec=7, magprec=6):
        if len(self.data.index):
            save_datafile(self.data, f"{self.data_dir}/{self.name}-{self.survey}.{ext}", hjdprec=hjdprec, magprec=magprec)

    def mk_phased(self, epoch, period):
        self.data = mk_phased(self.data, epoch, period, jdnam="hjd")


def read_crts_data(filename):
    """Read CRTS data"""
    columns = ["hjd", "mag", "magerr"]
    return pd.read_csv(filename, delim_whitespace=True, names=columns)


def read_gds_data(filename, filt, ra, dec, errlim=0.049):
    """Read GDS data"""
    columns = ["hjd", "mag", "magerr"]
    data_dir = "../data"
    data = pd.read_csv(f"{data_dir}/{filename}-g{filt}.dat", delim_whitespace=True, names=columns)
    data["hjd"] = mk_hjd_corr(data["hjd"], ra, dec, obs="Siding Springs")
    # data["hjd"] += JD_SHIFT
    data = data[data["magerr"] < errlim]
    return data


def read_ps_data(filename):
    """read data: mjd mag magerr"""
    return pd.read_csv(filename, delim_whitespace=True)


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
    elif obs == "Catalina":
        loc = CATALINA
    elif obs == "Siding Springs":
        loc = SIDINGSPRINGS
    elif obs == "Atlas-MLO":
        loc = ATLAS_MLO
    elif obs == "Atlas-CHL":
        loc = ATLAS_CHL
    elif obs == "Atlas-STH":
        loc = ATLAS_STH
    helio_time = time + time.light_travel_time(coord, "heliocentric", location=loc)
    return helio_time.mjd + JD_SHIFT


def mk_fsh(filtshift):
    """Generate text for plot labels"""
    if filtshift == 0:
        return ""
    elif float(filtshift) < 0:
        return f" +{abs(filtshift)}"
    else:
        return f" −{filtshift}"  # :+1.1f
