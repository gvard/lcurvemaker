import numpy as np
import pandas as pd

from astropy.io import ascii
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

from panstarrs import ps1cone, ps1search, addfilter


HALEAKALA = EarthLocation(lon=-156.257, lat=20.71, height=3048)
PALOMAR = EarthLocation(lon=-116.86194, lat=33.3581, height=1712)
JD_SHIFT = 2400000.5
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


def read_ps_data(filename):
    """read data: mjd mag magerr"""
    return pd.read_csv(filename, delim_whitespace=True)


def read_crts_data(filename):
    """Read CRTS data"""
    lcurve = pd.read_csv(filename, delim_whitespace=True)
    return lcurve.iloc[:, 0], lcurve.iloc[:, 1], lcurve.iloc[:, 2]


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
        "jd": (data["mjd"] + JD_SHIFT).round(6).astype("str").str.ljust(14, "0"),
        "mag": data["mag"].round(6).astype("str").str.ljust(9, "0"),
        "magerr": data["magerr"].round(6).astype("str").str.ljust(8, "0"),
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
