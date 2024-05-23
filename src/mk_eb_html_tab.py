"""Python script for generating an HTML table from the objects data
"""

import json

from astropy.coordinates import SkyCoord, get_constellation


objs = []
for mnum in (4, 5, 6, 8, 11, 15, 16, 17, 18, 19, 20, 21, 23, 24):
    objs.append(f"minkovskiy{mnum}")
objs.extend((
    "gusev4", "18abwvvvw", "cz1850", "c79-0", "19acdncga", "m1399", "nsv18979",
    "nsv13073", "nsv12175", "m1356", "m1366", "cz517", "c344-18", "18aaaawyw",
    "c210-16",
    "c257-37",
    "c111-36", "m1338", "nsv1409", "c180-34", "nsv18927",
    "dde2", "nsv61", "18acrvwcz", "a293-29", "g16ams", "18aagibwi",
    "u170", "m1398", "u248199", "navl02", "svkv70", "u248",
    "dde65", "u240", "u222", "n2fu0", "dde58", "dde60",
    ))

# TAB_HEAD = "<tr><th>#</th><th>Name</th><th>Other</th><th>RA (J2000)</th><th>DEC</th><th>Type</th><th>Max</th><th>Min</th><th>Sys</th><th>Period</th><th>Epoch (JD)</th><th>D, %</th><th>Sp</th><th>Comment</th><th>L.Curve</th><th>Find.Chart</th></tr>"
TAB_HEAD = "<tr><th>#</th><th>Name</th><th>Other</th><th>RA (J2000)</th><th>DEC</th><th>Con</th><th>Type</th><th>Max</th><th>Min</th><th>Amp</th><th>Sys</th><th>Period</th><th>D %</th><th>Prob</th><th>min II</th><th>Sp</th><th>Teff</th><th>logg</th><th>r</th><th>r ph</th><th>Exc</th><th>Gaia</th><th>L.Curve</th><th>Comment</th></tr>"
HTML_HEAD = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Eclipsing binaries data</title>
<script src="https://www.kryogenix.org/code/browser/sorttable/sorttable.js"></script>
<style>
body {
  background-color: #1b1f22;
  color: #dee2e6;
  font-family: Arial, Helvetica, sans-serif;
  font-size: 0.85em;
  margin: 0;
}
table, td {
  border-collapse: collapse;
  border: 1px solid #454b52;
  padding: 3px;
}
th {
  border: 2px solid #454b52;
}
.pz tr td:nth-of-type(2) {max-width: 119px;}
.pz tr td:nth-of-type(3) {max-width: 133px; overflow-wrap: break-word;}
.pz tr td:nth-of-type(4) {min-width: 81px;}
.pz tr td:nth-of-type(5) {min-width: 81px;}
.pz tr td:nth-of-type(8) {max-width: 50px; overflow-wrap: break-word;}
.pz tr td:nth-of-type(10) {max-width: 90px; overflow-wrap: break-word;}
.pz tr td:nth-of-type(12) {max-width: 100px; overflow-wrap: break-word;}
.pz tr td:nth-of-type(15) {max-width: 140px; overflow-wrap: break-word;}
table.sortable thead {
  background-color: inherit;
}
a {text-decoration: none; color: #0d6efd}
a:hover {text-decoration: underline;}
details > summary {
  list-style: none;
  cursor: pointer;
}
@media screen and (prefers-color-scheme: light) {
  body {background-color: white; color: black}
}
</style>
</head>
<body>
<table class="sortable pz">
"""


def hms_to_deg(ra="", dec=""):
    """Converting Between Decimal Degrees and hours, minutes, seconds.
    Snippet from http://www.bdnyc.org/2012/10/decimal-deg-to-hms/"""
    RA, DEC, rs, ds = "", "", 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == "-":
            ds, D = -1, abs(D)
        deg = D + (M / 60) + (S / 3600)
        DEC = round(deg * ds, 6)
    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == "-":
            rs, H = -1, abs(H)
        deg = (H * 15) + (M / 4) + (S / 240)
        RA = round(deg * rs, 5)
    if ra and dec:
        return RA, DEC
    return RA or DEC


k_prob = 0.8
vizier_gaia_lnk = '<a href="https://vizier.cds.unistra.fr/viz-bin/VizieR-6?-out.form=%2bH%2bm&-source=I/355/gaiadr3*&Source='
with open("../index.html", "w", encoding="utf8") as htmlfile:
    print(HTML_HEAD + TAB_HEAD, file=htmlfile)
    i = 0
    for obj in objs:
        i += 1
        with open(f"../objects/{obj}.json", encoding="utf-8") as json_data:
            data = json.loads(json_data.read())
        des = str(next(iter(data)))
        obj = data[des]
        nam = obj.get("name") if obj.get("name") else des
        obj_lnk = f'<a href="https://www.aavso.org/vsx/index.php?view=detail.top&amp;oid={obj.get("vsxoid")}" target="_blank">{nam}</a>'
        crds = obj.get("coord")
        ind = crds.replace(' ', ':', 2).find(" ")
        ra, dec = crds[:ind].strip(), crds[ind:].strip()
        coord = SkyCoord(ra=obj["coordeg"][0], dec=obj["coordeg"][1], unit="deg")
        con = get_constellation(coord, short_name=True)
        coord_ra = f'<a href="https://simbad.u-strasbg.fr/simbad/sim-coo?Coord={ra.replace(" ", "+")}+{dec.replace("+", "%2B").replace(" ", "+")}%09&CooEpoch=2000&CooEqui=2000&Radius=1&Radius.unit=arcmin&submit=submit+query" target="_blank">{ra}</a>'
        coord_dec = f"<details><summary>{dec}</summary>{hms_to_deg(ra=ra)} {hms_to_deg(dec=dec)}</details>"
        fclnk = obj.get("fclnk")
        lc = obj.get("lc")
        if isinstance(lc, list):
            lc = lc[0]
        lclink = f"https://github.com/gvard/lcurvemaker/raw/main/lc/{lc}" if lc else ""
        if obj.get("lclnk"):
            lc = "VSX plot"
            lclink = obj.get("lclnk")
        # if isinstance(fclnk, list):
        #     fclnk = fclnk[0]
        try:
            mina = float(obj.get("min").replace("<", "").replace(":", ""))
            min = obj.get("min").replace("<", "&lt;")
        except AttributeError:
            min = obj.get("min")
            mina = float(min)
        amp = obj.get("amp") if obj.get("amp") else round(mina - obj.get("max"), 2)
        d = round(obj.get("d") * 100, 2) if obj.get("d") else ""
        if d and d > 8.9:
            d = int(d)
        min2 = obj.get("min2") if obj.get("min2") else ""
        if isinstance(min2, str):
            min2 = min2.replace("<", "&lt;")
        if not obj.get("d") and obj.get("prob"):
            prob = obj["prob"]
        elif obj.get("d2"):
            prob = k_prob * (obj["d"] + obj["d2"])
        else:
            prob = k_prob * obj["d"]
        if "epphot" in obj:
            epphot = " EP✔️" if obj.get("epphot") else " EP❌"
        else:
            epphot = ""
        tr_content = [
            f"<td>{x}</td>"
            for x in (
                i,
                obj_lnk,
                obj.get("other"),
                coord_ra,
                coord_dec,
                con,
                obj.get("type"),
                obj.get("max"),
                min,
                amp,
                obj.get("system"),
                obj.get("period"),
                # obj.get("epoch"),
                d,
                round(100*prob, 1),
                min2,
                obj.get("sp") if obj.get("sp") else "",
                obj.get("teff") if obj.get("teff") else "",
                obj.get("logg") if obj.get("logg") else "",
                obj.get("r") if obj.get("r") else "",
                obj.get("rp") if obj.get("rp") else "",
                round(abs(0.5 - obj.get("2ndmin")), 3) if (obj.get("2ndmin") and obj.get("2ndmin") != 0.5) else "",
                # f'{href}{obj.get("gdr3")}" target="_blank">Gaia DR3</a>' if (obj.get("gdr3") and obj.get("epphot")) else "",
                # obj.get("gdr3") if obj.get("gdr3") else "",
                f'{vizier_gaia_lnk}{obj.get("gdr3")}" target="_blank">GDR3</a>{epphot}' if obj.get("gdr3") else "",
                f'<a href="{lclink}" target="_blank">{lc}</a>' if lc else "",
                obj.get("comm"),
                # f'<a href="{fclnk}" target="_blank">{nam.lower().replace(" ", "")}-fc</a>',
            )
        ]
        tr = f"<tr>{''.join(tr_content)}</tr>"
        print(
            tr,
            file=htmlfile,
        )
    print(
        "</table>",
        file=htmlfile,
    )
