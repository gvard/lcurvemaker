import json

HTML_HEAD = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Исследования переменных звезд в Московском дворце пионеров</title>
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
.sne tr td:nth-of-type(2), .sne tr td:nth-of-type(3) {white-space: nowrap;}
.merged tr td:nth-of-type(3), .merged tr td:nth-of-type(4) {white-space: nowrap;}
.generated tr td:nth-of-type(3), .generated tr td:nth-of-type(4) {white-space: nowrap;}
.pz tr td:nth-of-type(2) {max-width: 119px; white-space: nowrap;}
.pz tr td:nth-of-type(3) {max-width: 146px; overflow-wrap: break-word;}
.pz tr td:nth-of-type(4) {min-width: 81px;}
.pz tr td:nth-of-type(5) {min-width: 81px;}
.pz tr td:nth-of-type(8) {max-width: 50px; overflow-wrap: break-word;}
.pz tr td:nth-of-type(10) {max-width: 90px; overflow-wrap: break-word;}
.pz tr td:nth-of-type(12) {max-width: 39px; overflow-wrap: break-word;}
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
<tr><th>#</th><th>Name</th><th>Other</th><th>RA (J2000)</th><th>DEC</th><th>Type</th><th>Max</th><th>Min</th><th>Sys</th><th>Period</th><th>Epoch (JD)</th><th>D, %</th><th>Sp</th><th>Comment</th><th>L.Curve</th><th>Find.Chart</th></tr>"""


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


objs = []
for mnum in (4, 5, 6, 8, 11, 15, 16, 17, 18, 19, 20, 21, 23, 24):
    objs.append(f"minkovskiy{mnum}")
objs.append("gusev4")

with open("../index.html", "w", encoding="utf8") as htmlfile:
    print(HTML_HEAD, file=htmlfile)
    i = 0
    for obj in objs:
        i += 1
        with open(f"../objects/{obj}.json", encoding="utf-8") as json_data:
            obj = json.loads(json_data.read())
        nam = str(next(iter(obj)))
        obj = obj[nam]
        obj_lnk = f'<a href="https://www.aavso.org/vsx/index.php?view=detail.top&amp;oid={obj.get("vsxoid")}" target="_blank">{nam}</a>'
        ra, dec = obj.get("coord")[:12].strip(), obj.get("coord")[12:].strip()
        coord_ra = f'<a href="https://simbad.u-strasbg.fr/simbad/sim-coo?Coord={ra.replace(" ", "+")}+{dec.replace("+", "%2B").replace(" ", "+")}%09&CooEpoch=2000&CooEqui=2000&Radius=1&Radius.unit=arcmin&submit=submit+query" target="_blank">{ra}</a>'
        coord_dec = f"<details><summary>{dec}</summary>{hms_to_deg(ra=ra)} {hms_to_deg(dec=dec)}</details>"
        fclnk = obj.get("fclnk")
        lc = obj.get("lc")
        if isinstance(lc, list):
            lc = lc[0]
        lclink = f"https://github.com/gvard/lcurvemaker/raw/main/lc/{lc}" if lc else ""
        if isinstance(fclnk, list):
            fclnk = fclnk[0]
        try:
            min = obj.get("min").replace("<", "&lt;")
        except AttributeError:
            min = obj.get("min")
        tr_content = [
            f"<td>{x}</td>"
            for x in (
                i,
                obj_lnk,
                obj.get("other"),
                coord_ra,
                coord_dec,
                obj.get("type"),
                obj.get("max"),
                min,
                obj.get("system"),
                obj.get("period"),
                obj.get("epoch"),
                int(obj.get("d") * 100),
                obj.get("sp"),
                obj.get("comm"),
                f'<a href="{lclink}" target="_blank">{lc}</a>' if lc else "",
                f'<a href="{fclnk}" target="_blank">{nam.lower().replace(" ", "")}-fc</a>',
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
