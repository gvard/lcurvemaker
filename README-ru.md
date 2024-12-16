# lcurvemaker

Python код для работы с кривыми блеска переменных звезд

[![en](https://img.shields.io/badge/lang-en-red.svg)](README.md)
[![ru](https://img.shields.io/badge/lang-ru-green.svg)](README-ru.md)

## Зависимости

* [Matplotlib](https://matplotlib.org)
* [pandas](https://pandas.pydata.org)
* [Astropy](https://www.astropy.org)

## Установка зависимостей

```bash
pip install matplotlib pandas astropy requests
```

или

```bash
pip install -r requirements.txt
```

## Использование

`python plot_merged_phase_data.py [-v] [-l] [-s] [-c RA DEC] [-p PERIOD] [-e EPOCH] [-r MIN MAX] [-o] [-z] [-m] [-t PLOT] nickname [savedir]`

### Аргументы

`nickname` - обозначение объекта без пробелов и специальных символов, опционально с именем каталога.
Используется для поиска файлов и присвоения имен результирующим файлам.

`savedir` устанавливает директорию для сохранения иллюстраций (опционально).

Скрипт ищет файл настроек `nickname.json` в директории `objects`.

### Опции

* `-h, --help` показывает краткую справку
* `-v, --verbose` выводит больше информации
* `-l, --lines` нарисовать линии на кривых блеска, отмечающие эпоху, максимальное значение блеска,
  длительность затмения, фазы главного и вторичного (при наличии) минимумов
* `-s, --show` показать интерактивный график вместо сохранения изображения
* `-c RA DEC, --coord RA DEC` указать координаты объекта в градусах
* `-p PERIOD, --period PERIOD` указать период в днях для фазовой кривой блеска
* `-e EPOCH, --epoch EPOCH` указать эпоху в [HJD](https://en.wikipedia.org/wiki/Heliocentric_Julian_Day)
* `-r MIN MAX, --ztfran MIN MAX` удалить все данные ZTF вне указанного диапазона
* `-o, --localps` использовать локальные данные PS1 вместо их запроса через API
* `-z, --zoom` использовать настройки для увеличенного фрагмента иллюстрации
* `-m, --model` нарисовать простую модель кривой блеска по данным настроек
* `-t, --plot` выбрать, какие данные использовать для построения. Возвможные варианты: zt ps as at cs ga og gd

### Примеры

```bash
python plot_merged_phase_data.py gusev4
python plot_merged_phase_data.py minkovskiy24 -l
```

Результат:

* [Gusev 4](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227045) [settings](objects/gusev4.json), [phase plot](lc/gusev4-ps1-ztf-atlas-poss-ph.png), [light curve](lc/gusev4-ps1-ztf-atlas.png)
* [Minkovskiy 24](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2387050) [settings](objects/minkovskiy24.json), [phase plot](lc/minkovskiy24-ps1-ztf-ph.png), [light curve](lc/minkovskiy24-ps1-ztf.png)

Для работы с несколькими объектами можно использовать скрипты. Shell:

```bash
for nam in gusev4 minkovskiy17 minkovskiy24; do python plot_merged_phased_data.py $nam; done
```

PowerShell:

```powershell
ForEach ($nam in "gusev4", "minkovskiy17", "minkovskiy24") { python .\plot_merged_phased_data.py $nam }
```

## Файлы настроек

Настройки объекта находятся в файле [JSON](https://en.wikipedia.org/wiki/JSON).
[sample.json](objects/sample.json) содержит большую часть возможных настроек.

Основные настройки для объекта с фотометрическими данными

* [ZTF](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?projshort=ZTF),
* [PS1](https://ps1images.stsci.edu/ps1_dr2_api.html),
* [ASAS-SN](https://asas-sn.osu.edu/),
* [ASAS-3](https://www.astrouw.edu.pl/asas/?page=catalogues),
* [ATLAS](https://fallingstar.com/),
* [CRTS](http://nunuku.caltech.edu/cgi-bin/getcssconedb_priv.cgi),
* [Gaia DR3](https://ui.adsabs.harvard.edu/abs/2022yCat.1355....0G),
* [GDS](https://ui.adsabs.harvard.edu/abs/2015AN....336..590H),
* [SuperWASP](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblSearch/nph-tblSearchInit?app=ExoTbls&config=superwasptimeseries),
* [OGLE](https://ogledb.astrouw.edu.pl/~ogle/OCVS/catalog_query.php),
* [SDSS](https://skyserver.sdss.org/dr18/en/tools/search/radial.aspx),
* [SkyMapper Southern Sky Survey](https://skymapper.anu.edu.au/cone-search/),
* [Hipparcos](https://ui.adsabs.harvard.edu/abs/1997ESASP1200.....E),
* [TESS](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html),
* [CoRoT](https://cdsarc.cds.unistra.fr/viz-bin/cat/B/corot),
* [Kepler](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html),
* [K2](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)

```json
{
"sample-nickname": {
  "name": "ZTF19acdncga",
  "other": "USNO-B1.0 1462-0437198, 2MASS J22470763+5617523, GSC2.3 N1CQ180004",
  "coord": "22 47 07.629 +56 17 52.26",
  "coordeg": [341.7817888, 56.2978492],
  "max": 14.35,
  "min": 14.93,
  "min2": 14.78,
  "system": "g",
  "period": 24.8435,
  "epoch": 2460161.091,
  "2ndmin": 0.399,
  "d": 0.016,
  "d2": 0.019,
  "ztfran": [15.0, 14.2], "ztflim": 0.05,
  "atlasfnam": "sample-nickname-atlas.txt", "atlaslim": {"o": 0.013, "c": 0.017},
  "asasfnam": "sample-nickname-asas.csv", "asaslim": {"V": 0.05, "g": 0.06},
  "pslim": 0.1, "pslocal": false,
  "gaiadr3fnam": "sample-nickname-gdr3.dat",
  "crtsfnam": "sample-nickname-crts.csv", "crtslim": 0.06,
  "oglefnam": "sample-nickname-ogle.dat",
  "sdssfnam": "sample-nickname-sdss.dat",
  "zcatf": "ri",
  "curveshift": true,
  "clrshift": {
    "g": 0.3, "psg": 0.32, "r": 0, "psr": 0.01, "i": -0.1, "psi": -0.1, "I": 0.1,
    "psz": -0.2, "psy": -0.3, "o": -0.05, "c": 0.1, "V": 0.3, "asasg": 0.7,
    "G": 0.05, "CV": -0.1, "gr": 0.02, "gi": -0.1
  },
  "plot": {
    "ztffilt": "gri", "psfilt": "grizy", "atlasfilt": "oc", "asasfilt": "gV", "gdsfilt": "r",
    "atlasms": {"o": 2.8, "c": 2.8}, "atlaselw": 0.26, "atledgclr": "darkred",
    "asasms": {"V": 3.5, "g": 3.5}, "asaselw": 0.26,
    "zms": 4, "gms": 4, "psms": 6, "crtsms": 5,
    "xmal": 500, "xmil": 100, "ymal": 0.1, "ymil": 0.01, "leg": "lower left",
    "xlim": [-0.5, 1.0], "xedges": 90, "xlima": [56963, 60331],
    "ylim": [14.7, 13.94], "ylima": [14.88, 13.8]
  }
}}
```

## Литературные ссылки к файлам в директории `data`

* [Masci, F. J.; et al., 2019, The Zwicky Transient Facility: Data Processing, Products, and Archive](https://ui.adsabs.harvard.edu/abs/2019PASP..131a8003M), [ZTF DR22 query form](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd?catalog=ztf_objects_dr22)
* [Chambers, K. C.; et al., 2016, The Pan-STARRS1 Surveys](https://ui.adsabs.harvard.edu/abs/2016arXiv161205560C)
* [Kochanek, C. S.; et al., 2017, The All-Sky Automated Survey for Supernovae (ASAS-SN) Light Curve Server v1.0](https://ui.adsabs.harvard.edu/abs/2017PASP..129j4502K)
* [Pojmanski, G., 2002, The All Sky Automated Survey](https://ui.adsabs.harvard.edu/abs/2002AcA....52..397P)
* [Tonry, J. L.; et al., 2018, ATLAS: A High-cadence All-sky Survey System](https://ui.adsabs.harvard.edu/abs/2018PASP..130f4505T)
* [Drake, A. J.; et al., 2009, Catalina Real-time Transient Survey](https://ui.adsabs.harvard.edu/abs/2009ApJ...696..870D)
* [Gaia collaboration; et al., 2022, Gaia Data Release 3 (Gaia DR3) Part 1 Main source](https://ui.adsabs.harvard.edu/abs/2022yCat.1355....0G)
* [Gaia Photometric Science Alerts](http://gsaweb.ast.cam.ac.uk/alerts)
* [Hackstein, M.; et al., 2015, The Bochum Survey of the Southern Galactic Disk: II. Follow-up measurements and multi-filter photometry for 1323 square degrees monitored in 2010 - 2015](https://ui.adsabs.harvard.edu/abs/2015AN....336..590H)
* [Butters, O. W.; et al., 2010, The first WASP public data release](https://ui.adsabs.harvard.edu/abs/2010A%26A...520L..10B)

## Переменные звезды, исследованные при помощи lcurvemaker: таблицы с данными и кривые блеска

* [Затменные двойные](https://gvard.github.io/variability/eb/)

## Примеры кривых блеска (представлены в карточках объектов [VSX](https://aavso.org/vsx/))

### Затменные двойные

* [Gusev 4](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227045) [settings](objects/gusev4.json), [phase plot](https://www.aavso.org/vsx_docs/2227045/5780/gusev4-phased-ps1-ztf-atlas.png)
* [Minkovskiy 16](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2225441) [settings](objects/minkovskiy16.json), [phase plot](https://www.aavso.org/vsx_docs/2225441/5780/minkovskiy16-phased-ps1-ztf-atlas.png)
* [GSC 02150-01562](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359207) [phase plot](https://www.aavso.org/vsx_docs/359207/5780/a293-29-asas-gaia-ph_3.422652.png)
* [UCAC3 170-058819](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359210) [phase plot](https://www.aavso.org/vsx_docs/359210/5780/u170-ztf-gaia-ph_3.20341.png)
* [UCAC3 240-187355](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359214) [phase plot](https://www.aavso.org/vsx_docs/359214/5780/u240-ztf-ph_0.699445.png)
* [UCAC3 248-199991](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359217) [phase plot](https://www.aavso.org/vsx_docs/359217/5780/u248-199-ztf-gaia-ps1-ph_1.976275.png)
* [UCAC3 248-205306](https://www.aavso.org/vsx/index.php?view=detail.top&oid=359219) [phase plot](https://www.aavso.org/vsx_docs/359219/5780/u248-ztf-ps1-ph_1.328946.png)
* [USNO-B1.0 1534-0126575](https://www.aavso.org/vsx/index.php?view=detail.top&oid=256288) [phase plot](https://www.aavso.org/vsx_docs/256288/5780/u1534-0126575-ztf-ps1-gaia-ph_1.853704.png)
* [USNO-B1.0 1534-0125222](https://www.aavso.org/vsx/index.php?view=detail.top&oid=256287) [phase plot](https://www.aavso.org/vsx_docs/256287/5780/u1534-222-ps1-ztf-atlas-gaia-ph_0_1.664525), [light curve](https://www.aavso.org/vsx_docs/256287/5780/u1534-222-ps1-ztf-atlas-gaia.png)
* [ZTF19acdncga](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2344600) [phase plot](https://www.aavso.org/vsx_docs/2344600/5780/19acdncga-ztf-ph_24.8435.png)

### Катаклизмические переменные

* [NSV 14686](https://www.aavso.org/vsx/index.php?view=detail.top&oid=53310) [phase plot](https://www.aavso.org/vsx_docs/53310/5780/nsv14686-ztf-atlas-asas-gaia-ph_0_1.2101554), [light curve](https://www.aavso.org/vsx_docs/53310/5780/nsv14686-ztf-atlas-asas-gaia_1.png)

### Пульсирующие переменные

* [DAV V8](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227371) [phase plot](https://www.aavso.org/vsx_docs/2227371/5780/davv8-ps1-ztf-atlas-ogle-ph.png), [light curve](https://www.aavso.org/vsx_docs/2227371/5780/davv8-ps1-ztf-atlas-ogle.png)
* [DAV V10](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227460) [phase plot](https://www.aavso.org/vsx_docs/2227438/5780/davv10-ztf-atlas-ph_378_1_2.png), [light curve](https://www.aavso.org/vsx_docs/2227438/5780/davv10-ps1-ztf-atlas_1_2_3.png)
* [DAV V11](https://www.aavso.org/vsx/index.php?view=detail.top&oid=2227460) [phase plot](https://www.aavso.org/vsx_docs/2227460/5780/davv11-ps1-ztf-atlas-gaia-ph_360.png), [light curve](https://www.aavso.org/vsx_docs/2227460/5780/davv11-ps1-ztf-atlas-gaia.png)

## Применена оптимизация изображений

* [TinyPNG: WebP, PNG, JPEG optimization](https://tinypng.com/)
* [OptiPNG](https://optipng.sourceforge.net/), см. [guide to PNG optimization](https://optipng.sourceforge.net/pngtech/optipng.html)
* [Jpegoptim](https://www.kokkonen.net/tjko/projects.html), [для Windows](https://github.com/XhmikosR/jpegoptim-windows)
* [JPEGoptim + OptiPNG + TinyPNG - image optimization (на русском)](https://open-networks.ru/d/14-jpegoptim-optipng-tinypng-optimizaciya-izobrazenii)
