# Данные

## Фотометрические данные

* [ZTF (Zwicky Transient Facility, Установка для поиска транзиентов имени Цвикки): список доступных выпусков данных](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?projshort=ZTF)
* [PTF (Palomar Transient Factory, Паломарская фабрика транзиентов): форма поиска](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd)
* [Pan-STARRS1, PS1 DR2: примеры кода для запроса фотометрии](https://ps1images.stsci.edu/ps1_dr2_api.html)
* [ASAS-3 (All Sky Automated Survey)](https://www.astrouw.edu.pl/asas/?page=catalogues)
* [ASAS-SN (All-Sky Automated Survey for Supernovae)](https://asas-sn.osu.edu/)
* [ATLAS (Asteroid Terrestrial-impact Last Alert System)](https://fallingstar.com/)
* [The Catalina Surveys: CRTS (Catalina Real-time Transient Survey), форма поиска](http://nunuku.caltech.edu/cgi-bin/getcssconedb_priv.cgi)

* [Bochum Survey of the Southern Galactic Disk, GDS (Galactic Disk Survey)](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/AN/336/590),
  [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/AN/336/590): в поле ID ввести обозначение объекта, начинающееся с "GDS_J". В результирующей выдаче кликнуть по полю LC (Show the light curves), затем кликнуть по Data as a Table. Фотометрия в фильтрах r_s и i_s (при ее наличии) разделена пустой строкой. Для получения фотометрии программным способом необходимо сохранить файл [varlc.dat](https://cdsarc.cds.unistra.fr/ftp/J/AN/336/590/varlc.dat.gz) (размер файла 74 Мб). Обозначения объектов здесь содержатся в явном виде.

* [SuperWASP (Wide Angle Search for Planets)](https://en.wikipedia.org/wiki/Wide_Angle_Search_for_Planets), [форма поиска](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblSearch/nph-tblSearchInit?app=ExoTbls&config=superwasptimeseries)
* [OGLE (Optical Gravitational Lensing Experiment) catalogue of variable stars, форма поиска](https://ogledb.astrouw.edu.pl/~ogle/OCVS/catalog_query.php)
* [MASCARA (Multi-site All-Sky CAmeRA) variable star catalogue](https://home.strw.leidenuniv.nl/~burggraaff/MASCARA_variables/)
* [SDSS (Sloan Digital Sky Survey) DR19, форма поиска](https://skyserver.sdss.org/dr19/en/tools/search/radial.aspx)
* [SkyMapper Southern Sky Survey, форма поиска](https://skymapper.anu.edu.au/cone-search/)

* [Hipparcos](https://ui.adsabs.harvard.edu/abs/1997ESASP1200.....E),
  [VizieR](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/239/hip_main):
  в поле HIP ввести идентификатор объекта из каталога. В настройках (Preferences) слева установить флажок All columns, либо кликнуть по единице в первом столбце Full.
  В поле HIPep кликнуть по ссылке, затем кликнуть по Data as a Table. Данные представлены в виде таблицы из трех столбцов, разделенных пробелами: дата, зв. величина и ее погрешность. Для получения фотометрии программным способом можно скачать весь архив фотометрических данных миссии.  
  Также в поле SpType содержится спектральный класс.

* [Gaia DR3](https://ui.adsabs.harvard.edu/abs/2022yCat.1355....0G),
  [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/epphot):
  в поле Source ввести идентификатор из каталога. При наличии уведомления Result truncated to 50 rows в настройках (Preferences) слева увеличить максимальное количество строк (селектор max). Там же можно выбрать формат вывода.  
  Данные можно получить программно, [используя Astroquery](https://astroquery.readthedocs.io/en/latest/vizier/vizier.html) - пакет Python, координируемый проектом Astropy.

* [TESS, форма поиска портала MAST (Barbara A. Mikulski Archive for Space Telescopes)](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html):
  В поле поиска ввести координаты или имя объекта, слева выбрать необходимый тип данных. Например, Mission: HLSP (High Level Science Products), Provenance Name: [QLP](https://archive.stsci.edu/hlsp/qlp) (Quick‑Look Pipeline), TGLC (TESS‑Gaia Light Curve), TESS‑SPOC (TESS Light Curves From Full Frame Images, SPOC - Science Processing Operations Center). При необходимости уменьшить область поиска: ограничить расстояние в угловых секундах (Distance, arcsec).

* [CoRoT, VizieR](https://cdsarc.cds.unistra.fr/viz-bin/cat/B/corot). Данные поделены на 2 части: объекты, которые наблюдались в режиме ярких и тусклых звезд. В последнем случае кривые блеска доступны через VizieR в трех вариантах: Raw level, BARFILL level и Corrected level.
* [Kepler, форма поиска портала MAST](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)
* [K2, форма поиска портала MAST](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)

## Спектральный класс, эффективная температура и ускорение силы тяжести (log g)

* [Gaia DR3](https://ui.adsabs.harvard.edu/abs/2022yCat.1355....0G),
  [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/paramp):
  в поле Source ввести идентификатор из каталога Gaia DR3. В настройках (Preferences) слева установить флажок All columns, либо кликнуть по единице в первом столбце Full. Спектральный класс содержится в поле SpType-ELS, эффективная температура - в поле Teff, ускорение силы тяжести - в поле logg.

* [TESS Input Catalog (TIC) version 8.2](https://ui.adsabs.harvard.edu/abs/2021arXiv210804778P),
  [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=IV/39/tic82):
  в поле TIC ввести идентификатор TIC. Эффективная температура содержится в поле Teff, ускорение силы тяжести - в поле logg.  
  Следует помнить, что эффективная температура, когда это возможно, берется из списка спектроскопических каталогов и в порядке предпочтения, [указанном в таблице 1 данной статьи](https://ui.adsabs.harvard.edu/abs/2019AJ....158..138S).
  Однако значения log g не берутся из спектроскопических каталогов, а рассчитываются с использованием значений массы и радиуса, так что их следует использовать с осторожностью и лишь как грубую оценку.


## Список литературы

* [Masci, F. J.; et al., 2019, The Zwicky Transient Facility: Data Processing, Products, and Archive](https://ui.adsabs.harvard.edu/abs/2019PASP..131a8003M)
* [Chambers, K. C.; et al., 2016, The Pan-STARRS1 Surveys](https://ui.adsabs.harvard.edu/abs/2016arXiv161205560C)
* [Pojmanski, G., 2002, The All Sky Automated Survey](https://ui.adsabs.harvard.edu/abs/2002AcA....52..397P)
* [Kochanek, C. S.; et al., 2017, The All-Sky Automated Survey for Supernovae (ASAS-SN) Light Curve Server v1.0](https://ui.adsabs.harvard.edu/abs/2017PASP..129j4502K)
* [Tonry, J. L.; et al., 2018, ATLAS: A High-cadence All-sky Survey System](https://ui.adsabs.harvard.edu/abs/2018PASP..130f4505T)
* [Drake, A. J.; et al., 2009, Catalina Real-time Transient Survey](https://ui.adsabs.harvard.edu/abs/2009ApJ...696..870D)
* [Hackstein, M.; et al., 2015, The Bochum Survey of the Southern Galactic Disk: II. Follow-up measurements and multi-filter photometry for 1323 square degrees monitored in 2010 - 2015](https://ui.adsabs.harvard.edu/abs/2015AN....336..590H)
* [Butters, O. W.; et al., 2010, The first WASP public data release](https://ui.adsabs.harvard.edu/abs/2010A%26A...520L..10B)
* [Udalski, A.; et al., 2008, The Optical Gravitational Lensing Experiment. Final Reductions of the OGLE-III Data](https://ui.adsabs.harvard.edu/abs/2008AcA....58...69U)
* [Udalski, A.; et al., 2015, OGLE-IV: Fourth Phase of the Optical Gravitational Lensing Experiment](https://ui.adsabs.harvard.edu/abs/2015AcA....65....1U)
* [Burggraaff, O.; et al., 2018, Studying bright variable stars with the Multi-site All-Sky CAmeRA (MASCARA)](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..32B)

* [Perryman, M. A. C.; et al., 1997, The HIPPARCOS and TYCHO catalogues](https://ui.adsabs.harvard.edu/abs/1997HIP...C......0E)
* [Gaia collaboration; et al., 2022, Gaia Data Release 3 (Gaia DR3) Part 1 Main source](https://ui.adsabs.harvard.edu/abs/2022yCat.1355....0G)
* [Huang, C. X.; et al., 2020, Photometry of 10 Million Stars from the First Two Years of TESS Full Frame Images: Part I](https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..204H)
* [Han, T.; Brandt, T. D., 2023, TESS-Gaia Light Curve: A PSF-based TESS FFI Light-curve Product](https://ui.adsabs.harvard.edu/abs/2023AJ....165...71H)
* [Solano, E.; et al., 2009, The LAEX and NASA portals for CoRoT public data](https://ui.adsabs.harvard.edu/abs/2009A&A...506..455S)
* [Howell, S. B.; et al., 2014, The K2 Mission: Characterization and Early Results](https://ui.adsabs.harvard.edu/abs/2014PASP..126..398H)

* [Creevey, O. L.; et al., 2023, Gaia Data Release 3: Astrophysical parameters inference system (Apsis) I -- methods and content overview](https://ui.adsabs.harvard.edu/abs/2023A&A...674A..26C)
* [Stassun, K. G.; et al., 2019, The Revised TESS Input Catalog and Candidate Target List](https://ui.adsabs.harvard.edu/abs/2019AJ....158..138S)


См. также [Frequently Asked Questions сайта AAVSO VSX](https://vsx.aavso.org/index.php?view=about.faq), раздел References.