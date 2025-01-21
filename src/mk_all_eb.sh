lcdir="ebecc"
mkdir -p ../lc/$lcdir

fnamecc="g194 g523 g261 gsc02874 18abwvvvw nsv14401 g505 hd74588 m1396 cr0629 cr0640 sao25077"\
" g176 m1404 v1403cyg epic1649 cz1850 v481dra g182 cz701 cr0622 c79-0 19acdncga v1166tau bst0730 m1399"\
" nsv18979 nsv12175 m1356 m1366 cr1029 g163 bst21588 v338lac g1829 gds1549 c235-13 epic1643 g187 cr0630 c344-18"\
" 18aaaawyw gsc04267 g219 oxdel g179 c210-16 bst0731 gds0722 g460 cz1848 18aceoesf m1397 c257-37 cz1156"\
" hd180641 tyc3983 cz1153 g452 c132-17 c336-2 g217 18aabpshf 18aainzhi c88-0 g684 c111-36 g204 m1338 g220"\
" gsc0219 cz1847 g956 hd49241 18abuksln cr0849 nsv1409 g333 nsv2702 c311-3 c180-34 nsv18927 nsv18857 c256-24"\
" m1393 g434 g965 g021-005 gds280-4 m1349 nsv61 nsv19039 v1123cas gds0705 m1374 nsv18898 nsv14698 m1398"\
" nsv3760 nsv8389 nsv3749 nsv3481"

fnameccref=" kic9153621 epic2061 oblg254039 oblg138838 tulyn mgab282 oblg029297 oblg288795 oblg165652 wsp12-29 nsv14400 ewcru"
fnamkplr = " koi1008 koi436 koi3528 koi3155 koi3510 koi3511 koi3711 koi7141 koi7056 koi6782 koi7385"

for nam in $fnamecc $fnameccref $fnamkplr
do
  echo $nam $lcdir/
  python plot_merged_phased_data.py $nam $lcdir/ -loaz
  # python plot_hlsp_data.py $nam $lcdir/ -lz
done

lcdir="eb"
mkdir -p ../lc/$lcdir

fnameb="v2291oph gds0644 lwpup v698cyg epic0603 gds0852 hngem hd10220 cr3573 ztf356-58 nsv13073 v476lac cz517"\
" cz250 nsv17180 g109 oblg117530 v374cyg nsv17721 c350-33 as357-18 c65-1 17aablgmo gds0754 cz1826"\
" n171549 cz501 c349-7 17aadnujm cz257 g2613 n637-68 nsv18968 nsv5437 ptf5-6432 dde2 cz1225 cz127"\
" da240 cz2494 nsv4595 17aabgexm brh160 18abryisf c70-6 18acrvwcz hd277600 a293-29 g16ams c248-3 18aagibwi u170"\
" 18abqadhy da243 dde63 u248199 navl02 c86-5 svkv70 nsv1791 u248 dde65 c154-3 bmam440 sd1443 c66-7 z0308 c23-13"\
" u240 u222 nada99 u244-187 n2fu0 g986 z1954 sd0344 z1454 mg297 dde193 dde157 z0547 2m1358 z1927"

fnamdde=" dde58 dde62 dde60 dde203 dde64 dde78 dde74 dde80 dde176"

for nam in $fnameb $fnamdde
do
  python plot_merged_phased_data.py $nam $lcdir/ -loaz
  python plot_hlsp_data.py $nam $lcdir/ -lz
done

for mnum in 4 5 6 8 11 15 16 17 18 19 20 21 23 24
do
  python plot_merged_phased_data.py minkovskiy$mnum $lcdir/ -loaz
done
python plot_merged_phased_data.py gusev4 $lcdir/ -loa
