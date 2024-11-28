lcdir="ebecc"
mkdir -p ../lc/$lcdir

fnamecc="g194 g523 g261 gsc02874 18abwvvvw nsv14401 g505 m1396 cr0629 cr0640 sao25077"\
" g176 m1404 v1403cyg cz1850 g182 cz701 cr0622 c79-0 19acdncga v1166tau bst0730 m1399 nsv18979"\
" nsv12175 m1356 m1366 cr1029 g163 v338lac g1829 gds1549 c235-13 g187 cr0630 c344-18 18aaaawyw"\
" gsc04267 g219 oxdel g179 c210-16 bst0731 gds0722 g460 cz1848 18aceoesf m1397 c257-37 cz1156"\
" hd180641 cz115 cz1153 g452 c132-17 c336-2 g217 18aabpshf 18aainzhi c88-0 g684 c111-36"\
" g204 m1338 g220 gsc0219 cz1847 g956 hd49241 18abuksln cr0849 nsv1409 g333 nsv2702 c311-3 c180-34"\
" nsv18927 nsv18857 c256-24 m1393 g434 g965 g021-005 gds280-4 m1349 nsv61 nsv19039 gds0705 m1374"\
" nsv18898 nsv14698 m1398 nsv3760 nsv8389 nsv3749 nsv3481"

fnameccref=" koi7056 tulyn mgab282 wsp12-29 v1123cas nsv14400 ewcru"

for nam in $fnamecc $fnameccref
do
  echo $nam $lcdir/
  python plot_merged_phased_data.py $nam $lcdir/ -loaz
  # python plot_hlsp_data.py $nam $lcdir/ -lz
done


lcdir="eb"
mkdir -p ../lc/$lcdir

fnameb="gds0644 v698cyg hd10220 epic0603 gds0852 hngem ztf356-58 nsv13073 cr3573 v476lac cz517"\
" cz250 nsv17180 g109 v374cyg nsv17721 c350-33 as357-18 tyc3983 c65-1 17aablgmo gds0754 cz1826"\
" n171549 cz501 c349-7 17aadnujm cz257 g2613 n637-68 nsv18968 nsv5437 ptf5-6432 dde2 cz1225 cz127"\
" da240 cz2494 nsv4595 17aabgexm brh160 18abryisf 18acrvwcz a293-29 g16ams c248-3 18aagibwi u170"\
" 18abqadhy da243 u248199 navl02 svkv70 u248 dde65 sd1443 z0308 c23-13 u240 u222 nada99 n2fu0 g986 sd0344"

fnamdde=" dde58 dde60 dde62 dde78 dde80 dde74 dde176 dde203"

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
