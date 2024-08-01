lcdir="ebecc"
mkdir -p ../lc/$lcdir
for mnum in 4 5 6 8 11 15 16 17 18 19 20 21 23 24
do
  python plot_merged_phased_data.py minkovskiy$mnum $lcdir/ -loa
done
python plot_merged_phased_data.py gusev4 $lcdir/ -loa

fnamecc="g194 g523 g261 18abwvvvw nsv14401 g505 m1396 koi7056 g176 m1404 cz1850 g182"\
" c79-0 19acdncga m1399 nsv18979 nsv12175 m1356 m1366 cr1029 g163 v338lac g1829 c235-13"\
" g187 c344-18 18aaaawyw g219 oxdel g179 c210-16 g460 cz1848 18aceoesf m1397 c257-37"\
" cz1156 cz115 cz1153 g452 c132-17 c336-2 g217 18aabpshf 18aainzhi c88-0 g684 c111-36"\
" g204 m1338 g220 gsc0219 cz1847 g956 18abuksln nsv1409 g333 nsv2702 c311-3 c180-34"\
" nsv18927 nsv18857 c256-24 m1393 g434 g021-005 gds280-4 m1349 nsv61 nsv19039 m1374"\
" nsv18898 nsv14698 m1398 nsv3760 nsv8389 nsv3749 nsv3481"

for nam in $fnamecc
do
  python plot_merged_phased_data.py $nam $lcdir/ -loa
done

lcdir="eb"
mkdir -p ../lc/$lcdir
fnameb="v698cyg ztf356-58 nsv13073 v476lac cz517 cz250 g109 v374cyg c350-33 as357-18"\
" c65-1 cz1826 n171549 c349-7 17aadnujm cz501 g2613 n637-68 nsv18968 nsv5437 ptf5-6432 dde2"\
" cz1225 da240 cz2494 nsv4595 17aabgexm brh160 18abryisf 18acrvwcz a293-29 g16ams 18aagibwi"\
" u170 18abqadhy da243 u248199 navl02 svkv70 u248 dde65 c23-13 u240 u222 n2fu0"
for nam in $fnameb
do
  python plot_merged_phased_data.py $nam $lcdir/ -loa
done

fnamdde="dde58 dde60 dde62 dde78 dde80 dde74 dde176 dde203"
for nam in $fnamdde
do
  python plot_merged_phased_data.py $nam $lcdir/ -loa
done
