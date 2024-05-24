mkdir -p ../lc/eb/
for mnum in 4 5 6 8 11 15 16 17 18 19 20 21 23 24
do
  python plot_merged_phased_data.py minkovskiy$mnum eb/ -l
done
python plot_merged_phased_data.py gusev4 eb/ -l
fnams="18abwvvvw cz1850 c79-0 19acdncga m1399 nsv18979 nsv13073 nsv12175 m1356 m1366"\
" cz517 c344-18 18aaaawyw c210-16 c257-37 c111-36 m1338 nsv1409 c180-34 nsv18927"\
" dde2 nsv61 18acrvwcz a293-29 g16ams 18aagibwi u170 m1398 u248199 navl02 svkv70"\
" u248 dde65 u240 u222 n2fu0 dde58 dde60"
for nam in $fnams
do
  python plot_merged_phased_data.py $nam eb/ -l
done
