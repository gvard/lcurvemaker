lcdir="cv"
mkdir -p ../lc/$lcdir

fnamcv="brodskaya2 davv4 lanat2 davv16 shcheglov1 dizepa2 khrapov1 davv18 exuzyan2 exuzyan1 davv22 rybka1 bhanavira6 rsmr1 mem1 khrapov2"\
" mazepa2 davkld1 davv5 lisniak6 davv2 karachurin16 lisniak5 karachurin14 shcheglov3 bhanavira4 bhanavira20 karachurin12 lisniak2 lisniak4"\
" minkovskiy7 minkovskiy3 minkovskiy10 bhanavira12 druzhbina1 bhanavira19 minkovskiy14 minkovskiy22 bhanavira13 minkovskiy13 bhanavira7"\
" lisniak3 bhanavira9 bhanavira17 minkovskiy12 bhanavira1 kazakevichoa2 kazakevichoa4 kazakevichoa1 kazakevichoa3 kld1 druzhbina2 minkovskiy2"\
" minkovskiy1 bhanavira5 exuzyan4 taya2 koveshnikov1 apa-v1 mazepa1 lanat1 dizepa1 rybka2 davv17 karachurin15 exuzyan3"

for nam in $fnamcv
do
  echo $nam $lcdir/
  python plot_merged_phased_data.py $nam $lcdir/ -o
done
