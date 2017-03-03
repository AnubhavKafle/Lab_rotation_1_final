for i in $(ls *.txt) ; do awk -F'\t' '{print ($1,$9,tolower($16)) }' $i | grep -v -E "Silent|RNA" | sed 's/ /\t/g'  >> "Filtered_"$i; done 
for i in $(ls Filtered*) ; do  cut -d"-" -f1,2,3 $i >> "filtered_"$i; done
rm Filtered*
mv filtered* ./Mu
rename 's/filtered_Filtered_/Filtered_/' filtered*
for i in $(ls *.txt) ; do  sed 's/Tumor_Sample_Barcode/Patience_ID/g' $i >> "tested_"$i; done
rm Filtered*
for i in $(ls *.txt) ; do  sed 's/-/./g' $i >> "point_"$i; done
rm tested*
rename 's/point_//g' point_*
