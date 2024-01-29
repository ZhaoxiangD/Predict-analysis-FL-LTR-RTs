cat *.xls > all.xls
awk -F'\t' '{print $1,$2,$3,$5,$6,$14}' all.xls > all.xls_t
sed '/ID/d' all.xls_t >> all.xls_t_1
sed '/no/d' all.xls_t_1 >> all.xls_t_2
python3 ~/Wuhaotian/LTR_showing.py -i all.xls_t_2 -o all.xls_t_2_new
awk '{if(/ID/){print "ID-chr\tarm_start\tarm_end\tstrand"}else{print $1"-"$2"\t"$3"\t"$4"\t"$NF;print $1"-"$2"\t"$5"\t"$6"\t"$NF}}' all.xls_t_2_new >all.xls_t_2_new1
perl -pe  "s/_N/\tN/g" all.xls_t_2_new1 > all.xls_t_2_new2
sed 'y/-c/\tc/' all.xls_t_2_new2 >all.xls_t_2_new3
cat *.fna > all.fna
perl /home/yinglu/pl/liyalin/tools/getFragmentFromASeqBySite_batch.pl all.fna all.xls_t_2_new3 0 1 3 4 5 >z 2>zz
awk '{if(/>/){if($NF~/+/){$NF="p"}else{$NF="m"}print ">"$3"_"$2"_"$NF}else{print}}'  z >zz
awk '{if(/>/){if(dict[$1]=="y"){print $1"_"2}else{print $1"_"1;dict[$1]="y"}}else{print}}' zz >zzz
perl /home/yinglu/pl/linux_fasta2line.pl zzz

