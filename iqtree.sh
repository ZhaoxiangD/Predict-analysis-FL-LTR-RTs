set -e
set -u
set -o pipefail
name=$1
cat *_${name}.aa > all_${name}.aa
mafft all_${name}.aa > all_${name}.aa.aln
sed -i 's/:/_/g' all_${name}.aa.aln
/home/software/trimal-trimAl/source/trimal -in all_${name}.aa.aln -out all_${name}.aa.aln.trimed -automated1
iqtree -s all_${name}.aa.aln.trimed -m MFP -bb 1000 -bnni -redo > ${name}.log
