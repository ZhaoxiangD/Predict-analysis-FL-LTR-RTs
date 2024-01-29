set -e
set -u
set -o pipefail
fnaname=$(find . -name "c*.xls" | xargs basename -s '_result.xls' | sort -V)
fnaname=(${fnaname//\n/ })
head -1 "./${fnaname[0]}_result.xls" > all_result.xls
for file in ${fnaname[@]}
do
  tail -n+2 ./${file}_result.xls >> all_result.xls
done
python ~/Long_terminal_repeated_species_list/script/hmmstep.py > hmmresult.xls
sed -i 's///g' hmmresult.xls
