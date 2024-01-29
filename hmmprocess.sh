set -e
set -u
set -o pipefail
fnaname=$(find . -name "c*.fna" | xargs basename -s '.fna' | sort -V)
fnaname=(${fnaname//\n/ })
head -1 "./${fnaname[0]}_result.xls" > all_result.xls
for file in ${fnaname[@]}
do
  tail -n+2 ./${file}_result.xls >> all_result.xls
done
python3 ../script/hmmstep.py > hmmresult.xls
sed -i 's///g' hmmresult.xls
