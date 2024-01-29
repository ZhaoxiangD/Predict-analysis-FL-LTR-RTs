import pandas as pd
from pandas.core import frame
import sys
import subprocess

hmmfile = 'hmmresult.xls_time.txt'
#family = 'LTR/Gypsy '
#subfamily = sys.argv[1]
#subfamily = str(subfamily) + " "
hmm_list = pd.read_csv("{0}".format(hmmfile), sep= "\t")
labels, uniques = pd.factorize(hmm_list['subfamily '])
sub_list = uniques.tolist()

for subfamily in sub_list:
    hmm_list = pd.read_csv("{0}".format(hmmfile), sep= "\t")
    hmm_list = hmm_list.loc[hmm_list['subfamily '] == subfamily]
    hmm_list = hmm_list.loc[hmm_list['env_stop '] - hmm_list['env_start '] >= 200]
    for i, rows in hmm_list.iterrows():
        file_name = rows['proteinID '] + '*' + '.aa'
        file_name = file_name.replace(' ','')
        p = subprocess.Popen("grep -v '>' {0}".format(file_name), shell=True, stdout=subprocess.PIPE)
        rt_seq = p.stdout.read()
        rt_seq = rt_seq.decode()
        rt_seq = rt_seq.replace('\n', '')
        name = '>' + rows['proteinID '] + '_' + str(rows['insertion_time'])
        name = name.replace(' ','')
        print(name)
        subprocess.check_call("echo '{0}' >> {1}.time.aa && echo {2} >> {1}.time.aa".format(name, subfamily.strip(' '), rt_seq), shell = True)

print(hmm_list)
