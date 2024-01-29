import pandas as pd
from pandas.core import frame
import sys
import subprocess

time_file = sys.argv[1]
hmmfile = sys.argv[2]
subprocess.check_call( "grep '>' {0} | cut -f 1,3 > {0}.grep".format(time_file), shell =True)
time_list = pd.read_csv("{0}".format(time_file+".grep"), header=None, sep= "\t")
hmm_list = pd.read_csv("{0}".format(hmmfile), sep= "\t")
time_list['insertion_time'] = abs(time_list[1])/100/(2*3.71*10**-9)/1000000
time_list[0] = time_list[0].str.replace('>', '')
time_list['ID '] = time_list[0].str.split('_').apply(lambda x:x[2] + '_' + x[0] + '_' + x[1])
try:
	f1 = int(time_list['ID '].str.split('_').apply(lambda x:x[0])[0])
except:
	f1 = time_list['ID '].str.split('_').apply(lambda x:x[0])[0]
if isinstance(f1, str):
	time_list['ID '] = time_list[0].str.split('_').apply(lambda x:x[1] + '_' + x[0])
hmm_list['ID '] = hmm_list['ID '].str.strip(' ')
del time_list[0]
del time_list[1]
del hmm_list['Unnamed: 34']
print(time_list)
all_list = pd.merge(hmm_list, time_list, how='left', on = 'ID ')
all_list.to_csv(hmmfile+'_time' +'.txt', sep = '\t', index=False)

