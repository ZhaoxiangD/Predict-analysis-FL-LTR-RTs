import pandas as pd
from pandas.core import frame
import sys

list_name = sys.argv[1]
try:
	threshhold = float(sys.argv[2])
except:
	threshhold = 0.2

sylist = pd.read_csv("{0}".format(list_name), header = None, sep = "\t")
a_chrlist = []
b_chrlist = []
for i in range(len(sylist)):
	if sylist.iloc[i,0] not in a_chrlist:
		a_chrlist.append(sylist.iloc[i,0])
	if sylist.iloc[i,3] not in b_chrlist:
		b_chrlist.append(sylist.iloc[i,3])

dataf = pd.DataFrame(columns=['a','b','num'])
i=0
for aname in a_chrlist:
	alist = sylist[sylist.iloc[:, 0] == aname]
	alen = len(alist)
	for bname in b_chrlist:
		percent = len(alist[alist.iloc[:, 3] == bname])/alen
		dataf.loc[i] = [aname, bname, percent]
		i = i+1

dataf[dataf['num'] > threshhold][["a","b"]].to_csv("{0}.correspond_list".format(list_name), sep = '\t', header = False, index = False)
