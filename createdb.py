import pandas as pd
import subprocess
import sys

file_name = sys.argv[1]
sub_name = sys.argv[2]

hmm_list = pd.read_csv("{0}".format(file_name), sep="\t")
hmm_list = hmm_list[hmm_list['subfamily '].str.strip(" ")== sub_name]
hmm_list.reset_index(drop=True, inplace=True)
try:
	subprocess.check_call("rm {0}.db.aa".format(sub_name), shell = True)
except:
	a = 'a'
for index, dom in hmm_list.iterrows():
	subprocess.check_call("cat {0} >> {1}.db.aa".format(dom["proteinID "].strip(" ")+"*", sub_name), shell = True)
subprocess.check_call("diamond makedb --in {0}.db.aa -d {0}.db.aa.dia".format(sub_name), shell = True)
