import pandas as pd
from pandas.core import frame
import sys
import subprocess
import argparse

from pandas.core.reshape.concat import concat

def get_args():
    parser = argparse.ArgumentParser(description="find LTR in breakpoint")
    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', default="breakpoint.csv")
    parser.add_argument('-s', '--species_name', default= "zebrafish")
    parser.add_argument('-f', '--fasta_name', default= "all.fna")
    args = parser.parse_args()
    input_name = args.input_file
    output_name = args.output_file
    species = args.species_name
    faname = args.fasta_name
    return input_name, output_name, species, faname

input_name, output_name, species, faname = get_args()

sylist = pd.read_csv("{0}".format(input_name), header=None, sep= "\t")
a_chrlist=[] #dogfish, iloc=0
b_chrlist=[] #salmon, iloc=3
for i in range(len(sylist)):
    if sylist.iloc[i,0] not in a_chrlist:
        a_chrlist.append(sylist.iloc[i,0])
    if  sylist.iloc[i,3] not in b_chrlist:
        b_chrlist.append(sylist.iloc[i,3])

breakpointlist = pd.DataFrame(columns=['chra','chrb1', 'chrb2','break_start', 'break_stop', 'gene_start', 'gene_stop'])

for aname in a_chrlist:
    alist = sylist[sylist.iloc[:, 0] == aname]
    alist = alist.drop_duplicates([0,1,2],keep = False)
    for i in range(len(alist)-1):
        if alist.iloc[i,3] != alist.iloc[i+1, 3] and aname == alist.iloc[i+1,0]:
            breakpointlist = breakpointlist.append([{'chra': aname, 'chrb1': alist.iloc[i, 3], 'chrb2': alist.iloc[i+1,3],
            'break_start': alist.iloc[i,2], 'break_stop': alist.iloc[i+1, 1], 'gene_start': alist.iloc[i,1], 'gene_stop':alist.iloc[i,2]}], ignore_index=True)
#print(breakpointlist)

def list_process(breakpointlist, a_chrlist):
    breaklist = pd.DataFrame(columns=['chra','chrb1', 'chrb2','break_start', 'break_stop', 'reverse'])
    breaklist['reverse'] = 'False'
    num = 0
    for aname in a_chrlist:
        breaktmp = breakpointlist[breakpointlist.iloc[:,0] == aname]
        breaktmp.reset_index(drop = True, inplace = True)
        chr_3rd = True
        chralist, chrb1list, chrb2list, break_start_list, break_stop_list, break_num = [],[],[],[],[],[]
        for index, rows in breaktmp.iterrows():
            if index == 0:
                breaklist.loc[num,'chra'], breaklist.loc[num,'chrb1'], breaklist.loc[num,'break_start'], breaklist.loc[num,'break_stop'], breaklist.loc[num,'chrb2'] = rows['chra'], rows['chrb1'], rows['break_start'], rows['break_stop'], rows['chrb2']
                chralist.append(rows['chra'])
                chrb1list.append(rows['chrb1'])
                chrb2list.append(rows['chrb2'])
                break_start_list.append(rows['break_start'])
                break_stop_list.append(rows['break_stop'])
                break_num.append(num)
                num += 1
            elif index == len(breaktmp)-1 and chr_3rd:
                if rows['chrb1'] in chrb1list and rows['chrb2'] in chrb2list:
                    breaklist.loc[num-1,'break_stop'] = rows['break_stop']
                else:
                    breaklist.loc[num-1,'break_stop'] = breaktmp.iloc[index-1, 4]
                    breaklist.loc[num,'chra'], breaklist.loc[num,'chrb1'], breaklist.loc[num,'break_start'], breaklist.loc[num,'break_stop'], breaklist.loc[num,'chrb2'] = rows['chra'], rows['chrb1'], rows['break_start'], rows['break_stop'], rows['chrb2']
                    num += 1
            else:
                if rows['chrb2'] not in chrb1list and rows['chrb2'] not in chrb2list:
                    breaklist.loc[num-1,'break_stop'] = breaktmp.iloc[index-1, 4]
                    breaklist.loc[num,'chra'], breaklist.loc[num,'chrb1'], breaklist.loc[num,'break_start'], breaklist.loc[num,'break_stop'], breaklist.loc[num,'chrb2'] = rows['chra'], rows['chrb1'], rows['break_start'], rows['break_stop'], rows['chrb2']
                    chralist.append(rows['chra'])
                    chrb1list.append(rows['chrb1'])
                    chrb2list.append(rows['chrb2'])
                    break_start_list.append(rows['break_start'])
                    break_stop_list.append(rows['break_stop'])
                    break_num.append(num)
                    num += 1
                    chr_3rd = False
    return breaklist

def get_breakltr(breakpointlist,faname,speciename):
    for i, rows in breakpointlist.iterrows():
        if int(rows['break_start']) > int(rows['break_stop']):
            start = rows['break_start']
            stop = rows['break_stop']
            rows['break_start'], rows['break_stop'] = stop, start
            rows['reverse'] = 'True'
        subprocess.check_call("samtools faidx {0} {1}:{2}-{3} > {1}:{4}-{5}_{2}-{3}.fa".format(faname, rows['chra'],
        rows["break_start"], rows["break_stop"], rows["chrb1"], rows["chrb2"]), shell=True)
        subprocess.check_call("/home/software/RepeatMasker/RepeatMasker {1}:{4}-{5}_{2}-{3}.fa -pa 20 -species '{6}'".format(faname, rows['chra'],
        rows["break_start"], rows["break_stop"], rows["chrb1"], rows["chrb2"], speciename), shell=True, stdout = open('/dev/null','w'))
        print("/home/software/RepeatMasker/RepeatMasker {1}:{4}-{5}_{2}-{3}.fa -species '{6}'".format(faname, rows['chra'],
        rows["break_start"], rows["break_stop"], rows["chrb1"], rows["chrb2"], speciename))
        subprocess.check_call('grep -A 4 "LTR" {0}:{1}-{2}_{3}-{4}.fa.tbl > {0}:{1}-{2}_{3}-{4}.fa.tbl.grep'.format(rows['chra'], rows["chrb1"], rows["chrb2"], rows["break_start"], rows["break_stop"]), shell = True)
        tblname = '{0}:{1}-{2}_{3}-{4}.fa.tbl.grep'.format(rows['chra'], rows["chrb1"], rows["chrb2"], rows["break_start"], rows["break_stop"])
        ltr_tbl = read_tbl(tblname)
        breakpointlist = breakpointlist.append([{'LTR_num': ltr_tbl.loc[0,2], 'LTR_len': ltr_tbl.loc[0,3], 'LTR_percent':ltr_tbl.loc[0,5]
        , 'bel_num': ltr_tbl.loc[1,1], 'bel_len': ltr_tbl.loc[1,2], 'bel_percent': ltr_tbl.loc[1,4], 'copia_num':ltr_tbl.loc[2,1],
        'copia_len':ltr_tbl.loc[2,2], 'copia_percent':ltr_tbl.loc[2,4], 'gypsy_num':ltr_tbl.loc[3,1], 'gypsy_len':ltr_tbl.loc[3,2],
        'gypsy_percent':ltr_tbl.loc[3,4], 'erv_num':ltr_tbl.loc[4,1], 'erv_len':ltr_tbl.loc[4,2], 'erv_percent':ltr_tbl.loc[4,4]}],
        ignore_index = True)
    return breakpointlist
def read_tbl(tblname):
    ltr_tbl = pd.read_csv("{}".format(tblname), delim_whitespace=True, header=None)
    return ltr_tbl

breakpointlist = list_process(breakpointlist, a_chrlist)
print(breakpointlist)
breakpointlist = get_breakltr(breakpointlist, faname, species)
breakpointlist.to_csv('{}'.format(output_name), sep= '\t', header= True, index= False)
