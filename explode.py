from numpy import NaN
import pandas as pd
from pandas.core import frame
import sys

file_name = sys.argv[1]
t = sys.argv[2]
if t == 'r':
    c='rear'
elif t == 'f':
    c='foward'
else:
    c='caried'

gene_list = pd.read_csv("{0}".format(file_name), sep= "\t", dtype=object)
gene_list = gene_list.drop(labels='Unnamed: 34', axis=1)
gene_list = gene_list.dropna()
gene_list[c+'_gene_start'] = gene_list[c+'_gene_start'].str.split('/')
gene_list = gene_list.explode(c+'_gene_start')
gene_list[c+'_gene_end'] = gene_list[c+'_gene_end'].str.split('/')
gene_list = gene_list.explode(c+'_gene_end')
gene_list[c+'_gene_type'] = gene_list[c+'_gene_type'].str.split('/')
gene_list = gene_list.explode(c+'_gene_type')
gene_list[[c+'_gene_start',c+'_gene_end']] = gene_list[[c+'_gene_start',c+'_gene_end']].astype(str)
gene_list = gene_list.drop_duplicates([c+'_gene_start', c+'_gene_end',c+'_gene_type' ], keep = 'last')
gene_list[c+'_gene_attributes'] = gene_list[c+'_gene_attributes'].str.split('/')
gene_list = gene_list.explode(c+'_gene_attributes')
gene_list = gene_list.drop_duplicates([c+'_gene_start', c+'_gene_end', c+'_gene_type'], keep = 'last')
gene_list[c+'_gene_score'] = gene_list[c+'_gene_score'].str.split('/')
gene_list = gene_list.explode(c+'_gene_score')
gene_list = gene_list.drop_duplicates([c+'_gene_start', c+'_gene_end', c+'_gene_type'], keep = 'last')
gene_list[c+'_gene_strand'] = gene_list[c+'_gene_strand'].str.split('/')
gene_list = gene_list.explode(c+'_gene_strand')
gene_list = gene_list.drop_duplicates([c+'_gene_start', c+'_gene_end', c+'_gene_type'], keep = 'last')
gene_list[c+'_gene_phase'] = gene_list[c+'_gene_phase'].str.split('/')
gene_list = gene_list.explode(c+'_gene_phase')
gene_list = gene_list.dropna()
gene_list.reset_index(drop = True, inplace = True)
gene_list[[c+'_gene_type',c+'_gene_start',c+'_gene_end',c+'_gene_score',c+'_gene_strand',c+'_gene_phase',
c+'_gene_attributes']] = gene_list[[c+'_gene_type',c+'_gene_start',c+'_gene_end',c+'_gene_score',c+'_gene_strand',c+'_gene_phase',
c+'_gene_attributes']].astype(str)
gene_list = gene_list.drop_duplicates([c+'_gene_start', c+'_gene_end', c+'_gene_type'], keep = 'last')
gene_list.reset_index(drop = True, inplace = True)
print(gene_list)
gene_list.to_csv(file_name+'.exp', sep = '\t')