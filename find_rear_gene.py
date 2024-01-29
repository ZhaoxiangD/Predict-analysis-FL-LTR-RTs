import pandas as pd
from pandas.core import frame
import argparse
import subprocess
import time

def get_args():
    parser = argparse.ArgumentParser(description="find LTR rear gene")
    parser.add_argument('-i', '--input_file', default='hmmresult.xls')
    parser.add_argument('-o', '--output_file', default="rear_gene.csv")
    parser.add_argument('-g', '--gff_file', default= "GCF_006149115.1_Oner_1.0_genomic.sorted.gff")
    parser.add_argument('-s', '--synaty_file', default= "test.1.txt")
    parser.add_argument('-b', '--bp_num', default= "10000", type=int)
    args = parser.parse_args()
    input_name = args.input_file
    output_name = args.output_file
    gff_name = args.gff_file
    syn_name = args.synaty_file
    bp = args.bp_num
    return input_name, output_name, gff_name, syn_name, bp

def get_ltr_gene(input_name):
    ltr_list = pd.read_csv("{0}".format(input_name), sep= "\t")
    ltr_list['chr_code'] = ltr_list['{}'.format('ID ')].str.split('_',1).apply(lambda x:x[-1])
    ltr_list['ltr_start'] = ltr_list['location '].str.split('-', 1).apply(lambda x:x[0])
    ltr_list['ltr_end'] = ltr_list['location '].str.split('-', 1).apply(lambda x:x[-1])
    ltr_list[['rear_gene_type','rear_gene_start','rear_gene_end','rear_gene_score',
    'rear_gene_strand','rear_gene_phase','rear_gene_attributes']] = 'No gene exist in the chromosome'
    return ltr_list

def get_gff(gff_name):
    gff_list = pd.read_csv("{0}".format(gff_name), header= None, names=None , sep= "\t")
    return gff_list

def get_synaty(syn_name):
    syn_list = pd.read_csv("{}".format(syn_name), header=None, names=None, sep="\t")
    return syn_list

def find_loci(ltr_list, gff_list, syn_list, bp):
    for i, rows in ltr_list.iterrows():
        gff_list_filt = gff_list[gff_list[0] == rows['chr_code'].strip(' ')]
        gff_list_filt.reset_index(drop = True, inplace = True)
        ltr_foward_start = int(rows['ltr_start']) - bp
        if ltr_foward_start < 0:
            ltr_foward_start = 0
        ltr_rear_end = int(rows['ltr_end']) + bp
        mid = round(len(gff_list_filt)/2)
        if mid == 0:
            continue
        if int(rows['ltr_end']) < int(gff_list_filt.loc[mid,3]) and ltr_rear_end > int(gff_list_filt.loc[mid,4]):
            syn_result = syn_filt(rows['chr_code'], gff_list_filt.loc[mid,3], gff_list_filt.loc[mid,4], syn_list)
            if syn_result:
                ltr_list.loc[i,'rear_gene_type'] = gff_list_filt.loc[mid,2]
                ltr_list.loc[i,'rear_gene_start'] = gff_list_filt.loc[mid,3]
                ltr_list.loc[i,'rear_gene_end'] = gff_list_filt.loc[mid,4]
                ltr_list.loc[i,'rear_gene_score'] = gff_list_filt.loc[mid,5]
                ltr_list.loc[i,'rear_gene_strand'] = gff_list_filt.loc[mid,6]
                ltr_list.loc[i,'rear_gene_phase'] = gff_list_filt.loc[mid,7]
                ltr_list.loc[i,'rear_gene_attributes'] = gff_list_filt.loc[mid,8]
        elif ltr_rear_end < int(gff_list_filt.loc[mid,3]):
            g_type, g_start, g_end, g_score, g_strand, g_phase, g_attri = ['NA'],['NA'],['NA'],['NA'],['NA'],['NA'],['NA']
            gff_list_filt_half = gff_list_filt.loc[:mid, :]
            gff_list_filt_half.reset_index(drop = True, inplace = True)
            for gi, grows in gff_list_filt_half.iterrows():
                if int(rows['ltr_end']) < int(grows[3]) and ltr_rear_end > grows[4]:
                    syn_result = syn_filt(rows['chr_code'], grows[3], grows[4], syn_list)
                    if syn_result:
                        if g_type[0] == 'NA':
                            g_type, g_start, g_end, g_score, g_strand, g_phase, g_attri = [],[],[],[],[],[],[]
                        g_type.append(grows[2])
                        g_start.append(str(grows[3]))
                        g_end.append(str(grows[4]))
                        g_score.append(str(grows[5]))
                        g_strand.append(str(grows[6]))
                        g_phase.append(str(grows[7]))
                        g_attri.append(grows[8])
            ltr_list.loc[i,'rear_gene_type'] = "/".join(g_type)
            ltr_list.loc[i,'rear_gene_start'] = "/".join(g_start)
            ltr_list.loc[i,'rear_gene_end'] = "/".join(g_end)
            ltr_list.loc[i,'rear_gene_score'] = "/".join(g_score)
            ltr_list.loc[i,'rear_gene_strand'] = "/".join(g_strand)
            ltr_list.loc[i,'rear_gene_phase'] = "/".join(g_phase)
            ltr_list.loc[i,'rear_gene_attributes'] = "/".join(g_attri)
        elif int(rows['ltr_end']) > int(gff_list_filt.loc[mid,4]):
            g_type, g_start, g_end, g_score, g_strand, g_phase, g_attri = ['NA'],['NA'],['NA'],['NA'],['NA'],['NA'],['NA']
            gff_list_filt_half = gff_list_filt.loc[mid:, :]
            gff_list_filt_half.reset_index(drop = True, inplace = True)
            for gi, grows in gff_list_filt_half.iterrows():
                if int(rows['ltr_end']) < int(grows[3]) and ltr_rear_end > grows[4]:
                    syn_result = syn_filt(rows['chr_code'], grows[3], grows[4], syn_list)
                    if syn_result:
                        if g_type[0] == 'NA':
                            g_type, g_start, g_end, g_score, g_strand, g_phase, g_attri = [],[],[],[],[],[],[]
                        g_type.append(grows[2])
                        g_start.append(str(grows[3]))
                        g_end.append(str(grows[4]))
                        g_score.append(str(grows[5]))
                        g_strand.append(str(grows[6]))
                        g_phase.append(str(grows[7]))
                        g_attri.append(grows[8])
            ltr_list.loc[i,'rear_gene_type'] = "/".join(g_type)
            ltr_list.loc[i,'rear_gene_start'] = "/".join(g_start)
            ltr_list.loc[i,'rear_gene_end'] = "/".join(g_end)
            ltr_list.loc[i,'rear_gene_score'] = "/".join(g_score)
            ltr_list.loc[i,'rear_gene_strand'] = "/".join(g_strand)
            ltr_list.loc[i,'rear_gene_phase'] = "/".join(g_phase)
            ltr_list.loc[i,'rear_gene_attributes'] = "/".join(g_attri)
    return ltr_list

def syn_filt(ltr_chr, gene_start, gene_end, syn_list):
    syn_list_filt = syn_list[syn_list[0] == ltr_chr.strip(' ')]
    syn_list_filt.reset_index(drop = True, inplace = True)
    synaty_result = True
    for i, rows in syn_list_filt.iterrows():
        if gene_start == rows[1] and gene_end == rows[2]:
            synaty_result = False
    return synaty_result


def main():
    input_name, output_name, gff_name, syn_name, bp= get_args()
    ltr_list = get_ltr_gene(input_name)
    gff_list = get_gff(gff_name)
    syn_list = get_synaty(syn_name)
    #print(syn_list)
    time_start = time.time()
    ltr_list = find_loci(ltr_list,gff_list, syn_list, bp)
    time_stop = time.time()
    #print(ltr_list)
    print(time_stop - time_start)
    ltr_list.to_csv('{}'.format(output_name), sep= '\t', header= True, index= False)

main()