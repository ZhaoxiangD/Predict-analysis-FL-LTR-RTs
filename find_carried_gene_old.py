import pandas as pd
from pandas.core import frame
import argparse
import subprocess
import time

def get_args():
    parser = argparse.ArgumentParser(description="find LTR carried gene")
    parser.add_argument('-i', '--input_file', default='hmmresult.xls')
    parser.add_argument('-o', '--output_file', default="carried_gene.csv")
    parser.add_argument('-g', '--gff_file', default= "FISH_ST_Oncorhynchus_keta.sorted.gff")
    args = parser.parse_args()
    input_name = args.input_file
    output_name = args.output_file
    gff_name = args.gff_file
    return input_name, output_name, gff_name

def get_ltr_gene(input_name):
    ltr_list = pd.read_csv("{0}".format(input_name), sep= "\t")
    ltr_list['chr_code'] = ltr_list['{}'.format('ID ')].str.split('_',1).apply(lambda x:x[-1])
    ltr_list['ltr_start'] = ltr_list['location '].str.split('-', 1).apply(lambda x:x[0])
    ltr_list['ltr_end'] = ltr_list['location '].str.split('-', 1).apply(lambda x:x[-1])
    ltr_list[['caried_gene_type','caried_gene_start','caried_gene_end','caried_gene_score',
    'caried_gene_strand','caried_gene_phase','caried_gene_attributes']] = 'No gene exist in the chromosome'
    return ltr_list

def get_gff(gff_name):
    gff_list = pd.read_csv("{0}".format(gff_name), header= None, names=None , sep= "\t")
    return gff_list

def find_loci(ltr_list, gff_list):
    for i, rows in ltr_list.iterrows():
        gff_list_filt = gff_list[gff_list[0] == rows['chr_code'].strip(' ')]
        gff_list_filt.reset_index(drop = True, inplace = True)
        mid = round(len(gff_list_filt)/2)
        if mid == 0:
            continue
        if int(rows['ltr_start']) < int(gff_list_filt.loc[mid,3]) and int(rows['ltr_end']) > int(gff_list_filt.loc[mid,4]):
            ltr_list.loc[i,'caried_gene_type'] = gff_list_filt.loc[mid,2]
            ltr_list.loc[i,'caried_gene_start'] = gff_list_filt.loc[mid,3]
            ltr_list.loc[i,'caried_gene_end'] = gff_list_filt.loc[mid,4]
            ltr_list.loc[i,'caried_gene_score'] = gff_list_filt.loc[mid,5]
            ltr_list.loc[i,'caried_gene_strand'] = gff_list_filt.loc[mid,6]
            ltr_list.loc[i,'caried_gene_phase'] = gff_list_filt.loc[mid,7]
            ltr_list.loc[i,'caried_gene_attributes'] = gff_list_filt.loc[mid,8]
        elif int(rows['ltr_end']) < int(gff_list_filt.loc[mid,3]):
            g_type, g_start, g_end, g_score, g_strand, g_phase, g_attri = [],[],[],[],[],[],[]
            gff_list_filt_half = gff_list_filt.loc[:mid, :]
            gff_list_filt_half.reset_index(drop = True, inplace = True)
            for gi, grows in gff_list_filt_half.iterrows():
                if int(rows['ltr_start']) < int(grows[3]) and int(rows['ltr_end']) > grows[4]:
                    g_type.append(grows[2])
                    g_start.append(str(grows[3]))
                    g_end.append(str(grows[4]))
                    g_score.append(str(grows[5]))
                    g_strand.append(str(grows[6]))
                    g_phase.append(str(grows[7]))
                    g_attri.append(grows[8])
            ltr_list.loc[i,'caried_gene_type'] = "/".join(g_type)
            ltr_list.loc[i,'caried_gene_start'] = "/".join(g_start)
            ltr_list.loc[i,'caried_gene_end'] = "/".join(g_end)
            ltr_list.loc[i,'caried_gene_score'] = "/".join(g_score)
            ltr_list.loc[i,'caried_gene_strand'] = "/".join(g_strand)
            ltr_list.loc[i,'caried_gene_phase'] = "/".join(g_phase)
            ltr_list.loc[i,'caried_gene_attributes'] = "/".join(g_attri)
        elif int(rows['ltr_start']) > int(gff_list_filt.loc[mid,4]):
            g_type, g_start, g_end, g_score, g_strand, g_phase, g_attri = [],[],[],[],[],[],[]
            gff_list_filt_half = gff_list_filt.loc[mid:, :]
            gff_list_filt_half.reset_index(drop = True, inplace = True)
            for gi, grows in gff_list_filt_half.iterrows():
                if int(rows['ltr_start']) < int(grows[3]) and int(rows['ltr_end']) > grows[4]:
                    g_type.append(grows[2])
                    g_start.append(str(grows[3]))
                    g_end.append(str(grows[4]))
                    g_score.append(str(grows[5]))
                    g_strand.append(str(grows[6]))
                    g_phase.append(str(grows[7]))
                    g_attri.append(grows[8])
            ltr_list.loc[i,'caried_gene_type'] = "/".join(g_type)
            ltr_list.loc[i,'caried_gene_start'] = "/".join(g_start)
            ltr_list.loc[i,'caried_gene_end'] = "/".join(g_end)
            ltr_list.loc[i,'caried_gene_score'] = "/".join(g_score)
            ltr_list.loc[i,'caried_gene_strand'] = "/".join(g_strand)
            ltr_list.loc[i,'caried_gene_phase'] = "/".join(g_phase)
            ltr_list.loc[i,'caried_gene_attributes'] = "/".join(g_attri)
    return ltr_list

def main():
    input_name, output_name, gff_name = get_args()
    ltr_list = get_ltr_gene(input_name)
    gff_list = get_gff(gff_name)
    time_start = time.time()
    ltr_list = find_loci(ltr_list,gff_list)
    time_stop = time.time()
    print(ltr_list)
    print(time_stop - time_start)
    ltr_list.to_csv('{}'.format(output_name), sep= '\t', header= True, index= False)

main()