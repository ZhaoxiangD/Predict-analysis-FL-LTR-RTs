import pandas as pd
import subprocess
import sys
hmmfile = sys.argv[1]
classes  = sys.argv[2]
if classes == 'gypsy':
    family = 'LTR/Gypsy '
elif classes == 'copia':
    family = 'LTR/Copia '
elif classes == 'belpao':
    family = 'LTR/Pao '
elif classes == 'erv':
    family = 'LTR/ERV1 '
else:
    print('Wrong family name!')
    sys.exit(1)
ltr_list = pd.read_csv("{0}".format(hmmfile), sep= "\t")
ltr_list_filt = ltr_list[ltr_list['family '] == family]
ltr_list_filt_F = ltr_list_filt[ltr_list_filt['unclass '] == 'False ']
ltr_list_filt_T = ltr_list_filt[ltr_list_filt['unclass '] == 'True ']
try:
    subprocess.check_call('rm {0}_all.cds'.format(classes), shell = True)
except:
    a ='a' 
for index, rows in ltr_list_filt_F.iterrows():
    try:
        subprocess.check_call('samtools faidx {3}_all.fa.fgenesh.cds {0}:{1}-{2} >> {3}_all.cds'.format(rows['proteinID '].strip(' '), rows['env_start ']*3-2, rows['env_stop ']*3, classes), shell = True)
        subprocess.check_call("sed -i 's/{0}:{1}-{2}/{0}:{3}-{4}/' {5}_all.cds".format(rows['proteinID '].strip(' '), rows['env_start ']*3-2, rows['env_stop ']*3, rows['env_start '], rows['env_stop '], classes), shell = True)
    except:
        continue
for index, rows in ltr_list_filt_T.iterrows():
    try:
        subprocess.check_call('samtools faidx unclassified_all.fa.fgenesh.cds {0}:{1}-{2} >> {3}_all.cds'.format(rows['proteinID '].strip(' '), rows['env_start ']*3-2, rows['env_stop ']*3, classes), shell = True)
        subprocess.check_call("sed -i 's/{0}:{1}-{2}/{0}:{3}-{4}/' {5}_all.cds".format(rows['proteinID '].strip(' '), rows['env_start ']*3-2, rows['env_stop ']*3, rows['env_start '], rows['env_stop '], classes), shell = True)
    except:
        continue
#for index, rows in sup_list.iterrows():
    #subprocess.check_call('samtools faidx gypsy_all.fa.fgenesh.cds {0}:{1}-{2} >> gypsy_all.cds'.format(rows['proteinID'], (rows['env_start']*3)-2, (rows['env_stop']*3)-2), shell = True)
    #subprocess.check_call("sed -i 's/{0}:{1}-{2}/{0}:{3}-{4}/' gypsy_all.cds".format(rows['proteinID'], (rows['env_start']*3)-2, (rows['env_stop']*3)-2, rows['env_start'], rows['env_stop']), shell = True)
