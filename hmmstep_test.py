import time
import readxls #custom scirp to read previous step's xls file
import readdom # custom scrip to read hmmscan's dom file
import subprocess
xlsdict = readxls.xlsdict
'''
following files must be deleted if this script need to be runned again

files created : gypsy/copia/belpao/erv/unclassified_all.fa
*F/fgenesh*
*gypsy/copia/erv/belpao.aa

'''
#domdict = readdom.domdict
def getfgenesh(xlsdict):
    '''
    categorlize LTR's fasta file in to summary fasta file differentiated by superfamily name
    fgenesh summary fasta file
    
    summary file = <family>_all.fa
    aa file = <family>_all.fa.fgenesh.aa
    '''
    copia = 'False_1'
    gypsy = 'False_1'
    belpao = 'False_1'
    erv = 'False_1'
    unclass = 'False_1'
    for num in range(len(xlsdict)): # categorlize ltr's fasta file
        if xlsdict[num]['family'] == 'LTR/Gypsy':
            catgypsy = subprocess.Popen('cat {0}/{1}_{0}.fa >> gypsy_all.fa'.format(xlsdict[num]['chr'], xlsdict[num]['ID']), shell = True, stdout = open('/dev/null','w')) 
            catgypsy.wait()
            gypsy = 'True_gypsy' #True + family name
        elif xlsdict[num]['family'] == 'LTR/Copia':
            catcopia = subprocess.Popen('cat {0}/{1}_{0}.fa >> copia_all.fa'.format(xlsdict[num]['chr'], xlsdict[num]['ID']), shell = True, stdout = open('/dev/null','w'))
            catcopia.wait()
            copia = 'True_copia'
        elif xlsdict[num]['family'] == 'LTR/Pao':
            catbelpao = subprocess.Popen('cat {0}/{1}_{0}.fa >> belpao_all.fa'.format(xlsdict[num]['chr'], xlsdict[num]['ID']), shell = True, stdout = open('/dev/null','w'))
            catbelpao.wait()
            belpao = 'True_belpao'
        elif xlsdict[num]['family'] == 'LTR/ERV1':
            caterv = subprocess.Popen('cat {0}/{1}_{0}.fa >> erv_all.fa'.format(xlsdict[num]['chr'], xlsdict[num]['ID']), shell = True, stdout = open('/dev/null','w'))
            caterv.wait()
            erv = 'True_erv'
        else:
            catunclass = subprocess.Popen('cat {0}/{1}_{0}.fa >> unclassified_all.fa'.format(xlsdict[num]['chr'], xlsdict[num]['ID']), shell = True, stdout = open('/dev/null','w'))
            catunclass.wait()
            unclass = 'True_unclassified'
    familylist = [gypsy, copia, belpao, erv, unclass]
    for familyname in range(len(familylist)):
        if familylist[familyname].split('_')[0] == 'True':
            name = familylist[familyname].split('_')[1]
            fgenesh_step1 = subprocess.check_call('perl /home/yinglu/pl/run_FgeneSH_fish.pl {0}_all.fa {0}_all.fa.fgenesh seq & '.format(name), shell = True, stdout = open('/dev/null','w'))
            time.sleep(5)
    time.sleep(5)
    '''
    DO NOT CHANGE following function's sentences
    '''
    if gypsy == 'True_gypsy':
        fgenesh_step2 = subprocess.check_call('perl /home/yinglu/pl/linux_fgenesh_outseq_pick.pl gypsy_all.fa.fgenesh gypsy_all.fa.fgenesh.cds gypsy_all.fa.fgenesh.aa ', shell = True, stdout = open('/dev/null','w'))
        #fgenesh_step2.wait()
        time.sleep(5)
    if copia == 'True_copia':
        fgenesh_step2_copia = subprocess.check_call('perl /home/yinglu/pl/linux_fgenesh_outseq_pick.pl copia_all.fa.fgenesh copia_all.fa.fgenesh.cds copia_all.fa.fgenesh.aa', shell = True, stdout = open('/dev/null','w'))
        #fgenesh_step2_copia.wait()
        time.sleep(5)
    if belpao == 'True_belpao':
        fgenesh_step2_belpao = subprocess.check_call('perl /home/yinglu/pl/linux_fgenesh_outseq_pick.pl belpao_all.fa.fgenesh belpao_all.fa.fgenesh.cds belpao_all.fa.fgenesh.aa', shell =True, stdout = open('/dev/null','w'))
        #fgenesh_step2_belpao.wait()
        time.sleep(5)
    if erv == 'True_erv':
        fgenesh_step2_erv = subprocess.check_call('perl /home/yinglu/pl/linux_fgenesh_outseq_pick.pl erv_all.fa.fgenesh erv_all.fa.fgenesh.cds erv_all.fa.fgenesh.aa', shell =True, stdout = open('/dev/null','w'))
        #fgenesh_step2_erv.wait()
        time.sleep(5)
    if unclass == 'True_unclassified':
        fgenesh_step2_un = subprocess.check_call('perl /home/yinglu/pl/linux_fgenesh_outseq_pick.pl unclassified_all.fa.fgenesh unclassified_all.fa.fgenesh.cds unclassified_all.fa.fgenesh.aa', shell = True, stdout = open('/dev/null','w'))
        #fgenesh_step2_un.wait()
        time.sleep(10)
    return familylist

def gethmm(familylist):
    '''
    hmmscan summary fasta file
    domfile = <family>_all_RT.dom
    unclassfied file remain to be done
    '''
    for name in range(len(familylist)):
        if familylist[name].split('_')[0] == 'True' and familylist[name].split('_')[1] != 'unclassified':
            selectedname = familylist[name].split('_')[1]
            hmmstep = subprocess.check_call('hmmscan -E 1E-20 --incE 1E-20 --incdomE 1E-20 --domE 1E-20 --noali -o {0}_all_RT.txt --tblout {0}_all_RT.tbl --domtblout {0}_all_RT.dom --pfamtblout {0}_all_RT.pfa ../lib/RT_{0}.hmm {0}_all.fa.fgenesh.aa'.format(selectedname), shell = True)
        else:
            for i in range(len(familylist)-1):
                scan_name = familylist[i].split('_')[1]
                hmmstep = subprocess.check_call('hmmscan -E 1E-20 --incE 1E-20 --incdomE 1E-20 --domE 1E-20 --noali -o {0}_unclass_all_RT.txt --tblout {0}_unclass_all_RT.tbl --domtblout {0}_unclass_all_RT.dom --pfamtblout {0}_unclass_all_RT.pfa ../lib/RT_{0}.hmm unclassified_all.fa.fgenesh.aa'.format(scan_name), shell = True)
    return

def getdom(truefamily):
    domdict = []
    for i in range(len(truefamily)):
        domdict.append(readdom.readdom(truefamily[i]))
    return domdict

def getunclassdom():
    unclass_gypsy_domdict = []
    unclass_copia_domdict = []
    unclass_belpao_domdict = []
    unclass_erv_domdict = []
    unclass_gypsy_domdict.append(readdom.readdom('gypsy_unclass'))
    unclass_copia_domdict.append(readdom.readdom('copia_unclass'))
    unclass_belpao_domdict.append(readdom.readdom('belpao_unclass'))
    unclass_erv_domdict.append(readdom.readdom('erv_unclass'))
    unclass_gypsy_domdict = unclass_gypsy_domdict[0]
    unclass_copia_domdict = unclass_copia_domdict[0]
    unclass_belpao_domdict = unclass_belpao_domdict[0]
    unclass_erv_domdict = unclass_erv_domdict[0]
    return unclass_gypsy_domdict, unclass_copia_domdict, unclass_belpao_domdict, unclass_erv_domdict

def unclassdom_process(unclass_gypsy_domdict, unclass_copia_domdict, unclass_belpao_domdict, unclass_erv_domdict, selected_domdict):
    selected_gypsy_domdict = []
    selected_copia_domdict = []
    selected_belpao_domdict = []
    selected_erv_domdict = []
    ID_log_list = []
    tmpdict = {'domID': '0', 'sequence_Evalue': 10000}
    unclass_gypsy_exist = False if len(unclass_gypsy_domdict) == 0 else True
    unclass_copia_exist = False if len(unclass_copia_domdict) == 0 else True
    unclass_belpao_exist = False if len(unclass_belpao_domdict) == 0 else True
    unclass_erv_exist = False if len(unclass_erv_domdict) == 0 else True
    if unclass_gypsy_exist:
        selected_gypsy_domdict.append(unclass_gypsy_domdict[0])
        for gi in range(1,len(unclass_gypsy_domdict)):
            if unclass_gypsy_domdict[gi]['domID'] != unclass_gypsy_domdict[gi-1]['domID']:
                selected_gypsy_domdict.append(unclass_gypsy_domdict[gi])
    if unclass_copia_exist:
        selected_copia_domdict.append(unclass_copia_domdict[0])
        for ci in range(1,len(unclass_copia_domdict)):
            if unclass_copia_domdict[ci]['domID'] != unclass_copia_domdict[ci-1]['domID']:
                selected_copia_domdict.append(unclass_copia_domdict[gi])
    else:
        selected_copia_domdict.append(tmpdict)
    if unclass_belpao_exist:
        selected_belpao_domdict.append(unclass_belpao_domdict[0])
        for bi in range(1,len(unclass_belpao_domdict)):
            if unclass_belpao_domdict[bi]['domID'] != unclass_belpao_domdict[bi-1]['domID']:
                selected_belpao_domdict.append(unclass_belpao_domdict[bi])
    else:
        selected_belpao_domdict.append(tmpdict) 
    if unclass_erv_exist:
        selected_erv_domdict.append(unclass_erv_domdict[0])
        for ei in range(1,len(unclass_erv_domdict)):
            if unclass_erv_domdict[ei]['domID'] != unclass_erv_domdict[ei-1]['domID']:
                selected_erv_domdict.append(unclass_erv_domdict[ei])
    else:
        selected_erv_domdict.append(tmpdict)
    #selected_gypsy_domdict = selected_gypsy_domdict[0]
    #print(selected_gypsy_domdict)
    #selected_copia_domdict = selected_copia_domdict[0]
    #selected_belpao_domdict = selected_belpao_domdict[0]
    #selected_erv_domdict = selected_erv_domdict[0]
    #print(selected_erv_domdict)
    for g in range(len(selected_gypsy_domdict)):
        #g = int(g)
        GID = selected_gypsy_domdict[g]['domID']
        #print(GID)
        #print(type(selected_gypsy_domdict[g]))
        minvalue = float(selected_gypsy_domdict[g]['sequence_Evalue'])
        mindict = selected_gypsy_domdict[g]
        mindict['family'] = 'LTR/Gypsy'
        ID_log_list.append(GID)
        for ci in range(len(selected_copia_domdict)):
            if GID == selected_copia_domdict[ci]['domID']:
                if minvalue > float(selected_copia_domdict[ci]['sequence_Evalue']):
                    minvalue = float(selected_copia_domdict[ci]['sequence_Evalue'])
                    #minfamily = 'copia'
                    mindict = selected_copia_domdict[ci]
                    mindict['family'] = 'LTR/Copia'
        for bi in range(len(selected_belpao_domdict)):
            if GID == selected_belpao_domdict[bi]['domID']:
                if minvalue > float(selected_belpao_domdict[bi]['sequence_Evalue']):
                    minvalue = float(selected_belpao_domdict[bi]['sequence_Evalue'])
                    #minfamily = 'belpao'
                    mindict = selected_belpao_domdict[bi]
                    mindict['family'] = 'LTR/Pao'
        for ei in range(len(selected_copia_domdict)):
            if GID == selected_copia_domdict[ei]['domID']:
                if minvalue > float(selected_erv_domdict[ei]['sequence_Evalue']):
                    minvalue = float(selected_erv_domdict[ei]['sequence_Evalue'])
                    #minfamily = 'erv'
                    mindict = selected_erv_domdict[ei]
                    mindict['family'] = 'LTR/ERV1'
        mindict['unclass'] = 'True'
        #print(mindict)
        selected_domdict.append(mindict)
    #print(selected_domdict)
    for ci in range(len(selected_copia_domdict)):
        if selected_copia_domdict[ci]['domID'] in ID_log_list == False and selected_copia_domdict[ci]['domID'] != '0':
            ID = selected_copia_domdict[ci]['domID']
            minvalue = float(selected_copia_domdict[ci]['sequence_Evalue'])
            mindict = selected_copia_domdict[ci]
            mindict['family'] = 'LTR/Copia'
            ID_log_list.append(ID)
            for gi in range(len(selected_gypsy_domdict)):
                if ID == selected_gypsy_domdict[gi]['domID']:
                    if minvalue > float(selected_gypsy_domdict[gi]['sequence_Evalue']):
                        minvalue = float(selected_gypsy_domdict[gi]['sequence_Evalue'])
                        #minfamily = 'copia'
                        mindict = selected_gypsy_domdict[gi]
                        mindict['family'] = 'LTR/Gypsy'
            for bi in range(len(selected_belpao_domdict)):
                if ID == selected_belpao_domdict[bi]['domID']:
                    if minvalue > float(selected_belpao_domdict[bi]['sequence_Evalue']):
                        minvalue = float(selected_belpao_domdict[bi]['sequence_Evalue'])
                        #minfamily = 'belpao'
                        mindict = selected_belpao_domdict[bi]
                        mindict['family'] = 'LTR/Pao'
            for ei in range(len(selected_copia_domdict)):
                if ID == selected_copia_domdict[ei]['domID']:
                    if minvalue > float(selected_erv_domdict[ei]['sequence_Evalue']):
                        minvalue = float(selected_gypsy_domdict[gi]['sequence_Evalue'])
                        #minfamily = 'copia'
                        mindict = selected_gypsy_domdict[gi]
                        mindict['family'] = 'LTR/Gypsy'
            for bi in range(len(selected_belpao_domdict)):
                if ID == selected_belpao_domdict[bi]['domID']:
                    if minvalue > float(selected_belpao_domdict[bi]['sequence_Evalue']):
                        minvalue = float(selected_belpao_domdict[bi]['sequence_Evalue'])
                        #minfamily = 'belpao'
                        mindict = selected_belpao_domdict[bi]
                        mindict['family'] = 'LTR/Pao'
            for ei in range(len(selected_copia_domdict)):
                if ID == selected_copia_domdict[ei]['domID']:
                    if minvalue > float(selected_erv_domdict[ei]['sequence_Evalue']):
                        minvalue = float(selected_erv_domdict[ei]['sequence_Evalue'])
                        #minfamily = 'erv'
                        mindict = selected_erv_domdict[ei]
                        mindict['family'] = 'LTR/ERV1'
            mindict['unclass'] = 'True'
            selected_domdict.append(mindict)

    for bi in range(len(selected_belpao_domdict)):
        if selected_belpao_domdict[bi]['domID'] in ID_log_list == False and selected_belpao_domdict[bi]['domID'] != '0':
            ID = selected_belpao_domdict[bi]['domID']
            minvalue = float(selected_belpao_domdict[bi]['sequence_Evalue'])
            mindict = selected_belpao_domdict[bi]
            mindict['family'] = 'LTR/Pao'
            ID_log_list.append(ID)
            for gi in range(len(selected_gypsy_domdict)):
                if ID == selected_gypsy_domdict[gi]['domID']:
                    if minvalue > float(selected_gypsy_domdict[gi]['sequence_Evalue']):
                        minvalue = float(selected_gypsy_domdict[gi]['sequence_Evalue'])
                        #minfamily = 'copia'
                        mindict = selected_gypsy_domdict[gi]
                        mindict['family'] = 'LTR/Gypsy'
            for ci in range(len(selected_copia_domdict)):
                if ID == selected_copia_domdict[ci]['domID']:
                    if minvalue > float(selected_copia_domdict[ci]['sequence_Evalue']):
                        minvalue = float(selected_copia_domdict[ci]['sequence_Evalue'])
                        #minfamily = 'belpao'
                        mindict = selected_copia_domdict[ci]
                        mindict['family'] = 'LTR/Copia'
            for ei in range(len(selected_erv_domdict)):
                if ID == selected_erv_domdict[ei]['domID']:
                    if minvalue > float(selected_erv_domdict[ei]['sequence_Evalue']):
                        minvalue = float(selected_erv_domdict[ei]['sequence_Evalue'])
                        #minfamily = 'erv'
                        mindict = selected_erv_domdict[ei]
                        mindict['family'] = 'LTR/ERV1'
            mindict['unclass'] = 'True'
            selected_domdict.append(mindict)
    for ei in range(len(selected_erv_domdict)):
        #print(len(selected_erv_domdict))
        #print(ei)
        if selected_erv_domdict[ei]['domID'] in ID_log_list == False and selected_erv_domdict[ei]['domID'] != '0':
            ID = selected_erv_domdict[ei]['domID']
            minvalue = float(selected_erv_domdict[ei]['sequence_Evalue'])
            mindict = selected_erv_domdict[ei]
            mindict['family'] = 'LTR/ERV1'
            ID_log_list.append(ID)
            for gi in range(len(selected_gypsy_domdict)):
                if ID == selected_gypsy_domdict[gi]['domID']:
                    if minvalue > float(selected_gypsy_domdict[gi]['sequence_Evalue']):
                        minvalue = float(selected_gypsy_domdict[gi]['sequence_Evalue'])
                        #minfamily = 'copia'
                        mindict = selected_gypsy_domdict[gi]
                        mindict['family'] = 'LTR/Gypsy'
            for ci in range(len(selected_copia_domdict)):
                if ID == selected_copia_domdict[ci]['domID']:
                    if minvalue > float(selected_copia_domdict[ci]['sequence_Evalue']):
                        minvalue = float(selected_copia_domdict[ci]['sequence_Evalue'])
                        #minfamily = 'belpao'
                        mindict = selected_copia_domdict[ci]
                        mindict['family'] = 'LTR/Copia'
            for bi in range(len(selected_belpao_domdict)):
                if ID == selected_belpao_domdict[bi]['domID']:
                    if minvalue > float(selected_belpao_domdict[bi]['sequence_Evalue']):
                        minvalue = float(selected_belpao_domdict[bi]['sequence_Evalue'])
                        #minfamily = 'erv'
                        mindict = selected_belpao_domdict[bi]
                        mindict['family'] = 'LTR/Pao'
            mindict['unclass'] = 'True'
            selected_domdict.append(mindict)
    #print(selected_domdict)
    return selected_domdict

def domdict_process(domdict):
    selected_domdict = []
    for family in range(len(domdict)):
        selected_domdict.append(domdict[family][0])
        for i in range(1,len(domdict[family])):
            if domdict[family][i]['domID'] != domdict[family][i-1]['domID']:
                domdict[family][i]['unclass'] = 'False'
                selected_domdict.append(domdict[family][i])
    return selected_domdict

def ali_dict(selected_domdict,xlsdict):
    for xls in range(len(xlsdict)):
        for dom in range(len(selected_domdict)):
            if xlsdict[xls]['domID'] == selected_domdict[dom]['domID']:
                for key in selected_domdict[dom].keys():
                    xlsdict[xls][key] = selected_domdict[dom][key] 
            else:
                continue
            #elif xlsdict[xls]['domID'] == domdict[dom]['domID'] and domdict[dom]['unclass'] == 'True':
                #for key in domdict[dom].keys():
                    #xlsdict[xls][key] = domdict[dom][key] 
    #print(xlsdict)            
    return xlsdict

def getaafile(alidict):
    for i in range(len(alidict)):
        if alidict[i].get('unclass') == 'False':
            if alidict[i]['family'] == 'LTR/Gypsy' and 'subfamily' in alidict[i]:
                samtools = subprocess.check_call('samtools faidx gypsy_all.fa.fgenesh.aa {0}:{1}-{2} >> {0}_{3}_gypsy.aa'.format(alidict[i]['proteinID'], alidict[i]['env_start'], alidict[i]['env_stop'], alidict[i]['chr']), shell = True)
            elif alidict[i]['family'] == 'LTR/Copia' and 'subfamily' in alidict[i]:
                samtools1 = subprocess.check_call('samtools faidx copia_all.fa.fgenesh.aa {0}:{1}-{2} >> {0}_{3}_copia.aa'.format(alidict[i]['proteinID'], alidict[i]['env_start'], alidict[i]['env_stop'], alidict[i]['chr']), shell = True)
            elif alidict[i]['family'] == 'LTR/Pao' and 'subfamily' in alidict[i]:
                samtools2 = subprocess.check_call('samtools faidx belpao_all.fa.fgenesh.aa {0}:{1}-{2} >> {0}_{3}_belpao.aa'.format(alidict[i]['proteinID'], alidict[i]['env_start'], alidict[i]['env_stop'], alidict[i]['chr']), shell = True)
            elif alidict[i]['family'] == 'LTR/ERV1' and 'subfamily' in alidict[i]:
                samtools3 = subprocess.check_call('samtools faidx erv_all.fa.fgenesh.aa {0}:{1}-{2} >> {0}_{3}_erv.aa'.format(alidict[i]['proteinID'], alidict[i]['env_start'], alidict[i]['env_stop'], alidict[i]['chr']), shell = True)
        elif alidict[i].get('unclass') == 'True':  
            if alidict[i]['family'] == 'LTR/Gypsy' and 'subfamily' in alidict[i]:
                samtools = subprocess.check_call('samtools faidx unclassified_all.fa.fgenesh.aa {0}:{1}-{2} >> {0}_{3}_gypsy.aa'.format(alidict[i]['proteinID'], alidict[i]['env_start'], alidict[i]['env_stop'], alidict[i]['chr']), shell = True)
            elif alidict[i]['family'] == 'LTR/Copia' and 'subfamily' in alidict[i]:
                samtools1 = subprocess.check_call('samtools faidx unclassified_all.fa.fgenesh.aa {0}:{1}-{2} >> {0}_{3}_copia.aa'.format(alidict[i]['proteinID'], alidict[i]['env_start'], alidict[i]['env_stop'], alidict[i]['chr']), shell = True)
            elif alidict[i]['family'] == 'LTR/Pao' and 'subfamily' in alidict[i]:
                samtools2 = subprocess.check_call('samtools faidx unclassified_all.fa.fgenesh.aa {0}:{1}-{2} >> {0}_{3}_belpao.aa'.format(alidict[i]['proteinID'], alidict[i]['env_start'], alidict[i]['env_stop'], alidict[i]['chr']), shell = True)
            elif alidict[i]['family'] == 'LTR/ERV1' and 'subfamily' in alidict[i]:
                samtools3 = subprocess.check_call('samtools faidx unclassified_all.fa.fgenesh.aa {0}:{1}-{2} >> {0}_{3}_erv.aa'.format(alidict[i]['proteinID'], alidict[i]['env_start'], alidict[i]['env_stop'], alidict[i]['chr']), shell = True)
 
        else:
            continue
    return


familylist = getfgenesh(xlsdict)
#familylist = ['Ture_gypsy','Ture_copia','Ture_belpao','Ture_erv']
truefamily=[]
for name in range(len(familylist)):
    if familylist[name].split('_')[0] == 'True' and familylist[name].split('_')[1] != 'unclassified':
        truefamily.append(familylist[name].split('_')[1])
hmmact = gethmm(familylist)
#truefamily=['gypsy','belpao','copia','erv']
domdict = getdom(truefamily)
unclass_gypsy_domdict, unclass_copia_domdict, unclass_belpao_domdict, unclass_erv_domdict = getunclassdom()
selected_domdict = domdict_process(domdict)
selected_domdict = unclassdom_process(unclass_gypsy_domdict, unclass_copia_domdict, unclass_belpao_domdict, unclass_erv_domdict, selected_domdict)
alidict = ali_dict(selected_domdict, xlsdict)
aastep = getaafile(alidict)
#print(alidict[0]['family'] == 'LTR/Gypsy')
#print('subfamily' in alidict[0])
#print(alidict[0]['family'] == 'LTR/Gypsy' and 'subfamily' in alidict[0])
for i in range(len(alidict)):
    if 'unclass' in alidict[i]:
        for key in alidict[i].keys():
            print(key, '\t', end='\r')
        print('\n')
        break

for i in range(len(alidict)):
    print(alidict[i])
