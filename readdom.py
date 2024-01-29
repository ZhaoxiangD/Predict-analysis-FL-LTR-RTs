def readdom(name):
    with open ('{0}_all_RT.dom'.format(name),'r') as f:
        domfile=f.readlines()
    domfile=domfile[3:-10]
    domtmp = {}
    domdict = []
    for i in range(len(domfile)):
        domfile[i]=domfile[i].split(' ')
        for num in range(len(domfile[i])-1,0,-1):
            if domfile[i][num] == '':
                domfile[i].pop(num)        
        domtmp['protein_domain'] = domfile[i][0].split('_')[0]
        domtmp['subfamily'] = domfile[i][0].split('_')[1]
        domtmp['proteinID'] = domfile[i][3]
        if domfile[i][3].split('_')[0] != 'NC':
            domtmp['domID'] = domfile[i][3].split('_')[0]
        else:
            domtmp['domID'] = domfile[i][3].split('_')[0] + ':' + domfile[i][3].split('_')[1].split(':')[1]
        domtmp['qlen'] = domfile[i][5]
        domtmp['sequence_Evalue'] = domfile[i][6]
        domtmp['sequence_score'] = domfile[i][7]
        domtmp['sequence_bias'] = domfile[i][8]
        domtmp['domain_num'] = domfile[i][9]
        domtmp['domain_c-Evalue'] = domfile[i][11]
        domtmp['domain_i-Evalue'] = domfile[i][12]
        domtmp['domain_score'] = domfile[i][13]
        domtmp['domain_bias'] = domfile[i][14]
        domtmp['env_start'] = domfile[i][19]
        domtmp['env_stop'] = domfile[i][20]
        domtmp['acc'] = domfile[i][21]
        domdict.append(domtmp)
        domtmp = {}
    return domdict
#print(domdict)
