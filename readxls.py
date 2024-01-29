with open ('all_result.xls','r') as f:
    xlsfile=f.readlines()
xlshead = xlsfile[0].split('\t') #read the header 
xlsfile = xlsfile[1:] #read the content
xlstmp = {}
xlsdict = [] #xlsdict is a LIST, not a dict
for num in range(len(xlsfile)):
    xlsfile[num] = xlsfile[num].split('\t')
    for dictnum in range(len(xlshead)):
        xlstmp[xlshead[dictnum].strip('\n').strip(' ')] = xlsfile[num][dictnum].strip('\n').strip(' ')
    xlstmp['domID'] = xlsfile[num][0].split('_')[1].strip(' ') + ':' + xlsfile[num][2].strip(' ') #domID = ID + location
    xlsdict.append(xlstmp)
    xlstmp={}
#print(xlsdict)
