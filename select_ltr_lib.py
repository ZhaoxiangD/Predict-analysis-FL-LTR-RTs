import os
import subprocess
import sys
alldata=[]
result=[]
line=sys.stdin.read().strip('\n')
    #alldata.append(line)
data=line.split('\n')
#print(data)
for lines in range(len(data)-2,1,-2):
    dic={}
    #if data[lines][0]!='>' and data[lines-1][0]!='>':
        #data[lines-1]=data[lines-1]+data[lines]
        #data.pop[lines]
    datatmp=data[lines].split(' ')
    datatmp=datatmp[0].split('#')
    try:
        datatmp1=datatmp[1].split('/')
        if datatmp1[0]=='LTR':
            dic['ID']=data[lines]
            dic['sequence']=data[lines+1]
            result.append(dic)
    except:
        break
for i in range(len(result)):
    print(result[i]['ID'])
    print(result[i]['sequence'])
