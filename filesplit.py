import sys
import subprocess
raw=[]
result=[]
readtmp=subprocess.Popen('cat /home/software/RepeatMasker/Libraries/Dfam.hmm',shell=True,stdout=subprocess.PIPE)
#print (readtmp.stdout)
for line in readtmp.stdout:
    line=str(line,'utf-8')
    raw.append(line)
    #print(raw)
    if raw[:-1] == '//\n':
        result.append(raw)
        LTR=FALSE
        raw = []
        for i in range(len(result[-1])):
            if result[-1][i] == 'CC         Type: LTR\n':
                LTR=TURE
        if LTR == FALSE:
            result.pop(-1)
print(result)
