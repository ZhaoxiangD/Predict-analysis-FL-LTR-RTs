import sys
import os
import time
import subprocess
'''
----2021/3/20----
标准输入ltr-finder输出的数据
输出去除无ltr转座子的序列信息  #未启用
----2021/3/21----
新增delete_repeat_ltr函数，功能：删除重叠的ltr
新增make_dic函数，功能：
----2021/3/22----
新增get_samtools函数，功能：调用samtools获取对应的序列信息
----2021/3/23----
get_samtools函数追加功能：调用tandem repeats finder筛选ssr占比超过35%的假ltr
----2021/3/24----
增加chromosome变量
----2021/3/26----
增加get_RepeatMasker函数，功能：调用RepeatMasker确定所属家族
----2021/4/10----
create function: rm_noneltr  # no more error caused by 'No LTR Retrotransposons Found\n'
'''
chromosome=sys.stdin.read().strip('\n')

def delete_none_ltr(result):
    #去除无ltr转座子的序列信息
    for i1 in range(len(result)-1,0,-1): #倒叙防止出错
        for i2 in range(len(result[i1])):
            if result[i1][i2]=='No LTR Retrotransposons Found\n':
                result.pop(i1)
    return result

def make_dic(result,chromosome):
    #建立字典，便于搜寻
    result_dict={}
    result_processed=[]
    for i in range(3,len(result[1])):#第二个列表的第4行开始，前几行为注释
        line=result[1][i].split( ) #以空格为分隔符，可能要改为\t
        #for i1 in range(6,len(line)): 
            #if line[i1] == 'N-N':
                #line[i1] = 0
            #else:
                #line[i1] = 1
        if line[0] == '[' or line[0] == ']':
            line.pop(0)
            line[0] = '['+'0'+line[0]
        line[0]=line[0][1:-1]
        result_dict['ID'] = line[0] + '_' + line[1] #ID为index+SeqID
        result_dict['SeqID'] = line[1]
        result_dict['location'] = line[2]
        result_dict['LTR_len'] = line[3]
        result_dict['Inserted_element_len'] = line[4]
        result_dict['TSR'] = line[5]
        result_dict['PBS'] = line[6]
        result_dict['PPT'] = line[7]
        result_dict['RT'] = line[8]
        result_dict['IN(core)'] = line[9]
        result_dict['IN(c-term)'] = line[10]
        result_dict['RH'] = line[11]
        result_dict['Strand'] = line[12]
        result_dict['Score'] = line[13]
        result_dict['Sharpness'] = line[14]
        result_dict['Similarity'] = line[15]
        result_dict['chr'] = chromosome
        result_processed.append(result_dict)
        result_dict={}
    return result_processed

def delete_repeat_ltr(result_processed):
    #除去重复LTR
    for i in range(len(result_processed)-2,0,-1):
        #print(result_processed[i]['Similarity'])
        location=result_processed[i]['location'].split('-')
        location_previous=result_processed[i-1]['location'].split('-')
        location_after=result_processed[i+1]['location'].split('-')
        if int(location[0]) <= int(location_previous[-1]) and int(location[0]) >= int(location_previous[0]):
            if float(result_processed[i]['Score']) > float(result_processed[i-1]['Score']):
                result_processed.pop(i-1)
            elif float(result_processed[i]['Score']) < float(result_processed[i-1]['Score']):
                result_processed.pop(i)
            elif float(result_processed[i]['Score']) == float(result_processed[i-1]['Score']):
                if float(result_processed[i]['Similarity']) > float(result_processed[i-1]['Similarity']):#[:-2]除去/n
                    result_processed.pop(i-1)
                else:
                    result_processed.pop(i)
        elif int(location[-1]) <= int(location_after[-1]) and int(location[-1]) >= int(location_after[0]):
            if float(result_processed[i]['Score']) > float(result_processed[i+1]['Score']):
                result_processed.pop(i+1)
            elif float(result_processed[i]['Score']) < float(result_processed[i+1]['Score']):
                result_processed.pop(i)
            elif float(result_processed[i]['Score']) == float(result_processed[i+1]['Score']):
                if float(result_processed[i]['Similarity']) > float(result_processed[i+1]['Similarity']):
                    result_processed.pop(i+1)
                else:
                    result_processed.pop(i)
    return result_processed

def get_samtools(result_final):
    #调用samtools获取ltr序列信息，并筛选ssr比例小于35%的ltr
    for i in range(len(result_final)):
        sam = subprocess.Popen('samtools faidx ./{2}.fna {0}:{1}'.format(result_final[i]['SeqID'],result_final[i]['location'],result_final[i]['chr']), shell=True, stdout=subprocess.PIPE)#popen能返回信息
        sam.wait()
        tmp = sam.stdout.readlines()
        #subprocess.Popen('mkdir {0}'.format(result_final[i]['chr']), shell = True)
        sam_file=open('{0}_{1}.fa'.format(result_final[i]['ID'],result_final[i]['chr']),'w')
        for i1 in range(len(tmp)):
            tmp[i1]=str(tmp[i1],'utf-8')#标准输出结果是字节，转换成字符
            sam_file.write(tmp[i1])
        sam_file.close()
        trf=subprocess.Popen('trf {0}_{1}.fa 2 7 7 80 10 50 500 -h'.format(result_final[i]['ID'],result_final[i]['chr']), shell=True, stdout=subprocess.PIPE)#调用tandem repeats finder
        trf.wait()
        trf_waste = trf.stdout.readline()
        trf_waste = ''
        ssr = subprocess.Popen('cat {0}_{1}.fa.2.7.7.80.10.50.500.dat'.format(result_final[i]['ID'],result_final[i]['chr']),shell=True, stdout=subprocess.PIPE)
        ssr.wait()
        ssr_file = ssr.stdout.readlines()
        repeat_num = 0
        for i2 in range(15,len(ssr_file)):
            ssr_file[i2]=str(ssr_file[i2],'utf-8')#标准输出结果是字节，转换成字符
            ssr_file[i2]=ssr_file[i2].split( )
            #repeat_num += int(ssr_file[i2][2])
            repeat_num += float(ssr_file[i2][4])
            repeat_num = repeat_num * float(ssr_file[i2][3])#统计重复序列长度
        result_final[i]['ssr']=repeat_num
        #删除ssr占比超过35%的ltr ！！！已被替代！！！
        #if repeat_num/int(result_final[i]['Inserted_element_len']) > 0.35: #若ssr占比>35%认为是假ltr
            #os.system('rm {0}_chr4.fa'.format(result_final[i]['ID']))
            #os.system('rm {0}_chr4.fa.2.7.7.80.10.50.500.dat'.format(result_final[i]['ID']))
            #result_final.pop[i]
            #result_final[i]['ID']='Na'
        #else:
        os.system('rm {0}_{1}.fa.2.7.7.80.10.50.500.dat'.format(result_final[i]['ID'],result_final[i]['chr']))
    #删除ssr占比超过35%的ltr
    for i in range(len(result_final)-1,-1,-1):
        if int(result_final[i]['ssr']) / int(result_final[i]['Inserted_element_len']) >= 0.35:
            os.system('rm {0}_{1}.fa'.format(result_final[i]['ID'],result_final[i]['chr']))
            result_final.pop(i)
    return result_final

def get_RepeatMasker(result_final):
#调用repeatmasker获取ltr所属家族，并删除所有无所属家族的ltr？
    for i in range(len(result_final)-1,-1,-1):
        #调用RepeatMasker
        rm=subprocess.Popen('/home/software/RepeatMasker/RepeatMasker {0}_{1}.fa -low -lib /home/software/RepeatMasker/Libraries/RepeatMasker.lib '.format(result_final[i]['ID'],result_final[i]['chr']), shell=True, stdout=subprocess.PIPE)
        rm.wait()
        rmtmp = rm.stdout.readline()
        rmtmp = ''
        #调出repeatmasker生成的文件内容，out文件即为所需文件
        rm_file_tmp = subprocess.Popen('cat {0}_{1}.fa.out'.format(result_final[i]['ID'],result_final[i]['chr']), shell=True, stdout=subprocess.PIPE)
        rm_file_tmp.wait()
        rm_file=rm_file_tmp.stdout.readlines()
        result_final[i]['family']=''
        #确定ltr所属家族并删除所有无确定家族的ltr
        if len(rm_file) > 1:#文件行数小于1的为未找到所属家族的ltr
            for i1 in range(3,len(rm_file)):#从第三行开始才为有效信息
                rm_file[i1]=str(rm_file[i1],'utf-8')
                rm_file[i1]=rm_file[i1].split(' ')
                #以空格分割后会产生大量'',需去除以方便后续操作
                for i2 in range(len(rm_file[i1])-1,-1,-1):
                    if rm_file[i1][i2]=='':
                        rm_file[i1].pop(i2)
                #结果可能会出现重复，去除重复结果
                if result_final[i]['family'] != ' '+ rm_file[i1][-6]:
                    result_final[i]['family']=result_final[i]['family'] +' '+ rm_file[i1][-6]
            #删除repeatmasker产生的文件
            os.system('rm -r {0}_{1}.fa.*'.format(result_final[i]['ID'],result_final[i]['chr']))
        else:
            #删除无确定家族的ltr文件和dic
            #os.system('rm -r {0}_{1}.txt*'.format(result_final[i]['ID'],result_final[i]['chr']))
            #result_final.pop(i)
            #无结果的ltr的family定义为Unclassified
            result_final[i]['family']=' Unclassified'
            os.system('rm -r {0}_{1}.fa.*'.format(result_final[i]['ID'],result_final[i]['chr']))
    return result_final

def rm_repeat_family(result_final):
    for i in range(len(result_final)):
        family=''
        result_final[i]['family']=result_final[i]['family'].split(' ')
        result_final[i]['family']=sorted(result_final[i]['family'])[1:]
        if len(result_final[i]['family']) > 1:
            for i1 in range(len(result_final[i]['family'])-1,0,-1):
                if result_final[i]['family'][i1]==result_final[i]['family'][i1-1]:
                    result_final[i]['family'].pop(i1)
        for i2 in range(len(result_final[i]['family'])):
            family=family+' '+result_final[i]['family'][i2]
        result_final[i]['family']=family[1:]
    return result_final

def delete_dna(result_final):
    for i in range(len(result_final)-1,-1,-1):
        family=result_final[i]['family'].split(' ')
        if len(family)==1 and family[0]=='Unclassified': # delete unclassified canlender
            result_final.pop(i)
        else:
            for num in range(len(family)-1,-1,-1):
                family_tmp=family[num].split('/')
                if family_tmp[0] != 'LTR': # delete all family except LTR
                    family.pop(num)
            if len(family)<1 : # if family do not include ltr
                result_final.pop(i)
            else:
                result_final[i]['family']=' '.join(family)
    return result_final

def mv2dir(chromosome):
    md = subprocess.Popen('mkdir {0}'.format(chromosome), shell = True, stdout = open('/dev/null','w'))
    md.wait()
    mv = subprocess.Popen('mv *{0}.fa {0}'.format(chromosome), shell=True, stdout = open('/dev/null','w'))
    mv.wait()
    return

def rm_noneltr(result):
    ltr = 1
    if result[1][-1] == 'No LTR Retrotransposons Found\n':
        ltr=0
    return ltr
 
raw=[]
n=0
result=[]
readtmp=subprocess.Popen('cat {0}_finder_result.txt'.format(chromosome),shell=True,stdout=subprocess.PIPE)
for line in readtmp.stdout:
    line=str(line,'utf-8')
    raw.append(line)
    if raw[-1] == '\n': #以空行为分re隔符，将两个空行内的数据单独放入列表中
        result.append(raw[:-1])
        raw=[]
ltr_fact = rm_noneltr(result)
if ltr_fact == 0:
    print('no LTR detected')
else:
    result_processed=make_dic(result,chromosome)
    result_final=delete_repeat_ltr(result_processed)
    result_final=get_samtools(result_final)
    result_final=get_RepeatMasker(result_final)
    result_final=rm_repeat_family(result_final)
    result_final=delete_dna(result_final)
    mdir = mv2dir(chromosome)
    print('ID','\t','chr','\t','location','\t','family','\t','LTR_len','\t','Inserted_element_len','\t','TSR','\t','PBS','\t','PPT',\
        '\t','RT','\t','IN(core)','\t','IN(c-term)','\t','RH','\t','Strand','\t','Score','\t',\
            'Sharpness','\t','Similarity')
    for i in range(len(result_final)):
        print(result_final[i]['ID'],'\t',result_final[i]['chr'],'\t',result_final[i]['location'],'\t',result_final[i]['family'],'\t',result_final[i]['LTR_len'],'\t',result_final[i]['Inserted_element_len'],'\t',\
            result_final[i]['TSR'],'\t',result_final[i]['PBS'],'\t',result_final[i]['PPT'],'\t',\
                result_final[i]['RT'],'\t',result_final[i]['IN(core)'],'\t',result_final[i]['IN(c-term)'],\
                    '\t',result_final[i]['RH'],'\t',result_final[i]['Strand'],'\t',result_final[i]['Score'],\
                        '\t',result_final[i]['Sharpness'],'\t',result_final[i]['Similarity'])
