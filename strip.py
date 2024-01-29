import sys
filename = sys.argv[1]
with open(filename) as f:
    cds = f.readlines()
file1 =[]
for rows in cds:
    if rows[0] != ">":
        rows = rows.replace('.','')
    file1.append(rows)
f = open(filename+'.rep','w')
for row in file1:
    f.write(row)
f.close()
