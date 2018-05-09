# Find facotrisations
# Get first non-intermediate result
t1=[]
t2=[]
f=open('gecco.itfaa', 'r')
while True:
    line = f.readline().strip().split()
    result=line[0]
    if "STIN" not in result:
        break

f.close()

f=open("gecco.itfaa","r")
for line in f:
    words=line.split()
    if "STIN" not in words[0]:
        if words[0]==result:
            t1.append(words[2])
            t2.append(words[3])
        else:
            result=words[0]

f.close()

factorise=[]
for i in range(0, 1):
    for j in range(0, len(t1)):
        if i==j:
            continue
        if t1[i]!=t2[j]:
            print(t1[i])
            print(t2[j])
        
            factorise.append(t2[j])

print(factorise)
