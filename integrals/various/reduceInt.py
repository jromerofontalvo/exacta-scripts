import os
import numpy
import scipy
from sys import argv

orig_file = str(argv[1])
number= str(argv[2])
limit= str(argv[3])
lista=[]
pfilename=os.path.splitext(orig_file)[0]
filename=pfilename.split("-")[-1]
mol=pfilename.split("-")[-2]
print(mol)
print(orig_file)
print(filename)
working_file=os.path.join("",mol + "-" + filename + "-" + number + ".int")

for n in range(int(number),int(limit)):
  lista += [n]

print(lista)

with open(orig_file) as f, open(working_file,"w") as output:
    for line in f:
        opt = True
        for word in line.split():
            for n in lista :
              if str(word)==str(n):
                opt = False
        if opt == True:
            output.write(line)
            
            
