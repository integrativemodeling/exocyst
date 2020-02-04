#!/bin/python
import os,sys,string
chain="D"
chaino="X"
pdb=open(sys.argv[1],'r')
pdb_lines=[line.strip('\n') for line in pdb]
pdb.close()


of=open(sys.argv[2],'w')
for ln in pdb_lines:
        ln_fields=ln.split()
        res=ln_fields[3]+ " " + chaino
	res1=ln_fields[3]+ " " + chain
        print(res)
        ln=ln.replace(str(res),str(res1))
        of.write(ln+"\n")
	
of.close()


