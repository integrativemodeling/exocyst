#!/bin/python
import os,sys,string
chain="A"
pdb=open(sys.argv[1],'r')
pdb_lines=[line.strip('\n') for line in pdb]
pdb.close()


of=open(sys.argv[2],'w')
for ln in pdb_lines:
        ln_fields=ln.split()
        if ln_fields[4] in chain:
            of.write(ln+"\n")
	
of.close()
