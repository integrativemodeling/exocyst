import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.em
import IMP.core
import IMP.atom
import IMP.rmf
import RMF
import sys
import os
import glob
import numpy as np
from numpy.linalg import svd
import csv

###
#Reading invivo data and writing it in readable form
###
sel1=[];sel2=[];dist=[];sel3=[];sel4=[]
filename=sys.argv[1]
sf=open(filename,'r')
for i,ln in enumerate(sf.readlines()):
    line =ln.strip().split('\t')
    sel1.append(line[0].strip().split('-')[0])
    sel2.append(line[1].strip().split('_')[0])
    sel3.append(line[1].strip().split('_')[2])
    sel4.append(line[1].strip().split('_')[0]+line[1].strip().split('_')[2])
    dist.append(float(line[2])*10)

exocyst_components={'Sec03N':[1,600],'Sec03C':[601,1336],
                    'Sec05N':[1,485],'Sec05C':[486,971],
                    'Sec06N':[1,410],'Sec06C':[411,805],
                    'Sec08N':[1,473],'Sec08C':[474,1065],
                    'Sec10N':[1,490],'Sec10C':[491,871],
                    'Sec15N':[1,381],'Sec15C':[382,910],
                    'Exo70N':[1,310],'Exo70C':[311,623],
                    'Exo84N':[1,374],'Exo84C':[375,753]}


###
#Write information to file
###
tmxl1 = open('Invivo_satisfied.txt','w')
tmxl2 = open('Invivo_violated.txt','w')
tmxl3 = open('Invivo_dist_statistics.txt','w')

###
#functions
###
def get_center_of_mass(model,component):
    total_mass = 0.0
    #model_ps = [leaf.get_particle() for leaf in IMP.core.get_leaves(component)]
    com = IMP.algebra.Vector3D(0, 0, 0)
    total_mass=0.0
    for p in component:
        mass = IMP.atom.Mass(p).get_mass()
        pos = IMP.core.XYZ(p).get_coordinates()
        com += pos * mass
        total_mass += mass
    com /= total_mass
    return com

def get_center_of_mass_component(hier,component,resind1,resind2):
    total_mass = 0.0
    com = IMP.algebra.Vector3D(0, 0, 0)
    total_mass=0.0
    parts=IMP.atom.Selection(hier,molecule =component,
                                    residue_indexes=[resid1,resid2].get_selected_particles())
    for p in parts:
        mass = IMP.atom.Mass(p).get_mass()
        pos = IMP.core.XYZ(p).get_coordinates()
        com += pos * mass
        total_mass += mass
    com /= total_mass
    return com

def write_distance(sel1,sel2,sel3,sel4,com1,com2,dist,measure_dist,file):
    if measure_dist<dist:
        print >> tmxl1,i,file,sel1,sel2,sel3,sel4,measure_dist,dist
    else:
        print >> tmxl2,i,file,sel1,sel2,sel3,sel4,measure_dist,dist

def get_distance(com1,com2):
    x1=float(com1[0]);y1=float(com1[1]);z1=float(com1[2])
    x2=float(com2[0]);y2=float(com2[1]);z2=float(com2[2])
    measure_dist=((x1-x2)**2 + (y1-y2)**2 +(z1-z2)**2)**0.5
    return measure_dist

def write_statistics(dtot):
    for i,j in enumerate(dtot):
        print >> tmxl3,sel1[i],sel2[i],sel3[i],np.mean(dtot[i]),np.std(dtot[i]),min(dtot[i]),max(dtot[i]),dist[i]

def write_dtot(dtot):
    with open('Dtot_for_boxplot_invivo.csv','w') as f:
        f.write('#each row is data pertaining to 1 invivo XL from all models')
        for i in dtot:
            f.write(','.join([str(n) for n in i]))
            f.write('\n')

##

dtot=[[] for _ in range(0,len(sel1))]
for file in glob.glob("cluster.0/RMF/*.rmf3"):
    print (file)
    m = IMP.Model()
    inf = RMF.open_rmf_file_read_only(file)
    h = IMP.rmf.create_hierarchies(inf, m)[0]
    IMP.rmf.load_frame(inf, 0)
    for i in range(len(sel1)):
        p1 = IMP.atom.Selection(h,molecule =sel1[i]).get_selected_particles()
        com1=get_center_of_mass(m,p1)
        p2 = IMP.atom.Selection(h,molecule =sel2[i],residue_indexes=exocyst_components[sel4[i]]).get_selected_particles()
        com2=get_center_of_mass(m,p2)
        dtot[i].append(get_distance(com1,com2))
        write_distance(sel1[i],sel2[i],sel3[i],exocyst_components[sel4[i]],com1,com2,dist[i],get_distance(com1,com2),file)

write_statistics(dtot)
write_dtot(dtot)


