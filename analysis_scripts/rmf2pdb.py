import IMP
import IMP.atom
import IMP.core
import IMP.display
import IMP.rmf
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import RMF
import glob
import numpy as np



color_dict = {'Sec03':(0,0,0.5),
              'Sec05':(0,0.3,0.2),
              'Sec06':(0,0.9,0.1),
              'Sec08':(0,0.9,1),
              'Sec10':(0.5,0.2,0.3),
              'Sec15':(0.5,0.1,0.9),
              'Exo70':(0.9,0,0),
              'Exo84':(1,1,0)}


def color_rmf(rmf_file):
    rmf_file = rmf_file.split('.rmf3')[0]+'.rmf3'
    m = IMP.Model()
    prots = IMP.pmi.analysis.get_hiers_from_rmf(m, 0, rmf_file)
    print(rmf_file, prots)
    prots = prots[0]
    for p,col in color_dict.iteritems():
        s=IMP.atom.Selection(prots,
                            molecule=p,
                            resolution=IMP.atom.ALL_RESOLUTIONS)
        psel=s.get_selected_particles()
        color=IMP.display.Color(col[0],col[1],col[2])
        for part in psel:
            IMP.display.Colored(part).set_color(color)

    o=IMP.pmi.output.Output()
    o.init_rmf(rmf_file,[prots])
    o.write_rmf(rmf_file)
    o.close_rmf(rmf_file)
    return rmf_file

    

def write_pdb(rmf_file):
    pdb_file = rmf_file.split('.rmf3')[0]+'.pdb'
    mm = IMP.Model()
    prots = IMP.pmi.analysis.get_hiers_from_rmf(mm, 0, rmf_file)
    print(rmf_file, prots)
    prot=prots[0]
    for p,col in color_dict.iteritems():
        s=IMP.atom.Selection(prots,
                            molecule=p,
                            resolution=IMP.atom.ALL_RESOLUTIONS)
        psel=s.get_selected_particles()
        color=IMP.display.Color(col[0],col[1],col[2])
        for part in psel:
            IMP.display.Colored(part).set_color(color)

    if not prots:
        raise ValueError("Cannot read hiearchy from rmf")
    o = IMP.pmi.output.Output()
    o.init_pdb(pdb_file, prot)
    o.write_pdb(pdb_file)
    del o, mm

color_rmf('cluster_center_model.rmf3')
write_pdb(color_rmf('cluster_center_model.rmf3 '))





