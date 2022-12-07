import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import IMP.pmi.mmcif
import ihm
import IMP.bayesianem
import IMP.bayesianem.restraint
import ihm.cross_linkers
import ihm.location
import ihm.model
import RMF
import IMP.rmf
import os
import sys
import ihm.dumper
import ihm.format
import ihm.location
import ihm.representation
import ihm.startmodel
import ihm.dataset
import ihm.protocol
import ihm.analysis
import ihm.model
import ihm.restraint
import ihm.geometry


#---------------------------
# Define Input Files
#---------------------------
datadirectory = "../data/"
topology_file = datadirectory+"topology.txt"
target_gmm_file = datadirectory+'emd_21226.map.mrc.gmm.200.txt'

m = IMP.Model()

topology = IMP.pmi.topology.TopologyReader(topology_file,
                                  pdb_dir=datadirectory,
                                  fasta_dir=datadirectory,
                                  gmm_dir=datadirectory)
bs = IMP.pmi.macros.BuildSystem(m,resolutions=[1])

if '--mmcif' in sys.argv:
    po = IMP.pmi.mmcif.ProtocolOutput()
    bs.system.add_protocol_output(po)
    po.system.title = "Integrative structure and function of the yeast exocyst complex"
    po.system.citations.append(ihm.Citation.from_pubmed_id(32239688))

bs.dry_run = '--dry-run' in sys.argv

bs.add_state(topology)#keep_chain_id=True)

root_hier, dof = bs.execute_macro(max_rb_trans=4.0,
                                  max_rb_rot=0.3,
                                  max_bead_trans=4.0,
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)

fixed_particles=[]
for prot in ["Sec03","Sec05","Sec06","Sec08","Sec10","Sec15","Exo70","Exo84"]:
    fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot).get_selected_particles()

fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover,
                                          IMP.pmi.TransformMover])

IMP.pmi.tools.shuffle_configuration(root_hier,
                                    excluded_rigid_bodies=fixed_rbs,
                                    max_translation=50,
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)

outputobjects = [] 

mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname=mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
        mol, label=molname)
    cr.add_to_model()
    outputobjects.append(cr)

ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=1)
ev.add_to_model()
outputobjects.append(ev)

kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_unique_id_key("id")
kw.set_protein1_key("protein1")
kw.set_protein2_key("protein2")
kw.set_residue1_key("residue1")
kw.set_residue2_key("residue2")
kw.set_id_score_key(None)

xldb_dss = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
xldb_dss.create_set_from_file(datadirectory+'Exo.DSS.csv')

xls_dss = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=root_hier, database=xldb_dss,
        length=21, label="XLS_dss", resolution=1.0, slope=0.02,
        linker=ihm.cross_linkers.dss)

xls_dss.add_to_model()
outputobjects.append(xls_dss)
dof.get_nuisances_from_restraint(xls_dss)
xls_dss.set_psi_is_sampled(True)
psi_dss = xls_dss.psi_dictionary["PSI"][0]
psi_dss.set_scale(0.1)
dof.get_nuisances_from_restraint(xls_dss) 


xldb_edc= IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
xldb_edc.create_set_from_file(datadirectory+'Exo.EDC.csv')

xls_edc = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=root_hier, database=xldb_edc,
        length=16, label="XLS_edc", resolution=1.0, slope=0.02,
        linker=ihm.cross_linkers.edc)

xls_edc.add_to_model()
outputobjects.append(xls_edc)
dof.get_nuisances_from_restraint(xls_edc)
xls_edc.set_psi_is_sampled(True)
psi_edc = xls_edc.psi_dictionary["PSI"][0]
psi_edc.set_scale(0.1)
dof.get_nuisances_from_restraint(xls_edc)

densities = IMP.atom.Selection(root_hier,representation_type=IMP.atom.DENSITIES).get_selected_particles()
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,datadirectory+target_gmm_file,slope=0.0000001,target_radii_scale=3.0,scale_target_to_mass=True)

gem.set_label("EM")
gem.add_to_model()
outputobjects.append(gem)
num_frames = 1
num_mc_steps = 10

# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    root_hier=root_hier,
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    output_objects=outputobjects,
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=True,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=2.5,
                                    simulated_annealing_minimum_temperature_nframes=200,
                                    simulated_annealing_maximum_temperature_nframes=20,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=10,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames,
                                    global_output_directory="output",
                                    test_mode=bs.dry_run)

# Start Sampling
mc1.execute_macro()
po.finalize()
s = po.system
import ihm.dumper
with open('initial.cif', 'w') as fh:
    ihm.dumper.write(fh, [s])
for r in s.restraints:
    if isinstance(r, ihm.restraint.CrossLinkRestraint):
        print("XL-MS dataset at:", r.dataset.location.path)
        print("Details:", r.dataset.location.details)
em, = [r for r in s.restraints
       if isinstance(r, ihm.restraint.EM3DRestraint)]
d = em.dataset
dss, edc = [r for r in s.restraints
               if isinstance(r, ihm.restraint.CrossLinkRestraint)]
dss.linker = ihm.cross_linkers.dss
edc.linker = ihm.cross_linkers.edc

last_step = s.orphan_protocols[-1].steps[-1]
#print(last_step.num_models_end)
last_step.num_models_end = 200000

protocol = po.system.orphan_protocols[-1]
analysis = ihm.analysis.Analysis()
protocol.analyses.append(analysis)
analysis.steps.append(ihm.analysis.ClusterStep(
                      feature='RMSD', num_models_begin=2000000,
                      num_models_end=9669))
mg = ihm.model.ModelGroup(name="Cluster 0")

po.system.state_groups[-1][-1].append(mg)

e = ihm.model.Ensemble(model_group=mg,
                       num_models=9669,
                       post_process=analysis.steps[-1],
                       name="Cluster 0",
                       clustering_method='Other',
                       details='Density based threshold-clustering',
                       clustering_feature='RMSD',
                       precision='38'
                       )
po.system.ensembles.append(e)

Uniprot={'Sec03.0':'P33332',
         'Sec05.0':'P89102',
         'Sec06.0':'P32844',
         'Sec08.0':'P32855',
         'Sec10.0':'Q06245',
         'Sec15.0':'P22224',
         'Exo70.0':'P19658',
         'Exo84.0':'P38261'}
for prot, entry in Uniprot.items():
     ref = ihm.reference.UniProtSequence.from_accession(entry)
     po.asym_units[prot].entity.references.append(ref)

m1 = IMP.Model()
inf1 = RMF.open_rmf_file_read_only('../results/models/cluster_center_model.rmf3')
h1 = IMP.rmf.create_hierarchies(inf1, m1)[0]

for state in h1.get_children():
    comp={}
    for component in state.get_children():
            part1={}
            for i,leaf in enumerate(IMP.core.get_leaves(component)):
                p=IMP.core.XYZ(leaf.get_particle())
                part1[p.get_name()]=p.get_coordinates()
            comp[component.get_name()]=part1
del h1

for state in root_hier.get_children():
    #comp2={}
    for component in state.get_children():
            #part2={}
            for i,leaf in enumerate(IMP.core.get_leaves(component)):
                p=IMP.core.XYZ(leaf.get_particle())
                if 'Residue' in p.get_name():
                    name=p.get_name().split('_')[1]
                else:
                    name=p.get_name()
                p.set_coordinates(comp[component.get_name()][name])
                #part2[name]=p.get_coordinates()
            #comp2[component.get_name()]=part2

model = po.add_model(e.model_group)
print (e.model_group)

densityrepos = ihm.location.Repository(
    doi="10.5281/zenodo.3951752",
    url="https://zenodo.org/record/3951752/files/Localization_density.zip")

for i in ["Sec03","Sec05","Sec06","Sec08","Sec10","Sec15","Exo70","Exo84"]: 
    prot2=i+'.0'
    asym = po.asym_units[prot2]
    loc = ihm.location.OutputFileLocation('Localization_density/LPD_'+i+'.mrc',
            repo=densityrepos)
    den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
    e.densities.append(den)

dss, edc = [r for r in s.restraints
               if isinstance(r, ihm.restraint.CrossLinkRestraint)]
d = dss.dataset

scriptrepos = ihm.location.Repository(
    doi="10.5281/zenodo.3951752",
    root=".",
    top_directory="modeling",
    url="https://zenodo.org/record/3951752/files/modeling.zip")

datarepos = ihm.location.Repository(
    doi="10.5281/zenodo.3951752",
    root="../data",
    top_directory="data",
    url="https://zenodo.org/record/3951752/files/data.zip")

dcdrepos = ihm.location.Repository(
    doi="10.5281/zenodo.3951752",
    url="https://zenodo.org/record/3951752/files/cluster.0.dcd")
dcd_location = ihm.location.OutputFileLocation(path=None, repo=dcdrepos,
                                            details="cluster ensemble")

ss = ihm.model.IndependentSubsample(name='Cluster 0 subsample',num_models='9741',file=dcd_location)
e.subsamples.append(ss)

po.system.update_locations_in_repositories([scriptrepos, datarepos])

po.finalize()
with open('exocyst.cif', 'w') as fh:
    ihm.dumper.write(fh, [po.system])

import ihm.reader
with open('exocyst.cif') as fh:
    s, = ihm.reader.read(fh)
print(s.title, s.restraints, s.ensembles, s.state_groups)

