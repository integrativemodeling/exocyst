from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
import sys, os

env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl = complete_pdb(env, sys.argv[1])
zscore = mdl.assess_normalized_dope()
