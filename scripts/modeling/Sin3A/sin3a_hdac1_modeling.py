#################### Modeling the coREST Complex ############################
############## ------>"May the Force serve u well..." <------################
#############################################################################

############# One above all #############
##-------------------------------------##
from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
#import IMP.pmi.restraints.saxs
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.dof
import IMP.atom
#import IMP.saxs
import IMP.desmosome
import numpy as np
import os
import sys
import time

runType = sys.argv[1] # Specify test or prod
runID = sys.argv[2]   # Specify the number of runs
run_output_dir = "Output_" + str(runID)

if runType == "test":
    num_frames = 10
elif runType == "prod":
    num_frames = 30000

max_shuffle = 50
rex_max_temp = 5

# Identify data files
xls_data = "../../../input/Sin3A/XL/Sin3A_Mapped_XLs_HDAC1.csv"
pdb_dir = "../../../input/Sin3A/PDB/"
fasta_dir = "../../../input/Sin3A/Fasta/"

# Restraint weights
xl_weight = 20
mpdbr_weight = 10
mpdbr_res = 10
# em_weight = 60

# Topology File
topology_file = "../../../input/Sin3A/topology_Sin3A_HDAC1.txt"


# FUNCTIONS AND WRAPPERS -----------------------------------------------
# wrapper for the MPDBR Restraint
class MinimumPairDistanceBindingRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist1, plist2, x0=0, kappa=1, label=None, weight=1):
        name = 'MinimumPairDistanceBindingRestraint%1%'
        super(MinimumPairDistanceBindingRestraint, self).__init__(model, name=name, label=label, weight=weight)
        l1 = IMP.container.ListSingletonContainer(mdl)
        l1.add(plist1)
        l2 = IMP.container.ListSingletonContainer(mdl)
        l2.add(plist2)
        bipartite_container = IMP.container.AllBipartitePairContainer(l1, l2)
        score = IMP.core.HarmonicSphereDistancePairScore(x0, kappa)
        res_main = IMP.container.MinimumPairRestraint(score, bipartite_container, 1)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added MPDBR on particles", len(plist1), ":", len(plist2), "at x0:kappa", str(x0), ":",
              str(kappa))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All IMP systems start out with a Model
mdl = IMP.Model()

# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file)
t = IMP.pmi.topology.TopologyReader(topology_file,
                                  pdb_dir=pdb_dir,
                                  fasta_dir=fasta_dir)
                                  # gmm_dir=gmm_dir)


# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)

# executing the macro will return the root hierarchy and degrees of freedom (dof) objects
root_hier, dof = bs.execute_macro(max_rb_trans = 1,
                                  max_rb_rot = 0.2,
                                  max_bead_trans = 4.5,
                                  max_srb_trans = 1,
                                  max_srb_rot = 0.05)


# It's useful to have a list of the molecules.
molecules = t.get_components()


##### Uncomment the following lines to get test.rmf file to visualise the system representation
#
# # Uncomment this line for verbose output of the representation
# IMP.atom.show_with_representations(root_hier)
# # output to RMF
# fname = 'test.rmf'
# rh = RMF.create_rmf_file(fname)
# IMP.rmf.add_hierarchy(rh, root_hier)
# IMP.rmf.save_frame(rh)




#####################################################
##################### RESTRAINTS ####################
#####################################################

# Restraints define functions that score the model based on
# input information.
#
# Restraint objects are first created in the definition.
# To be evaluated, the restraint object must be add_to_model().
#
# In some cases, sampled parameters for restraints must be added to the DOF
# object

# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.
output_objects = []

# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT
#
# Restrains residues/particles that are connected in sequence
# This should be used for any system without an atomic force field (e.g. CHARMM)
# We apply the restraint to each molecule

for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

print("Connectivity restraint applied")


# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of all molecule chains to included_objects to apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=[root_hier],
                                            resolution=1000)
output_objects.append(evr)

print("Excluded volume restraint applied")


# -------------------------
# %%%%% CROSSLINKING RESTRAINT
#
# Restrains two particles via a distance restraint based on
# an observed crosslink.
#
# First, create the crosslinking database from the input file
# The "standard keys" correspond to a crosslink csv file of the form:
#
# Protein1,Residue1,Protein2,Residue2
# A,18,G,24
# A,18,G,146
# A,50,G,146
# A,50,G,171
# A,50,G,189
#
# This restraint allows for ambiguity in the crosslinked residues,
# a confidence metric for each crosslink and multiple states.
# See the PMI documentation or the MMB book chapter for a
# full discussion of implementing crosslinking restraints.

# This first step is used to translate the crosslinking data file.
# The KeywordsConverter maps a column label from the xl data file
# to the value that PMI understands.
# Here, we just use the standard keys.
# One can define custom keywords using the syntax below.
# For example if the Protein1 column header is "prot_1"
# xldbkc["Protein1"]="prot_1"

# The CrossLinkDataBase translates and stores the crosslink information
# from the file "xl_data" using the KeywordsConverter.

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

xl = IMP.pmi.io.crosslink.CrossLinkDataBase()
xl.create_set_from_file( file_name = xls_data,
                              converter = xldbkc )
xls = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier = root_hier,    # Must pass the root hierarchy to the system
                database = xl,            # The crosslink database.
                length = 25,              # The crosslinker plus side chain length
                resolution = 1,           # The resolution at which to evaluate the crosslink
                slope = 0.0001,           # This adds a linear term to the scoring function
                label = "xl_Sin3A_hdac1",                        #   to bias crosslinks towards each other
                weight = xl_weight)       # Scaling factor for the restraint score.

output_objects.append(xls)

print("Cross-linking restraint applied")


# -------------------------------------------------
# %%%%% MIN PAIR DISTANCE BINDING RESTRAINT (MPDBR)
# SAP30-HDAC1 interacting residues.

plist1 = []
for res1 in [84, 86, 88, 110]:
    selection_tuple = (res1, res1, 'SAP30', 0, None)
    p1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, mpdbr_res)
    if p1 not in plist1:
        plist1 = np.append( plist1, p1 )

plist2 = []
for res2 in [31, 33, 270, 306]:
    selection_tuple = (res2, res2, 'HDAC1', 0, None)
    p2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, mpdbr_res)
    if p2 not in plist2:
        plist2 = np.append( plist2, p2 )

mpdbr = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, x0 = 0, kappa = 1, label = f"SAP30.0-HDAC1.0",
                                                      weight = mpdbr_weight)
output_objects.append( mpdbr )
print( f"Min Pair Distance Binding Restraint Applied --> SAP30.0-HDAC1.0" )


# SIN3A-HDAC1 interacting residues.
selection_tuple = ( 687, 829, 'SIN3A', 0, None )
plist1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, mpdbr_res)

selection_tuple = ( 8, 376, 'HDAC1', 0, None )
plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, mpdbr_res)

mpdbr = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, x0 = 0, kappa = 1, label = f"SIN3A.0-HDAC1.0",
                                                      weight = mpdbr_weight)
output_objects.append( mpdbr )
print( f"Min Pair Distance Binding Restraint Applied --> SIN3A.0-HDAC1.0" )


#####################################################
###################### SAMPLING #####################
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.
print("The type of run is: " + str(runType))
print("Number of sampling frames: " + str(num_frames))
# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation

IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation = max_shuffle)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(1500)


# Now, add all of the other restraints to the scoring function to start sampling
for o in output_objects:
    o.add_to_model()

print( "Replica Exchange Maximum Temperature : " + str( rex_max_temp ) )

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier = root_hier,                    # pass the root hierarchy
        crosslink_restraints = xls,
        # This allows viewing the crosslinks in Chimera. Also, there is not inter-protein ADH crosslink available. Hence it is not mentioned in this list
        monte_carlo_temperature = 1.0,
        replica_exchange_minimum_temperature = 1.0,
        replica_exchange_maximum_temperature = rex_max_temp,
	    monte_carlo_sample_objects = dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
        global_output_directory = run_output_dir,      # The output directory for this sampling run.
        output_objects = output_objects,          # Items in output_objects write information to the stat file.
        monte_carlo_steps = 10,                   # Number of MC steps between writing frames
        number_of_best_scoring_models = 0,        # set >0 to store best PDB files (but this is slow)
        number_of_frames = num_frames )            # Total number of frames to run / write to the RMF file.
        # test_mode=test_mode)                    # (Ignore this) Run in test mode (don't write anything)

# Ok, now we finally do the sampling!
rex.execute_macro()
