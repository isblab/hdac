# Homology modeling by the automodel class
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

log.level(1)     # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
#env.io.atom_files_directory = ['.', '../atom_files']

a = automodel(env,
              alnfile  = 'SUDS3_205-228-SIN3A_604-728.ali',     # alignment filename
              knowns   = ('2n2h_AB',),              # codes of the templates
              sequence = 'SUDS3_205-228-SIN3A_604-728',
	      assess_methods=assess.normalized_dope)              # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 20                 # index of the last model
                                    # (determines how many models to calculate)
a.make()                            # do the actual homology modeling

