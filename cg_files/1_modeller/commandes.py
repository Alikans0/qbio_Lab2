# Very fast homology modeling by the automodel class
from modeller import *
from modeller.automodel import *    # Load the automodel class

# Redefine the special_patches routine to include the additional disulfides
# (this routine is empty by default):
#class mymodel(automodel):
    #def special_patches(self, aln):
         #A disulfide between residues 115 and 197:
        #self.patch(residue_type='DISU', residues=(self.residues['91'],
                                                  #self.residues['174']))

log.verbose()
env = environ()
# directories for input atom files
env.io.atom_files_directory = './:../atom_files'
#my remplace auto
a = automodel(env,
              alnfile='ali.ali',      # alignment filename
              knowns='7f9z',                # codes of the templates
              sequence='ghsr',
	      assess_methods=(assess.DOPE,                 
                              assess.GA341))              # code of the target
              
a.starting_model = 1
a.ending_model = 20
# a.final_malign3d = True

a.make()                            # make the homology model

