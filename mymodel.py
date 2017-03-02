#! -*- coding: utf-8 -*-

#KC#CG# packages import
from modeller import *
from modeller.automodel import *

#CG# set all Modeller log output levels
log.level(output=0, notes=0, warnings=0, errors=0, memory=0)

#KC#CG# own files import
from fonctions import *
from modelisation import *

########################################################################
#########################                       ########################
#########################        MyModel        ########################
#########################                       ########################
######################################################################## --> CG

class MyModel(automodel):
    
    def __init__(self,env, alnfile, knowns, sequence, assess_methods, parent, terminus_folding=False):

        automodel.__init__(self, env=env, alnfile=alnfile, knowns=knowns, sequence=sequence, assess_methods=assess_methods)
        self.parent=parent
        self.terminus_folding=terminus_folding
                                       
    def special_restraints(self, aln):
        
        #CG# restraint format: [type, distance, stdev, atom1, atom2]
        rsr=self.restraints
        at=self.atoms
        restraints=self.parent.restraints

        for restraint in restraints:
            restraint_type=restraint[0]
            restraint_atom1=restraint[1]
            restraint_atom2=restraint[2]
            if restraint_type == "distance":
                restraint_distance=float(restraint[3])
                restraint_stdev=float(restraint[4])
                rsr.add(forms.gaussian(group=physical.xy_distance, feature=features.distance(at[restraint_atom1], at[restraint_atom2]), mean=restraint_distance, stdev=restraint_stdev))
            elif restraint_type == "upperdistance":
                restraint_distance=float(restraint[3])
                restraint_stdev=float(restraint[4])
                rsr.add(forms.upper_bound(group=physical.xy_distance, feature=features.distance(at[restraint_atom1], at[restraint_atom2]), mean=restraint_distance, stdev=restraint_stdev))
            elif restraint_type == "patch_DISU":
                self.patch(residue_type='DISU', residues=(self.residues[restraint_atom1], self.residues[restraint_atom2]))
            elif restraint_type == "helix":
                rsr.add(secondary_structure.alpha(self.residue_range(restraint_atom1, restraint_atom2)))
            elif restraint_type == "sheet":
                rsr.add(secondary_structure.sheet(self.residue_range(restraint_atom1, restraint_atom2)))
            elif restraint_type == "strand":
                rsr.add(secondary_structure.strand(self.residue_range(restraint_atom1, restraint_atom2)))
            else :
                print  "WARNING: Unknown restraint type"
        
        #CG# restraints for terminus added
        if self.terminus_folding:
            center=pseudo_atom.gravity_center(self.atoms)
            top_center=pseudo_atom.gravity_center([atom for atom in self.atoms if atom.residue.index in self.parent.blank_model.top])
            down_center=pseudo_atom.gravity_center([atom for atom in self.atoms if atom.residue.index in self.parent.blank_model.down])
            self.restraints.pseudo_atoms.append(center) 
            self.restraints.pseudo_atoms.append(down_center) 
            self.restraints.pseudo_atoms.append(top_center)    
            
            for atom in self.atoms:
                #CG# restraints added
                if atom.residue.index in self.parent.blank_model.C_ter:
                    rsr.add(forms.upper_bound(group=physical.xy_distance,feature=features.distance(atom,center),mean=self.parent.blank_model.z_max*1.1, stdev=1))
                                                
                for i,idx in enumerate(self.parent.blank_model.N_ter):
                    if idx==atom.residue.index:
                        rsr.add(forms.upper_bound(group=physical.xy_distance,feature=features.distance(atom,top_center),mean=10, stdev=1))
