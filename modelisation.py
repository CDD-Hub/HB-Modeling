#! -*- coding: utf-8 -*-

#KC#CG# packages import
import os
import shutil
import random as rd
from modeller import *
from modeller.automodel import *

#CG# set all Modeller log output levels
log.level(output=0, notes=0, warnings=0, errors=0, memory=0)

#KC#CG# own files import
from model import *
from mymodel import *
from fonctions import *

#KC# working directory
workdir=os.getcwd()

########################################################################
######################                            ######################
######################        modelisation        ######################
######################                            ######################
######################################################################## --> CG

class modelisation:
    """ KC - class for modelisation """
    
    #KC#CG# class "constructor"
    def __init__(self, job_name, sequence_name, pdb_template=None, seq_model=None, seq_template=None, pdb_ref=None, n_models=1, sequence_filename=None, restraints=[], actions=[], blank_model=None, terminus_folding=False):
        
        #KC#CG# parameters tests        
        assert isinstance(sequence_name, str), "ERROR: sequence_name have to be a string"
        assert isinstance(sequence_filename, str), "ERROR: sequence_filename have to be a string"
        assert isinstance(pdb_template, str), "ERROR: pdb_template have to be a string"
        assert (isinstance(pdb_ref, str) or pdb_ref is None), "ERROR: pdb_ref have to be a string"
        assert (isinstance(n_models, int) or n_models>1), "ERROR: n_models have to be an integer > 1"
        assert isinstance(restraints, list), "ERROR: restraints have to be a list"
        assert isinstance(actions, list), "ERROR: actions have to be a list"
        assert os.path.isfile(sequence_filename), 'ERROR: sequence file '+sequence_filename+' not exist'
        assert os.path.isfile(pdb_template+'.pdb'), 'ERROR: sequence file '+pdb_template+' not exist'
        if pdb_ref is not None:
            assert os.path.isfile(pdb_ref+'.pdb'), 'ERROR: sequence file '+pdb_ref+' not exist'

        #KC#CG# variables declaration
        self.blank_model=blank_model
        self.sequence_name=sequence_name
        self.sequence_filename=sequence_filename
        self.pdb_ref=pdb_ref
        self.pdb_template=pdb_template
        self.n_models=n_models
        self.job_name=job_name
        self.actions=actions #CG# actions before modelling saved
        self.terminus_folding=terminus_folding
        self.restraints=restraints #CG# tab with all restraints to apply in the modelisation [[type, dist, stddev, residue/atom1, residue/atom2], ['SS', 2.2, 0.01, 45, 123], [...], [...], ...]
        self.seq_model=seq_model
        self.seq_template=seq_template
        
    ##################################################
    ####################### run ######################
    ################################################## --> CG
        
    def run(self):
        """ KC - run of a modelisation """
        
        print "JOB =  "+self.job_name
        self.models=[]
        
        #KC#CG# "job name" folder created and alignment file moved into
        name=self.job_name
        os.chdir(name)
        shutil.copyfile("../res_final.ali", "./resultat.ali") #KC# alignment between template and blank
        
        #KC#CG# spaces replaced in string name
        assert isinstance(name, str), 'ERROR: name argument in _modelisation have to be a string'
        name=name.replace(' ', '_')
        
        #CG# random seed
        env=environ(rand_seed=rd.randrange(-50000,-1))
        #KC#CG# directories for input atom files
        env.io.atom_files_directory='./'+name+'/.'
        #KC#CG# MyModel launching wich add restraints to the default ones
        a=MyModel(env, alnfile='resultat.ali', knowns="blank_model", sequence=self.sequence_name, assess_methods=(assess.DOPE, assess.GA341), parent=self, terminus_folding=self.terminus_folding)
        a.starting_model=1
        a.ending_model=self.n_models
        a.make()

        #CG# get the best model from the n realized by Modeller
        for model in a.outputs:
            if model['failure'] is None:
                self.models.append(my_model(filename=os.getcwd()+"/"+model["name"], seq_model=self.seq_model, seq_template=self.seq_template, sequence_name=self.sequence_name, \
                    list_helix=self.blank_model.res_list_helix, DOPE=model["DOPE score"], molpdf=model["molpdf"], GA341=model["GA341 score"], ASA_per=self.blank_model.ASA_per))
        clean_tempfile()
        os.chdir('../.')
        
        print "End of \""+self.job_name+"\" job"
        
        return self.models