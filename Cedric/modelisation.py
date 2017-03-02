#! -*- coding: utf-8 -*-

#KC#CG# packages import
import os
import shutil
import random as rd
from pylab import *
from modeller import *
from modeller.automodel import *
from modeller.parallel import *

#CG# set all Modeller log output levels
log.level(output=0, notes=0, warnings=0, errors=0, memory=0)

#KC#CG# own files import
from model import *
from mymodel import *
from fonctions import *

########################################################################
######################                            ######################
######################        modelisation        ######################
######################                            ######################
########################################################################

class modelisation:
    """ KC - class for modelisation """
    
    #KC#CG# class "constructor"
    def __init__(self, job_name, sequence_name, pdb_template, seq_model, seq_template, pdb_ref=None, n_models=10, sequence_filename='seq.ali', restraints=[], actions=[], blank_model=None,terminus_folding=False):
        
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
        self.restraints=restraints #CG# tab with all restraints to apply in the modelisation [ [type, dist, stddev, residue/atom1, residue/atom2], ['SS', 2.2, 0.01, 45, 123], [...], [...], ...]
        self.seq_model=seq_model
        self.seq_template=seq_template
        self.models=None
        self.atom_models_list=None
        self.residue_models_list=None

        #KC#CG# functions launching
        self.run()                   

    ##################################################
    ############ get and set functions ###############
    ##################################################
    
    def get_pdb_template(self):
        return self.pdb_template
    
    def get_n_models(self):
        return self.n_models
        
    def get_restraints(self):
        return self.restraints
                
    def set_pdb_template(self, pdb_template):
        if isinstance(pdb_template, str):
            self.color_conservation(self).pdb_template=pdb_template
        else:
            print 'ERROR: set_pdb_template argument have to be a string'
    
    def set_n_models(self, n_models):
        if isinstance(n_models, int):
            if n_models>1:
                self.n_models=models
            else:
                print 'ERROR: set_n_models argument have to be an integer > 1'
        else:
            print 'ERROR: set_n_models argument have to be an integer > 1'
            
    def set_restraints(self, restraints_list):
        if isinstance(restraints_list, list):
            self.restraints=restraints_list
        else:
                print 'ERROR: set_restraints argument have to be a list'

    ##################################################
    ####################### run ######################
    ##################################################
        
    def run(self):
        """ KC - run of a modelisation """
        
        print "JOB =  "+self.job_name
        
        #KC#CG# "job name" folder created and alignment file moved into
        name=self.job_name
        os.chdir(name)
        shutil.copyfile("../res_final.ali", "./resultat.ali") #KC# alignment between template and blank
        
        #KC#CG# spaces replaced in string name
        assert isinstance(name, str), 'ERROR: name argument in _modelisation have to be a string'
        name=name.replace(' ', '_')
        
        env=environ(rand_seed=rd.randrange(-50000,-1)) #CG# random seed
        env.io.atom_files_directory='./'+name+'/.'
        
        #CG# j = job()
        #CG# j.append(local_slave())
        #CG# j.append(local_slave())
        #CG# j.append(local_slave())
        #CG# j.append(local_slave())
        #CG# j.append(local_slave())
        
        a=MyModel(env, alnfile='resultat.ali', knowns="blank_model", sequence=self.sequence_name, assess_methods=(assess.DOPE, assess.GA341), parent=self, terminus_folding=self.terminus_folding)
        a.starting_model=1
        a.ending_model=self.n_models
        
        #CG# Very thorough VTFM optimization:
        #CG# a.library_schedule = autosched.slow
        #CG# a.max_var_iterations = 400

        #CG#  Thorough MD optimization:
        #CG# a.md_level = refine.very_slow

        #CG#  Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
        #CG# a.repeat_optimization = 2
        #CG# a.max_molpdf = 1e6
        
        #CG# a.use_parallel_job(j)

        a.make()

        #CG# get the DOPE score
        dic={}
        """ CG - model example
        {'molpdf': 1024.6275634765625, 'name': '5HT2B.B99990001.pdb',
        'pdfterms': physical.values(bond=11.132513, angle=220.248734, dihedral=406.517822,
        improper=32.589676, soft_sphere=6.329496, lennard_jones=0.000000, coulomb=0.000000,
        h_bond=0.000000, ca_distance=52.086723, n_o_distance=33.452244, phi_dihedral=0.000000,
        psi_dihedral=0.000000, omega_dihedral=66.401932, chi1_dihedral=33.396786, chi2_dihedral=62.984882,
        chi3_dihedral=31.286905, chi4_dihedral=13.019203, disulfide_distance=0.009089,
        disulfide_angle=0.666614, disulfide_dihedral=1.290730, lower_distance=0.000000,
        upper_distance=0.000000, sd_mn_distance=109.917717, chi5_dihedral=0.000000,
        phi_psi_dihedral=-171.499817, sd_sd_distance=114.796341, xy_distance=0.000000,
        nmr_distance=0.000000, nmr_distance2=0.000000, min_distance=0.000000, nonbond_spline=0.000000,
        accessibility=0.000000, density=0.000000, absposition=0.000000, dihedral_diff=0.000000,
        gbsa=0.000000, em_density=0.000000, saxs=0.000000, symmetry=0.000000),
        'GA341 score': [1.0, 0.2325633019208908, -165.8031005859375, 15.074830055236816,
        -1.0837429761886597, -4.53900146484375, -2.0401101112365723, -4.477954864501953],
        'failure': None, 'num': 1, 'DOPE score': -37603.42578125}
        """
        
        for model in a.outputs:
            if model['failure'] is None:
                dic[model["name"]]=model

        self.models=[]
        for name in dic:
            pdb_filename=name
            #KC#CG# pdb_code=pdb_filename[:-4]
            self.models.append(my_model(os.getcwd()+"/"+pdb_filename, self.seq_model, self.seq_template, self.sequence_name, self.blank_model.res_list_helix, DOPE=dic[name]["DOPE score"], molpdf=dic[name]["molpdf"]))
        self.structure_alignement()
        self.models_analysis()
        os.chdir('../.')
        
        print "End of \""+self.job_name+"\" job"

    ##################################################
    ############## structure alignement ##############
    ##################################################
                
    def structure_alignement(self):
        """ KC - loop on models to do structure alignment """
        
        #CG# alignement based on helix
        for model in self.models[1:]:
            model.structure_alignement(self.models[0])
            
    ##################################################
    ################# models analysis ################
    ##################################################
            
    def models_analysis(self):
        """ KC - analyse of models """
        
        self.ASA_calcul()
        self.dihedral_calcul()
        self.RMSD_calcul()
        self.color_conservation()
        self.analyse_contacts()

    ##################################################
    ################### ASA calcul ###################
    ##################################################
            
    def ASA_calcul(self):
        """ KC - loop on models to calculate ASA """
        
        for model in self.models:
            model.ASA_calcul()
    
    ##################################################
    ################# dihedral calcul ################
    ##################################################
    
    def dihedral_calcul(self):
        """ KC - loop on models to calculate dihedral angles """
        
        self.models[0].dihedral_calcul(self.models[1:])
        
    ##################################################
    ################## RMSD calcul ###################
    ##################################################
    
    def RMSD_calcul(self):
        """ KC - loop on models to calculate RMSD """
        
        self.models[0].RMSD_calcul(self.models[1:])
        
    ##################################################
    ############### color conservation ###############
    ##################################################
        
    def color_conservation(self):
        """ KC - loop on models to do color conservation """
        
        for model in self.models:
            model.color_conservation()
    
    ##################################################
    ################ analyse contacts ################
    ##################################################

    def analyse_contacts(self):
        """ KC - loop on models to analyse contacts """
        
        for model in self.models:
            model.analyse_contacts()
    
    ##################################################
    ################## get results ###################
    ##################################################
            
    def get_results(self):
        
        result_list=[]
        for model in self.models:
            dico=model.get_model_interaction()
            dico["DOPE score"]=model.DOPE
            dico["molpdf"]=model.molpdf
            dico["RMSD"]=model.get_RMSD_2(self.blank_model)
            dico["actions"]=self.actions
            result_list.append(dico)
    
        return result_list    
    
    ##################################################
    #################### get hbond ###################
    ##################################################
            
    def get_hbond(self):
        
        hbond_list=[]
        for model in self.models:
            hbond_list+=model.get_hbond()
    
        return hbond_list
        
    ##################################################
    ################## get hbond bis #################
    ##################################################
            
    def get_hbond_bis(self):
        
        hbond_list=[]
        for model in self.models:
            hbond_list+=model.get_hbond_bis()
    
        return hbond_list
