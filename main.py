#! -*- coding: utf-8 -*-

#KC#CG# packages import
import gc
import os
import time
import shutil
import Bio.PDB
import numpy as np
import random as rd
from multiprocessing import Pool
from modeller import *

#CG# no log output display
log.none()

#KC#CG# own files import
from fonctions import *
from modelisation import *
from analysis import *

#KC#CG# garbage collector activation
gc.enable()

#KC# working directory
workdir=os.getcwd()

MAKE=True

########################################################################
##############                                            ##############
##############        GPCR homology modeling (MAIN)       ##############
##############                                            ##############
######################################################################## --> CG

class GPCR_homology_modeling:
    """ KC - class for GPCR homology modeling (MAIN)  """
    
    #KC#CG# class "constructor"
    def __init__(self):
        
        #KC#CG# variables declaration
        self.filein=open("parametres.txt", 'r')
        self.TMs={}
        self.actions=[]
        self.n_proc=1
        self.n_models=10
        self.N_models=1
        self.pdb_template=None
        self.sequence_name=None
        self.sequence_filename=None
        self.alignment=True
        self.pdb_query_structure=None

        #KC#CG# functions launching
        self.parsing_parameter_file()
        if os.path.isdir('blank_model'):
            self.archive(_type="error")
            if self.GPCR_structures_test==True:
                self.archive(_type="error", _test=True)
        if self.GPCR_structures_test:
            self.GPCR_test()
            self.archive(_test=True)            
        else:
            if self.alignment:
                self.seq_alignement()
            self.pool=Pool(int(self.n_proc-1)) #CG# multiprocessing reservation
            self.seq_template,self.seq_model=get_sequences('res_align.ali', self.pdb_template, self.sequence_name) #KC#CG# sequences collected
            self.blank_model_modelisation()
            if self.pdb_query_structure != None : 
                self.blank_model.structure_alignement(self.query_structure_model)
            self.set_actions()
            try :
                self.modeling.run()
            except:
                pass

            if self.opt_boolean:
                self.opt="H_bonds"
                self.optimization()
                #self.selection()

        clean_tempfile()
        if not self.GPCR_structures_test:
            self.archive()
        
    ##################################################
    ################### archive ######################
    ################################################## --> CG
    
    def archive(self, _type=None, _test=False):
        """ KC - archive creation """
        
        #KC#CG# folder creation with unique name
        if _type=="error":
            dirname="error_archive_"+time.strftime('%d-%m-%y_%HH%Mmin%S',time.localtime())
        else :
            dirname="archive_"+time.strftime('%d-%m-%y_%HH%Mmin%S',time.localtime())
        
        #KC#CG# files moved
        if _test==False:
            if _type!="error":
                os.mkdir("template_model")
                os.system("mv "+self.template_model.filename[:-4]+".hbnds "+self.template_model.filename[:-4]+".psa "+self.template_model.filename[:-4]+"*.png template_model/. > /dev/null 2>&1")
                if self.query_structure_model != None:                
                    os.mkdir("query_model")
                    os.system("mv "+self.query_structure_model.filename[:-4]+".hbnds "+self.query_structure_model.filename[:-4]+".psa "+self.query_structure_model.filename[:-4]+"*.png query_model/. > /dev/null 2>&1")
            os.mkdir(dirname)
            os.system("mv blank_model op* sel* query* template* *.csv *.hbnds *conservation* *aligned* *angles* *network* *.psa *asa* *core* *access* *matrix* probas* with_* "+dirname+"/. > /dev/null 2>&1")
        elif _test==True :
            os.chdir(self.GPCR_structures_path)
            os.system("rm *.sol "+self.GPCR_structures_path+"/"+dirname+"/. > /dev/null 2>&1")
            os.mkdir(dirname)
            os.system("mv *results* *.png *.txt "+self.GPCR_structures_path+"/"+dirname+"/. > /dev/null 2>&1")
            os.chdir(workdir)
            
    ##################################################
    ########### parsing parameter file ###############
    ################################################## --> CG
    
    def parsing_parameter_file(self):
        """ KC - parameters file loaded """
        
        #KC#CG# file browse line by line
        for line in self.filein:
            line=line.strip().replace('\t',' ')
            if line !="" and not line.startswith("#"):
                tab=line.split(" ")
                tab=[ele for ele in tab if ele!=""]
                
                #KC#CG# parameters added to lists
                if tab[0]=="action":
                    self.actions.append(tab[1:])
                elif tab[0].startswith("TM"):
                    self.TMs[tab[0]]=tab[1:]
                elif tab[0]=="sequence_name":
                    self.sequence_name=tab[1]
                elif tab[0]=="pdb_template":
                    self.pdb_template=tab[1]
                    if self.pdb_template.endswith(".pdb"):
                        self.pdb_template=self.pdb_template[:-4]
                elif tab[0]=="pdb_query_structure":
                    self.pdb_query_structure=tab[1]
                    if self.pdb_query_structure.endswith(".pdb"):
                        self.pdb_query_structure=self.pdb_query_structure[:-4]
                elif tab[0]=="n_models":
                    self.n_models=int(tab[1])
                elif tab[0]=="n_proc":
                    self.n_proc=int(tab[1]) 
                elif tab[0]=="sequence_filename":
                    self.sequence_filename=tab[1]
                elif tab[0]=="alignment":
                    if tab[1].lower()=="true":
                        self.alignment=True
                    elif tab[1].lower()=="false":        
                        self.alignment=False
                    else:
                        print "unsuported value for alignement in parameter file"
                elif tab[0]=="optimization":
                    if tab[1].lower()=="true":
                        self.opt_boolean=True
                    elif tab[1].lower()=="false":        
                        self.opt_boolean=False
                    else:
                        print "unsuported value for alignement in parameter file"
                elif tab[0]=="p_min":
                    self.p_min=int(tab[1])
                elif tab[0]=="p_max":
                    self.p_max=int(tab[1])
                elif tab[0]=="p_step":
                    self.p_step=int(tab[1])
                elif tab[0]=="ASA_per":
                    self.ASA_per=float(tab[1])
                elif tab[0]=="selection_criteria":
                    self.selection_criteria=tab[1]
                elif tab[0]=="GPCR_structures_test":
                    if tab[1].lower()=="true":
                        self.GPCR_structures_test=True
                    elif tab[1].lower()=="false":        
                        self.GPCR_structures_test=False
                    else:
                        print "unsuported value for alignement in parameter file"
                elif tab[0]=="N_models":
                    self.N_models=int(tab[1])
                elif tab[0]=="GPCR_structures_path":
                    self.GPCR_structures_path=tab[1]
                else:
                    print "unsuported keyword in parameter file"
                    print line

        #parameters tests
        assert isinstance(self.sequence_name, str), "ERROR: sequence_name have to be a string"
        assert isinstance(self.sequence_filename, str), "ERROR: sequence_filename have to be a string"
        assert isinstance(self.pdb_template, str), "ERROR: pdb_template have to be a string"
        assert (isinstance(self.n_models, int) or self.n_models>1), "ERROR: n_models have to be an integer > 1"
        assert os.path.isfile(self.sequence_filename), 'ERROR: sequence file '+self.sequence_filename+' not exist'
        assert os.path.isfile(self.pdb_template+'.pdb'), 'ERROR: sequence file '+self.pdb_template+' not exist'
        
        self.filein.close()


    ##################################################
    ################### GPCR test ####################
    ################################################## --> KC

    def GPCR_test(self):
        """ KC - analyse of each GPCR structure """
        
        os.system('ls '+self.GPCR_structures_path+'> GPCR_structures_name.txt')
        self.filein_GPCR=open("GPCR_structures_name.txt", 'r')
        self.GPCR_structures_name=[]
        for line in self.filein_GPCR:
            line=line.strip()
            if line.endswith(".pdb"):
                self.GPCR_structures_name.append(line)
        self.filein_GPCR.close()
        os.system('rm GPCR_structures_name.txt > /dev/null 2>&1')
        for structure_name in self.GPCR_structures_name:
            structure_name_path=self.GPCR_structures_path+"/"+structure_name
            my_model(structure_name_path, sequence_name=structure_name[:-4], _test=True, ASA_per=self.ASA_per)
            os.chdir(self.GPCR_structures_path)
            folder_name=structure_name[:-4]+"_results"
            os.mkdir(folder_name)
            os.system("mv *.hbnds *.psa *Hbond* *stacking* *access* "+folder_name+"/. > /dev/null 2>&1")
            os.chdir(workdir)
        os.system("./R_scripts/GPCR_one_graph.R " + workdir)
        os.system("mv *_test* *_all* "+self.GPCR_structures_path+"/. > /dev/null 2>&1")

    ##################################################
    ########### blank model modelisation #############
    ################################################## --> CG

    def blank_model_modelisation(self):
        """ KC - modelisation of a blank model """
        
        #KC#CG# MODELLER automodel realisation
        env=environ()
        env.io.hetatm=True #CG# read in HETATM records from template PDBs
        a=automodel(env, alnfile='res_align.ali', knowns=self.pdb_template, sequence=self.sequence_name, assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model=1
        a.ending_model=1
        a.make()
        
        #KC#CG# scores collected
        self.molpdf_blank=a.outputs[0]["molpdf"]
        self.DOPE_blank=a.outputs[0]["DOPE score"]
        self.GA341_blank=a.outputs[0]["GA341 score"]
            
        #KC#CG# folder creation, file moved into and file renamed
        os.mkdir("blank_model")
        shutil.copyfile(self.sequence_name+".B99990001.pdb", "blank_model/blank.pdb")
        shutil.move(self.sequence_name+".B99990001.pdb", "blank.pdb")

        #KC# filename of blank pdb definition
        self.pdb_blank="blank_model/blank"
        self.pdb_filename_blank = "%s.pdb" % self.pdb_blank

        #KC#CG# models (blank and template) realisation and pdb file colored by conservation
        self.blank_model=my_model(self.pdb_filename_blank, seq_model=self.seq_model, seq_template=self.seq_template, sequence_name=self.sequence_name, DOPE=self.DOPE_blank, molpdf=self.molpdf_blank, GA341=self.GA341_blank, _type="blank", ASA_per=self.ASA_per)
        self.template_model=my_model(self.pdb_template+".pdb", seq_model=self.seq_template, seq_template=self.seq_model, sequence_name=self.pdb_template, _type="xray", ASA_per=self.ASA_per)
        if self.pdb_query_structure != None :        
            self.query_structure_model=my_model(self.pdb_query_structure+".pdb", sequence_name=self.pdb_query_structure, _type="xray", ASA_per=self.ASA_per)
        else :
            self.query_structure_model=None

        #KC#CG# alignment between blank model structure and orexin sequence for next modelisation
        aln=alignment(env)
        mdl=model(env, file='blank')
        aln.append_model(mdl, align_codes='blank_model', atom_files='blank.pdb')
        aln.append(file=self.sequence_filename, align_codes=self.sequence_name)
        aln.align2d()
        aln.write(file='res_final.ali', alignment_format='PIR')

    ############################################
    ############## seq alignement ##############
    ############################################ --> CG
                
    def seq_alignement(self):
        """ KC - alignment 2D between model/query and result file written (res_align.ali) """
        
        env=environ()
        env.io.hetatm = True #CG# read in HETATM records from template PDBs
        aln=alignment(env)
        mdl=model(env, file=self.pdb_template) #KC#CG# model definition
        aln.append_model(mdl, align_codes=self.pdb_template, atom_files=self.pdb_template+'.pdb')
        aln.append(file=self.sequence_filename, align_codes=self.sequence_name)
        aln.align2d()
        aln.write(file='res_align.ali', alignment_format='PIR')
        
    ##################################################
    ################# set actions ####################
    ################################################## --> CG
    
    def set_actions(self):
        """ KC - actions in parameters file realised """

        structure=Bio.PDB.PDBParser().get_structure("blank", "blank.pdb") #KC#GC# PDB structure of blank model collected
        model=structure[0]
        actions_temp=[] #CG# action_temp to remind what's happen in each modelisation
        restraints=[]

        #KC#CG# loop on the actions of the parameters file or of the optimization
        for action in self.actions:
            actions_temp.append(action)

            #KC#CG# rotation ou translation
            if action[0] in ["rotation", "translation"]:                
                #CG# search the first point
                if action[1].startswith("TM"):
                    #CG# search begin and end of TM
                    debut=int(self.TMs[action[1][:3]][0])-int(self.TMs[action[1][:3]][1])
                    fin = int(self.TMs[action[1][:3]][0])+int(self.TMs[action[1][:3]][2])
                elif action[1].startswith("res"):
                    #GG# search begin and end of TM
                    debut=int(action[1].split(":")[1])
                    fin = int(action[1].split(":")[2])
                elif action[1]=="protein":
                    debut=-1
                    fin="Inf"
                elif action[1]=="ligand":
                    print "To be implemented..."
                
                #CG# generate some data that lies along a line
                selection=[atom for atom in model.get_atoms() if debut<=int(atom.get_full_id()[3][1])<=fin]
                data=np.array([atom.get_coord() for atom in selection])
                #CG# calculate the mean of the points, i.e. the 'center' of the cloud
                axe_point_1=data.mean(axis=0) #KC# axis=0 to take the mean of each col
                
                #CG# define axis
                if len(action[2].split(','))==3:
                    axe=[float(var) for var in action[2].split(',')]
                elif action[2] in ["x","y","z"]:
                    dico={"x":[1,0,0], "y":[0,1,0], "z":[0,0,1]}
                    axe=dico[action[2]]
                #CG# search the second point
                else:
                    if action[2].startswith("TM"):
                        #CG# search begin and end of TM
                        debut=int(self.TMs[action[2][:3]][0])-int(self.TMs[action[2][:3]][1])
                        fin = int(self.TMs[action[2][:3]][0])+int(self.TMs[action[2][:3]][2])
                    elif action[2].startswith("res"):
                        #CG# search begin and end of TM
                        debut=int(action[2].split(":")[1])
                        fin = int(action[2].split(":")[2])
                    elif action[2]=="protein":
                        debut=-1
                        fin="Inf"
                    elif action[2]=="ligand":
                        print "To be implemented..."
                    
                    #CG# generate some data that lies along a line
                    data=np.array([atom.get_coord() for atom in model.get_atoms() if debut<=int(atom.get_full_id()[3][1])<=fin])
                    #CG# calculate the mean of the points, i.e. the 'center' of the cloud
                    axe_point_2=data.mean(axis=0) #KC#CG# axis=0 to take the mean of each col
                    
                    #KC#CG# define axis
                    if np.all(axe_point_2==axe_point_1):
                        axe=determine_axe(selection, None)
                    elif action[2].endswith("_perp"):
                        axe=determine_axe(selection, axe_point_2, _type="perp")
                    else:
                        axe=determine_axe(selection, axe_point_2, _type="normal")
                
                #CG# movement
                if action[0] == "rotation":
                    pivot_tab=action[4].split(":")
                    pivot=model[pivot_tab[2]][int(pivot_tab[1])][pivot_tab[0]].get_coord() #KC#CG# coordinates of given atom collected
                    #KC#CG# rotation for atoms selected
                    for atom in selection:
                        atom.set_coord(rotation(vecteur_coord=atom.get_coord(), axe=np.array(axe), angle=2*np.pi/360*float(action[3]), centre=pivot))
                    io=Bio.PDB.PDBIO()
                    io.set_structure(structure)
                    io.save("blank_model/blank_rotation.pdb")
                
                elif action[0] == "translation":
                    #KC#CG# translation for atoms selected
                    for atom in selection:
                        atom.set_coord(translation(vecteur_coord=atom.get_coord(),valeur=float(action[3]), axe=axe))
                    io=Bio.PDB.PDBIO()
                    io.set_structure(structure)
                    io.save("blank_model/blank_translation.pdb")
            
            #KC#CG# modelisations
            if action[0]=="modelisation":
                #KC#CG# folder creation if possible
                try:
                    os.mkdir(action[1])
                except:
                    pass
                #KC#CG# blank structure saved
                io=Bio.PDB.PDBIO()
                io.set_structure(structure)
                io.save(action[1]+"/blank.pdb")
                #KC#CG# modelisations added to the list
                self.modeling=modelisation(job_name=action[1], sequence_name=self.sequence_name, pdb_template=self.pdb_template, seq_model=self.seq_model, seq_template=self.seq_template, \
                n_models=self.n_models, sequence_filename=self.sequence_filename, restraints=restraints, actions=actions_temp, blank_model=self.blank_model)
            
            #KC#CG# restraints added to the list
            elif action[0]=="restraint":
                restraints.append(action[1:])
                
            #KC#CG# restraints removed
            elif action[0]=="clear_restraint":
                restraints=[]
                
            #KC#CG# template initialization
            elif action[0]=="init_template":
                structure=Bio.PDB.PDBParser().get_structure("blank", "blank.pdb")
                model=structure[0]
                restraints=[]
                actions_temp=[]

    ##################################################
    ################# optimization ###################
    ################################################## --> KC

    def optimization(self):
        """ KC - optimization of number of core hydrogen bonds """

        #KC# list of donors/acceptors with parametres.py file
        if self.opt=="H_bonds":
            self.list_donor=[(donor.get_full_id()[4][0])+":"+str(donor.get_full_id()[3][1])+":"+str(donor.get_full_id()[2]) for donor in self.blank_model.get_donors()]
            self.list_acceptor=[(acceptor.get_full_id()[4][0]+":"+str(acceptor.get_full_id()[3][1])+":"+str(acceptor.get_full_id()[2])) for acceptor in self.blank_model.get_acceptors()]
        elif self.opt=="ionic_bonds":
            self.list_donor=[(donor.get_full_id()[4][0])+":"+str(donor.get_full_id()[3][1])+":"+str(donor.get_full_id()[2]) for donor in self.blank_model.get_ionic_donors() if donor.get_full_id()[4][0]!="N"]
            self.list_acceptor=[(acceptor.get_full_id()[4][0]+":"+str(acceptor.get_full_id()[3][1])+":"+str(acceptor.get_full_id()[2])) for acceptor in self.blank_model.get_ionic_acceptors() if acceptor.get_full_id()[4][0]!="O"]
        
        #KC# matrix initialization
        self.contact_matrix_probas=[[self.p_min for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]
        self.contact_matrix_blank=[[0 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]
        self.contact_matrix_query=[[0 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]
        self.contact_matrix_best_model=[[0 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]

        #KC# negative weight applied in contact matrix
        self.zones_constraints()

        for (atom1,atom2) in self.blank_model.interactions.get_bonds(self.opt):
            self.contact_matrix_blank[self.list_donor.index(atom1)][self.list_acceptor.index(atom2)] = 1
        save_matrix(self.contact_matrix_blank, "contact_matrix_blank", self.list_acceptor, self.list_donor)
        
        """self.num_query,self.num_model=get_aligned_residues('structure_sequence.ali')
        dictionary=dict(zip(self.num_query, self.num_model))
        for (atom1,atom2) in self.query_structure_model.interactions.get_bonds(self.opt):
            atom1=atom1.split(":")
            atom2=atom2.split(":")
            chain1=[res.get_full_id()[2] for res in self.blank_model.res if res.get_full_id()[3][1]==dictionary[int(atom1[1])]]
            chain2=[res.get_full_id()[2] for res in self.blank_model.res if res.get_full_id()[3][1]==dictionary[int(atom2[1])]]
            atom1=atom1[0]+":"+str(dictionary[int(atom1[1])])+":"+chain1[0]
            atom2=atom2[0]+":"+str(dictionary[int(atom2[1])])+":"+str(chain2[0])
            self.contact_matrix_query[self.list_donor.index(atom1)][self.list_acceptor.index(atom2)] = 1"""

        save_matrix(self.contact_matrix_query, "contact_matrix_query", self.list_acceptor, self.list_donor)

        self.all_models=[]
        self.list_tanimoto=[]
        self.sum_probas=[]
        self.sum_Hbonds=[]
        self.sum_restraints=[]
        self.blank_Hbonds_number=len(self.blank_model.interactions.get_bonds(self.opt))
        for step in xrange(500):   
            self.actions=[]
            for donor in rd.sample(self.list_donor, len(self.list_donor)/3):
                acceptors=[]
                for acceptor in self.list_acceptor:
                    if self.contact_matrix_probas[self.list_donor.index(donor)][self.list_acceptor.index(acceptor)]>0:
                            acceptors.append(acceptor)
                for acceptor in rd.sample(acceptors, len(acceptors)/2):
                    random_number=int(rd.random()*100)
                    if random_number>= 0 and random_number<self.contact_matrix_probas[self.list_donor.index(donor)][self.list_acceptor.index(acceptor)]:
                        self.actions.append(['restraint', "upperdistance", donor, acceptor, 6, 3])
            self.actions.append(["modelisation", "optimization"+str(step+1)])
            self.set_actions()
            self.modelisation_run("optimization")
            if (step+1) in [1, 50, 100, 200, 300, 500]:
                for model in self.modeling.models:
                    model.interactions.Hbonds_dist_model(self.list_donor, step+1)
                if self.query_structure_model != None:
                    for model in self.modeling.models:
                        analysis(None, None, None, None, self.query_structure_model, self.list_donor, self.list_acceptor, None, None, None, self.opt, \
                            None, None, None, None, None, None, None, True, model, step+1)     
        save_matrix(self.contact_matrix_probas, "contact_matrix_probas", self.list_acceptor, self.list_donor)
        #print self.contact_matrix_probas

        #KC# models analysis
        analysis(self.all_models, self.n_models, self.blank_model, self.template_model, self.query_structure_model, self.list_donor, self.list_acceptor, \
            self.contact_matrix_probas, self.contact_matrix_best_model, self.contact_matrix_query, self.opt, self.p_min, self.p_max, self.p_step, self.list_tanimoto, \
            self.sum_probas, self.sum_Hbonds, self.sum_restraints, False, None, None)

        os.mkdir("optimization")
        os.system("mv *matrix* opt* models* probas_dist* Hbonds_dist* optimization/. > /dev/null 2>&1")      
        
    ##################################################
    #################### selection ###################
    ################################################## --> KC

    def selection(self):
        """ KC - selection according to selection criteria """       

        self.best_model=self.blank_model
        
        #KC# loop with arbitrary number
        for step in xrange(self.N_models):   
            self.actions=[]
            #CG# random choice of restraints which will be done
            for donor in self.list_donor:
                acceptors=[acceptor for acceptor in self.list_acceptor if self.contact_matrix_probas[self.list_donor.index(donor)][self.list_acceptor.index(acceptor)] > 0]
                for acceptor in rd.sample(acceptors, len(acceptors)/2):
                    random_number=int(rd.random()*100)
                    if random_number>= 0 and random_number<self.contact_matrix_probas[self.list_donor.index(donor)][self.list_acceptor.index(acceptor)]:
                        self.actions.append(['restraint', "upperdistance", donor, acceptor, 6, 3])
            self.actions.append(["modelisation", "selection"+str(step+1)])
            self.set_actions()
            self.modelisation_run("selection")
            for model in self.modeling.models:
                model.interactions.Hbonds_dist_model(self.list_donor, step+1)
            
        if self.best_model != self.blank_model:
            print "\n Best model :"
            print "filename : "+str(self.best_model.filename)
            print "GA341 score : "+str(self.best_model.GA341[0])
            print "molpdf score : "+str(self.best_model.molpdf)
            print "DOPE score : "+str(self.best_model.DOPE)
            print self.selection_criteria+" : "+str(self.best_selection_criteria)
            for restraint in self.best_modelisation.restraints:
                if [restraint[1],restraint[2]] in self.best_model.interactions.get_bonds(self.opt):
                    self.contact_matrix_best_model[self.list_donor.index(restraint[1])][self.list_acceptor.index(restraint[2])]=1
            save_matrix(self.contact_matrix_best_model, "contact_matrix_best_model", self.list_acceptor, self.list_donor)
            analysis(None, None, None, None, self.query_structure_model, self.list_donor, self.list_acceptor, None, None, None, self.opt, \
                None, None, None, None, None, None, None, True, self.best_model, "final")
            self.best_model.structure_alignement(self.query_structure_model)
            self.best_model.interactions.network_distribution("Hbond")
            self.best_model.interactions.accessibility_VS_Hbond()
        else :
            print "No better model than blank found !"
            print "filename : "+str(self.blank_model.filename)
            print "GA341 score : "+str(self.blank_model.GA341[0])
            print "molpdf score : "+str(self.blank_model.molpdf)
            print "DOPE score : "+str(self.blank_model.DOPE)
            print self.selection_criteria+" : "+str(-(self.best_selection_criteria))
            analysis(None, None, None, None, self.query_structure_model, self.list_donor, self.list_acceptor, None, None, None, self.opt, \
                None, None, None, None, None, None, None, True, self.blank_model, "final")
            
        os.mkdir("selection")
        os.system("mv *matrix* sel* Hbonds_dist* selection/. > /dev/null 2>&1") 
        
    ##################################################
    ################ zones constraints ###############
    ################################################## --> KC

    def zones_constraints(self):
        """ KC - constraints applied according to protein zones """
        
        #KC# list of donors/acceptors to use "get_parent()"
        if self.opt=="H_bonds":
            LD=[donor for donor in self.blank_model.get_donors()]
            LA=[acceptor for acceptor in self.blank_model.get_acceptors()]
        elif self.opt=="ionic_bonds":
            LD=[donor for donor in self.blank_model.get_ionic_donors() if donor.get_full_id()[4][0]!="N"]
            LA=[acceptor for acceptor in self.blank_model.get_ionic_acceptors() if acceptor.get_full_id()[4][0]!="O"]
            
        #KC# negative weight for same donor/acceptor
        for donor in LD:
            for acceptor in LA:
                res_acceptor = acceptor.get_parent()
                res_donor = donor.get_parent()
                #KC# negative weight for same donor/acceptor
                if acceptor == donor :
                    self.contact_matrix_probas[LD.index(donor)][LA.index(acceptor)] = -10000
                elif (res_acceptor in self.blank_model.extra and res_donor in self.blank_model.intra) or (res_acceptor in self.blank_model.intra and res_donor in self.blank_model.extra) :
                    self.contact_matrix_probas[LD.index(donor)][LA.index(acceptor)] = -10000

    ##################################################
    ################ modelisation run ################
    ################################################## --> KC

    def modelisation_run(self, part):
        """ KC - run of modelisations and probability in matrix changed """
        
        if part=="optimization":
            
            self.all_models.append(self.modeling.run())
            for model in self.modeling.models:
                list_tanimoto_temp=[]
                model_Hbond=model.interactions.get_bonds(self.opt)
                self.sum_Hbonds.append(len(model_Hbond))
                if len(model_Hbond) > self.blank_Hbonds_number:
                    for restraint in self.modeling.restraints:
                        if [restraint[1],restraint[2]] in model_Hbond:
                            if restraint not in list_tanimoto_temp:
                                list_tanimoto_temp.append(restraint)
                            if self.contact_matrix_probas[self.list_donor.index(restraint[1])][self.list_acceptor.index(restraint[2])]<self.p_max:
                                self.contact_matrix_probas[self.list_donor.index(restraint[1])][self.list_acceptor.index(restraint[2])] += self.p_step
            self.sum_restraints.append(len(list_tanimoto_temp))
            self.list_tanimoto.append(list_tanimoto_temp)
            sum_probability=0
            for donor in self.list_donor:
                for acceptor in self.list_acceptor:
                    if self.contact_matrix_probas[self.list_donor.index(donor)][self.list_acceptor.index(acceptor)] > 0:
                        sum_probability+=self.contact_matrix_probas[self.list_donor.index(donor)][self.list_acceptor.index(acceptor)]
            self.sum_probas.append(sum_probability)
        
        elif part =="selection":
            
            self.modeling.run()
            for model in self.modeling.models:
                if self.selection_criteria == "RMSD_VS_query_structure":
                    self.best_selection_criteria=float(-(self.best_model.get_RMSD_query_structure(self.query_structure_model)))
                    selection_criteria=float(-(model.get_RMSD_query_structure(self.query_structure_model)))
                elif self.selection_criteria == "RMSD_VS_template":
                    self.best_selection_criteria=float(-(self.best_model.get_RMSD_template(self.template_model)))
                    selection_criteria=float(-(model.get_RMSD_template(self.template_model)))
                elif self.selection_criteria == "number_Hbonds_core":
                    self.best_selection_criteria=len(self.best_model.interactions.search_H_bond_core(self.best_model.res_core))
                    selection_criteria=float(model.interactions.search_H_bond_core(model.res_core))
                elif self.selection_criteria == "length_network_core":
                    self.best_selection_criteria=len(max(network for network in self.best_model.interactions.H_bonds_networks))
                    selection_criteria=len(max(network for network in model.interactions.H_bonds_networks))
                elif self.selection_criteria == "number_residues_core":
                    self.best_selection_criteria=len(self.best_model.res_core)
                    selection_criteria=len(model.res_core)
                if (selection_criteria>self.best_selection_criteria) :
                        print "YES"
                        self.best_model=model
                        self.best_modelisation=self.modeling
                        if self.selection_criteria in ["RMSD_VS_query_structure", "RMSD_VS_template"]:
                            self.best_selection_criteria=-(self.best_selection_criteria)

        
if __name__ == '__main__':
    GPCR_homology_modeling()
 
