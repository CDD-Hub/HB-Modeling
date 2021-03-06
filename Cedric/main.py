#! -*- coding: utf-8 -*-

#KC#CG# packages import
import gc
import os
import time
import shutil
import random as rd
from multiprocessing import Pool
import Bio.PDB
from Bio.PDB.Vector import *
from pylab import *
from modeller import *

#CG# no log output display
log.none()

#KC#CG# own files import
from fonctions import *
from modelisation import *

#KC#CG# garbage collector activation
gc.enable()

MAKE=True

########################################################################
########################                        ########################
########################        functions       ########################
########################                        ########################
########################################################################

def run_multiprocessing(tab):
    """ KC - run of modelisation in the list """
    
    modelisation,start_analysis,end_analysis=tab
    modelisation.run()
    for model in modelisation.models:
        model.interactions.run(start=start_analysis, end=end_analysis)
    return modelisation

########################################################################
##############                                            ##############
##############        GPCR homology modeling (MAIN)       ##############
##############                                            ##############
########################################################################

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
        self.pdb_ref=None
        self.pdb_template=None
        self.sequence_name=None
        self.sequence_filename=None
        self.list_modelisations=None
        self.alignment=True
        
        #KC#CG# functions launching
        self.archive()
        self.parsing_parameter_file()
        if self.alignment:
            self.seq_alignement()
        self.pool=Pool(int(self.n_proc-1)) #CG# multiprocessing reservation
        self.seq_template,self.seq_model=get_sequences('res_align.ali', self.pdb_template, self.sequence_name) #KC#CG# sequences collected
        self.conservation=[AA1==AA2 for AA1,AA2 in zip(self.seq_template,self.seq_model)] #KC#CG# list of boolean (same amino acid or not)
        self.blank_model_modelisation()
        #KC#CG# self.set_actions()
        #KC#CG# self.run()
        #CG# self.hbond_optimisation()
        #CG# self.helix_optimisation()
        #CG# self.helix_hbond_optimisation()       
        
    ##################################################
    ################### archive ######################
    ##################################################
    
    def archive(self):
        """ KC - archive creation """
        
        #KC#CG# folder creation with unique name
        dirname="archive_"+time.strftime('%d-%m-%y_%HH%Mmin%S',time.localtime())
        os.mkdir(dirname)
        #KC#CG# files moved
        os.system("mv blank_model ?_? ??_? ???_? ????_? ?????_? "+dirname+"/. > /dev/null 2>&1")
        
    ##################################################
    ########### parsing parameter file ###############
    ##################################################
    
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
                elif tab[0]=="n_models":
                    self.n_models=int(tab[1])
                elif tab[0]=="n_proc":
                    self.n_proc=int(tab[1])
                elif tab[0]=="pdb_ref":
                    if tab[1].lower()=="none":
                        self.pdb_ref=None
                    else:
                        self.pdb_ref=tab[1]
                elif tab[0]=="sequence_filename":
                    self.sequence_filename=tab[1]
                elif tab[0]=="alignment":
                    if tab[1].lower()=="true":
                        self.alignment=True
                    elif tab[1].lower()=="false":        
                        self.alignment=False
                    else:
                        print "unsuported value for alignement in parameter file"
                else:
                    print "unsuported keyword in parameter file"
                    print line

        #parameters tests
        assert isinstance(self.sequence_name, str), "ERROR: sequence_name have to be a string"
        assert isinstance(self.sequence_filename, str), "ERROR: sequence_filename have to be a string"
        assert isinstance(self.pdb_template, str), "ERROR: pdb_template have to be a string"
        assert (isinstance(self.pdb_ref, str) or self.pdb_ref is None), "ERROR: pdb_ref have to be a string"
        assert (isinstance(self.n_models, int) or self.n_models>1), "ERROR: n_models have to be an integer > 1"
        assert os.path.isfile(self.sequence_filename), 'ERROR: sequence file '+self.sequence_filename+' not exist'
        assert os.path.isfile(self.pdb_template+'.pdb'), 'ERROR: sequence file '+self.pdb_template+' not exist'
        if self.pdb_ref is not None:
            assert os.path.isfile(self.pdb_ref+'.pdb'), 'ERROR: sequence file '+self.pdb_ref+' not exist'
        
        self.filein.close()

    ##################################################
    ########### blank model modelisation #############
    ##################################################

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
        self.DOPE_blank=a.outputs[0]["DOPE score"]
        self.molpdf_blank=a.outputs[0]["molpdf"]
        
        #KC# last folder removed
        if os.path.isdir('blank_model'):
            shutil.rmtree("blank_model")
            
        #KC#CG# folder creation, file moved into and file renamed
        os.mkdir("blank_model")
        shutil.copyfile(self.sequence_name+".B99990001.pdb", "blank_model/blank.pdb")
        shutil.move(self.sequence_name+".B99990001.pdb", "blank.pdb")

        #KC#CG# filename of blank pdb definition
        self.pdb_blank="blank_model/blank"
        self.pdb_filename_blank = "%s.pdb" % self.pdb_blank

        #KC#CG# models (blank and template) realisation and pdb file colored by conservation
        self.blank_model=my_model(self.pdb_filename_blank, self.seq_model, self.seq_template, self.sequence_name, None)
        self.blank_model.color_conservation()
        self.blank_model.ASA_calcul(sequence_name="blank")
        self.template_model=my_model(self.pdb_template+".pdb", self.seq_template, self.seq_model, self.pdb_template, None, _type="xray")
        self.template_model.color_conservation()
        self.template_model.ASA_calcul()
        
        #KC#CG# files with .B99, .png, colored or conservation moved
        move_files("blank_model")
        clean_tempfile()
        
        #KC#CG# alignment between blank model structure and orexin sequence for next modelisation
        env.io.hetatm = True #CG# read in HETATM records from template PDBs
        aln=alignment(env)
        mdl=model(env, file='blank')
        aln.append_model(mdl, align_codes='blank_model', atom_files='blank.pdb')
        aln.append(file=self.sequence_filename, align_codes=self.sequence_name)
        aln.align2d()
        aln.write(file='res_final.ali', alignment_format='PIR')

    ############################################
    ############## seq alignement ##############
    ############################################
                
    def seq_alignement(self):
        """ KC - alignment 2D between model/query and result file written (res_align.ali) """
        
        env=environ()
        env.io.hetatm = True #read in HETATM records from template PDBs
        aln=alignment(env)
        mdl=model(env, file=self.pdb_template) #KC#CG# model definition
        aln.append_model(mdl, align_codes=self.pdb_template, atom_files=self.pdb_template+'.pdb')
        aln.append(file=self.sequence_filename, align_codes=self.sequence_name)
        aln.align2d()
        aln.write(file='res_align.ali', alignment_format='PIR')
        #CG# aln.to_profile().write(file='alntoprof.prf', profile_format='TEXT')
        
    ##################################################
    ################# set actions ####################
    ##################################################
    
    def set_actions(self):
        """ KC - actions in parameters file realised """
        
        #CG# pdb_code = self.pdb_template
        #CG# pdb_filename = "%s.pdb" % pdb_code

        structure=Bio.PDB.PDBParser().get_structure("blank", "blank.pdb") #KC# PDB structure of blank model collected
        model=structure[0]

        self.list_modelisations=[] #CG# list_modelisations to run and analyses
        actions_temp=[] #CG# action_temp to remind what's happen in each modelisation
        restraints=[]

        #KC#CG# loop on the actions of the parameters file
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
                        atom.set_coord(rotation(vecteur_coord=atom.get_coord(), axe=np.array(axe), angle=2*pi/360*float(action[3]), centre=pivot))
                    io=Bio.PDB.PDBIO()
                    io.set_structure(structure)
                    io.save(self.sequence_name+"_rotation.pdb")
                
                elif action[0] == "translation":
                    #KC#CG# translation for atoms selected
                    for atom in selection:
                        atom.set_coord(translation(vecteur_coord=atom.get_coord(),valeur=float(action[3]), axe=axe))
                    io=Bio.PDB.PDBIO()
                    io.set_structure(structure)
                    io.save(self.sequence_name+"_translation.pdb")
            
            #KC#CG# modelisations
            elif action[0]=="modelisation":
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
                self.list_modelisations.append(modelisation(job_name=action[1],
                                                            sequence_name=self.sequence_name,
                                                            pdb_template=self.pdb_template,
                                                            seq_model=self.seq_model,
                                                            seq_template=self.seq_template,
                                                            n_models=self.n_models,
                                                            sequence_filename=self.sequence_filename,
                                                            restraints=restraints,
                                                            actions=actions_temp, 
                                                            blank_model=self.blank_model))
            
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
    ###################### run #######################
    ##################################################

    def run(self, start_analysis=0, end_analysis=-1):
        """ KC - run of modelisations in the list together """
        
        """ CG - for modelisation in self.list_modelisations:
                modelisation.run()
                for model in modelisation.models:
                model.interactions.run(start=start_analysis, end=end_analysis)"""

        self.list_modelisations=self.pool.map(run_multiprocessing, [[modelisation, start_analysis, end_analysis] for modelisation in self.list_modelisations] )
        return    
        
    ##################################################
    ############## helix optimisation ################
    ##################################################

    def helix_optimisation(self):
        """ KC - random of helix restraint (action restraint helix 11 21) """
        
        """ CG -
        The helix optimisation follow this steps:
        - optimisation in 3D of the helix
        - optimisation in 3D of the half-helix around the pivot
        """
        
        for step in xrange(2):
            self.actions=[]
            self.actions.append(["modelisation","ref"])
            
            for i in xrange(self.n_proc-2): #CG# for each model
                for change in xrange(5): #CG# apply 5 change in model
                    self.actions.append(["init_template"])
                    
                    #KC#CG# random choice of helix
                    helix=rd.choice([TM for TM in self.TMs])
                    part=rd.choice(["all","N","C"])
                    direction=rd.choice(["x", "y", "z", helix])
                    print helix, part, direction
                    
                    #KC#CG# helix parameters
                    pivot = int(self.TMs[helix][0])
                    start = pivot-int(self.TMs[helix][1])
                    end = pivot+int(self.TMs[helix][2])
                    
                    #KC#CG# translation or rotation according to the part of helix
                    if part=="all":
                        mvt=rd.choice(["translation", "rotation"])
                    else:
                        mvt="rotation"
                        if part=="C":
                            helix="res:"+str(start)+":"+str(pivot)
                        else:
                            helix="res:"+str(pivot)+":"+str(end)
                            
                    #KC#CG# actions added
                    if mvt=="translation":
                        d=rd.uniform(-2,2)
                        self.actions.append(['translation', helix, direction, str(d)])
                    elif mvt=="rotation":
                        if direction == helix and part=="all":
                            a=rd.uniform(-20,20)
                        else:
                            a=rd.uniform(-10,10)
                            self.actions.append(['rotation', helix, direction, str(a), "CA:"+str(pivot)+": "])
                
                self.actions.append(["modelisation",str(step)+"_"+str(i+1)])
            self.set_actions()
            self.run(start_analysis=start, end_analysis=end)
            self.run()
            
            #CG# new blank moved
            name=self.get_results()
            if name=="":
                shutil.copy("blank.pdb", "blank_model/blank_"+str(step)+".pdb")
            else:
                os.rename("blank.pdb", "blank_model/blank_"+str(step)+".pdb")
                shutil.copy(name, "blank.pdb")
            
            print "#####################################################"
            print name
            print "#####################################################"
            
        shutil.copy("blank.pdb", "blank_model/blank_"+str(step)+".pdb")

    ##################################################
    ############## hbond optimisation ################
    ##################################################

    def hbond_optimisation(self):
        """ KC - random of distance restraint between donor and acceptor (action restraint distance CA:100:A CA:200:A 10 0.1) """
        
        #KC#CG# list of donors/acceptors with parametres.py file
        self.list_donor=[donor.get_full_id()[4][0]+":"+str(donor.get_full_id()[3][1])+":"+donor.get_full_id()[2] for donor in self.blank_model.get_donors() if donor.get_full_id()[4][0]!="N"]
        self.list_acceptor=[acceptor.get_full_id()[4][0]+":"+str(acceptor.get_full_id()[3][1])+":"+acceptor.get_full_id()[2] for acceptor in self.blank_model.get_acceptors() if acceptor.get_full_id()[4][0]!="O"]

        #CG# dico to easy access
        self.dico_donor={name:number for number, name in enumerate(self.list_donor)}
        self.dico_acceptor={name:number for number, name in enumerate(self.list_acceptor)}
        
        #KC#CG# matrix initialisation
        matrix_ref=[[0 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]
        
        #KC#CG# loop with arbitrary number
        for step in xrange(1):
            self.actions=[]
            for n_model in xrange(self.n_proc-1):
                self.actions.append(["init_template"])
                
                #CG# random choice of restraints which will be done
                for indice_donor, line in enumerate(matrix_ref):
                    #CG# [2,5,1]==>[0,0,1,1,1,1,1,2]
                    value_weight=[]
                    for value,n in enumerate(line):
                        outset=100
                        value_weight+=[value for i in xrange(n+outset)]
                    #CG# random choice
                    for indice_acceptor in rd.sample(value_weight,5):
                        #CG# self.actions.append(['restraint', "distance", self.list_donor[indice_donor], self.list_acceptor[indice_acceptor], 3, max(2-0.001*line[indice_acceptor], 1.2, 100)])
                        self.actions.append(['restraint', "distance", self.list_donor[indice_donor], self.list_acceptor[indice_acceptor], 3, 1.5-0.005*line[indice_acceptor]]) # ???
                        
                self.actions.append(["modelisation", str(step+1)+"_"+str(n_model+1)])
            self.set_actions()
            #KC#CG# self.run()

            matrix_ref=self.get_hbond(matrix_ref)
            
            #CG# hbond matrix saved
            with open("result_hbond.csv", "w") as fileout:
                for line in matrix_ref:
                    fileout.write("\t".join([str(i) for i in line])+"\n")
                    
    ##################################################
    ######### helix and hbond Optimisation ###########
    ##################################################

    def helix_hbond_optimisation(self):
        """ KC - helix and hbond optimisation together """
        
        self.list_donor=[donor.get_full_id()[4][0]+":"+str(donor.get_full_id()[3][1])+":"+donor.get_full_id()[2] for donor in self.blank_model.get_donors() if donor.get_full_id()[4][0]!="N"]
        self.list_acceptor=[acceptor.get_full_id()[4][0]+":"+str(acceptor.get_full_id()[3][1])+":"+acceptor.get_full_id()[2] for acceptor in self.blank_model.get_acceptors() if acceptor.get_full_id()[4][0]!="O"]

        self.dico_donor={name:number for number, name in enumerate(self.list_donor)}
        self.dico_acceptor={name:number for number, name in enumerate(self.list_acceptor)}
        
        matrix_ref=[[0 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]

        for step in xrange(1,2,1):       
            self.actions=[]
            self.actions.append(["modelisation","ref"])
            for i in xrange(self.n_proc-2):
                for change in xrange(5): 
                    self.actions.append(["init_template"])
                    helix=rd.choice([TM for TM in self.TMs])
                    part=rd.choice(["all","N","C"])
                    direction=rd.choice(["x", "y", "z", helix])
                    pivot = int(self.TMs[helix][0])
                    start = pivot-int(self.TMs[helix][1])
                    end   = pivot+int(self.TMs[helix][2])
                    if part=="all":
                        mvt=rd.choice(["translation", "rotation"])
                    else:
                        mvt="rotation"
                        if part=="C":
                            helix="res:"+str(start)+":"+str(pivot)
                        else:
                            helix="res:"+str(pivot)+":"+str(end)
                    if mvt=="translation":
                        d=rd.uniform(-2,2)
                        self.actions.append(['translation', helix, direction, str(d)])
                    elif mvt=="rotation":
                        if direction == helix and part=="all":
                            a=rd.uniform(-20,20)
                        else:
                            a=rd.uniform(-10,10)
                            self.actions.append(['rotation', helix, direction, str(a), "CA:"+str(pivot)+": "])
                
                #CG# hbond restraints added
                for indice_donor, line in enumerate(matrix_ref):
                    value_weight=[]
                    for value,n in enumerate(line):
                        #CG# outset=500
                        #CG# value_weight+=[value for z in xrange(n+outset)]
                        value_weight+=[value]
                    for indice_acceptor in rd.sample(value_weight,3):
                        #CG# self.actions.append(['restraint', "distance", self.list_donor[indice_donor], self.list_acceptor[indice_acceptor], 3, max(2-0.001*line[indice_acceptor], 1.2, 100)])
                        self.actions.append(['restraint', "upperdistance", self.list_donor[indice_donor], self.list_acceptor[indice_acceptor], 6, 3-0.005*line[indice_acceptor]])
                        pass
                self.actions.append(["modelisation",str(step)+"_"+str(i+1)])
            self.set_actions()
            self.run(start_analysis=start, end_analysis=end)
            self.run()
            matrix_ref=self.get_hbond(matrix_ref)

            with open("result_hbond.csv", "w") as fileout:
                for line in matrix_ref:
                    fileout.write("\t".join([str(i) for i in line])+"\n")
            
            name=self.get_results()
            if name=="":
                shutil.copy("blank.pdb", "blank_model/blank_"+str(step)+".pdb")
            else:
                os.rename("blank.pdb", "blank_model/blank_"+str(step)+".pdb")
                shutil.copy(name, "blank.pdb")
                
            print "#####################################################"
            print name
            print "#####################################################"
            
        os.rename("blank.pdb", "blank_model/blank_"+str(step)+".pdb")
        
    ##################################################
    ################## get results ###################
    ##################################################
    
    def get_results(self):
        """ KC - analyse of results """
        
        name=""
        value=-1
        result_list=[]
        for modelisation in self.list_modelisations:
            result_list+=modelisation.get_results()
                    
        name=""
        DOPE="Inf"
        for res in result_list:            
            #CG# if res["DOPE score"]<0.9*self.DOPE_blank and res["molpdf"]<1.5*self.molpdf_blank:
            print "RMSD: ", res["RMSD"]
            print res["DOPE score"], 0.9*self.DOPE_blank
            if res["DOPE score"]<0.9*self.DOPE_blank and res["RMSD"]<3:
                print "ok"
                ok=False
                if res["hbond size"]>value:
                    ok=True
                elif res["hbond size"]==value and res["DOPE score"]<DOPE:
                    ok=True
                elif res["hbond size"]==value and res["DOPE score"]==DOPE and "ref" in res["filename"]:
                    ok=True
                    
                if ok:
                    name=res["filename"]
                    value=res["hbond size"]
                    DOPE=res["DOPE score"]
        return name        
                     
    ##################################################
    ################## get hbonds ####################
    ##################################################
    
    def get_hbond(self, matrix_ref):
        """ KC - get H bond """
        
        for modelisation in self.list_modelisations:
            hbond_list=modelisation.get_hbond()
            #CG# matrix_ref=[[1 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]
            for (acceptor, donor) in hbond_list:
                if not (donor.name=="N" or acceptor.name=="O"):
                    indice_donor=self.dico_donor[donor.get_full_id()[4][0]+":"+str(donor.get_full_id()[3][1])+":"+donor.get_full_id()[2]]
                    indice_acceptor=self.dico_acceptor[acceptor.get_full_id()[4][0]+":"+str(acceptor.get_full_id()[3][1])+":"+acceptor.get_full_id()[2]]
                    matrix_ref[indice_donor][indice_acceptor]+=1
                    
        return matrix_ref

if __name__ == '__main__':
    GPCR_homology_modeling()
 
