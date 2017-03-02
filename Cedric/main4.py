#! -*- coding: utf-8 -*-
import Bio.PDB
from Bio.PDB.Vector import *
from modeller import *
from modeller.automodel import *
## Set all Modeller log output levels
##log.level(output=0, notes=0, warnings=0, errors=0, memory=0)
log.none()
import numpy
from pylab import *
import os
import shutil
import math
import sys
import random as rd
from multiprocessing import Pool
import time

from fonctions import *
import parametres
from modelisation import *

actual_path=os.getcwd()
import gc
gc.enable()


MAKE=True

########################################################################
########################                        ########################
########################        FONCTIONS       ########################
########################                        ########################
########################################################################
def run_multiprocessing(tab):
    
    modelisation,start_analysis,end_analysis=tab
    modelisation.run()
    for model in modelisation.models:
        model.interactions.run(start=start_analysis, end=end_analysis)
    return modelisation

   
    
    

########################################################################
########################                        ########################
########################          MAIN          ########################
########################                        ########################
########################################################################
class GPCR_homology_modeling:
    
    def __init__(self):
        self.archive()
        self.filein = open("parametres", 'r')
        self.actions=[]
        self.TMs={}
        self.n_models=10
        self.pdb_ref=None
        self.sequence_name=None
        self.sequence_filename=None
        self.pdb_template=None
        self.alignment=True
        self.n_proc=1
        self.list_modelisations=None
        
        self.parsing_parameter_file()
        
        ## multiprocessing reservation
        self.pool = Pool(int(self.n_proc-1))        
        
        if self.alignment:
            self.seq_alignement()

        self.seq_template,self.seq_model=get_sequences('resultat.ali', self.sequence_name, self.pdb_template)
        self.conservation=[AA1==AA2 for AA1,AA2 in zip(self.seq_template,self.seq_model)]
        
        self.blank_model_modelisation()
        self.set_actions()
        self.run()
        
        
        ##self.hbond_optimisation()
        ##self.helix_optimisation()
        ##self.helix_hbond_optimisation()
            
        
        
        
    ##################################################
    ###########         archive        ###############
    ##################################################
    
    def archive(self):
        dirname="archive_"+time.strftime('%d-%m-%y_%HH%Mmin%S',time.localtime())
        os.mkdir(dirname)
        os.system("mv blank_model ?_? ??_? ???_? ????_? ?????_? "+dirname+"/. > /dev/null 2>&1")
        
    ##################################################
    ########### parsing parameter file ###############
    ##################################################
    
    def parsing_parameter_file(self):
        
        for line in self.filein:
            line=line.strip().replace('\t',' ')
            if line !="" and not line.startswith("#"):
                ## Split and remove empty element in tab
                tab=line.split(" ")
                tab=[ele for ele in tab if ele!=""]
                
                if tab[0]=="action":
                    self.actions.append(tab[1:])
                elif tab[0].startswith("TM"):
                    self.TMs[tab[0]]=tab[1:]
                elif tab[0]=="sequence":
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


        assert isinstance(self.sequence_name, str), "ERROR: sequence_name have to be a string"
        assert isinstance(self.sequence_filename, str), "ERROR: sequence_filename have to be a string"
        assert isinstance(self.pdb_template, str), "ERROR: pdb_template have to be a string"
        assert (isinstance(self.pdb_ref, str) or self.pdb_ref is None), "ERROR: pdb_ref have to be a string"
        assert (isinstance(self.n_models, int) or self.n_models>1), "ERROR: n_models have to be an integer > 1"
        assert os.path.isfile(self.sequence_filename), 'ERROR: sequence file '+self.sequence_filename+' not exist'
        assert os.path.isfile(self.pdb_template+'.pdb'), 'ERROR: sequence file '+self.pdb_template+' not exist'
        if self.pdb_ref is not None:
            assert os.path.isfile(self.pdb_ref+'.pdb'), 'ERROR: sequence file '+self.pdb_ref+' not exist'


    ##################################################
    ################# Blank model ####################
    ##################################################
    def blank_model_modelisation(self):

        print "### Blank model ###"
        env = environ()
        env.io.hetatm = True ## accept BLK residues and ligand 
        a = automodel(env, alnfile='resultat.ali',
                knowns=self.pdb_template, sequence=self.sequence_name,
                assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model = 1
        a.ending_model = 1
        a.make()
        
        self.DOPE_blank=a.outputs[0]["DOPE score"]
        self.molpdf_blank=a.outputs[0]["molpdf"]
        
        os.mkdir("blank_model")
        shutil.copyfile(self.sequence_name+".B99990001.pdb", "blank_model/blank.pdb")
        shutil.move(self.sequence_name+".B99990001.pdb", "blank.pdb")

        self.pdb_blank="blank_model/blank"
        self.pdb_filename_blank = "%s.pdb" % self.pdb_blank

        self.blank_model=my_model(self.pdb_filename_blank, self.seq_model, self.seq_template, self.sequence_name, None)
        self.blank_model.color_conservation() ## make the pdb file colored by conservation

        ## Same for template
        my_model(self.pdb_template+".pdb", self.seq_template, self.seq_model, self.pdb_template, None, _type="xray").color_conservation()
        

        move_files("blank_model")
        #clean_tempfile()
        
        print "INNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        
        env.io.hetatm = True ## accept BLK residues and ligand 
        ## alignment between blank model structure and blank model sequence for next modelisation
        aln = alignment(env)
        mdl = model(env, file='blank.pdb')
        aln.append_model(mdl, align_codes='blank_model', atom_files='blank.pdb')
        aln.append(file=self.sequence_filename, align_codes=self.sequence_name)
        aln.align2d()
        aln.write(file='resultat_final.ali', alignment_format='PIR')


        
        
    ############## seq_alignement ##############    
                
    def seq_alignement(self):
        print "### Alignment ###"
        env = environ()
        env.io.hetatm = True ## accept BLK residues and ligand 
        #-- alignement
        aln = alignment(env)
        mdl = model(env, file=self.pdb_template)
        aln.append_model(mdl, align_codes=self.pdb_template, atom_files=self.pdb_template+'.pdb')
        aln.append(file=self.sequence_filename, align_codes=self.sequence_name)
        aln.align2d()
        aln.write(file='resultat.ali', alignment_format='PIR')
        ## aln.to_profile().write(file='alntoprof.prf', profile_format='TEXT')
        
    
    ##################################################
    ################### ACTIONS ######################
    ##################################################
    def set_actions(self):
        # Open the pdb file
        pdb_code = self.pdb_template
        pdb_filename = "%s.pdb" % pdb_code
        structure=Bio.PDB.PDBParser().get_structure("blank", "blank.pdb")
        model=structure[0]

        restraints=[]
        ## list_modelisations to run and analyses
        self.list_modelisations=[]
        ## action_temp to remind what's happen in each modelisation
        actions_temp=[]

        for action in self.actions:
            actions_temp.append(action)
            if action[0] in ["rotation", "translation"]:                
                # Search the first point
                if action[1].startswith("TM"):
                    ## Search begin and end of TM
                    debut=int(self.TMs[action[1][:3]][0])-int(self.TMs[action[1][:3]][1])
                    fin = int(self.TMs[action[1][:3]][0])+int(self.TMs[action[1][:3]][2])
                elif action[1].startswith("res"):
                    ## Search begin and end of TM
                    debut=int(action[1].split(":")[1])
                    fin = int(action[1].split(":")[2])
                elif action[1]=="protein":
                    debut=-1
                    fin="Inf"
                elif action[1]=="ligand":
                    print "a implementer"
                
                # Generate some data that lies along a line
                selection=[atom for atom in model.get_atoms() if debut<=int(atom.get_full_id()[3][1])<=fin]
                data = np.array([atom.get_coord() for atom in selection])
                # Calculate the mean of the points, i.e. the 'center' of the cloud
                axe_point_1 = data.mean(axis=0)
                
                # Define axe
                # Search the first point
                if len(action[2].split(','))==3:
                    axe=[float(var) for var in action[2].split(',')]
                elif action[2] in ["x","y","z"]:
                    dico={"x":[1,0,0], "y":[0,1,0], "z":[0,0,1]}
                    axe=dico[action[2]]
                else:
                    if action[2].startswith("TM"):
                        ## Search begin and end of TM
                        debut=int(self.TMs[action[2][:3]][0])-int(self.TMs[action[2][:3]][1])
                        fin = int(self.TMs[action[2][:3]][0])+int(self.TMs[action[2][:3]][2])
                    elif action[2].startswith("res"):
                        ## Search begin and end of TM
                        debut=int(action[2].split(":")[1])
                        fin = int(action[2].split(":")[2])
                    elif action[2]=="protein":
                        debut=-1
                        fin="Inf"
                    elif action[2]=="ligand":
                        print "a implementer"
                    
                    # Generate some data that lies along a line
                    data = np.array([atom.get_coord() for atom in model.get_atoms() if debut<=int(atom.get_full_id()[3][1])<=fin])
                    # Calculate the mean of the points, i.e. the 'center' of the cloud
                    axe_point_2 = data.mean(axis=0)
                    
                    
                    if np.all(axe_point_2==axe_point_1):
                        axe=determine_axe(selection, None)
                    elif action[2].endswith("_perp"):
                        axe=determine_axe(selection, axe_point_2, _type="perp")
                    
                    else:
                        axe=determine_axe(selection, axe_point_2, _type="normal")
                
                ## Mouvement
                if action[0] == "rotation":
                    pivot_tab=action[4].split(":")
                    pivot=model[pivot_tab[2]][int(pivot_tab[1])][pivot_tab[0]].get_coord()
                    for atom in selection:
                        atom.set_coord(rotation(vecteur_coord=atom.get_coord(), axe=np.array(axe), angle=2*pi/360*float(action[3]), centre=pivot))
                
                elif action[0] == "translation":
                    for atom in selection:
                        atom.set_coord(translation(vecteur_coord=atom.get_coord(),valeur=float(action[3]), axe=axe))
            
            elif action[0]=="modelisation":
                try:
                    os.mkdir(action[1])
                except:
                    pass
                io=Bio.PDB.PDBIO()
                io.set_structure(structure)
                io.save(action[1]+"/blank.pdb")
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
            
            
            elif action[0]=="restraint":
                restraints.append(action[1:])
                    
            elif action[0]=="clear_restraint":
                restraints=[]
                    
            elif action[0]=="init_template":
                structure=Bio.PDB.PDBParser().get_structure("blank", "blank.pdb")
                model=structure[0]
                restraints=[]
                actions_temp=[]




    ##################################################
    ################### RUN JOBS #####################
    ##################################################
    def run(self, start_analysis=0, end_analysis=-1):
        """
        for modelisation in self.list_modelisations:
            modelisation.run()
            for model in modelisation.models:
                model.interactions.run(start=start_analysis, end=end_analysis)
            

        """
        self.list_modelisations=self.pool.map(run_multiprocessing, [[modelisation, start_analysis, end_analysis] for modelisation in self.list_modelisations] )


        
        return
        

    ##################################################
    ############## helix Optimisation ################
    ##################################################


    def helix_optimisation(self):
        """
        The helix optimisation follow this steps:
        - optimisation in 3D of the helix
        - optimisation in 3D of the half-helix around the pivot
        """
        
        stop=False
        step=0 
        
        while not stop:       
            self.actions=[]
            step+=1
            self.actions.append(["modelisation","ref"])
            for i in xrange(self.n_proc-2):

                for change in xrange(5):
                    self.actions.append(["init_template"])
                    helix=rd.choice([TM for TM in self.TMs])
                    pivot = int(self.TMs[helix][0])
                    start = pivot-int(self.TMs[helix][1])
                    end   = pivot+int(self.TMs[helix][2])
                    part=rd.choice(["all","N","C"])
                    print helix, part
                    if part=="all":
                        mvt=rd.choice(["translation", "rotation"])
                    else:
                        mvt="rotation"
                        if part=="C":
                            helix="res:"+str(start)+":"+str(pivot)
                        else:
                            helix="res:"+str(pivot)+":"+str(end)
                            
                        
                    direction=rd.choice(["x", "y", "z", helix])
                    
                    if mvt=="translation":
                        d=rd.uniform(-2,2)
                        self.actions.append(['translation', helix, direction, str(d)])
                    elif mvt=="rotation":
                        if direction == helix and part=="all":
                            a=rd.uniform(-20,20)
                        else:
                            a=rd.uniform(-10,10)
                            self.actions.append(['rotation', helix, direction, str(a), "CA:"+str(pivot)+": "])
                
                ##self.actions.append(["restraint", "helix", "11", "21"])
                self.actions.append(["modelisation",str(step)+"_"+str(i+1)])
                    
                    
            ## echantillonage
            if step=="Inf":
                stop=True
            self.set_actions()
            ##self.run(start_analysis=start, end_analysis=end)
            self.run()
            
            
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
    ############## hbond Optimisation ################
    ##################################################

    def hbond_optimisation(self):

        #action  restraint		distance	CA:100:A	CA:200:A	10	0.1
        self.list_donor=[donor.get_full_id()[4][0]+":"+str(donor.get_full_id()[3][1])+":"+donor.get_full_id()[2] for donor in self.blank_model.get_donors() if donor.get_full_id()[4][0]!="N"]
        self.list_acceptor=[acceptor.get_full_id()[4][0]+":"+str(acceptor.get_full_id()[3][1])+":"+acceptor.get_full_id()[2] for acceptor in self.blank_model.get_acceptors() if acceptor.get_full_id()[4][0]!="O"]
                
        ## dico to easy acces
        self.dico_donor={name:number for number, name in enumerate(self.list_donor)}
        self.dico_acceptor={name:number for number, name in enumerate(self.list_acceptor)}
        
        
        matrix_ref=[[0 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]
        
        
        for step in xrange(1000):
            self.actions=[]
            for n_model in xrange(self.n_proc-1):
                self.actions.append(["init_template"])
                ## random choice of restraintes which will be done
                for indice_donor, line in enumerate(matrix_ref):
                    ##[2,5,1]==>[0,0,1,1,1,1,1,2]
                    value_weight=[]
                    for value,n in enumerate(line):
                        outset=100
                        value_weight+=[value for i in xrange(n+outset)]
                    ##random choice
                    for indice_acceptor in rd.sample(value_weight,5):
                        ##self.actions.append(['restraint', "distance", self.list_donor[indice_donor], self.list_acceptor[indice_acceptor], 3, max(2-0.001*line[indice_acceptor], 1.2, 100)])
                        self.actions.append(['restraint', "distance", self.list_donor[indice_donor], self.list_acceptor[indice_acceptor], 3, 1.5-0.005*line[indice_acceptor]])
                self.actions.append(["modelisation", str(step+1)+"_"+str(n_model+1)])
            
            self.set_actions()
            self.run()

            matrix_ref=self.get_hbond(matrix_ref)
            
            ## save hbond matrix
            with open("result_hbond.csv", "w") as fileout:
                for line in matrix_ref:
                    fileout.write("\t".join([str(i) for i in line])+"\n")            
            
            
           
        
                    
                    
                    
                    
                    
    ##################################################
    ######### helix and hbond Optimisation ###########
    ##################################################


    def helix_hbond_optimisation(self):
        """
        The helix optimisation follow this steps:
        - optimisation in 3D of the helix
        - optimisation in 3D of the half-helix around the pivot
        """
        
        ### INITIALISATION ###
        
        stop=False
        step=0 
        
        self.list_donor=[donor.get_full_id()[4][0]+":"+str(donor.get_full_id()[3][1])+":"+donor.get_full_id()[2] for donor in self.blank_model.get_donors() if donor.get_full_id()[4][0]!="N"]
        self.list_acceptor=[acceptor.get_full_id()[4][0]+":"+str(acceptor.get_full_id()[3][1])+":"+acceptor.get_full_id()[2] for acceptor in self.blank_model.get_acceptors() if acceptor.get_full_id()[4][0]!="O"]
                
        ## dico to easy acces
        self.dico_donor={name:number for number, name in enumerate(self.list_donor)}
        self.dico_acceptor={name:number for number, name in enumerate(self.list_acceptor)}
        
        matrix_ref=[[0 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]

        
        for step in xrange(1,101,1):       
            self.actions=[]
            self.actions.append(["modelisation","ref"])
            for i in xrange(self.n_proc-2):## for each model
                ## apply 5 change in model
                for change in xrange(5): 
                    self.actions.append(["init_template"])
                    helix=rd.choice([TM for TM in self.TMs])
                    pivot = int(self.TMs[helix][0])
                    start = pivot-int(self.TMs[helix][1])
                    end   = pivot+int(self.TMs[helix][2])
                    part=rd.choice(["all","N","C"])
                    if part=="all":
                        mvt=rd.choice(["translation", "rotation"])
                    else:
                        mvt="rotation"
                        if part=="C":
                            helix="res:"+str(start)+":"+str(pivot)
                        else:
                            helix="res:"+str(pivot)+":"+str(end)
                            
                        
                    direction=rd.choice(["x", "y", "z", helix])
                    
                    if mvt=="translation":
                        d=rd.uniform(-2,2)
                        self.actions.append(['translation', helix, direction, str(d)])
                    elif mvt=="rotation":
                        if direction == helix and part=="all":
                            a=rd.uniform(-20,20)
                        else:
                            a=rd.uniform(-10,10)
                            self.actions.append(['rotation', helix, direction, str(a), "CA:"+str(pivot)+": "])
                
                ## add Hbond restraints
                for indice_donor, line in enumerate(matrix_ref):
                    ##[2,5,1]==>[0,0,1,1,1,1,1,2]
                    value_weight=[]
                    for value,n in enumerate(line):
                        #outset=500
                        #value_weight+=[value for z in xrange(n+outset)]
                        value_weight+=[value]
                    ##random choice
                    for indice_acceptor in rd.sample(value_weight,3):
                        ##self.actions.append(['restraint', "distance", self.list_donor[indice_donor], self.list_acceptor[indice_acceptor], 3, max(2-0.001*line[indice_acceptor], 1.2, 100)])
                        self.actions.append(['restraint', "upperdistance", self.list_donor[indice_donor], self.list_acceptor[indice_acceptor], 6, 3-0.005*line[indice_acceptor]])
                        pass
                
                
                
                ## modelisation
                self.actions.append(["modelisation",str(step)+"_"+str(i+1)])
                    
                    
            ## echantillonage
            if step=="Inf":
                stop=True
            self.set_actions()
            ##self.run(start_analysis=start, end_analysis=end)
            self.run()
            matrix_ref=self.get_hbond(matrix_ref)
            
            ## save hbond matrix
            with open("result_hbond.csv", "w") as fileout:
                for line in matrix_ref:
                    fileout.write("\t".join([str(i) for i in line])+"\n")
            
            ## move new blank
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
    ################## Get Results ###################
    ##################################################
    
    def get_results(self):
        name=""
        value=-1
        result_list=[]
        for modelisation in self.list_modelisations:
            result_list+=modelisation.get_results()
                    
        name=""
        DOPE="Inf"
        for res in result_list:            
            ##if res["DOPE score"]<0.9*self.DOPE_blank and res["molpdf"]<1.5*self.molpdf_blank:
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
    ################## Get hbonds ####################
    ##################################################
    
    def get_hbond(self, matrix_ref):
        for modelisation in self.list_modelisations:
            hbond_list=modelisation.get_hbond()
            ##matrix_ref=[[1 for i in xrange(len(self.list_acceptor))] for j in xrange(len(self.list_donor))]
            
            for (acceptor, donor) in hbond_list:
                if not (donor.name=="N" or acceptor.name=="O"):
                    indice_donor=self.dico_donor[donor.get_full_id()[4][0]+":"+str(donor.get_full_id()[3][1])+":"+donor.get_full_id()[2]]
                    indice_acceptor=self.dico_acceptor[acceptor.get_full_id()[4][0]+":"+str(acceptor.get_full_id()[3][1])+":"+acceptor.get_full_id()[2]]
                    matrix_ref[indice_donor][indice_acceptor]+=1
                    
        return matrix_ref



if __name__ == '__main__':
    GPCR_homology_modeling()

