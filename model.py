#! -*- coding: utf-8 -*-

#KC# packages import
import os
import Bio.PDB
import numpy as np
from modeller import *
from modeller.automodel import *

#CG# set all Modeller log output levels
log.level(output=0, notes=0, warnings=0, errors=0, memory=0)

#KC# own files import
import parametres
from fonctions import *
from interactions import *

#KC# working directory
workdir=os.getcwd()

########################################################################
########################                        ########################
########################        my model        ########################
########################                        ########################
######################################################################## --> CG

class my_model:
    
    #KC#CG# class constructor
    def __init__(self, filename, seq_model=None, seq_template=None, sequence_name="Unknown", list_helix=None, DOPE=None, molpdf=None, GA341=None, _type=None, _test=False, ASA_per=float(30)):
        
        print "MODEL = "+sequence_name
        
        #KC#CG# parameters tests
        assert _type in [None,"blank","xray"], "_type ERROR"
        
        #KC#CG# variables declaration
        self.DOPE=DOPE
        self.molpdf=molpdf
        self.GA341=GA341
        self._type=_type
        self._test=_test
        self.filename=filename
        self.seq_model=seq_model
        self.seq_template=seq_template
        self.sequence_name=sequence_name
        pdb_code=self.filename[:-4]
        pdb_filename = "%s.pdb" % pdb_code
        self.ASA_per=ASA_per
        self.RMSD=None #CG# tab with RMSD on n_models for each atom
        self.asa=None #CG# tab with ASA for each residues for each models
        self.angle_dihedres=None #CG# tab with difference of dihedral angle for each residues for each models
        
        #KC#CG# PDB structure, model, atoms, alpha carbons and residus collected
        self.structure=Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename) 
        self._Model=self.structure[0]
        self.atoms=[atom for atom in self._Model.get_atoms()]
        self.CA=[atom for atom in self._Model.get_atoms() if atom.get_name()=="CA"]
        self.res=[res for res in self._Model.get_residues() if res.get_resname() not in parametres.list_lig_name]

        #KC# functions launching
        if self._test==False:
            if self._type=="xray":
                self.renum_dict = self.dict_renumbering()
                self.find_don_acc()
                self.find_ionic_don_acc()   
                self.ASA_calcul(sequence_name=self._type)
                self.interactions=interactions(self, renum_dict=self.renum_dict)
                self.interactions.Hbonds_dist(self.donors, self.filename[:-4])
            elif self._type=="blank":
                if (list_helix==None):
                    self.find_secondary_structure()
                else:
                    self.res_list_helix=list_helix
                self.orient_protein()
                self.calc_box_size()
                self.color_conservation()
                self.search_zones()
                self.color_zones()
                self.find_don_acc()
                self.find_ionic_don_acc()  
                self.ASA_calcul(sequence_name=self._type)
                self.color_core()
                self.interactions=interactions(self)
                self.interactions.Hbonds_dist(self.donors, self.filename[:-4])
            else:
                if (list_helix==None):
                    self.find_secondary_structure()
                else:
                    self.res_list_helix=list_helix
                self.orient_protein()
                self.calc_box_size()
                self.search_zones()
                self.find_don_acc()
                self.find_ionic_don_acc()  
                self.ASA_calcul(sequence_name=self._type)
                self.interactions=interactions(self, _color=False, _graph=False)
        else:
            self.renum_dict = self.dict_renumbering()
            self.find_don_acc()
            self.find_ionic_don_acc()   
            self.ASA_calcul(sequence_name=self._type)
            self.interactions=interactions(self, renum_dict=self.renum_dict, _color=False, _test=self._test)
            self.interactions.Hbonds_dist(self.donors, self.filename[:-4])
            
    ##################################################
    ############# get and set functions ##############
    ################################################## --> CG

    def get_CA(self):
        return self.CA
        
    def get_res(self):
        return self.res
        
    def get_atoms(self):
        return self.atoms
        
    def get_filename(self):
        return self.filename
        
    def get_donors(self):
        return self.donors
    
    def get_acceptors(self):
        return self.acceptors
    
    def get_ionic_donors(self):
        return self.ionic_donors
    
    def get_ionic_acceptors(self):
        return self.ionic_acceptors

    ##################################################
    ################ dict renumbering ################
    ################################################## --> KC
    
    def dict_renumbering(self):
        """ KC - dictionary for correspondance between numbers of atoms """

        liste_num = os.popen("cat "+self.filename+" | tr -s ' ' | grep ^ATOM | cut -d ' ' -f 6 | uniq", "r").read().split()
        liste_num = [int(i) for i in liste_num]
        liste_renum = list(range(len(self.res)))
        return dict(zip(liste_num, liste_renum))
    
    ##################################################
    ########### find secondary structure #############
    ################################################## --> CG
    
    def find_secondary_structure(self):
        """ KC - executable stride used to find secondary structure """

        os.system('./stride '+self.filename+'>stride.out')
        filein=open('stride.out','r')
        SEQ=''
        self.STR=''
        
        #KC#CG# sequence and secondary structure prediction collected
        for line in filein:
            if line.startswith('SEQ'):
                SEQ+=line[10:60].strip()
            if line.startswith('STR'):
                self.STR+=line[10:10+len(SEQ)-len(self.STR)]

        #KC#CG# list of the position of helix in the sequence
        self.res_list_helix=[]
        for i,S in enumerate(self.STR):
            if S=='H':
                self.res_list_helix.append(i)
        
        os.remove('stride.out')
        
        #KC#CG# C-ter and N-ter parts collected
        self.C_ter=[]
        self.N_ter=[]
        i=0
        while self.STR[i]!="H":
            self.N_ter.append(i)
            i+=1
        i=len(self.STR)-1
        while self.STR[i]!="H":
            self.C_ter.append(i)
            i-=1
        
        return

    ##################################################
    ################ orient protein ##################
    ################################################## --> CG
           
    def orient_protein(self):
        """ KC - protein right up according to z and centered on 0 """
        
        #CG# list with atoms in helix
        list_atoms=[atom for atom in self.atoms if atom.get_full_id()[3][1] in self.res_list_helix]
        axe_from=determine_axe(list_atoms, None)
        axe_to=[0,0,1]
        angle=calcul_angle(axe_from, [0,0,0], axe_to)
        
        #CG# search the axe between witch one the protein have to rotate
        x1,y1,z1=axe_from
        x2,y2,z2=axe_to
        x3=1        
        #CG# development after calculating scalar products
        y3=(z1/z2*x2*x3-x1*x3)/(y1-z1/z2*y2)
        z3=-(x2*x3+y2*y3)/z2
        axe3=[x3,y3,z3]
        
        #CG# calculate the mean of the points, i.e. the 'center' of the cloud
        data=np.array([atom.get_coord() for atom in list_atoms])
        pivot=data.mean(axis=0)

        #KC#CG# rotation and translation
        for atom in self.atoms:
            atom.set_coord(rotation(vecteur_coord=atom.get_coord(), axe=np.array(axe3), angle=-2*np.pi/360*angle, centre=pivot))
            atom.set_coord(translation(vecteur_coord=atom.get_coord(), valeur=distance(pivot, [0,0,0]), axe=np.array(pivot)))
                
        if self.atoms[0].get_coord()[2]<0:
            for atom in self.atoms:
                atom.set_coord(rotation(vecteur_coord=atom.get_coord(), axe=np.array([1,0,0]), angle=np.pi, centre=[0,0,0]))
        
        return
        
    ##################################################
    ################# calc box size ##################
    ################################################## --> CG
        
    def calc_box_size(self, _type="helix"):
        """ KC - size of the box with helix calculated """
        
        if _type == "helix":
            tab=[atom.get_coord() for atom in self.atoms if atom.get_full_id()[3][1] in self.res_list_helix]
        elif _type == "protein":
            tab=[atom.get_coord() for atom in self.atoms]

        #KC#CG# min/max selected
        self.x_min=min([ele[0] for ele in tab])
        self.x_max=max([ele[0] for ele in tab])
        self.y_min=min([ele[1] for ele in tab])
        self.y_max=max([ele[1] for ele in tab])
        self.z_min=min([ele[2] for ele in tab])
        self.z_max=max([ele[2] for ele in tab])

        return 

    ##################################################
    ################ search top down #################
    ################################################## --> CG (not used)
    
    def search_top_down(self):
        """ KC - links top and down of helix search """
        
        self.top=[]
        self.down=[]
        
        for (i,_type),res in zip(enumerate(self.STR),self.res):
            if _type not in "H" and i not in self.C_ter+self.N_ter:
                if res['CA'].get_coord()[2]<0.5*self.z_min:
                    self.down.append(i)
                elif res['CA'].get_coord()[2]>0.5*self.z_max:
                    self.top.append(i)
                    
        return

    ##################################################
    ############## color conservation ################
    ################################################## --> CG
    
    def color_conservation(self):
        """ KC - file conservation created """
        
        i=0
        for aa1, aa2 in zip(self.seq_model,self.seq_template):
            if aa1 in ['-','/','*'] and aa2 in ['-','/','*']:
                pass
            elif aa1 in ['-','/']: #CG# deletion in the template sequence
                pass
            elif aa2 in ['-','/']: #CG# insertion in the model sequence
                for atom in self.get_res()[i]:
                    if self._type=="xray":
                        atom.set_bfactor(25)
                    else:
                        atom.set_bfactor(75)
                i+=1
            else:
                assert aa1 == Bio.PDB.Polypeptide.three_to_one(self.get_res()[i].resname), "\n\nPlease check your alignment file...\nCombinations with special characters ( - / * ) and letters ( L ) :\n-   /   *   -   /   -   /   L   L\n/   -   *   -   /   L   L   -   /\nFile %s_conservation.pdb not saved !\n" % self.sequence_name
                if aa1==aa2: #CG# identity
                    for atom in self.get_res()[i]:
                        atom.set_bfactor(0)
                elif (aa1,aa2) in [('F','Y'),('I','L'),('S','T'),('D','E'),('R','K'),('A','G')] or (aa2,aa1) in [('F','Y'),('I','L'),('S','T'),('D','E'),('R','K'),('A','G')]: #KC# similarity
                    for atom in self.get_res()[i]:
                        atom.set_bfactor(50)
                else: #CG# difference
                    for atom in self.get_res()[i]:
                        atom.set_bfactor(100)
                i+=1
        io=Bio.PDB.PDBIO()
        io.set_structure(self.structure)
        io.save(self.filename[:-4]+"_conservation.pdb")
        
    ##################################################
    ################## search zones ##################
    ################################################## --> KC
    
    def search_zones(self):
        """ KC - 3 zones search """
        
        self.intra=[]
        self.mb=[]
        self.extra=[]
        
        #KC# protein separated into 3 parts
        limit_max = self.z_max-(self.z_max-self.z_min)/3
        limit_min = self.z_min+(self.z_max-self.z_min)/3
        
        for res in self.res:
            if res['CA'].get_coord()[2]<limit_min:
                self.intra.append(res)
            elif res['CA'].get_coord()[2]>limit_max:
                self.extra.append(res)
            else :
                self.mb.append(res)

        return

    ##################################################
    ############## color conservation ################
    ################################################## --> KC
    
    def color_zones(self):
        """ KC - color by zones file created """
        
        for atom in self.atoms :
            if atom.get_parent() in self.intra :
                atom.set_bfactor(0)
            elif atom.get_parent() in self.mb:
                atom.set_bfactor(50)
            else :
                atom.set_bfactor(100)
                
        io=Bio.PDB.PDBIO()
        io.set_structure(self.structure)
        io.save(self.filename[:-4]+"_zones.pdb")
        
    ##################################################
    ################## find don acc ##################
    ################################################## --> CG
    
    def find_don_acc(self):
        """ KC - loop on the residue to find donors ans acceptors  """
    
        self.acceptors=[]
        self.donors=[]
        
        #KC#CG# donors or acceptors atoms added for main and side chains
        for res in self.res:
            if res.resname in parametres.donors_dict:
                for atom_name in parametres.donors_dict[res.resname]:
                    try:
                        self.donors.append(res[atom_name])
                    except:
                        pass
            if res.resname in parametres.acceptors_dict:
                for atom_name in parametres.acceptors_dict[res.resname]:
                    try:
                        self.acceptors.append(res[atom_name])
                    except:
                        pass

    ##################################################
    ############### find ionic don acc ###############
    ################################################## --> KC
    
    def find_ionic_don_acc(self):
        """ KC - loop on the residue to find donors ans acceptors  """
    
        self.ionic_acceptors=[]
        self.ionic_donors=[]
        
        #KC# ionic donors or acceptors atoms added for main and side chains
        for res in self.res:
            if res.resname in parametres.ionic_donors_dict:
                for atom_name in parametres.ionic_donors_dict[res.resname]:
                    try:
                        self.ionic_donors.append(res[atom_name])
                    except:
                        pass
            if res.resname in parametres.ionic_acceptors_dict:
                for atom_name in parametres.ionic_acceptors_dict[res.resname]:
                    try:
                        self.ionic_acceptors.append(res[atom_name])
                    except:
                        pass
        
    ##################################################
    ################### ASA calcul ###################
    ################################################## --> CG/KC
    
    def ASA_calcul(self, sequence_name=None):
        """ KC - model ASA calcul """
        
        ######################################################################################
        #KC#CG# previous method with Naccess --> to be compared
        """os.system('/home/ckaren/Documents/Codes/Naccess/naccess '+self.filename)
        pdb_code=self.filename[:-4]
        pdb_out_filename = "%s_asa.pdb" % pdb_code
        if sequence_name==None:
            pdb_filename = "%s.asa" % self.sequence_name
        else:
            pdb_filename = "%s.asa" % sequence_name
    
        filein=open(pdb_filename,'r')
        self.asa=[float(line[64:69].strip()) for line in filein if line[12:16].strip()!='OXT']    
        filein.close()
        
        maxi=max(self.asa)
        mini=min(self.asa)
        self.asa=[(bfactor-mini)*100/(maxi-mini) for bfactor in self.asa]
    
        filein=open(pdb_filename,'r')
        fileout=open(pdb_out_filename,'w')

        for bfactor,line in zip(self.asa, filein):
            fileout.write(line[:64]+str(bfactor)[:4]+'\n')
        
        fileout.close()
        filein.close()"""
        ######################################################################################
        
        #KC# new method with modeller function
        #KC# get topology library for radii and the model without waters and HETATMs:
        env = environ()
        env.io.atom_files_directory = [os.getcwd()]
        env.io.hetatm = False
        env.io.water = False
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        mdl = model(env, file=self.filename)
        
        #KC# calculate and write PSA score to file
        myedat = energy_data()
        myedat.radii_factor = 1.0 #KC# default = 0.82 (for soft-sphere restraints)
        mdl.write_data(file=self.filename[:-4], edat=myedat, output='PSA')
        os.system('rm *.sol > /dev/null 2>&1')
        self.search_core()
        
    ##################################################
    ################### search core ##################
    ################################################## --> KC

    def search_core(self):
        """ KC - core research in blank and template model """

        #KC# file reading to get access res in list
        self.liste_access = os.popen("cat "+self.filename[:-4]+".psa | tr -s ' ' | grep -v \# | cut -d ' ' -f 2,6,10,14", "r").read()
        self.liste_access = self.liste_access.split("\n")
        del self.liste_access[-1]
        self.res_core=[]

        #KC# each line of the file analysed
        for line in self.liste_access :
            line = line.split()
            if float(line[2]) < self.ASA_per :
                self.res_core.append(line[0])
        return
        
    ##################################################
    ################### color core ###################
    ################################################## --> KC

    def color_core(self):
        """ KC - core coloration in blank and template model """
        
        #KC# bfactor changed according to the position of the residue (in the core or not)
        for atom in self.atoms :
            if str(atom.get_parent().get_full_id()[3][1]) in self.res_core:
                atom.set_bfactor(0)
            else :
                atom.set_bfactor(100)
                
        #KC# structure saved
        io=Bio.PDB.PDBIO()
        io.set_structure(self.structure)
        io.save(self.filename[:-4]+"_core.pdb")
        return
        
    ##################################################
    ################ get RMSD template ###############
    ################################################## --> CG/KC
    
    def get_RMSD_template(self, template):
        """ KC - RMSD score calculated between model and template """
        
        self.atoms_model1=[]
        self.atoms_template=[]
        
        #KC# aligned residues selected
        num_template,num_model=get_aligned_residues('res_align.ali')
        for residue in self.res:
            if residue.get_full_id()[3][1] in num_model:
                for atom in residue:
                    self.atoms_model1.append(atom)  
        for residue in template.res:
            if residue.get_full_id()[3][1] in num_template:
                for atom in residue:
                    self.atoms_template.append(atom)

        #KC# length test
        assert len([atom for atom in self.atoms_model1 if atom.name=="CA"])==len([atom for atom in self.atoms_template if atom.name=="CA"]), "ERROR length CA atoms"
            
        #KC# superposition of aligned residues and RMSD calculation
        super_imposer=Bio.PDB.Superimposer()
        super_imposer.set_atoms([atom for atom in self.atoms_template if atom.name=="CA"], [atom for atom in self.atoms_model1 if atom.name=="CA"])
        #super_imposer.apply(self.atoms)
        
        return (float(super_imposer.rms))

    ##################################################
    ############# get RMSD template core #############
    ################################################## --> KC
    
    def get_RMSD_template_core(self, template):
        """ KC - RMSD score calculated between model and template in the protein core """
        
        at_model=[]
        at_template=[]

        for (at1, at2) in zip([atom for atom in self.atoms_model1 if atom.name=="CA"], [atom for atom in self.atoms_template if atom.name=="CA"]):
            if str(at1.get_parent().get_full_id()[3][1]) in self.res_core:
                at_model.append(at1)
                at_template.append(at2)

        #KC# length test
        assert len(at_model)==len(at_template), "ERROR length CA atoms"
            
        #KC# superposition of aligned residues and RMSD calculation
        super_imposer=Bio.PDB.Superimposer()
        super_imposer.set_atoms(at_template, at_model)
        #super_imposer.apply(self.atoms)
        
        return (float(super_imposer.rms))

    ##################################################
    ########### get RMSD template helices ############
    ################################################## --> KC
    
    def get_RMSD_template_helices(self, template):
        """ KC - RMSD score calculated between model and template in the protein core """
        
        at_model=[]
        at_template=[]

        for (at1, at2) in zip([atom for atom in self.atoms_model1 if atom.name=="CA"], [atom for atom in self.atoms_template if atom.name=="CA"]):
            if at1.get_parent().get_full_id()[3][1] in self.res_list_helix:
                at_model.append(at1)
                at_template.append(at2)

        #KC# length test
        assert len(at_model)==len(at_template), "ERROR length CA atoms"
            
        #KC# superposition of aligned residues and RMSD calculation
        super_imposer=Bio.PDB.Superimposer()
        super_imposer.set_atoms(at_template, at_model)
        #super_imposer.apply(self.atoms)
        
        return (float(super_imposer.rms))
        
    ##################################################
    ############ get RMSD query structure ###########
    ################################################## --> CG/KC
    
    def get_RMSD_query_structure(self, query_structure):
        """ KC - RMSD score calculated between model and query_structure """

        #KC#CG# alignment between crytal structure and orexin sequence
        """env=environ()
        env.io.hetatm=True #CG# read in HETATM records from template PDBs
        aln=alignment(env)
        mdl=model(env, file='4S0V')
        aln.append_model(mdl, align_codes='4S0V', atom_files='4S0V.pdb')
        aln.append(file=self.sequence_name+".fasta", align_codes=self.sequence_name)
        aln.align2d()
        aln.write(file='structure_sequence.ali', alignment_format='PIR')"""

        self.atoms_model2=[]
        self.atoms_query=[]
        
        #KC# aligned residues selected
        num_query,num_model=get_aligned_residues('structure_sequence.ali')
        for residue in self.res:
            if residue.get_full_id()[3][1] in num_model:
                for atom in residue:
                    self.atoms_model2.append(atom)  
        for residue in query_structure.res:
            if residue.get_full_id()[3][1] in num_query:
                for atom in residue:
                    self.atoms_query.append(atom)

        #KC# length test
        assert len([atom for atom in self.atoms_model2 if atom.name=="CA"])==len([atom for atom in self.atoms_query if atom.name=="CA"]), "ERROR length CA atoms"

        #KC# superposition of aligned residues and RMSD calculation
        super_imposer=Bio.PDB.Superimposer()
        super_imposer.set_atoms([atom for atom in self.atoms_query if atom.name=="CA"], [atom for atom in self.atoms_model2 if atom.name=="CA"])
        #super_imposer.apply(self.atoms)
        return (float(super_imposer.rms))

    ##################################################
    ######### get RMSD query structure core ##########
    ################################################## --> KC
    
    def get_RMSD_query_structure_core(self, query_structure):
        """ KC - RMSD score calculated between model and query structure in the protein core """
        
        at_model=[]
        at_query=[]
        
        for (at1, at2) in zip([atom for atom in self.atoms_model2 if atom.name=="CA"], [atom for atom in self.atoms_query if atom.name=="CA"]):
            if str(at1.get_parent().get_full_id()[3][1]) in self.res_core:
                at_model.append(at1)
                at_query.append(at2)

        #KC# length test
        assert len(at_model)==len(at_query), "ERROR length CA atoms"
            
        #KC# superposition of aligned residues and RMSD calculation
        super_imposer=Bio.PDB.Superimposer()
        super_imposer.set_atoms(at_query, at_model)
        #super_imposer.apply(self.atoms)
        
        return (float(super_imposer.rms))

    ##################################################
    ####### get RMSD query structure helices #########
    ################################################## --> KC
    
    def get_RMSD_query_structure_helices(self, query_structure):
        """ KC - RMSD score calculated between model and query structure in the protein core """
        
        at_model=[]
        at_query=[]

        for (at1, at2) in zip([atom for atom in self.atoms_model2 if atom.name=="CA"], [atom for atom in self.atoms_query if atom.name=="CA"]):
            if at1.get_parent().get_full_id()[3][1] in self.res_list_helix:
                at_model.append(at1)
                at_query.append(at2)

        #KC# length test
        assert len(at_model)==len(at_query), "ERROR length CA atoms"
            
        #KC# superposition of aligned residues and RMSD calculation
        super_imposer=Bio.PDB.Superimposer()
        super_imposer.set_atoms(at_query, at_model)
        #super_imposer.apply(self.atoms)
        
        return (float(super_imposer.rms))
        
    ##################################################
    ############## structure alignment ###############
    ################################################## --> CG
    
    def structure_alignement(self, model_ref):
        """ KC - structure alignment between model and reference model """
        
        self.res_model=[]
        self.res_model_ref=[]
        
        #KC# aligned residues selected
        num_query,num_model=get_aligned_residues('structure_sequence.ali')
        for residue in self.res:
            if residue.get_full_id()[3][1] in num_model:
                self.res_model.append(residue)
        for residue in model_ref.res:
            if residue.get_full_id()[3][1] in num_query:
                self.res_model_ref.append(residue)

        if model_ref.get_filename()!=self.get_filename():
            self.CA_list_helix=[]
            CA_list_ref_helix=[]
            for res1,res2 in zip(self.res_model, self.res_model_ref):
                if res1.get_full_id()[3][1] in self.res_list_helix:
                    for atom in res1 :
                        if atom.name=="CA":
                            self.CA_list_helix.append(atom)
                    for atom in res2 :
                        if atom.name=="CA":
                            CA_list_ref_helix.append(atom)
                
            assert len(CA_list_ref_helix)==len(self.CA_list_helix), "ERROR: Fixed and moving atom lists differ in size: "+str(len(CA_list_ref_helix))+' '+str(len(self.CA_list_helix))
                    
            super_imposer=Bio.PDB.Superimposer()
            super_imposer.set_atoms(CA_list_ref_helix, self.CA_list_helix)
            super_imposer.apply(self.get_atoms())

        io=Bio.PDB.PDBIO()
        io.set_structure(self.structure)
        io.save(self.filename[:-4]+"_aligned.pdb")