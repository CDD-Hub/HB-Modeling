#! -*- coding: utf-8 -*-

#KC#CG# packages import
import os
import math
import Bio.PDB
from Bio.PDB.Vector import *
from pylab import *
from modeller import *
from modeller.automodel import *

#CG# set all Modeller log output levels
log.level(output=0, notes=0, warnings=0, errors=0, memory=0)

#KC#CG# own files import
from fonctions import *
from interactions import *
import parametres

########################################################################
########################                        ########################
########################        my model        ########################
########################                        ########################
########################################################################

class my_model:
    
    #KC#CG# class constructor
    def __init__(self, filename, seq_model, seq_template, sequence_name, list_helix, DOPE=None, molpdf=None, _type=None):
        
        print "MODEL = "+sequence_name
        
        #KC#CG# parameters tests
        assert _type in [None,"blank","xray"], "_type ERROR"
        
        #KC#CG# variables declaration
        self.DOPE=DOPE
        self.molpdf=molpdf
        self._type=_type
        self.filename=filename
        self.seq_model=seq_model
        self.seq_template=seq_template
        self.sequence_name=sequence_name
        pdb_code=self.filename[:-4]
        pdb_filename = "%s.pdb" % pdb_code
        self.RMSD=None #CG# tab with RMSD on n_models for each atom
        self.asa=None #CG# tab with ASA for each residues for each models
        self.angle_dihedres=None #CG# tab with difference of dihedral angle for each residues for each models
        
        #KC#CG# PDB structure, model, atoms, alpha carbons and residus collected
        self.structure=Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename) 
        self._Model=self.structure[0]
        self.atoms=[atom for atom in self._Model.get_atoms()]
        self.CA=[atom for atom in self._Model.get_atoms() if atom.get_name()=="CA"]
        self.res=[res for res in self._Model.get_residues() if res.get_resname() not in parametres.liste_lig_name]

        #KC#CG# functions launching
        """if (self._type!="xray") and (list_helix!=None):
            self.find_don_acc()
        if (self._type!="xray") and (list_helix==None): #KC# not sure for the second condition
            self.find_secondary_structure()
        else:
            self.res_list_helix=list_helix
        if self._type != "xray":
            self.orient_protein()
            #KC# self.calc_box_size()
            if list_helix==None:
                #KC# self.search_top_down()    
                #KC# self.find_SS_bridges()
                self.find_don_acc()
        if self._type != "xray":
            self.interactions=interactions(self)"""
        
        if (list_helix==None):
            self.find_secondary_structure()
        else:
            self.res_list_helix=list_helix
        if (self._type!="xray"):
            self.orient_protein()
        self.find_don_acc()
        self.interactions=interactions(self)
        
    ##################################################
    ############# get and set functions ##############
    ##################################################    

    def get_CA(self):
        return self.CA
        
    def set_CA(self, value):
        self.CA=value
        
    def get_res(self):
        return self.res
        
    def set_res(self, value):
        self.res=value
        
    def get_atoms(self):
        return self.atoms
        
    def set_atoms(self, value):
        self.atoms=value
        
    def get_filename(self):
        return self.filename
        
    def set_filename(self, value):
        self.filename=value
        
    def get_RMSD(self):
        return self.RMSD
        
    def set_RMSD(self, value):
        self.RMSD=value
        
    def get_asa(self):
        return self.asa
        
    def set_asa(self, value):
        self.asa=value
        
    def get_angles_dihedres(self):
        return self.angles_dihedres
        
    def set_angles_dihedres(self, value):
        self.angles_dihedres=value
    
    def get_box_size(self):
        return self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max
    
    def get_C_ter(self):
        return self.C_ter
           
    def get_N_ter(self):
        return self.N_ter 
         
    def get_list_helix(self):
        return self.res_list_helix
        
    def get_model_interaction(self):
        return self.interactions.get_results()
        
    def get_donors(self):
        return self.donors
    
    def get_acceptors(self):
        return self.acceptors
    
    def get_hbond(self):
        return self.interactions.get_hbond()
    
    def get_hbond_bis(self):
        return self.interactions.get_hbond_bis()
        
    ##################################################
    ################ orient protein ##################
    ################################################## 
           
    def orient_protein(self):
        """ KC - protein right up and centered on 0 """
        
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
            atom.set_coord(rotation(vecteur_coord=atom.get_coord(), axe=np.array(axe3), angle=-2*pi/360*angle, centre=pivot))
            atom.set_coord(translation(vecteur_coord=atom.get_coord(), valeur=distance(pivot, [0,0,0]), axe=np.array(pivot)))
                
        if self.atoms[0].get_coord()[2]<0:
            for atom in self.atoms:
                atom.set_coord(rotation(vecteur_coord=atom.get_coord(), axe=np.array([1,0,0]), angle=pi, centre=[0,0,0]))
        
        return
    
    ##################################################
    ################# calc box size ##################
    ################################################## 
        
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
    ########### find secondary structure #############
    ################################################## 
    
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
    ################ search top down #################
    ################################################## 
    
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
    ############## create bfactor file ###############
    ################################################## 
    
    def create_bfactor_file(self, bfactor, extension):
        """ Kc - file with bfactor created """
        
        #KC#CG# bfactor for all residues
        if len(bfactor)==len(self.get_res()):
            for res,value in zip(self.get_res(), bfactor):
                for atom in res:
                    atom.set_bfactor(value)
        
        #KC#CG# bfactor for all atoms
        elif len(bfactor)==(self.get_atoms()):
            for atom,value in zip(self.get_atoms(), bfactor):
                atom.set_bfactor(value)
        
        #KC#CG# length problem
        else:
            assert False, "ERROR: Wrong bfactor tab size"+str(len(bfactor))+" "+str(len(self.res))+" "+str(len(self.atoms))+""

        #KC#CG# new file saved
        io=Bio.PDB.PDBIO()
        io.set_structure(self.structure)
        io.save(self.filename[:-4]+extension+".pdb")
        print "File "+self.filename[:-4]+extension+".pdb saved"

    ##################################################
    ############## color conservation ################
    ################################################## 
    
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
        io.save(self.sequence_name+"_conservation.pdb")

    ##################################################
    ################### ASA calcul ###################
    ################################################## 
    
    def ASA_calcul(self, sequence_name=None):
        """ KC - model ASA calcul """

        #KC#CG# previous method
        os.system('/home/ckaren/Documents/Codes/Naccess/naccess '+self.filename)
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
        filein.close()
        
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

    ##################################################
    ############## structure alignment ###############
    ################################################## 
    
    def structure_alignement(self, model_ref):
        """ KC - structure alignment between model and reference model """
        
        if model_ref.get_filename()!=self.get_filename():
            self.CA_list_helix=[]
            for i in self.res_list_helix:
                self.CA_list_helix.append(self.CA[i])
                
            CA_list_ref_helix=[]
            for i in self.res_list_helix:
                CA_list_ref_helix.append(model_ref.get_CA()[i])
        
            assert len(CA_list_ref_helix)==len(self.res_list_helix), "ERROR: Fixed and moving atom lists differ in size: "+len(CA_list_ref_helix)+' '+len(self.res_list_helix)
            assert len(self.CA_list_helix)==len(self.res_list_helix), "ERROR: Fixed and moving atom lists differ in size: "+len(self.CA_list_helix)+' '+len(self.res_list_helix)
                    
            super_imposer=Bio.PDB.Superimposer()
            super_imposer.set_atoms(CA_list_ref_helix, self.CA_list_helix)
            super_imposer.apply(self.get_atoms())

        io=Bio.PDB.PDBIO()
        io.set_structure(self.structure)
        io.save(self.filename[:-4]+"_aligned.pdb")

    ##################################################
    ################### RMSD calcul ##################
    ################################################## 
    
    def RMSD_calcul(self, others):
        """ KC - RMSD score calculated """
        
        self.RMSD=[]
        all=others
        all.append(self)
        
        for n in xrange(len(self.CA)):
            var=RMSD_calc([all[i].get_CA()[n] for i in xrange(len(all))])
            self.RMSD.append(var)

        plot([x+1 for x in xrange(len(self.RMSD))], self.RMSD)
        savefig("RMSD.png")
        clf()

        #CG# normalisation
        mini=min(self.RMSD)
        maxi=max(self.RMSD)
        bfactor_RMSD=[(var-mini)*100/(maxi-mini) for var in self.RMSD]

        for model in others:
            model.set_RMSD(self.RMSD)
            model.create_bfactor_file(bfactor_RMSD, "_RMSD")
            
    ##################################################
    ################ dihedral calcul #################
    ################################################## 
        
    def dihedral_calcul(self, others):
        """ KC - dihedral score calculated """
        
        self.angle_dihedres=[]
        all=others
        all.append(self)

        for n in xrange(len(self.res)):
            try:
                ecart_type=stat_ecart_type([calc_dihedral(Vector(all[i].get_res()[n]['C'].get_coord()),
                                                          Vector(all[i].get_res()[n]['CA'].get_coord()),
                                                          Vector(all[i].get_res()[n]['CB'].get_coord()),
                                                          Vector(all[i].get_res()[n]['CG'].get_coord()))/math.pi*180
                                            for i in xrange(len(all))])
            except:
                ecart_type=0.01
            self.angle_dihedres.append(ecart_type)

        plot([x+1 for x in xrange(len(self.angle_dihedres))], self.angle_dihedres)
        savefig("angle_dihedres.png")
        clf()
        
        mini=min(self.angle_dihedres)
        maxi=max(self.angle_dihedres)
        bfactor_angle=[(var-mini)*100/(maxi-mini) for var in self.angle_dihedres]
        assert max(bfactor_angle)<=100, "maximum de bfactor trop haut apres normalisation: "+str(max(bfactor_angle))
        assert min(bfactor_angle)>=0, "minimum de bfactor trop bas apres normalisation: "+str(min(bfactor_angle))

        for model in others:
            model.set_angles_dihedres(self.angle_dihedres)
            model.create_bfactor_file(bfactor_angle, "_dihedre")

    ##################################################
    ############### analyse contacts #################
    ################################################## 
    
    def analyse_contacts(self):
        """ KC - analyse of contacts """
        
        self.find_neighbor()
        self.calc_fraction_neighbor()
            
    ##################################################
    ################ find SS bridges #################
    ################################################## 
    
    def find_SS_bridges(self):
        """ KC - SS bridges search (parameters = between 1 and 3 A) """
        
        self.SS_bridges=[]
        liste_cys=[res for res in self.res if res.get_resname()=="CYS"]
        for cys1 in liste_cys:
            for cys2 in liste_cys:
                if cys1.get_id()!=cys2.get_id():
                    if parametres.d_min_SS<distance(cys1['SG'].get_coord(),cys2['SG'].get_coord())<parametres.d_max_SS:
                        self.SS_bridges.append((cys1.get_id(),cys2.get_id()))   

    ##################################################
    ################## find don acc ##################
    ################################################## 
    
    def find_don_acc(self):
        """ KC - loop on the residue to find donors ans acceptors  """
        
        self.acceptors=[]
        self.donors=[]
        for res in self.res:
            #KC#CG# "N" or "O" atoms added for main chain
            try:
                self.donors.append(res["N"])
            except:
                pass
            try:
                self.acceptors.append(res["O"])
            except:
                pass
            
            #KC#CG# donors or acceptors atoms added for side chain
            if res.resname in parametres.donors:
                for atom_name in parametres.donors[res.resname]:
                    try:
                        self.donors.append(res[atom_name])
                    except:
                        pass
            if res.resname in parametres.acceptors:
                for atom_name in parametres.acceptors[res.resname]:
                    try:
                        self.acceptors.append(res[atom_name])
                    except:
                        pass
        
    ##################################################
    #################### get RMSD ####################
    ################################################## 
    
    def get_RMSD_2(self, blank):
        for atom1, atom2 in zip(self.atoms, blank.atoms):
            assert atom1.name==atom2.name, "ERROR atom type"
            
        super_imposer=Bio.PDB.Superimposer()
        super_imposer.set_atoms([atom for atom in self.atoms if atom.name=="CA"], [atom for atom in blank.atoms if atom.name=="CA"])
        super_imposer.apply([atom for atom in self.atoms if atom.name=="CA"])
            
        return super_imposer.rms
            
    ##################################################
    ################# find neighbor ##################
    ################################################## 

    def find_neighbor(self):
        pass

    ##################################################
    ############ calc fraction neighbor ##############
    ################################################## 

    def calc_fraction_neighbor(self):
        #KC# print "cherche contacts    "
        for res in self.res:
            dico_contact={"CC":0,
                          "OO":0,
                          "NN":0,
                          "NO":0,
                          "ON":0,
                          "CO":0,
                          "OC":0,
                          "CN":0,
                          "NC":0,
                          "SS":0,
                          "SO":0,
                          "OS":0,
                          "SN":0,
                          "NS":0,
                          "SC":0,
                          "CS":0}
                          
            for atom in res:
                for neighbor in self.neighbor[parametres.cutoffs_neighbor[-1]][atom.get_full_id()]:
                    dico_contact[atom.element+neighbor[4][0][0]]+=1

