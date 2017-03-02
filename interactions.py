#! -*- coding: utf-8 -*-

#KC# packages import
import os
import Bio.PDB
from modeller import *

#KC# own files import
import parametres
from fonctions import *
from geometrie import *

#KC# working directory
workdir=os.getcwd()

########################################################################
########################                        ########################
########################        functions       ########################
########################                        ########################
######################################################################## --> CG

def is_stacking(res1, res2):
    """ KC - return boolean for stacking """

    return cycle(res1).is_stacking(cycle(res2))
    
########################################################################
######################                          ########################
######################        interactions      ########################
######################                          ########################
######################################################################## --> CG

class interactions:
    
    #KC#CG class constructor
    def __init__(self, model, renum_dict=None, _color=False, _graph=True, _test=False):
        
        self.model=model
        self.renum_dict=renum_dict
        self._color=_color
        self._graph=_graph
        self._test=_test
        #self.find_SS_bridges()
        #self.search_stacking()
        self.search_H_bond()
        if self._graph:
            self.accessibility_VS_Hbond()

    ##################################################
    ################### is H bond ####################
    ################################################## --> CG
    
    def is_H_bond(self, atom_A, atom_D):
        """ KC - H bond research between 2 atoms """
    
        """ CG -
        acceptor:atom1
        donors:atom2
        source RDOCK
        X-A - - - D-Y
        is Hbond if:
        2.5 < dist (A-D) < 3.5
        90<angle(YDA)<150
        """
             
        #KC# result file written to compare with modeller function
        fileout=open(self.model.filename[:-4]+"_angles.txt",'a')
        
        #KC#CG# H bond research
        if atom_A.get_parent()!=atom_D.get_parent():
            NS=Bio.PDB.NeighborSearch([atom for atom in self.model.atoms if atom_D!=atom])
            is_hbond=False
            for atom_Y in NS.search(center=atom_D.get_coord(), radius=parametres.cutoff_bonds):
                angle=calcul_angle(atom_Y.get_coord(), atom_D.get_coord(), atom_A.get_coord())
                angle=min(angle, 360-angle)
                if 90 < angle < 150 and 2.5 < atom_D-atom_A < 3.5:
                    fileout.write(str(atom_D.get_parent())+'\t'+str(atom_D)+'\t'+str(atom_A.get_parent())+'\t'+str(atom_A)+'\t'+str(angle)+'\n')
                    is_hbond=True
                    break
            if is_hbond:
                return True
        fileout.close()
        
        return False

    ##################################################
    ################ find SS bridges #################
    ################################################## --> CG
    
    def find_SS_bridges(self):
        """ KC - SS bridges search (parameters = between 1 and 3 A) """
        
        self.SS_bridges=[]
        liste_cys=[res for res in self.model.res if res.get_resname()=="CYS"]
        for cys1 in liste_cys:
            for cys2 in liste_cys:
                if cys1.get_id()!=cys2.get_id():
                    if distance(cys1['SG'].get_coord(),cys2['SG'].get_coord())<parametres.dist_SS:
                        if (cys2.get_id(),cys1.get_id()) not in self.SS_bridges:
                            self.SS_bridges.append((cys1.get_id(),cys2.get_id()))

    ##################################################
    ################# search H bond ##################
    ################################################## --> KC

    def search_H_bond(self):
        """ KC - H bond research in model """

        ##############################################################################
        #KC#CG# previous method with function is_H_bond --> to be compared
        """acceptors=self.model.get_acceptors()
        donors=self.model.get_donors()

        NS=Bio.PDB.NeighborSearch(donors)
        self.HB=[]
        
        #CG# for each acceptors, find all the donors with Hbond
        for atom_A in acceptors:
            for atom_D in NS.search(center=atom_A.get_coord(), radius=3.5, level="A"):
                if self.is_H_bond(atom_A, atom_D):
                    self.HB.append((atom_A, atom_D))"""
        ##############################################################################

        #KC# new method with modeller function
        #KC# environment definition
        env = environ()
        #KC# localisation of PDB file of the model
        env.io.atom_files_directory = [os.getcwd()]
        #KC# waters and heteroatoms excluded
        env.io.hetatm = False
        env.io.water = False
        #KC# get topology library for radii
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        #KC# model class recorded in mdl
        mdl = model(env, file=self.model.filename)
        #KC# function write.data in model class launched
        mdl.write_data(file=self.model.filename[:-4], output='HBONDS')

        #KC# file reading to get hbond in list
        liste_hbonds = os.popen("cat "+self.model.filename[:-4]+".hbnds | tr -s ' ' | tr -s '-' ' ' | grep -v \# | cut -d ' ' -f 2-11", "r").read()
        liste_hbonds = liste_hbonds.split("\n")
        del liste_hbonds[-1]
        self.H_bonds=[]
        self.H_bonds_atoms=[]
        self.ionic_bonds_atoms=[]
        #KC# each line of the file analysed
        for line in liste_hbonds :
            res1=[]
            res2=[]
            atom1=""
            atom2=""
            H_bond=[]
            H_bond_atoms=[]
            ionic_bond_atoms=[]
            line = line.split()
            #KC# if chain written in the file (otherwise problems in columns number)
            if line[1] in ["A", "B", "C", "D"] :
                atom1=line[3]+":"+str(line[0])+":"+str(line[1])
                atom2=line[7]+":"+str(line[4])+":"+str(line[5])
                del line[1]
                del line[4]
            else :
                print "list of Hbond atoms not done, chain identification missing !"
            #KC# model sequence which start at 1
            if self.renum_dict==None :
                residue1 = self.model.res[int(line[0])-1]
                residue2 = self.model.res[int(line[3])-1]
            #KC# other sequence which maybe don't start at 1
            else :
                residue1 = self.model.res[self.renum_dict[int(line[0])]]
                residue2 = self.model.res[self.renum_dict[int(line[3])]]
            #KC# both atoms in Hbond added to the H_bond_atom_list
            try:
                if (line[1] in ("LYS", "ARG", "HIS")) and (line[4] in ("ASP", "GLU")) and (line[6]=="SC") and (line[7]=="SC"):
                    ionic_bond_atoms.append(atom1)
                    ionic_bond_atoms.append(atom2)
                    self.ionic_bonds_atoms.append(ionic_bond_atoms)
                H_bond_atoms.append(atom1)
                H_bond_atoms.append(atom2)
                self.H_bonds_atoms.append(H_bond_atoms)
            except:
                pass
            #KC# both residues in hbond added to the H_bond list
            res1.append(residue1)
            res1.append(line[6])
            res2.append(residue2)
            res2.append(line[7])
            H_bond.append(res1)
            H_bond.append(res2)
            
            #KC# H_bond added in the whole list of H_bonds
            self.H_bonds.append(H_bond)

        #KC# network research
        self.H_bonds_networks=self.search_network(self.H_bonds, _type="Hbond")
        if self._graph:
            self.network_distribution(_type="Hbond")
        
        return

    ##################################################
    #################### get bonds ###################
    ################################################## --> KC

    def get_bonds(self, opt):
        """ KC - return ionic bonds or Hbonds """

        if opt=="H_bonds":
            return self.H_bonds_atoms
        elif opt=="ionic_bonds":
            return self.ionic_bonds_atoms

    ##################################################
    ############### search H bond core ###############
    ################################################## --> KC

    def search_H_bond_core(self, res_core):
        """ KC - H bond in hydrophobic core research in model """

        self.H_bonds_atoms_core=[]
        for (res1, res2) in self.H_bonds_atoms:
            H_bond_atoms_core=[]
            res1 = res1.split(":")
            res2 = res2.split(":")
            if (res1[1] in res_core) and (res2[1] in res_core):
                H_bond_atoms_core.append(res1)
                H_bond_atoms_core.append(res2)
                self.H_bonds_atoms_core.append(H_bond_atoms_core)

        return self.H_bonds_atoms_core
        
    ##################################################
    ################ search stacking #################
    ################################################## --> CG
        
    def search_stacking(self):
        """ KC - stacking research in model """
        
        #KC#CG# list of aromatics residues atoms
        atom_for_search=[]
        for res in self.model.res:
            if res.resname in ['PHE','TYR','TRP','HIS']:
                for name in parametres.dico_cycles[res.resname]["cycle"]:
                    try:
                        atom_for_search.append(res[name])
                    except:
                        pass

        #KC#CG# list of residues pairs that have atoms pairs within 6A radius
        self.aromatic=Bio.PDB.NeighborSearch(atom_for_search).search_all(6, level='R')

        #KC#CG# search of stacking
        self.stacking=[]
        for res1,res2 in self.aromatic:
            if is_stacking(res1,res2):
                if res1.get_id()[1]<res2.get_id()[1]:
                    self.stacking.append((res1,res2))
                else:
                    self.stacking.append((res2,res1))

        #KC# network research
        self.stacking_networks=self.search_network(self.stacking, _type="stacking")
        if self._graph:
            self.network_distribution(_type="stacking")
        
        return

    ##################################################
    ################ search network ##################
    ################################################## --> CG/KC
    
    def search_network(self, list_couple, _type):
        """ KC - network research in model """

        list_unique=[]
        for (res1,res2) in list_couple:
            if res1 not in list_unique:
                list_unique.append(res1)
            if res2 not in list_unique:
                list_unique.append(res2)
                
        list_unique_core=[]
        for res in list_unique:
            for res_core in self.model.res_core:
                if res_core==str(res[0].get_full_id()[3][1]):
                    list_unique_core.append(res)
                    
        list_following_res=[]
        for residue in list_unique_core :
            is_in_network=False
            for element in list_following_res:
                if residue in element:
                    is_in_network=True
            #CG# not already in one network
            if not is_in_network:
                following_res=[]
                for res in list_unique_core:
                    previous_res=[]
                    next_res=[]
                    if (_type=="Hbond") and (res[1]=="MC"):
                        if res[0].get_full_id()[3][1]==(residue[0].get_full_id()[3][1]-1):
                            previous_res.append(res)
                            previous_res.append('MC')
                        if res[0].get_full_id()[3][1]==(residue[0].get_full_id()[3][1]+1):
                            next_res.append(res)
                            next_res.append('MC')
                    elif (_type=="stacking"):
                        if res.get_full_id()[3][1]==(residue.get_full_id()[3][1]-1):
                            previous_res.append(res)
                        if res.get_full_id()[3][1]==(residue.get_full_id()[3][1]+1):
                            next_res.append(res)
                    try :
                        if previous_res in list_unique_core:
                            following_res.append(previous_res)
                            following_res.append(residue)
                            list_following_res.append(following_res)
                            following_res=[]
                    except :
                        pass
                    try :
                        if next_res in list_unique_core:
                            following_res.append(residue)
                            following_res.append(next_res)
                            list_following_res.append(following_res)
                    except :
                        pass

        list_couple = list_couple + list_following_res

        list_unique_core=[]
        for (res1,res2) in list_couple:
            if res1 not in list_unique_core:
                list_unique_core.append(res1)
            if res2 not in list_unique_core:
                list_unique_core.append(res2)
                
        networks=[]
        for residue in list_unique_core:
            #CG# already in one network
            is_in_network=False
            for net in networks:
                if residue in net:
                    is_in_network=True
            #CG# not already in one network
            if not is_in_network:
                network=[residue]
                CONTINUE=True
                #CG# continue if we add new members
                while CONTINUE:
                    CONTINUE=False
                    #CG# for each couple, if one element is in the network, add the other one
                    for (res1,res2) in list_couple:
                        if res1 in network and res2 not in network:
                            network.append(res2)
                            CONTINUE=True
                        elif res2 in network and res1 not in network:
                            network.append(res1)
                            CONTINUE=True
                networks.append(network)                

        if not self._color:
            self.select_to_color_network(networks, _type)
        
        return networks
        
    ##################################################
    ############## network distribution ##############
    ################################################## --> KC

    def network_distribution(self, _type):
        """ KC - plot of H bond network distribution realized by R script """

        #KC# initialization
        net_vector=[]
        net_vector_norm=[]
        vector1=""
        vector2=""
        
        #KC# number of D/A for normalisation
        nbAD = len(self.model.get_acceptors()) + len(self.model.get_donors())        
        
        #KC# length of each network collected
        if _type == "Hbond":
            for net in self.H_bonds_networks:
                net_vector.append(len(net))
                net_vector_norm.append(len(net)*1000/nbAD)
        elif _type == "stacking":
            for net in self.stacking_networks:
                net_vector.append(len(net))
                net_vector_norm.append(len(net)*1000/nbAD)
            
        #KC# list sorted in descending order
        net_vector.sort(reverse=True)
        net_vector_norm.sort(reverse=True)
        
        #KC# vector written as one string for easy access by R
        for nb in net_vector:
            vector1=vector1+str(nb)+"_"
        for nb in net_vector_norm:
            vector2=vector2+str(nb)+"_"  
                
        #KC# graph(s) realisation
        if self._test:
            fileout1=open("network_distribution_"+_type+"_test.txt",'a')
            fileout2=open("network_distribution_"+_type+"_test_norm.txt",'a')
            fileout1.write(self.model.filename.replace('_','-')+"_"+str(vector1)+'\n')
            fileout2.write(self.model.filename.replace('_','-')+"_"+str(vector2)+'\n')
            fileout1.close()
            fileout2.close()
        
        #KC# script R launching
        path1 = str(self.model.filename[:-4]) + "_" + _type + "_network.png"
        path2 = str(self.model.filename[:-4]) + "_" + _type + "_network_norm.png" 
        os.system(workdir+"/R_scripts/network_distribution.R " + vector1 + " " + path1 + " " + _type)
        os.system(workdir+"/R_scripts/network_distribution.R " + vector2 + " " + path2 + " " + _type)
        
        return
        
    ##################################################
    ############# select to color network ############
    ################################################## --> KC

    def select_to_color_network(self, networks, _type):
        """ KC - color residue of the network """

        for typ in parametres.list_network_color:
            
            if typ=="all":
                self.color_network(networks, _type, typ)

            net_to_color=[]
            if typ=="first":
                max_length=max(len(net) for net in networks)
                for network in networks:
                    if len(network)==max_length :
                        net_to_color.append(network)
                self.color_network(net_to_color, _type, typ)    
            if typ=="last":
                min_length=min(len(net) for net in networks)
                for network in networks:
                    if len(network)==min_length :
                        net_to_color.append(network)
                self.color_network(net_to_color, _type, typ)
        
        return
        
    ##################################################
    ################## color network #################
    ################################################## --> KC

    def color_network(self, net_to_color, _type, typ):
        """ KC - color residue of the network """

        for rang,network in enumerate(net_to_color):
            for residue in self.model.res :
                if residue in network :
                    for atom in residue :
                        atom.set_bfactor(0)
                else :
                    for atom in residue :
                        atom.set_bfactor(100)
            io=Bio.PDB.PDBIO()
            io.set_structure(self.model.structure)
            io.save(self.model.filename[:-4] + "_" + _type + "_network_" + typ + str(rang+1) + "_" + str(len(network)) + "res.pdb")
        
        return
        
    ##################################################
    ############# accessibility VS Hbond #############
    ################################################## --> KC

    def accessibility_VS_Hbond(self):
        """ KC - analysis of chain (main or polar side) exposure VS number of Hbond """
        
        for typ in ["donor", "acceptor"]:
            Hbond_count, access_MC, access_PSC = [], [], []
            if typ == "donor":
                #KC# donor Hbonds counted for each donor and residue accessibility added
                for donor in self.model.donors :
                    Hbond_nb=0
                    for Hbond in self.H_bonds_atoms:
                        if (donor.get_full_id()[4][0])+":"+str(donor.get_full_id()[3][1])+":"+str(donor.get_full_id()[2])==Hbond[0]:
                            Hbond_nb+=1
                    Hbond_count.append(Hbond_nb)
                    for line in self.model.liste_access:
                        line=line.split()
                        if line[0]==str(donor.get_full_id()[3][1]):
                            access_MC.append(line[3])
                            access_PSC.append(line[2])
            elif typ == "acceptor":
                #KC# donor Hbonds counted for each donor and residue accessibility added
                for acceptor in self.model.acceptors :
                    Hbond_nb=0
                    for Hbond in self.H_bonds_atoms:
                        if (acceptor.get_full_id()[4][0])+":"+str(acceptor.get_full_id()[3][1])+":"+str(acceptor.get_full_id()[2])==Hbond[1]:
                            Hbond_nb+=1
                    Hbond_count.append(Hbond_nb)
                    for line in self.model.liste_access:
                        line=line.split()
                        if line[0]==str(acceptor.get_full_id()[3][1]):
                            access_MC.append(line[3])
                            access_PSC.append(line[2])

            #KC# vector written as one string for easy access by R
            vec_access_MC, vec_access_PSC, vec_Hbond_count = "", "", ""
            for nb in access_MC:
                vec_access_MC=vec_access_MC+str(nb)+"_"
            for nb in access_PSC:
                vec_access_PSC=vec_access_PSC+str(nb)+"_"
            for nb in Hbond_count:
                vec_Hbond_count=vec_Hbond_count+str(nb)+"_"

            #KC# graph(s) realisation
            if self._test:
                fileout1=open("MC_access_VS_Hbond_"+str(typ)+"_test.txt",'a')
                fileout2=open("PSC_access_VS_Hbond_"+str(typ)+"_test.txt",'a')
                fileout3=open("Hbond_count_"+str(typ)+"_test.txt",'a')
                fileout1.write(str(vec_access_MC)+'\n')
                fileout2.write(str(vec_access_PSC)+'\n')
                fileout3.write(str(vec_Hbond_count)+'\n')
                fileout1.close()
                fileout2.close()
                fileout3.close()

            #KC# script R launching
            os.system(workdir+"/R_scripts/access_VS_Hbond.R " + vec_Hbond_count + " " + vec_access_MC + " " + self.model.filename[:-4]+"_MC_access_VS_Hbond_"+str(typ)+".png MC "+typ+" "+str(self.model.ASA_per))
            os.system(workdir+"/R_scripts/access_VS_Hbond.R " + vec_Hbond_count + " " + vec_access_PSC + " " + self.model.filename[:-4]+"_PSC_access_VS_Hbond_"+str(typ)+".png PSC "+typ+" "+str(self.model.ASA_per))

        return
        
    ##################################################
    ################ Hbonds dist model ###############
    ################################################## --> KC

    def Hbonds_dist_model(self, donors, step):
        """ KC - distribution of Hbonds for each donor """
        
        vector=""
        for donor in donors:
            count_donor=0
            for Hbond in self.H_bonds_atoms:
                if donor==Hbond[0]:
                    count_donor+=1
            vector=vector+str(count_donor)+"_"

        os.system(workdir + "/R_scripts/Hbonds_dist.R " + vector + " Hbonds_dist_model_" + str(step) + ".png")            
        
    ##################################################
    #################### Hbonds dist #################
    ################################################## --> KC

    def Hbonds_dist(self, donors, filename):
        """ KC - distribution of Hbonds for each donor """
        
        vector=""
        for donor in donors:
            count_donor=0
            for Hbond in self.H_bonds_atoms:
                if (donor.get_full_id()[4][0])+":"+str(donor.get_full_id()[3][1])+":"+str(donor.get_full_id()[2])==Hbond[0]:
                    count_donor+=1
            vector=vector+str(count_donor)+"_"
            
        if self._test:
            fileout=open("Hbonds_dist_test.txt",'a')
            fileout.write(str(vector)+'\n')
            fileout.close()

        os.system(workdir + "/R_scripts/Hbonds_dist.R " + vector + " " + filename + "_Hbonds_dist.png")