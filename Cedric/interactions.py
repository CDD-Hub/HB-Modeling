#! -*- coding: utf-8 -*-

#KC#CG# packages import
import pylab as pl
from Bio.PDB.NeighborSearch import *

#KC#CG# own files import
import parametres
from fonctions import *
from geometrie import *

########################################################################
########################                        ########################
########################        functions       ########################
########################                        ########################
########################################################################

def is_stacking(res1, res2):
    """ KC - return boolean for stacking """

    return cycle(res1).is_stacking(cycle(res2))
    
def is_H_bond_old(atom_A, atom_D):
    """
    CG -
    acceptor:atom1
    donors:atom2
    source RDOCK
    X-A - - - D-Y
    is Hbond if:
    dist (A-D) <3.5
    100<angle(XAD)<180
    60<angle(YDA)<120
    """
    
    if atom_A.get_parent()!=atom_D.get_parent():
        NS=NeighborSearch([atom for atom in atom_A.get_parent() if atom_A!=atom])
        is_hbond=False
        for atom_X in NS.search(center=atom_A.get_coord(), radius=parametres.cutoff_bonds):
            angle=calcul_angle(atom_X.get_coord(), atom_A.get_coord(), atom_D.get_coord())
            angle=min(angle, 360-angle)
            if 100<angle<180:
                is_hbond=True
                break
        
        if is_hbond:
            NS=NeighborSearch([atom for atom in atom_D.get_parent() if atom_D!=atom])
            is_hbond=False
            for atom_Y in NS.search(center=atom_D.get_coord(), radius=parametres.cutoff_bonds):
                angle=calcul_angle(atom_Y.get_coord(), atom_D.get_coord(), atom_A.get_coord())
                angle=min(angle, 360-angle)
                if 60<angle<120:
                    is_hbond=True
                    break
            if is_hbond:
                return True
        
    return False
    
########################################################################
######################                          ########################
######################        interactions      ########################
######################                          ########################
########################################################################

class interactions:
    
    #KC#CG class constructor
    def __init__(self, model):
        self.model=model
        self.search_stacking()
        self.search_H_bond()
        self.search_H_bond_bis()
        #CG# for atom in self.model.atoms:
            #CG# self.search_connectivity(atom, 3)

    ##################################################
    ############# get and set functions ##############
    ##################################################

    def get_results(self):
        return self.results
        
    def get_hbond(self):
        return self.H_bonds
        
    def get_hbond_bis(self):
        return self.H_bonds_bis

    ##################################################
    ###################### run #######################
    ##################################################
    
    def run(self, start=0, end=-1):
        """ KC - run of analyse """
        
        #CG# plot Tx satisfaction O,N
        N=0
        satisf=0
        
        selection_test=[]
        for res in self.model.res[start:end]:
            for a in res:
                selection_test.append(a)
        
        #CG# for atom in self.model.atoms:
        for atom in selection_test:
            if atom.name in ["N","O"]:
				N+=1
				for network in self.H_bonds_networks:
					if atom in network:
							satisf+=1		
                
        #CG# H bond number
        count=0
        for a,b in self.H_bonds:
            if a in selection_test or b in selection_test:
                count+=1
        
        networks_tmp=[]
        for network in self.H_bonds_networks:
            for a in network:
                if a in selection_test:
                    networks_tmp.append(network)
                    break
                
        N_Hbond=count
        
        #CG# H bond length
        
        size_Hbond=sum([len(network) for network in networks_tmp])
        mean_Hbond=float(size_Hbond)/len(networks_tmp)
        
        N_reseaux=len(self.H_bonds_networks)
        
        outfile=open("./res.csv","a")
        outfile.write(self.model.filename+" "+str(N)+" "+str(satisf)+" "+str(float(satisf)/N)+" "+str(N_Hbond)+" "+str(size_Hbond)+" "+str(mean_Hbond)+" "+str(N_reseaux)+"\n")
        self.results={"filename":self.model.filename,"hbond size":size_Hbond, "hbond_bis": self.H_bonds_bis}
        outfile.close()
        return
        
        
        #CG# plot score=f(residue)
        
        liste1=[]
        liste2=[]
        liste3=[] 
        
        list_H_main=[]
        list_H_side=[]
        
        dic={True:1, False:0}
        for res in self.model.res:
            is_stack=False
            is_H_bond=0
            is_H_bond_main=0
            is_H_bond_side=0
                        
            l_H_network=0
            for network in self.stacking_networks:
                if res in network:
                    is_stack=True
                    break
            
            for network in self.H_bonds_networks:
                for atom in network:
                    if atom.get_parent()==res:
                        if len(network)>l_H_network:
                            l_H_network=len(network)
                        is_H_bond+=1
                        
                        if atom.name in ['N','O'] and len(network)>is_H_bond_main:
                            is_H_bond_main=len(network)
                        elif atom.name not in ['N','O'] and len(network)>is_H_bond_side:
                            is_H_bond_side=len(network)
                
            if res.get_id()[1] not in self.model.res_list_helix:
                liste1.append(is_H_bond)
                liste2.append(dic[is_stack])
                liste3.append(l_H_network)
            else:
                liste1.append(0)
                liste2.append(0)
                liste3.append(0)
            
            list_H_main.append(is_H_bond_main)
            list_H_side.append(is_H_bond_side)
            
        print float(sum([1 for i in liste1 if i!=0]))/len(liste1)
        print self.model.filename
        
        pl.clf()
        pl.plot(liste1, label="liaison is Hb")
        pl.plot(liste2, label="is Stack")
        pl.plot([a+b for a,b in zip(liste1,liste2)], label="is stack + is Hb")
        pl.plot(liste3, label="liaison l Hb")
        legend(loc=2)
        pl.savefig(self.model.filename[:-4]+"all.png")
        
        pl.clf()
        pl.plot(list_H_main, label="liaison H main")
        pl.plot(list_H_side, label="liaison H side")
        legend(loc=2)
        pl.savefig(self.model.filename[:-4]+"H_bond.png")

    ##################################################
    ################ search stacking #################
    ##################################################
        
    def search_stacking(self):
        """ KC - stacking research in model """
        
        #KC#CG# list of aromatics residues atoms
        atom_for_search=[]
        for res in self.model.res:
            if res.resname in ['PHE','TYR','TRP','HIS']:
                for name in parametres.dico_cycles[res.resname]["cycle"]:
                    atom_for_search.append(res[name])

        #KC#CG# list of residues pairs that have atoms pairs within 6A radius
        self.aromatic_less_5=NeighborSearch(atom_for_search).search_all(6, level='R')

        #KC#CG# search of stacking
        self.stacking=[]
        for res1,res2 in self.aromatic_less_5:
            if is_stacking(res1,res2):
                if res1.get_id()[1]<res2.get_id()[1]:
                    self.stacking.append((res1,res2))
                else:
                    self.stacking.append((res2,res1))

        self.stacking_networks=self.search_network(self.stacking, _type='stacking')

    ##################################################
    ################ search network ##################
    ##################################################
    
    def search_network(self, list_couple_tmp, _type):
        """ KC - network research in model """
        
        list_couple=[ele for ele in list_couple_tmp]

        list_unique=[]
        for (ele1,ele2) in list_couple:
            if ele1 not in list_unique:
                list_unique.append(ele1)
            if ele2 not in list_unique:
                list_unique.append(ele2)

        #CG# if H_bonds we add the couple by groups
        if _type=="H_bonds":
            for ele1 in list_unique:
                for ele2 in list_unique:
                    if ele1!=ele2:
                        parent1=ele1.get_parent()
                        if parent1 == ele2.get_parent() and parent1.resname in parametres.dico_groups_H_bonds:
                            if ele1.get_name() in parametres.dico_groups_H_bonds[parent1.resname] and ele2.get_name() in parametres.dico_groups_H_bonds[parent1.resname]:
                                list_couple.append((ele1,ele2))
            
        networks=[]
        for ele in list_unique:
            #CG# already in one network
            is_in_network=False
            for network in networks:
                if ele in network:
                    is_in_network=True
            #CG# not already in one network
            if not is_in_network:
                network_temp=[ele]
                CONTINUE=True
                #CG# continue if we add new members
                while CONTINUE:
                    CONTINUE=False
                    #CG# for each couple, if one element is in the network, add the other one
                    for (ele1,ele2) in list_couple:
                        if ele1 in network_temp and ele2 not in network_temp:
                            network_temp.append(ele2)
                            CONTINUE=True
                        elif ele2 in network_temp and ele1 not in network_temp:
                            network_temp.append(ele1)
                            CONTINUE=True
                
                networks.append(network_temp)

        return networks

    ##################################################
    ################# search H bond ##################
    ##################################################

    def search_H_bond(self):
        """ KC - H bond research in model (parameter = 3,5) """
        
        acceptors=self.model.get_acceptors()
        donors=self.model.get_donors()

        NS=NeighborSearch(donors)
        self.H_bonds=[]
        self.H_bonds_networks=[]
        
        #CG# for each acceptors, find all the donors with Hbond
        for atom_A in acceptors:
            for atom_D in NS.search(center=atom_A.get_coord(), radius=3.5, level="A"):
                if self.is_H_bond(atom_A, atom_D):
                    self.H_bonds.append((atom_A, atom_D))
                    
        self.H_bonds_networks=self.search_network(self.H_bonds, _type='H_bonds')
        
        return

    ##################################################
    ############### search H bond bis ################
    ##################################################

    def search_H_bond_bis(self):
        """ KC - H bond research in model (parameter = 6) """
        
        acceptors=self.model.get_acceptors()
        donors=self.model.get_donors()
        NS=NeighborSearch(donors)
        self.H_bonds_bis=[]
        #CG# for each acceptors, find all the donors with Hbond
        for atom_A in acceptors:
            for atom_D in NS.search(center=atom_A.get_coord(), radius=6, level="A"):
                if self.is_H_bond(atom_A, atom_D):
                    self.H_bonds_bis.append((atom_A, atom_D))
                    
        self.H_bonds_networks_bis=self.search_network(self.H_bonds_bis, _type='H_bonds')
        
        return
   
    ##################################################
    ################# search dipoles #################
    ##################################################

    def search_dipoles(self):
        pass
    
    ##################################################
    ############## search connectivity ###############
    ##################################################
        
    def search_connectivity(self, atom_center, rank, NS=None, previous=None):
        
        if NS is None:
            NS=NeighborSearch(self.model.atoms)
        
        resultat=[atom for atom in NS.search(center=atom_center.get_coord(), radius=parametres.cutoff_bonds, level="A") if atom not in [atom_center,previous]]
        
        if rank == 1:
            return resultat
        else:
            connectivity={}
            for atom in resultat:
                resultat_temp=self.search_connectivity(atom_center=atom, rank=rank-1, NS=NS, previous=atom_center)
                if len(resultat_temp)==0:
                    connectivity[atom]=None
                else:
                    connectivity[atom]=resultat_temp
            return connectivity
            
    ##################################################
    ################### is H bond ####################
    ##################################################
    
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
            NS=NeighborSearch([atom for atom in self.model.atoms if atom_D!=atom])
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
