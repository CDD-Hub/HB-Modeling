#! -*- coding: utf-8 -*-

#KC# packages import
import os

############################
##### used in model.py #####
############################

#KC#CG# list of ligands which will not be analysed as residues in the model
list_lig_name=["UNK","BLK"]

#KC# creation of donors and acceptors dictionaries from modeller librairy (modeller_donors_acceptors_lib.txt = /nfs/modeller/modlib/donor_acceptor.lib)
for name in ("DONOR", "ACCEPTOR"):
    liste_keys = os.popen("cat modeller_donors_acceptors_lib.txt | tr -s ' ' | grep "+name+" | cut -d ' ' -f 1 | uniq", "r").read().split()    
    liste = os.popen("cat modeller_donors_acceptors_lib.txt | tr -s ' ' | grep "+name+" | cut -d ' ' -f 1,4", "r").read().split("\n")
    del liste[-1]
    
    liste_tup = []
    for key in liste_keys :
        tup = []
        for line in liste :
            line = line.split()
            if line[0] == key :
                tup.append(line[1])
        liste_tup.append(tup)
    if name == "DONOR" :
        donors_dict = dict(zip(liste_keys, liste_tup))
    else :
        acceptors_dict = dict(zip(liste_keys, liste_tup))

#KC# creation of ionic donors and acceptors dictionaries
ionic_donors_dict={}
for key in donors_dict :
    if key in ("LYS", "ARG", "HIS"):
        if key not in ionic_donors_dict:
            ionic_donors_dict[key]=[]
        for atom in donors_dict[key]:
            ionic_donors_dict[key].append(atom)
ionic_acceptors_dict={}
for key in acceptors_dict :
    if key in ("ASP", "GLU"):
        if key not in ionic_acceptors_dict:
            ionic_acceptors_dict[key]=[]
        for atom in acceptors_dict[key]:
            ionic_acceptors_dict[key].append(atom)

#KC#CG# parameter used for SS bond research
dist_SS=7

#KC# list of network which will be color ("first", "last" or "all" possible)
list_network_color=[]

###################################
##### used in interactions.py #####
###################################

#CG# cutoff for the research of bond
cutoff_bonds=1.9

#KC#CG# creation of cycles dictionary for stacking research
dico_cycles={'PHE':{"cycle":["CG","CD1","CD2","CE1","CE2","CZ"], "atom_A":"CG", "atom_B":"CZ", "atom_C":"CD1"},
             'TYR':{"cycle":["CG","CD1","CD2","CE1","CE2","CZ"], "atom_A":"CG", "atom_B":"CZ", "atom_C":"CD1"},
             'HIS':{"cycle":["CG","ND1","CE1","NE2","CD2"], "atom_A":"CG", "atom_B":"NE2", "atom_C":"CE1"},
             'TRP':{"cycle":["CG","CD2","CE3","CZ3","CH2","CZ2","CE2","NE1","CD1"], "atom_A":"CE3", "atom_B":"CZ2", "atom_C":"CD1"}}