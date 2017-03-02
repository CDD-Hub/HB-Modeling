##### model.py #####

liste_lig_name=["UNK","BLK"]

donors={     "TYR": {"OH":1},
             "TRP": {"NE1":1}, 
             "LYS": {"NZ":3}, 
             "ARG": {"NE":1, "NH1":2, "NH2":2},
             "HIS": {"ND1":1, "NE2":1},
             "SER": {"OG":1},
             "THR": {"OG1":1},
             "ASN": {"ND2":2},
             "GLN": {"NE2":2}}

             #CG# A VOIR EN FONCTION DE LA CHARGE
             #CG# "ASP": {"OD1":"CG","OD2":"CG"},
             #CG# "GLU": {"OE1":"CD","OE2":"CD"},}
 
acceptors={  "TYR": {"OH":1},
             "HIS": {"NE2":1, "ND1":1},
             "ASP": {"OD1":2},
             "GLU": {"OE1":2,"OE2":2},
             "SER": {"OG":2},
             "THR": {"OG1":2},
             "ASN": {"OD1":2},
             "GLN": {"OE1":2}}

            #CG# A VOIR EN FONCTION DE LA CHARGE
            #CG# "LYS": {"NZ"},  
            #CG# "ARG": {"NE", "NH1", "NH2"}, 

d_min_SS=1
d_max_SS=3

cutoffs_neighbor=[3,3.25,3.5,3.75,4,4.25,4.5,4.75]

##### interactions.py #####

dico_cycles={'PHE':{"cycle":["CG","CD1","CD2","CE1","CE2","CZ"], "atom_A":"CG", "atom_B":"CZ", "atom_C":"CD1"},
             'TYR':{"cycle":["CG","CD1","CD2","CE1","CE2","CZ"], "atom_A":"CG", "atom_B":"CZ", "atom_C":"CD1"},
             'HIS':{"cycle":["CG","ND1","CE1","NE2","CD2"], "atom_A":"CG", "atom_B":"NE2", "atom_C":"CE1"},
             'TRP':{"cycle":["CG","CD2","CE3","CZ3","CH2","CZ2","CE2","NE1","CD1"], "atom_A":"CE3", "atom_B":"CZ2", "atom_C":"CD1"}}

cutoff_bonds=1.9 #CG# cutoff for the research of bond

dico_groups_H_bonds={}                          
for res in acceptors:
    if res not in dico_groups_H_bonds:
        dico_groups_H_bonds[res]=[]
    for atom_name in acceptors[res]:
        dico_groups_H_bonds[res].append(atom_name)
for res in donors:
    if res not in dico_groups_H_bonds:
        dico_groups_H_bonds[res]=[]
    for atom_name in donors[res]:
        if atom_name not in dico_groups_H_bonds[res]:
            dico_groups_H_bonds[res].append(atom_name)

key_to_remove=[]
for key in dico_groups_H_bonds:
    if len(dico_groups_H_bonds[key])==1:
        key_to_remove.append(key)

for key in key_to_remove:
    del dico_groups_H_bonds[key]

#KC#CG# dico_hydrophobiques={}
