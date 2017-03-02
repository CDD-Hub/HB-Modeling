#! -*- coding: utf-8 -*-

#KC#CG# packages import
import os
import numpy as np
from modeller import *

#CG# set all Modeller log output levels
log.level(output=0, notes=0, warnings=0, errors=0, memory=0)

########################################################################
########################                        ########################
########################        functions       ########################
########################                        ########################
######################################################################## --> CG

#KC#CG# function used in fonctions.py
def distance(coord1, coord2):
    """
    calculate distance between 2 points
    """
    
    return np.sqrt((coord1[0]-coord2[0])**2+
                (coord1[1]-coord2[1])**2+
                (coord1[2]-coord2[2])**2)

#KC#CG# function used in model.py
def calcul_angle(point1, point2, point3):
    """
    calculate angle between 3 points
    """
    
    x1,y1,z1=point1
    x2,y2,z2=point2
    x3,y3,z3=point3
        
    vec1=[x1-x2, y1-y2, z1-z2]
    vec2=[x3-x2, y3-y2, z3-z2]

    return calcul_angle_vector(vec1, vec2)

#KC#CG# function used in model.py and geometrie.py
def calcul_angle_vector(vec1, vec2):
    """
    calculate angle between 2 vectors
    """
    
    try:
        div=(vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])/(distance(vec1,[0,0,0])*distance(vec2,[0,0,0]))
        if div>1:
            div=1
        if div<-1:
            div=-1
        #KC# tranlation to degrees
        angle=180/np.pi*np.arccos(div)
    except:
        print vec1
        print vec2
        print (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])/(distance(vec1,[0,0,0])*distance(vec2,[0,0,0]))
    return angle

#KC#CG# function used in main.py
def translation(vecteur_coord, valeur, axe):

    size=np.sqrt(axe[0]*axe[0]+axe[1]*axe[1]+axe[2]*axe[2])
    u_vecteur=axe/size

    vecteur_coord[0]=vecteur_coord[0]-valeur*u_vecteur[0]
    vecteur_coord[1]=vecteur_coord[1]-valeur*u_vecteur[1]
    vecteur_coord[2]=vecteur_coord[2]-valeur*u_vecteur[2]
    return vecteur_coord
    
#KC#CG# function used in main.py
def rotation(vecteur_coord, axe, angle, centre=None):
    #CG# rotation around an axis --> http://fr.wikipedia.org/wiki/Matrice_de_rotation
    
    if centre is not None:
        vecteur_coord[0]=vecteur_coord[0]-centre[0]
        vecteur_coord[1]=vecteur_coord[1]-centre[1]
        vecteur_coord[2]=vecteur_coord[2]-centre[2]

    #CG# vector transformed into unit vector
    size=np.sqrt(axe[0]*axe[0]+axe[1]*axe[1]+axe[2]*axe[2])
    u_vecteur=axe/size
    ux=u_vecteur[0]
    uy=u_vecteur[1]
    uz=u_vecteur[2]
    
    #CG# calculs
    #CG# [R = P +(I-P)cos(angle) + Q*sin(angle)]
    P=np.matrix([[ux*ux, ux*uy, ux*uz],
                 [ux*uy, uy*uy, uy*uz],
                 [ux*uz, uy*uz, uz*uz]])
    I=np.matrix([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]])
    Q=np.matrix([[ 0, -uz,  uy],
                 [uz,   0, -ux],
                 [-uy, ux,  0]])  
    R=P+(I-P)*np.cos(angle)+Q*np.sin(angle)
    vecteur_coord_modif=np.array(np.dot(R,vecteur_coord))[0]
    
    if centre is not None:
        vecteur_coord_modif[0]=vecteur_coord_modif[0]+centre[0]
        vecteur_coord_modif[1]=vecteur_coord_modif[1]+centre[1]
        vecteur_coord_modif[2]=vecteur_coord_modif[2]+centre[2]

    return vecteur_coord_modif

#KC#CG# function used in model.py and main.py
def determine_axe(selection, direction, _type="normal"):

    assert _type in ["normal","perp"], "ERROR"
    #CG# SVD (singular value decomposition) on the mean-centered data
    data=np.array([atom.get_coord() for atom in selection])
    datamean=data.mean(axis=0)
    uu,dd,vv=np.linalg.svd(data-datamean)
    #KC#CG# vv contains axis of SVD
    axe1=vv[0]
    
    #KC#CG# no direction given
    if direction==None:
        return axe1
    
    #KC#CG# direction given
    x2=datamean[0]-direction[0]
    y2=datamean[1]-direction[1]
    z2=datamean[2]-direction[2]
    axe2=[x2,y2,z2]
        
    #KC#CG# perp type
    if _type=="perp":
        return axe2

    #KC#CG# normal type
    #CG# perpendicular axis determination
    x1=axe1[0]
    y1=axe1[1]
    z1=axe1[2]
    x3=1
    #CG# development after calculating scalar products
    y3=(z1/z2*x2*x3-x1*x3)/(y1-z1/z2*y2)
    z3=-(x2*x3+y2*y3)/z2

    axe3=[x3,y3,z3]
    return axe3

#KC#CG# function used in main.py
def clean_tempfile():
    
    os.system('rm *.ini *.sch *.rsr *.D000* *.V999* *.asa *.log *.rsa *.rsa+ *.pyc *.sol blank.pdb res_final.ali blank_model/*.sol > /dev/null 2>&1')
            
#KC#CG# function used in main.py
def get_sequences(filename, pdb_template, sequence_name):
        
    filein=open(filename,'r')
    i=0
    lines=filein.readlines()
    while ">" not in lines[i]:
        i+=1
    name1=lines[i].split(';')[1].strip()
    i+=2
    seq1=""
    while "*" not in seq1:
        seq1+=lines[i].strip()
        i+=1
            
    while ">" not in lines[i]:
        i+=1
    name2=lines[i].split(';')[1].strip()
    i+=2
    seq2=""
    while "*" not in seq2:
        seq2+=lines[i].strip()
        i+=1
    filein.close()
        
    assert len(seq1)==len(seq2), "Sequence length differe between sequence and template"

    if (name1,name2)==(pdb_template,sequence_name):
        return seq1, seq2
    elif (name1,name2)==(sequence_name,pdb_template):
        return seq2,seq1
    else:
        print "ERROR: Sequences cannot be found"
        print name1,name2
        print pdb_template,sequence_name
        return '',''
        
def get_aligned_residues(aln_file):
    
    num_structure=[]
    num_query=[]
    
    env=environ()
    aln=alignment(env, file=aln_file)

    #KC# unaligned residues not selected
    for res_str in aln[0].residues:
        if res_str.get_aligned_residue(aln[1]) != None:
            num_structure.append(int(res_str.get_aligned_residue(aln[0]).num))
    for res_que in aln[1].residues:
        if res_que.get_aligned_residue(aln[0]) != None:
            num_query.append(int(res_que.get_aligned_residue(aln[1]).index))
            
    return num_structure, num_query

def save_matrix(matrix, matrix_name, list_acceptor, list_donor):
    
    with open(matrix_name+".csv", "w") as fileout:       
        fileout.write("\t"+"\t".join([acceptor for acceptor in list_acceptor])+"\n")
        for cpt,line in enumerate(matrix):
            fileout.write(str(list_donor[cpt])+"\t"+"\t".join([str(val) for val in line])+"\n")
    fileout.close()
    
#KC#CG# function not used
def clear_restraints(self):
    self.restraints=[]

#KC#CG# function not used
def add_restraint(self, restraint):
    if isinstance(restraint, list):
        if restraint not in self.restraints:
            self.restraints.append(restraint)
        else:
            'WARNING: add_restraints argument already in restraints list\n         Use get_restraints to see restraints list'
    else:
            print 'ERROR: add_restraints argument have to be a list'

#KC#CG# function not used
def del_restraint(self, restraint):
    if isinstance(restraint, list):
        if restraint in self.restraints:
            self.restraints.remove(restraint)
        else:
            'WARNING: del_restraints argument not in restraints list\n         Use get_restraints to see restraints list'
    else:
            print 'ERROR: del_restraints argument have to be a list'