#! -*- coding: utf-8 -*-

#KC#CG# packages import
import os
import math
import shutil
from Bio.PDB.Vector import *
from pylab import *

#KC#CG# own files import
from geometrie import *

########################################################################
########################                        ########################
########################        functions       ########################
########################                        ########################
########################################################################

#KC#CG# function not used
def distance_droite_point(point_A, point_B, point_M):
    """
    calculate distance between point M and a droite passing by the point A and the point B
    H is the projected of M on AB
    H=A+k*(A-B)
    """
    xa,ya,za=point_A
    xb,yb,zb=point_B
    xm,ym,zm=point_M
  
    xab=xa-xb
    yab=ya-yb
    zab=za-zb

    xma=xm-xa
    yma=ym-ya
    zma=zm-za
    
    if (xab*xma + yab*yma + zab*zma)!=0:
        k=(xab*xma + yab*yma + zab*zma)/float((xab*xab + yab*yab + zab*zab))
        
        xh=xa + k*(xa-xb)
        yh=ya + k*(ya-yb)
        zh=za + k*(za-zb)
        
        point_H=[xh, yh, zh]
        
        return distance(point_H, point_M)
    else:
        return distance(point_A, point_M)

#KC#CG# function used in fonctions.py
def stat_moyenne(echantillon):
    """
    calculate average of a sample
    """
    
    taille=len(echantillon)
    moyenne=sum(echantillon)/taille
    return moyenne

#KC#CG# function used in fonctions.py
def stat_variance(echantillon):
    """
    calculate variance of a sample
    """
    
    n=len(echantillon) #CG# size
    mq=stat_moyenne(echantillon)**2
    s=sum([x**2 for x in echantillon])
    variance=s/n-mq
    return variance

#KC#CG# function used in model.py
def stat_ecart_type(echantillon):
    """
    calculate standard deviation of a sample
    """
    
    variance=stat_variance(echantillon)
    ecart_type=math.sqrt(variance)
    return ecart_type

#KC#CG# function used in model.py
def RMSD_calc(tab_atoms, reference=None):
    """
    calculate RMSD score
    """
    
    n_atoms=len(tab_atoms)
    if reference == None:
        barycenter=0.0,0.0,0.0
        for i in xrange(1, n_atoms):
            barycenter+=tab_atoms[i].get_coord()
        for i in [0,1,2]:
            barycenter = list(barycenter)
            barycenter[i]=barycenter[i]/n_atoms
            barycenter = tuple(barycenter)
    else:
        barycenter=reference.get_coord()
    
    total=0
    for i in xrange(n_atoms):
        total+=(distance(tab_atoms[i].get_coord(), barycenter))**2
    
    return sqrt(total/n_atoms)

#KC#CG# function used in fonctions.py
def distance(coord1, coord2):
    """
    calculate distance between 2 points
    """
    
    return sqrt((coord1[0]-coord2[0])**2+
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
        #KC#CG# tranlation to degrees
        angle=180/math.pi*math.acos(div)
    except:
        print vec1
        print vec2
        print (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])/(distance(vec1,[0,0,0])*distance(vec2,[0,0,0]))
    return angle

#########################
#########################
#### To be completed ####

#KC#CG# function not used
def find_site_actif(model):
    
    return [1,2,3,4,5,6,7,8,9]
    
#KC#CG# function not used
def def_donn_acc(model):
    
    return [1,2,3,4,5,6,7,8,9],[1,2,3,4,5,6,7,8,9]
    
#KC#CG# function not used
def def_ions(model):
    
    return [1,2,3,4,5,6,7,8,9]

#KC#CG# function not used
def def_aromatics(model):
    
    return [1,2,3,4,5,6,7,8,9]

#########################
#########################

#KC#CG# function used in main.py
def translation(vecteur_coord, valeur, axe):

    size=sqrt(axe[0]*axe[0]+axe[1]*axe[1]+axe[2]*axe[2])
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
    size=sqrt(axe[0]*axe[0]+axe[1]*axe[1]+axe[2]*axe[2])
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
    R=P+(I-P)*cos(angle)+Q*sin(angle)
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
    os.system('rm *.ini *.sch *.rsr *.D000* *.V999* *.asa *.log *.rsa+ > /dev/null 2>&1')
    
#KC#CG# function used in main.py
def move_files(dir_name):
    try:
        os.mkdir(dir_name)
    except:
        pass
    for filename in os.listdir('.'):
        if ".B99" in filename or '.png' in filename or 'colored' in filename or 'conservation' in filename:
            shutil.move(filename, dir_name)
            
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
