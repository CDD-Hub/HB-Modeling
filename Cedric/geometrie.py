#! -*- coding: utf-8 -*-

#KC#CG# packages import
import math
import numpy as np
from Bio.PDB.Vector import *

#KC# own files import
import fonctions
import parametres

########################################################################
##########################                     #########################
##########################        cycle        #########################
##########################                     #########################
########################################################################

class cycle:
    
    #KC#CG# class "constructor"
    def __init__(self, residue):
        
        self.residue=residue
        list_points=[residue[parametres.dico_cycles[residue.resname]["atom_A"]].get_coord(),
                     residue[parametres.dico_cycles[residue.resname]["atom_B"]].get_coord(),
                     residue[parametres.dico_cycles[residue.resname]["atom_C"]].get_coord()]
        self.plan=plan(list_points=list_points)

        data=np.array([residue[name_atom].get_coord() for name_atom in parametres.dico_cycles[residue.resname]["cycle"]])
        self.center=data.mean(axis=0)

        return
        
    #######################-CG-#######################
    ################## is stacking ###################
    ##################################################
        
    def is_stacking(self, other):
        """ KC - stacking determination between 2 residues """
        
        distance=float("inf")
                    
        for name_atom1 in parametres.dico_cycles[self.residue.resname]["cycle"]:
            for name_atom2 in parametres.dico_cycles[other.residue.resname]["cycle"]:
                d=fonctions.distance(self.residue[name_atom1].get_coord(), other.residue[name_atom2].get_coord())
                if d<distance:
                    distance=d

        angle=fonctions.calcul_angle_vector(other.plan.normale, self.plan.normale)
        angle=min(angle, 180-angle) 

        if 3<distance<4+angle*2/90:
            return True
        else:
            return False

########################################################################
##########################                    ##########################
##########################        plan        ##########################
##########################                    ##########################
########################################################################

class plan:
    
    #KC#CG# class "constructor"
    def __init__(self, list_points=None):
        
        #KC#CG# parameters tests
        if list_points:
            assert len(list_points)>=3, "ERROR: length of list_points avec to be 3 or more"
            for points in list_points:
                assert len(list_points)==3, "ERROR: length of one point in list_points have to be 3"
        
        #CG# plan and normale calculated
        self.points=list_points
        self._calc_plan()
        
    #######################-CG-#######################
    ################### calc plan ####################
    ##################################################
        
    def _calc_plan(self):
        """ KC - calculate equation of plan """
        
        #CG# plan of equation ax+by+cz+d=0
        if len(self.points)==3:
            self._calc_normale()
            self.a, self.b, self.c=self.normale
            self.d=-self.a*self.points[0][0]-self.b*self.points[0][1]-self.c*self.points[0][2]
        else:
            pass
        return

    #######################-CG-#######################
    ################## calc normale ##################
    ##################################################

    def _calc_normale(self):
        """ KC - calculate normale from plan """
        
        #KC#CG# vectors of plan
        vecAB=[i-j for i,j in zip(self.points[0],self.points[1])]
        vecBC=[i-j for i,j in zip(self.points[1],self.points[2])]
        x1,y1,z1=vecAB
        x2,y2,z2=vecBC
        
        #KC#CG# coordinates of normal line
        x3=1
        y3=(z1/z2*x2*x3-x1*x3)/(y1-z1/z2*y2)
        z3=-(x2*x3+y2*y3)/z2
        self.normale=[x3,y3,z3]
        
        return

    #######################-CG-#######################
    #################### distance ####################
    ##################################################

    def distance(self, point_M):
        """ KC - calculate distance point from plan """
        
        #CG# http://fr.wikipedia.org/wiki/Distance_d%27un_point_%C3%A0_un_plan
        xm,ym,zm=point_M
        d=abs(self.a*xm+self.b*ym+self.c*zm)/math.sqrt(self.a*self.a+self.b*self.b+self.c*self.c)
        
        return d