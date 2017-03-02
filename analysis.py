#! -*- coding: utf-8 -*-

#KC#CG# packages import
import os
import numpy as np
from modeller import *

#CG# set all Modeller log output levels
log.level(output=0, notes=0, warnings=0, errors=0, memory=0)

#KC#CG# own files import
from model import *

########################################################################
######################                            ######################
######################          analysis          ######################
######################                            ######################
######################################################################## --> KC

class analysis:
    """ KC - class for analysis """
    
    #KC#CG# class "constructor"
    def __init__(self, all_models, n_models, blank_model, template_model, query_structure_model, list_donor, list_acceptor, \
        contact_matrix_probas, contact_matrix_best_model, contact_matrix_query, opt, p_min, p_max, p_step, list_tanimoto, \
        sum_probas, sum_Hbonds, sum_restraints, confusion, model, step):
              
        #KC#CG# variables declaration
        self.all_models=all_models
        self.n_models=n_models
        self.blank_model=blank_model
        self.template_model=template_model
        self.query_structure_model=query_structure_model
        self.list_donor=list_donor
        self.list_acceptor=list_acceptor
        self.contact_matrix_probas=contact_matrix_probas
        self.contact_matrix_best_model=contact_matrix_best_model
        self.contact_matrix_query=contact_matrix_query
        self.opt=opt
        self.p_min=p_min
        self.p_max=p_max
        self.p_step=p_step
        self.list_tanimoto=list_tanimoto
        self.sum_probas=sum_probas
        self.sum_Hbonds=sum_Hbonds
        self.sum_restraints=sum_restraints
        self.confusion=confusion
        self.model=model
        
        #KC# functions launching
        if self.confusion:
            if self.query_structure_model!=None:
                self.confusion_matrix(model, step)
        else:
            self.probas_graph()
            self.optimization_analysis()
            self.tanimoto_restraints()

    ##################################################
    ################ models analysis #################
    ################################################## --> KC
          
    def probas_graph(self):
        """ KC - graph of probabilities """

        #KC# vector written as one string for easy access by R
        vector1, vector2 ="", ""
        for x in np.arange(self.p_min, self.p_max+self.p_step, self.p_step):
            count=0
            for donor in self.list_donor:
                for acceptor in self.list_acceptor:
                    #print self.contact_matrix_probas[self.list_donor.index(donor)][self.list_acceptor.index(acceptor)]
                    if self.contact_matrix_probas[self.list_donor.index(donor)][self.list_acceptor.index(acceptor)]==x:
                        count+=1
            vector1=vector1+str(x)+"_"
            vector2=vector2+str(count)+"_"

        #KC# script R launching        
        os.system(workdir + "/R_scripts/probas_dist.R " + vector1 + " " + vector2 + " probas_dist.png")
        
        return
    
    ##################################################
    ################## stop criteria #################
    ################################################## --> KC
          
    def optimization_analysis(self):
        """ KC - analysis of stop criteria """

        vector1, vector2, vector3 = "", "", ""
            
        #KC# vectors written as one string for easy access by R
        for nb in self.sum_probas:
            vector1=vector1+str(nb)+"_"
        for nb in self.sum_Hbonds:
            vector2=vector2+str(nb)+"_"
        for nb in self.sum_restraints:
            vector3=vector3+str(nb)+"_"

        #KC# script R launching        
        os.system(workdir + "/R_scripts/optimization_scores.R " + vector1 + " models_sum_probas.png sum_probas")
        os.system(workdir + "/R_scripts/optimization_scores.R " + vector2 + " models_sum_Hbonds.png sum_Hbonds")
        os.system(workdir + "/R_scripts/optimization_scores.R " + vector3 + " models_sum_restraints.png sum_restraints")

        #KC# average of the n_models done by MODELLER and scores added to the list for each modeling
        self.molpdf_all_models, self.DOPE_all_models, self.GA341_all_models, self.RMSD_template_all_models, self.RMSD_query_structure_all_models, self.RMSD_template_helices_all_models, self.RMSD_query_structure_helices_all_models = [], [], [], [], [], [], []
        vector1, vector2, vector3, vector4, vector5, vector6, vector7 = "", "", "", "", "", "", ""
        for models in self.all_models:
            molpdf_models, DOPE_models, GA341_models, RMSD_template_models, RMSD_query_structure_models, RMSD_template_helices_models, RMSD_query_structure_helices_models = 0, 0, 0, 0, 0, 0, 0
            for model in models:
                molpdf_models+=model.molpdf
                DOPE_models+=model.DOPE
                GA341_models+=model.GA341[0]
                RMSD_template_models+=model.get_RMSD_template(self.template_model)
                RMSD_template_helices_models+=model.get_RMSD_template_helices(self.template_model)
                if self.query_structure_model != None:
                    RMSD_query_structure_models+=model.get_RMSD_query_structure(self.query_structure_model)
                    RMSD_query_structure_helices_models+=model.get_RMSD_query_structure_helices(self.query_structure_model)
            self.molpdf_all_models.append(molpdf_models/self.n_models)
            self.DOPE_all_models.append(DOPE_models/self.n_models)
            self.GA341_all_models.append(GA341_models/self.n_models)
            self.RMSD_template_all_models.append(RMSD_template_models/self.n_models)
            self.RMSD_template_helices_all_models.append(RMSD_template_helices_models/self.n_models)
            self.RMSD_query_structure_all_models.append(RMSD_query_structure_models/self.n_models)
            self.RMSD_query_structure_helices_all_models.append(RMSD_query_structure_helices_models/self.n_models)
            
        #KC# vectors written as one string for easy access by R
        for nb in self.molpdf_all_models:
            vector1=vector1+str(nb)+"_"
        for nb in self.DOPE_all_models:
            vector2=vector2+str(nb)+"_"
        for nb in self.GA341_all_models:
            vector3=vector3+str(nb)+"_"
        for nb in self.RMSD_template_all_models:
            vector4=vector4+str(nb)+"_"
        for nb in self.RMSD_query_structure_all_models:
            vector5=vector5+str(nb)+"_"
        for nb in self.RMSD_template_helices_all_models:
            vector6=vector6+str(nb)+"_"
        for nb in self.RMSD_query_structure_helices_all_models:
            vector7=vector7+str(nb)+"_"

        #KC# script R launching        
        os.system(workdir + "/R_scripts/optimization_scores.R " + vector1 + " models_molpdf_score.png molpdf")
        os.system(workdir + "/R_scripts/optimization_scores.R " + vector2 + " models_DOPE_score.png DOPE")
        os.system(workdir + "/R_scripts/optimization_scores.R " + vector3 + " models_GA341_score.png GA341")
        os.system(workdir + "/R_scripts/optimization_scores.R " + vector4 + " models_RMSD_template_score.png RMSD_template " + vector6)
        os.system(workdir + "/R_scripts/optimization_scores.R " + vector5 + " models_RMSD_query_structure_score.png RMSD_query_structure " + vector7)

        self.save_results()
        
        return
        
    ##################################################
    ################## save results ##################
    ################################################## --> KC
            
    def save_results(self):
        """ KC - saveguard of results """
        
        fileout1=open("models_molpdf_score.txt", "w")
        fileout2=open("models_DOPE_score.txt", "w")
        fileout3=open("models_GA341_score.txt", "w")
        fileout4=open("models_RMSD_template_score.txt", "w")
        fileout5=open("models_RMSD_query_structure_score.txt", "w")
        fileout6=open("models_sum_probas.txt", "w")
        fileout7=open("models_sum_Hbonds.txt", "w")
        fileout8=open("models_sum_restraints.txt", "w")
        fileout9=open("models_RMSD_template_helices_score.txt", "w")
        fileout10=open("models_RMSD_query_structure_helices_score.txt", "w")

        fileout1.write("\t".join([str(molpdf) for molpdf in self.molpdf_all_models])+"\n")
        fileout2.write("\t".join([str(DOPE) for DOPE in self.DOPE_all_models])+"\n")
        fileout3.write("\t".join([str(GA341) for GA341 in self.GA341_all_models])+"\n")
        fileout4.write("\t".join([str(RMSD_temp) for RMSD_temp in self.RMSD_template_all_models])+"\n")
        fileout5.write("\t".join([str(RMSD_query) for RMSD_query in self.RMSD_query_structure_all_models])+"\n")
        fileout6.write("\t".join([str(sum_p) for sum_p in self.sum_probas])+"\n")
        fileout7.write("\t".join([str(sum_H) for sum_H in self.sum_Hbonds])+"\n")
        fileout8.write("\t".join([str(sum_r) for sum_r in self.sum_restraints])+"\n")
        fileout9.write("\t".join([str(RMSD_temp_helices) for RMSD_temp_helices in self.RMSD_template_helices_all_models])+"\n")
        fileout10.write("\t".join([str(RMSD_query_helices) for RMSD_query_helices in self.RMSD_query_structure_helices_all_models])+"\n")
            
        fileout1.close()
        fileout2.close()
        fileout3.close()
        fileout4.close()
        fileout5.close()
        fileout6.close()
        fileout7.close()
        fileout8.close()
        
        return
        
    ##################################################
    ############## tanimoto restraints ###############
    ################################################## --> KC
            
    def tanimoto_restraints(self):
        """ KC - analyse of scoring model """
        
        Tanimoto_matrix=[[0 for i in xrange(len(self.list_tanimoto))] for j in xrange(len(self.list_tanimoto))]
        if self.list_tanimoto == []:
            print "Tanimoto list empty !"
        else :
            for rg1, list1 in enumerate(self.list_tanimoto):
                for rg2, list2 in enumerate(self.list_tanimoto):
                    if rg2 > rg1 :
                        intersecting_set=0
                        for restraint in list1:
                            if restraint in list2:
                                intersecting_set+=1
                        A=float(len(list1))
                        B=float(len(list2))
                        C=float(intersecting_set)
                        if (A+B-C) == 0:
                            Tanimoto_matrix[rg2][rg1]=0
                        else :
                            Tanimoto_matrix[rg2][rg1]=(C/(A+B-C))
            
            #KC# Tanimoto matrix saved
            with open("Tanimoto_matrix.csv", "w") as fileout:       
                for line in Tanimoto_matrix:
                    fileout.write("\t".join([str(val) for val in line])+"\n")
            fileout.close()
        
        return
        
    ##################################################
    ################ confusion matrix ################
    ################################################## --> KC
            
    def confusion_matrix(self, model, step):
        """ KC - confusion matrix creation """
        
        TP, TN, FP, FN = 0, 0, 0, 0
        self.confusion_matrix={}
        for donor in self.list_donor:
            for acceptor in self.list_acceptor:
                if [donor,acceptor] in self.model.interactions.get_bonds(self.opt) and [donor,acceptor] in self.query_structure_model.interactions.get_bonds(self.opt) :
                    TP+=1
                elif [donor,acceptor] not in self.model.interactions.get_bonds(self.opt) and [donor,acceptor] not in self.query_structure_model.interactions.get_bonds(self.opt) :
                    TN+=1
                elif [donor,acceptor] in self.model.interactions.get_bonds(self.opt) and [donor,acceptor] not in self.query_structure_model.interactions.get_bonds(self.opt) :
                    FP+=1
                else:
                    FN+=1
                    
        self.confusion_matrix["True +"]=TP
        self.confusion_matrix["True -"]=TN
        self.confusion_matrix["False +"]=FP
        self.confusion_matrix["False -"]=FN
        
        with open("confusion_matrix_" + str(step) + ".txt", "w") as fileout:       
            fileout.write(str(self.confusion_matrix))
        fileout.close()
        
        return
        