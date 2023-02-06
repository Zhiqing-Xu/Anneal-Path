# ==================================================================================== #
#!/usr/bin/python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit import DataStructs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem
# ==================================================================================== #
# header
import os, os.path
from sys import platform
if os.name == 'nt' or platform == 'win32':
    os.chdir(os.path.dirname(__file__))
import copy
import pickle
import ast
import numpy as np
import pandas as pd
import scipy.io
import subprocess
from random import shuffle
# ==================================================================================== #
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import pylab as pl
import matplotlib.mlab as mlab
# ==================================================================================== #
#-------------------- (1) --------------------#
def similarity_score(smiles_a, smiles_b, parameter_1="ECFP", parameter_2=2): # Return the similarity of two compounds
    #print smiles_a, smiles_b
    try:
        # parameter_1 is similarity metric selected
        cmpd_a=Chem.MolFromSmiles(str(smiles_a))
        cmpd_b=Chem.MolFromSmiles(str(smiles_b))
        if (parameter_1=="top"):
            fp_a=FingerprintMols.FingerprintMol(cmpd_a)
            fp_b=FingerprintMols.FingerprintMol(cmpd_b)  
            similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
        elif (parameter_1=="MACCS"):
            fp_a=MACCSkeys.GenMACCSKeys(cmpd_a)
            fp_b=MACCSkeys.GenMACCSKeys(cmpd_b)
            similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
        elif (parameter_1=="atom_pairs"):
            fp_a=Pairs.GetAtomPairFingerprint(cmpd_a)
            fp_b=Pairs.GetAtomPairFingerprint(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="vec_pairs"):
            fp_a=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_a)
            fp_b=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="torsions"):
            fp_a=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_a)
            fp_b=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="FCFP"):
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2,useFeatures=True)
            fp_b=AllChem.GetMorganFingerprint(cmpd_b,parameter_2,useFeatures=True)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        else: #ECFP
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2)
            fp_b=AllChem.GetMorganFingerprint(cmpd_b,parameter_2)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)

    except Exception:
        print "ERROR: similarity score"
        print "ERROR: similarity score"
        print "ERROR: similarity score"
        print "ERROR: similarity score"
        print "ERROR: similarity score"
        similarity=0

    return similarity
#-------------------- (2) --------------------#
def maxsimscore(list_a,list_b,fptype,parameter_2=2):
    score_list=[]
    for smiles_a in list_a:
        for smiles_b in list_b:
            score_list.append(similarity_score(smiles_a, smiles_b, fptype, parameter_2=2))
    return max(score_list)
#-------------------- (3) --------------------#
def replace_n(str1, n, str2):
    letters = (
    str2 if i == n else char
        for i, char in enumerate(str1)
    )
    return ''.join(letters)

# ==================================================================================== #
def main():
    #-------------------- (1) --------------------#
    # Open the Step04 temp file and retrieve useful info.
    pickle_in1=open("./Step07_SMILESPair_Distancedict_5000","rb")
    SMILESPair_Distancedict=pickle.load(pickle_in1)
    pickle_in1.close()
    print "Done importing"

    paired_cmpds_list=[] #[   [ { fr{}, fr{} },d ],   [ { fr{}, fr{} },d ],  [{{},{}},d], ....  ]
    for one_pair_SMILES in SMILESPair_Distancedict.keys():
        #print one_pair_SMILES
        paired_cmpds_list.append([{frozenset([one_pair_SMILES[0]]),frozenset([one_pair_SMILES[1]])},SMILESPair_Distancedict[one_pair_SMILES]])
    print len(paired_cmpds_list)
    #-------------------- (2) --------------------#  

    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats

    plt.figure()
    metric_list=["MACCS","atom_pairs","torsions","ECFP"]
    color_list=["red","gold","blue","green"]
    metric_list_2=["MACCS","Atom Pairs","Topological Torsions","ECFP"]
    for i in [0,1,2,3]:
        plt.subplot(2, 2, i+1)
        dis_list=[]
        sim_list=[]
        score_list=[]
        count_x=0
        sim_dist_1=[]
        sim_dist_2=[]
        sim_dist_3=[]
        sim_dist_4=[]
        
        for one_pair_info in paired_cmpds_list: #[   [ { fr{}, fr{} },d ],   [ { fr{}, fr{} },d ],  [{{},{}},d], ....  ]
            count_x+=1
            print count_x
            if count_x>=10002:
                break
            if len(list(list(one_pair_info[0])))==2:
                score=maxsimscore(list(list(one_pair_info[0])[0]),list(list(one_pair_info[0])[1]),metric_list[i])
                distance=one_pair_info[1]
                sim=1/(distance+1.)
                dis_list.append(distance)
                sim_list.append(sim)
                score_list.append(score)
                if distance==1:
                    sim_dist_1.append(score)
                if distance==2:
                    sim_dist_2.append(score)
                if distance==3:
                    sim_dist_3.append(score)
                if distance==4:
                    sim_dist_4.append(score)



        print dis_list # real distance
        print sim_list # similarity correspoding to real distance
        print score_list # similarity scores
        #real_value_list
        #predicted_value
        
        '''
        plt.scatter(sim_list[1:10000],score_list[1:10000],marker='.',s=50,alpha=0.01, color=color_list[i])
        font = {'family' : 'normal', 'weight' : 'bold', 'size': 10}
        plt.rc('font', **font)
        plt.rc('font', size=10)          # controls default text sizes
        plt.rc('axes', titlesize=10)     # fontsize of the axes title
        plt.rc('axes', labelsize=10)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
        plt.rc('legend', fontsize=10)    # legend fontsize
        plt.rc('figure', titlesize=20)  # fontsize of the figure title
        plt.legend([metric_list_2[i]+", R = "+str(round((np.corrcoef(sim_list,score_list)[1,0]),3))+"  "])
        plt.ylim(0,1)
        plt.xlabel('Similarity showing the right distance',**font)
        plt.ylabel('Similarity Score',**font)
        plt.title('Similarity Score vs. Estimated Real Similarity',**font)
        print np.corrcoef(sim_list,score_list)[1,0]
        '''
        
        x1 = np.array(sim_dist_1)
        sns.distplot(x1, hist=False, rug=False, color="red")
        x2 = np.array(sim_dist_2)
        sns.distplot(x2, hist=False, rug=False, color="gold")     
        x3 = np.array(sim_dist_3)
        sns.distplot(x3, hist=False, rug=False, color="green")     
        x4 = np.array(sim_dist_4)
        sns.distplot(x4, hist=False, rug=False, color="blue") 

        font = {'family' : 'normal', 'weight' : 'bold', 'size': 10}
        plt.rc('font', **font)
        plt.rc('font', size=10)          # controls default text sizes
        plt.rc('axes', titlesize=10)     # fontsize of the axes title
        plt.rc('axes', labelsize=10)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
        plt.rc('legend', fontsize=10)    # legend fontsize
        plt.rc('figure', titlesize=20)  # fontsize of the figure title
        plt.legend([metric_list_2[i]+", R = "+str(round((np.corrcoef(sim_list,score_list)[1,0]),3))+"  "])
        plt.xlabel('Estimated Similarity Scores',**font)
        plt.ylabel('PDF',**font)
        plt.title('Distribution Plot',**font)

    plt.show()
    #--------------------------------------------------#
    

   
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':
    main()
