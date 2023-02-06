# ==================================================================================== #
#!/usr/bin/python
from rdkit import Chem
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

# ==================================================================================== #
def main():
    #-------------------- (1) --------------------#
    # Open the Step01 temp file and retrieve useful info.
    
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict1","rb")
    cmpd_mnxid_smiles_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict2","rb")
    cmpd_mnxid_smiles_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict3","rb")
    cmpd_mnxid_smiles_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict4","rb")
    cmpd_mnxid_smiles_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict5","rb")
    cmpd_mnxid_smiles_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict6","rb")
    cmpd_mnxid_smiles_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict7","rb")
    cmpd_mnxid_smiles_dict7=pickle.load(pickle_in1)
    pickle_in1.close()
    cmpd_mnxid_smiles_dict=cmpd_mnxid_smiles_dict1.copy()
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict2)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict3)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict4)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict5)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict6)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict7)
    
    #-------------------- (2) --------------------#
    # Open the Step06 temp file and retrieve useful info.
    pickle_in1=open("./Step06_RXN_Network_5_step_dict","rb")
    RXN_Network_5_step_dict=pickle.load(pickle_in1)
    pickle_in1.close()

    #-------------------- (3) --------------------#


    #-------------------- (4) --------------------#
    # Again, use brute-force here since it wont take too long for this amount of data
    print "Start working on CmpdPair_Distance_dict", len(RXN_Network_5_step_dict)
    CmpdPair_Distance_dict=dict([])
    count_x=0
    for one_cmpd in RXN_Network_5_step_dict.keys():
        count_x+=1
        print count_x
        if count_x>5000:
            break
        for i in range(4):
            for another_cmpd in RXN_Network_5_step_dict[one_cmpd][i]:
                #cmpd_pairs=(one_cmpd,another_cmpd,)
                cmpd_pairs=(cmpd_mnxid_smiles_dict[one_cmpd],cmpd_mnxid_smiles_dict[another_cmpd],)
                if cmpd_pairs in CmpdPair_Distance_dict.keys():
                    CmpdPair_Distance_dict[cmpd_pairs]=min(CmpdPair_Distance_dict[cmpd_pairs],i+1)
                else:
                    CmpdPair_Distance_dict[cmpd_pairs]=i+1
    '''
    for i in CmpdPair_Distance_dict:
        print i, CmpdPair_Distance_dict[i]
        '''

    print "Done CmpdPair_Distance_dict"

    pickle_out1=open("./Step07_SMILESPair_Distancedict_5000","wb")
    pickle.dump(CmpdPair_Distance_dict, pickle_out1)
    pickle_out1.close()

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':

    main()
