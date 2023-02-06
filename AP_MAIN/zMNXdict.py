#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Microsoft VS header
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
#if os.name == 'nt' or platform == 'win32':
#    print("Running on Windows")
#    if 'ptvsd' in sys.modules:
#        print("Running in Visual Studio")
#        try:
#            os.chdir(os.path.dirname(__file__))
#            print('CurrentDir: ', os.getcwd())
#        except:
#            pass
##--------------------------------------------------#
#    else:
#        print("Running outside Visual Studio")
#        try:
#            if not 'workbookDir' in globals():
#                workbookDir = os.getcwd()
#                print('workbookDir: ' + workbookDir)
#                os.chdir(workbookDir)
#        except:
#            pass
#--------------------------------------------------#
if os.name != 'nt' and platform != 'win32':
    print("Not Running on Windows")
#--------------------------------------------------#
from rdkit import Chem
from rdkit.Chem import AllChem
#--------------------------------------------------#
# header
import os, os.path
from sys import platform
if os.name == 'nt' or platform == 'win32':
    os.chdir(os.path.dirname(__file__))
import pickle
import ast
import numpy as np
import scipy.io
import subprocess
from random import shuffle
from AP_convert import *
from AP_convert import unique_canonical_smiles_AP as unis
# ==================================================================================== #
def MNX_dict_1():

    pickle_in1=open("./MNXScrapSavings/Test01_smiles_id_nme_dict1","rb")
    smiles_id_nme_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_smiles_id_nme_dict2","rb")
    smiles_id_nme_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_smiles_id_nme_dict3","rb")
    smiles_id_nme_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_smiles_id_nme_dict4","rb")
    smiles_id_nme_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_smiles_id_nme_dict5","rb")
    smiles_id_nme_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_smiles_id_nme_dict6","rb")
    smiles_id_nme_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_smiles_id_nme_dict7","rb")
    smiles_id_nme_dict7=pickle.load(pickle_in1)
    pickle_in1.close()

    smiles_id_nme_dict=smiles_id_nme_dict1.copy()
    smiles_id_nme_dict.update(smiles_id_nme_dict2)
    smiles_id_nme_dict.update(smiles_id_nme_dict3)
    smiles_id_nme_dict.update(smiles_id_nme_dict4)
    smiles_id_nme_dict.update(smiles_id_nme_dict5)
    smiles_id_nme_dict.update(smiles_id_nme_dict6)
    smiles_id_nme_dict.update(smiles_id_nme_dict7)

    # Add meaningless keys to the dict to avoid errors
    smiles_id_nme_dict["BIOMASS"]=["",""]
    smiles_id_nme_dict["MNXM0"]=["",""]
    # ==================================================================================== #
    # Open all pickles and retrieve compound/reaction info.
    pickle_in1=open("./MNXScrapSavings/Test01_nme_smiles_id_dict1","rb")
    nme_smiles_id_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_nme_smiles_id_dict2","rb")
    nme_smiles_id_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_nme_smiles_id_dict3","rb")
    nme_smiles_id_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_nme_smiles_id_dict4","rb")
    nme_smiles_id_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_nme_smiles_id_dict5","rb")
    nme_smiles_id_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_nme_smiles_id_dict6","rb")
    nme_smiles_id_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/Test01_nme_smiles_id_dict7","rb")
    nme_smiles_id_dict7=pickle.load(pickle_in1)
    pickle_in1.close()

    nme_smiles_id_dict=nme_smiles_id_dict1.copy()
    nme_smiles_id_dict.update(nme_smiles_id_dict2)
    nme_smiles_id_dict.update(nme_smiles_id_dict3)
    nme_smiles_id_dict.update(nme_smiles_id_dict4)
    nme_smiles_id_dict.update(nme_smiles_id_dict5)
    nme_smiles_id_dict.update(nme_smiles_id_dict6)
    nme_smiles_id_dict.update(nme_smiles_id_dict7)

    return smiles_id_nme_dict, nme_smiles_id_dict

# ==================================================================================== #
# ==================================================================================== #
def test():
    def read_smiles(dictionary, one_smiles):
        try:    
            print [unis(one_smiles),] + dictionary[ unis( one_smiles ) ]," ,"
        except Exception:
            print ("!!!!! Cannot Find : ", unis(one_smiles), " !!!!!" )
        return 
    # ==================================================================================== #
    # ==================================================================================== #
    smiles_id_nme_dict, nme_smiles_id_dict=MNX_dict_1()
    d1=smiles_id_nme_dict
    # ==================================================================================== #
    # CoA
    #(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)

    '''
    print "#2"
    read_smiles(d1, "OCCCCO")
    read_smiles(d1, "O=CCCCO")
    read_smiles(d1, unis( "C(=O)(CCCO)(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)" ) )
    read_smiles(d1, unis( "O=C(O)CCCO" ) )
    read_smiles(d1, unis( "O=CCCC(=O)O" ) )

    '''
    # ==================================================================================== #
    print ("#3")
    read_smiles(d1, "[H]OC([H])([H])C([H])(N([H])[H])C([H])([H])O[H]")
    read_smiles(d1, "[H]OC([H])([H])C([H])(N([H])[H])C([H])([H])OP(=O)(O[H])O[H]")
    read_smiles(d1, "[H]OC(=O)C([H])([H])C([H])([H])C([H])(C(=O)O[H])N([H])[H]")
    read_smiles(d1, "[H]OC([H])([H])C(=O)C([H])([H])OP(=O)(O[H])O[H]")
    print (nme_smiles_id_dict["2-amino-1,3-propanediol"])
    # ==================================================================================== #
    
    
    # ==================================================================================== #
    print("#8")
    read_smiles(d1, "NCCCCC(N)C(=O)O")
    read_smiles(d1, "N=C(O)CCCCN")   
    read_smiles(d1, "NCCCCC(=O)O")
    read_smiles(d1, "O=CCCCC(=O)O")
    read_smiles(d1, "O=C(O)CCCC(=O)O")
    # ==================================================================================== #
    print("#11")
    read_smiles(d1, "NC(Cc1ccc(O)cc1)C(=O)O")
    read_smiles(d1, "O=C(O)C=Cc1ccc(O)cc1")   
    read_smiles(d1, "C(=O)(C=Cc1ccc(O)cc1)(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
    #read_smiles(d1, "CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(O)=NCCC(O)=NCCSC(=O)C=Cc1ccc(O)cc1")
    read_smiles(d1, "O=C(C=Cc1ccc(O)cc1)c1c(O)cc(O)cc1O")
    read_smiles(d1, "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c21")   
    # ==================================================================================== #
    print("#17")
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O")
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC=C(C)C")   
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC1OC1(C)C")
    read_smiles(d1, "CC(C)=CCCC(C)(O)C1CCC2(C)C1CCC1C3(C)CCC(O)C(C)(C)C3CCC12C")
    read_smiles(d1, "CC(C)=CCCC(C)(O)C1CCC2(C)C1C(O)CC1C3(C)CCC(O)C(C)(C)C3CCC12C") 
    # ==================================================================================== #
    print("#18")
    read_smiles(d1, "Cc1ccc(C)cc1")
    read_smiles(d1, "Cc1ccc(CO)cc1")   
    read_smiles(d1, "Cc1ccc(C=O)cc1")
    read_smiles(d1, "Cc1ccc(C(=O)O)cc1")
    read_smiles(d1, "O=C(O)c1ccc(CO)cc1")   
    read_smiles(d1, "O=Cc1ccc(C(=O)O)cc1")
    read_smiles(d1, "O=C(O)c1ccc(C(=O)O)cc1")

    # ==================================================================================== #
    print("#20")
    read_smiles(d1, "NC(Cc1c[nH]c2ccccc12)C(=O)O")
    read_smiles(d1, "N=C(Cc1c[nH]c2ccccc12)C(=O)O")    
    read_smiles(d1, "N=C(C(=O)O)C(c1c[nH]c2ccccc12)C(C(=N)C(=O)O)c1c[nH]c2ccccc12")
    read_smiles(d1, "O=C(O)c1[nH]c(-c2c[nH]c3ccccc23)cc1-c1c[nH]c2ccccc12")
    read_smiles(d1, "O=C(O)c1[nH]c(-c2c[nH]c3ccc(O)cc23)cc1-c1c[nH]c2ccccc12")    
    read_smiles(d1, "O=C(O)c1[nH]c(-c2c[nH]c3ccc(O)cc23)cc1-c1c(O)[nH]c2ccccc12")
    read_smiles(d1, "O=C1NC(c2c[nH]c3ccc(O)cc23)=CC1=C1C(=O)Nc2ccccc21")

    # ==================================================================================== #
    read_smiles(d1, "NCCN")

    # ==================================================================================== #
    print ("-"*20)
    print (nme_smiles_id_dict["acetyl-coa"])
    print (nme_smiles_id_dict["Malonyl-CoA".lower()])
    print (nme_smiles_id_dict["Succinyl-CoA".lower()])
    print (smiles_id_nme_dict[ unis( "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O" ) ])
    print (nme_smiles_id_dict["2-oxoglutarate".lower()])
# ==================================================================================== #
# ==================================================================================== #
if (__name__ == '__main__'):
    test()









