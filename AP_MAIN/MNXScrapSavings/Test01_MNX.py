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
import pickle
import ast
import numpy as np
import scipy.io
import subprocess
from random import shuffle
from chemconvert import *
from chemconvert import unique_canonical_smiles_zx as unis
# ==================================================================================== #
def readlines_chem_prop(file, mark="", skip_header=False):
    count_x=0
    smiles_id_nme_dict=dict([])
    # skip lines not containing information
    # There are 388 lines of header (trivial information about the DB)
    if skip_header==True:
        for i in range(388): 
            line = file.readline()

    # Process 10^6 at a time in case there are problematic lines in the data file.
    # It is easy to combine the dictionaries with a few lines of code.
    for i in range(100000): 
        count_x+=1
        line = file.readline()
        if line != "" and line != "\n" :
            cmpd_info_list=line.split("\t")
            print cmpd_info_list[0], count_x+1

            if cmpd_info_list[1].find("compound")==-1:
                one_smiles=cmpd_info_list[6]
                smiles_processed=max(one_smiles.split("."), key = len)
                #smiles_processed=smiles_processed.replace("([*])", "")
                #smiles_processed=smiles_processed.replace("[*]", "")
                smiles_processed=smiles_processed.replace("[O-]", "O")
                smiles_processed=smiles_processed.replace("[N+]", "N")
                smiles_processed=smiles_processed.replace("[NH+]", "N")
                smiles_processed=smiles_processed.replace("[NH2+]", "N")
                smiles_processed=smiles_processed.replace("[NH3+]", "N")
                smiles_processed=smiles_processed.replace("[N-]", "N")
                #smiles_processed=smiles_processed.replace("*-", "N")
                #smiles_processed=smiles_processed.replace("-*", "N")
                #smiles_processed=smiles_processed.replace("(-*)", "N")
                smiles_processed=unique_canonical_smiles_zx(smiles_processed)

                smiles_id_nme_dict[smiles_processed]=[cmpd_info_list[1], cmpd_info_list[0],]
        else:
            return smiles_id_nme_dict
    return smiles_id_nme_dict

# ==================================================================================== #
def readlines_chem_prop2(file, mark="", skip_header=False):
    count_x=0
    smiles_id_nme_dict=dict([])
    # skip lines not containing information
    # There are 388 lines of header (trivial information about the DB)
    if skip_header==True:
        for i in range(388): 
            line = file.readline()

    # Process 10^6 at a time in case there are problematic lines in the data file.
    # It is easy to combine the dictionaries with a few lines of code.
    for i in range(100000): 
        count_x+=1
        line = file.readline()
        if line != "" and line != "\n" :
            cmpd_info_list=line.split("\t")
            print cmpd_info_list[0], count_x+1

            if cmpd_info_list[1].find("compound")==-1:
                one_smiles=cmpd_info_list[6]
                smiles_processed=max(one_smiles.split("."), key = len)
                #smiles_processed=smiles_processed.replace("([*])", "")
                #smiles_processed=smiles_processed.replace("[*]", "")
                smiles_processed=smiles_processed.replace("[O-]", "O")
                smiles_processed=smiles_processed.replace("[N+]", "N")
                smiles_processed=smiles_processed.replace("[NH+]", "N")
                smiles_processed=smiles_processed.replace("[NH2+]", "N")
                smiles_processed=smiles_processed.replace("[N-]", "N")
                #smiles_processed=smiles_processed.replace("*-", "N")
                #smiles_processed=smiles_processed.replace("-*", "N")
                #smiles_processed=smiles_processed.replace("(-*)", "N")
                smiles_processed=unique_canonical_smiles_zx(smiles_processed)

                smiles_id_nme_dict[cmpd_info_list[1].lower()]=[smiles_processed, cmpd_info_list[0],]
        else:
            return smiles_id_nme_dict
    return smiles_id_nme_dict
# ==================================================================================== #
def main():
    # ==================================================================================== #
    # ==================================================================================== #
    # Process 10^6 at a time in case there are problematic lines in the data file.
    # Write the compound info into 7 dictionaries, each info of with up to 10^6 compounds
    # It is easy to combine the dictionaries with a few lines of code.
    '''
    file_MNX_cmpd_address="./chem_prop.tsv"
    file_MNX_cmpd=open(file_MNX_cmpd_address)
    smiles_id_nme_dict1=readlines_chem_prop(file_MNX_cmpd,"",True)
    smiles_id_nme_dict2=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict3=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict4=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict5=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict6=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict7=readlines_chem_prop(file_MNX_cmpd,"",False)

    pickle_out1=open("./Test01_smiles_id_nme_dict1","wb")
    pickle.dump(smiles_id_nme_dict1, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_smiles_id_nme_dict2","wb")
    pickle.dump(smiles_id_nme_dict2, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_smiles_id_nme_dict3","wb")
    pickle.dump(smiles_id_nme_dict3, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_smiles_id_nme_dict4","wb")
    pickle.dump(smiles_id_nme_dict4, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_smiles_id_nme_dict5","wb")
    pickle.dump(smiles_id_nme_dict5, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_smiles_id_nme_dict6","wb")
    pickle.dump(smiles_id_nme_dict6, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_smiles_id_nme_dict7","wb")
    pickle.dump(smiles_id_nme_dict7, pickle_out1)
    pickle_out1.close()
    '''
    # ==================================================================================== #
    '''
    file_MNX_cmpd_address="./chem_prop.tsv"
    file_MNX_cmpd=open(file_MNX_cmpd_address)
    nme_smiles_id_dict1=readlines_chem_prop2(file_MNX_cmpd,"",True)
    nme_smiles_id_dict2=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict3=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict4=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict5=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict6=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict7=readlines_chem_prop2(file_MNX_cmpd,"",False)

    pickle_out1=open("./Test01_nme_smiles_id_dict1","wb")
    pickle.dump(nme_smiles_id_dict1, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_nme_smiles_id_dict2","wb")
    pickle.dump(nme_smiles_id_dict2, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_nme_smiles_id_dict3","wb")
    pickle.dump(nme_smiles_id_dict3, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_nme_smiles_id_dict4","wb")
    pickle.dump(nme_smiles_id_dict4, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_nme_smiles_id_dict5","wb")
    pickle.dump(nme_smiles_id_dict5, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_nme_smiles_id_dict6","wb")
    pickle.dump(nme_smiles_id_dict6, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Test01_nme_smiles_id_dict7","wb")
    pickle.dump(nme_smiles_id_dict7, pickle_out1)
    pickle_out1.close()    
    '''
    # ==================================================================================== #
    # ==================================================================================== #
    # Open all pickles and retrieve compound/reaction info.
    pickle_in1=open("./Test01_smiles_id_nme_dict1","rb")
    smiles_id_nme_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_smiles_id_nme_dict2","rb")
    smiles_id_nme_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_smiles_id_nme_dict3","rb")
    smiles_id_nme_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_smiles_id_nme_dict4","rb")
    smiles_id_nme_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_smiles_id_nme_dict5","rb")
    smiles_id_nme_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_smiles_id_nme_dict6","rb")
    smiles_id_nme_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_smiles_id_nme_dict7","rb")
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
    pickle_in1=open("./Test01_nme_smiles_id_dict1","rb")
    nme_smiles_id_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_nme_smiles_id_dict2","rb")
    nme_smiles_id_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_nme_smiles_id_dict3","rb")
    nme_smiles_id_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_nme_smiles_id_dict4","rb")
    nme_smiles_id_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_nme_smiles_id_dict5","rb")
    nme_smiles_id_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_nme_smiles_id_dict6","rb")
    nme_smiles_id_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Test01_nme_smiles_id_dict7","rb")
    nme_smiles_id_dict7=pickle.load(pickle_in1)
    pickle_in1.close()

    nme_smiles_id_dict=nme_smiles_id_dict1.copy()
    nme_smiles_id_dict.update(nme_smiles_id_dict2)
    nme_smiles_id_dict.update(nme_smiles_id_dict3)
    nme_smiles_id_dict.update(nme_smiles_id_dict4)
    nme_smiles_id_dict.update(nme_smiles_id_dict5)
    nme_smiles_id_dict.update(nme_smiles_id_dict6)
    nme_smiles_id_dict.update(nme_smiles_id_dict7)


    # ==================================================================================== #
    # ==================================================================================== #
    print "Done Merge!"
    print smiles_id_nme_dict[ unis( "CC=O" ) ]
    print nme_smiles_id_dict["acetyl-coa"]
    print
    print
    print nme_smiles_id_dict["Malonyl-CoA".lower()]
    print nme_smiles_id_dict["Succinyl-CoA".lower()]
    print smiles_id_nme_dict[ unis( "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O" ) ]

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':
    main()




