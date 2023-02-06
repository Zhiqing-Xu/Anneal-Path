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
# ==================================================================================== #
def readlines_chem_prop(file, mark="", skip_header=False):
    count_x=0
    cmpd_mnxid_smiles_dict=dict([])
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
            cmpd_mnxid_smiles_dict[cmpd_info_list[0]]=cmpd_info_list[6]
        else:
            return cmpd_mnxid_smiles_dict
    return cmpd_mnxid_smiles_dict
# ==================================================================================== #
def main():

    # Process 10^6 at a time in case there are problematic lines in the data file.
    # Write the compound info into 7 dictionaries, each info of with up to 10^6 compounds
    # It is easy to combine the dictionaries with a few lines of code.
    file_MNX_cmpd_address="./chem_prop.tsv"
    file_MNX_cmpd=open(file_MNX_cmpd_address)
    cmpd_mnxid_smiles_dict1=readlines_chem_prop(file_MNX_cmpd,"",True)
    cmpd_mnxid_smiles_dict2=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict3=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict4=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict5=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict6=readlines_chem_prop(file_MNX_cmpd,"",False)
    cmpd_mnxid_smiles_dict7=readlines_chem_prop(file_MNX_cmpd,"",False)

    pickle_out1=open("./Step01_cmpd_mnxid_smiles_dict1","wb")
    pickle.dump(cmpd_mnxid_smiles_dict1, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Step01_cmpd_mnxid_smiles_dict2","wb")
    pickle.dump(cmpd_mnxid_smiles_dict2, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Step01_cmpd_mnxid_smiles_dict3","wb")
    pickle.dump(cmpd_mnxid_smiles_dict3, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Step01_cmpd_mnxid_smiles_dict4","wb")
    pickle.dump(cmpd_mnxid_smiles_dict4, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Step01_cmpd_mnxid_smiles_dict5","wb")
    pickle.dump(cmpd_mnxid_smiles_dict5, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Step01_cmpd_mnxid_smiles_dict6","wb")
    pickle.dump(cmpd_mnxid_smiles_dict6, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./Step01_cmpd_mnxid_smiles_dict7","wb")
    pickle.dump(cmpd_mnxid_smiles_dict7, pickle_out1)
    pickle_out1.close()


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':
    main()




