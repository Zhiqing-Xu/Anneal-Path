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
def retrieve_cmpd_id(extended_mnx_id):
# Remove the useless part of the compounds id strings
# Example: MNXM89588@MNXD1 -> MNXM89588
    if extended_mnx_id.find("@")==-1:
        print "Unexpected Case", extended_mnx_id 
        return extended_mnx_id
    return extended_mnx_id.split("@")[0]
# ==================================================================================== #
def rxn_str_to_cmpds(one_str):
# Define this function to retrieve reactants/products as well as 
# stoichiometric coefficients from a reaction string.
    coef_list=[]
    cmpds_list=[]
    if one_str.find(" + ") == -1:
        coef_list=[one_str.split(" ")[0]]
        cmpds_list=[retrieve_cmpd_id(one_str.split(" ")[1])]
    else:
        split_temp_list=one_str.split(" + ")
        for one_cmpd in split_temp_list:
            coef_list.append(one_cmpd.split(" ")[0])
            cmpds_list.append(retrieve_cmpd_id(one_cmpd.split(" ")[1]))
    return [coef_list,cmpds_list]
# ==================================================================================== #
def readlines_reac_prop(file, mark=""):
    count_x=0
    mnxid_rxn_list=[]
    # skip lines not containing information
    # There are 388 lines of header (trivial information about the DB)

    for i in range(386): 
        line = file.readline()

    for i in range(50000): 
        count_x+=1
        print count_x 
        line = file.readline()

        if line != "" and line != "\n" :
            one_rxn_info_list=line.split("\t")
            rxn_id=one_rxn_info_list[0]
            rctt_str=one_rxn_info_list[1].split(" = ")[0]
            prod_str=one_rxn_info_list[1].split(" = ")[1]
            print rxn_str_to_cmpds(rctt_str)+rxn_str_to_cmpds(prod_str)
            mnxid_rxn_list.append(rxn_str_to_cmpds(rctt_str)+rxn_str_to_cmpds(prod_str))

        else:
            return mnxid_rxn_list
    return mnxid_rxn_list
# ==================================================================================== #  
def main():
    # From MetaNetX Database (reac_prop.tsv), 
    # obtain reaction information (stoichiometric coefficient and MNXid) 
    # for ~40000 reactions and store the data into a python list.
    file_MNX_rxn_address="./reac_prop.tsv"
    file_MNX_rxn=open(file_MNX_rxn_address)
    mnxid_rxn_list=readlines_reac_prop(file_MNX_rxn,"")

    pickle_out1=open("./Step02_mnxid_rxn_list","wb")
    pickle.dump(mnxid_rxn_list, pickle_out1)
    pickle_out1.close()


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':
    main()




