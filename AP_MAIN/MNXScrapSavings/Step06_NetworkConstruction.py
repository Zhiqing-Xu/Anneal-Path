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
def expand_rxn_tree(one_MNXid, RXN_Network_5_step_dict, RXN_dict, max_len=4 ):
    # A brute-force algorithm is used here (INEFFICIENT). 
    # Computation does not take long for this easy problem
    # Efficient graph search algorithm can be used here if constructing a large network
    print one_MNXid
    RXN_tree_5_step=[[],[],[],[]]
    for i in range(max_len):
        if i==0:
            RXN_tree_5_step[i]=list(set(RXN_dict[one_MNXid]))
        else:
            cmpds_tb_expanded=RXN_tree_5_step[i-1]
            next_lv_cmpds=[]
            duplicates=set([one_MNXid,])
            for j in range(i):
                k=i-j-1
                duplicates=duplicates.union(set(RXN_tree_5_step[k]))
            for one_cmpd in cmpds_tb_expanded:
                next_lv_cmpds=[item for item in set(RXN_dict[one_cmpd]) if item not in duplicates]
            RXN_tree_5_step[i]=next_lv_cmpds
    print RXN_tree_5_step
    RXN_Network_5_step_dict[one_MNXid]=RXN_tree_5_step
    return
    
def main():
    #-------------------- (1) --------------------#
    # Open all pickles and retrieve compound/reaction info.
    pickle_in1=open("./Step03_screened_mnxid_rxn_list","rb")
    screened_mnxid_rxn_list=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step03_Used_MNXid_set","rb")
    Used_MNXid_set=pickle.load(pickle_in1)
    pickle_in1.close()

    #-------------------- (2) --------------------#
    # Format the reactions into a dictionary, 
    # values are a list of MNXid's that are linked to the key with one step reaction
    RXN_dict=dict([])
    for one_MNXid in Used_MNXid_set:
        prod_list=[]
        for one_pair in screened_mnxid_rxn_list:
            if one_MNXid in one_pair:
                temp_pair=list(one_pair)
                temp_pair.remove(one_MNXid)
                prod_list.append(temp_pair[0])
        RXN_dict[one_MNXid]=prod_list

    #-------------------- (3) --------------------#
    RXN_Network_5_step_dict=dict([])
    count_x=0
    for one_MNXid in Used_MNXid_set:
        count_x+=1
        print count_x
        expand_rxn_tree(one_MNXid, RXN_Network_5_step_dict, RXN_dict)


    pickle_out1=open("./Step06_RXN_Network_5_step_dict","wb")
    pickle.dump(RXN_Network_5_step_dict, pickle_out1)
    pickle_out1.close()

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':


    main()

