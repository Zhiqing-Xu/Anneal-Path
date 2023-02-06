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
import scipy.io
import subprocess
from random import shuffle
# ==================================================================================== #
# ==================================================================================== #
def main():
    #-------------------- (1) --------------------#
    # Open all pickles required.
    pickle_in1=open("./Step03_screened_mnxid_rxn_list","rb")
    screened_mnxid_rxn_list=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step03_Used_SMILES_set","rb")
    Used_SMILES_set=pickle.load(pickle_in1)
    pickle_in1.close()
    
    #-------------------- (2) --------------------#
    # Format the data into two sets, nodes_set and edges_set
    count_x=0
    nodes_set=set([])
    edges_list=[]
    for i in screened_mnxid_rxn_list:
        count_x+=1
        # The plot cannot include all 10000 reactions (plot 1/10 to make it not too messy)
        if count_x>2000:
            break
        nodes_set=nodes_set.union(i)
        edges_list.append([set([list(i)[0]]),set([list(i)[1]])])

    nodes_list=list(nodes_set)

    #-------------------- (3) --------------------#
    # Write a json file based on the formated data for making the graph in D3.js
    text_file = open("./Step05_graph.json", "w")
    text_file.write("{")
    text_file.write("\n")
    text_file.write("  \"nodes\": [")
    text_file.write("\n")

    for i in nodes_list:
        if nodes_list.index(i) == len(nodes_list)-1:
            text_file.write("    {\"id\": \""+i+"\", \"group\": "+ "0" +"}" )
        else:
            text_file.write("    {\"id\": \""+i+"\", \"group\": "+ "0" +"}," )
        text_file.write("\n")
    text_file.write("  ],")
    text_file.write("\n")
    text_file.write("  \"links\": [")
    text_file.write("\n")

    for i in edges_list:
        str1="\"source\": [ \""
        for j in i[1]: # source
            str1=str1+j+"\", \""
        str1=str1[0:-3]+" ]"

        str2="\"target\": [ \""
        for j in i[0]: # target
            str2=str2+j+"\", \""
        str2=str2[0:-3]+" ]"

        str3="\"value\": 0"
        if edges_list.index(i)==len(edges_list)-1:
            text_file.write("    {" + str1 + "," + str2 + "," + str3 + "}")
        else:
            text_file.write("    {" + str1 + "," + str2 + "," + str3 + "},")
        text_file.write("\n")

    text_file.write("  ]")
    text_file.write("\n")
    text_file.write("}")
    text_file.write("\n")
    print "ALL DONE!"

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':
    main()

