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
from rdkit.Chem import Draw


import sys

import time
import itertools
import pickle


from AP_convert import *
from copy import deepcopy

global bkgd_cmpd_list; bkgd_cmpd_list=['O','[H]O[H]', 'CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(C(O)C1OP(=O)(O)O)N1:C:N:C2:C(N):N:C:N:C:2:1)C(O)C(=O)NCCC(=O)NCCS', 'O=P(O)(O)O', 'O=C=O', 'N']
from AP_funcs import *
from pprint import pprint

global bad_ss_dict; bad_ss_dict=dict([]) # bad smiles smarts dict

def format_rxn_rules(one_rxn_rule):
    one_rxn=one_rxn_rule.split(">>")
    return one_rxn[0][1:-1]+">>"+one_rxn[1][1:-1]
#####----------parse reaction rules in csv files
import csv
print ("Start parsing csv files and loading reaction rules")
rxn_dict=dict([])
rxn_id_list=[]
rxn_str_list=[]
rxn_ec_num_list=[]
rxn_order_list=[]
rxn_diameter_list=[]
rxn_score_list=[]
#with open('rules/knime-ready-rules_mnx-all-forward_ECOLI-iJO1366.csv', 'r') as csvfile:
with open('rules/APrules.csv', 'r') as csvfile:
    spamreader = csv.reader(csvfile)
    yicixing=0 # for not reading the first line in the csv file
    for row in spamreader:
        if yicixing==0:
            yicixing=1
            continue

        '''
        one_rule_id=row[0]
        #one_rule_id=row[0]+'-'+row[4]
        one_rule_str=row[1]
        one_rule_ec_num=row[2]
        one_rule_rxn_order=row[1].split(">>")[0].count(".")+1
        one_rule_diameter=row[4]
        one_rule_score=row[5]
        one_rxn_tag="forward"
        '''

        one_rule_id=row[0]
        #one_rule_id=row[0]+'-'+row[4]
        one_rule_str=row[1]
        one_rule_ec_num=row[2]
        one_rule_rxn_order=row[1].split(">>")[0].count(".")+1
        one_rule_diameter=1
        one_rule_score=1
        one_rxn_tag="forward"

        rxn_dict[one_rule_id]=[format_rxn_rules(one_rule_str),one_rule_ec_num,one_rule_rxn_order,one_rule_diameter,one_rule_score,one_rxn_tag]
        rxn_id_list.append(one_rule_id)
        rxn_str_list.append(format_rxn_rules(one_rule_str))
        rxn_ec_num_list.append(one_rule_ec_num)
        rxn_order_list.append(one_rule_rxn_order)
        rxn_diameter_list.append(one_rule_diameter)
        rxn_score_list.append(one_rule_score)


j="a"
for i in sorted(rxn_dict):
    print (i)
    print (rxn_dict[i][0].split(">>")[0],rxn_dict[i][0].split(">>")[1])
    ms=[Chem.MolFromSmiles(rxn_dict[i][0].split(">>")[0]),Chem.MolFromSmiles(rxn_dict[i][0].split(">>")[1])]
    #ms=[Chem.MolFromSmarts(rxn_dict[i][0].split(">>")[0]),Chem.MolFromSmarts(rxn_dict[i][0].split(">>")[1])]
    img=Draw.MolsToGridImage(ms,molsPerRow=2,subImgSize=(800,800),legends=["caonima" for x in ms])
    #img.save('rules/MNX_fwd/'+i+'.png')
    img.save('rules/APrules/'+i+'.png')
    #img.save('rules/BNICE_mo/'+i+'.png')
    j=i
print (j)

