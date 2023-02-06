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
import sys
import time
import pickle
import itertools
from tqdm import tqdm
from pprint import pprint
from copy import deepcopy
#--------------------------------------------------#
from AP_funcs import *
from AP_convert import *
from AP_convert import unique_canonical_smiles_zx as unis
#--------------------------------------------------#
global bad_ss_dict; bad_ss_dict=dict([]) # bad smiles smarts dict
global bkgd_cmpd_list; bkgd_cmpd_list=bkgd_cmpd_list_func()
global CoA_cmpd_list; CoA_cmpd_list=CoA_cmpd_list_func()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
global fwd_CoA_rxn; fwd_CoA_rxn=[( ("Acetyl-CoA",) , (unis("O=C(C(=O)O)C"),) , "1.2.4.a.1.f" ),
                                 ( (unis("O=C(C(=O)O)C"),) , ("Acetyl-CoA",) , "1.2.4.a.1.r" ),
                                 ( ("Acetyl-CoA",) , (unis("CC(=O)O"),) , "3.1.2.a.1.f" ),
                                 ( (unis("CC(=O)O"),) , ("Acetyl-CoA",) , "3.1.2.a.1.r" ),
                                 ( ( "Acetyl-CoA", unis("C(=O)O") ) , ( unis("O=C(C(=O)O)C"),) , "2.3.1.a.1.f" ),
                                 ( ( unis("O=C(C(=O)O)C"),) , ( "Acetyl-CoA", unis("C(=O)O") ) , "2.3.1.a.1.r" ),
                                 ( ("Formyl-CoA",) , (unis("C(=O)"),) , "1.2.1.a.1.r" ),
                                 ( (unis("C(=O)"),) , ("Formyl-CoA",) , "1.2.1.a.1.f" ),

                                 ] # For test only


global bwd_CoA_rxn; bwd_CoA_rxn=[( ("Acetyl-CoA",) , (unis("O=C(C(=O)O)C"),) , "1.2.4.a.1.f" ),
                                 ( (unis("O=C(C(=O)O)C"),) , ("Acetyl-CoA",) , "1.2.4.a.1.r" ),
                                 ( ("Acetyl-CoA",) , (unis("CC(=O)O"),) , "3.1.2.a.1.f" ),
                                 ( (unis("CC(=O)O"),) , ("Acetyl-CoA",) , "3.1.2.a.1.r" ),
                                 ( ( "Acetyl-CoA", unis("C(=O)O") ) , ( unis("O=C(C(=O)O)C"),) , "2.3.1.a.1.f" ),
                                 ( ( unis("O=C(C(=O)O)C"),) , ( "Acetyl-CoA", unis("C(=O)O") ) , "2.3.1.a.1.r" ),
                                 ( ("Formyl-CoA",) , (unis("C(=O)"),) , "1.2.1.a.1.r" ),
                                 ( (unis("C(=O)"),) , ("Formyl-CoA",) , "1.2.1.a.1.f" ),

                                 ] # For test only
                                
fwd_CoA_rxn=[] # For test only
bwd_CoA_rxn=[] # For test only

#######################################################################################################################################
#######################################################################################################################################
global unfavorable_rxn_list; unfavorable_rxn_list=[
                                                 ( (unis("CC(=O)I"),) , (unis("CC(=O)O"),) ),
                                                 ( (unis("CC(=O)O"),) , (unis("CC(=O)I"),) ),

                                                 ] # For test only

def unfavorable_rxn_list_func():
    return unfavorable_rxn_list # For test only



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
class AP_Reactor:
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def __init__(self, rxn_diameter=2, rxn_score_lb=0,rxn_db="AP"):
        def format_rxn_rules(one_rxn_rule):
            one_rxn=one_rxn_rule.split(">>")
            return one_rxn[0][1:-1]+">>"+one_rxn[1][1:-1]
        #####----------parse reaction rules in csv files
        import csv
        #print ("Start parsing csv files and loading reaction rules")
        rxn_dict=dict([])
        rxn_id_list=[]
        rxn_str_list=[]
        rxn_ec_num_list=[]
        rxn_order_list=[]
        rxn_diameter_list=[]
        rxn_score_list=[]
        # ============================================================================================================================ #
        #####---------- (1) RetroPath
        if rxn_db=="retropath":
            #####---------- # 1. parse forward
            with open('../rxn_rule/knime-ready-rules_mnx-all-forward_ECOLI-iJO1366.csv', 'r') as csvfile:
                spamreader = csv.reader(csvfile)
                for row in spamreader:
                    one_rule_id=row[0]+"-"+row[4]+"-"+"fwd"
                    one_rule_str=row[1]
                    one_rule_ec_num=row[2]
                    one_rule_rxn_order=row[1].split(">>")[0].count(".")+1
                    one_rule_diameter=row[4]
                    one_rule_score=row[5]
                    one_rxn_tag="forward"
                    if one_rule_id != "Rule ID-Diameter-fwd":
                        rxn_dict[one_rule_id]=[format_rxn_rules(one_rule_str),one_rule_ec_num,one_rule_rxn_order,one_rule_diameter,one_rule_score,one_rxn_tag]
                        rxn_id_list.append(one_rule_id)
                        rxn_str_list.append(format_rxn_rules(one_rule_str))
                        rxn_ec_num_list.append(one_rule_ec_num)
                        rxn_order_list.append(one_rule_rxn_order)
                        rxn_diameter_list.append(one_rule_diameter)
                        rxn_score_list.append(one_rule_score)
            #####---------- # 2. parse backward
            with open('../rxn_rule/knime-ready-rules_mnx-all-reverse_ECOLI-iJO1366.csv', 'r') as csvfile:
                spamreader = csv.reader(csvfile)
                for row in spamreader:

                    one_rule_id=row[0]+"-"+row[4]+"-"+"bwd"
                    one_rule_str=row[1]
                    one_rule_ec_num=row[2]
                    one_rule_rxn_order=row[1].split(">>")[0].count(".")+1
                    one_rule_diameter=row[4]
                    one_rule_score=row[5]
                    one_rxn_tag="reverse"
                    if one_rule_id != "Rule ID-Diameter-bwd":
                        rxn_dict[one_rule_id]=[format_rxn_rules(one_rule_str),one_rule_ec_num,one_rule_rxn_order,one_rule_diameter,one_rule_score,one_rxn_tag]
                        rxn_id_list.append(one_rule_id)
                        rxn_str_list.append(format_rxn_rules(one_rule_str))
                        rxn_ec_num_list.append(one_rule_ec_num)
                        rxn_order_list.append(one_rule_rxn_order)
                        rxn_diameter_list.append(one_rule_diameter)
                        rxn_score_list.append(one_rule_score)
            print ("Finish loading reaction rules")
            #####----------Select Reaction Rules
            print ("rxn_diameter: ", rxn_diameter)
            print ("rxn_score_lb: ", rxn_score_lb)
            self.selected_rxn_id_list=[]
            self.selected_rxn_dict=dict([])
            for one_rule_id in rxn_id_list:
                if rxn_dict[one_rule_id][3]==str(rxn_diameter) and float(rxn_dict[one_rule_id][4])>=rxn_score_lb:
                    self.selected_rxn_id_list.append(one_rule_id)
                    self.selected_rxn_dict[one_rule_id]=rxn_dict[one_rule_id]
                
        # ============================================================================================================================ #
        #####---------- (2) APrules
        if rxn_db=="APrules":
            with open('../rxn_rule/APrules.csv', 'r') as csvfile:
                spamreader = csv.reader(csvfile)
                yicixing=0 # for not reading the first line in the csv file
                for row in spamreader:
                    if yicixing==0:
                        yicixing=1
                        continue

                    one_rule_id=row[0]
                    one_rule_str=row[1]
                    one_rule_ec_num=row[2]
                    one_rule_rxn_order=row[1].split(">>")[0].count(".")+1
                    one_rule_diameter=row[4]
                    one_rule_score=row[5]
                    one_rxn_tag=row[5]

                    rxn_dict[one_rule_id]=[format_rxn_rules(one_rule_str),one_rule_ec_num,one_rule_rxn_order,one_rule_diameter,one_rule_score,one_rxn_tag]
                    rxn_id_list.append(one_rule_id)
                    rxn_str_list.append(format_rxn_rules(one_rule_str))
                    rxn_ec_num_list.append(one_rule_ec_num)
                    rxn_order_list.append(one_rule_rxn_order)
                    rxn_diameter_list.append(one_rule_diameter)
                    rxn_score_list.append(one_rule_score)
                self.selected_rxn_id_list=rxn_id_list
                self.selected_rxn_dict=rxn_dict
            
        # ============================================================================================================================ #
        if rxn_db != "APrules" and rxn_db != "retropath":
            with open('../rxn_rule/knime-ready-rules_BNICE.csv', 'r') as csvfile:
                spamreader = csv.reader(csvfile)
                yicixing=0 # for not reading the first line in the csv file
                for row in spamreader:
                    if yicixing==0:
                        yicixing=1
                        continue
                    one_rule_id=row[0]
                    one_rule_str=remove_hydrogen_nodes_in_rule(row[1])
                    one_rule_ec_num=row[3]
                    one_rule_rxn_order=row[1].split(">>")[0].count(".")+1
                    one_rule_diameter=-1
                    one_rule_score=row[2]
                    one_rxn_tag="BNICE"

                    rxn_dict[one_rule_id]=[format_rxn_rules(one_rule_str),one_rule_ec_num,one_rule_rxn_order,one_rule_diameter,one_rule_score,one_rxn_tag]
                    rxn_id_list.append(one_rule_id)
                    rxn_str_list.append(format_rxn_rules(one_rule_str))
                    rxn_ec_num_list.append(one_rule_ec_num)
                    rxn_order_list.append(one_rule_rxn_order)
                    rxn_diameter_list.append(one_rule_diameter)
                    rxn_score_list.append(one_rule_score)
                self.selected_rxn_id_list=rxn_id_list
                self.selected_rxn_dict=rxn_dict


#######################################################################################################################################
#######################################################################################################################################
    def pattern_matching_zx(self, cmpd_smiles, rxn_rule_subgraph):
        #####----------return boolean variable (if the compound match the subgraph)
        molecule = MolFromSmiles_zx(cmpd_smiles, bad_ss_dict)
        try:
            pattern_matching = molecule.HasSubstructMatch(Chem.MolFromSmarts(rxn_rule_subgraph))
        except Exception:
            pattern_matching = False
        return pattern_matching

    # ============================================================================================================================ #
    def generate_new_reactions(self, reactants_smiles_list, reaction_id):
        #####----------return a list of all possible products sets, each set is a group of possible products
        r_fmt_rxn=AllChem.ReactionFromSmarts(self.selected_rxn_dict[reaction_id][0])
        r_fmt_reactants_list=[]
        for one_reactant_smiles in reactants_smiles_list:
            r_fmt_reactants_list.append(MolFromSmiles_zx(one_reactant_smiles, bad_ss_dict))
        r_fmt_reactants_tuple=tuple(r_fmt_reactants_list)
        #####----------try/except to avoid rdkit error
        try:
            r_fmt_products_lists = r_fmt_rxn.RunReactants(r_fmt_reactants_tuple)
        except Exception:
            r_fmt_products_lists=[]
        products_lists=[]
        for r_fmt_one_set_products in r_fmt_products_lists:
            one_set_products=[]
            for r_fmt_one_product in r_fmt_one_set_products:
                one_set_products.append(MolToSmiles_zx(r_fmt_one_product, bad_ss_dict))
            products_lists.append(one_set_products)
        return products_lists
    # ============================================================================================================================ #
    def generate_new_reactions_reverse(self, reactants_smiles_list, reaction_id): # probably not need this one for now
        #####----------
        def swap_one_reaction_rule(one_reaction_rule):
            return one_reaction_rule.split(">>")[1]+'>>'+one_reaction_rule.split(">>")[0]
        #####----------return a list of all possible products sets, each set is a group of possible products
        r_fmt_rxn=AllChem.ReactionFromSmarts(swap_one_reaction_rule(self.selected_rxn_dict[reaction_id][0]))
        r_fmt_reactants_list=[]
        for one_reactant_smiles in reactants_smiles_list:
            r_fmt_reactants_list.append(MolFromSmiles_zx(one_reactant_smiles, bad_ss_dict))
        r_fmt_reactants_tuple=tuple(r_fmt_reactants_list)
        r_fmt_products_lists = r_fmt_rxn.RunReactants(r_fmt_reactants_tuple)

        products_lists=[]
        for r_fmt_one_set_products in r_fmt_products_lists:
            one_set_products=[]
            for r_fmt_one_product in r_fmt_one_set_products:
                one_set_products.append(MolToSmiles_zx(r_fmt_one_product, bad_ss_dict))
            products_lists.append(one_set_products)
        return products_lists

#######################################################################################################################################
#######################################################################################################################################
    def apply_enzymes(self,fwd_smiles_set,fwd_smiles_set_level_list,subs_level,max_C_num,max_O_num):
    # Return only all the reactions found, [new_reactions_list], would be enough, the lists used will be updated in zxpathfinder.py
        print ("Reaction Simulation ... ")
        #============================================================================================================================#
        #print fwd_smiles_set
        #print fwd_smiles_set_level_list[subs_level]
        #####----------Step 0 Process the inputs:

        #============================================================================================================================#
        #####----------Step 1 Generate new reactions and store in the new_reaction_list (in one loop):
        new_reactions_list=[]
        for one_rule_id in tqdm(self.selected_rxn_id_list):
            #print "############################################################"
            #print 'Sub-Side-Depth', subs_level+1, '-',
            #print one_rule_id
                
            #begin_time = time.time()
            # There can be one reactant or more than one reactant, including: (1)reaction coef>1  (2)two reactants  (3)one reactant contains all subs.)
            #============================================================================================================================#
            #####----------Step 1-1 Obtain enzyme data: 
            one_rule_str=self.selected_rxn_dict[one_rule_id][0]
            one_rule_rxn_order=self.selected_rxn_dict[one_rule_id][2]
            one_rule_tag=self.selected_rxn_dict[one_rule_id][5]
            if one_rule_tag =="reverse": # skip if not forward (not applicable to APrules)
                continue
            #============================================================================================================================#
            #####----------Step 1-2 Type One Enzyme, one reactant: 
            if one_rule_rxn_order==1 :
                #print 'one reactant'
                for one_cmpd_smiles in fwd_smiles_set_level_list[subs_level]:
                    #if one_rule_id == "1.1.1.a.1.f" and subs_level==3:
                        #print one_cmpd_smiles
                        #continue
                    if one_cmpd_smiles.find("CoA")!=-1:
                        continue # CoA cannot be processed.
                    products_list_duplicates=[]
                    products_lists=self.generate_new_reactions([one_cmpd_smiles,], one_rule_id)
                    for one_products_list in products_lists:
                        if one_products_list not in products_list_duplicates:
                            products_list_duplicates.append(one_products_list)
                            unique_smiles_products_list=[]
                            for one_smiles in unique_canonical_smiles_list_zx(one_products_list):
                                if one_smiles not in bkgd_cmpd_list:
                                    unique_smiles_products_list.append(one_smiles)
                            if unique_smiles_products_list != []:
                                new_reactions_list.append(   ( tuple(unique_smiles_products_list) , (one_cmpd_smiles,) , one_rule_id )   )
            #============================================================================================================================#
            #####----------Step 1-3 Type Two Enzyme, reaction coef>1 OR two reactants:
            else:
                #print 'multiple reactants'
                #####----------Step 1-3-1 Initiate templates of enzyme reactants:
                reactant_substructure_list=one_rule_str.split(">>")[0].split(".")
                numreactants=len(reactant_substructure_list)
                numproducts=one_rule_str.split(">>")[1].count(".")+1
                #####----------Step 1-3-2 Initiate compounds combinations list:
                # Initiate good_cmpds_list:
                good_compounds_list=[]
                good_compounds_set_list=[]
                for i in range(numreactants):
                    good_compounds_list.append([])
                #####----------Step 1-3-3 Generate compounds groups that match the reactants' templates:
                # Generate a list of compounds for each reactants:
                for i in range(numreactants):
                    for one_cmpd_smiles in fwd_smiles_set.union(set(bkgd_cmpd_list)):
                        if one_cmpd_smiles.find("CoA")!=-1:
                            continue # CoA cannot be processed.
                        if self.pattern_matching_zx(one_cmpd_smiles, reactant_substructure_list[i]):
                                good_compounds_list[i].append(one_cmpd_smiles)
                num_combinations=get_num_combinations(good_compounds_list)
                #print num_combinations

                num_group=10000                                     ##### if num_combinations is larger than 10^9 than will cause memory error (on windows)
                group_size=int(num_combinations/num_group+1)             ##### if num_combinations is larger than 10^9 than will cause memory error (on windows)

                if numproducts==1: # IF THERE IS ONLY ONE PRODUCT, CAN DO PRE_SCREENING !!! (EFFECTIVELY SAVE RUNTIME!!!)

                    for group_i in range(num_group+1):              ##### if num_combinations is larger than 10^9 than will cause memory error (on windows)
                        for num_group_i in range(group_size):       
                            num_i=group_i*group_size+num_group_i    
                            if num_i>=num_combinations:             
                                continue                            

                            one_set_reactants_list=get_ith_combination(good_compounds_list, num_i)
                            # Pre-screening 
                            if (set(one_set_reactants_list) - set(fwd_smiles_set_level_list[subs_level-1]).union(bkgd_cmpd_list))==set([]) and subs_level != 0:
                                continue
                            if str(one_set_reactants_list).count('C')>max_C_num+1 or str(one_set_reactants_list).count('O')>max_O_num+1:
                                continue
                            products_list_duplicates=[]
                            products_lists=self.generate_new_reactions(one_set_reactants_list, one_rule_id)
                            for one_products_list in products_lists:
                                if one_products_list not in products_list_duplicates:
                                    products_list_duplicates.append(one_products_list)
                                    unique_smiles_products_list=[]
                                    no_bkgd_reactants_list=[]
                                    for one_smiles in unique_canonical_smiles_list_zx(one_products_list):
                                        if one_smiles not in bkgd_cmpd_list:
                                            unique_smiles_products_list.append(one_smiles)
                                    for one_smiles in one_set_reactants_list:
                                        if one_smiles not in bkgd_cmpd_list:
                                            no_bkgd_reactants_list.append(one_smiles)
                                    if unique_smiles_products_list!=[] and no_bkgd_reactants_list!=[]:
                                        new_reactions_list.append(   ( tuple(unique_smiles_products_list) , tuple(no_bkgd_reactants_list) , one_rule_id )   )
                else:

                    for group_i in range(num_group+1):              ##### if num_combinations is larger than 10^9 than will cause memory error (on windows)
                        for num_group_i in range(group_size):       
                            num_i=group_i*group_size+num_group_i    
                            if num_i>=num_combinations:             
                                continue                            
                        
                            one_set_reactants_list=get_ith_combination(good_compounds_list, num_i)
                            # Pre-screening 
                            if (set(one_set_reactants_list) - set(fwd_smiles_set_level_list[subs_level-1]).union(bkgd_cmpd_list))==set([]) and subs_level != 0:
                                continue
                            products_list_duplicates=[]
                            products_lists=self.generate_new_reactions(one_set_reactants_list, one_rule_id)
                            for one_products_list in products_lists:
                                if one_products_list not in products_list_duplicates:
                                    products_list_duplicates.append(one_products_list)
                                    unique_smiles_products_list=[]
                                    no_bkgd_reactants_list=[]
                                    for one_smiles in unique_canonical_smiles_list_zx(one_products_list):
                                        if one_smiles not in bkgd_cmpd_list:
                                            unique_smiles_products_list.append(one_smiles)
                                    for one_smiles in one_set_reactants_list:
                                        if one_smiles not in bkgd_cmpd_list:
                                            no_bkgd_reactants_list.append(one_smiles)
                                    if unique_smiles_products_list!=[] and no_bkgd_reactants_list!=[]:
                                        new_reactions_list.append(   ( tuple(unique_smiles_products_list) , tuple(no_bkgd_reactants_list) , one_rule_id )   )
        #####----------Step 2-0 CoA rxn:
        for one_CoA_rxn in fwd_CoA_rxn:
            if set(one_CoA_rxn[1]).issubset(fwd_smiles_set_level_list[subs_level]):
                print (("CoA rxn found."))
                new_reactions_list.append(one_CoA_rxn)
        #print bad_ss_dict
        print ("Reaction Simulation Finished. ", end="")
        return new_reactions_list

#######################################################################################################################################
#######################################################################################################################################
    def bwd_apply_enzymes_zx(self,cmpds_need_expand,prod_level):
        print ("Reaction Simulation ... ")

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 1 Generate new reactions and store in the new_reaction_list (in one loop):
        new_reactions_list=[]

        for one_rule_id in tqdm(self.selected_rxn_id_list):
            #print ('Prod-Side-Depth', prod_level+1, '-',) #cmpds_need_expand
            #print one_rule_id
            # For backward direction prediction, need to find those rules with only one product and 
            # ==================================================================================== #
            #####----------Step 1-1 Obtain enzyme datsa: 
            one_rule_str=self.selected_rxn_dict[one_rule_id][0]
            one_rule_ec_num=self.selected_rxn_dict[one_rule_id][1]
            one_rule_rxn_order=self.selected_rxn_dict[one_rule_id][2]
            one_rule_diameter=self.selected_rxn_dict[one_rule_id][3]
            one_rule_score=self.selected_rxn_dict[one_rule_id][4]
            one_rule_tag=self.selected_rxn_dict[one_rule_id][5]
            if one_rule_tag =="forward": # skip if not backward (not applicable to APrules)
                continue
            # ==================================================================================== #
            #####----------Step 1-2 Type One Enzyme, one reactant: 
            if one_rule_str.split(">>")[0].count(".")==0: # There is only one product for this reaction rule
                #print one_rule_str
                for one_cmpd_smiles in cmpds_need_expand:
                    if one_cmpd_smiles.find("CoA")!=-1:
                        continue # CoA cannot be processed
                    products_list_duplicates=[]
                    products_lists=self.generate_new_reactions([one_cmpd_smiles,], one_rule_id)
                    for one_products_list in products_lists:
                        if one_products_list not in products_list_duplicates:
                            products_list_duplicates.append(one_products_list)
                            unique_smiles_products_list=[]
                            for one_smiles in unique_canonical_smiles_list_zx(one_products_list):
                                if one_smiles not in bkgd_cmpd_list:
                                    unique_smiles_products_list.append(one_smiles)
                            if unique_smiles_products_list != []:
                                new_reactions_list.append(   ( tuple(unique_smiles_products_list) , (one_cmpd_smiles,) , one_rule_id )   )
            else:
            #####----------Step 1-3 Type Two Enzyme, reaction coef>1 OR two reactants:
                #print 'multiple reactants'
                #####----------Step 1-3-1 Initiate templates of enzyme reactants:
                reactant_substructure_list=one_rule_str.split(">>")[0].split(".")
                numreactants=len(reactant_substructure_list)
                #####----------Step 1-3-2 Initiate compounds combinations list:
                # Initiate good_cmpds_list:
                good_compounds_list=[]
                good_compounds_set_list=[]
                for i in range(numreactants):
                    good_compounds_list.append([])
                #####----------Step 1-3-3 Generate compounds groups that match the reactants' templates:
                # Generate a list of compounds for each reactants:
                #for one_cmpd_need_expand in cmpds_need_expand:
                for i in range(numreactants):
                    for one_cmpd_smiles in set(bkgd_cmpd_list).union(set(cmpds_need_expand)):
                        if one_cmpd_smiles.find("CoA")!=-1:
                            continue # CoA cannot be processed
                        if self.pattern_matching_zx(one_cmpd_smiles, reactant_substructure_list[i]):
                                good_compounds_list[i].append(one_cmpd_smiles)
                num_combinations=get_num_combinations(good_compounds_list)
                #print num_combinations
                for num_i in range(num_combinations):
                    one_set_reactants_list=get_ith_combination(good_compounds_list, num_i)
                    if len(set(one_set_reactants_list)-set(bkgd_cmpd_list))!=1:
                        continue
                    one_cmpd_need_expand=list(set(one_set_reactants_list)-set(bkgd_cmpd_list))[0]
                    products_list_duplicates=[]
                    products_lists=self.generate_new_reactions(one_set_reactants_list, one_rule_id)
                    for one_products_list in products_lists:
                        if one_products_list not in products_list_duplicates:
                            products_list_duplicates.append(one_products_list)
                            unique_smiles_products_list=[]

                            for one_smiles in unique_canonical_smiles_list_zx(one_products_list):
                                if one_smiles not in bkgd_cmpd_list:
                                    unique_smiles_products_list.append(one_smiles)
                            if unique_smiles_products_list!=[]:
                                new_reactions_list.append(   ( tuple(unique_smiles_products_list) , (one_cmpd_need_expand,) , one_rule_id )   ) 

        #####----------Step 2-0 CoA rxn:
        for one_CoA_rxn in bwd_CoA_rxn:
            if len(one_CoA_rxn[1])==1 and one_CoA_rxn[1][0] in cmpds_need_expand:
                print ("CoA rxn found.")
                new_reactions_list.append(one_CoA_rxn)
        #print new_reactions_list
        print ("Reaction Simulation Finished. ", end="")
        return new_reactions_list


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def test_x(self):
        for rxn_id in self.selected_rxn_id_list:
            if self.selected_rxn_dict[rxn_id][2]==1:
                result=self.generate_new_reactions(["O[C@H]1[C@H](O)[C@@H](COP(O)(O)=O)OC(O)[C@@H]1O",],rxn_id)
                if result  !=[]:
                    print (result)
        return


#######################################################################################################################################
#######################################################################################################################################
    def test_z(self):

        fwd_smiles_set=set(['O=CC(O)COP(=O)(O)O',])
        fwd_smiles_set_level_list=[list(fwd_smiles_set),[]]
        max_C_num=6
        max_O_num=6
        subs_level=0
        new_rct_list= self.apply_enzymes(fwd_smiles_set,fwd_smiles_set_level_list,subs_level,3,3)
        #new_rct_list= self.bwd_apply_enzymes_zx(fwd_smiles_set,subs_level)
        print (new_rct_list)
        print (len(set(new_rct_list)))

        return

#######################################################################################################################################
#######################################################################################################################################

    def test_w(self):
        cmpd_a='C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N'  # CC(C)(=O)C(O)C(O)C(=O)O
        rule_a='MNXR55960_MNXM960-2-fwd'
        rule_str=AP_reactor.selected_rxn_dict['MNXR55960_MNXM960-2-fwd'][0]
        print (rule_str)
        rxn = AllChem.ReactionFromSmarts(rule_str)
        rf_cmpd= Chem.MolFromSmiles(cmpd_a)
        ps = rxn.RunReactants((rf_cmpd,))
        rf_prod=ps[0][0]
        print (rf_prod )
        bad_ss_dict=dict([])
        r_fmt_cmpd_x=rf_prod
        print ("start!")
        print (MolToSmiles_zx(r_fmt_cmpd_x, bad_ss_dict))
        print (bad_ss_dict)
        print ("finished!")

        print (Chem.MolToSmiles(rf_prod))
        print (Chem.MolToSmarts(rf_prod))
        rf_cmpd= Chem.MolFromSmarts(cmpd_a)
    
        Draw.MolToFile(rf_prod,'cdk2_mol1.o.png',size=(500,500)) # depending on molecule size, shall adjust the img size input
        Draw.MolToFile(rf_cmpd,'cdk2_mol2.o.png',size=(500,500)) # depending on molecule size, shall adjust the img size input
        print ("done")
        return


    def test_ww(self):
        rule="([#16]-[#6:2]([C:4])=[O:1])>>([H][#6:2]([C:4])=[O:1])"
        smiles_a="[H][#6:2]([C:4])=[O:1]"
        print (AP_reactor.pattern_matching_zx("CCO","[CH,CH2,CH3:2][C:4][OH:1]"))
        print (AP_reactor.pattern_matching_zx("O","[O:1]-[C:2]"))


#######################################################################################################################################
#######################################################################################################################################
    def test_y(self):
        a='O=C(O)C(=O)CC(=O)O'
        b="O[C@H]1[C@H](O)[C@@H](COP(O)(O)=O)OC(O)[C@@H]1O"
        for rxn_id in self.selected_rxn_id_list:
            if self.selected_rxn_dict[rxn_id][0].split(">>")[1].count(".")==0:
                if self.pattern_matching_zx(a, self.selected_rxn_dict[rxn_id][0].split(">>")[1]):
                    print (self.selected_rxn_dict[rxn_id][5])
                result=self.generate_new_reactions_reverse([a,],rxn_id)
                if result  !=[]:
                    print (result)

            if self.selected_rxn_dict[rxn_id][0].split(">>")[0].count(".")==0:
                if self.pattern_matching_zx(a, self.selected_rxn_dict[rxn_id][0].split(">>")[0]):
                    print (self.selected_rxn_dict[rxn_id][5])
                result=self.generate_new_reactions([a,],rxn_id)
                if result  !=[]:
                    print (result)

        return
#######################################################################################################################################
#######################################################################################################################################

if __name__ == '__main__':

    AP_reactor = AP_Reactor()
    print (AP_reactor.pattern_matching_zx("CC(=O)O","[CH3:1][C:2](=[O:3])[OH:4]"))
    print (AP_reactor.pattern_matching_zx("CC(=O)OP(O)(O)=O","[C:2][O:6][P:7]([OH:8])([OH:9])=[O:10]"))

    #AP_reactor.test_ww()

    print ("done")

