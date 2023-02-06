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
import math
import itertools
import operator
import time
import random
import subprocess
import pymysql.cursors
#--------------------------------------------------#
from tqdm import tqdm
from copy import copy, deepcopy
#--------------------------------------------------#
from AP_funcs import *
from AP_convert import * 
from AP_convert import unique_canonical_smiles_AP as unis
from AP_reactor import *
from AP_output import *
#--------------------------------------------------#
global bkgd_cmpd_list; bkgd_cmpd_list=bkgd_cmpd_list_func()
global CoA_cmpd_list; CoA_cmpd_list=CoA_cmpd_list_func()
global unfavorable_rxn_list; unfavorable_rxn_list=unfavorable_rxn_list_func()
global XI_global; XI_global=0.2
global T0_global; T0_global=10

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
class AnnealPath: 
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def __init__(self, outstream=sys.stderr, rxn_diameter=2, rxn_score_lb=0, rxn_db="AP"):
        self.outstream = outstream
        self.rxn_db = rxn_db
        self.reactor = AP_Reactor(rxn_diameter=rxn_diameter, rxn_score_lb=rxn_score_lb, rxn_db=rxn_db)
    # ============================================================================================================================ #
    def pathway_found_textart(self):
        thumb_up=['        $$$$', '       $$  $', '       $   $$', '       $   $$', '       $$   $$', '        $    $$', '        $$    $$$', '         $$     $$', '         $$      $$', '          $       $$', '    $$$$$$$        $$', '  $$$               $$$$$$', ' $$    $$$$            $$$', ' $   $$$  $$$            $$', ' $$        $$$            $', '  $$    $$$$$$            $', '  $$$$$$$    $$           $', '  $$       $$$$           $', '   $$$$$$$$$  $$         $$', '    $        $$$$     $$$$', '    $$    $$$$$$    $$$$$$', '     $$$$$$    $$  $$', '       $     $$$ $$$', '        $$$$$$$$$$']
        for line in thumb_up:
            print(line)
        return
    # ============================================================================================================================ #
    def KEGG_nme_canonical_SMILES_dict(self):
        pickle_in1=open("./KEGGScrapSavings/KEGG_nme_canonical_SMILES_dict.pickle","rb")
        KEGG_nme_canonical_SMILES_dict=pickle.load(pickle_in1)
        pickle_in1.close()
        return KEGG_nme_canonical_SMILES_dict

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def pathway_construction(self, all_lists):
        #--------------------------------------------------#
        # Information need:
        # cmplt_pathways_list
        # target_pathways_list
        # subs_side_compoundshash_list
        # reactions_list
        # fwd_smiles_list_level_list
        #--------------------------------------------------#
        # Will have to deal with three situations: 
        # 1. pathways already found on subs side searching
        # 2. pathways already found on searching prod side
        # 3. combine the searching results from both searching direction and check the existence of pathways
        #######################################################################################################################################
        #######################################################################################################################################
        def subs_side_pathways(target_pathways_list, subs_side_compoundshash_list, reactions_list): 
            all_pathways_list=[]
            for i in range(subs_level-target_compound_level+1):#1
                pathway_length=target_compound_level+i
                print ("pathway length: ")
                print (pathway_length)
                pathways_list=[]
                # now have pathway_length and target_compound_level
                # give the very first reactions to the pathways_list and start to expand
                for j in range(len(reactions_list)):#1.1
                    if (target_compound in reactions_list[j][0]) and (iftuplestrinlist(reactions_list[j][1],fwd_smiles_list_level_list[pathway_length-1])):
                        pathways_list.append([(pathway_length,reactions_list[j][0],reactions_list[j][1],reactions_list[j][2]),])
                for k in range(pathway_length-1):#1.2
                    l=pathway_length-1-k
                    temp_pathways_list=[]
                    # ============================================================================================================================ #
                    for m in range(len(pathways_list)):#1.2.1
                        temp_pathways_list.append(pathways_list[m])
                    # ============================================================================================================================ #
                    for n in range(len(temp_pathways_list)):#1.2.2
                        pathways_list.remove(temp_pathways_list[n])
                        reactantstobeexpanded=[]
                        for o in range(len(temp_pathways_list[n])):#1.2.2.1
                            if temp_pathways_list[n][o][0]==l+1:
                                for p in range(len(temp_pathways_list[n][o][2])):#1.2.2.1.1
                                    if temp_pathways_list[n][o][2][p] not in fwd_smiles_list_level_list[0]:
                                        reactantstobeexpanded.append(temp_pathways_list[n][o][2][p])
                        #----------------------------------------------------------------------------------------------------#
                        next_level_reactions=[]
                        #----------------------------------------------------------------------------------------------------#
                        if len(reactantstobeexpanded)>0:
                            for q in range(len(reactantstobeexpanded)):#1.2.2.2
                                next_level_reactions.append([])
                                for r in range(len(reactions_list)):#1.2.2.2.1
                                    if (reactantstobeexpanded[q] in reactions_list[r][0]) and iftuplestrinlist(reactions_list[r][1],fwd_smiles_list_level_list[l-1]):
                                        next_level_reactions[q].append((l,reactions_list[r][0],reactions_list[r][1],reactions_list[r][2]))
                            next_level_reactions_set=cart_prod(next_level_reactions) ###
                            for u in range(len(next_level_reactions_set)):#1.2.2.3
                                one_pathway_list=deepcopy(temp_pathways_list[n])
                                for v in range(len(next_level_reactions_set[u])):#1.2.2.3.1
                                    if next_level_reactions_set[u][v] not in one_pathway_list:
                                        one_pathway_list.append(next_level_reactions_set[u][v])
                                if one_pathway_list not in pathways_list:
                                    pathways_list.append(one_pathway_list)
                        #----------------------------------------------------------------------------------------------------#
                        if len(reactantstobeexpanded)==0:
                            one_pathway_list=temp_pathways_list[n]
                            length_shorten=one_pathway_list[-1][0]-1
                            for a in range(len(one_pathway_list)):#1.2.2.4
                                temp_thing=list(one_pathway_list[a])
                                temp_thing[0]=temp_thing[0]-length_shorten
                                one_pathway_list[a]=tuple(temp_thing)
                            if one_pathway_list not in pathways_list:
                                pathways_list.insert(0,one_pathway_list)
                # ============================================================================================================================ #
                for r in range(len(pathways_list)):#1.3
                    all_pathways_list.append(pathways_list[r])
            return all_pathways_list

        #######################################################################################################################################
        #######################################################################################################################################
        def convert_bwd_pathways(cmplt_pathways_list,number_to_add):
            # Used to return a complete pathway found in backward searching alone (prod_side_pathways)
            # Or to return an incomplete pathway (prod side completed) for subs side construction
            result_list=[]
            if number_to_add==0:
                for i in range(len(cmplt_pathways_list)):
                    one_pathway_list=[]
                    number_tb_add=1-cmplt_pathways_list[i][-1][0]
                    for j in range(len(cmplt_pathways_list[i])):
                        one_pathway_list.append((cmplt_pathways_list[i][j][0]+number_tb_add,cmplt_pathways_list[i][j][1],cmplt_pathways_list[i][j][2],cmplt_pathways_list[i][j][3]))
                    result_list.append(one_pathway_list)
            else:
                for i in range(len(cmplt_pathways_list)):
                    one_pathway_list=[]
                    number_tb_add=number_to_add
                    for j in range(len(cmplt_pathways_list[i])):
                        one_pathway_list.append((cmplt_pathways_list[i][j][0]+number_tb_add,cmplt_pathways_list[i][j][1],cmplt_pathways_list[i][j][2],cmplt_pathways_list[i][j][3]))
                    result_list.append(one_pathway_list)
            return result_list

        #######################################################################################################################################
        #######################################################################################################################################
        def let_the_pathway_form(target_pathways_list,reactions_list,fwd_smiles_list_level_list):
            all_pathways_list=[]
            #print prod_level
            for i in range(prod_level): # i helps get different length of pathways, pathway_length=i+1+subs_level
                pathways_list=convert_bwd_pathways(target_pathways_list[i],i+2+subs_level)
                # ============================================================================================================================ #
                for j in range(subs_level):
                    print("subs_level: ", j, "prod_level: ", i)
                    l=subs_level-j # l is the level of cmpds to be searched now
                    temp_pathways_list=copy(pathways_list)
                    all_rctts_tb_expdd=set([]) # all_rctts_tb_expdd is the reactants to be substituted by lower level compounds
                    reactants_expansion_dict=dict([])
                    #----------------------------------------------------------------------------------------------------#
                    for one_icmplt_pathway in temp_pathways_list: # n means processing (doing substitution to) the n^th pathway in the pathway list
                        pathways_list.remove(one_icmplt_pathway)

                        for one_itmdt_rxn in one_icmplt_pathway: # o means the o^th reaction in the reaction list of the n^th pathway
                            if one_itmdt_rxn[0]==l+1:
                                for one_rctt_tb_expdd in one_itmdt_rxn[2]:#1.2.2.1.1
                                    if one_rctt_tb_expdd not in fwd_smiles_list_level_list[0]:
                                        all_rctts_tb_expdd.add(one_rctt_tb_expdd)
                    #print all_rctts_tb_expdd
                    #----------------------------------------------------------------------------------------------------#
                    for one_rctt_tb_expdd in all_rctts_tb_expdd: #1.2.2.2
                        reactants_expansion_dict[one_rctt_tb_expdd]=[]
                        for one_rxn in reactions_list: #1.2.2.2.1
                            if (one_rctt_tb_expdd in one_rxn[0]) and \
                                iftuplestrinlist(one_rxn[1],fwd_smiles_list_level_list[l-1]):
                                    reactants_expansion_dict[one_rctt_tb_expdd].append((l,one_rxn[0],one_rxn[1],one_rxn[2]))
                    #----------------------------------------------------------------------------------------------------#
                    for one_icmplt_pathway in temp_pathways_list: # n means processing (doing substitution to) the n^th pathway in the pathway list
                        reactants_found_already=[] # reactants_found_already is checked to erase meaningless pathways
                        reactants_found_already_1=set([])
                        for one_itmdt_rxn in one_icmplt_pathway: # o means the o^th reaction in the reaction list of the n^th pathway
                            reactants_found_already.append(tuple(set(one_itmdt_rxn[1])))
                            #reactants_found_already_1=reactants_found_already_1.union(set(one_itmdt_rxn[1]))
                            #reactants_found_already_1=reactants_found_already_1.union(set(one_itmdt_rxn[2]))
                        #print reactants_found_already
                        #--------------------------------------------------#
                        reactantstobeexpanded=[] # reactantstobeexpanded is the reactants to be substituted by lower level compounds
                        for one_itmdt_rxn in one_icmplt_pathway: # o means the o^th reaction in the reaction list of the n^th pathway
                            if one_itmdt_rxn[0]==l+1:
                                for one_rctt_tb_expdd in one_itmdt_rxn[2]:#1.2.2.1.1
                                    if one_rctt_tb_expdd not in fwd_smiles_list_level_list[0]:
                                        reactantstobeexpanded.append(one_rctt_tb_expdd)
                        #print reactantstobeexpanded
                        reactantstobeexpanded=list(set(reactantstobeexpanded))
                        #--------------------------------------------------#
                        if reactantstobeexpanded!=[]:
                            next_level_reactions=[]
                            for one_rctt_tb_expdd in reactantstobeexpanded: #1.2.2.2
                                one_set_next_level_reactions=[]
                                for one_possible_reaction in reactants_expansion_dict[one_rctt_tb_expdd]:
                                    if (tuple(set(one_possible_reaction[2])) not in reactants_found_already):
                                        one_set_next_level_reactions.append(one_possible_reaction)
                                next_level_reactions.append(one_set_next_level_reactions)
                            next_level_reactions_set=cart_prod(next_level_reactions)
                            #print next_level_reactions
                            #print next_level_reactions_set
                            for one_next_lv_rxn_set in next_level_reactions_set:#1.2.2.4
                                one_pathway_list=copy(one_icmplt_pathway)
                                for one_next_lv_rxn in one_next_lv_rxn_set:#1.2.2.4.1
                                    if one_next_lv_rxn not in one_pathway_list:
                                        one_pathway_list.append(one_next_lv_rxn)
                                if one_pathway_list not in pathways_list:
                                    pathways_list.append(one_pathway_list)
                            #print pathways_list
                        #--------------------------------------------------#
                        else:
                            one_pathway_list=one_icmplt_pathway
                            length_shorten=one_pathway_list[-1][0]-1
                            for a in range(len(one_pathway_list)):
                                temp_thing=list(one_pathway_list[a])
                                temp_thing[0]=temp_thing[0]-length_shorten
                                one_pathway_list[a]=tuple(temp_thing)
                            if one_pathway_list not in pathways_list:
                                pathways_list.insert(0,one_pathway_list)
                # ============================================================================================================================ #
                for w in range(len(pathways_list)):
                    all_pathways_list.append(pathways_list[w])
            #print 'pathways constructed by results of both fwd and bwd search: '
            #print all_pathways_list
            return all_pathways_list
            
        #######################################################################################################################################
        #######################################################################################################################################
        # Main part of pathway_construction():
        # Get the search results
        cmplt_pathways_list=all_lists[0]
        target_pathways_list=all_lists[1]
        subs_side_compoundshash_list=all_lists[2]
        reactions_list=all_lists[3]
        fwd_smiles_list_level_list=all_lists[4]
        smilesfirstacceptlevel=all_lists[5]
        '''
        print "all_lists: "
        print " complete_pathways_list: ", all_lists[0]
        print " target_pathways_list: ", all_lists[1]
        print " list(fwd_smiles_set): ", all_lists[2]
        print " list(fwd_rxn_set): ", all_lists[3]
        print " fwd_smiles_list_level_list: ", all_lists[4]
        print " smilesfirstacceptlevel: ", all_lists[5]
        '''
        # ============================================================================================================================ #
        final_pathways_list=[]
        # Three situations:
        if found_on_subs_side==True:
            final_pathways_list=subs_side_pathways(target_pathways_list, subs_side_compoundshash_list, reactions_list)
        print ("Any pwys already found through fwd search alone: ")
        print (final_pathways_list)
        print ("Any pwys already found through bwd search alone: ")
        print (cmplt_pathways_list)
        # ============================================================================================================================ #
        if cmplt_pathways_list!=[]:
            prod_side_pathways=convert_bwd_pathways(cmplt_pathways_list,0)
            for i in range(len(prod_side_pathways)):
                if prod_side_pathways[i] not in final_pathways_list:
                    final_pathways_list.append(prod_side_pathways[i])
        # ============================================================================================================================ #
        print("Now construct all complete pathways through bridge compounds ...")
        both_side_pathways=let_the_pathway_form(target_pathways_list,reactions_list,fwd_smiles_list_level_list)
        print ("Final pathways constructed by double direction search: ")
        # ============================================================================================================================ #
        for j in range(len(both_side_pathways)):
            if both_side_pathways[j] not in final_pathways_list:
                final_pathways_list.append(both_side_pathways[j])
                #print both_side_pathways[j]
        #print final_pathways_list
        return final_pathways_list

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def expand_reaction_tree(self,fwd_smiles_set,fwd_rxn_set,fwd_smiles_set_level_list,trfm_odict,ALL_RXN_SET):
        # Expand the following objects at each level AND UPDATE THE fwd_probabilities_dict:
        # fwd_smiles_set            : intermediate compounds that are ACCEPTED (FOUND and SELECTED), used as reactants at next level
        # fwd_rxn_set             : expands at each level, all reactions FOUND are stored in it
        # processed_rxn_list      : used in reactor.apply_enzymes(), all reactants set are stored in it.
        # fwd_smiles_set_level_list     : all compounds ACCEPTED grouped in different levels, 
        # Could have used the global variables. Those inputs are just for tests.
        #######################################################################################################################################
        #######################################################################################################################################
        def fwd_probability_two_cmpds(target_cmpd,reactant,prod):
            # 1. Get similarity score for the reactant and product.
            reactant_score=similarity_score(target_cmpd,reactant,fp_type_global)
            prod_score=similarity_score(target_cmpd,prod,fp_type_global)
            # 2. Calculate probability based on the two scores generated:
            if prod_score>=reactant_score:
                probability=1
            else:
                probability=math.exp((prod_score-reactant_score)/T0_global/(XI_global**subs_level))
            #print probability, 
            return probability
	    # ============================================================================================================================ #
        def probability_select(fwd_smiles_set_level_list,cmpds_list,rxn_list,target_cmpd,num_cmpds): # Top ranked similarity improvement
            #print ("e^(-delta(sim_score)/(T0^subs_level))", "subs_level=", subs_level)
            # Return a list of compounds according to probability rankings, in other words, similarity improvement
            # 1. Update the probability_dict and the global fwd_probabilities_dict directly:
            probability_dict={}
            top_similarity_improvement=[] # Output list.
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        for reactant in rxn[1]:
                            if reactant in fwd_smiles_set_level_list[subs_level]:
                                probability=fwd_probability_two_cmpds(target_cmpd,reactant,prod)

                                if (prod in fwd_probabilities_dict.keys() and probability > fwd_probabilities_dict[prod]) or (prod not in fwd_probabilities_dict.keys()):
                                    fwd_probabilities_dict[prod]=probability
            return top_similarity_improvement,probability_dict
	    # ============================================================================================================================ #
        def similarity_improvement_select_7(fwd_smiles_set_level_list,cmpds_list,rxn_list,target_cmpd,num_cmpds): # Top ranked similarity improvement difference (half number)
            similarity_improvement_diff_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        for reactant in rxn[1]:
                            if reactant in fwd_smiles_set_level_list[subs_level]:
                                similarity_improvement_diff=similarity_score(target_cmpd,prod,fp_type_global)-similarity_score(target_cmpd,reactant,fp_type_global)

                                if (prod in similarity_improvement_diff_dict.keys() and similarity_improvement_diff > similarity_improvement_diff_dict[prod]) \
                                    or (prod not in similarity_improvement_diff_dict.keys()):
                                    similarity_improvement_diff_dict[prod] = round(similarity_improvement_diff, 3)

            sorted_similarity_improvement_diff_list = sorted(similarity_improvement_diff_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print ("sorted_similarity_improvement_diff_list: ")
            #print (sorted_similarity_improvement_diff_list)
            top_similarity_improvement_diff=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_similarity_improvement_diff_list) else len(sorted_similarity_improvement_diff_list)
            for j in range(num_cmpds):
                if sorted_similarity_improvement_diff_list[j][1] != -1:
                    top_similarity_improvement_diff.append(sorted_similarity_improvement_diff_list[j][0])
                else:
                    break
            # 3. Output list of compounds hashes
            return top_similarity_improvement_diff
	    # ============================================================================================================================ #
        def similarity_improvement_select_8(fwd_smiles_set_level_list,cmpds_list,rxn_list,target_cmpd,num_cmpds): # Top ranked similarity improvement ratio (half number)
            similarity_improvement_ratio_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        for reactant in rxn[1]:
                            if reactant in fwd_smiles_set_level_list[subs_level]:
                                similarity_improvement_ratio=similarity_score(target_cmpd,prod,fp_type_global)/(similarity_score(target_cmpd,reactant,fp_type_global)+0.00001)

                                if (prod in similarity_improvement_ratio_dict.keys() and similarity_improvement_ratio > similarity_improvement_ratio_dict[prod]) \
                                    or (prod not in similarity_improvement_ratio_dict.keys()):
                                    similarity_improvement_ratio_dict[prod]=similarity_improvement_ratio

            sorted_similarity_improvement_ratio_list = sorted(similarity_improvement_ratio_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print ("sorted_similarity_improvement_ratio_list: ")
            #print (sorted_similarity_improvement_ratio_list)
            top_similarity_improvement_ratio=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_similarity_improvement_ratio_list) else len(sorted_similarity_improvement_ratio_list)
            for j in range(num_cmpds):
                if sorted_similarity_improvement_ratio_list[j][1] != 0:
                    top_similarity_improvement_ratio.append(sorted_similarity_improvement_ratio_list[j][0])
                else:
                    break
            # 3. Output list of compounds hashes
            return top_similarity_improvement_ratio
        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 0-0 Generate new compounds
        new_rxn_list=self.reactor.apply_enzymes(fwd_smiles_set,fwd_smiles_set_level_list,subs_level,max_C_num,max_O_num)
        #####----------Step 0-1 Update trfm_rxn_dict & trfm_rule_dict.  Modify new_rxn_list. (new version)
        modified_new_rxn_set=set([])
        fwd_current_smiles_set=set([])
        fwd_current_accepted_smiles_set=set([])
        for one_rxn in new_rxn_list:
            if tuple([one_rxn[0],one_rxn[1],]) not in unfavorable_rxn_list:
                one_trfm=(one_rxn[0],one_rxn[1])
                if trfm_odict[one_trfm]==0:
                    trfm_odict[one_trfm]=set([one_rxn[2],])
                    trfm_id=len(trfm_odict)-1
                else:
                    trfm_odict[one_trfm].add(one_rxn[2])
                    trfm_id=list(trfm_odict.keys()).index(one_trfm)
                one_modified_rxn=(one_rxn[0],one_rxn[1],trfm_id)

                modified_new_rxn_set.add(one_modified_rxn)
                fwd_rxn_set.add(one_modified_rxn)
                fwd_current_smiles_set=fwd_current_smiles_set.union(set(one_modified_rxn[0]))
        new_rxn_list=list(modified_new_rxn_set)
        print ("Number of new reactions: " + str(len(new_rxn_list)) +". ", end="")
        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 0-2 Select a number of compounds based on SI and SA (assuming only one target compound here)
        #####----------Step 0-2-0 If last subs level, no need for screening or selecting compounds for next level
        if subs_level ==( int((max_levels_global+1)/2) - 1 ): # If last subs level
            print ("last level, skip compound selection")
            for one_smiles in fwd_current_smiles_set:
                if one_smiles not in fwd_smiles_set:
                    fwd_current_accepted_smiles_set.add(one_smiles)
                    fwd_smiles_set.add(one_smiles)
        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 0-2-1 If NOT last subs level, Screen and Select.
        if subs_level !=( int((max_levels_global+1)/2) - 1 ): # If NOT last subs level, Screen.
            screened_smiles_list=[]
            # ============================================================================================================================ #
            #####----------Step 0-2-1-0 Pre-Screening: Remove too large intermediate compounds (have a lot more C or a lot more O atoms than the target compound)
            if P0_global==[0,0,0]:
                for one_smiles in fwd_current_smiles_set:
                    one_smiles_expand=one_smiles
                    if one_smiles.find("I")!=-1:
                        one_smiles_expand=one_smiles.replace("(I)", "(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
                        one_smiles_expand=one_smiles_expand.replace("I", "(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
                        one_smiles_expand=unis(one_smiles_expand)
                    
                    if one_smiles_expand in KEGG_nme_canonical_SMILES_dict.keys() and one_smiles_expand not in bkgd_cmpd_list:
                        #print ("KEGG cmpd selected, ", KEGG_nme_canonical_SMILES_dict[one_smiles_expand], one_smiles_expand)
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                    if one_smiles not in bkgd_cmpd_list and one_smiles not in bad_ss_dict.keys():
                        if one_smiles.count("C")<=max_C_num and \
                            one_smiles.count("O")<=max_O_num:
                            screened_smiles_list.append(one_smiles)
            else:
                for one_smiles in fwd_current_smiles_set:
                    if one_smiles not in bkgd_cmpd_list and one_smiles not in bad_ss_dict.keys():
                        if one_smiles.count("C")<=max_C_num and \
                            one_smiles.count("O")<=max_O_num:
                            screened_smiles_list.append(one_smiles)
            #print "number of compounds screened on searching this level (subs side): " + str(len(screened_smiles_list))
            # ============================================================================================================================ #
            #####----------Step 0-2-1-1 Select
            print ("Number of new compounds: " + str(len(screened_smiles_list)) + ".")
            #####----------Step 0-2-1-1-0 If first subs level, no need to select compounds
            if subs_level==0:
                for one_smiles in fwd_current_smiles_set:
                    if one_smiles not in fwd_smiles_set:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
            #####----------Step 0-2-1-1-1 If NOT first subs level, select compounds
            else:
                #----------------------------------------------------------------------------------------------------#
                # 0. Determine the number of compounds to accept based on the P0_global value
                if P0_global==[0,0,0]: # Select half number of compounds from three groups
                    SISAaccept=len(screened_smiles_list)
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                elif P0_global==[2,2,2]: # Select all compounds
                    SISAaccept=len(screened_smiles_list)
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[-8,-8,-8]: # Select half number of compounds with higher similarity improvement ratios #(used to be x0.5)
                    SISAaccept=len(screened_smiles_list)
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),0,0]
                elif P0_global==[-9,-9,-9]: # Select half number of compounds from three groups #(used to be x0.167)
                    SISAaccept=len(screened_smiles_list)
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                else: # Select certain numbers of compounds from the three groups (User-specified)
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[min(int(SISAaccept*bin_adj_global[0]),P0_global[0]),\
                                            min(int(SISAaccept*bin_adj_global[1]),P0_global[1]),\
                                            min(int(SISAaccept*bin_adj_global[2]),P0_global[2])]
                #----------------------------------------------------------------------------------------------------#
                if P0_global==[0,0,0]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    # Since this is the major bin for selecting compounds, KEGG cmpds acceptance are placed here.
                    count_x=num_screening_global[0]
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=similarity_improvement_select_7(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,count_x)
                    print ("Similarity improvement: ", len(top_similarity_improvement))
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)


                    # Probability Select (select according to the probability dict)
                    (_,_)=probability_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    Probabilities_select=[]
                    random_order_list=randomList(list(fwd_probabilities_dict.keys()))
                    count_y=num_screening_global[1]
                    
                    for one_smiles in random_order_list*10:
                        if count_y==0:
                            break
                        if fwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in fwd_smiles_set:
                            count_y=count_y-1
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print ("Probabilities Select: ", len(Probabilities_select))
                    # Similarity Select (based on similarity scores)
                    Similarities_select=[]
                    taofactors_dict=similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print ("Similarities Select: ", len(Similarities_select))
                #----------------------------------------------------------------------------------------------------#
                elif P0_global==[2,2,2]: # Select all compounds
                    for one_smiles in screened_smiles_list:
                        if one_smiles not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                    print ("Accepted All:", len(screened_smiles_list))
                #----------------------------------------------------------------------------------------------------#
                elif P0_global==[-8,-8,-8]:
                    # Similarity Improvement ( # Select half number of compounds with higher similarity improvements #(***) )
                    top_similarity_improvement=similarity_improvement_select_8(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print ("Similarity improvement: ", len(top_similarity_improvement))
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                #----------------------------------------------------------------------------------------------------#
                elif P0_global==[-9,-9,-9]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=similarity_improvement_select_7(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print ("Similarity improvement: ", len(top_similarity_improvement))
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)

                    # Probability Select (select according to the probability dict)
                    (top_similarity_improvement,probability_dict)=probability_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    Probabilities_select=[]
                    random_order_list=randomList(list(fwd_probabilities_dict.keys()))
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if fwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in fwd_smiles_set:
                            count_x=count_x-1
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print ("Probabilities Select: ", len(Probabilities_select))
                    # Similarity Select (based on similarity scores)
                    Similarities_select=[]
                    taofactors_dict=similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print ("Similarities Select: ", len(Similarities_select))
                #----------------------------------------------------------------------------------------------------#
                else:
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement,probability_dict)=probability_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print ("Similarity improvement: ", len(top_similarity_improvement))
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)

                    # Probability Select (select according to the probability dict)
                    Probabilities_select=[]
                    random_order_list=randomList(list(fwd_probabilities_dict.keys()))
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if fwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in fwd_smiles_set:
                            count_x=count_x-1
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print ("Probabilities Select: ", len(Probabilities_select))
                    # Similarity Select (based on similarity scores)

                    Similarities_select=[]
                    taofactors_dict=similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print ("Similarities Select: ", len(Similarities_select))
        #######################################################################################################################################
        #######################################################################################################################################
        fwd_smiles_set_level_list[subs_level+1] = fwd_current_accepted_smiles_set
        return

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def bwd_expand_reaction_tree(self,bwd_pathways_list,cmplt_pathways_list,bwd_rxn_list,bwd_level_dict,trfm_odict,ALL_RXN_SET):
        # Expand (perform backward substitution to) the following objects at each level AND UPDATE THE bwd_probabilities_dict:
        # bwd_pathways_list  : shown in the long comments below
        # cmplt_pathways_lis : ???
        # bwd_rxn_list       : expands at each level, all reactions FOUND are stored in it
        # bwd_level_dict     : ???
        # bwd_pathways_list:
            #                [
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme') ],
            #                        [ ], [ ], ... [ ]
            #                    ],
            #
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), ( ), ... ( )  ], [ ], ... [ ],
            #                        [ (-2, ('target_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme') ], [ ], ... [ ]
            #                    ],
            #
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'),
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), ( ), ... ( ), 
            #                                (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), 
            #                                (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), ( ), ... ( ) ], [ ], ... [ ],
            #                        [ (-2, ('target_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'),
            #                            (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), 
            #                            (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), ( ), ... ( ) ], [ ], ... [ ],
            #                        [ (-3, ('target_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme') ], [ ], ... [ ]
            #                    ],
            #                    ......
            #                ]
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        def bwd_SA_probability(hash_a,rxn_list,starting_cmpds,bwd_level_dict):
            # Return a simulated probability for hash_a, which is a compound found on current level AND UPDATE THE GLOBAL bwd_probabilities_dict
            # 1. Use the following initiated objects to get one similarity score for the reactants (can be many reactants) of one FOUND compound.
            reactants_list=[] # For each reaction contains hash_a as a product, get the reactants.
            simindex_list1=[] # For all reactants in reactants_list, get the similarity score, select the minimum one
            simindex_list2=[] # Get similarity score between hash_a and all starting compounds
            final_reactant_score=0
            product_score=0
            for i in range(len(rxn_list)):
                if hash_a in rxn_list[i][0]:
                    if rxn_list[i][1][0] not in reactants_list: # There shall be only one reactant (for all reactions in rxn_list).
                        reactants_list.append(rxn_list[i][1][0])
            for one_reactant in reactants_list: # Reactants here actually means the products of a reaction in a pathway
                for one_starting_cmpd in starting_cmpds:
                    simindex_list1.append(similarity_score(one_starting_cmpd,one_reactant,fp_type_global))
            final_reactant_score=max(simindex_list1)
            for one_starting_cmpd in starting_cmpds:
                simindex_list2.append(similarity_score(one_starting_cmpd,hash_a,fp_type_global))
            product_score=max(simindex_list2)
            # 2. Calculate probability based on the two scores generated:
            if product_score>=final_reactant_score:
                probability=1
            else:
                probability=math.exp((product_score-final_reactant_score)/T0_global/(XI_global**prod_level))
            # 3. Update the global bwd_probabilities_dict:
            if bwd_probabilities_dict.has_key(hash_a):
                if probability>bwd_probabilities_dict[hash_a]:
                    bwd_probabilities_dict[hash_a]=probability
            else:
                bwd_probabilities_dict[hash_a]=probability
                #bwd_level_dict[hash_a]=prod_level+1
            return round(probability,3)
        # ==================================================================================== #
        def bwd_probability_two_cmpds(starting_cmpds,reactant,prod):
            # 1. Use the following initiated objects to get similarity score for the reactant and product.
            simindex_list1=[] # All scores for the reactant (reactant tb predicted, 'products' in the pathways)
            simindex_list2=[] # All scores for the prod (prod tb expanded, 'reactants' in the pathways)
            for one_starting_cmpd in starting_cmpds:
                simindex_list1.append(similarity_score(one_starting_cmpd,reactant,fp_type_global))
                simindex_list2.append(similarity_score(one_starting_cmpd,prod,fp_type_global))
            reactant_score=max(simindex_list1)
            prod_score=max(simindex_list2)
            # 2. Calculate probability based on the two scores generated:
            if prod_score>=reactant_score:
                probability=1
            else:
                probability=math.exp((prod_score-reactant_score)/T0_global/(XI_global**prod_level))
            return probability
        # ==================================================================================== #
        def bwd_probability_select(cmpds_list,rxn_list,starting_cmpds,num_cmpds,bwd_level_dict):
            #print ("e^(-delta(sim_score)/(T0^prod_level))", ", prod_level=", prod_level)
            # Return a list of compounds according to probability rankings, in other words, similarity improvement
            probability_dict={}
            top_similarity_improvement=[] # Output list.
            # 1. Update the probability_dict and the global bwd_probabilities_dict directly:

            probability_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        probability=bwd_probability_two_cmpds(starting_cmpds,rxn[1][0],prod)
                        if (prod in bwd_probabilities_dict.keys() and probability > bwd_probabilities_dict[prod]) or (prod not in bwd_probabilities_dict.keys()):
                            bwd_probabilities_dict[prod]=probability
            return top_similarity_improvement, probability_dict

        # ==================================================================================== #
        def bwd_similarity_improvement_select_7(cmpds_list,rxn_list,starting_cmpds,num_cmpds,bwd_level_dict):
            # Used for [-7,-7,-7] alone.
            similarity_improvement_diff_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        similarity_improvement_diff=max([similarity_score(prod,one_starting_cmpd,fp_type_global) for one_starting_cmpd in starting_cmpds])- \
                                                        max([similarity_score(rxn[1][0],one_starting_cmpd,fp_type_global) for one_starting_cmpd in starting_cmpds])
                        if (prod in similarity_improvement_diff_dict.keys() and similarity_improvement_diff > similarity_improvement_diff_dict[prod]) or (prod not in similarity_improvement_diff_dict.keys()):
                            similarity_improvement_diff_dict[prod] = round(similarity_improvement_diff, 3)
            sorted_similarity_improvement_diff_list = sorted(similarity_improvement_diff_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print ("sorted_similarity_improvement_diff_list: ")
            #print (sorted_similarity_improvement_diff_list)
            top_similarity_improvement=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_similarity_improvement_diff_list) else len(sorted_similarity_improvement_diff_list)
            for j in range(num_cmpds):
                if sorted_similarity_improvement_diff_list[j][1] !=-1:
                    top_similarity_improvement.append(sorted_similarity_improvement_diff_list[j][0])
                else:
                    break
            return top_similarity_improvement
        # ==================================================================================== #
        def bwd_similarity_improvement_select_8(cmpds_list,rxn_list,starting_cmpds,num_cmpds,bwd_level_dict):
            # Used for [-8,-8,-8] alone.
            similarity_improvement_ratio_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        similarity_improvement_ratio=max([similarity_score(prod,one_starting_cmpd,fp_type_global) for one_starting_cmpd in starting_cmpds])/ \
                                                        max([similarity_score(rxn[1][0],one_starting_cmpd,fp_type_global) for one_starting_cmpd in starting_cmpds]+[0.00001])
                        if (prod in similarity_improvement_ratio_dict.keys() and similarity_improvement_ratio > similarity_improvement_ratio_dict[prod]) or (prod not in similarity_improvement_ratio_dict.keys()):
                            similarity_improvement_ratio_dict[prod]=similarity_improvement_ratio
            sorted_similarity_improvement_ratio_list = sorted(similarity_improvement_ratio_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print ("sorted_similarity_improvement_ratio_list: ")
            #print (sorted_similarity_improvement_ratio_list)
            top_similarity_improvement=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_similarity_improvement_ratio_list) else len(sorted_similarity_improvement_ratio_list)
            for j in range(num_cmpds):
                if sorted_similarity_improvement_ratio_list[j][1] !=0:
                    top_similarity_improvement.append(sorted_similarity_improvement_ratio_list[j][0])
                else:
                    break
            return top_similarity_improvement
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 0-0 Draw all reactants to be substitued and prepare for backward reaction search
        cmpds_need_expand=[] # A list of compounds input into reactor.bwd_apply_enzymes()
        for i in range(len(bwd_pathways_list[-1])):
            for j in range(len(bwd_pathways_list[-1][i])):
                if bwd_pathways_list[-1][i][j][0] == 0 - prod_level:
                    for k in range(len(bwd_pathways_list[-1][i][j][2])):
                        if bwd_pathways_list[-1][i][j][2][k] not in subs_smiles and bwd_pathways_list[-1][i][j][2][k] not in bkgd_cmpd_list:
                            cmpds_need_expand.append(bwd_pathways_list[-1][i][j][2][k])
        screened_cmpds_need_expand_list=[]
        for one_smiles in cmpds_need_expand:
            if one_smiles not in bkgd_cmpd_list:
                if one_smiles.count("C")<=max_C_num and \
                    one_smiles.count("O")<=max_O_num:
                    screened_cmpds_need_expand_list.append(one_smiles)

        #####----------Step 0-1 Generate new compounds and update reaction list
        # bwd_new_rxn_list : A list of compounds contains outputs of reactor.bwd_apply_enzymes()
        bwd_new_rxn_list=self.reactor.bwd_apply_enzymes_AP(screened_cmpds_need_expand_list,prod_level)
        #####----------Step 0-2 Update trfm_rxn_dict & trfm_rule_dict.  Modify bwd_new_rxn_list. (new version)
        
        modified_bwd_new_rxn_set=set([])
        for one_rxn in bwd_new_rxn_list:
            one_trfm=(one_rxn[1],one_rxn[0])
            if trfm_odict[one_trfm]==0:
                trfm_odict[one_trfm]=set([one_rxn[2],])
                trfm_id=len(trfm_odict)-1
            else:
                trfm_odict[one_trfm].add(one_rxn[2])
                trfm_id=list(trfm_odict.keys()).index(one_trfm)
            modified_bwd_new_rxn_set.add((one_rxn[0],one_rxn[1],trfm_id))
        bwd_new_rxn_list=list(modified_bwd_new_rxn_set)
        print ("Number of new reactions: " + str(len(bwd_new_rxn_list)) +". ", end="")

        #####----------Step 0-2 Update bwd_rxn_set.  Initialize bwd_current_smiles_set & bwd_current_accepted_smiles_set 
        bwd_current_smiles_set=set([])
        bwd_current_accepted_smiles_set=set([])
        for i in range(len(bwd_new_rxn_list)):
            if (bwd_new_rxn_list[i][0],bwd_new_rxn_list[i][1],) not in unfavorable_rxn_list:
                bwd_rxn_list.append(bwd_new_rxn_list[i])
                bwd_current_smiles_set=bwd_current_smiles_set.union(set(bwd_new_rxn_list[i][0]))

        #####----------Step 0-3 Draw all products from predicted reactions and get ready for SI and SA selection
        #print "number of compounds found on searching this level (prod side): " + str(len(bwd_current_smiles_set))
        for one_smiles in bwd_current_smiles_set:
            if one_smiles not in bwd_level_dict.keys():
                bwd_level_dict[one_smiles]=prod_level+1
        #####----------Step 0-4 Select (accept) a number of compounds based on SI and SA
        # -1. If last prod level, no need for screening or selecting compounds for next level 
        if prod_level==((max_levels_global)/2) -1:
            print ("last level, skip compound selection.")
            print ("="*100)
            bwd_current_accepted_smiles_set=bwd_current_accepted_smiles_set.union(bwd_current_smiles_set)
        else:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#
        # 0. Pre-screening: actually don't need to remove large compounds for backward search! So could have remove this step here!
            screened_smiles_list=[]
            if P0_global==[0,0,0]:
                for one_smiles in bwd_current_smiles_set:
                    one_smiles_expand=one_smiles
                    if one_smiles.find("I")!=-1:
                        one_smiles_expand=one_smiles.replace("(I)", "(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
                        one_smiles_expand=one_smiles_expand.replace("I", "(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
                        one_smiles_expand=unis(one_smiles_expand)

                    if one_smiles_expand in KEGG_nme_canonical_SMILES_dict.keys() and one_smiles_expand not in bkgd_cmpd_list:
                        #print ("KEGG cmpd selected, ", KEGG_nme_canonical_SMILES_dict[one_smiles_expand], one_smiles_expand)
                        bwd_current_accepted_smiles_set.add(one_smiles)

                    if one_smiles not in bkgd_cmpd_list:
                        if one_smiles.count("C")<=max_C_num and \
                            one_smiles.count("O")<=max_O_num:
                            screened_smiles_list.append(one_smiles)
            else:
                for one_smiles in bwd_current_smiles_set:
                        if one_smiles.count("C")<=max_C_num and \
                            one_smiles.count("O")<=max_O_num:
                            screened_smiles_list.append(one_smiles)

            #print "number of compounds screened on searching this level (subs side): " + str(len(screened_smiles_list))
            print ("Number of new compounds: " + str(len(screened_smiles_list)) + ".")
            if prod_level==0:
                for one_smiles in bwd_current_smiles_set:
                    bwd_current_accepted_smiles_set.add(one_smiles)
            else:
                # 0. Determine the number of compounds to accept based on the P0_global value
                if P0_global==[0,0,0]: # Select half number of compounds from three groups
                    SISAaccept=len(screened_smiles_list)
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                elif P0_global==[2,2,2]: # Select all compounds
                    SISAaccept=len(screened_smiles_list)
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[-8,-8,-8]: # Select half number of compounds with higher similarity improvement ratios #(used to be x0.5)
                    SISAaccept=len(screened_smiles_list)
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),0,0]
                elif P0_global==[-9,-9,-9]: # Select half number of compounds from three groups #(used to be x0.167)
                    SISAaccept=len(screened_smiles_list)
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                else: # Select certain numbers of compounds from the three groups (User-specified)
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[min(int(SISAaccept*bin_adj_global[0]),P0_global[0]),\
                                            min(int(SISAaccept*bin_adj_global[1]),P0_global[1]),\
                                            min(int(SISAaccept*bin_adj_global[2]),P0_global[2])]

                if P0_global==[0,0,0]: 
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=bwd_similarity_improvement_select_7(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print ("Similarity improvement: ", len(top_similarity_improvement))
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)

                    # Probability Select (select according to the probability dict)
                    # The bwd_probabilities_dict used here contains the probabilities of accepting compounds (ALL generated).
                    (_,_)=bwd_probability_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    Probabilities_select=[]
                    random_order_list=randomList(list(bwd_probabilities_dict.keys()))
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if bwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            count_x=count_x-1
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print ("Probabilities Select: ", len(Probabilities_select))
                    # Similarity Select (based on similarity scores)
                    Similarities_select=[]
                    taofactors_dict=similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print ("Similarities Select: ", len(Similarities_select))

                elif P0_global==[2,2,2]: # Select all compounds
                    bwd_current_accepted_smiles_set=screened_smiles_list
                    print ("Accepted All:", len(screened_smiles_list))
                elif P0_global==[-8,-8,-8]:
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=bwd_similarity_improvement_select_8(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print ("Similarity improvement: ", len(top_similarity_improvement))
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)

                elif P0_global==[-9,-9,-9]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=bwd_similarity_improvement_select_7(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print ("Similarity improvement: ", len(top_similarity_improvement))
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)


                    # Probability Select (select according to the probability dict)
                    # The bwd_probabilities_dict used here contains the probabilities of accepting compounds (ALL generated).
                    (top_similarity_improvement, probability_dict)=bwd_probability_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    Probabilities_select=[]
                    random_order_list=randomList(list(bwd_probabilities_dict.keys()))
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if bwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            count_x=count_x-1
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print ("Probabilities Select: ", len(Probabilities_select))
                    # Similarity Select (based on similarity scores)
                    Similarities_select=[]
                    taofactors_dict=similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print ("Similarities Select: ", len(Similarities_select))

                else:
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement, probability_dict)=bwd_probability_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print ("Similarity improvement: ", len(top_similarity_improvement))
                    Similarity_improvement=[]
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                        Similarity_improvement.append(one_smiles)

                    # Probability Select (select according to the probability dict)
                    # The bwd_probabilities_dict used here contains the probabilities of accepting compounds (ALL generated).
                    Probabilities_select=[]
                    random_order_list=randomList(list(bwd_probabilities_dict.keys()))
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if bwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            count_x=count_x-1
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print ("Probabilities Select: ", len(Probabilities_select))
                    # Similarity Select (based on similarity scores)
                    Similarities_select=[]
                    taofactors_dict=similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print ("Similarities Select: ", len(Similarities_select))
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 1 Take the newest group of incomplete pathways from bwd_pathways_list (in last [])
        tb_expanded_list=bwd_pathways_list[-1] # get the last list in bwd_pathways_list, that is the list to be expanded upon
        bwd_pathways_list.append([]) # Prepare for adding a new group of incomplete pathways of current level

        #####----------Step 2 Expand 2 pathway lists for backward direction search(, FIRST add reactions with compounds found on THIS LEVEL).
        if bwd_pathways_list[0][0][0][0]==0:
            #####----------Step 2-0 
            # At the first backward search step, initiate the bwd_pathways_list with the first backward level search result
            bwd_pathways_list=[[]] # 
            for n in range(len(bwd_new_rxn_list)):
                for o in range(len(bwd_new_rxn_list[n][0])):
                    if bwd_new_rxn_list[n][0][o] in bwd_current_accepted_smiles_set:
                        one_pathway_list=[]
                        one_pathway_list.append((-1,bwd_new_rxn_list[n][1],bwd_new_rxn_list[n][0],bwd_new_rxn_list[n][2]))
                        bwd_pathways_list[-1].append(one_pathway_list)
                        del one_pathway_list
                        break

        else:
            #####----------Step 2-1 For both full-length and non-full-length pathways, add reactions with compounds found on THIS LEVEL.
            # After the first backward search step, bwd_pathways_list now has been initiated and are ready for backward substitution.
            print("Start backward substitution, number of compounds to be expanded: ", len(tb_expanded_list))

            for n in tqdm(range(len(tb_expanded_list))):
                # 0. Obtain all compounds along this one certain pathway, store in reactants_found_already.
                reactants_found_already=set([])
                for o in range(len(tb_expanded_list[n])):
                    for p in range(len(tb_expanded_list[n][o])):
                        if p==1 or p==2:
                            reactants_found_already = reactants_found_already | set(tb_expanded_list[n][o][p])

                # 1. Obtain all compounds to be substituted by predicted reactants(, duplicates and starting compounds removed).
                reactantstobeexpanded=[]
                for o in range(len(tb_expanded_list[n])):
                    if tb_expanded_list[n][o][0] == 0 - prod_level:
                        for p in range(len(tb_expanded_list[n][o][2])):
                            if (tb_expanded_list[n][o][2][p] not in reactantstobeexpanded) and tb_expanded_list[n][o][2][p] not in subs_smiles:
                                reactantstobeexpanded.append(tb_expanded_list[n][o][2][p])

                # 2. Normally reactantstobeexpanded shall not be empty list.
                if len(reactantstobeexpanded)==0:
                    print ("completed pathway found in incomplete pathway list")
                    continue

                # 3. For each compounds to be substituted, get all reactions (from this step of search) that can be used to substitute this compound.
                next_level_reactions=[]
                for q in range(len(reactantstobeexpanded)):
                    next_level_reactions.append([]) # Each reactants has a corresponding list in next_level_reactions, which contains all reactions that can be used to substitute this reactant.
                    for r in range(len(bwd_new_rxn_list)):
                        for s in range(len(bwd_new_rxn_list[r][0])):
                            if bwd_new_rxn_list[r][0][s] in bwd_current_accepted_smiles_set and \
                                reactantstobeexpanded[q] in bwd_new_rxn_list[r][1] and \
                                bwd_new_rxn_list[r][0][s] not in reactants_found_already and \
                                bwd_level_dict[bwd_new_rxn_list[r][0][s]] == prod_level+1:# ! ! ! ! ! ! ! ! : This line makes sure only compounds found on this level is substituted.
                                next_level_reactions[q].append((0 - prod_level -1, bwd_new_rxn_list[r][1],bwd_new_rxn_list[r][0],bwd_new_rxn_list[r][2]))
                                break # Once a reaction is accepted, the other reactants (, or backwardly products) dont need to be checked.
                

                next_level_reactions=[]
                for q in range(len(reactantstobeexpanded)):
                    next_level_reactions.append([]) # Each reactants has a corresponding list in next_level_reactions, which contains all reactions that can be used to substitute this reactant.
                    for r in range(len(bwd_new_rxn_list)):
                        for s in range(len(bwd_new_rxn_list[r][0])):
                            if bwd_new_rxn_list[r][0][s] in bwd_current_accepted_smiles_set and \
                                reactantstobeexpanded[q] in bwd_new_rxn_list[r][1] and \
                                bwd_new_rxn_list[r][0][s] not in reactants_found_already and \
                                bwd_level_dict[bwd_new_rxn_list[r][0][s]] == prod_level+1:# ! ! ! ! ! ! ! ! : This line makes sure only compounds found on this level is substituted.
                                next_level_reactions[q].append((0 - prod_level -1, bwd_new_rxn_list[r][1],bwd_new_rxn_list[r][0],bwd_new_rxn_list[r][2]))
                                break # Once a reaction is accepted, the other reactants (, or backwardly products) dont need to be checked.


                # 4. Use cart_prod to get all combinations of those reactions that is all possibilities how compounds can be substituted
                next_level_reactions_set=cart_prod(next_level_reactions)
                # 5. Expand one pathway and check if it is completed
                for u in range(len(next_level_reactions_set)):
                    one_pathway_list=deepcopy(tb_expanded_list[n])
                    for v in range(len(next_level_reactions_set[u])):
                        if next_level_reactions_set[u][v] not in one_pathway_list: # Prevent adding two identical reactions into one pathway (, doesnt do anything to avoid bad pathways)
                            one_pathway_list.append(next_level_reactions_set[u][v])
                    if one_pathway_list not in bwd_pathways_list[prod_level] and one_pathway_list not in cmplt_pathways_list:
                        # 5.1 Now check if the pathway is already completed.
                        reactantswillbeexpanded=[] # Predict all compounds to be substituted at next level, if none, pathway is completed.
                        for w in range(len(one_pathway_list)):
                            if one_pathway_list[w][0] == 0 - prod_level - 1:
                                for x in range(len(one_pathway_list[w][2])):
                                    if (one_pathway_list[w][2][x] not in reactantswillbeexpanded) and one_pathway_list[w][2][x] not in subs_smiles:
                                        reactantswillbeexpanded.append(one_pathway_list[w][2][x])
                        # 5.2 If a complete pathway, add to another list not the bwd_pahways_list
                        if reactantswillbeexpanded!=[]:
                            bwd_pathways_list[prod_level].append(one_pathway_list)
                        else:
                            cmplt_pathways_list.append(one_pathway_list)
            

        # ! ! ! ! ! ! ! ! : Step 3 has almost been rewritten, so may have mistakes unfound.
        #####----------Step 3 Expand 2 pathway lists for backward direction search (, NOW add reactions with compounds found on PREVIOUS LEVELS).

        for z in range(prod_level):
            #print ("prod_level", z)
            cmpd_level=z+1 # cmpd_level is the level of the compounds to be substituted.
            #####----------Step 3-0 Compounds FOUND at FIRST step get DROPPED and then ACCEPTED at this current level
            if cmpd_level==1: # First deal with compounds FOUND at FIRST step (, this means another new pathway needs to be initiated).
                for n in range(len(bwd_rxn_list)):
                    for o in range(len(bwd_rxn_list[n][0])):
                        if (bwd_rxn_list[n][0][o] in bwd_current_accepted_smiles_set and \
                            bwd_level_dict[bwd_rxn_list[n][0][o]]==1 and \
                            bwd_rxn_list[n][1][0]==target_compound): # ! ! ! ! ! ! ! ! : This line makes sure only compounds found on FIRST level is substituted.
                            one_pathway_list=[]
                            one_pathway_list.append((-1-prod_level,bwd_rxn_list[n][1],bwd_rxn_list[n][0],bwd_rxn_list[n][2]))
                            bwd_pathways_list[-1].append(one_pathway_list)
                            break
                continue # Use a continue here so that dont need a long else part below.
            #print("done step #1.-2")
            #####----------Step 3-1 Compounds FOUND AFTER 1st step and before current step get DROPPED and then ACCEPTED at this current level
            # Prepare for expanding the bwd_pathways_list:
            tb_expanded_list=[]
            for m in range(len(bwd_pathways_list[-1-prod_level+cmpd_level-2])):
                # bwd_pathways_list contains all incomplete pathways histories at each level, stored in a different list.
                # Need to go to one earlier incomplete pathways list, so that subsitution can be made to add compounds found earlier but just accepted.
                # Need to adjust the earlier incomplete pathways to make substitution.
                adjusted_earlier_pathway=[]
                for n in range(len(bwd_pathways_list[-1-prod_level+cmpd_level-2][m])):
                    # adjust the reaction level index to prepare for substitution (use a,b below):
                    a=list(bwd_pathways_list[-1-prod_level+cmpd_level-2][m][n])
                    a[0]=a[0]-prod_level-1+cmpd_level
                    b=tuple(a)
                    adjusted_earlier_pathway.append(b)
                tb_expanded_list.append(adjusted_earlier_pathway)
            #print("done step #1.-1")

            #print "len(tb_expanded_list)", len(tb_expanded_list)
            for n in range(len(tb_expanded_list)):
                #print n
                # 0. Obtain all compounds along this one certain pathway, store in reactants_found_already.
                reactants_found_already=set([])
                for o in range(len(tb_expanded_list[n])):
                    for p in range(len(tb_expanded_list[n][o])):
                        if p==1 or p==2:
                            reactants_found_already = reactants_found_already | set(tb_expanded_list[n][o][p])

                # 1. Obtain all compounds to be substituted by predicted reactants(, duplicates and starting compounds removed).
                reactantstobeexpanded=[]
                for o in range(len(tb_expanded_list[n])):
                    if tb_expanded_list[n][o][0] == 0 - prod_level:
                        for p in range(len(tb_expanded_list[n][o][2])):
                            if (tb_expanded_list[n][o][2][p] not in reactantstobeexpanded) and tb_expanded_list[n][o][2][p] not in subs_smiles:
                                reactantstobeexpanded.append(tb_expanded_list[n][o][2][p])
                # 2. Normally reactantstobeexpanded shall not be empty list.
                if len(reactantstobeexpanded)==0:
                    print ("completed pathway found in incomplete pathway list2: ")
                    #print tb_expanded_list[n]
                    continue
                # 3. For each compounds to be substituted, get all reactions (from this step of search) that can be used to substitute this compound.
                next_level_reactions=[]
                for q in range(len(reactantstobeexpanded)):
                    next_level_reactions.append([]) # Each reactants has a corresponding list in next_level_reactions, which contains all reactions that can be used to substitute this reactant.
                    for r in range(len(bwd_rxn_list)):
                        for s in range(len(bwd_rxn_list[r][0])):
                            if bwd_rxn_list[r][0][s] in bwd_current_accepted_smiles_set and \
                                reactantstobeexpanded[q] in bwd_rxn_list[r][1] and \
                                bwd_rxn_list[r][0][s] not in reactants_found_already and \
                                bwd_level_dict[bwd_rxn_list[r][0][s]] == cmpd_level:# ! ! ! ! ! ! ! ! : This line makes sure only compounds found on that earlier level is substituted.
                                next_level_reactions[q].append((0 - prod_level -1, bwd_rxn_list[r][1],bwd_rxn_list[r][0],bwd_rxn_list[r][2]))
                                break # Once a reaction is accepted, the other reactants (, or backwardly products) dont need to be checked.
                # 4. Use cart_prod to get all combinations of those reactions that is all possibilities how compounds can be substituted
                next_level_reactions_set=cart_prod(next_level_reactions)
                # 5. Expand one pathway and check if it is completed
                for u in range(len(next_level_reactions_set)):
                    one_pathway_list=deepcopy(tb_expanded_list[n])
                    for v in range(len(next_level_reactions_set[u])):
                        if next_level_reactions_set[u][v] not in one_pathway_list:
                            one_pathway_list.append(next_level_reactions_set[u][v])
                    if one_pathway_list not in bwd_pathways_list[prod_level] and one_pathway_list not in cmplt_pathways_list:
                        # 5.1 Now check if the pathway is already completed.
                        reactantswillbeexpanded=[] # Predict all compounds to be substituted at next level, if none, pathway is completed.
                        for w in range(len(one_pathway_list)):
                            if one_pathway_list[w][0] == 0 - prod_level - 1:
                                for x in range(len(one_pathway_list[w][2])):
                                    if (one_pathway_list[w][2][x] not in reactantswillbeexpanded) and one_pathway_list[w][2][x] not in subs_smiles:
                                        reactantswillbeexpanded.append(one_pathway_list[w][2][x])
                        # 5.2 If a complete pathway, add to another list not the bwd_pahways_list
                        if reactantswillbeexpanded!=[]:
                            bwd_pathways_list[prod_level].append(one_pathway_list)
                        else:
                            cmplt_pathways_list.append(one_pathway_list)

        return bwd_pathways_list

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def pathway_searching(self, starting_cmpds_smiles_list, target_cmpd_smiles_list, max_levels, max_value):
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 0-0 Initialization of variables used in forward search
        # The following 3 objects are used (expanded at each level) in forward direction search. 
        # fwd_smiles_set                  : { 'smiles_A' , 'smiles_B' , 'smiles_C' , ...} 
        # fwd_rxn_set                     : { (('smiles_B','smiles_C'),('smiles_A'),'trfm_id') , () , () , ...} 
        # fwd_smiles_set_level_list       : A list of sets, the [i]th list contains the compounds first found on level i, the levels of cmpds will be used to reconstruct the pathways later.
        fwd_smiles_set=set(starting_cmpds_smiles_list)
        fwd_rxn_set=set([])
        fwd_smiles_set_level_list=[]
        for i in range( int((max_levels+1)/2+1) ):
            fwd_smiles_set_level_list.append(set([]))
        fwd_smiles_set_level_list[0]=deepcopy(fwd_smiles_set)
        #####----------Step 0-1 Initialization of variables used in forward search
        # The following 2 objects are used (expanded at each level) in both forward and backward direction search.
        # trfm_odict
        trfm_odict=OrderedCounter([])
        # ALL_RXN_SET
        ALL_RXN_SET=set([])
        #####----------Step 0-2 Initialization of variables used in backward search
        # The following 4 objects are used (expanded at each level) in backward direction search.
        # target_pathways_list    : shown in the long comments below
            #                [
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme') ],
            #                        [ ], [ ], ... [ ]
            #                    ],
            #
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), ( ), ... ( )  ], [ ], ... [ ],
            #                        [ (-2, ('target_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme') ], [ ], ... [ ]
            #                    ],
            #
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'),
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), ( ), ... ( ), 
            #                                (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), 
            #                                (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), ( ), ... ( ) ], [ ], ... [ ],
            #                        [ (-2, ('target_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'),
            #                            (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), 
            #                            (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), ( ), ... ( ) ], [ ], ... [ ],
            #                        [ (-3, ('target_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme') ], [ ], ... [ ]
            #                    ],
            #                    ......
            #                ]
        # complete_pathways_list  : 
        # backward_reactions_list : [ (('hash_B','hash_C'),('hash_A'),'enzyme_id') , () , () , ...]
        # backward_level_dict     : {'cmpd':lv, 'cmpd':lv, 'cmpd':lv, ...}
        # CAUTION: Here assume that there is only one target product!
        target_pathways_list=[[[(0,(target_cmpd_smiles_list[0],),(target_cmpd_smiles_list[0],),'trfm_id'),]]]
        complete_pathways_list=[]
        backward_reactions_list=[]
        backward_level_dict=dict([])
        #####----------Step 0-3 Global variables: 
        # Initialize global vars, target_smiles is a str, subs_smiles is a list of strs
        global searching_level; searching_level=0
        global max_levels_global; max_levels_global=max_levels
        global subs_level; subs_level=0
        global prod_level; prod_level=0
        global subs_smiles; subs_smiles=starting_cmpds_smiles_list
        global target_smiles; target_smiles=target_cmpd_smiles_list
        global target_compound; target_compound=target_cmpd_smiles_list[0]
        global target_compound_level; target_compound_level=-1 # target compound FIRST found on this level
        global found_on_subs_side; found_on_subs_side=False
        global KEGG_nme_canonical_SMILES_dict; KEGG_nme_canonical_SMILES_dict=self.KEGG_nme_canonical_SMILES_dict()
        global fwd_probabilities_dict; fwd_probabilities_dict={}
        global bwd_probabilities_dict; bwd_probabilities_dict={}
        global max_C_num; max_C_num=max(2,subs_smiles[0].count("C")+max_value[0],target_smiles[0].count("C")+max_value[1])
        global max_O_num; max_O_num=max(2,subs_smiles[0].count("O")+max_value[2],target_smiles[0].count("O")+max_value[3])
        print ("max_O_num",max_O_num)
        print ("max_C_num",max_C_num)
        #####----------Step 0-4 System Output: 
        # Print the beginning and targeting compounds
        print ("Substrate:")
        for one_smiles in starting_cmpds_smiles_list:
            print (one_smiles)
        print ("Products:")
        for one_smiles in target_cmpd_smiles_list:
            print (one_smiles)
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 1 Search
        # Expand 3 + 2 + 4 objects:
        # fwd_smiles_set
        # fwd_rxn_set
        # fwd_smiles_set_level_list
        #####
        # trfm_odict
        #####
        # target_pathways_list
        # complete_pathways_list
        # backward_reactions_list
        # backward_level_dict
        while (subs_level + prod_level < max_levels):
            if subs_level<=prod_level:
                # Forward direction search and expand
                print ("\n"+">"*25 + " Level #%d, searching on the subs side " % (searching_level + 1) + "<"*25 )
                self.expand_reaction_tree(fwd_smiles_set,fwd_rxn_set,fwd_smiles_set_level_list,trfm_odict,ALL_RXN_SET)

                searching_level += 1
                subs_level += 1
                # Print 
                print ("Compounds currently accepted: " , len(fwd_smiles_set))
                # Check if target has been found:
                if target_cmpd_smiles_list[0] in fwd_smiles_set and found_on_subs_side==False:
                    print ("Target product first found on: forward search level num %d," % (subs_level))
                    found_on_subs_side=True
                    self.pathway_found_textart()
                if (found_on_subs_side==True and target_compound_level==-1):
                    target_compound_level=subs_level # target compound FIRST found on this level
            else:
                # Backward direction search and expand
                print ("\n"+">"*25 + " Level #%d, searching on the prod side " % (searching_level + 1) + "<"*25 )
                target_pathways_list=self.bwd_expand_reaction_tree(target_pathways_list,complete_pathways_list,backward_reactions_list,backward_level_dict,trfm_odict,ALL_RXN_SET)
                searching_level += 1
                prod_level += 1
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #print ("Bad Rules")
        #print trfm_odict[list(trfm_odict.keys())[297]]   
        
        #####----------Step 2 Modify the bwd search results (target_pathways_list and complete_pathways_list)
        # 1. Modify target_pathways_list: 
        adj_target_pathways_list=[]
        for i in range(prod_level):
            adj_target_pathways_list.append([])
        for i in range(len(target_pathways_list)):
            for j in range(len(target_pathways_list[i])):
                if target_pathways_list[i][j][0][0]!=-1:
                    difference_value=-1-target_pathways_list[i][j][0][0]
                    adj_pathway=[]
                    for k in range(len(target_pathways_list[i][j])):
                        adj_reaction=list(target_pathways_list[i][j][k])
                        adj_reaction[0]=adj_reaction[0]+difference_value
                        adj_pathway.append(tuple(adj_reaction))
                    adj_target_pathways_list[adj_pathway[-1][0]*(-1)-1].append(adj_pathway)
                else:
                    adj_target_pathways_list[i].append(target_pathways_list[i][j])
        target_pathways_list=adj_target_pathways_list

        # 2. Modify complete_pathways_list:
        adj_complete_pathways_list=[]
        for i in range(len(complete_pathways_list)):
            if complete_pathways_list[i][0][0]!=-1:
                print ("Pathway found by iteration process") # ? ? ? ? ? ? ? ? : 
                difference_value=-1-complete_pathways_list[i][0][0]
                adj_pathway=[]
                for k in range(len(complete_pathways_list[i])):
                    adj_reaction=list(complete_pathways_list[i][k])
                    adj_reaction[0]=adj_reaction[0]+difference_value
                    adj_pathway.append(tuple(adj_reaction))
                adj_complete_pathways_list.append(adj_pathway)
            else:
                adj_complete_pathways_list.append(complete_pathways_list[i])
        complete_pathways_list=adj_complete_pathways_list
        #print complete_pathways_list

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 3 Prepare everything for next step, constructing pathways
        # Rename the fwd_smiles_set_level_list
        smilesfirstacceptlevel=deepcopy(fwd_smiles_set_level_list) # ? ? ? ? ? ? ? ? : the list isnt used in the program.
        # Modify the compoundshashlevel (add compounds in lower levels to higher levels) for substitute pathway construction:
        for h in range(len(fwd_smiles_set_level_list)-1):
            fwd_smiles_set_level_list[h+1]=fwd_smiles_set_level_list[h+1].union(fwd_smiles_set_level_list[h])
        fwd_smiles_list_level_list=[]
        for one_set in fwd_smiles_set_level_list:
            fwd_smiles_list_level_list.append(list(one_set))
        #print 'fwd_smiles_list_level_list: '
        #print fwd_smiles_list_level_list

        # Return all lists
        return [complete_pathways_list, target_pathways_list, list(fwd_smiles_set), list(fwd_rxn_set), fwd_smiles_list_level_list, smilesfirstacceptlevel, trfm_odict]



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def network_construction(self, all_lists, pathway_name):
        # Information need:
        # cmplt_pathways_list, target_pathways_list, subs_side_compoundshash_list,reactions_list,fwd_smiles_list_level_list
        # will have to deal with three situations: 
        # 1. pathways already found on subs side searching
        # 2. pathways already found on searching prod side
        # 3. combine the searching results from both searching direction and check the existence of pathways
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        def subs_side_networks_generator(target_pathways_list, fwd_smiles_list_level_list, reactions_list, temp_target): 
            all_pathways_list=[]
            pathway_length=subs_level
            pathways_list=[]

            # give the very first reactions to the pathways_list and start to expand
            for j in range(len(reactions_list)):#1.1
                #print (temp_target in reactions_list[j][0]), (iftuplestrinlist(reactions_list[j][1],fwd_smiles_list_level_list[pathway_length-1]))
                if (temp_target in reactions_list[j][0]) and (iftuplestrinlist(reactions_list[j][1],fwd_smiles_list_level_list[pathway_length-1])):
                    pathways_list.append([(pathway_length,reactions_list[j][0],reactions_list[j][1],reactions_list[j][2]),])
            for k in range(pathway_length-1):#1.2
                l=pathway_length-1-k
                temp_pathways_list=[]
                for m in range(len(pathways_list)):#1.2.1
                    temp_pathways_list.append(pathways_list[m])

                for n in range(len(temp_pathways_list)):#1.2.2
                    pathways_list.remove(temp_pathways_list[n])
                    reactantstobeexpanded=[]
                    for o in range(len(temp_pathways_list[n])):#1.2.2.1
                        if temp_pathways_list[n][o][0]==l+1:
                            for p in range(len(temp_pathways_list[n][o][2])):#1.2.2.1.1
                                if temp_pathways_list[n][o][2][p] not in fwd_smiles_list_level_list[0]:
                                    reactantstobeexpanded.append(temp_pathways_list[n][o][2][p])
                    next_level_reactions=[]
                    if len(reactantstobeexpanded)>0:
                        for q in range(len(reactantstobeexpanded)):#1.2.2.2
                            next_level_reactions.append([])
                            for r in range(len(reactions_list)):#1.2.2.2.1
                                if (reactantstobeexpanded[q] in reactions_list[r][0]) and iftuplestrinlist(reactions_list[r][1],fwd_smiles_list_level_list[l-1]):
                                    next_level_reactions[q].append((l,reactions_list[r][0],reactions_list[r][1],reactions_list[r][2]))
                        next_level_reactions_set=cart_prod(next_level_reactions) ###
                        for u in range(len(next_level_reactions_set)):#1.2.2.4
                            one_pathway_list=deepcopy(temp_pathways_list[n])
                            for v in range(len(next_level_reactions_set[u])):#1.2.2.4.1
                                if next_level_reactions_set[u][v] not in one_pathway_list:
                                    one_pathway_list.append(next_level_reactions_set[u][v])
                            # one_pathway_list.append(("subs_side_pathway",)) # 
                            if one_pathway_list not in pathways_list:
                                pathways_list.append(one_pathway_list)
                    else:
                        one_pathway_list=temp_pathways_list[n]
                        length_shorten=one_pathway_list[-1][0]-1
                        for a in range(len(one_pathway_list)):
                            temp_thing=list(one_pathway_list[a])
                            temp_thing[0]=temp_thing[0]-length_shorten
                            one_pathway_list[a]=tuple(temp_thing)
                       
                        if one_pathway_list not in pathways_list:
                            pathways_list.insert(0,one_pathway_list)

            for r in range(len(pathways_list)):
                all_pathways_list.append(pathways_list[r])
            return all_pathways_list
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        def convert_bwd_pathways(cmplt_pathways_list,number_to_add):
            # Used to return a complete pathway found in backward searching alone (prod_side_pathways)
            # Or to return an incomplete pathway (prod side completed) for subs side construction
            result_list=[]
            if number_to_add==0:
                for i in range(len(cmplt_pathways_list)):
                    one_pathway_list=[]
                    number_tb_add=1-cmplt_pathways_list[i][-1][0]
                    for j in range(len(cmplt_pathways_list[i])):
                        one_pathway_list.append((cmplt_pathways_list[i][j][0]+number_tb_add,cmplt_pathways_list[i][j][1],cmplt_pathways_list[i][j][2],cmplt_pathways_list[i][j][3]))
                    result_list.append(one_pathway_list)
            else:
                for i in range(len(cmplt_pathways_list)):
                    one_pathway_list=[]
                    number_tb_add=number_to_add
                    for j in range(len(cmplt_pathways_list[i])):
                        one_pathway_list.append((cmplt_pathways_list[i][j][0]+number_tb_add,cmplt_pathways_list[i][j][1],cmplt_pathways_list[i][j][2],cmplt_pathways_list[i][j][3]))
                    result_list.append(one_pathway_list)
            return result_list

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #def let_the_pathway_form(target_pathways_list,reactions_list,fwd_smiles_list_level_list):
            #...
        
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        # Main part of pathway_construction():
        # Get the search results
        cmplt_pathways_list=all_lists[0]
        target_pathways_list=all_lists[1]
        subs_side_compoundshash_list=all_lists[2]
        reactions_list=all_lists[3]
        fwd_smiles_list_level_list=all_lists[4]
        smilesfirstacceptlevel=all_lists[5]
        final_pathways_list=[]
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        dir="../results/" + pathway_name + "/"
        # Bwd_network
        prod_side_networks=convert_bwd_pathways(target_pathways_list[prod_level-1],prod_level+1)
        print ("Outputting Bwd_network")
        #print prod_side_networks
        pickle_out1=open(dir+"networks_bwd_00","wb")
        pickle.dump(prod_side_networks, pickle_out1)
        pickle_out1.close()
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        # Fwd_network
        print ("Outputted Bwd_network")
        print ("Outputting Fwd_network")
        subs_side_networks=[]
        print (smilesfirstacceptlevel[0])
        print (len(smilesfirstacceptlevel[subs_level]))
        count_x=0
        for one_subs_side_temp_target in smilesfirstacceptlevel[subs_level]:
            count_x+=1
            #print count_x
            subs_side_networks_part=subs_side_networks_generator(target_pathways_list, fwd_smiles_list_level_list, reactions_list, one_subs_side_temp_target)
            #print subs_side_networks_part
            for i in subs_side_networks_part:
                if i not in subs_side_networks:
                    subs_side_networks.append(i)
                    #print i
        print (len(subs_side_networks))
        pickle_out1=open(dir+"networks_fwd_00","wb")
        pickle.dump(subs_side_networks, pickle_out1)
        pickle_out1.close()
        print ("Outputted Fwd_network")
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        return final_pathways_list

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def result_pathways_conversion(self, final_pathways_list):
        converted_pwy_list=[]
        CoA_smiles_dict=dict([])
        CoA_string="(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)"
        for one_pathway in final_pathways_list:
            converted_pwy=[]
            for one_reaction in one_pathway:
                converted_rctt=tuple([unis((one_rctt.replace("(I)",CoA_string)).replace("I", CoA_string)) for one_rctt in one_reaction[2]])
                converted_prod=tuple([unis((one_prod.replace("(I)",CoA_string)).replace("I", CoA_string)) for one_prod in one_reaction[1]])
                for one_rctt in one_reaction[2]:
                    if one_rctt.find("I")!=-1:
                        one_rctt=one_rctt.replace("I", "CoA")
                        one_rctt_smiles=one_rctt.replace("(CoA)", "(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
                        one_rctt_smiles=one_rctt_smiles.replace("CoA", "(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
                        CoA_smiles_dict[one_rctt_smiles]=one_rctt
                for one_prod in one_reaction[1]:
                    if one_prod.find("I")!=-1:
                        one_prod=one_prod.replace("I", "CoA")
                        one_prod_smiles=one_prod.replace("(CoA)", "(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
                        one_prod_smiles=one_prod_smiles.replace("CoA", "(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
                        CoA_smiles_dict[one_prod_smiles]=one_prod

                converted_rxn=tuple([one_reaction[0], converted_prod, converted_rctt, one_reaction[3]])
                converted_pwy.append(converted_rxn)
            converted_pwy_list.append(converted_pwy)
        return converted_pwy_list, CoA_smiles_dict

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
    def solve_pathway(self, original_subs_string, original_prod_string, pathway_name, pathway_prefix, max_levels, rxn_db, pruning_method, fp_type, max_value, bin_adj,plot_type):
        # Subs in the inputs is a sting telling the starting compounds, which is converted to starting_cmpds_smiles_list (a list containing all the graphs of starting compounds)
        # Prod in the inputs is a string telling the target compounds, which is converted to target_cmpd_smiles_list
        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 0-1 Process the inputs:
        # Inputs are SMILES
        original_subs_string=original_subs_string.replace('[CoA]', 'CoA')
        original_prod_string=original_prod_string.replace('[CoA]', 'CoA')
        
        starting_cmpds_smiles_list = []
        if original_subs_string.find('.') == -1:
            if original_subs_string.find('CoA') == -1:
                starting_cmpds_smiles_list = [ unis(original_subs_string) ]
            else:
                starting_cmpds_smiles_list = [ unis(original_subs_string.replace('CoA', 'I')), ]
        else:
            for one_smiles_string in original_subs_string.split('.'):
                if one_smiles_string.find('CoA') == -1:
                    starting_cmpds_smiles_list.append(unis(one_smiles_string))
                else:
                    starting_cmpds_smiles_list.append(unis(one_smiles_string.replace('CoA', 'I')))

        target_cmpd_smiles_list = []
        if original_prod_string.find('.') == -1:
            if original_prod_string.find('CoA') == -1:
                target_cmpd_smiles_list = [ unis(original_prod_string) ]
            else:
                target_cmpd_smiles_list = [ unis(original_prod_string.replace('CoA', 'I')), ]
        else:
            for one_smiles_string in original_prod_string.split('.'):
                if one_smiles_string.find('CoA') == -1:
                    target_cmpd_smiles_list.append(unis(one_smiles_string))
                else:
                    target_cmpd_smiles_list.append(unis(one_smiles_string.replace('CoA', 'I')))


        #####----------Step 0-3 Global variables:
        global P0_global; P0_global=pruning_method
        global fp_type_global; fp_type_global=fp_type
        global max_value_global; max_value_global=max_value
        global bin_adj_global; bin_adj_global=bin_adj

        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 1 pathway_searching():  
        # Call the pathway_searching function to search pathways
        print ("="*100 + '\nStep 1: pathway searching\n')
        time0 = time.time()
        all_lists=[]
        all_lists=self.pathway_searching(starting_cmpds_smiles_list, target_cmpd_smiles_list, max_levels, max_value)
        trfm_odict=all_lists[6]
        
        dir="../results/" + pathway_name + "/"

        time1 = time.time()-time0
        print ( "Time Elapsed", time.time()-time0)
        print ('\nStep 1 completed\n' + "="*100 + "\n")

        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 2 pathway_construction():  
        # Call the pathway_construction function to construct pathways
        print ("*"*100  +'\nStep 2: pathway construction\n')
        time0 = time.time()
        result_pathways=[]
        result_pathways=self.pathway_construction(all_lists)
        time2 = time.time()-time0
        print ( "Time Elapsed", time.time()-time0)
        print ('\nStep 2 completed\n' + "*"*100 + "\n")
        result_pathways, CoA_smiles_dict =self.result_pathways_conversion(result_pathways)
        len_result_pathways=len(result_pathways)
        print("Number of pathways found: ", len_result_pathways)
        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 3 network_construction():  
        '''
        # Call the network_construction function to construct networks
        print "="*80 +'\nStep 3: network construction\n' + "="*80
        time0 = time.time()
        result_networks=self.network_construction(all_lists,pathway_name)
        time2 = time.time()-time0
        print time.time()-time0
        print "="*80 +'\nStep 2 completed\n' + "="*80
        '''

        print ("="*100 +'\nStep 3: network construction\n')
        time0 = time.time()

        print("Skipped for examples.")
        pickle_out1=open(dir+"all_lists_00","wb")
        pickle.dump(all_lists, pickle_out1)
        pickle_out1.close()
        
        time3 = time.time()-time0
        print ( "Time Elapsed", time.time()-time0)
        print ('\nStep 3 completed\n' + "="*100 + "\n")

        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 4 wrtie text files: 
        #####---------- 0. Write pathways.txt
        # Output the pathways in a txt file:
        text_file = open("../results/" + pathway_name + "/pathways.txt", "w")
        shortest_len=max_levels
        for one_pwy in result_pathways:
            shortest_len = min(one_pwy[0][0],shortest_len)

        shortest_count=0
        for one_pwy in result_pathways:
            if one_pwy[0][0]==shortest_len:
                shortest_count+=1
            for one_rxn in one_pwy:
                text_file.write(str(one_rxn))
                text_file.write("\n")
            text_file.write("\n")
        text_file.close()
        #os.system('xdg-open "../results/" + pathway_name + "/pathways.txt"')
        #subprocess.Popen(["notepad.exe", "../results/" + pathway_name + "/pathways.txt"])

        #####---------- 1. Write trfm_odict.txt
        text_file = open("../results/" + pathway_name + "/trfm_odict.txt", "w")
        for one_index in range(len(trfm_odict.keys())):

            text_file.write(str(one_index))
            text_file.write(" , ")
            text_file.write(str(list(trfm_odict[list(trfm_odict.keys())[one_index]])))
            text_file.write(" , ")
            text_file.write(str(list(trfm_odict.keys())[one_index]))
            text_file.write("\n")
        text_file.close()


        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 5 Plot: 
        print ("*"*100 +'\nStep 5: pathway visulization\n')
        time0 = time.time()

        num_selected_pwys,zpathways,zpathways_selected=output_zpathways(result_pathways, KEGG_nme_canonical_SMILES_dict, trfm_odict, CoA_cmpd_list, pathway_name, CoA_smiles_dict)
        run_r_plot(zpathways, zpathways_selected, pathway_name, plot_type)

        time4 = time.time()-time0
        print ( "Time Elapsed", time.time()-time0)
        print ('\nStep 5 completed\n' + "*"*100+ "\n")



        #####---------- 3. Write search_log.txt
        write_search_log(pathway_name,fp_type_global,T0_global,XI_global,P0_global,max_levels,max_value,bin_adj,subs_smiles,target_smiles,time1,time2,len_result_pathways,shortest_len,shortest_count,num_selected_pwys)
        #######################################################################################################################################
        #######################################################################################################################################
        #####----------Step 6 update database:
        text_file = open("../results/" + "/search_log.txt", "a+")
        text_file.write(str(pathway_name)+" "+str(fp_type_global)+" "+str(max_levels)+" "+str(T0_global)+" "+str(XI_global)+" "+str(P0_global)\
                                    +" "+str(subs_smiles)+" "+str(target_smiles)+" "+str(max_value)\
                                    +" "+str(bin_adj)+" "+str(len_result_pathways)+" "+str(shortest_len)+" "+str(shortest_count)+" "+str(time1)+" "+str(time2))
        text_file.write("\n")

        
        # Connect to the database
        '''
        print "connect to mysql db"
        connection = pymysql.connect(host='localhost',
                                     user='xuzhiqin',
                                     password='xxxxxxxx',
                                     db='annealpath_20190501_for_volatility',
                                     cursorclass=pymysql.cursors.DictCursor)
        # Insert the search result into DB
        cursor=connection.cursor()
        try:
            with connection.cursor() as cursor:
                sql_query = "INSERT INTO `search_result` (pwy_name,sim_mtrc,max_len,ini_temp,cool_coef,dsim_sel,simu_sel,simi_sel,\
                                subs_sms,prod_sms,max_C_subs,max_C_prod,max_O_prod,max_O_subs,dsim_adj,simu_adj,simi_adj,rxn_time,pwy_time,pwy_num,shtst,shtst_num\
                                ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
                val_tb_inserted=(pathway_name,str(fp_type_global),max_levels,T0_global,XI_global,P0_global[0],P0_global[1],P0_global[2],\
                                    str(subs_smiles),str(target_smiles),max_value[0],max_value[1],max_value[2],max_value[3],\
                                    bin_adj[0],bin_adj[1],bin_adj[2],time1,time2,len_result_pathways,shortest_len,shortest_count)
                cursor.execute(sql_query,val_tb_inserted)
            # connection is not autocommit by default. So you must commit to save
            # your changes.
            connection.commit()
        finally:
            connection.close()
        '''
        #######################################################################################################################################
        #######################################################################################################################################
        print ("$"*100)
        if result_pathways != []:
            print ("Success !!!",)
            #print ("%d pathways found" % (len(result_pathways)))
        else:
            print ("Failure !!!")
        print ("$"*100)



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def test_3():

    PF=AnnealPath()
    #pickle_in1=open("../zz_test_savings/all_lists_001","rb")
    #all_lists=pickle.load(pickle_in1)
    #pickle_in1.close()

    global subs_level; subs_level=3
    global prod_level; prod_level=3

    #test_result=PF.network_construction(all_lists, "test_pathway_name")



if __name__ == '__main__':
    test_3()
    PF=AnnealPath()
    PF.pathway_found_textart()

