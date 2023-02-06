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
import os
import random
import pathlib
import subprocess
random.seed(42)
#--------------------------------------------------#
from AP_funcs import *
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def output_zpathways(result_pathways, KEGG_nme_canonical_SMILES_dict, trfm_odict, CoA_cmpd_list, pathway_name, CoA_smiles_dict):
    #####-----Step 01 add CoA_cmpd_list to KEGG_nme_canonical_SMILES_dict:
    for one_cmpd in CoA_cmpd_list:
        KEGG_nme_canonical_SMILES_dict[one_cmpd]=[one_cmpd, "None"]
    for one_CoA_smiles in CoA_smiles_dict.keys():
        if one_CoA_smiles in KEGG_nme_canonical_SMILES_dict.keys():
            KEGG_nme_canonical_SMILES_dict[one_CoA_smiles]=[CoA_smiles_dict[one_CoA_smiles], "CoA_cmpds"]
    #####-----Step 02 Modify result pathways, screen zpathways and zpathways_selected, and write text outputs.
    text_file = open("../results/" + pathway_name + "/zpathways.txt", "w")
    zpathways=[]
    zpathways_selected=[]
    for one_pwy in result_pathways:
        modified_pwy=[]
        bool_x=1
        for one_rxn in one_pwy:
            identified_prod_list=[]
            identified_rctt_list=[]
            identified_rules=[]
            modified_rxn_text=[]
            for one_prod in one_rxn[1]:
                if one_prod not in KEGG_nme_canonical_SMILES_dict.keys() : 
                    identified_prod_list.append(one_prod)
                    bool_x=0
                else:
                    identified_prod_list.append(KEGG_nme_canonical_SMILES_dict[one_prod][0])
            for one_rctt in one_rxn[2]:
                if one_rctt not in KEGG_nme_canonical_SMILES_dict.keys() : 
                    identified_rctt_list.append(one_rctt)
                    bool_x=0
                else:
                    identified_rctt_list.append(KEGG_nme_canonical_SMILES_dict[one_rctt][0])
            modified_rxn_text=tuple( [one_rxn[0], tuple(identified_prod_list), tuple(identified_rctt_list), tuple( trfm_odict[ list(trfm_odict.keys())[one_rxn[3]] ] ) ])
            text_file.write(str(tuple(modified_rxn_text)))
            text_file.write("\n")
            modified_pwy.append(modified_rxn_text)
        zpathways.append(modified_pwy)
        if bool_x==1:
            zpathways_selected.append(modified_pwy)
        text_file.write("\n")
    text_file.close()

    text_file = open("../results/" + pathway_name + "/zpwys_sel.txt", "w")
    for one_pwy in zpathways_selected:
        for one_rxn in one_pwy:
            text_file.write(str(one_rxn))
            text_file.write("\n")
        text_file.write("\n")
    text_file.close()
    return len(zpathways_selected),zpathways , zpathways_selected

#######################################################################################################################################
#######################################################################################################################################
def write_search_log(pathway_name,fp_type_global,T0_global,XI_global,P0_global,max_levels,max_value,bin_adj,subs_smiles,target_smiles,time1,time2,len_result_pathways,shortest_len,shortest_count,num_selected_pwys):
    text_file = open("../results/" + pathway_name + "/search_log.txt", "a+")
    text_file.write("pathway_name                    :" + pathway_name)
    text_file.write("\n")
    text_file.write("parameters                      :" + str(fp_type_global) + ',' + str(T0_global) + ',' + str(XI_global) + ','+ str(P0_global))
    text_file.write("\n")
    text_file.write("max_levels                      :" + str(max_levels))
    text_file.write("\n")
    text_file.write("max_values                      :" + str(max_value))
    text_file.write("\n")
    text_file.write("bin_adjs                        :" + str(bin_adj))
    text_file.write("\n")
    text_file.write("sub_smiles                      :" + str(subs_smiles))
    text_file.write("\n")
    text_file.write("target_smiles                   :" + str(target_smiles))
    text_file.write("\n")
    text_file.write("time of reaction computation    :" + str(time1))
    text_file.write("\n")
    text_file.write("time of pathway reconstruction  :" + str(time2))
    text_file.write("\n")
    text_file.write("number of pathways found        :" + str(len_result_pathways))
    text_file.write("\n")
    text_file.write("shortest length of pathways     :" + str(shortest_len))
    text_file.write("\n")
    text_file.write("number of shortest pathways     :" + str(shortest_count))
    text_file.write("\n")
    text_file.write("number of selected pathways     :" + str(num_selected_pwys))
    text_file.close()
    #os.system('xdg-open "../results/" + pathway_name + "/search_log.txt"')
    #subprocess.Popen(["notepad.exe", "../results/" + pathway_name + "/search_log.txt"])

#######################################################################################################################################
#######################################################################################################################################
def run_r_plot(zpathways, zpathways_selected, pathway_name, plot_type):

    if plot_type== "selected_unlimited" or plot_type == "all" or plot_type == "selected" or plot_type == "novel":
        if plot_type == "selected_unlimited":
            r_plot_pwy(zpathways_selected, pathway_name, "pwyr", 2)
            for one_pwy in zpathways_selected:
                r_plot_pwy([one_pwy], pathway_name, "pwyr", 1)

        elif plot_type == "selected":
            r_plot_pwy(zpathways_selected, pathway_name, "pwyr", 2)
            count_x=0
            for one_pwy in zpathways_selected:
                count_x+=1
                if count_x<=20:
                    r_plot_pwy([one_pwy], pathway_name, "pwyr", 1)

        elif plot_type == "novel":
            r_plot_pwy(zpathways, pathway_name, "pwyr", 2)
            count_x=0
            for one_pwy in zpathways:
                count_x+=1
                if count_x<=120:
                    r_plot_pwy([one_pwy], pathway_name, "pwyr", 1)
        else:
            r_plot_pwy(zpathways, pathway_name, "pwyr", 2)
            for one_pwy in zpathways:
                r_plot_pwy([one_pwy], pathway_name, "pwyr", 1)
    #print pathways_list
    return

#######################################################################################################################################
#######################################################################################################################################
def r_plot_pwy(pathways_list, pathway_name, r_plot_name, r_plot_type):
    #--------------------------------------------------#
    nodes_list=[] # actually 
    node_value_dict=dict([])
    link_set_list=[]
    link_name_list=[]
    link_set_s_t_list=[]
    #random.shuffle(pathways_list,)

    for one_pathway in pathways_list:
        for i in one_pathway[-1][2]:
            if i not in nodes_list:
                nodes_list.append(i)
                node_value_dict[i]=0
        for one_rxn in one_pathway:
            for i in one_rxn[1]:
                if i not in nodes_list:
                    nodes_list.append(i)
                    node_value_dict[i]=one_rxn[0]
                else:
                    node_value_dict[i]=min(one_rxn[0],node_value_dict[i])
            one_rxn_1=[]
            one_rxn_2=[]
            for i in one_rxn[1]:
                one_rxn_1.append(i)
            for i in one_rxn[2]:
                one_rxn_2.append(i)
            if [set(one_rxn_1),set(one_rxn_2)] not in link_set_s_t_list:
                link_set_list.append(set(one_rxn_1+one_rxn_2))
                link_set_s_t_list.append([set(one_rxn_1),set(one_rxn_2)])
                EC_id="EC"+one_rxn[3][0].split(".")[0]+"."+one_rxn[3][0].split(".")[1]+"."+one_rxn[3][0].split(".")[2]+"_0"
                while EC_id in link_name_list:
                    EC_id=EC_id.split("_")[0]+"_"+str(int(EC_id.split("_")[1])+1)
                link_name_list.append(EC_id)
    # ============================================================================================================================ #
    #print node_value_dict
    #print nodes_list
    #print link_set_list
    #print link_set_s_t_list
    #print link_name_list
    #print "first step done"
    # ============================================================================================================================ #
    #print "Starting write R script"
    # Reprocess the nodes
    nodes_C_dict=dict([])
    nodes_C_dict_rev=dict([])
    for i in range(len(nodes_list)):
        node_C_name=nodes_list[i] 
        nodes_C_dict[node_C_name]=nodes_list[i]
        nodes_C_dict_rev[nodes_list[i]]=node_C_name
    # ============================================================================================================================ #
    nodes_R_dict=dict([])
    for i in range(len(link_set_s_t_list)):
        node_R_name=link_name_list[i]
        nodes_R_dict[node_R_name]=link_set_s_t_list[i]
    # ============================================================================================================================ #
    #print nodes_R_dict
    #print nodes_C_dict
    C_nodes_num=len(nodes_C_dict)
    R_nodes_num=len(nodes_R_dict)
    # ============================================================================================================================ #
    str4="0"
    while os.path.isfile("../results/" + pathway_name + "/pwy_r/pathways/pathway" + str4 + ".png")==True:
        str4=str(int(str4)+1)
    #print ("ploting", str4)
    # ============================================================================================================================ #
    text_file = open("../results/" + pathway_name + "/pwy_r/" + r_plot_name + str4 + ".r", "w")
    text_file.write("if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n")
    text_file.write("    install.packages(\"BiocManager\")\n")
    text_file.write("BiocManager::install(\"hyperdraw\")\n")
    text_file.write("BiocManager::install(\"hypergraph\")\n")
    text_file.write("#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#\n")
    text_file.write("#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#\n")
    text_file.write("library(hypergraph)\n")
    text_file.write("library(hyperdraw)\n")
    text_file.write("\n")
    # ============================================================================================================================ #
    DHE_list=[]
    node_C_final_list=[]
    for i in range(R_nodes_num):
        DHE_name=replace_n("DHE"+str(1000+i),3,"0")
        node_R_name=link_name_list[i]
        edge_info=nodes_R_dict[node_R_name]
        str1="c("
        #--------------------------------------------------#
        for one_C in edge_info[1]:
            str1=str1+"\""+one_C+"\""+","
            if one_C not in node_C_final_list:
                node_C_final_list.append(one_C)
        str1=str1[:-1]
        str1=str1+")"+", "+ "c("
        #--------------------------------------------------#
        for one_C in edge_info[0]:
            str1=str1+"\""+one_C+"\""+","
            if one_C not in node_C_final_list:
                node_C_final_list.append(one_C)
        str1=str1[:-1]
        str1=str1+")"+", "+ "\""+node_R_name+"\""+")"
        text_file.write(DHE_name + " <- DirectedHyperedge(" + str1 + "\n")
        DHE_list.append(DHE_name)
    #--------------------------------------------------#
    str2="Cnodes <- c("
    for one_C in node_C_final_list:
        str2=str2+"\""+one_C+"\""+","
    str2=str2[:-1] if len(node_C_final_list)!=0 else str2+"\"C\""
    str2=str2+")"+"\n"
    text_file.write(str2)
    #--------------------------------------------------#
    str3="Rnodes <- list("
    for one_DHE in DHE_list:
        str3=str3+one_DHE+","
    str3=str3[:-1] if len(DHE_list)!=0 else str3+"\"E\""
    str3=str3+")"+"\n"
    text_file.write(str3)
    #--------------------------------------------------#
    text_file.write("hg <- Hypergraph(Cnodes,Rnodes)\n")
    text_file.write("hgbph <- graphBPH(hg)\n")
    text_file.write("testrabph <- graphLayout(hgbph)\n")
    text_file.write("edgeDataDefaults(testrabph, \"lwd\") <- 1\n")
    text_file.write("edgeDataDefaults(testrabph, \"color\") <- \"black\"\n")
    text_file.write("nodeDataDefaults(testrabph, \"margin\") <- 'unit(1, \"mm\")'\n")
    text_file.write("nodeDataDefaults(testrabph, \"shape\") <- \"box\"\n")
    text_file.write("nodeDataDefaults(testrabph, \"color\") <- \"black\"\n")
    text_file.write("#plot(testrabph)\n")
    #--------------------------------------------------#
    text_file.write("initial.options <- commandArgs(trailingOnly = FALSE)\n")
    text_file.write("file.arg.name <- \"--file=\"\n")
    text_file.write("script.name <- sub(file.arg.name, \"\", initial.options[grep(file.arg.name, initial.options)])\n")
    text_file.write("script.basename <- dirname(script.name)\n")
    text_file.write("other.name <- file.path(script.basename, \"\")\n")
    text_file.write("save_path<-paste(other.name," + "\"/pathways/pathway\"" + ",\"" + str4 + ".png\",sep=\"\")\n")
    #--------------------------------------------------#
    text_file.write("png(file = save_path, width = " + str(1920*r_plot_type) + ", height = " + str(1080*r_plot_type) + ", res = 300)\n")
    text_file.write("plot(testrabph)\n")
    text_file.write("dev.off()\n")
    #--------------------------------------------------#
    text_file.close()
    # ============================================================================================================================ #
    R_path=["C:\\R\\R-4.1.1\\bin\\Rscript.exe", "/usr/lib/R/bin/exec/R"][0]
    R_script_path=str(pathlib.Path(__file__).parent.parent.absolute())+"\\results\\"+pathway_name+"\\pwy_r\\"+ r_plot_name + str4 + ".r"
    subprocess_R_output = subprocess.check_output("\""+R_path+"\" \""+R_script_path+"\"")
    # ============================================================================================================================ #
    return


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
if (__name__ == '__main__'):
    #run_r_plot()

    print (os.path.dirname(os.path.abspath(__file__)))

    print (pathlib.Path(__file__).parent.parent.absolute())

    pathway_name = "test_pathway_name"
    print (os.path.isdir("../results/" + pathway_name + "/pwy_r/pathways/pathway" + str4 + ".png")==True)

