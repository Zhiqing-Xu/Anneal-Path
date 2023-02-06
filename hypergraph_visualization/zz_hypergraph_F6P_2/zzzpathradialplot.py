#!/usr/bin/python
import os, os.path
from sys import platform
if os.name == 'nt' or platform == 'win32':
    os.chdir(os.path.dirname(__file__))
import ast
import pickle
from copy import copy, deepcopy

def replace_n(str1, n, str2):
    letters = (
        str2 if i == n else char
        for i, char in enumerate(str1)
    )
    return ''.join(letters)

def readlines_until_empty_line(file, mark=""):
    count_x=0
    reactions = []
    while True:
        line = file.readline()
        
        if (line == "" and count_x==0):
            return None
        one_reaction = line.rstrip()
        if (line == "\n" and count_x != 0):
            return reactions
        reactions.append(one_reaction)
        count_x+=1

def read_pathway_txt(pathways_file):
    pathways_list=[]
    while True:
        one_pathway_text_list=readlines_until_empty_line(pathways_file, mark="")
        one_pathway_list=[]
        if one_pathway_text_list==None:
            break
        for one_reaction_str in one_pathway_text_list:
            one_reaction=ast.literal_eval(one_reaction_str)

            #if one_reaction[2]==('CH@3,C,CH2,CH2,CH2,CH2,NH2,NH2,OH,O~100001100100100000001100000001000000020000000',):
                #one_reaction=tuple([one_reaction[0],one_reaction[1],('CH,C,CH2,CH2,CH2,CH2,NH2,NH2,OH,O~100001100100100000001100000001000000020000000',),one_reaction[3]])
            
            one_pathway_list.append(one_reaction)
        pathways_list.append(one_pathway_list)
    return pathways_list

def find_new_tb_sub(one_pathway,tb_sub,current_lv,distance):

    tb_sub_copy=copy.copy(tb_sub)
    new_tb_sub=[]
    for reaction in one_pathway:
        for cmpd_tb_sub in tb_sub:
            if reaction[0]==current_lv-distance and cmpd_tb_sub in reaction[1] and len(reaction[1])==1:
                if cmpd_tb_sub in tb_sub_copy:
                    tb_sub_copy.remove(cmpd_tb_sub)
                    for rctt in reaction[2]:
                        new_tb_sub.append(rctt)
                    


    if tb_sub_copy!=[]:
        return []
    else:
        return new_tb_sub


def replace_node(node_name):
    if node_name=="O=C(O)C(=O)CC(O)C(O)C(O)COP(=O)(O)O":
        return "DAHP"
    if node_name=="O=P(O)(O)OCC1OC(O)(CO)C(O)C1O":
        return "F6P"
    if node_name=="C=C(OP(=O)(O)O)C(=O)O":
        return "PEP"
    if node_name=="O=CC(O)C(O)COP(=O)(O)O":
        return "E4P"
    if node_name=="O=C(O)C(CO)OP(=O)(O)O":
        return "2PG"
    if node_name=="O=C(O)C(O)COP(=O)(O)O":
        return "3PG"
    if node_name=="C(C(C=O)O)OP(=O)(O)O":
        return "GA3P"
    if node_name=="O=C(CO)C(O)C(O)COP(=O)(O)O":
        return "X5P"
    if node_name=="O=P(O)(O)OCC1OC(O)(COP(=O)(O)O)C(O)C1O":
        return "F16bP"
    if node_name=="O=C(O)C(O)COP(=O)(O)O":
        return "DHAP"
    return node_name


def replace_f6p(node_name):
    if node_name=="O=P(O)(O)OCC1OC(O)(CO)C(O)C1O":
        return "F6P"
    return node_name


def main2():
    #print pathways_list

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



    nodes_list=[] # actually 
    node_value_dict=dict([])
    link_set_list=[]
    link_set_s_t_list=[]

    '''
    pathways_file=open("pathways_F6P3.txt")
    pathways_list=read_pathway_txt(pathways_file)
    print pathways_list
    for one_pathway in pathways_list:
        for i in one_pathway[-1][2]:
            if replace_node(i) not in nodes_list:
                nodes_list.append(replace_node(i))
                node_value_dict[replace_node(i)]=0
        for one_rxn in one_pathway:
            for i in one_rxn[1]:
                if replace_node(i) not in nodes_list:
                    nodes_list.append(replace_node(i))
                    node_value_dict[replace_node(i)]=one_rxn[0]
                else:
                    node_value_dict[replace_node(i)]=min(one_rxn[0],node_value_dict[replace_node(i)])
            one_rxn_1=[]
            one_rxn_2=[]
            for i in one_rxn[1]:
                one_rxn_1.append(replace_node(i))
            for i in one_rxn[2]:
                one_rxn_2.append(replace_node(i))
            if set(one_rxn_1+one_rxn_2) not in link_set_list:
                link_set_list.append(set(one_rxn_1+one_rxn_2))
                link_set_s_t_list.append([set(one_rxn_1),set(one_rxn_2)])

    '''
    #print node_value_dict
    #print nodes_list
    #print link_set_list
    #print link_set_s_t_list
    #print "first step done"
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    
    pickle_in1=open("../zz_test_savings/networks_bwd_F6P","rb")
    networks_bwd=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open("../zz_test_savings/networks_fwd_F6P","rb")
    networks_fwd=pickle.load(pickle_in2)
    pickle_in2.close()
    print (len(networks_bwd))

    #--------------------------------------------------#
    main_nodes_list=deepcopy(nodes_list)
    count_x=0
    for one_pathway in networks_fwd:
        count_x+=1
        if count_x>150:
            break
        for i in one_pathway[-1][2]:
            if replace_node(i) not in nodes_list:
                nodes_list.append(replace_node(i))
            if replace_node(i) not in main_nodes_list:
                node_value_dict[replace_node(i)]=7
        for one_rxn in one_pathway:
            for i in one_rxn[1]:
                if replace_node(i) not in nodes_list:
                    nodes_list.append(replace_node(i))
                if replace_node(i) not in main_nodes_list:
                    node_value_dict[replace_node(i)]=7
            one_rxn_1=[]
            one_rxn_2=[]
            for i in one_rxn[1]:
                one_rxn_1.append(replace_node(i))
            for i in one_rxn[2]:
                one_rxn_2.append(replace_node(i))
            if set(one_rxn_1+one_rxn_2) not in link_set_list:
                link_set_list.append(set(one_rxn_1+one_rxn_2))
                link_set_s_t_list.append([set(one_rxn_1),set(one_rxn_2)])


    #--------------------------------------------------#
    #--------------------------------------------------#
    # I am writing this for the fwd network alone, since i only need to show linear paths in the fwd network for now.
    nodes_list_l=[] # actually 
    node_value_dict_l=dict([])
    link_set_list_l=[]
    link_set_s_t_list_l=[]


    linear_pwys=[]
    count_x=0
    for one_pathway in networks_fwd:
        count_x+=1
        if count_x>150:
            break
        pwy_len=one_pathway[0][0]
        if pwy_len==len(one_pathway):
            count_tt=0
            for i in range(pwy_len):
                if len(one_pathway[i][1])==1:
                    count_tt+=1
            if count_tt==pwy_len:
                linear_pwys.append(one_pathway)



    for one_pathway in linear_pwys:
        for i in one_pathway[-1][2]:
            if replace_node(i) not in nodes_list_l:
                nodes_list_l.append(replace_node(i))
            if replace_node(i) not in main_nodes_list:
                node_value_dict_l[replace_node(i)]=7
        for one_rxn in one_pathway:
            for i in one_rxn[1]:
                if replace_node(i) not in nodes_list_l:
                    nodes_list_l.append(replace_node(i))
                if replace_node(i) not in main_nodes_list:
                    node_value_dict_l[replace_node(i)]=7
            one_rxn_1=[]
            one_rxn_2=[]
            for i in one_rxn[1]:
                one_rxn_1.append(replace_node(i))
            for i in one_rxn[2]:
                one_rxn_2.append(replace_node(i))
            if set(one_rxn_1+one_rxn_2) not in link_set_list_l:
                link_set_list_l.append(set(one_rxn_1+one_rxn_2))
                link_set_s_t_list_l.append([set(one_rxn_1),set(one_rxn_2)])
    print ("link_set_s_t_list_l, "), link_set_s_t_list_l


    #--------------------------------------------------#
    #--------------------------------------------------#


    count_x=0
    for one_pathway in networks_bwd:
        count_x+=1
        if count_x>0:
            break
        for i in one_pathway[-1][2]:
            if replace_node(i) not in nodes_list:
                nodes_list.append(replace_node(i))
            if replace_node(i) not in main_nodes_list:
                node_value_dict[replace_node(i)]=8
        for one_rxn in one_pathway:
            for i in one_rxn[1]:
                if replace_node(i) not in nodes_list:
                    nodes_list.append(replace_node(i))
                if replace_node(i) not in main_nodes_list:
                    node_value_dict[replace_node(i)]=8
            one_rxn_1=[]
            one_rxn_2=[]
            for i in one_rxn[1]:
                one_rxn_1.append(replace_node(i))
            for i in one_rxn[2]:
                one_rxn_2.append(replace_node(i))
            if set(one_rxn_1+one_rxn_2) not in link_set_list:
                link_set_list.append(set(one_rxn_1+one_rxn_2))
                link_set_s_t_list.append([set(one_rxn_1),set(one_rxn_2)])
                

    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#
    #--------------------------------------------------##--------------------------------------------------#                
    
    #--------------------------------------------------#

    print (node_value_dict)

    print (nodes_list)
    print (link_set_s_t_list)
    print ("starting write R script")

    #--------------------------------------------------#
    # Reprocess the nodes
    nodes_C_dict=dict([])
    nodes_C_dict_rev=dict([])
    for i in range(len(nodes_list)):
        if nodes_list[i] not in ['F6P', 'DAHP', 'PEP', '2PG', 'X5P', 'E4P', '3PG', 'GA3P', 'DHAP', 'F16bP']:
            node_C_name=replace_n("C"+str(1000+i),1,"0")
        else:
            node_C_name=nodes_list[i] 
        nodes_C_dict[node_C_name]=nodes_list[i]
        nodes_C_dict_rev[nodes_list[i]]=node_C_name

    link_set_list_rewrite=[]
    for one_list in link_set_s_t_list:
        new_first_set=set([])
        for one_cmpd in one_list[0]:
            new_first_set.add(nodes_C_dict_rev[one_cmpd])
        new_second_set=set([])
        for one_cmpd in one_list[1]:
            new_second_set.add(nodes_C_dict_rev[one_cmpd])
        new_list_of_set=[new_first_set,new_second_set]
        if new_list_of_set not in link_set_list_rewrite:
            link_set_list_rewrite.append([new_first_set,new_second_set])


#####################################
    link_set_list_rewrite_l=[]
    for one_list in link_set_s_t_list_l:
        new_first_set=set([])
        for one_cmpd in one_list[0]:
            new_first_set.add(nodes_C_dict_rev[one_cmpd])
        new_second_set=set([])
        for one_cmpd in one_list[1]:
            new_second_set.add(nodes_C_dict_rev[one_cmpd])
        new_list_of_set=[new_first_set,new_second_set]
        if new_list_of_set not in link_set_list_rewrite_l:
            link_set_list_rewrite_l.append([new_first_set,new_second_set])
    print ("link_set_list_rewrite_l, "), link_set_list_rewrite_l
#####################################

    nodes_R_dict=dict([])
    nodes_R_dict_l=dict([])
    for i in range(len(link_set_list_rewrite)):
        node_R_name=replace_n("R"+str(1000+i),1,"0")
        nodes_R_dict[node_R_name]=link_set_list_rewrite[i]
        if link_set_list_rewrite[i] in link_set_list_rewrite_l:
            nodes_R_dict_l[node_R_name]=link_set_list_rewrite[i]


    print ("nodes_R_dict, "), nodes_R_dict
    print ("nodes_R_dict_l, "), nodes_R_dict_l
    print ("nodes_C_dict, "), nodes_C_dict





    C_nodes_num=len(nodes_C_dict)
    R_nodes_num=len(nodes_R_dict)

    #--------------------------------------------------#

    text_file = open("R_hypergraph_fwd_8.r", "w")
    text_file.write("#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#")
    text_file.write("\n")
    text_file.write("#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#")
    text_file.write("\n")
    text_file.write("library(hypergraph)")
    text_file.write("\n")
    text_file.write("library(hyperdraw)")
    text_file.write("\n")
    text_file.write("\n")
    text_file.write("\n")

    DHE_list=[]
    node_C_final_list=[]
    for i in range(R_nodes_num):
        DHE_name=replace_n("DHE"+str(1000+i),3,"0")
        node_R_name=replace_n("R"+str(1000+i),1,"0")
        edge_info=nodes_R_dict[node_R_name]
        str1="c("
        
        for one_C in edge_info[1]:
            str1=str1+"\""+one_C+"\""+","
            if one_C not in node_C_final_list:
                node_C_final_list.append(one_C)
        str1=str1[:-1]
        str1=str1+")"+", "+ "c("
        for one_C in edge_info[0]:
            str1=str1+"\""+one_C+"\""+","
            if one_C not in node_C_final_list:
                node_C_final_list.append(one_C)
        str1=str1[:-1]
        str1=str1+")"+", "+ "\""+node_R_name+"\""+")"
        
        text_file.write(DHE_name + " <- DirectedHyperedge(" + str1)
        text_file.write("\n")
        DHE_list.append(DHE_name)

    str3="Cnodes <- c("
    for one_C in node_C_final_list:
        str3=str3+"\""+one_C+"\""+","
    str3=str3[:-1]
    str3=str3+")"
    text_file.write(str3)
    text_file.write("\n")

    str4="Rnodes <- list("
    for one_DHE in DHE_list:
        str4=str4+one_DHE+","
    str4=str4[:-1]
    str4=str4+")"
    text_file.write(str4)
    text_file.write("\n")



    str2="hg <- Hypergraph(Cnodes,Rnodes)"
    text_file.write(str2)
    text_file.write("\n")

    text_file.write("hgbph <- graphBPH(hg)")
    text_file.write("\n")
    text_file.write("#plot(hgbph)")
    text_file.write("\n")
    text_file.write("\n")
    text_file.write("\n")
    text_file.write("\n")

    for one_edge in nodes_R_dict.keys():
        edge_info_l=nodes_R_dict[one_edge]
        if one_edge in nodes_R_dict_l.keys():
            str_1="edgeData(testrabph, c("
            str_1=str_1 + "\"" + list(edge_info_l[1])[0]+ "\"" + ", "
            str_1=str_1 + "\"" + one_edge + "\"" + "), c("
            str_1=str_1 + "\"" + one_edge + "\"" + ", "
            str_1=str_1 + "\"" + list(edge_info_l[0])[0]+ "\"" + "), "
            str_1=str_1+ "\"lwd\") <- c(\"3\", \"3\")"
            text_file.write(str_1)
            text_file.write("\n")
            str_1="edgeData(testrabph, c("
            str_1=str_1 + "\"" + list(edge_info_l[1])[0]+ "\"" + ", "
            str_1=str_1 + "\"" + one_edge + "\"" + "), c("
            str_1=str_1 + "\"" + one_edge + "\"" + ", "
            str_1=str_1 + "\"" + list(edge_info_l[0])[0]+ "\"" + "), "
            str_1=str_1+ "\"color\") <- \"blue\""
            text_file.write(str_1)
            text_file.write("\n")
            text_file.write("\n")
        else:
            str_1="edgeData(testrabph, c("
            str_1=str_1 + "\"" + list(edge_info_l[1])[0]+ "\"" + ", "
            str_1=str_1 + "\"" + one_edge + "\"" + "), c("
            str_1=str_1 + "\"" + one_edge + "\"" + ", "
            str_1=str_1 + "\"" + list(edge_info_l[0])[0]+ "\"" + "), "
            str_1=str_1+ "\"lwd\") <- c(\"5\", \"5\")"
            text_file.write(str_1)
            text_file.write("\n")
            str_1="edgeData(testrabph, c("
            str_1=str_1 + "\"" + list(edge_info_l[1])[0]+ "\"" + ", "
            str_1=str_1 + "\"" + one_edge + "\"" + "), c("
            str_1=str_1 + "\"" + one_edge + "\"" + ", "
            str_1=str_1 + "\"" + list(edge_info_l[0])[0]+ "\"" + "), "
            str_1=str_1+ "\"color\") <- \"firebrick\""
            text_file.write(str_1)
            text_file.write("\n")
            text_file.write("\n")

            if len(edge_info_l[1]) == 2:
                str_1="edgeData(testrabph, c("
                str_1=str_1 + "\"" + list(edge_info_l[1])[1]+ "\"" + ", "
                str_1=str_1 + "\"" + one_edge + "\"" + "), c("
                str_1=str_1 + "\"" + one_edge + "\"" + ", "
                str_1=str_1 + "\"" + list(edge_info_l[0])[0]+ "\"" + "), "
                str_1=str_1+ "\"lwd\") <- c(\"5\", \"5\")"
                text_file.write(str_1)
                text_file.write("\n")
                str_1="edgeData(testrabph, c("
                str_1=str_1 + "\"" + list(edge_info_l[1])[1]+ "\"" + ", "
                str_1=str_1 + "\"" + one_edge + "\"" + "), c("
                str_1=str_1 + "\"" + one_edge + "\"" + ", "
                str_1=str_1 + "\"" + list(edge_info_l[0])[0]+ "\"" + "), "
                str_1=str_1+ "\"color\") <- \"firebrick\""
                text_file.write(str_1)
                text_file.write("\n")
                text_file.write("\n")

            if len(edge_info_l[0]) == 2:
                str_1="edgeData(testrabph, c("
                str_1=str_1 + "\"" + list(edge_info_l[1])[0]+ "\"" + ", "
                str_1=str_1 + "\"" + one_edge + "\"" + "), c("
                str_1=str_1 + "\"" + one_edge + "\"" + ", "
                str_1=str_1 + "\"" + list(edge_info_l[0])[1]+ "\"" + "), "
                str_1=str_1+ "\"lwd\") <- c(\"5\", \"5\")"
                text_file.write(str_1)
                text_file.write("\n")
                str_1="edgeData(testrabph, c("
                str_1=str_1 + "\"" + list(edge_info_l[1])[0]+ "\"" + ", "
                str_1=str_1 + "\"" + one_edge + "\"" + "), c("
                str_1=str_1 + "\"" + one_edge + "\"" + ", "
                str_1=str_1 + "\"" + list(edge_info_l[0])[1]+ "\"" + "), "
                str_1=str_1+ "\"color\") <- \"firebrick\""
                text_file.write(str_1)
                text_file.write("\n")
                text_file.write("\n")

            if len(edge_info_l[1]) == 2 and len(edge_info_l[0]) == 2:
                str_1="edgeData(testrabph, c("
                str_1=str_1 + "\"" + list(edge_info_l[1])[1]+ "\"" + ", "
                str_1=str_1 + "\"" + one_edge + "\"" + "), c("
                str_1=str_1 + "\"" + one_edge + "\"" + ", "
                str_1=str_1 + "\"" + list(edge_info_l[0])[1]+ "\"" + "), "
                str_1=str_1+ "\"lwd\") <- c(\"5\", \"5\")"
                text_file.write(str_1)
                text_file.write("\n")
                str_1="edgeData(testrabph, c("
                str_1=str_1 + "\"" + list(edge_info_l[1])[1]+ "\"" + ", "
                str_1=str_1 + "\"" + one_edge + "\"" + "), c("
                str_1=str_1 + "\"" + one_edge + "\"" + ", "
                str_1=str_1 + "\"" + list(edge_info_l[0])[1]+ "\"" + "), "
                str_1=str_1+ "\"color\") <- \"firebrick\""
                text_file.write(str_1)
                text_file.write("\n")
                text_file.write("\n")





if __name__ == '__main__':
    #main1()
    main2()
