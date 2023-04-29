#!/usr/bin/python
import os, os.path
from sys import platform
if os.name == 'nt' or platform == 'win32':
    os.chdir(os.path.dirname(__file__))
import ast
import pickle
from copy import deepcopy
import random
rand_list=[]
for i in range(2000):
    random.seed(i)
    rand_list.append(random.random())
print rand_list

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
            
            one_pathway_list.append(one_reaction)
        pathways_list.append(one_pathway_list)
    return pathways_list

def find_new_tb_sub(one_pathway,tb_sub,current_lv,distance):

    tb_sub_copy=copy(tb_sub)
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

#--------------------------------------------------#
    nodes_list=[] # actually 
    node_value_dict=dict([])
    link_set_list=[]
    link_set_s_t_list=[]

    #print node_value_dict
    #print nodes_list
    #print link_set_list
    #print link_set_s_t_list
    #print "first step done"
#--------------------------------------------------#
    pwys=[]
    rct1=[]
    rct2=[]
    rct3=[]
    rcts=[rct1,rct2,rct3]
    lv0=["A","B"]
    Rules=["X","Y"]
    lv1=[]
    lv2=[]
    lv3=[]
    lv=[lv0,lv1,lv2,lv3]
    for i in range(len(lv)-1): # i = 0,1,2
        if i<2: 
            for j in range(len(lv[i])-1): # j = 0,1,2...
                for k in range(len(lv[i])-1-j): # k = 0,1,2,..
                    l=len(lv[i])-1-k
                    rctt_pair=lv[i][j]+lv[i][j-l]
                    for one_rule in Rules:
                        for prod_no in ["1","2"]:
                            one_prod=(rctt_pair+one_rule+prod_no)
                            lv[i+1].append(one_prod)
                        one_rct=(i+1 ,(rctt_pair+one_rule+"1",rctt_pair+one_rule+"2"),(lv[i][j],lv[i][j-l]),one_rule)
                        print one_rct
                        rcts[i].append(one_rct)
        else:
            count_xx=0
            for j in range(len(lv[i])-1): # j = 0,1,2...
                for k in range(len(lv[i])-1-j): # k = 0,1,2,..
                    l=len(lv[i])-1-k
                    rctt_pair=lv[i][j]+lv[i][j-l]
                    for one_rule in Rules:
                        for prod_no in ["1","2"]:
                            one_prod=(rctt_pair+one_rule+prod_no)
                            lv[i+1].append(one_prod)
                        one_rct=(i+1 ,(rctt_pair+one_rule+"1",rctt_pair+one_rule+"2"),(lv[i][j],lv[i][j-l]),one_rule)
                        count_xx+=1
                        if rand_list[count_xx]<0.05:
                            rcts[i].append(one_rct)
                            #print one_rct

    hypernet1=[[(1, ('ABX1', 'ABX2'), ('A', 'B'), 'X')],[(1, ('ABY1', 'ABY2'), ('A', 'B'), 'Y')]]


    hypernet2=[]
    for one_rct_2 in rct2:
        selected_rct1=[]
        tb_epd_cmpd=one_rct_2[2]
        for one_cmpd in tb_epd_cmpd:
            for one_rct_1 in rct1:
                if one_cmpd in one_rct_1[1]:
                    if one_rct_1 not in selected_rct1:
                        selected_rct1.append(one_rct_1)
        if len(selected_rct1)==2:
            hypernet2.append([one_rct_2,selected_rct1[0],selected_rct1[1]])
        if len(selected_rct1)==1:
            hypernet2.append([one_rct_2,selected_rct1[0]])       


    hypernet3=[]
    for one_rct_3 in rct3:
        selected_rct2=[]
        selected_rct1=[]
        tb_epd_cmpd2=one_rct_3[2]
        for one_cmpd in tb_epd_cmpd2:
            for one_rct_2 in rct2:
                if one_cmpd in one_rct_2[1]:
                    if one_rct_2 not in selected_rct2:
                        selected_rct2.append(one_rct_2)
        if len(selected_rct2)==2:
            one_pwy=[one_rct_3,selected_rct2[0],selected_rct2[1]]
            tb_epd_cmpd1=selected_rct2[0][2]+selected_rct2[1][2]
        if len(selected_rct2)==1:
            one_pwy=[one_rct_3,selected_rct2[0]]
            tb_epd_cmpd1=selected_rct2[0][2]

        for one_cmpd in tb_epd_cmpd1:
            for one_rct_1 in rct1:
                if one_cmpd in one_rct_1[1]:
                    if one_rct_1 not in selected_rct1:
                        selected_rct1.append(one_rct_1)
        if len(selected_rct1)==2:
            hypernet3.append(one_pwy+[selected_rct1[0],selected_rct1[1]])
        if len(selected_rct1)==1:
            hypernet3.append(one_pwy+[selected_rct1[0],])               



    print hypernet1
    print hypernet2
    print hypernet3


    pwys=[]
    rct1=[]
    rct2=[]
    rct3=[]
    rcts=[rct1,rct2,rct3]
    lv0=["A","B"]
    Rules=["X","Y"]
    lv1=[]
    lv2=[]
    lv3=[]
    lv=[lv0,lv1,lv2,lv3]

    
    for i in lv0:
        for j in Rules:
            for k in Rules:
                for l in Rules:
                    lv3.append(i+j+k+l)
    print lv3
    net_3=[]
    for i in lv3:
        net_3.append([(3,(i,),(i[0:-1],),i[-1]),(2,(i[0:-1],),(i[0:-1][0:-1],),i[0:-1][-1]),(1,(i[0:-1][0:-1],),(i[0:-1][0:-1][0:-1],),i[0:-1][0:-1][-1])])

    for i in net_3:
        print i
#--------------------------------------------------#
#--------------------------------------------------#
#--------------------------------------------------#
#--------------------------------------------------#
#--------------------------------------------------#
    main_nodes_list=deepcopy(nodes_list)
    count_x=0
    for one_pathway in net_3:
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
    print "link_set_s_t_list_l, ", link_set_s_t_list_l


#--------------------------------------------------#
#--------------------------------------------------#
    print node_value_dict

    print nodes_list
    print link_set_s_t_list
    print "starting write R script"

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
    print "link_set_list_rewrite_l, ", link_set_list_rewrite_l
#####################################

    nodes_R_dict=dict([])
    nodes_R_dict_l=dict([])
    for i in range(len(link_set_list_rewrite)):
        node_R_name=replace_n("R"+str(1000+i),1,"0")
        nodes_R_dict[node_R_name]=link_set_list_rewrite[i]
        if link_set_list_rewrite[i] in link_set_list_rewrite_l:
            nodes_R_dict_l[node_R_name]=link_set_list_rewrite[i]


    print "nodes_R_dict, ", nodes_R_dict
    print "nodes_R_dict_l, ", nodes_R_dict_l
    print "nodes_C_dict, ", nodes_C_dict





    C_nodes_num=len(nodes_C_dict)
    R_nodes_num=len(nodes_R_dict)

    #--------------------------------------------------#

    text_file = open("ZX_R_graph_figure_Bi.r", "w")
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




    #--------------------------------------------------#

    '''
    for one_edge in nodes_R_dict_l.keys():
        edge_info_l=nodes_R_dict_l[one_edge]
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
        '''

if __name__ == '__main__':
    #main1()
    main2()
