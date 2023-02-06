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
from numpy import *
import scipy.io
import subprocess

from zpathfinder import *


def humanbytes(B):
   'Return the given bytes as a human friendly KB, MB, GB, or TB string'
   B = float(B)
   KB = float(1024)
   MB = float(KB ** 2) # 1,048,576
   GB = float(KB ** 3) # 1,073,741,824
   TB = float(KB ** 4) # 1,099,511,627,776

   if B < KB:
      return '{0} {1}'.format(B,'Bytes' if 0 == B > 1 else 'Byte')
   elif KB <= B < MB:
      return '{0:.2f} KB'.format(B/KB)
   elif MB <= B < GB:
      return '{0:.2f} MB'.format(B/MB)
   elif GB <= B < TB:
      return '{0:.2f} GB'.format(B/GB)
   elif TB <= B:
      return '{0:.2f} TB'.format(B/TB)





# ==================================================================================== #
# parse pathway.txt
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

def read_pathway_txt(pathways_file, which_10_5): # from_which_10_5 is used to read only 10^5 pathways from the text 
    pathways_list=[]
    count_read=0
    count_parse=0
    while (True and count_parse<=100000):
        count_read+=1
        #print count_read
        one_pathway_text_list=readlines_until_empty_line(pathways_file, mark="")
        if count_read>which_10_5*100000:
            #print count_read 
            count_parse+=1
            one_pathway_list=[]
            if one_pathway_text_list==None:
                print "one_pathway_text_list==None, break!"
                break
            for one_reaction_str in one_pathway_text_list:
                one_reaction=ast.literal_eval(one_reaction_str)
                one_pathway_list.append(one_reaction)
            pathways_list.append(one_pathway_list)

    return (pathways_list, count_parse-1)


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #

def update_paired_cmpds_list(distance,paired_cmpds_set,paired_cmpds_list,all_pair_list):
    # Used in parse_one_pathway(one_pathway,paired_cmpds_list,all_pair_list)
    if paired_cmpds_set not in all_pair_list:
        all_pair_list.append(paired_cmpds_set)
        paired_cmpds_list.append([paired_cmpds_set,distance])
    else: 
        for i in range(len(paired_cmpds_list)):
            if paired_cmpds_list[i][0]==paired_cmpds_set and paired_cmpds_list[i][1]>distance: # means we have a update in distance 
                #print "update for shorter distance found"
                paired_cmpds_list[i][1]=distance
                
'''
def update_paired_cmpds_list(distance,paired_cmpds_set,paired_cmpds_list,all_pair_list):
    temp_all_paired_set=set(all_pair_list)
    current_size=len(temp_all_paired_set)
    # Used in parse_one_pathway(one_pathway,paired_cmpds_list,all_pair_list)
    if paired_cmpds_set not in temp_all_paired_set:
        all_pair_list.append(paired_cmpds_set)
        paired_cmpds_list.append([paired_cmpds_set,distance])
    else: 
        for i in range(len(paired_cmpds_list)):
            if paired_cmpds_list[i][0]==paired_cmpds_set and paired_cmpds_list[i][1]>distance: # means we have a update in distance 
                #print "update for shorter distance found"
                paired_cmpds_list[i][1]=distance
                '''
def parse_one_pathway(one_pathway,paired_cmpds_list,all_pair_list, full_pwy_bool=True):
    #--------------------------------------------------#
    # This function updats the paired_cmpds_list and all_pair_list
    # one_pathway has one target compounds and one or more reactants 
    # 1. Initialize
    pathway_length=one_pathway[0][0]
    target_cmpd=[]
    starting_cmpds=[]
    rcns_list=[]
    all_cmpds_list=[]
    cmpds_lv_list=[]
    rcns_lv_list=[]
    for i in range(pathway_length+1):
        cmpds_lv_list.append([])
        rcns_lv_list.append([])

    # 2. Change the data structure. Get rcns_list, rcns_lv_list, cmpds_lv_list, all_cmpds_list
    rcns_list=copy(one_pathway)
    for one_reaction in one_pathway:
        rcns_lv_list[one_reaction[0]].append(one_reaction)
        for one_product in one_reaction[1]: # all the intermediate compounds and target compound
            if one_product not in cmpds_lv_list[one_reaction[0]]:
                cmpds_lv_list[one_reaction[0]].append(one_product)
            if one_product not in all_cmpds_list:
                all_cmpds_list.append(one_product)
        if one_reaction[0]==1: # all the starting compounds
            for one_starting_cmpd in one_reaction[2]:
                if one_starting_cmpd not in cmpds_lv_list[0]:
                    cmpds_lv_list[0].append(one_starting_cmpd)
    starting_cmpds= cmpds_lv_list[0]
    target_cmpd=cmpds_lv_list[-1]

    #--------------------------------------------------#
    # Preparation for the compounds pairing
    cmpds_lv_number_list=[]
    for i in range(pathway_length+1):
        cmpds_lv_number_list.append(pathway_length-i)

    #--------------------------------------------------#
    # main part of the pwy to distance algorithm
    for lv_num in cmpds_lv_number_list: 
    # go down level by level, from the highest level(, that is, the target compound level)
        for cmpd_tb_paired in cmpds_lv_list[lv_num]: 
        # go through each compound in one level
            current_substitute=cmpd_tb_paired # at least one compound can be substituted
            paired_lv_num=lv_num
            cont_bool=True
            while cont_bool==True and paired_lv_num>=1:
                num_rcn_checked=0
                for one_rcn in rcns_lv_list[paired_lv_num]:
                    num_rcn_checked+=1
                    if current_substitute in one_rcn[1] and len(one_rcn[2])==1 and len(one_rcn[1])==1:

                        distance=lv_num-paired_lv_num+1
                        paired_cmpds=[(cmpd_tb_paired,),(one_rcn[2][0],)]
                        paired_cmpds_set=set(frozenset(i) for i in paired_cmpds) # set?????
                        update_paired_cmpds_list(distance,paired_cmpds_set,paired_cmpds_list,all_pair_list)

                        
                        paired_lv_num-=1
                        current_substitute=one_rcn[2][0]
                        break
                    if current_substitute in one_rcn[1] and len(one_rcn[2])>1 and len(one_rcn[1])==1:

                        distance=lv_num-paired_lv_num+1
                        paired_cmpds=[(cmpd_tb_paired,),one_rcn[2]]
                        paired_cmpds_set=set(frozenset(i) for i in paired_cmpds) # set?????
                        update_paired_cmpds_list(distance,paired_cmpds_set,paired_cmpds_list,all_pair_list)
                        
                        cont_bool=False
                        break

                    if current_substitute in one_rcn[1] and len(one_rcn[1])>1:
                        cont_bool=False
                        break
                    
                    if num_rcn_checked==len(rcns_lv_list[paired_lv_num]):
                    # In case the reaction list doesnt include the compound (which causes infinite loop)
                        cont_bool=False
                        break
    #--------------------------------------------------#
    if full_pwy_bool==True:
        num_lvs=len(cmpds_lv_list)
        half_num_lvs=num_lvs/2
        for i in range(half_num_lvs):
            for j in range(half_num_lvs):
                starting_side_lv=i
                target_side_lv=num_lvs-j-1
                if len(cmpds_lv_list[starting_side_lv])==1 and len(cmpds_lv_list[target_side_lv])==1:
                    # print cmpds_lv_list[starting_side_lv][0], starting_side_lv
                    next_lv_rctts=[] ### ensure that it is the only next level reactants
                    for one_rcn in rcns_lv_list[starting_side_lv+1]: ###
                        next_lv_rctts=next_lv_rctts+list(one_rcn[2]) ###
                    if len(set(next_lv_rctts))==1 and next_lv_rctts[0]==cmpds_lv_list[starting_side_lv][0]: ###
                        distance=target_side_lv-starting_side_lv
                        paired_cmpds=[cmpds_lv_list[starting_side_lv],cmpds_lv_list[target_side_lv]]
                        paired_cmpds_set=set(frozenset(i) for i in paired_cmpds) # set?????
                        update_paired_cmpds_list(distance,paired_cmpds_set,paired_cmpds_list,all_pair_list)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #

def rd_txt_sv_pk(read_text,save_pickle_paired_cmpds_list,save_pickle_all_pair_list):
    # read_text is the pathway file from which info is retrieved
    # save_pickle_paired_cmpds_list
    # save_pickle_all_pair_list

    #--------------------------------------------------#
    paired_cmpds_list=[]
    all_pair_list=[]
    #--------------------------------------------------#
    print "reading pathway text ..."
    pathways_file=open(read_text)
    for num_10_5 in range(30):
        print "num_10_5, ", num_10_5
        (pathways_list,pathway_numbers)=read_pathway_txt(pathways_file,num_10_5)
        if pathways_list ==[]:
            break
        count_x=0
        for one_pathway in pathways_list:
            count_x+=1
            print count_x, '/', pathway_numbers
            parse_one_pathway(one_pathway, paired_cmpds_list, all_pair_list)
    print "finished reading pathway text !!!"
    #--------------------------------------------------#
    pickle_out1=open(save_pickle_paired_cmpds_list,"wb")
    pickle.dump(paired_cmpds_list, pickle_out1)
    pickle_out1.close()
    pickle_out2=open(save_pickle_all_pair_list,"wb")
    pickle.dump(all_pair_list, pickle_out2)
    pickle_out2.close()


def op_pk_rd_txt_sv_pk(open_pickle_paired_cmpds_list,open_pickle_all_pair_list,\
    read_text,save_pickle_paired_cmpds_list,save_pickle_all_pair_list):
    #--------------------------------------------------#
    print "loading data, don't even click this window"
    pickle_in1=open(open_pickle_paired_cmpds_list,"rb")
    paired_cmpds_list=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open(open_pickle_all_pair_list,"rb")
    all_pair_list=pickle.load(pickle_in2)
    pickle_in2.close()
    print "pickle data loaded, you may click the window"
    #--------------------------------------------------#
    print "reading pathway text ..."
    pathways_file=open(read_text)
    for num_10_5 in range(30):
        print "num_10_5, ", num_10_5
        (pathways_list,pathway_numbers)=read_pathway_txt(pathways_file,num_10_5)
        if pathways_list ==[]:
            print "pathways_list==[], break!"
            break
        count_x=0
        for one_pathway in pathways_list:
            count_x+=1
            print count_x, '/', pathway_numbers
            parse_one_pathway(one_pathway, paired_cmpds_list, all_pair_list)
    print "finished reading pathway text !!!"
    #--------------------------------------------------#
    pickle_out1=open(save_pickle_paired_cmpds_list,"wb")
    pickle.dump(paired_cmpds_list, pickle_out1)
    pickle_out1.close()
    pickle_out2=open(save_pickle_all_pair_list,"wb")
    pickle.dump(all_pair_list, pickle_out2)
    pickle_out2.close()


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
    def sim_vs_dis(): # Seems not used in this _.py
    # After converting all the pathways, the sim vs dis plot can be made (not related to the main goal of this program)
    # Not finished yet 
        dis_list=[]
        sim_list=[]
        score_list=[]
        for one_pair_info in paired_cmpds_list: #[   [ { fr{}, fr{} },d ],   [ { fr{}, fr{} },d ],  [{{},{}},d], ....  ]
            if len(list(list(one_pair_info[0])))==2:
                score=maxsimscore(list(list(one_pair_info[0])[0]),list(list(one_pair_info[0])[1]),"FP4")
                distance=one_pair_info[1]
                sim=1/(distance+1.)
                dis_list.append(distance)
                sim_list.append(sim)
                score_list.append(score)

        print dis_list
        print sim_list
        print score_list


        print corrcoef(sim_list,score_list)[1,0]


        '''
        all_cmpds=[]
        for i in pair_list:
            for k in i:
                for j in k:
                    if j not in all_cmpds:
                        all_cmpds.append(j)
        print all_cmpds
        '''

        return


def open_pickles(open_pickle_paired_cmpds_list=""): # Seems not used in this _.py
    open_pickle_paired_cmpds_list= "../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle" # Comment out if needed
    pickle_in1=open(open_pickle_paired_cmpds_list,"rb")
    paired_cmpds_list=pickle.load(pickle_in1)
    pickle_in1.close()
    print len(paired_cmpds_list)
    for i in paired_cmpds_list:
        print i
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #



def main_in_main(read_pwy_txt_address): # Obtain paired compounds and distance information
    #--------------------------------------------------##--------------------------------------------------#
    read_text=read_pwy_txt_address
    #read_text="F://20181020 dt results//AKT to ARG0024//pathways.txt"
    save_pickle_paired_cmpds_list="../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"
    save_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"
    #--------------------------------------------------##--------------------------------------------------#
    open_pickle_paired_cmpds_list= "../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"# open_pickle is the pickle with saved paired compounds info
    open_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"
    #--------------------------------------------------##--------------------------------------------------#
    #rd_txt_sv_pk(read_text, save_pickle_paired_cmpds_list, save_pickle_all_pair_list)
    op_pk_rd_txt_sv_pk(open_pickle_paired_cmpds_list,open_pickle_all_pair_list,read_text,save_pickle_paired_cmpds_list,save_pickle_all_pair_list)



def main_add_fwd_bwd(network_address):
    #network_address="../zz_test_savings/networks_fwd_01"

    '''
    pickle_in1=open("../zz_test_savings/networks_bwd_001","rb")
    networks_bwd=pickle.load(pickle_in1)
    pickle_in1.close()
    '''

    pickle_in2=open(network_address,"rb")
    networks_fwd=pickle.load(pickle_in2)
    pickle_in2.close()

    save_pickle_paired_cmpds_list="../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"
    save_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"
    
    #--------------------------------------------------##--------------------------------------------------#
    open_pickle_paired_cmpds_list= "../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"# open_pickle is the pickle with saved paired compounds info
    open_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"

    #--------------------------------------------------##--------------------------------------------------#
    print "loading data, don't even click this window"
    pickle_in1=open(open_pickle_paired_cmpds_list,"rb")
    paired_cmpds_list=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open(open_pickle_all_pair_list,"rb")
    all_pair_list=pickle.load(pickle_in2)
    pickle_in2.close()
    print "pickle data loaded, you may click the window"
    #--------------------------------------------------##--------------------------------------------------#
    '''
    print "reading pathway text ..."
    count_x=0
    for one_pathway in networks_bwd:
        print count_x
        count_x+=1
        parse_one_pathway(one_pathway, paired_cmpds_list, all_pair_list)
    print "finished reading pathway text !!!"
    '''


    print "reading pathway text ..."
    count_x=0
    len_net=len(networks_fwd)
    for one_pathway in networks_fwd:
        print count_x, "/", len_net
        count_x+=1
        parse_one_pathway(one_pathway, paired_cmpds_list, all_pair_list)
    print "finished reading pathway text !!!"
    #--------------------------------------------------##--------------------------------------------------#
    pickle_out1=open(save_pickle_paired_cmpds_list,"wb")
    pickle.dump(paired_cmpds_list, pickle_out1)
    pickle_out1.close()
    pickle_out2=open(save_pickle_all_pair_list,"wb")
    pickle.dump(all_pair_list, pickle_out2)
    pickle_out2.close()

def generate_data():
    def update_paired_cmpds_list(distance,paired_cmpds_set,paired_cmpds_list,all_pair_list):
        # Used in parse_one_pathway(one_pathway,paired_cmpds_list,all_pair_list)
        if paired_cmpds_set not in all_pair_list:
            all_pair_list.append(paired_cmpds_set)
            paired_cmpds_list.append([paired_cmpds_set,distance])
        else: 
            for i in range(len(paired_cmpds_list)):
                if paired_cmpds_list[i][0]==paired_cmpds_set and paired_cmpds_list[i][1]>distance: # means we have a update in distance 
                    #print "update for shorter distance found"
                    paired_cmpds_list[i][1]=distance
        return
    #--------------------------------------------------##--------------------------------------------------#
    save_pickle_paired_cmpds_list="../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"
    save_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"
    #--------------------------------------------------##--------------------------------------------------#
    open_pickle_paired_cmpds_list= "../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"# open_pickle is the pickle with saved paired compounds info
    open_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"
    #--------------------------------------------------##--------------------------------------------------#
    print "loading data, don't even click this window"
    pickle_in1=open(open_pickle_paired_cmpds_list,"rb")
    paired_cmpds_list=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open(open_pickle_all_pair_list,"rb")
    all_pair_list=pickle.load(pickle_in2)
    pickle_in2.close()
    print "pickle data loaded, you may click the window"
    #--------------------------------------------------##--------------------------------------------------#

    from random import shuffle
    SOURCE=['NC(N)=NCCCC(N)C(=O)O',
        'NC(=O)CC(N)C(=O)O',
        'NC(CC(=O)O)C(=O)O',
        'O=C(O)C=CC(=O)O',
        'NC(=O)CCC(N)C(=O)O',
        'NC(CCC(=O)O)C(=O)O']

    for j in range(len(SOURCE)):

        all_set=set([])
        source1=set([])

        for one_pair in all_pair_list:
            #print one_pair
            all_set=all_set.union(one_pair)
            if frozenset([SOURCE[j],]) in one_pair:
                print "666"
                source1.union(one_pair)


        group_1=all_set-source1
        list_1=list(group_1)
        shuffle(list_1)
        for i in range(1000):
            paired_cmpds_set=set([frozenset([SOURCE[j],]), list_1[i]])
            print paired_cmpds_set
            update_paired_cmpds_list(6,paired_cmpds_set,paired_cmpds_list,all_pair_list)
    pickle_out1=open(save_pickle_paired_cmpds_list,"wb")
    pickle.dump(paired_cmpds_list, pickle_out1)
    pickle_out1.close()
    pickle_out2=open(save_pickle_all_pair_list,"wb")
    pickle.dump(all_pair_list, pickle_out2)
    pickle_out2.close()



def main():
    #main_in_main("../zz_metric_learn_pathways/pathways1.txt")

    
    pathways_text_num=20
    pathways_text_list=[]
    for i in range(pathways_text_num):
        print "../zz_metric_learn_pathways/pathways"+str(i)+".txt"
        pathways_text_list.append("../zz_metric_learn_pathways/pathways"+str(i)+".txt")
    for j in pathways_text_list:
        print j
        print 
        if os.path.isfile(j)==True:
            print humanbytes(os.stat(j).st_size)
            main_in_main(j)
            

    '''
    net_num=20
    net_text_list=[]
    for i in range(net_num):
        print "../zz_test_savings/networks_fwd_"+str(i)
        net_text_list.append("../zz_test_savings/networks_fwd_"+str(i))
    for j in net_text_list:
        print j
        print 
        if os.path.isfile(j)==True:
            print humanbytes(os.stat(j).st_size)
            main_add_fwd_bwd(j)
            '''




# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':

    #main_add_fwd_bwd()
    #main()
    #function_test_1()
    open_pickles()
    #main1()
    #generate_data()

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #

def function_test_1():
    # ==================================================================================== #
    #print len( get_full_ecfp("O=P(O)(O)OCC1OC(O)C(O)C(O)C1O","ECFP2",1) )
    # ==================================================================================== #
    paired_cmpds_list_test=[]
    all_pair_list_test=[]
    pathway=[(6, ('CC(C)C(N)C(=O)O',), ('CC(C)=C(N)C(=O)O',), 19),
        (5, ('CC(C)=C(N)C(=O)O',), ('CC(C)(O)C(N)C(=O)O',), 1064),
        (4, ('CC(C)(O)C(N)C(=O)O',), ('CC(C)=O', 'NCC(=O)O'), 81783),
        (3, ('CC(C)=O',), ('CC(N)CO',), 76239),
        (3, ('NC(CO)C(=O)C=O', 'NCC(=O)O'), ('NC(CO)C(=O)C(O)C(N)C(=O)O',), 62234),
        (2, ('CC(N)CO',), ('CC(=O)CO',), 85),
        (2, ('NC(CO)C(=O)C(O)C(N)C(=O)O',), ('NC(CO)C(=O)O', 'NC(CO)C(=O)O'), 409),
        (1, ('CC(=O)CO',), ('CC(=O)C(=O)O',), 10), 
        (1, ('NC(CO)C(=O)O',), ('CC(=O)C(=O)O',), 16)]
    parse_one_pathway(pathway, paired_cmpds_list_test,all_pair_list_test)
    
    for i in paired_cmpds_list_test+all_pair_list_test:
        print i
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #







