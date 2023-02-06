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
from zpathfinder import *

# ==================================================================================== #
# all_cmpds     :  list( ["X","X",...]         ) # list
# all_ecfps     :  set ( ["ecfp", "ecfp", ...] ) # set
# all_pairs     :  [{{},{}}, {{},{}}, {{},{}},... ]
# all_info      :  [   [  { fr{}, fr{} }, d  ],   [  { fr{}, fr{} }, d  ],  [  { fr{}, fr{} }, d  ], ....  ]
# ==================================================================================== #


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# ==================================================================================== #
def similarity_score(smiles_a, smiles_b, parameter_1="ECFP", parameter_2=2): # Return the similarity of two compounds
    try:
        # parameter_1 is similarity metric selected
        cmpd_a=Chem.MolFromSmiles(str(smiles_a))
        cmpd_b=Chem.MolFromSmiles(str(smiles_b))
        if (parameter_1=="top"):
            fp_a=FingerprintMols.FingerprintMol(cmpd_a)
            fp_b=FingerprintMols.FingerprintMol(cmpd_b)  
            similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
        elif (parameter_1=="MACCS"):
            fp_a=MACCSkeys.GenMACCSKeys(cmpd_a)
            fp_b=MACCSkeys.GenMACCSKeys(cmpd_b)
            similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
        elif (parameter_1=="atom_pairs"):
            fp_a=Pairs.GetAtomPairFingerprint(cmpd_a)
            fp_b=Pairs.GetAtomPairFingerprint(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="vec_pairs"):
            fp_a=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_a)
            fp_b=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="torsions"):
            fp_a=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_a)
            fp_b=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_b)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        elif (parameter_1=="FCFP"):
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2,useFeatures=True)
            fp_b=AllChem.GetMorganFingerprint(cmpd_b,parameter_2,useFeatures=True)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
        else: #ECFP
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2)
            fp_b=AllChem.GetMorganFingerprint(cmpd_b,parameter_2)
            similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    except Exception:
        print "ERROR: similarity score"
        similarity=0
    return similarity

def maxsimscore(list_a,list_b,fptype,parameter_2=2):
    score_list=[]
    for smiles_a in list_a:
        for smiles_b in list_b:
            score_list.append(similarity_score(smiles_a, smiles_b, fptype, parameter_2=2))
    return max(score_list)

def replace_n(str1, n, str2):
    letters = (
    str2 if i == n else char
        for i, char in enumerate(str1)
    )
    return ''.join(letters)


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# ==================================================================================== #
# ECFP from CDK java file
def CDK_ECFP(smiles_str,ecfp_type,iteration_number):    
    # Use java file CDKImpl class to get ECFP from cmd line
    query_str1='java -cp .;cdk-2.2.jar CDKImpl ' + smiles_str + ' ' + ecfp_type + ' ' + str(iteration_number)
    query_result = subprocess.check_output(query_str1, shell=True)
    query_result=query_result.replace('[','')
    query_result=query_result.replace(']','')
    query_result=query_result.replace(' ','')
    query_result=query_result.replace('\n','')
    query_result=query_result.replace('\r','')
    if query_result!="":
        if query_result[-1]==',':
            query_result=query_result[0:-1]
        list_of_ecfp=query_result.split(",")
    else:
        list_of_ecfp=[]
    return list_of_ecfp 

def get_full_ecfp(smiles_str,ecfp_type,iteration_number):   
    # ECFP4 + itr2 or ECFP2 + itr1
    full_ecfp_list=[]
    for i in range(iteration_number+1):
        full_ecfp_list=full_ecfp_list+CDK_ECFP(smiles_str,ecfp_type,i)
    return full_ecfp_list

def generate_all_ECFPs(list_cmpds,ecfp_type="ECFP2",iteration_number=1):
# return a list of ECFPs of all depth for a list of compounds (UNIQUE!!!)
    all_ecfps=set([])
    for smiles_a in list_cmpds:
        discriptors = get_full_ecfp(smiles_a,ecfp_type,iteration_number)
        print smiles_a
        all_ecfps=all_ecfps.union(set(discriptors))
    return all_ecfps

def update_all_ECFPs(new_list_cmpds,list_cmpds,all_ecfps,ecfp_type="ECFP2",iteration_number=1):
    for smiles_a in new_list_cmpds:
        if smiles_a not in list_cmpds:
            discriptors = get_full_ecfp(smiles_a,ecfp_type,iteration_number)
            print smiles_a
            all_ecfps=all_ecfps.union(set(discriptors))
    return all_ecfps


def generate_all_cmpds_ecfps_dict(list_cmpds,ecfp_type="ECFP2",iteration_number=1):
    all_cmpds_ecfps_dict=dict([])
    for smiles_a in list_cmpds:
        print smiles_a
        all_cmpds_ecfps_dict[smiles_a]=get_full_ecfp(smiles_a,ecfp_type,iteration_number)
    return all_cmpds_ecfps_dict

def update_all_cmpds_ecfps_dict(new_list_cmpds,list_cmpds,all_cmpds_ecfps_dict,ecfp_type="ECFP2",iteration_number=1):
    for smiles_a in new_list_cmpds:
        if smiles_a not in list_cmpds:
            print smiles_a
            all_cmpds_ecfps_dict[smiles_a]=get_full_ecfp(smiles_a,ecfp_type,iteration_number)
    return all_cmpds_ecfps_dict

def generate_list_and_dict(list_cmpds,ecfp_type="ECFP2",iteration_number=1):
    all_ecfps=set([])
    all_cmpds_ecfps_dict=dict([])
    for smiles_a in list_cmpds:
        discriptors = get_full_ecfp(smiles_a,ecfp_type,iteration_number)
        print smiles_a
        all_cmpds_ecfps_dict[smiles_a]=discriptors
        all_ecfps=all_ecfps.union(set(discriptors))
    return all_ecfps,all_cmpds_ecfps_dict



def update_all_ecfps_and_all_cmpds_ecfps_dict(new_list_cmpds,list_cmpds,all_ecfps,all_cmpds_ecfps_dict,ecfp_type="ECFP2",iteration_number=1):
    count_x=0
    for smiles_a in new_list_cmpds:
        if smiles_a not in list_cmpds:
            discriptors = get_full_ecfp(smiles_a,ecfp_type,iteration_number)
            print count_x
            count_x+=1
            all_cmpds_ecfps_dict[smiles_a]=discriptors
            all_ecfps=all_ecfps.union(set(discriptors))
    return all_ecfps,all_cmpds_ecfps_dict





# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
def open_pickles(open_pickle_paired_cmpds_list,open_pickle_all_pair_list):
    
    print "loading data, don't even click this window"
    pickle_in1=open(open_pickle_paired_cmpds_list,"rb")
    paired_cmpds_list=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open(open_pickle_all_pair_list,"rb")
    all_pair_list=pickle.load(pickle_in2)
    pickle_in2.close()
    print "pickle data loaded, you may click the window"
    return [paired_cmpds_list,all_pair_list]

def get_all_cmpds(paired_cmpds_list, all_pair_list):
    all_cmpds=[]
    for i in all_pair_list:
        for k in i:
            for j in k:
                if j not in all_cmpds:
                    all_cmpds.append(j)
    return all_cmpds


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
def Initialize_all_cmpds_all_ecfps():
    #--------------------------------------------------#
    open_pickle_paired_cmpds_list= "../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"
    open_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"
    open_lists=open_pickles(open_pickle_paired_cmpds_list,open_pickle_all_pair_list)
    paired_cmpds_list=open_lists[0]
    all_pair_list=open_lists[1]
    
    all_cmpds = get_all_cmpds(paired_cmpds_list, all_pair_list)
    print len(all_cmpds)
    #all_ecfps = generate_all_ECFPs(all_cmpds,ecfp_type="ECFP2",iteration_number=1)
    #all_cmpds_ecfps_dict=generate_all_cmpds_ecfps_dict(all_cmpds,ecfp_type="ECFP2",iteration_number=1)
    (all_ecfps,all_cmpds_ecfps_dict)=generate_list_and_dict(all_cmpds,ecfp_type="ECFP2",iteration_number=1)
    #--------------------------------------------------#
    pickle_out1=open("../zz_metric_learn_savings/zz1_all_cmpds_000","wb")
    pickle.dump(all_cmpds, pickle_out1)
    pickle_out1.close()
    pickle_out2=open("../zz_metric_learn_savings/zz1_all_ecfps_000","wb")
    pickle.dump(all_ecfps, pickle_out2)
    pickle_out2.close()
    pickle_out3=open("../zz_metric_learn_savings/zz1_all_cmpds_ecfps_dict_000","wb")
    pickle.dump(all_cmpds_ecfps_dict, pickle_out3)
    pickle_out3.close()
    print "DONE!!!"

def update_all_cmpds_all_ecfps():
    print "start update_all_cmpds_all_ecfps ..."
    open_pickle_paired_cmpds_list= "../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"
    open_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"
    open_lists=open_pickles(open_pickle_paired_cmpds_list,open_pickle_all_pair_list)
    paired_cmpds_list=open_lists[0]
    all_pair_list=open_lists[1]
    new_all_cmpds = get_all_cmpds(paired_cmpds_list, all_pair_list)
    print "total compounds number: ", len(new_all_cmpds)
    #--------------------------------------------------#
    pickle_in1=open("../zz_metric_learn_savings/zz1_all_cmpds_000","rb")
    all_cmpds=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open("../zz_metric_learn_savings/zz1_all_ecfps_000","rb")
    all_ecfps=pickle.load(pickle_in2)
    pickle_in2.close()
    pickle_in3=open("../zz_metric_learn_savings/zz1_all_cmpds_ecfps_dict_000","rb")
    all_cmpds_ecfps_dict=pickle.load(pickle_in3)
    pickle_in3.close()
    #--------------------------------------------------#
    #all_ecfps=update_all_ECFPs(new_all_cmpds,all_cmpds,all_ecfps,ecfp_type="ECFP2",iteration_number=1)
    #all_cmpds_ecfps_dict=update_all_cmpds_ecfps_dict(new_all_cmpds,all_cmpds,all_cmpds_ecfps_dict,ecfp_type="ECFP2",iteration_number=1)
    (all_ecfps,all_cmpds_ecfps_dict)=update_all_ecfps_and_all_cmpds_ecfps_dict(new_all_cmpds,all_cmpds,all_ecfps,all_cmpds_ecfps_dict,ecfp_type="ECFP2",iteration_number=1)
    #--------------------------------------------------#
    pickle_out1=open("../zz_metric_learn_savings/zz1_all_cmpds_000","wb")
    pickle.dump(new_all_cmpds, pickle_out1)
    pickle_out1.close()
    pickle_out2=open("../zz_metric_learn_savings/zz1_all_ecfps_000","wb")
    pickle.dump(all_ecfps, pickle_out2)
    pickle_out2.close()
    pickle_out3=open("../zz_metric_learn_savings/zz1_all_cmpds_ecfps_dict_000","wb")
    pickle.dump(all_cmpds_ecfps_dict, pickle_out3)
    pickle_out3.close()

    print("DONE UPDATE update_all_cmpds_all_ecfps")






# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
def open_all_all_pickles():
    open_pickle_paired_cmpds_list= "../zz_metric_learn_savings/zz0_paired_cmpds_list_000.pickle"
    open_pickle_all_pair_list="../zz_metric_learn_savings/zz0_all_pair_list_000.pickle"
    open_lists=open_pickles(open_pickle_paired_cmpds_list,open_pickle_all_pair_list)
    paired_cmpds_list=open_lists[0]
    all_pair_list=open_lists[1]
    #--------------------------------------------------#
    pickle_in1=open("../zz_metric_learn_savings/zz1_all_cmpds_000","rb")
    all_cmpds=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open("../zz_metric_learn_savings/zz1_all_ecfps_000","rb")
    all_ecfps=pickle.load(pickle_in2)
    pickle_in2.close()
    pickle_in3=open("../zz_metric_learn_savings/zz1_all_cmpds_ecfps_dict_000","rb")
    all_cmpds_ecfps_dict=pickle.load(pickle_in3)
    pickle_in3.close()
    print len(paired_cmpds_list)
    #shuffle(paired_cmpds_list)

    dis_list=[]
    sim_list=[]
    score_list=[]
    count_x=0
    for one_pair_info in paired_cmpds_list: #[   [ { fr{}, fr{} },d ],   [ { fr{}, fr{} },d ],  [{{},{}},d], ....  ]
        count_x+=1
        #print count_x
        if count_x>=10002:
            break
        if len(list(list(one_pair_info[0])))==2:
            score=maxsimscore(list(list(one_pair_info[0])[0]),list(list(one_pair_info[0])[1]),"ECFP")
            if score==0:
                continue
            distance=one_pair_info[1]
            sim=1/(distance+1.)
            dis_list.append(distance)
            sim_list.append(sim)
            score_list.append(score)

    print dis_list
    print sim_list
    print score_list

    sim_dist_1=[]
    sim_dist_2=[]
    sim_dist_3=[]
    sim_dist_4=[]
    sim_dist_5=[]
    sim_dist_6=[]
    sim_dist_7=[]
    sim_dist_8=[]
    for i in range(6031):
        #print i
        distance=dis_list[i]
        prediction=score_list[i]
        if distance==1:
            sim_dist_1.append(prediction)
        if distance==2:
            sim_dist_2.append(prediction)
        if distance==3:
            sim_dist_3.append(prediction)
        if distance==4:
            sim_dist_4.append(prediction)
        if distance==5:
            sim_dist_5.append(prediction)
        if distance==6:
            sim_dist_6.append(prediction)
            print "6 : ", prediction
        if distance==7:
            sim_dist_7.append(prediction)
            print "7 : ", prediction
        if distance==8:
            sim_dist_8.append(prediction)
            print "8 : ", prediction

    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats
    x=[]
    x.append(np.array(sim_dist_1))
    x.append(np.array(sim_dist_2))
    x.append(np.array(sim_dist_3))
    x.append(np.array(sim_dist_4))
    x.append(np.array(sim_dist_5))
    x.append(np.array(sim_dist_6))
    x.append(np.array(sim_dist_7))

    col_list=["darkred","orange","goldenrod","darkgreen","darkblue","purple","grey"]
    plt.figure()
    for i in [0,1,2,3,4,5,6]:
        plt.subplot(7, 1, i+1)
        #sns.kdeplot(x[i], bw = 0.01 , color="darkred")
        sns.distplot(x[i], hist = 1, bins = 40, rug=False, kde=True, color=col_list[i])
        plt.xlim((0,1))





    '''
    sns.kdeplot(x1,bw= 0.1 , color="darkred")
    sns.distplot(x1, hist=1, bins=40, rug=False, kde=False, color="red")
    sns.distplot(x1, hist=1, rug=False, kde=False, hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1, "color": "red"})
    '''


    '''
    font = {'family' : 'normal', 'weight' : 'bold', 'size': 20}
    plt.rc('font', **font)
    plt.rc('font', size=20)          # controls default text sizes
    plt.rc('axes', titlesize=20)     # fontsize of the axes title
    plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
    plt.rc('legend', fontsize=20)    # legend fontsize
    plt.rc('figure', titlesize=30)  # fontsize of the figure title
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel('Similarity Scores',**font)
    plt.ylabel('PDF',**font)
    plt.title('Distribution Plot',**font)
    '''
        


    plt.show()


    from numpy import * 
    print corrcoef(sim_list,score_list)[1,0]







def main():
    #--------------------------------------------------#
    from zz0_pwy_to_dis import *

    #--------------------------------------------------#
    # Update the zz0 pickles from pathways.txt(functions are all in zz0_pwy_to_dis.py )
    '''
    pathways_text_num=20
    pathways_text_list=[]
    for i in range(pathways_text_num):
        pathways_text_list.append("d:/zz_metric_learn_pathways/pathways"+str(i)+".txt")
    for j in pathways_text_list:
        if os.path.isfile(j)==true:
            print j 
            print humanbytes(os.stat(j).st_size)
            main_in_main(j)
            '''

    update_all_cmpds_all_ecfps()


def test1():



    ecfp_type="ECFP2"
    smiles_str="CCCCC"
    iteration_number=0
    query_str1='java -cp .;cdk-2.2.jar CDKImpl ' + smiles_str + ' ' + ecfp_type + ' ' + str(iteration_number)
    print(query_str1)
    print(os.path.dirname(__file__))
    query_result = os.system(query_str1)
    print(query_result)
    query_result = subprocess.Popen(query_str1, shell=True ,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)

    for lines in query_result.stdout:
        print lines.decode("gb2312")

if __name__ == '__main__':

    #main()
    #open_all_all_pickles()
    test1()
    #Initialize_all_cmpds_all_ecfps()
    #update_all_cmpds_all_ecfps()


    '''
    smiles_a="C=C(OP(=O)(O)O)C(=O)O"
    smiles_b="O=C(O)C(=O)CC(O)C(O)C(O)COP(=O)(O)O"
    fptype="ECFP"
    print similarity_score(smiles_a, smiles_b, fptype, parameter_2=2)
    '''
    '''
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
    '''