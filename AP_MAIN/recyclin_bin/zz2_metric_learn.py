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
import matplotlib.pyplot as plt
from sklearn import tree

from random import shuffle
from zpathfinder import *
from chemfuncs import *

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




# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #

def list_smiles_to_ecfp_through_dict(smiles_list, all_cmpds_ecfps_dict):
    ecfp_list=[]
    for one_smiles in smiles_list:
        ecfp_list=ecfp_list+all_cmpds_ecfps_dict[one_smiles]
    return ecfp_list


def parse_one_pair_info( one_pair_info,all_ecfps, all_cmpds_ecfps_dict):
    dimension=len(all_ecfps)
    X1i=[0]*dimension
    X2i=[0]*dimension
    X1i_ecfp_list=list_smiles_to_ecfp_through_dict(list(list(one_pair_info[0])[0]),all_cmpds_ecfps_dict)
    X2i_ecfp_list=list_smiles_to_ecfp_through_dict(list(list(one_pair_info[0])[1]),all_cmpds_ecfps_dict)
    #print X1i_MNA_list
    #print X2i_MNA_list
    distance=one_pair_info[1]
    for one_ecfp in X1i_ecfp_list:
        X1i[all_ecfps.index(one_ecfp)]=X1i_ecfp_list.count(one_ecfp)
    for one_ecfp in X2i_ecfp_list:
        X2i[all_ecfps.index(one_ecfp)]=X2i_ecfp_list.count(one_ecfp)
    Yi=distance
    return (X1i,X2i,Yi)

def list_subtract(list_a,list_b):
    list_out=[]
    for i in range(len(list_a)):
        list_out.append(list_a[i]-list_b[i])
    return list_out





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
    print len(all_ecfps)

def main():
    
    #--------------------------------------------------#
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
    #--------------------------------------------------#
    shuffle(paired_cmpds_list)
    ##----------------------------------##
    for one_pair_info in paired_cmpds_list:
        if len(one_pair_info[0])!=2:
            print one_pair_info[0]
            print "wtf?"
            paired_cmpds_list.remove(one_pair_info)
    print "shuffled!"
    ##----------------------------------##
    #--------------------------------------------------#
    all_ecfps=list(all_ecfps)


    #--------------------------------------------------#
    '''
    X1MX=[]
    X2MX=[]
    X_Diff_MX=[]
    YMX=[]
    print "paired_cmpds_list",len(paired_cmpds_list)
    count_x=0
    for one_pair_info in paired_cmpds_list:
        count_x+=1
        #print count_x
        (X1i,X2i,Yi)=parse_one_pair_info(one_pair_info,all_ecfps,all_cmpds_ecfps_dict)
        X1MX.append(X1i)
        X2MX.append(X2i)
        X_Diff_MX.append(list_subtract(X1i,X2i))
        YMX.append(Yi)
    #print X_Diff_MX
    #print YMX
    #--------------------------------------------------#
    import copy
    YYY=copy.copy(YMX)
    YMX=[]
    print YYY
    for i in YYY:
        if i >9:
            YMX.append(5)
        else:
            YMX.append(i)
            
    #--------------------------------------------------#
    
    pickle_out1=open("..\X_Diff_MX_2.pickle","wb")
    pickle.dump(X_Diff_MX, pickle_out1)
    pickle_out1.close()
    pickle_out2=open("..\YMX_2.pickle","wb")
    pickle.dump(YMX, pickle_out2)
    pickle_out2.close()
    '''

    #--------------------------------------------------#
    
    pickle_in1=open("..\zz_metric_learn_savings\X_Diff_MX_2.pickle","rb")
    X_Diff_MX=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open("..\zz_metric_learn_savings\YMX_2.pickle","rb")
    YMX=pickle.load(pickle_in2)
    pickle_in2.close()
    
    #--------------------------------------------------#

    print "all_pair_list",len(all_pair_list)

    print "start learning"

    #--------------------------------------------------#
    # Step 0
    '''
    print "starting PCA ..."
    from sklearn import decomposition
    pca=decomposition.PCA(n_components=300)
    pca.fit(X_Diff_MX)
    X_Diff_MX=pca.transform(X_Diff_MX)
    print "End PCA."
    '''
    #--------------------------------------------------#
    '''
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn import decomposition
    from sklearn import datasets
    np.random.seed(5)
    centers = [[1, 1], [-1, -1], [1, -1]]
    X = np.asarray(X_Diff_MX)
    y = np.asarray(YMX)
    #--------------------------------------------------#
    fig = plt.figure(1, figsize=(4, 3))
    plt.clf()
    ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
    plt.cla()
    #--------------------------------------------------#
    #pca = decomposition.PCA(n_components=3)
    #pca.fit(X)
    #X = pca.transform(X)
    #--------------------------------------------------#
    print "Start tSNE ..."
    from sklearn.manifold import TSNE
    tsne = TSNE(n_components=3, random_state=0)
    X = tsne.fit_transform(X)
    print "End tSNE."
    #--------------------------------------------------#
    #for name, label in [("1", 1), ("2", 2), ("3", 3)]:
    #    ax.text3D(X[y == label, 0].mean(),
    #              X[y == label, 1].mean() + 1.5,
    #              X[y == label, 2].mean(), name,
    #              horizontalalignment='center',
    #              bbox=dict(alpha=.5, edgecolor='w', facecolor='w'))


    print type(y)
    y = np.choose(y, [0,1,2,3,4,5,6,7]).astype(np.float)
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap=plt.cm.hsv,
               edgecolor='k',s=10,lw=0)

    
    ax.w_xaxis.set_ticklabels([])
    ax.w_yaxis.set_ticklabels([])
    ax.w_zaxis.set_ticklabels([])

    plt.show()
    '''
    #--------------------------------------------------#
    # Step 1
    print "len(YMX)"
    print len(YMX)
    from sklearn import tree
    X=[]
    Y=[]

    for i in range(len(X_Diff_MX)):
        if i <= 30000 and i >=32000 and YMX[i]>=4:
            X.append(X_Diff_MX[i])
            Y.append(YMX[i])
    print len(X)
    

    X = X + X_Diff_MX[1:20000]
    Y = Y + YMX[1:20000]

    clf = tree.DecisionTreeClassifier()
    #clf = tree.DecisionTreeRegressor()
    clf = clf.fit(X, Y)
    
    #print clf.decision_path(X_Diff_MX[13001])

    #--------------------------------------------------#
    '''
    from sklearn.discriminant_analysis import *
    X = X_Diff_MX[1:13000]
    Y = YMX[1:13000]
    clf=LinearDiscriminantAnalysis()
    clf.fit(X,Y)
    '''
    #--------------------------------------------------#
    #--------------------------------------------------#
    '''
    from sklearn.ensemble import RandomForestClassifier
    X = X_Diff_MX[1:13000]
    Y = YMX[1:13000]
    clf = RandomForestClassifier(n_estimators=800, max_depth=8,random_state=2)
    clf.fit(X, Y)
    '''

    #--------------------------------------------------#


    #--------------------------------------------------#
    
    count_x=0
    predicted_value=[]
    real_value_list=[]

    sim_dist_1=[]
    sim_dist_2=[]
    sim_dist_3=[]
    sim_dist_4=[]
    sim_dist_5=[]
    sim_dist_6=[]
    sim_dist_7=[]
    sim_dist_8=[]

    for i in range(30000,32000):
        prediction=clf.predict([X_Diff_MX[i]])[0]
        predicted_value.append(prediction)
        real_value=YMX[i]
        real_value_list.append(real_value)
        if prediction==real_value :
            count_x+=1

        list_a=str([prediction,real_value])
        str_a=replace_n(list_a, 3, "")

        distance=real_value
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
            print prediction
        if distance==6:
            sim_dist_6.append(prediction)
            print prediction
        if distance==7:
            sim_dist_7.append(prediction)
            print prediction
        if distance==8:
            sim_dist_8.append(prediction)
            print prediction

    print count_x/2000.
    list_a=[]
    for i in range(len(real_value_list)):
        list_a.append([real_value_list[i],predicted_value[i]]) 

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
    
    print x[3]
    print x[4]
    print x[5]
    print x[6]

    col_list=["darkred","orange","goldenrod","darkgreen","darkblue","purple","grey"]
    plt.figure()
    for i in [0,1,2,3,4,5,6]:
        plt.subplot(7, 1, i+1)
        #sns.kdeplot(x[i], bw = 0.01 , color="darkred")
        sns.distplot(x[i], hist = 1, bins = 40, rug=False, kde=True, color=col_list[i])
        plt.xlim((0,8))

    plt.show()

    come_some_music(3)

    return



if __name__ == '__main__':

    #test1()
    main()
    #open_all_all_pickles()

    #Initialize_all_cmpds_all_ecfps()
    #update_all_cmpds_all_ecfps()
