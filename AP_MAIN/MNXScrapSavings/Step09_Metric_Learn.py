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
import pandas as pd
import scipy.io
import subprocess
from random import shuffle
# ==================================================================================== #
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import pylab as pl
import matplotlib.mlab as mlab
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
def parse_one_pair_info(one_pair_tuple, smiles_VAEVEC_dict,smiles_pairs_distance_dict):
    if one_pair_tuple[0] not in smiles_VAEVEC_dict.keys() or one_pair_tuple[1] not in smiles_VAEVEC_dict.keys():
        return None
    X1i=smiles_VAEVEC_dict[one_pair_tuple[0]]
    X2i=smiles_VAEVEC_dict[one_pair_tuple[1]]
    Yi=smiles_pairs_distance_dict[one_pair_tuple]
    return (X1i,X2i,Yi) 
# ==================================================================================== #
def replace_n(str1, n, str2):
    letters = (
    str2 if i == n else char
        for i, char in enumerate(str1)
    )
    return ''.join(letters)
# ==================================================================================== #
def list_subtract(list_a,list_b):
    list_out=[]
    for i in range(len(list_a)):
        list_out.append(list_a[i]-list_b[i])
    return list_out
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
def main_1():
    #-------------------- (1) --------------------#
    # Open the Step04 temp file and retrieve useful info.
    pickle_in1=open(".\Step08_X_Diff_MX","rb")
    X_Diff_MX=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open(".\Step08_YMX","rb")
    YMX=pickle.load(pickle_in2)
    pickle_in2.close()
    print "Done importing data points"
    #-------------------- (2) --------------------#
    from itertools import cycle
    from sklearn import datasets
    from sklearn.metrics import roc_curve, auc
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import label_binarize
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.ensemble import RandomForestClassifier
    from scipy import interp
    X = np.array(X_Diff_MX)
    y = np.array(YMX)
    n_samples, n_features = X.shape
    # Run classifier with cross-validation and plot ROC curves
    # Binarize the output
    y = label_binarize(y, classes=[0,1,2,3,4,5])
    n_classes = y.shape[1]
    # shuffle and split training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,random_state=0)
    n_estimators = [1, 2, 4, 8, 16, 32, 64, 100, 200]
    train_results=[]
    test_results=[]
    #-------------------- (3) --------------------#
    for estimator in n_estimators:
        # Learn to predict each class against the other
        classifier = OneVsRestClassifier(RandomForestClassifier(n_estimators=estimator, n_jobs=-1))
        y_score = classifier.fit(X_train, y_train).predict(X_test)
        y_train_pred = classifier.fit(X_train, y_train).predict(X_train)
        # Compute ROC curve and ROC area for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])
        print "Computing n_estimator = ", estimator
        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
        roc_auc_1=(roc_auc[1]+roc_auc[2]+roc_auc[3]+roc_auc[4])/4.0
        test_results.append(roc_auc_1)

        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_train[:, i], y_train_pred[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])
        fpr["micro"], tpr["micro"], _ = roc_curve(y_train.ravel(), y_train_pred.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
        roc_auc_2=(roc_auc[1]+roc_auc[2]+roc_auc[3]+roc_auc[4])/4.0
        train_results.append(roc_auc_2)
    import matplotlib.pyplot as plt
    plt.plot(n_estimators, train_results, 'b', label='Test AUC')
    plt.plot(n_estimators, test_results, 'r', label='Train AUC')
    plt.ylabel('AUC score')
    plt.xlabel('n_estimators')
    plt.show()
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
def main_2():
    #-------------------- (0) --------------------#
    pickle_in1=open(".\Step08_X_Diff_MX","rb")
    X_Diff_MX=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in2=open(".\Step08_YMX","rb")
    YMX=pickle.load(pickle_in2)
    pickle_in2.close()
    print "Done importing data points"
    #-------------------- (1) --------------------#
    from sklearn.ensemble import RandomForestClassifier
    X = X_Diff_MX[1:9000]
    Y = YMX[1:9000]
    clf = RandomForestClassifier(n_estimators=50, max_depth=8, random_state=2)
    clf.fit(X, Y)
    #-------------------- (1) --------------------#
    count_x=0
    predicted_value=[]
    real_value_list=[]

    sim_dist_1=[]
    sim_dist_2=[]
    sim_dist_3=[]
    sim_dist_4=[]

    for i in range(9001,11000):
        prediction=clf.predict([X_Diff_MX[i]])[0]
        predicted_value.append(prediction)
        real_value=YMX[i]
        real_value_list.append(real_value)
        if prediction==real_value and ((prediction-real_value)<=1 and (prediction-real_value)>=-1):
            count_x+=1

        distance=real_value
        if distance==1:
            sim_dist_1.append(prediction)
        if distance==2:
            sim_dist_2.append(prediction)
        if distance==3:
            sim_dist_3.append(prediction)
        if distance==4:
            sim_dist_4.append(prediction)

    print count_x/2000.
    list_a=[]
    for i in range(len(real_value_list)):
        list_a.append([real_value_list[i],predicted_value[i]]) 

    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats

    x1 = np.array(sim_dist_1)
    sns.distplot(x1, hist=False, rug=False, color="red")
    x2 = np.array(sim_dist_2)
    sns.distplot(x2, hist=False, rug=False, color="gold")     
    x3 = np.array(sim_dist_3)
    sns.distplot(x3, hist=False, rug=False, color="green")     
    x4 = np.array(sim_dist_4)
    sns.distplot(x4, hist=False, rug=False, color="blue") 
    font = {'family' : 'normal', 'weight' : 'bold', 'size': 10}
    plt.rc('font', **font)
    plt.rc('font', size=10)          # controls default text sizes
    plt.rc('axes', titlesize=10)     # fontsize of the axes title
    plt.rc('axes', labelsize=10)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
    plt.rc('legend', fontsize=10)    # legend fontsize
    plt.rc('figure', titlesize=20)  # fontsize of the figure title
    plt.xlabel('Estimated Similarity Scores',**font)
    plt.ylabel('PDF',**font)
    plt.title('Distribution Plot',**font)
    plt.show()
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
def main_3():
    #-------------------- (0) --------------------#
    # Classification and ROC analysis
    import numpy as np
    import matplotlib.pyplot as plt
    from itertools import cycle
    from sklearn import svm, datasets
    from sklearn.metrics import roc_curve, auc
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import label_binarize
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.ensemble import RandomForestClassifier
    from scipy import interp
    X = np.array(X_Diff_MX)
    y = np.array(YMX)
    n_samples, n_features = X.shape
    # Run classifier with cross-validation and plot ROC curves
    # Binarize the output
    y = label_binarize(y, classes=[0,1,2,3,4,5])
    n_classes = y.shape[1]

    # shuffle and split training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,
                                                        random_state=0)
    n_estimators = [1, 2, 4, 8, 16, 32, 64, 100, 200]
    train_results=[]
    test_results=[]
    
    for estimator in n_estimators: # Use sample code for plotting OneVsRestClassifier AUC from sklearn website
        # Learn to predict each class against the other
        classifier = OneVsRestClassifier(RandomForestClassifier(n_estimators=estimator, n_jobs=-1))
        #classifier = OneVsRestClassifier(svm.SVC(kernel='linear', probability=True,random_state=0))
        y_score = classifier.fit(X_train, y_train).predict(X_test)
        y_train_pred = classifier.fit(X_train, y_train).predict(X_train)
        # Compute ROC curve and ROC area for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])
        print "1"    # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
        roc_auc_1=(roc_auc[1]+roc_auc[2]+roc_auc[3]+roc_auc[4])/4.0
        test_results.append(roc_auc_1)

        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_train[:, i], y_train_pred[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])
        print "1"    # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_train.ravel(), y_train_pred.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
        roc_auc_2=(roc_auc[1]+roc_auc[2]+roc_auc[3]+roc_auc[4])/4.0
        train_results.append(roc_auc_2)
    
    plt.figure()
    lw = 2
    plt.plot(fpr[1], tpr[1], color='red',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[1])
    plt.plot(fpr[2], tpr[2], color='gold',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
    plt.plot(fpr[3], tpr[3], color='green',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[3])
    plt.plot(fpr[4], tpr[4], color='blue',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[4])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':

    # metric learn based on random forest (tuned parameter)
    main_1()

    # AUC score plot over num extimators
    # Other parameters are also plotted vs. AUC scores to tune the model
    main_2()

    # AUC-ROC
    main_3()
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
def distance_metric_learn(): # Bad Results
   
    print "starting PCA ..."
    from sklearn import decomposition
    pca=decomposition.PCA(n_components=100)
    pca.fit(X_Diff_MX)
    X_Diff_MX=pca.transform(X_Diff_MX)
    print "End PCA."
    

    
    print "start learning"
    from metric_learn import LMNN
    from metric_learn import LSML_Supervised
    from metric_learn import NCA
    from metric_learn import SDML_Supervised
    from metric_learn import ITML_Supervised
    from metric_learn import RCA_Supervised
    from metric_learn import LFDA
    X = np.asarray(X_Diff_MX[1:9000],dtype=np.float32)
    Y = np.asarray(YMX[1:9000],dtype=np.int32)
    print Y
    #clf = LMNN()                  # LMNN: array is too big
    clf = LSML_Supervised(num_constraints=200)        # NCA : singular matrix, the algorithms need to inverse the matrix
    #clf = NCA(max_iter=1000, learning_rate=0.01)      # LSML: iterator too large
    #clf = SDML_Supervised(num_constraints=200)        # SDML: system ill-conditioned
    #clf = ITML_Supervised(num_constraints=100)        # ITML: performace not good at all
    #clf = RCA_Supervised(num_chunks=20, chunk_size=20)# RCA : training failed all values NaN
    #clf = LFDA(k=10)      # LFDA: error in calculation

    clf.fit(X, Y)
    metric_learned=clf.metric()
    print metric_learned

    for i in range(9001,11000):
        X_Diff=np.asarray([X_Diff_MX[i]],dtype=np.float32)
        print sqrt(X_Diff.dot(metric_learned).dot(X_Diff.transpose())),
        print YMX[i]
    count_x=0
    predicted_value=[]
    real_value_list=[]
    sim_dist_1=[]
    sim_dist_2=[]
    sim_dist_3=[]
    sim_dist_4=[]
    for i in range(9001,11000):
        X_Diff=np.asarray([X_Diff_MX[i]],dtype=np.float32)
        prediction=sqrt(X_Diff.dot(metric_learned).dot(X_Diff.transpose()))
        #prediction=clf.predict([X_Diff_MX[i]])[0]
        predicted_value.append(prediction)
        real_value=YMX[i]
        real_value_list.append(real_value)
        if prediction==real_value or ((prediction-real_value)<=1 and (prediction-real_value)>=-1):
            count_x+=1

        list_a=str([prediction,real_value])
        str_a=replace_n(list_a, 3, "")
        #print str_a, ';',

        distance=real_value
        if distance==1:
            sim_dist_1.append(prediction)
        if distance==2:
            sim_dist_2.append(prediction)
        if distance==3:
            sim_dist_3.append(prediction)
        if distance==4:
            sim_dist_4.append(prediction)
    print count_x/2000.
    list_a=[]
    for i in range(len(real_value_list)):
        list_a.append([real_value_list[i],predicted_value[i]]) 
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats
    x1 = np.array(sim_dist_1)
    sns.distplot(x1, hist=False, rug=False, color="red")
    x2 = np.array(sim_dist_2)
    sns.distplot(x2, hist=False, rug=False, color="gold")     
    x3 = np.array(sim_dist_3)
    sns.distplot(x3, hist=False, rug=False, color="green")     
    x4 = np.array(sim_dist_4)
    sns.distplot(x4, hist=False, rug=False, color="blue") 
    font = {'family' : 'normal', 'weight' : 'bold', 'size': 10}
    plt.rc('font', **font)
    plt.rc('font', size=10)          # controls default text sizes
    plt.rc('axes', titlesize=10)     # fontsize of the axes title
    plt.rc('axes', labelsize=10)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
    plt.rc('legend', fontsize=10)    # legend fontsize
    plt.rc('figure', titlesize=20)  # fontsize of the figure title
    plt.xlabel('Estimated Similarity Scores',**font)
    plt.ylabel('PDF',**font)
    plt.title('Distribution Plot',**font)
    plt.show()
    