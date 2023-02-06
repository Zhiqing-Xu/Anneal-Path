# ==================================================================================== #
#!/usr/bin/python
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit import DataStructs

from rdkit.Chem.AtomPairs import Torsions
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
import scipy.io
import subprocess
from random import shuffle
# ==================================================================================== #
# ==================================================================================== #
def similarity_metric_select(fp_a,fp_b,parameter_1,parameter=2):
    if (parameter_1=="top"):
        similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
    elif (parameter_1=="MACCS"):
        similarity=DataStructs.FingerprintSimilarity(fp_a,fp_b)
    elif (parameter_1=="atom_pairs"):
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    elif (parameter_1=="vec_pairs"):
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    elif (parameter_1=="torsions"):
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    elif (parameter_1=="FCFP"):
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    else: # ECFP
        similarity=DataStructs.DiceSimilarity(fp_a,fp_b)
    return similarity


def generate_fingerprint(smiles_a, parameter_1, parameter_2=2):
    try: 
        cmpd_a=Chem.MolFromSmiles(str(smiles_a))
        if (parameter_1=="top"):
            fp_a=FingerprintMols.FingerprintMol(cmpd_a)
        elif (parameter_1=="MACCS"):
            fp_a=MACCSkeys.GenMACCSKeys(cmpd_a)
        elif (parameter_1=="atom_pairs"):
            fp_a=Pairs.GetAtomPairFingerprint(cmpd_a)
        elif (parameter_1=="vec_pairs"):
            fp_a=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_a)
        elif (parameter_1=="torsions"):
            fp_a=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_a)
        elif (parameter_1=="FCFP"):
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2,useFeatures=True)
        else: #ECFP
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2)
    except Exception:
        print "ERROR: generate fingerprint"
        cmpd_a=Chem.MolFromSmiles(str('O'))
        if (parameter_1=="top"):
            fp_a=FingerprintMols.FingerprintMol(cmpd_a)
        elif (parameter_1=="MACCS"):
            fp_a=MACCSkeys.GenMACCSKeys(cmpd_a)
        elif (parameter_1=="atom_pairs"):
            fp_a=Pairs.GetAtomPairFingerprint(cmpd_a)
        elif (parameter_1=="vec_pairs"):
            fp_a=Pairs.GetAtomPairFingerprintAsBitVect(cmpd_a)
        elif (parameter_1=="torsions"):
            fp_a=Torsions.GetTopologicalTorsionFingerprintAsIntVect(cmpd_a)
        elif (parameter_1=="FCFP"):
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2,useFeatures=True)
        else: #ECFP
            fp_a=AllChem.GetMorganFingerprint(cmpd_a,parameter_2)
    return fp_a

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



# ==================================================================================== #
def remove_cofactors(one_rxn_info):
    # More cofactor id might be updated later
    cofactor_mnx_id=["MNXM01","MNXM1","MNXM2","MNXM3","MNXM4","MNXM5","MNXM6","MNXM7","MNXM8",\
                    "MNXM9","MNXM10","MNXM11","MNXM12","MNXM13","MNXM14","MNXM15","MNXM15"]
    for one_cofactor in cofactor_mnx_id:
        if one_cofactor in one_rxn_info[1]:
            one_rxn_info[0].remove(one_rxn_info[0][one_rxn_info[1].index(one_cofactor)])
            one_rxn_info[1].remove(one_cofactor)
        if one_cofactor in one_rxn_info[3]:
            one_rxn_info[2].remove(one_rxn_info[2][one_rxn_info[3].index(one_cofactor)])
            one_rxn_info[3].remove(one_cofactor)     
    return one_rxn_info

def remove_trivial_cmpds(one_rxn_info,cmpd_mnxid_smiles_dict):
    # Remove metal ions, small compounds, NA SMILES strings and empty strings from the list
    # The algorithm here is not perfected, probably need to take into account more conditions
    #print one_rxn_info
    rxn_info=[[],[],[],[]]
    for i in range(len(one_rxn_info[1])):
        one_cmpd=one_rxn_info[1][i]
        cmpd_smiles=cmpd_mnxid_smiles_dict[one_cmpd]
        if len(cmpd_smiles)>6 and cmpd_smiles.find(".")==-1:
            rxn_info[0].append(one_rxn_info[0][i])
            rxn_info[1].append(one_rxn_info[1][i])        

    for i in range(len(one_rxn_info[3])):
        one_cmpd=one_rxn_info[3][i]
        cmpd_smiles=cmpd_mnxid_smiles_dict[one_cmpd]
        if len(cmpd_smiles)>6 and cmpd_smiles.find(".")==-1:
            rxn_info[2].append(one_rxn_info[2][i])
            rxn_info[3].append(one_rxn_info[3][i])        
    #print rxn_info
    return rxn_info
# ==================================================================================== #
def main():
    #-------------------- (1) --------------------#
    # Open all pickles and retrieve compound/reaction info.
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict1","rb")
    cmpd_mnxid_smiles_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict2","rb")
    cmpd_mnxid_smiles_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict3","rb")
    cmpd_mnxid_smiles_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict4","rb")
    cmpd_mnxid_smiles_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict5","rb")
    cmpd_mnxid_smiles_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict6","rb")
    cmpd_mnxid_smiles_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./Step01_cmpd_mnxid_smiles_dict7","rb")
    cmpd_mnxid_smiles_dict7=pickle.load(pickle_in1)
    pickle_in1.close()

    cmpd_mnxid_smiles_dict=cmpd_mnxid_smiles_dict1.copy()
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict2)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict3)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict4)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict5)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict6)
    cmpd_mnxid_smiles_dict.update(cmpd_mnxid_smiles_dict7)

    # Add meaningless keys to the dict to avoid errors
    cmpd_mnxid_smiles_dict["BIOMASS"]=""
    cmpd_mnxid_smiles_dict["MNXM0"]=""

    # Could save another pickle for the merged dictionary.
    # However, the runtime here is negligible.
    print "done merge" 
    pickle_in1=open("./Step02_mnxid_rxn_list","rb")
    mnxid_rxn_list=pickle.load(pickle_in1)
    pickle_in1.close()
    #-------------------- (2) --------------------#
    # 1st Screening
    screened_mnxid_rxn_list_1=[]
    for one_rxn_info in mnxid_rxn_list:
        one_rxn_info=remove_cofactors(one_rxn_info)
        one_rxn_info=remove_trivial_cmpds(one_rxn_info,cmpd_mnxid_smiles_dict)

        if len(one_rxn_info[0])!=0 and len(one_rxn_info[2])!=0 \
            and len(one_rxn_info[0]) == len(one_rxn_info[1]) \
            and len(one_rxn_info[2]) == len(one_rxn_info[3]):
            screened_mnxid_rxn_list_1.append(one_rxn_info)
    print "Done 1st Screening"
    del(mnxid_rxn_list)
    mnxid_rxn_list=copy.copy(screened_mnxid_rxn_list_1)
    #-------------------- (3) --------------------#
    # 2nd Screening, get reactions in the form (A -> B) 
    # since we only reconstructs linear pathways here.
    screened_mnxid_rxn_list_2=[]
    count_x=0
    count_y=0
    count_u=0
    count_t=0

    pair_AB=[]
    pair_multi_rct=[]
    pair_multi_prt=[]
    pair_multi_r_p=[]


    parameter_1="ECFP"
    for one_rxn_info in mnxid_rxn_list:
        #print one_rxn_info
        
        if len(one_rxn_info[1])==1 and len(one_rxn_info[3])==1 \
            and one_rxn_info[1][0]!=one_rxn_info[3][0] \
            and set([one_rxn_info[1][0],one_rxn_info[3][0]]) not in screened_mnxid_rxn_list_2: # Consider simple reactions here (A->B)
            screened_mnxid_rxn_list_2.append(set([one_rxn_info[1][0],one_rxn_info[3][0]]))
            count_x+=1
            ssscore=similarity_score( cmpd_mnxid_smiles_dict[one_rxn_info[1][0]] , cmpd_mnxid_smiles_dict[one_rxn_info[3][0]] , parameter_1, parameter_2=2 )
            if ssscore!=1.0:
                pair_AB.append(         ssscore )
        elif len(one_rxn_info[1])==1 and len(one_rxn_info[3])>1:
            count_y+=1
            ssscore=similarity_score( cmpd_mnxid_smiles_dict[one_rxn_info[1][0]] , cmpd_mnxid_smiles_dict[one_rxn_info[3][0]] , parameter_1, parameter_2=2 )
            if ssscore!=1.0:
                pair_multi_rct.append(  ssscore )
        elif len(one_rxn_info[1])>1 and len(one_rxn_info[3])==1:
            count_u+=1
            ssscore=similarity_score( cmpd_mnxid_smiles_dict[one_rxn_info[1][0]] , cmpd_mnxid_smiles_dict[one_rxn_info[3][0]] , parameter_1, parameter_2=2 )
            if ssscore!=1.0:
                pair_multi_prt.append(  ssscore )

        elif len(one_rxn_info[1])>1 and len(one_rxn_info[3])>1:
            count_t+=1
            ssscore=similarity_score( cmpd_mnxid_smiles_dict[one_rxn_info[1][0]] , cmpd_mnxid_smiles_dict[one_rxn_info[3][0]] , parameter_1, parameter_2=2 )
            if ssscore!=1.0:
                pair_multi_r_p.append(  ssscore )

        else:
            aaa=1

    print count_x,count_y,count_u,count_t




    import statistics

    print pair_AB
    print pair_multi_rct
    print pair_multi_prt
    print pair_multi_r_p


    print statistics.mean(pair_AB)
    print statistics.mean(pair_multi_rct)
    print statistics.mean(pair_multi_prt)
    print statistics.mean(pair_multi_r_p)



    print "Done 2nd Screening"
    print len(screened_mnxid_rxn_list_2)
    #-------------------- (4) --------------------#
    # Get all SMILES needed for further analysis
    # No need to convert all SMILES (through VAE) to latent vectors
    Used_MNXid_set=set([])
    for i in screened_mnxid_rxn_list_2:
        Used_MNXid_set=Used_MNXid_set.union(i)
    print Used_MNXid_set

    Used_SMILES_set=set([])
    for one_mnxid in Used_MNXid_set:
        Used_SMILES_set.add(cmpd_mnxid_smiles_dict[one_mnxid])
    print Used_SMILES_set

    pickle_out1=open("./Step03_screened_mnxid_rxn_list","wb")
    pickle.dump(screened_mnxid_rxn_list_2, pickle_out1)
    pickle_out1.close()

    pickle_out1=open("./Step03_Used_MNXid_set","wb")
    pickle.dump(Used_MNXid_set, pickle_out1)
    pickle_out1.close()

    pickle_out1=open("./Step03_Used_SMILES_set","wb")
    pickle.dump(Used_SMILES_set, pickle_out1)
    pickle_out1.close()

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
if __name__ == '__main__':
    main()



