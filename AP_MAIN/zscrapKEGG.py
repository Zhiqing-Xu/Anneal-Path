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
if os.name == 'nt' or platform == 'win32':
    print("Running on Windows")
    if 'ptvsd' in sys.modules:
        print("Running in Visual Studio")
        try:
            os.chdir(os.path.dirname(__file__))
            print('CurrentDir: ', os.getcwd())
        except:
            pass
#--------------------------------------------------#
    else:
        print("Running outside Visual Studio")
        try:
            if not 'workbookDir' in globals():
                workbookDir = os.getcwd()
                print('workbookDir: ' + workbookDir)
                os.chdir(workbookDir)
        except:
            pass
#--------------------------------------------------#
if os.name != 'nt' and platform != 'win32':
    print("Not Running on Windows")
#--------------------------------------------------#

# ==================================================================================== #
import pickle
# ==================================================================================== #
# bioservices.KEGG
from bioservices import * # import KEGG
# pubchempy
import pubchempy
# chemspipy
from chemspipy import *
cs = ChemSpider("85cfd898-cc63-4347-9dec-0b4964a387c6")
# cirpy
import cirpy
# ==================================================================================== #
from AP_convert import *
import csv
from AP_convert import canonical_smiles_AP as unis
from bioservices import KEGG


# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
# 001 - Use saved KEGG_DBlinks, format: dict("KEGG_id" : ("CAS", "ChEBI", "PubChem"))
pickle_in1=open("./KEGGScrapSavings/KEGG_DBlinks_dict.pickle","rb")
KEGG_DBlinks_dict1=pickle.load(pickle_in1)
pickle_in1.close()
print (KEGG_DBlinks_dict1)

# ==================================================================================== #
# 002 - bioservices.ChEBI(), converting ChEBI to SMILES
def ChEBI_2_SMILES(ChEBI_a):
    ch=ChEBI()
    res = ch.getCompleteEntity("CHEBI:"+ChEBI_a)
    #print res.smiles
    return res.smiles

# ==================================================================================== #
# 003 - pubchempy, converting smiles to PubChem_id???
def smiles_2_PubChem_id(SMILES_a):
    # Example SMILES_a='Cn1cnc2n(C)c(=O)n(C)c(=O)c12'
    a=pubchempy.get_compounds(SMILES_a, 'smiles')
    return a[0]

# ==================================================================================== #
# 004 - pubchempy, converting PubChem_id to molecule info (including SMILES)
# properties_str_list=['bond_stereo_count', 'defined_atom_stereo_count', 'mmff94_partial_charges_3d', 'inchi', 'h_bond_donor_count', 'feature_selfoverlap_3d', 'canonical_smiles', 'shape_fingerprint_3d', 'isotope_atom_count', 'molecular_weight', 'coordinate_type', 'charge', 'conformer_rmsd_3d', 'isomeric_smiles', 'exact_mass', 'rotatable_bond_count', 'xlogp', 'defined_bond_stereo_count', 'iupac_name', 'monoisotopic_mass', 'tpsa', 'volume_3d', 'inchikey', 'elements', 'bonds', 'mmff94_energy_3d', 'conformer_id_3d', 'atoms', 'fingerprint', 'covalent_unit_count', 'shape_selfoverlap_3d', 'undefined_atom_stereo_count', 'cid', 'cactvs_fingerprint', 'pharmacophore_features_3d', 'effective_rotor_count_3d', 'record', 'complexity', 'heavy_atom_count', 'undefined_bond_stereo_count', 'h_bond_acceptor_count', 'molecular_formula', 'atom_stereo_count', 'multipoles_3d']
def PubChem_id_2_info_dict(PubChem_id):
    c = pubchempy.Compound.from_cid(PubChem_id)
    return c.to_dict().keys()
    # return c.to_dict(properties=properties_str.lower().split(', ')) # cannot remember why return this before

# ==================================================================================== #
# 005 - pubchempy, converting PubChem_id to SMILES
# properties_str_list=['bond_stereo_count', 'defined_atom_stereo_count', 'mmff94_partial_charges_3d', 'inchi', 'h_bond_donor_count', 'feature_selfoverlap_3d', 'canonical_smiles', 'shape_fingerprint_3d', 'isotope_atom_count', 'molecular_weight', 'coordinate_type', 'charge', 'conformer_rmsd_3d', 'isomeric_smiles', 'exact_mass', 'rotatable_bond_count', 'xlogp', 'defined_bond_stereo_count', 'iupac_name', 'monoisotopic_mass', 'tpsa', 'volume_3d', 'inchikey', 'elements', 'bonds', 'mmff94_energy_3d', 'conformer_id_3d', 'atoms', 'fingerprint', 'covalent_unit_count', 'shape_selfoverlap_3d', 'undefined_atom_stereo_count', 'cid', 'cactvs_fingerprint', 'pharmacophore_features_3d', 'effective_rotor_count_3d', 'record', 'complexity', 'heavy_atom_count', 'undefined_bond_stereo_count', 'h_bond_acceptor_count', 'molecular_formula', 'atom_stereo_count', 'multipoles_3d']
def PubChem_id_2_SMILES(PubChem_id):
    c = pubchempy.Compound.from_cid(PubChem_id)
    info_dict=c.to_dict()
    return info_dict['canonical_smiles']
    # return c.to_dict(properties=properties_str.lower().split(', ')) # cannot remember why return this before

#print PubChem_id_2_SMILES('1234')
# ==================================================================================== #
# 006 - cirpy, converting CAS to SMILES
def CAS_2_SMILES(CAS_a):
    return cirpy.resolve(CAS_a, 'smiles')

# ==================================================================================== #
# 006 - pubchempy, converting SMILES to canonical, non-big-pi SMILES
def SMILES_2_canonical_SMILES(SMILES_a):
    temp_a = pubchempy.get_properties('CanonicalSMILES', SMILES_a , 'smiles')
    a= temp_a[0]['CanonicalSMILES']
    return a.encode('ascii','ignore')
    
# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
# Use 001 and 002 to convert KEGG_id to SMILES
def KEGG_id_2_SMILES_1(KEGG_id_a):
    ChEBI_a=KEGG_DBlinks_dict1[KEGG_id_a][1][0]
    return ChEBI_2_SMILES(ChEBI_a)

# Use 001 and 006 to convert KEGG_id to SMILES
def KEGG_id_2_SMILES_2(KEGG_id_a):
    
    CAS_a=KEGG_DBlinks_dict1[KEGG_id_a][0][0]
    return CAS_2_SMILES(CAS_a)


# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
# Script used to scrap KEGG_DBlinks_dict (do not run!!)
def parse_KEGG_text_API(KEGG_id, attribute): 
# to (mainly) get the CAS and ChEBI of the compound with known KEGG_id
# attribute can only be one string type variable in the list below,
# [FORMULA,EXACT_MASS,MOL_WEIGHT,CAS,PubChem,ChEBI,ChEMBL,KNApSAcK,PDB-CCD,3DMET,NIKKAJI,ATOM]
# NAME,REACTION,PATHWAY,ENZYME will not be parsed for now

    KEGG_API_url="http://rest.kegg.jp/get/"+KEGG_id
    html=urlopen(KEGG_API_url)
    KEGG_cmpd_info=html.read()
    if KEGG_cmpd_info=='':
        print ('Cannot find this KEGG_id!')
        return 'Not Exist'
    KEGG_cmpd_info_list=KEGG_cmpd_info.split('\n')

    if attribute not in ['FORMULA','EXACT_MASS','MOL_WEIGHT','CAS','PubChem','ChEBI','ChEMBL','KNApSAcK','PDB-CCD','3DMET','NIKKAJI','ATOM']:
        print ('Attribute ERROR!')
        return
    if attribute == 'CAS':
        attribute = 'DBLINKS     CAS:'
    if attribute in ['PubChem','ChEBI','ChEMBL','KNApSAcK','PDB-CCD','3DMET','NIKKAJI']:
        attribute=attribute+':'


    for str_info in KEGG_cmpd_info_list:
        if str_info.find(attribute)!=-1:
            str_info=str_info.replace(attribute,'')
            str_info=str_info.replace('DBLINKS','')
            str_info_set=set(str_info.split(' '))
            str_info_set.discard('')
            str_info_list=tuple(str_info_set)
                

            return str_info_list
    #print 'Cannot find this attribute!'
    return 'None'

# ==================================================================================== #
def generate_KEGG_DBlinks_dict():
    def parse_KEGG_DBlinks(KEGG_id):
        DBlinks=[]
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'CAS'))
        if DBlinks[0]=='Not Exist':
            return DBlinks
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'ChEBI'))
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'PubChem'))
        return DBlinks

    '''
    All_KEGG_IDs_list=initialize_KEGG_id()
    KEGG_DBlinks_dict=dict([])
    for KEGG_id in All_KEGG_IDs_list:
        try:
            DBlinks=parse_KEGG_DBlinks(KEGG_id)
        except Exception as error_message:
            print 'Error raised when parsing', KEGG_id, ': ', error_message, '!!!'
        if DBlinks==['Not Exist']:
            continue
        KEGG_DBlinks_dict[KEGG_id]=DBlinks
        print KEGG_id, DBlinks
    print KEGG_DBlinks_dict

    pickle_out1=open("KEGG_DBlinks_dict.pickle","wb")
    pickle.dump(KEGG_DBlinks_dict, pickle_out1)
    pickle_out1.close()
    '''
    pickle_in1=open("KEGG_DBlinks_dict.pickle","rb")
    KEGG_DBlinks_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    return KEGG_DBlinks_dict1

# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
def KEGG_id_SMILES_txt(): # Step 01

    print(KEGG_id_2_SMILES_2("C00001"))

    # 1. Write trfm_odict.txt
    text_file = open("./KEGGScrapSavings/KEGG_id_SMILES.csv", "w")


    print(len(KEGG_DBlinks_dict1.keys()))
    bad_smiles_1=[]
    for one_KEGGid in KEGG_DBlinks_dict1.keys():
        '''
        try:
            smiles_1=KEGG_id_2_SMILES_1(one_KEGGid)
        except Exception:
            smiles_1="NAN"
            '''

        try:
            smiles_2=KEGG_id_2_SMILES_2(one_KEGGid)
        except Exception:
            smiles_2="NAN"
        print (one_KEGGid, smiles_2)

        text_file.write(one_KEGGid)
        #text_file.write(" , ")
        #text_file.write(smiles_1)
        text_file.write(" , ")
        text_file.write(str(smiles_2))
        text_file.write("\n")
    text_file.close()
    return

# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
def KEGG_id_smiles_canonical(): # Step 02 (ignorable runtime)
    import csv

    KEGG_id_list=[]
    smiles_canoncical_list=[]

    with open('./KEGGScrapSavings/KEGG_id_SMILES_1.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            KEGG_id=row[0][:-1]
            smiles_string=row[1][1:]
            print (KEGG_id)
            KEGG_id_list.append(KEGG_id)
            smiles_canoncial=unique_canonical_smiles_AP(smiles_string)
            smiles_processed=smiles_processed=smiles_canoncial.split(".")[0]
            smiles_processed=smiles_processed.replace("([*])", "")
            smiles_processed=smiles_processed.replace("[*]", "")
            smiles_processed=smiles_processed.replace("[O-]", "O")
            smiles_processed=smiles_processed.replace("[N+]", "N")
            smiles_processed=smiles_processed.replace("[NH+]", "N")
            smiles_processed=smiles_processed.replace("[NH2+]", "N")
            smiles_processed=smiles_processed.replace("[N-]", "N")
            smiles_processed=unique_canonical_smiles_AP(smiles_processed)
            smiles_canoncical_list.append(smiles_processed)
            print (smiles_processed)

    text_file = open("./KEGGScrapSavings/KEGG_id_canonical_SMILES.csv", "w")
    for i in range(len(KEGG_id_list)):
        text_file.write(KEGG_id_list[i])
        #text_file.write(" , ")
        #text_file.write(smiles_1)
        text_file.write(",")
        text_file.write(smiles_canoncical_list[i])
        text_file.write("\n")

    text_file.close()

    return
# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
def KEGGid_nme(): # Step 03

    s = KEGG()

    #print str(s.get("C05382"))

    All_KEGG_IDs_list=initialize_KEGG_id()
    text_file = open("./KEGGScrapSavings/KEGG_id_nme_1.csv", "w")

    for one_id in All_KEGG_IDs_list:
        KEGG_cmpd_info_bioservices=str(s.get(one_id))
        kcib_list=KEGG_cmpd_info_bioservices.split('\n')

        if KEGG_cmpd_info_bioservices.find("NAME")!=-1:
            for one_row in kcib_list:
                if one_row.find("NAME")!=-1:
                    str_info=one_row.replace("NAME",'')
                    str_info=str_info.replace("        ",'')
                    nme = str_info.replace(";",'')
                    break
        else:
            nme = "NAN"

        print (one_id, nme)
        text_file.write(one_id)
        text_file.write(",")
        text_file.write(str(nme))
        text_file.write("\n")
    text_file.close()

    return
# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
def KEGG_nme_smiles(): # Step 04 include step 02


    KEGG_id_list=[]
    smiles_canoncical_list=[]

    with open('./KEGGScrapSavings/KEGG_id_SMILES_1.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            KEGG_id=row[0][:-1]
            smiles_string=row[1][1:]
            KEGG_id_list.append(KEGG_id)
            smiles_canoncial=unique_canonical_smiles_AP(smiles_string)
            smiles_processed=smiles_processed=smiles_canoncial.split(".")[0]
            smiles_processed=smiles_processed.replace("([*])", "")
            smiles_processed=smiles_processed.replace("[*]", "")
            smiles_processed=smiles_processed.replace("[O-]", "O")
            smiles_processed=smiles_processed.replace("[N+]", "N")
            smiles_processed=smiles_processed.replace("[NH+]", "N")
            smiles_processed=smiles_processed.replace("[NH2+]", "N")
            smiles_processed=smiles_processed.replace("[N-]", "N")
            smiles_processed=smiles_processed.replace("*-", "N")
            smiles_processed=smiles_processed.replace("-*", "N")
            smiles_processed=smiles_processed.replace("(-*)", "N")
            smiles_processed=unique_canonical_smiles_AP(smiles_processed)
            smiles_canoncical_list.append(smiles_processed)


    KEGG_id_nme_dict=dict([])
    KEGG_nme_id_dict=dict([])
    with open('./KEGGScrapSavings/KEGG_id_nme_1.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile)
        
        for row in spamreader:
            row_str=row[0]
            for i in range(len(row)-1):
                row_str=row_str+","+row[i+1]
            row_str=row_str.replace(",",";",1)

            KEGG_id=row_str.split(";")[0]
            nme_string=row_str.split(";")[1]
            KEGG_id_nme_dict[KEGG_id]=nme_string
            KEGG_nme_id_dict[nme_string]=KEGG_id


    KEGG_nme_canonical_SMILES_dict=dict([]) # {key: smiles; value: nme}
    text_file = open("./KEGGScrapSavings/z4_KEGG_nme_canonical_SMILES.csv", "w") #z4 means 4th step 

    improper_cmpd_nmes={"Ketone", "Ester", "Fatty acid", "Calcium formate",}
    #"D-Xylulose 5-phosphate", "L-Xylulose 5-phosphate"
    for i in range(len(KEGG_id_list)):
        print (KEGG_id_list[i], smiles_canoncical_list[i], KEGG_id_nme_dict[KEGG_id_list[i]])
        if smiles_canoncical_list[i]!="NAN" and KEGG_id_list[i] in KEGG_id_nme_dict.keys():
            KEGG_nme=KEGG_id_nme_dict[KEGG_id_list[i]]
            KEGG_smiles_processed=smiles_canoncical_list[i]
            if KEGG_smiles_processed not in KEGG_nme_canonical_SMILES_dict.keys():
                KEGG_nme_canonical_SMILES_dict[KEGG_smiles_processed]=[KEGG_nme, KEGG_id_list[i]]
            else:
                for one_improper_nme in improper_cmpd_nmes:
                    if KEGG_nme_canonical_SMILES_dict[KEGG_smiles_processed][0] == one_improper_nme and KEGG_nme != one_improper_nme:
                        KEGG_nme_canonical_SMILES_dict[KEGG_smiles_processed]=[KEGG_nme, KEGG_nme_id_dict[KEGG_nme]]

                if KEGG_nme not in improper_cmpd_nmes:
                    if len(KEGG_nme)<len(KEGG_nme_canonical_SMILES_dict[KEGG_smiles_processed][0]):
                        KEGG_nme_canonical_SMILES_dict[KEGG_smiles_processed]=[KEGG_nme, KEGG_nme_id_dict[KEGG_nme]]

            text_file.write(KEGG_smiles_processed)
            text_file.write(",")
            text_file.write(KEGG_nme.replace(",",";"))
            text_file.write(",")
            text_file.write(KEGG_id_list[i])
            text_file.write(",")
            text_file.write("!" + KEGG_nme.replace(",",";") + "!")
            text_file.write("\n")

            #print KEGG_smiles_processed, KEGG_nme, KEGG_id_list[i]

    text_file.close()


    pickle_out1=open("./KEGGScrapSavings/KEGG_nme_canonical_SMILES_dict.pickle","wb")
    pickle.dump(KEGG_nme_canonical_SMILES_dict, pickle_out1)
    pickle_out1.close()

    return


# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
def Modify_z4_dict(): # Step 05 

    KEGG_id_nme_dict=dict([])
    KEGG_nme_id_dict=dict([])
    with open('./KEGGScrapSavings/KEGG_id_nme_1.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile)
        
        for row in spamreader:
            row_str=row[0]
            for i in range(len(row)-1):
                row_str=row_str+","+row[i+1]
            row_str=row_str.replace(",",";",1)
            KEGG_id=row_str.split(";")[0]
            nme_string=row_str.split(";")[1]
            KEGG_id_nme_dict[KEGG_id]=nme_string
            KEGG_nme_id_dict[nme_string]=KEGG_id

    # ==================================================================================== #
    pickle_in1=open("./KEGGScrapSavings/KEGG_nme_canonical_SMILES_dict.pickle","rb")
    KEGG_nme_canonical_SMILES_dict=pickle.load(pickle_in1)
    pickle_in1.close()

    # if smiles is in the dict but name in the dict is wrong, use this list to fix.
    smiles_nme_list=[["O=CCO","Glycolaldehyde"],
                     ["CC(=O)O","Acetate"],
                     ["O=C(CO)C(O)C(O)COP(=O)(O)O", "D-Xylulose 5-phosphate"],
                     ["OCCO", "Ethylene glycol"],
                     ["CN", "Methylamine"], 
                     ["O=P(O)(O)OCC1OC(O)(CO)C(O)C1O", "beta-D-Fructose 6-phosphate"],
                     ["O=CC(=O)O", "Glyoxylate"],
                     ["NC(Cc1c[nH]c2ccccc12)C(=O)O", "L-Tryptophan"],
                     ["CC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O","2 trans,trans-Farnesyl diphosphate"],
                     ["OCC1OC(O)C(O)C(O)C1O","D-Glucose"],
                     ]


    for smiles_nme in smiles_nme_list:
        one_smiles=smiles_nme[0]
        one_nme=smiles_nme[1]
        KEGG_nme_canonical_SMILES_dict[one_smiles]=[one_nme, KEGG_nme_canonical_SMILES_dict[one_smiles][1] if one_nme not in KEGG_nme_id_dict.keys() else KEGG_nme_id_dict[one_nme]]
    # ==================================================================================== #
    non_KEGG_list=[["smiles","cmpd_nme", "non_KEGG_id"],

                     ]


    # if smiles is not in the dict, use this list to add smiles:nme,id to the dict
    non_KEGG_list=[
                    ["OCCCCO",'butane-1,4-diol', 'MNXM16407']  , #2
                    ["O=CCCCO","4-hydroxybutanal", "MNXM37518"]  ,
                    ["CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCSC(=O)CCCO",'4-hydroxybutanoyl-CoA', 'MNXM2334' ]  ,
                    ["O=C(O)CCCO",'Sodium oxybate (USAN)', 'MNXM161191']  ,
                    ["O=CCCC(=O)O",'succinate semialdehyde', 'MNXM172']  ,

                    ['NC(CO)CO', '2-amino-1,3-propanediol', 'MNXM21729']  , #3
                    ['NC(CO)COP(=O)(O)O', 'serinol phosphate', 'MNXM12933']  ,
                    ['NC(CCC(=O)O)C(=O)O', 'disodium L-glutamate', 'MNXM511523']  ,
                    ['O=C(CO)COP(=O)(O)O', 'dihydroxyacetone phosphate', 'MNXM77']  ,

                    ['N=C(O)CCCCN', '5-Aminopentanamide', 'PubChemPy']  , #8

                    ['CC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O', '2-trans,6-cis-farnesyl diphosphate', 'MNXM35463']  , #16
                    ['CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC=C(C)C', 'squalene', 'MNXM292']  ,
                    ['CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC1OC1(C)C', '2,3-epoxysqualene', 'MNXM468658']  ,
                    ['CC(C)=CCCC(C)(O)C1CCC2(C)C1CCC1C3(C)CCC(O)C(C)(C)C3CCC12C', 'Dammarenediol-I', 'MNXM724543']  ,
                    ['CC(C)=CCCC(C)(O)C1CCC2(C)C1C(O)CC1C3(C)CCC(O)C(C)(C)C3CCC12C', '(20R)-protopanaxadiol', 'MNXM723118']  ,

                    ["O=C(O)c1ccc(CO)cc1","4-(Hydroxymethyl)benzoate","PubChemPy"]  , #18
                    ['O=Cc1ccc(C(=O)O)cc1', '4-carboxybenzaldehyde', 'MNXM10176']  ,

                    ['NC(Cc1c[nH]c2ccccc12)C(=O)O', 'L-Tryptophan', 'C00078']  , #20
                    ['N=C(Cc1c[nH]c2ccccc12)C(=O)O', 'IPA imine', 'C21124']  ,
                    ['O=C(O)c1[nH]c(-c2c[nH]c3ccccc23)cc1-c1c[nH]c2ccccc12', 'protodeoxyviolaceinate', 'C21131']  ,
                    ['O=C(O)c1[nH]c(-c2c[nH]c3ccc(O)cc23)cc1-c1c[nH]c2ccccc12', 'protoviolaceinate', 'C21134']  ,
                    ['O=C(O)c1[nH]c(-c2c[nH]c3ccc(O)cc23)cc1-c1c(O)[nH]c2ccccc12', 'violaceinate', 'C21135']  ,
                    ['O=C1NC(c2c[nH]c3ccc(O)cc23)=CC1=C1C(=O)Nc2ccccc21', 'violacein', 'C21136']  ,

                     ]


    
    for one_list in non_KEGG_list:
        if one_list[0] not in KEGG_nme_canonical_SMILES_dict.keys():
            KEGG_nme_canonical_SMILES_dict[one_list[0]]=[one_list[1], one_list[2]]


    # ==================================================================================== #
    pickle_out1=open("./KEGGScrapSavings/KEGG_nme_canonical_SMILES_dict.pickle","wb")
    pickle.dump(KEGG_nme_canonical_SMILES_dict, pickle_out1)
    pickle_out1.close()
    print ("done")

def Test_01():
    def read_smiles(dictionary, one_smiles):
        try:    
            print ([unis(one_smiles),] + dictionary[ unis( one_smiles ) ]," ,")
        except Exception:
            print ("!!!!! Cannot Find : ", unis(one_smiles), " !!!!!" )
        return 

    pickle_in1=open("./KEGGScrapSavings/KEGG_nme_canonical_SMILES_dict.pickle","rb")
    KEGG_nme_canonical_SMILES_dict=pickle.load(pickle_in1)
    pickle_in1.close()
    d1=KEGG_nme_canonical_SMILES_dict
    #print unis("CC(=O)I")
    
    print ("start read smiles")
    # ==================================================================================== #
    #2
    print("#2")
    print (KEGG_nme_canonical_SMILES_dict[ unis("OCCCCO") ])
    print (KEGG_nme_canonical_SMILES_dict[ unis("O=CCCCO") ])
    print (KEGG_nme_canonical_SMILES_dict[ unis("CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCSC(=O)CCCO") ])
    print (KEGG_nme_canonical_SMILES_dict[ unis("O=C(O)CCCO") ])
    print (KEGG_nme_canonical_SMILES_dict[ unis("O=CCCC(=O)O") ])
    
    # ==================================================================================== #
    #3
    print("#3")
    read_smiles(d1, "[H]OC([H])([H])C([H])(N([H])[H])C([H])([H])O[H]")
    read_smiles(d1, "[H]OC([H])([H])C([H])(N([H])[H])C([H])([H])OP(=O)(O[H])O[H]")
    read_smiles(d1, "[H]OC(=O)C([H])([H])C([H])([H])C([H])(C(=O)O[H])N([H])[H]")
    read_smiles(d1, "[H]OC([H])([H])C(=O)C([H])([H])OP(=O)(O[H])O[H]")
    read_smiles(d1, "O=C(O)CCC(=O)C(=O)O")
    
    # ==================================================================================== #
    print("#4")
    read_smiles(d1, "C=C(OC1C=CC=C(C(=O)O)C1O)C(=O)O")
    read_smiles(d1, "O=C(O)c1ccccc1O")
    read_smiles(d1, "O=C(O)c1cc(O)ccc1O")
    # ==================================================================================== #
    print("#5")
    read_smiles(d1, "[H]OC(=O)C(=O)C([H])([H])C1(C(=O)O[H])C([H])=C([H])C([H])(O[H])C([H])=C1[H]")
    read_smiles(d1, "O=C(O)C(=O)Cc1ccccc1")
    read_smiles(d1, "O=C(O)C(O)c1ccccc1")
    read_smiles(d1, "O=C(O)C(=O)c1ccccc1")
    read_smiles(d1, "O=Cc1ccccc1")
    read_smiles(d1, "OCc1ccccc1")
    # ==================================================================================== #
    print("#6")
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O")
    read_smiles(d1, "C=C(C)CCOP(=O)(O)OP(=O)(O)O")
    
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O")
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CC=CC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C")
    read_smiles(d1, "CC(C)=CCCC(C)=CC=CC(C)=CC=CC(C)=CC=CC=C(C)C=CC=C(C)C=CC=C(C)CCC=C(C)C")
    read_smiles(d1, "CC(C)=CCCC(C)=CC=CC(C)=CC=CC(C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CCCC1(C)C")
    read_smiles(d1, "CC(C=CC=C(C)C=CC1=C(C)CCCC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CCCC1(C)C")
    # ==================================================================================== #
    print("#7")
    read_smiles(d1, "C=C(OC1C=CC=C(C(=O)O)C1O)C(=O)O")
    read_smiles(d1, "O=C(O)c1ccccc1O")   
    read_smiles(d1, "Oc1ccccc1O")
    read_smiles(d1, "O=C(O)C=CC=CC(=O)O")
    # ==================================================================================== #
    print("#8")
    read_smiles(d1, "NCCCCC(N)C(=O)O")
    read_smiles(d1, "N=C(O)CCCCN")   
    read_smiles(d1, "NCCCCC(=O)O")
    read_smiles(d1, "O=CCCCC(=O)O")
    read_smiles(d1, "O=C(O)CCCC(=O)O")
    # ==================================================================================== #
    print("#9")
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O")
    read_smiles(d1, "C=C(C)CCOP(=O)(O)OP(=O)(O)O")   
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O")
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CC=CC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C")
    read_smiles(d1, "CC(C)=CCCC(C)=CC=CC(C)=CC=CC(C)=CC=CC=C(C)C=CC=C(C)C=CC=C(C)CCC=C(C)C")
    # ==================================================================================== #
    print("#10")
    read_smiles(d1, "NC(CCC(=O)O)C(=O)O")
    read_smiles(d1, "CC(C(=O)O)C(N)C(=O)O")   
    read_smiles(d1, "CC(=CC(=O)O)C(=O)O")
    # ==================================================================================== #
    print("#11")
    read_smiles(d1, "NC(Cc1ccc(O)cc1)C(=O)O")
    read_smiles(d1, "O=C(O)C=Cc1ccc(O)cc1")   
    read_smiles(d1, "C(=O)(C=Cc1ccc(O)cc1)(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
    #read_smiles(d1, "CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(O)=NCCC(O)=NCCSC(=O)C=Cc1ccc(O)cc1")
    read_smiles(d1, "O=C(C=Cc1ccc(O)cc1)c1c(O)cc(O)cc1O")
    read_smiles(d1, "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c21")
    # ==================================================================================== #
    print("#12")
    read_smiles(d1, "NCCCCN")
    read_smiles(d1, "CNCCCCN")   
    read_smiles(d1, "CNCCCC=O")
    read_smiles(d1, "CN1=CCCC1")
    # ==================================================================================== #
    print("#13")
    read_smiles(d1, "NC(Cc1ccc(O)cc1)C(=O)O")
    read_smiles(d1, "O=C(O)C=Cc1ccc(O)cc1")   
    read_smiles(d1, "C=Cc1ccc(O)cc1")
    # ==================================================================================== #
    print("#14")
    read_smiles(d1, "NC(Cc1ccc(O)cc1)C(=O)O")
    read_smiles(d1, "O=C(O)C=Cc1ccc(O)cc1")   
    read_smiles(d1, "C(=O)(C=Cc1ccc(O)cc1)(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
    read_smiles(d1, "Oc1ccc(C=Cc2cc(O)cc(O)c2)cc1")
    read_smiles(d1, "Oc1cc(O)cc(C=Cc2ccc(O)c(O)c2)c1")  

    # ==================================================================================== #
    print("#15")
    read_smiles(d1, "NC(Cc1ccccc1)C(=O)O")
    read_smiles(d1, "O=C(O)C=Cc1ccccc1")   
    read_smiles(d1, "C(=O)(C=Cc1ccccc1)(SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)")
    read_smiles(d1, "O=C(C=Cc1ccccc1)c1c(O)cc(O)cc1O")
    read_smiles(d1, "O=C1CC(c2ccccc2)Oc2cc(O)cc(O)c21") 
    # ==================================================================================== #
    print("#16")
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O")
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC=C(C)C")   
    read_smiles(d1, "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC1OC1(C)C")
    read_smiles(d1, "CC(C)=CCCC(C)(O)C1CCC2(C)C1CCC1C3(C)CCC(O)C(C)(C)C3CCC12C")
    read_smiles(d1, "CC(C)=CCCC(C)(O)C1CCC2(C)C1C(O)CC1C3(C)CCC(O)C(C)(C)C3CCC12C")  

    # ==================================================================================== #
    print("#17")
    read_smiles(d1, "NC(Cc1ccccc1)C(=O)O")
    read_smiles(d1, "O=C(O)C=Cc1ccccc1")   
    read_smiles(d1, "C=Cc1ccccc1")
    # ==================================================================================== #
    print("#18")
    read_smiles(d1, "Cc1ccc(C)cc1")
    read_smiles(d1, "Cc1ccc(CO)cc1")   
    read_smiles(d1, "Cc1ccc(C=O)cc1")
    read_smiles(d1, "Cc1ccc(C(=O)O)cc1")
    read_smiles(d1, "O=C(O)c1ccc(CO)cc1")   
    read_smiles(d1, "O=Cc1ccc(C(=O)O)cc1")
    read_smiles(d1, "O=C(O)c1ccc(C(=O)O)cc1")
    # ==================================================================================== #
    print("#19")
    read_smiles(d1, "O=C(O)C1=CC(=O)C(O)C(O)C1")
    read_smiles(d1, "O=C(O)c1ccc(O)c(O)c1")   
    read_smiles(d1, "O=Cc1ccc(O)c(O)c1")
    read_smiles(d1, "COc1cc(C=O)ccc1O")
    # ==================================================================================== #
    print("#20")
    read_smiles(d1, "NC(Cc1c[nH]c2ccccc12)C(=O)O") # NC(Cc1c[nH]c2ccccc12)C(=O)O
    read_smiles(d1, "N=C(Cc1c[nH]c2ccccc12)C(=O)O") # NC(=Cc1c[nH]c2ccccc12)C(=O)O
    #read_smiles(d1, "N=C(C(=O)O)C(c1c[nH]c2ccccc12)C(C(=N)C(=O)O)c1c[nH]c2ccccc12")
    read_smiles(d1, "O=C(O)c1[nH]c(-c2c[nH]c3ccccc23)cc1-c1c[nH]c2ccccc12")
    read_smiles(d1, "O=C(O)c1[nH]c(-c2c[nH]c3ccc(O)cc23)cc1-c1c[nH]c2ccccc12")    
    read_smiles(d1, "O=C(O)c1[nH]c(-c2c[nH]c3ccc(O)cc23)cc1-c1c(O)[nH]c2ccccc12")
    read_smiles(d1, "O=C1NC(c2c[nH]c3ccc(O)cc23)=CC1=C1C(=O)Nc2ccccc21")



    read_smiles(d1, "NCCN")



if __name__ == '__main__':
    #KEGG_nme_smiles()
    Modify_z4_dict()
    Test_01()
    #generate_KEGG_DBlinks_dict()



# ==================================================================================== #
# ==================================================================================== #
# ==================================================================================== #
# 004 - ChemSpider, converting CS_id to mol files (2D & 3D)
'''
cs = ChemSpider("PW2qJbvGaAK1gAthkR6n4jDxoUJ3xw2D")
c1 = cs.get_compound(2424)  # Specify compound by ChemSpider ID
c2 = cs.search('benzene')  # S

info = cs.get_extended_compound_info(2424)
print info
mol=cs.get_record_mol(2424, calc3d=True)
print mol
'''
# ==================================================================================== #
'''
c3=cs.search('CCCCC')
print 'caonima'
print c3
a=[]
b=[]
for result in c3:
    a.append(result.mol2D)
    b.append(result.mol3D)
print a
print b
'''
# ==================================================================================== #
# 004 - 
'''
from bioservices import *
s = ChemSpider("PW2qJbvGaAK1gAthkR6n4jDxoUJ3xw2D")
s.find("Pyridine")
results = s.GetExtendedCompoundInfo(1020)
print results['averagemass']
'''
# ==================================================================================== #
# 004 - 
'''
from chemspipy import *
cs = ChemSpider("PW2qJbvGaAK1gAthkR6n4jDxoUJ3xw2D")
c1 = cs.get_compound(2424)  # Specify compound by ChemSpider ID
c2 = cs.search('benzene')  # S
info = cs.get_extended_compound_info(2424)
print info
mol=cs.get_record_mol(2424, calc3d=True)
print mol
c3=cs.search('Cn1cnc2n(C)c(=O)n(C)c(=O)c12')
print 'caonima'
for result in c3:
    print(result.csid)
    '''
