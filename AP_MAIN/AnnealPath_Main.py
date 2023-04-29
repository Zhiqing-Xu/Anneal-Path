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
if __name__ == "__main__":
    print("="*80)
    if os.name == 'nt' or platform == 'win32':
        print("Running on Windows")
        if 'ptvsd' in sys.modules:
            print("Running in Visual Studio")
#--------------------------------------------------#
    if os.name != 'nt' and platform != 'win32':
        print("Not Running on Windows")
#--------------------------------------------------#
    if "__file__" in globals().keys():
        try:
            os.chdir(os.path.dirname(__file__))
            print('CurrentDir: ', os.getcwd())
        except:
            print("Problems with navigating to the file dir.")
    else:
        print("Running in python jupyter notebook.")
        try:
            if not 'workbookDir' in globals():
                workbookDir = os.getcwd()
                print('workbookDir: ' + workbookDir)
                os.chdir(workbookDir)
        except:
            print("Problems with navigating to the workbook dir.")
#--------------------------------------------------#
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
#--------------------------------------------------#
from AP_Solver import * 
from AP_funcs import multipletrialnames
#######################################################################################################################################
#######################################################################################################################################
# rxn_db            - reaction rule sets used to expand pathway network, "APrules", "BNICE", "retropath", "Simpheny".
# pathway_name      - name of the pathway defined by user.
# substrates        - a list of substrates as starting compounds (multiple compounds input eg.: ["CCO.CCC",]).
# products          - a list of products as target compounds (only one target compound in this example).
# max_levels        - maximum search levels for expanding the reaction network (pathway length upper bound).
# pruning_method    - different pruning methods.
#                   - [0,0,0],    - For each search stage, accept half number of compounds (recommended).
#                   - [2,2,2],    - For each search stage, accept all compounds found.
#                   - [x,x,x],    - Accept specified numbers of compounds from three bins (SI, delta(SI) and SA)
# fp_type           - fingerprint type for SimIndex or depth for MNA.
# fp_type           - "top" , "MACCS" , "atom_pairs" , "vec_pairs" , "torsions" , "ECFP" , "FCFP".
# max_value         - upper bound of size of intermediate compounds (used to screen too large compounds)
# bin_adj           - multipliers adjust the number of three bins
# plot_type         - plot only pathways with known compounds / plot all pathways
#######################################################################################################################################
#######################################################################################################################################
def main_AP(    rxn_db         =  "APrules",
                pathway_name   =  "sample pathways",
                substrates     =  ["C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O"], # G6P
                products       =  ["O=C(C(=O)O)C"], #PYR
                max_levels     =  6,
                pruning_method =  [0,0,0],
                fp_type        =  "ECFP",
                max_value      =  [0,0,0,0],
                bin_adj        =  [0.0,0.0,0.0],
                plot_type      =  "selected"
            ):
    #============================================================================================================================#
    APfunc = AnnealPath(outstream=sys.stderr, rxn_diameter=2, rxn_score_lb=0, rxn_db=rxn_db)
    print( "$"*35 + "  Computing pathway %s  " % pathway_name + "$"*35 + '\n')
    #============================================================================================================================#
    pathway_name = multipletrialnames(pathway_name, num_digits=4)
    pathway_r_files = "pwy_r" # R is used to efficiently plot mulitple hypergraphs representing pathways.
    #============================================================================================================================#
    # ! ! ! ! ! ! ! ! !: Go to AP_Solver.py and search! #
    APfunc.solve_pathway(substrates[0], products[0], 
                         pathway_name, pathway_r_files, 
                         max_levels, rxn_db, pruning_method, 
                         fp_type, max_value, bin_adj, plot_type)
    print('Done solving pathway!')
    #============================================================================================================================#
    return pathway_name
#######################################################################################################################################
#######################################################################################################################################
def main_test(  rxn_db         =  "APrules",
                pathway_name   =  "test_pathway",
                substrates     =  ["C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O"], # G6P
                products       =  ["O=C(C(=O)O)C"], #PYR
                max_levels     =  5,
                pruning_method =  [1,1,1],
                fp_type        =  "ECFP",
                max_value      =  [0,0,0,0],
                bin_adj        =  [0.0,0.0,0.0],
                plot_type      =  "selected",
                count_trial    =  1
             ):
    for one_trial in range(count_trial):
        main_AP(rxn_db,pathway_name,substrates,products,max_levels,pruning_method,fp_type,max_value,bin_adj,plot_type)
    return
#######################################################################################################################################
#######################################################################################################################################
def main():

    #============================================================================================================================#
    #main_test(  rxn_db = "retropath",
    #            pathway_name = "F6P to DAHP",
    #            substrates = ["C(C(C(C(C(=O)CO)O)O)O)OP(=O)(O)O"], # F6P
    #            products = ["C(C(C(C(COP(=O)(O)O)O)O)O)C(=O)C(=O)O"], #DAHP
    #            max_levels = 7,
    #            pruning_method = [0,0,0],
    #            fp_type = "ECFP",
    #            max_value = [0,0,0,0],
    #            bin_adj = [0.1,0.1,0.1],
    #            plot_type = "selected",
    #            count_trial = 1
    #         )
    #============================================================================================================================#
    pathway_result_2 = \
            main_AP(    rxn_db          =  "APrules",
                        pathway_name    =  "AKG to LYS",
                        substrates      =  ['O=C(O)C(=O)CCC(=O)O'],
                        products        =  ['NCCCCC(N)C(=O)O'],
                        max_levels      =  8,
                        pruning_method  =  [0,0,0],
                        fp_type         =  "ECFP",
                        max_value       =  [1,1,1,1],
                        bin_adj         =  [0.2,0.2,0.2],
                        plot_type       =  "selected"
                        )
    plot_result_path = Path("../results/" + pathway_result_2 + "/pwy_r/pathways/")
    Image(plot_result_path / 'pathway0.png')
    Image(plot_result_path / 'pathway4.png')


    #============================================================================================================================#
    #pathway_result_2 = \
    #main_AP(    rxn_db          =  "APrules",
    #            pathway_name    =  "EtOHBtOH to HxOH",
    #            substrates      =  ["CCO.CCCCO"],
    #            products        =  ["CCCCCCO"],
    #            max_levels      =  6,
    #            pruning_method  =  [2,2,2],
    #            fp_type         =  "ECFP",
    #            max_value       =  [3,3,3,3],
    #            bin_adj         =  [0,0,0],
    #            plot_type       =  "novel"
    #            )
    #============================================================================================================================#
    #main_AP(    rxn_db         =  "APrules",
    #            pathway_name   =  'Glucose to AcCoA',
    #            substrates     =  ['C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O'],
    #            products       =  ['CC(=O)CoA'],
    #            max_levels     =  6,
    #            pruning_method =  [0,0,0],
    #            fp_type        =  "ECFP",
    #            max_value      =  [0,0,0,0],
    #            bin_adj        =  [0.2,0.2,0.2],
    #            plot_type      =  "selected",
    #            )
    #============================================================================================================================#
    #main_AP(    rxn_db         =  "APrules",
    #            pathway_name   =  'AKG to LYS ',
    #            substrates     =  ['O=C(O)C(=O)CCC(=O)O'],
    #            products       =  ['NCCCCC(N)C(=O)O'],
    #            max_levels     =  8,
    #            pruning_method =  [2,2,2],
    #            fp_type        =  "ECFP",
    #            max_value      =  [2,2,2,2],
    #            bin_adj        =  [0.1,0.1,0.1],
    #            plot_type      =  "selected",
    #            )
    #============================================================================================================================#
    return

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
if __name__ == '__main__':
    main()

#######################################################################################################################################
#######################################################################################################################################
# Amino Acids
    # ARG                                  C(C(C(=O)O)N)C(=O)N
    # HIS                                  C1=C(NC=N1)CC(C(=O)O)N
    # LYS                                  C(CCN)CC(C(=O)O)N
    # ASP                                  O=C(O)CC(N)C(=O)O
    # GLU                                  C(CC(=O)O)C(C(=O)O)N
    # SER                                  C(C(C(=O)O)N)O
    # THR                                  CC(C(C(=O)O)N)O
    # ASN                                  C(C(C(=O)O)N)C(=O)N
    # GLN                                  C(CC(=O)N)C(C(=O)O)N
    # PRO                                  C1CC(NC1)C(=O)O
    # alanine                              CC(C(=O)O)N
    # valine                               CC(C)C(C(=O)O)N
    # isoleucine                           CCC(C)C(C(=O)O)N
    # leucine                              CC(C)CC(C(=O)O)N
    # glutamate                            C(CC(=O)O)C(C(=O)O)N
    # glutamine                            C(CC(=O)N)C(C(=O)O)N
    # proline                              C1CC(NC1)C(=O)O
    # arginine                             C(CC(C(=O)O)N)CN=C(N)N
    # oxaloacetate                         O=C(O)C(=O)CC(=O)O
    # aspartate                            O=C(O)CC(N)C(=O)O
    # asparagine                           O=C(N)CC(N)C(=O)O
    # threonine                            CC(C(C(=O)O)N)O
    # lysine                               C(CCN)CC(C(=O)O)N
    # glutamate                            C(CC(=O)O)C(C(=O)O)N
    # aspartate                            O=C(O)CC(N)C(=O)O
    # histidine                            C1=C(NC=N1)CC(C(=O)O)N
    # serine                               C(C(C(=O)O)N)O
    # glycine                              C(C(=O)O)N
    # methionine                           CSCCC(C(=O)O)N
    # cysteine                             C(C(C(=O)O)N)S
    # phenylalanine                        C1=CC=C(C=C1)CC(C(=O)O)N
    # Tyrosine                             C1=CC(=CC=C1CC(C(=O)O)N)O
    # Tryptophan                           C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N
#######################################################################################################################################
#######################################################################################################################################
# Common Central Carbon Metabolites
    # oxaloacetate (oxaloacetic acid)      O=C(O)C(=O)CC(=O)O
    # pyruvate                             O=C(C(=O)O)C
    # alpha-ketoglutarate                  O=C(O)C(=O)CCC(=O)O
    # G6P                                  C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O
    # PEP                                  O=C(O)C(OP(=O)(O)O)=C
    # E4P                                  C(C(C(C=O)O)O)OP(=O)(O)O
    # R5P                                  C(C(C(C(C=O)O)O)O)OP(=O)(O)O
    # F6P                                  C(C(C(C(C(=O)CO)O)O)O)OP(=O)(O)O
    # DAHP                                 C(C(C(C(COP(=O)(O)O)O)O)O)C(=O)C(=O)O
    # PYU                                  O=C(C(=O)O)C                            pyruvate
    # OXA                                  O=C(O)C(=O)CC(=O)O                      oxaloacetate
    # SUC                                  C(CC(=O)O)C(=O)O                        succinate
    # AKG                                  O=C(O)C(=O)CCC(=O)O                     alpha-ketoglutarate, 2-OG
#######################################################################################################################################
#######################################################################################################################################