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
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import *

import os, os.path
from sys import platform
if os.name == 'nt' or platform == 'win32':
    os.chdir(os.path.dirname(__file__))


#import matplotlib.pyplot as plt
from PIL import Image
import PIL

def replace68(str_x):
    str_x=str_x.replace("#6","C")
    return str_x.replace("#8","O")

'''
a="[#1]-[#6]=[#6]"
molecule_a=Chem.MolFromSmiles(a)

b="[#1]-[#6]=[#6]"
molecule_b=Chem.MolFromSmarts(b)
'''



rxn="([#6](=[#8])-[#6]-[#8])>>([#6](-[#8])-[#6]=[#8])"

rxn_1 = replace68(rxn)
print (rxn_1)
print (rxn_1.split(">>")[1]+">>"+rxn_1.split(">>")[0])



#rxn="([O:2]-[C:3](=[O:4])-[C:5](=[O:6])-[C:7](-[O:9])-[C:11](-[C:15](=[O:16])-[C:17]-[O:20])-[O:13])>>([C:5](-[C:3](-[O:2])=[O:4])(-[C:7](-[O:9])-[C:11](-[O:13])-[C:15](-[C:17]-[O:20])=[O:16])-[O:6])"


#rxn="([CH:5](-[C:3](-[O:2])=[O:4])(-[C:7](-[O:9])-[C:11](-[O:13])-[C:15](-[C:17]-[O:20])=[O:16])-[OH:6])>>([O:2]-[C:3](=[O:4])-[C:5](=[O:6])-[C:7](-[O:9])-[C:11](-[C:15](=[O:16])-[C:17]-[O:20])-[O:13])"



a=rxn.split(">>")[0][1:-1]
molecule_a=Chem.MolFromSmarts(a)
b=rxn.split(">>")[1][1:-1]
molecule_b=Chem.MolFromSmarts(b)


#img=Draw.MolsToGridImage(molecule,ImgSize=(200,200))
ms=[Chem.MolFromSmarts(a),Chem.MolFromSmarts(b)]
img=Draw.MolsToGridImage(ms,molsPerRow=2,subImgSize=(500,500),legends=[Chem.MolToSmarts(x) for x in ms])
#img.save('rules/MNX_fwd/'+i+'.png')
#img.save('rules/APrules/'+i+'.png')
img.save('temp_ab'+'.png')


'''
Chem.Draw.MolToFile(molecule_a,'cdk2_mol1.o.png',size=(500,500)) # depending on molecule size, shall adjust the img size input
'''

os.system('temp_ab'+'.png')
print (Chem.MolToSmarts(molecule_a))
print (Chem.MolToSmarts(molecule_b))



#rxn="([C:2](-[C:5](-[C:7](-[O:8])=[O:9])=[O:6])(-[O:4])-[C:3])>>([C:5](-[C:7](-[O:8])=[O:9])(-[C:2](-[O:4])-[C:3])-[O:6])"
#rxn="([CH:5](-[C:7](-[O:8])=[O:9])(-[C:2](-[O:4])-[C:3])-[OH:6])>>([C:2](-[C:5](-[C:7](-[O:8])=[O:9])=[O:6])(-[O:4])-[C:3])"



