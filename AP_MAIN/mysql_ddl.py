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

import pymysql.cursors
import time
# Connect to the database
connection = pymysql.connect(host='localhost',
                             user='xuzhiqin',
                             password='8'*8,
                             cursorclass=pymysql.cursors.DictCursor)
print ("connect to mysql db")
time.sleep(5)
print ("preparing")





cursor=connection.cursor()
DB_NAME = 'annealpath_20190501_for_volatility'
cursor.execute(
    "DROP DATABASE IF EXISTS "+DB_NAME+";")
cursor.execute(
    "CREATE DATABASE "+DB_NAME)
cursor.execute(
    "USE "+DB_NAME+";")
    

sql_line_1=(
    "CREATE TABLE IF NOT EXISTS search_result ("
    "  pwy_name varchar(50) NOT NULL,"
    "  sim_mtrc varchar(10) NOT NULL,"
    "  max_len int(11) NOT NULL,"
    "  ini_temp int(11) NOT NULL,"
    "  cool_coef double NOT NULL,"
    "  dsim_sel int(11) NOT NULL,"
    "  simu_sel int(11) NOT NULL,"
    "  simi_sel int(11) NOT NULL,"
    "  subs_sms varchar(200) NOT NULL,"
    "  prod_sms varchar(50) NOT NULL,"
    "  max_C_subs int(11) NOT NULL,"
    "  max_C_prod int(11) NOT NULL,"
    "  max_O_prod int(11) NOT NULL,"
    "  max_O_subs int(11) NOT NULL,"
    "  dsim_adj double NOT NULL,"
    "  simu_adj double NOT NULL,"
    "  simi_adj double NOT NULL,"
    "  rxn_time double NOT NULL,"
    "  pwy_time double NOT NULL,"
    "  pwy_num int(11) NOT NULL,"
    "  shtst int(11) NOT NULL,"
    "  shtst_num int(11) NOT NULL,"
    "  PRIMARY KEY (pwy_name)"
    ") ENGINE=InnoDB")

'''
try:
    with connection.cursor() as cursor:
        cursor.execute(sql_line_1,)
    connection.commit()

finally:
    connection.close()
    '''

try:
    with connection.cursor() as cursor:
        # Create a new record
        sql_line_2 = "INSERT INTO `search_result` (pwy_name,sim_mtrc,max_len,ini_temp,cool_coef,dsim_sel,simu_sel,simi_sel,\
                        subs_sms,prod_sms,max_C_subs,max_C_prod,max_O_prod,max_O_subs,dsim_adj,simu_adj,simi_adj,rxn_time,pwy_time,pwy_num,shtst,shtst_num\
                        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
        val_2=("test_pwy","ECFP4",4,1000,0.8,1,1,1,"CCO.CCCCO","CCCCCCO",2,2,2,2,0.5,0.5,2,9876.5,432.1,66666)
        cursor.execute(sql_line_1,)
        #cursor.execute(sql_line_2,val_2)
    # connection is not autocommit by default. So you must commit to save
    # your changes.
    connection.commit()

finally:
    connection.close()
    

