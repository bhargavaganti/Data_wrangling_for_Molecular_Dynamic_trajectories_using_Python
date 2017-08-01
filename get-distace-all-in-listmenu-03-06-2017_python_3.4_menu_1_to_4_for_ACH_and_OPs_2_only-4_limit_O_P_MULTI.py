## This program is indent to get coordinates of specified atoms and calculate
## their inter atomic distances. The out put will be stored in a separate file.


## import necessary packages
import os
import glob

import math
import re
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
#import multiprocessing
#from multiprocessing import 
from multiprocessing import Process, Queue, Pool
import struct
import gc

gc.collect()
## global variables
#file_names = []
#from biopandas.pdb import PandasPDB

############################################################################################
## Function to select main manu items to implement
def menu_list():

	select1 = False
	select2 = False
	select3 = False
	select4 = False
	select5 = False

	menu_itm = input('''Main menu..
		1: Write ligand protein and ligand water distances to .csv files seperately
		2: Give important distances vs score
		3: Analyze important distance vs score
		4: Read and analyze MD trajectories

Enter your selection number/s (comma, tab or space seperation
		''')

	#menu_itm = input(")")
	items =re.findall(r"[\w']+", menu_itm) ## seperate the input based on deliminators

	for itm in items:
		if itm == '1':
			select1 = True
		elif itm == '2':
			select2 = True
			'''

			print(Enter the atom pairs to measure distances:
			default distances:

			ser_O = "3038"
			his_N1 = "6619"  ## N originally without H
			his_N2 = "6616"  ## N originally with H
			ser_H = " 3043"
			his_H = "6626"
			And closest water distances
			)
			atm_list = input(For additional distance enter the atom numbers for pair.)
			'''

		elif itm == '3':
			select3 = True
		elif itm == '4':
			select4 = True
		elif itm == '5':
			select5 = True
		else:
			print ("entry "+itm+" is not a valid entry")
	to_return = [select1, select2,select3,select4, select5]
	return to_return
############################################################################################


###########################################################################################
## calculate distances if you give the coordinates of two atoms
def cal_distance3(row,liginf ):
	row1 = row.values
	t1 = datetime.now()

	lig = liginf.values
	## create a name to distance
	dis_nm = lig[8]+'/'+lig[0]+'/'+lig[1]+'/'+lig[2]+'/'+lig[3]+'_'+row1[8]+'/'+row1[4]+'/'+row1[5]+\
	'/'+row1[6]+'/'+row1[7]
	#print (dis_nm)
	## assigning coordinates to x, y ,z
	x1 = float(row1[0])
	y1 = float(row1[1])
	z1 = float(row1[2])
	x2 = float(lig[4])
	y2 = float(lig[5])
	z2 = float(lig[6])
	#row[['chID','res_nm','resID','atm_type']]

	## calculate distance
	calcd_distns = round(math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2),3)

	#print("cal_distance3: pandas" , datetime.now()-t1)
	return (lig[7], row1[3],dis_nm, calcd_distns)

###########################################################################################

###########################################################################################
## calculate distances if you give the coordinates of two atoms
def cal_distance1(row,lig ):
	row1 = row.values
	#t1 = datetime.now()

	## assigning coordinates to x, y ,z
	x1 = float(row1[0])
	y1 = float(row1[1])
	z1 = float(row1[2])
	x2 = lig[0]
	y2 = lig[1]
	z2 = lig[2]
	#row[['chID','res_nm','resID','atm_type']]

	## calculate distance
	calcd_distns = round(math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2),3)

	#print("cal_distance2: pandas" , datetime.now()-t1)
	return calcd_distns
###########################################################################################

###########################################################################################
## Function to separate protein atoms from .pdb file and give out atom info in dataframeformat
## total score is also given
## need .pdb file and ligand or protein (HETATM or ATOM)
## returns a list contains protein and lignd info as dataframes and file name and score.
def get_pdb_info_into_df(pdb_file, fil_name, lignd_id = None):
	t1 = datetime.now()
	#print ("now at function get_pdb_info_in_lists")
	res_lst = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS','ILE',\
 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SEC','SER', 'THR', 'TRP', 'TYR', 'VAL']
	op_res_list = ['DCT','DCV','EDF','ETP', 'MEP', 'MVP', 'ODM', 'OMT', 'TMT','VMT', 'ACH']
	ach_lst = ['ACH']
	prot_found = False
	liga_found = False
	score_found = False
	other_lig_found = False
	prot_atm_data_list = list() ## to store all data corresponding to all atoms
	lig_atm_data_list = list()
	lig_lst = []
	sep_ligs = []
	prot_ligs_all = []
	score = None
	interface_delta_X = None
	lette2atoms = ['CL', 'FE', 'MG'] ## this is two letter atoms list. update this appropiately
	#prot_df =pd.DataFrame()
	score = 'None'
	columns = ['atmID', 'atm_type', 'res_nm', 'chID', 'resID','x_cor', 'y_cor', 'z_cor', 'atm_sym', 'pro_lig', 'seg_id']
	for line in pdb_file:
		line = line.rstrip() # read line and take off white spaces at the end
		if line == '': continue # check for empty line (to get rid of split error)
		atm_data = list() ## store single line atom data
		if line.startswith("ATOM") or line.startswith("HETATM"):
			atmID     = line[6:11].strip()  ## list index 0 atom number
			atm_typ   = line[12:16].strip() ## list index 1 atom type or name e.g. NE2
			res_nam   = line[17:21].strip() ## list index 2 residue type

			chID      = line[21:22].strip() ## list index 3 chain ID
			resID     = line[22:26].strip() ## list index 4 residue number
			x_cor     = line[30:38].strip() ## list index 5 x coordinate
			y_cor     = line[38:46].strip() ## list index 6 y coordinate
			z_cor     = line[46:54].strip() ## list index 7 z coordinate
			atm_sym   = line[77:78].strip() ## list index 8 atom symbol e.g. H, C
			prt_o_lig = line[0:6].strip()   ## list index 9 ATOM or HETATM
			seg_id    = line[71:76].strip() ## segment ID

			###### Carry out some alterations####
			## renaming histadine
			if res_nam == 'HSE':
				res_nam = 'HIS'
				
			## only for actylcholine	
			if res_nam == 'ACHO':
				## To MD and Rosetta atom names compatibility 9Acetylcholine)
				## Rosetta and MD atom names are diffrent C3=C1, O1 =OM, C6=C6, O2=O, N1=N
				ach_atoms = {'OM':'O1', 'C': 'C6', 'O': 'O2', 'N':'N1', 'C1': 'C3', 'C3': 'C1'}
				if atm_typ in ach_atoms.keys():
					#print(ach_atoms[atm_typ])
					atm_typ = ach_atoms[atm_typ]
				res_nam = 'ACH'
				
			## give atom symbols
			if len(atm_sym.strip()) == 0 or atm_sym == 'X':
				if atm_typ[0:2] in lette2atoms:
					atm_sym = atm_typ[0:2]
					#print(atm_typ[0:2])
				else:
					atm_sym = atm_typ[0]
			## Identify retro atoms
			if res_nam not in res_lst:

				prt_o_lig = 'HETATM'
				temp_chID = chID
				chID = 'X'
				temp_resID = resID
				resID = '1'
				if res_nam in ['TIP', 'TIP3','WAT']:
					res_nam = 'WAT'
					chID = 'W'
					resID = temp_resID
					
				if res_nam in ['SOD']:
					chID = temp_chID
					resID = temp_resID
					atm_sym = 'NA'
				
			if res_nam in res_lst:  ## change chain ID to A for all standard aminiacids
				chID = 'A'
			#print (atm_sym)

			atm_data.append(atmID)      ## list index 0 atom number
			atm_data.append(atm_typ)    ## list index 1 atom type or name e.g. NE2
			atm_data.append(res_nam)    ## list index 2 residue type
			atm_data.append(chID)       ## list index 3 chain ID
			atm_data.append(resID)      ## list index 4 residue number
			atm_data.append(x_cor)      ## list index 5 x coordinate
			atm_data.append(y_cor)      ## list index 6 y coordinate
			atm_data.append(z_cor)      ## list index 7 z coordinate
			atm_data.append(atm_sym )   ## list index 8 atom symbol e.g. H, C
			atm_data.append(prt_o_lig)  ## list index 9 ATOM or HETATM
			atm_data.append(seg_id)     ## list index 10 segment ID

			prot_ligs_all.append(atm_data)


		if line.startswith("total_score"):  ## change this to interface_delta_X
			score = line[12:]
			score_found = True		
		
		if line.startswith("interface_delta_X"):  ## change this to interface_delta_X
			interface_delta_X = line[18:]
			#score_found = True
		scores = (score, interface_delta_X)
	## conver list to dataframe
	prot_ligs_all_df = pd.DataFrame(prot_ligs_all, columns = columns)
	#pd.to_numeric(prot_ligs_all_df['y_cor'])
	#pd.to_numeric(prot_ligs_all_df['z_cor'])
	## select only proteins
	prot_df = prot_ligs_all_df[prot_ligs_all_df['pro_lig']=='ATOM'].copy()
	if len(prot_df) != 0:
		prot_found = True
	## get ligand name list
	lig_name_list = prot_ligs_all_df['res_nm'][prot_ligs_all_df['pro_lig']=='HETATM'].drop_duplicates().values
	if len(lig_name_list) != 0:
		liga_found = True
		if len(lig_name_list)>1:
			other_lig_found = True
	lig_dic = {}
	## seperate different ligands into seperate df in a list
	#print('Ligands found....: ', lig_name_list)
	for lig in lig_name_list:
		lig_df= prot_ligs_all_df[prot_ligs_all_df['res_nm']==lig].copy()
		## stote OP ligand info
		if lig in op_res_list:
			lig_dic['OP'] = lig_df
		elif lig == 'WAT':
			lig_dic['WAT'] = lig_df
		else:
			lig_dic[lig] = lig_df
	## if no water found assgn dummy water data
	if 'WAT' not in lig_name_list:	
		fake_wat = [['0', 'O', 'WAT', 'W', '0', '0.00', '0.00', '0.00', 'O', 'HETATM', 'WT1'],\
		['0', 'H', 'WAT', 'W', '0', '0.00', '0.00', '0.00', 'H', 'HETATM', 'WT1']]	
		fake_wat_df = pd.DataFrame(fake_wat, columns = columns)
		lig_dic['WAT'] = fake_wat_df


	if prot_found == False:
		print ("No protein ATOM list found!")

	if other_lig_found == False and liga_found == False:
		print ("No ligand found!")
	#if score_found == False:
	#	print (" No score record found!")

	## lig_dic can have several ligand info dataframes
	#pdb_info = [prot_atm_data_list,lig_dic, score, fil_name]
	pdb_info = [prot_df,lig_dic, scores, fil_name,prot_ligs_all_df]
	#print("get_pdb_info_in_lists--pandas----: ", datetime.now()-t1)
	return pdb_info


###########################################################################################

def get_selec_distances2(to_df, from_df):

	#totstarttime = datetime.now()
	## to_df is the df to get distance from atoms from from_df
	colums = ['ligAtmID', 'OtherAtmID','dist_name', 'dist']
	lig_dis_df = pd.DataFrame(columns =colums)
	all_dist_df = pd.DataFrame(columns =colums)
	df_tuples = pd.DataFrame()
	df=from_df[['atmID','x_cor','y_cor','z_cor', 'atm_type','res_nm','resID','chID', 'seg_id']].copy()
	#print(df.ix[:,['atmID','x_cor','y_cor','z_cor', 'atm_type','res_nm','resID','chID']])
	to_df1 = to_df[['x_cor','y_cor', 'z_cor','atmID','chID', 'res_nm', 'resID', 'atm_type', 'seg_id']].copy()
	#print(to_df1)
	if len(df) != 0 and len(to_df1) != 0:
		for ligline in df.index:
			lig_atm_nm = df.ix[ligline, ['chID','res_nm','resID','atm_type','x_cor','y_cor','z_cor','atmID', 'seg_id']]
			#print(lig_atm_nm)
			df_tuples['data'] = to_df1.apply(lambda row :cal_distance3(row,lig_atm_nm), axis =1)#\
	
			new_col_list = ['ligAtmID', 'OtherAtmID', 'dist_name','dist']
			for n,col in enumerate(new_col_list):
				lig_dis_df[col] = df_tuples['data'].apply(lambda data: data[n])
			if len(lig_dis_df) == 0: continue
	
	
			lig_dis_df = lig_dis_df[lig_dis_df['dist']<20]  ## can change distance cutoff here
			all_dist_df = all_dist_df.append(lig_dis_df, ignore_index=True)
		#print(all_dist_df)

	#print("get_selec_distances pandas2: ", datetime.now()-totstarttime)
	return all_dist_df
###########################################################################################

def get_selec_distances_x(t_df, q):
	to_df, from_df = t_df[0], t_df[1]
	
	#totstarttime = datetime.now()
	## to_df is the df to get distance from atoms from from_df
	colums = ['ligAtmID', 'OtherAtmID','dist_name', 'dist']
	lig_dis_df = pd.DataFrame(columns =colums)
	all_dist_df = pd.DataFrame(columns =colums)
	df_tuples = pd.DataFrame()
	df=from_df[['atmID','x_cor','y_cor','z_cor', 'atm_type','res_nm','resID','chID']].copy()
	#print(df.ix[:,['atmID','x_cor','y_cor','z_cor', 'atm_type','res_nm','resID','chID']])
	to_df1 = to_df[['x_cor','y_cor', 'z_cor','atmID','chID', 'res_nm', 'resID', 'atm_type']].copy()
	#print(to_df1)
	for ligline in df.index:
		lig_atm_nm = df.ix[ligline, ['chID','res_nm','resID','atm_type','x_cor','y_cor','z_cor','atmID']]
		#print(lig_atm_nm)
		df_tuples['data'] = to_df1.apply(lambda row :cal_distance3(row,lig_atm_nm), axis =1)#\

		new_col_list = ['ligAtmID', 'OtherAtmID', 'dist_name','dist']
		for n,col in enumerate(new_col_list):
			lig_dis_df[col] = df_tuples['data'].apply(lambda data: data[n])
		if len(lig_dis_df) == 0: continue


		lig_dis_df = lig_dis_df[lig_dis_df['dist']<20]  ## can change distance cutoff here
		all_dist_df = all_dist_df.append(lig_dis_df, ignore_index=True)
	#print(all_dist_df)
	q.put(all_dist_df)
	#print("get_selec_distances pandasx: ", datetime.now()-totstarttime)
	return all_dist_df






###########################################################################################

def write_to_csv(df, file_name,file_n,is_new,score='None'):


	if is_new ==True:
		with open(file_name, 'w') as f:
			f.write(file_n+','+score+'\n')

			df.to_csv(f, header=True,index=False)

	elif is_new == False:
		with open(file_name, 'a') as f:
			f.write(file_n+','+score+'\n')

			df.to_csv(f, header=False,index=False)

###########################################################################################

## Thsi function get approximate centre of mass of a molecule based on the coordinates
def get_aprox_centerof_molecule(df):

	## for a molecule. (multiple atoms)
	if len(df) > 1:
		x_min = df['x_cor'].min()
		x_max = df['x_cor'].max()
		y_min = df['y_cor'].min()
		y_max = df['y_cor'].max()
		z_min = df['z_cor'].min()
		z_max = df['z_cor'].max()
		x= (float(x_min) + float(x_max))/2
		y= (float(y_min) + float(y_max))/2
		z= (float(z_min) + float(z_max))/2
	## for single atom
	elif len(df)==1:
		x = df['x_cor']
		y = df['y_cor']
		z = df['z_cor']
	else:
		print('No atoms found for calculating center of mass')

	return (x,y,z)
###########################################################################################
def get_atoms_arount_actvsite(df,centre,threshold = 10.0):
	to_df2 = df[['x_cor','y_cor', 'z_cor','atmID','chID', 'res_nm', 'resID', 'atm_type', 'pro_lig','atm_sym']].copy()
	## calculate the dictance to an atom from a given centre
	to_df2['dist_from_centre'] = to_df2.apply(lambda row :cal_distance1(row,centre), axis =1)

	df_new = to_df2[to_df2['dist_from_centre']<threshold].copy()
	df_new.drop('dist_from_centre', axis=1, inplace=True)

	return df_new

###########################################################################################
def get_connectivity(df):
	pass

###########################################################################################

def reorder_colmns(df):

	#df6 =df.copy()
	#w1 = datetime.now()
	## get a list of columns
	cols = list(df)
	## move the column to head of list using index, pop and insert
	cols.insert(0, cols.pop(cols.index('Molecule_name')))
	cols.insert(1, cols.pop(cols.index('Score')))
	## use ix to reorder
	df = df.ix[:, cols]

	#print("col reorder 1: ", datetime.now()-w1)

	'''
	## this method also working
	w2 = datetime.now()
	mol_nm = df6['Molecule_name']
	scr = df6['Score']
	df6.drop(labels=['Molecule_name'], axis=1,inplace = True)
	df6.insert(0, 'Molecule_name', mol_nm)
	df6.drop(labels=['Score'], axis=1,inplace = True)
	df6.insert(1, 'Score', scr)

	print("col reorder 2: ", datetime.now()-w2)
	'''
	return df
###########################################################################################
def cal_rms(df_row):
	## to calculate RMS, only numeric values of DF should be passed to this function
	#print (df_row)
	#print('DONE')
	sqr = 0
	for val in df_row:

		#print ('val: ', val)
		sqr = sqr + val**2
		#print ('sqr: ', sqr)
	row_rms = round(math.sqrt(sqr/(len(df_row))),3)

	#print('rms: ', row_rms)
	#print(df_row[0],df_row[1] ,df_row[2],row_rms)

	return row_rms
###########################################################################################

def impnt_dist_vs_score(prt_inf_,lig1_, lig2_):
	#def impnt_dist_vs_score(score_, file_,prt_inf_,lig1_, lig2_,select_dist_df_, file_count_, time_st_ ):
	#t1 = datetime.now()
		
	'''
	his437 HE2 6626, his NE2 6619 his ND1 6616
	ser 199,OG3038, HG 3043
	glu 198 OE1 3025, OE2 3026
	'''

	resi_list = ['437','199' ,'198']
	## FOR SOME OP
	#atm_lst = ['OG', 'HG','HG1', 'OE1', 'OE2', 'ND1','ND2','NE1', 'NE2', 'HD1','HE2']
	## FOR ACH
	atm_lst = ['OG', 'HG','HG1', 'OE1', 'OE2', 'ND1','ND2','NE1', 'NE2', 'HD1','HE2'] ##***************** modify loc1
	hevy_atoms = ['P','S','O','N','Cl','CL','C','F']
	### To select appropriate His atoms
	only_impotant_atms = prt_inf_[prt_inf_['resID'].isin(resi_list) & prt_inf_['atm_type'].isin(atm_lst)]
	#prt_inf_ = prt_inf_2[prt_inf_2['resID'].isin(resi_list)]
	#only_impotant_atms = prt_inf_[prt_inf_['atm_type'].isin(atm_lst)]
	#imp_atm_list = ['6626','6619','6616', '3038','3043','3025', '3026']
	## Selected protein atoms
	#only_impotant_atms = prt_inf_.apply(lambda row :select_atoms(row), axis =1).dropna()
	
	
	
	## change this part accordingly	
	#if lig1_.ix[lig1_.index[1],'res_nm'] in ['DCV', 'DCT', 'EDF', 'ETP', 'MEP', 'MVP', 'ODM', 'OMT', 'VMT']:
	#	## selected ligand atoms for some OP
	#	lig_inf = lig1_[lig1_['atm_sym'].isin(['P','S','O','N','Cl', 'CL'])]
	#	
	#if lig1_.ix[lig1_.index[1],'res_nm'] in ['TMT']:
	#	## selected ligand atoms for some OP
	#	lig_inf = lig1_[lig1_['atm_sym'].isin(['P','S','O','N','Cl', 'CL'])]
	
	prot_all = prt_inf_[prt_inf_['atm_sym'].isin(hevy_atoms)]
	
	if lig1_.ix[lig1_.index[1],'res_nm'] in ['ACH']:
		## selected ligand atoms for ACH
		lig_inf = lig1_[lig1_['atm_type'].isin(['O1','O2','C6','N1','C1', 'C2','C4'])] 
		lig_inf_all = lig1_[lig1_['atm_sym'].isin(hevy_atoms)]
	else:
		
		lig_inf = lig1_[lig1_['atm_sym'].isin(['P','S','O','N','Cl', 'CL'])]
		## lig_inf2 only for water distances calculation, just to lower number of distances
		lig_inf2 = lig1_[lig1_['atm_type'].isin(['P1','S1','O1'])]
		lig_inf_all = lig1_[lig1_['atm_sym'].isin(hevy_atoms)]
	

	## get distances form serine HG to histidine ND ans NE
	## get important his atoms
	#his_atms = only_impotant_atms[only_impotant_atms['atmID'].isin(['6616','6619','6626'])]
	his_atms = only_impotant_atms[only_impotant_atms['resID'] == '437']
	## get important ser atoms
	ser_H_O =  only_impotant_atms[only_impotant_atms['resID']=='199']
	## get important GLU atoms
	#glu_O =  only_impotant_atms[only_impotant_atms['atmID'].isin(['3026','3025'])]
	glu_O =  only_impotant_atms[only_impotant_atms['resID'] == '198']
	## combine his and ser atom info
	his_ser_OP_atms = lig_inf.append([his_atms,ser_H_O], ignore_index=True)
	## get distances to waters from ligand1 selected atoms and protein selected atioms
	lig2_only_water_O = lig2_[lig2_['atm_sym'].isin(['O','OH2' ])]
	#lig_inf2 = lig1[lig1['atm_type'].isin(['P1','O1'])]
	

	
	
	'''
	his_ser_dis_df = get_selec_distances2(ser_H_O, his_atms[his_atms['atm_type'].isin(['NE2','ND1', 'NE1', 'ND2'])])	
		
	## get only selected important distances with ligand and protein (important residues)
	select_impt_dis_df = get_selec_distances2(only_impotant_atms, lig_inf)
	## calculate the distances between impotant glu atoms to selected his atoms
	his_glu_dis_df = get_selec_distances2(glu_O, his_atms)

	## calculte distances to water molecules from his and ser selected atoms
	lig_water_O_dis_df = get_selec_distances2(lig2_only_water_O, his_ser_OP_atms)
	'''
	lig_pro_hevy_atm = get_selec_distances2(prot_all, lig_inf_all)
	
	
	lig_water_dist = get_selec_distances2(lig2_only_water_O, lig_inf) ## note lig_inf2 vs lig_inf
	his_water_dist = get_selec_distances2(lig2_only_water_O, his_atms[his_atms['atm_type'].isin(['NE2','ND1', 'NE1', 'ND2'])])
	ser_water_dist = get_selec_distances2(lig2_only_water_O, ser_H_O)
	
	## calculate the distances between impotant glu atoms to selected his atoms
	his_glu_dis_df = get_selec_distances2(glu_O, his_atms)
	#print('####################################',len(his_atms[his_atms['atm_type'].isin(['NE2','ND1', 'NE1', 'ND2'])]), len(ser_H_O))
	his_ser_dis_df = get_selec_distances2(ser_H_O, his_atms[his_atms['atm_type'].isin(['NE2','ND1', 'NE1', 'ND2'])])
		
	## get only selected important distances with ligand and protein (important residues)
	select_impt_dis_df = get_selec_distances2(only_impotant_atms, lig_inf)
	
	return [lig_water_dist,his_water_dist, ser_water_dist, his_glu_dis_df, his_ser_dis_df, select_impt_dis_df, lig_pro_hevy_atm]
###########################################################################################	
###########################################################################################	
def df_append_transpose(main_df, df_list_to_append, file__, score__, time_st__, file_count__):
	## main_df is ths df to append other dfs
	## df_list_to_append a list og data frmes
	## append both prot-lig and prot-prot distances
	filenm_scor_df = pd.DataFrame({'Molecule_name':file__, 'Score': score__[0], 'interface_delta_X': score__[1], 'Time_step':time_st__ }, index = [0])
	select_impt_dis_df_ = main_df 
	select_impt_dis_df_ = (select_impt_dis_df_.append(df_list_to_append, ignore_index=True))[['dist_name', 'dist']].copy()
	select_impt_dis_df_.set_index('dist_name', inplace=True)

	select_impt_dis_df_ = select_impt_dis_df_.transpose()
	select_impt_dis_df_.set_index([[0]], inplace=True)

	## concatenate molecule name and score in to a  dataframe which contains selected distance informatio
	tempdf = pd.concat([filenm_scor_df,select_impt_dis_df_], join_axes  = [select_impt_dis_df_.index], axis=1)
	#print(file_count__)
	
	#print("impnt_dist_vs_score pandas2: ", datetime.now()-t1)
	return tempdf
	
	#return df_concatenate(select_dist_df_, tempdf)
	
###########################################################################################
## this dataframe receive cumulative didtance data frame and,
## one line transpose data frame to concatenate with cumulative df
def df_concatenate(select_dist_df__, tempdf__):
	
	select_dist_df__.drop_duplicates(inplace=True)
	#tempdf__.drop_duplicates(inplace=True)
	if  True: #len(tempdf__) != 0:#	len(select_dist_df__) !=0 or
		## merge temporary dataframe to accumulation df
		select_dist_df__ = pd.concat([select_dist_df__, tempdf__ ],axis=0,ignore_index =True, join='outer')
		#select_dist_df = select_dist_df.combine(tempdf, func)
		##select_dist_df__ = pd.merge(select_dist_df__, tempdf__ , how='outer', sort=False)
		#print('Not zero   iiiiiiiiiiiiiiiiiiii')
		#print('Column lengtth: ', len(select_dist_df__.columns))
	
		## re arrange column sequence
		select_dist_df__ = reorder_colmns(select_dist_df__)
	
	return select_dist_df__
	#print(select_dist_df.columns)
	#select_dist_df_.to_csv(file_[10:30]+'impotant_dist.csv', index=False)
	
###########################################################################################
def plot_score_vs_RMS(dat_frm):

	#print(prf_dst_df['RMS'])
	#pref_dist1 = pref_dist[(pref_dist['X/'+resi_name+'/1/P1_A/SER/199/OG'] < 5.0) & (pref_dist['A/HIS/437/ND1_A/SER/199/HG'] < 4.0 ) & (pref_dist['X/'+resi_name+'/1/O1_A/HIS/437/HD1'] < 4.0)]
	#print(pref_dist4)#['Molecule_name'].values)
	pref_dis_rms = dat_frm['RMS']
	
	if len(pref_dis_rms) > 100:
		
		try:
			scor = dat_frm['interface_delta_X']
		except:
			scor = dat_frm['Score']
		#print(dat_frm)
		
		
		
		
		fig_rms =plt.figure(figsize=(12,8))
		ax3 = fig_rms.add_subplot(111)
		ax3.scatter(pref_dis_rms,scor,color = 'g', s = 20, marker ='o' )
		ax3.set_xlim(2.5, 7.5)
		#ax3.set_ylim(-450, -250)
		#ax.axis(xmin = 2.75, xmax = 5.25) # this works too
		###ax3.set_title('RMS vs Score (for MD, score = nonbonded energy (kcal/mol); for Rosetta docking, score = Rosetta score)') ## just added some letex
		ax3.set_xlabel(('RMS: '+str(get_rms_cols(dat_frm))))
		ax3.set_ylabel('E(Non bonded)')
		ax3.grid()
		plt.show()
	else:
		print ('Less than 100 occurances')



	dat_frm = pd.DataFrame()
###########################################################################################
def plot_score_vs_interface_delta(dat_frm):
	#print(prf_dst_df['RMS'])
	#pref_dist1 = pref_dist[(pref_dist['X/'+resi_name+'/1/P1_A/SER/199/OG'] < 5.0) & (pref_dist['A/HIS/437/ND1_A/SER/199/HG'] < 4.0 ) & (pref_dist['X/'+resi_name+'/1/O1_A/HIS/437/HD1'] < 4.0)]
	#print(pref_dist4)#['Molecule_name'].values)
	#pref_dis_rms = dat_frm['RMS']
	
	if len(dat_frm.dropna()) > 100:
		#scor = dat_frm['Score']
		inf_delta = dat_frm['interface_delta_X']
		score = dat_frm['Score']
		
		
		fig_rms =plt.figure(figsize=(12,8))
		ax3 = fig_rms.add_subplot(111)
		ax3.scatter(inf_delta,score,color = 'g', s = 20, marker ='o' )
		#ax3.set_xlim(3.0, 7.5)
		#ax3.set_ylim(-450, -250)
		#ax.axis(xmin = 2.75, xmax = 5.25) # this works too
		###ax3.set_title('RMS vs Score (for MD, score = nonbonded energy (kcal/mol); for Rosetta docking, score = Rosetta score)') ## just added some letex
		ax3.set_xlabel('E(interface_delta_X)')
		ax3.set_ylabel('E(Score)')
		ax3.grid()
		plt.show()
	else:
		print ('Less than 100 occurances')



	dat_frm = pd.DataFrame()
	
###########################################################################################	
def plot_dist_vs_TS(dat_f, dis_nm):
	
	dat_frm = dat_f.copy()
	dat_frm.sort_values(by='Time_step', ascending=True, inplace=True )

	if len(dat_frm['Time_step'].dropna()) > 10:
		#print(prf_dst_df['RMS'])
		#pref_dist1 = pref_dist[(pref_dist['X/'+resi_name+'/1/P1_A/SER/199/OG'] < 5.0) & (pref_dist['A/HIS/437/ND1_A/SER/199/HG'] < 4.0 ) & (pref_dist['X/'+resi_name+'/1/O1_A/HIS/437/HD1'] < 4.0)]
		#print(pref_dist4)#['Molecule_name'].values)
	
		scor = dat_frm['Score']
		rms = dat_frm['RMS']
		time_step =  dat_frm['Time_step']
		
		#print (dat_frm.head())
		
		
		y_item_list = [scor, rms]
		
		for y_item in y_item_list:
			#print(y_item)
		
		
			fig_rms =plt.figure(figsize=(12,8))
			ax3 = fig_rms.add_subplot(111)
			#ax3.scatter(time_step,y_item,color = 'g', s = 20, marker ='o' )
			ax3.plot(time_step,y_item,'.g-' )
			ax3.set_xlim(0.0, 10000.0)
			#ax3.set_ylim(-500, -250)
			#ax.axis(xmin = 2.75, xmax = 5.25) # this works too
			####ax3.set_title('Time_step vs Distance (in A)') ## just add some letex
			#axis_tit = dat_frm.columns.values
			ax3.set_xlabel('Time_step (x 100 fs)')
			#ax3.set_xlabel(('RMS: '+str(axis_tit[2:(len(axis_tit)-1)])))
			
			if y_item.name == 'RMS':
				name = y_item.name+'\n'+str(dis_nm)
			else:
				name = y_item.name
			ax3.set_ylabel(name)
			ax3.grid()
			plt.show()
	else:
		print ('Less than 100 occurances or no occurances available')
	
		
	dat_frm = pd.DataFrame()
	
	
	
	
	
###########################################################################################	
## This function plots time step dat like ts vs RMS and TS vs Score
def plot_dist_vs_TS2(dat_f_lst):
	import matplotlib.mlab as mlab
	num_of_dfs = len(dat_f_lst) ## get number of dataframes 
	
	# start with a rectangular Figure
	#fig = plt.figure(1, figsize=(30,30),dpi=100)
	fig = plt.figure(1, dpi=300, figsize = (num_of_dfs*2+4,num_of_dfs*2+4), facecolor='#f0f0f0')
	
	
	
	ax = [[] for i in range(num_of_dfs)]
	ax_his = [[] for i in range(num_of_dfs)]
	hist_line = [[] for i in range(num_of_dfs)]
	row_len = num_of_dfs*2
	
	gs = gridspec.GridSpec(row_len, 7)
	
	roBegin = 0
	row_end = 0
	for i, df in enumerate(dat_f_lst):
		dat_frm = df[0]
		
		col_txt = ''
		for itm in df[1]:
			col_txt = col_txt+itm+' : '
		col_txt = col_txt.rstrip(' : ')

		
		

		dat_frm.sort_values(by='Time_step', ascending=True, inplace=True )
		
		scor = dat_frm['Score']  ## y 
		scor_min, scor_max = scor.min(), scor.max()
		rms = dat_frm['RMS']  ## y 
		rms_min, rms_max, rms_mean, rms_std = rms.min(), rms.max(), rms.mean(), rms.std()
		time_step =  dat_frm['Time_step']  ## x 
		ts_min, ts_max = time_step.min(), time_step.max()
		

		
		if i == 0: ## set reference plot
			roBegin = i
			row_end = roBegin +2
			ax[i] = plt.subplot(gs[roBegin:row_end, 0:5])
			ax[i].plot(time_step, rms)		
			ax[i].set_xlim(ts_min, ts_max)
			ax[i].get_xaxis().set_visible(False)
			ax[i].set_ylabel('RMS', fontsize=16, fontproperties='Times New Roman')
			#ax[i].legend(loc=9, ncol=1, prop={'size': 11})
			#ax[i].set_title(col_nam)
			ax[i].text(.5,.9, 'RMS = '+col_txt, horizontalalignment='center',transform=ax[i].transAxes,\
				fontsize=10, fontproperties='Times New Roman')
			
			ax_his[i] = plt.subplot(gs[roBegin:row_end, 5:],sharey=ax[0])
			
			n, bins, patches = ax_his[i].hist(rms, bins=10,orientation='horizontal', color='g',normed=1, facecolor='green', alpha=0.5)
			y = mlab.normpdf(bins, rms_mean, rms_std)			
			hist_line[i] = plt.plot(y, bins, 'r-')
			
			ax_his[i].get_yaxis().set_visible(False)
			ax_his[i].get_xaxis().set_visible(False)
			
		else:
			roBegin = row_end
			row_end = roBegin +2
			ax[i] = plt.subplot( gs[roBegin:row_end, 0:5], sharex=ax[0])
			ax[i].plot(time_step, rms)		
			ax[i].set_xlim(ts_min, ts_max)
			plt.subplots_adjust(wspace=0, hspace=0.2)
			ax[i].set_ylabel('RMS', fontsize=16, fontproperties='Times New Roman')
			ax[i].text(.5,.9, 'RMS = '+col_txt, horizontalalignment='center',transform=ax[i].transAxes,\
				fontsize=10, fontproperties='Times New Roman')
			
			ax_his[i] = plt.subplot(gs[roBegin:row_end, 5:],sharey=ax[0], sharex=ax_his[0])
			n, bins, patches = ax_his[i].hist(rms, bins=10,orientation='horizontal', color='g',normed=1, facecolor='green', alpha=0.5)
			y = mlab.normpdf(bins, rms_mean, rms_std )			
			hist_line[i] = plt.plot(y, bins, 'r-')
			
			ax_his[i].get_yaxis().set_visible(False)
			
			
			if i != num_of_dfs -1:
				ax[i].get_xaxis().set_visible(False)
				ax_his[i].get_xaxis().set_visible(False)
				
			elif i == num_of_dfs -1:
				
				ax[i].set_xlabel('Time step', fontsize=16, fontproperties='Times New Roman')
				ax_his[i].set_xlabel('Count%', fontsize=16, fontproperties='Times New Roman')
			#ax1_hist.set_yticklabels([])
			#fig.add_subplot(ax[i])
		

	plt.show()
	fig.savefig('TS_vs_RMS.png', facecolor=fig.get_facecolor())
	
		
	dat_frm = pd.DataFrame()
###########################################################################################

	
def read_huge_file(big_fname):

	with open(big_fname, 'rb') as big_file:##closes file after all the lines have been processed
		for line in big_file: ## not using readlines(), because it store lines in memory
			data = line.decode('ascii')
			#lene2 =struct.unpack('c',line)
			#line2 =decode(line,'utf-8')
			print (data)

###########################################################################################
def menu_0(prt_inf0, lig10, lig20, score0, lig_centr_coord0 ,is_new_prot_run0, is_new_water_run0, f_nam0):
	## get all the atom info from a given given treshold around a given centre
	prot_atoms_arond_act_site = get_atoms_arount_actvsite(prt_inf0,lig_centr_coord0 ,14.0)

	## get distances from ligand1 atoms to ligand2 atoms
	#if lig20 != None:
	lig_water_dis_df = get_selec_distances2(lig20, lig10)
	## write water ligand distances to csv
	write_to_csv(lig_water_dis_df, f_nam0[10:30]+'_lig_water.csv',f_nam0,is_new_water_run0,score0)
	## get distances from ligand1 atom to protein atoms
	lig_prot_dis_df = get_selec_distances2(prot_atoms_arond_act_site, lig10)
	#lig_prot_dis_df = get_selec_distances2(prt_inf0, lig10)

	## write peotein ligand distances to csv

	write_to_csv(lig_prot_dis_df, f_nam0[10:30]+'_lig_prot.csv',f_nam0,is_new_prot_run0,score0)
	#is_new_prot_run0 = False

	#is_new_water_run0 = False

	#return [is_new_prot_run, is_new_water_run]

###########################################################################################
## get selected column values from a data frame
def get_rms_cols(rms_df):
	dis_cols = re.compile(r'\w/\w{3}/\d{1,3}/.{1,3}_\w/\w{3}/\d{3}/.{1,3}')
	rms_cols = []
	for col in rms_df.columns:
		if len(re.findall(dis_cols, col))>0:
			rms_cols.append(col)
	return rms_cols
	## list comprehension way
	#rms_cols = [re.findall(dis_cols, col) for col in prf_dst_df.columns if len(re.findall(dis_cols, col))>0]



###########################################################################################
def menu_2():

	csv_dist_file_name = input('''Enter file name and path to distanes and score CSV file
	Default: 'dynamics_best_comple_impotant_dist.csv
	''')

	if csv_dist_file_name == '':
		csv_dist_file_name = 'dynamics_best_comple_impotant_dist.csv'
	impt_dis_df = pd.read_csv(csv_dist_file_name,low_memory=False)
	colums = impt_dis_df.columns.values
	#print (colums)

	## get residue name
	resi_select = re.compile('X/\w{3}/.*')
	resi_name = [m.group(0) for l in colums for m in [resi_select.search(l)] if m][0][2:5]

	#impt_dis_df['hisHs'] = pd.concat([impt_dis_df['X/'+resi_name+'/1/O1_A/HIS/437/HE2'].dropna(), impt_dis_df['X/'+resi_name+'/1/O1_A/HIS/437/HD1'].dropna()]).reindex_like(impt_dis_df)

	#print(impt_dis_df['hisHs'])
	#####-----------start HISTIDINE-------------------------------------------------------######
	##************************************************************************************************************	
	##without H involved
	
	if resi_name == 'ACH':
		dis_list1 = ['X/'+resi_name+'/1/C6_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/OG','A/HIS/437/NE2_A/GLU/198/OE2']
		dis_list2 = ['X/'+resi_name+'/1/C6_A/SER/199/OG','A/HIS/437/NE2_A/SER/199/OG','A/HIS/437/NE2_A/GLU/198/OE2']
		dis_list3 = ['X/'+resi_name+'/1/C6_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/OG','A/HIS/437/ND1_A/GLU/198/OE2']
		dis_list4 = ['X/'+resi_name+'/1/C6_A/SER/199/OG','A/HIS/437/NE2_A/SER/199/OG','A/HIS/437/ND1_A/GLU/198/OE2']
		dis_list5 = ['X/'+resi_name+'/1/C6_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/OG','A/HIS/437/NE2_A/GLU/198/OE1']
		dis_list6 = ['X/'+resi_name+'/1/C6_A/SER/199/OG','A/HIS/437/NE2_A/SER/199/OG','A/HIS/437/NE2_A/GLU/198/OE1']
		dis_list7 = ['X/'+resi_name+'/1/C6_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/OG','A/HIS/437/ND1_A/GLU/198/OE1']
		dis_list8 = ['X/'+resi_name+'/1/C6_A/SER/199/OG','A/HIS/437/NE2_A/SER/199/OG','A/HIS/437/ND1_A/GLU/198/OE1']

	elif resi_name == 'TMT':
		dis_list1 = ['X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/OG','X/'+resi_name+'/1/S1_A/HIS/437/ND1']
		dis_list2 = ['X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/NE2_A/SER/199/OG','X/'+resi_name+'/1/S1_A/HIS/437/NE2']
		dis_list3 = ['X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/OG','X/'+resi_name+'/1/S1_A/HIS/437/NE2']
		dis_list4 = ['X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/NE2_A/SER/199/OG','X/'+resi_name+'/1/S1_A/HIS/437/ND1']
	else:
		dis_list1 = ['X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/OG','X/'+resi_name+'/1/O1_A/HIS/437/ND1']
		dis_list2 = ['X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/NE2_A/SER/199/OG','X/'+resi_name+'/1/O1_A/HIS/437/NE2']
		dis_list3 = ['X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/OG','X/'+resi_name+'/1/O1_A/HIS/437/NE2']
		dis_list4 = ['X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/NE2_A/SER/199/OG','X/'+resi_name+'/1/O1_A/HIS/437/ND1']
	##************************************************************************************************************
	non_bonded = 'interface_delta_X'
	if len(impt_dis_df['interface_delta_X'].dropna())==0:		
		non_bonded = 'Score'
	
	if resi_name == 'ACH':
		pref_dist1 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list1].dropna().copy()
		pref_dist2 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list2].dropna().copy()
		pref_dist3 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list3].dropna().copy()
		pref_dist4 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list4].dropna().copy()
		pref_dist5 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list5].dropna().copy()
		pref_dist6 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list6].dropna().copy()
		pref_dist7 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list7].dropna().copy()
		pref_dist8 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list8].dropna().copy()
	else:
		pref_dist1 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list1].dropna().copy()
		pref_dist2 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list2].dropna().copy()
		pref_dist3 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list3].dropna().copy()
		pref_dist4 = impt_dis_df[['Molecule_name', non_bonded, 'Time_step'] + dis_list4].dropna().copy()
		
	
	try:
		score_vs_interface_delta_X = impt_dis_df[['interface_delta_X', 'Score']].dropna().copy()
		plot_score_vs_interface_delta(score_vs_interface_delta_X)
	except:
		print('Not a Rossetta calculation')
	


	
	#############################################################################################
	#############################################################################################
	if resi_name == 'ACH':
		prf_dst_df_lst = [pref_dist1, pref_dist2, pref_dist3, pref_dist4, pref_dist5, pref_dist6,\
		 pref_dist7, pref_dist8]
	else:
		prf_dst_df_lst = [pref_dist1, pref_dist2, pref_dist3, pref_dist4]

	is_new_his_df =True

	print('######################################################################################')
	print('Histadine involment')
	
	plot_ts_lst = []
	for prf_dst_df in prf_dst_df_lst:
		if len(prf_dst_df) == 0: 
			continue
		

		rms_col = get_rms_cols(prf_dst_df)
		prf_dst_df['RMS'] = prf_dst_df.ix[:,rms_col].apply(lambda row :cal_rms(row), axis =1)

		
		#prf_dst_df_unsorted = prf_dst_df.sort_values(by=prf_dst_df.index, ascending=True).copy()

		

		prf_dst_df.sort_values(by=non_bonded, ascending=True, inplace=True )
		#####un comment_following
		plot_score_vs_RMS(prf_dst_df)
		
		if prf_dst_df['Time_step'].sum() != 0:
			
			#plot_dist_vs_TS(prf_dst_df, rms_col)
			
			plot_ts_lst.append((prf_dst_df, rms_col))
		
		

		## write dfs into one file

		df_colums = prf_dst_df.columns.values

		data_title = str(df_colums[2:(len(df_colums)-1)])
		#data_title2 = str(df_colums[0:(len(df_colums)-1)])
		#print (data_title2)
		write_to_csv(prf_dst_df.head(20), 'top_20_rms_vs_score.csv',data_title,is_new_his_df)
		#write_to_csv(prf_dst_df_unsorted, 'top_20_rms_vs_score_unsort.csv',data_title2,is_new_his_df)

		is_new_his_df = False


	plot_dist_vs_TS2(plot_ts_lst)
	
	
	
	
	
	
	#####----------END HISTIDINE---------------------------------------------------------######
"""
	#####----------start WATER------------------------------------------------------######

	## search and select a components in a given list
	## select phophinyl O ans water O
	po_O_wat_regex = re.compile('X/'+resi_name+'/1/O1_W/WAT/.*')
	phosphynil_O_wat = [m.group(0) for l in colums for m in [po_O_wat_regex.search(l)] if m]  ## m.group(0) gives all part of the selection

	## select serinr O ans water O
	ser_O_wat_regex = re.compile('A/SER/199/[O]G_W/WAT/.*')
	ser_O_wat = [m.group(0) for l in colums for m in [ser_O_wat_regex.search(l)] if m]

	## select serinr H ans water O
	ser_H_wat_regex = re.compile('A/SER/199/[H]G\d{0,1}_W/WAT/.*')
	ser_H_wat = [m.group(0) for l in colums for m in [ser_H_wat_regex.search(l)] if m]


	for pO_wat in phosphynil_O_wat:
		wat_pO = pO_wat.split('_')[1]
		
		print('######################################################################################')
		print('Water HETEROATOMS involment')

		## arange data to plot serine O and water O
		for serO_wat in ser_O_wat:
			wat_serO = serO_wat.split('_')[1]
			


			if wat_pO == wat_serO:
				##******************************************** P1 to C6
				wat_assist_df = impt_dis_df[['Molecule_name', 'Time_step', 'Score','X/'+resi_name+'/1/C6_A/SER/199/OG', pO_wat, serO_wat]].dropna()

				if len(wat_assist_df) == 0: continue ## do not pass empty dataframe to plot
				wat_assist_df['RMS'] = wat_assist_df.ix[:, 3:].apply(lambda row :cal_rms(row), axis =1)
				plot_dist_vs_TS(wat_assist_df)
				wat_assist_df.sort_values(by=['RMS','Score'], ascending=[True,True], inplace=True )
				plot_score_vs_RMS(wat_assist_df)
		#############################################################################################
		#############################################################################################
		## To add water H uncomment following block
		'''
		print('######################################################################################')
		print('Water HYDROGEN involment')
		## arange data to plot serine H and water O
		for serH_wat in ser_H_wat:
			wat_serH = serH_wat.split('_')[1]

			if wat_pO == wat_serH:


				wat_assist_df2 = impt_dis_df[['Molecule_name', 'Time_step', 'Score','X/'+resi_name+'/1/P1_A/SER/199/OG', pO_wat, serH_wat]].dropna()

				if len(wat_assist_df2) == 0: continue ## do not pass empty dataframe to plot
				wat_assist_df2['RMS'] = wat_assist_df2.ix[:, 3:].apply(lambda row :cal_rms(row), axis =1)
				plot_dist_vs_TS(wat_assist_df2)
				wat_assist_df2.sort_values(by=['RMS','Score'], ascending=[True,True], inplace=True )
				plot_score_vs_RMS(wat_assist_df2)
				
		'''
		#############################################################################################
		#############################################################################################
	#####----------end WATER------------------------------------------------------######

	############## examples for regular expression #########################
	#water_df = impt_dis_df.filter(regex=('X/'+resi_name+'/1/O1_W/WAT/.*')) ## select useing reg ex
	#regexp = re.compile(r'\b[A-Z]{3,}\b')
	# use only one of the following lines, whichever you prefer           A/SER/199/[H,O]G_W/WAT/.*']+
	#filtered = filter(lambda i: not regex.search(i), full)
	#filtered = [i for i in full if not regex.search(i)]

	#rpt[rpt['STK_ID'].str.contains(r'^600[0-9]{3}$')] # ^ means start of string # [0-9]{3} means any three digits
	# $ means end of string
	#endstrings = ['01$', '02$', '05$']
	#rpt[rpt['STK_ID'].str.contains('|'.join(endstrings)]
	########################################################################
	'''
	pref_dist = impt_dis_df[['Molecule_name', 'Score','X/'+resi_name+'/1/P1_A/SER/199/OG','A/HIS/437/ND1_A/SER/199/HG','X/'+resi_name+'/1/O1_A/HIS/437/HD1']].dropna()
	pref_dist = pref_dist[(pref_dist['X/'+resi_name+'/1/P1_A/SER/199/OG'] < 5.0) & (pref_dist['A/HIS/437/ND1_A/SER/199/HG'] < 4.0 ) & (pref_dist['X/'+resi_name+'/1/O1_A/HIS/437/HD1'] < 4.0)]
	print(pref_dist['Molecule_name'].values)
	'''
	'''
	fig1 =plt.figure(figsize=(12,8))
	ax = fig1.add_subplot(111)
	ax.scatter(x,z,color = 'g', s = 20, marker ='o' )
	ax.set_xlim(2.75, 5.1)
	ax.set_ylim(-500, -250)
	#ax.axis(xmin = 2.75, xmax = 5.25) # this works too
	ax.set_title('O__P dist vs Score $y=x^2$') ## just added some letex
	ax.set_xlabel('O__P dist')
	ax.set_ylabel('Score')
	ax.grid()
	plt.show()


	fig2 =plt.figure(figsize=(12,8))
	ax2 = fig2.add_subplot(111, projection = '3d')
	#ax2 = fig.gca(projection='3d')
	ax2.scatter(x,y,z)
	ax2.set_xlim3d(2.75, 5.1)
	#ax.set_ylim3d(0, 10)
	#ax.set_zlim3d(0, 10)

	ax2.set_xlabel('O__P dist')
	ax2.set_ylabel('hisN__serH')
	ax2.set_zlabel('Score')
	#ax2.set_xslim()
	plt.show()

	'''
"""
###########################################################################################
def menu_x(ad1, ad2, ad3, ad4, score1, file1, prt_inf1, lig11, lig21, file_count1, time_st = '0'):
	
	
	
	dt_fms = impnt_dist_vs_score(prt_inf1,lig11, lig21)
	
	lig_water_dist_ = dt_fms[0]
	his_water_dist_ = dt_fms[1]
	ser_water_dist_ = dt_fms[2]
	his_glu_dis_df_ = dt_fms[3]
	his_ser_dis_df_ = dt_fms[4]
	select_impt_dis_df_ = dt_fms[5] ## Ligand to H_S_G
	hevy_lig_dis_df_ = dt_fms[6]
	
	lig_wat_nam = file1[:61]+'_lig_wt_impt_dis.csv'
	his_wat_nam = file1[:61]+'_his_wt_impt_dis.csv'
	ser_wat_nam = file1[:61]+'_ser_wt_impt_dis.csv'
	lig_ser__his_glu_nam = file1[:61]+'_lig_S_H_G_wt_impt_dis.csv'
	
	emp_df =pd.DataFrame()
	
	accumulative_df_L_W = ad1
	accumulative_df_H_W = ad2
	accumulative_df_S_W = ad3
	accumulative_df_L_H_S_G = ad4
	
	
	#### some of following parallen jobs are working but slower
	#args_to_pass = (
	#		[lig_water_dist_, [emp_df], file1, lig_wat_nam, score1, time_st, file_count1, accumulative_df_L_W],\
	#		[his_water_dist_, [emp_df], file1,	his_wat_nam, score1, time_st, file_count1, accumulative_df_H_W],\
	#		[ser_water_dist_, [emp_df], file1,	ser_wat_nam, score1, time_st, file_count1, accumulative_df_S_W],\
	#		[select_impt_dis_df_, [his_glu_dis_df_, his_ser_dis_df_], file1, lig_ser__his_glu_nam,\
	#		score1, time_st, file_count1, accumulative_df_L_H_S_G])
	'''
	p = multiprocessing.Pool(2)
	data = p.starmap(seperate_process, args_to_pass)
	p.close()
	p.terminate()
	p.join()
	
	#print(type(data), len(data))
	
	accumulative_df_L_W = data[0]
	accumulative_df_H_W = data[1]
	accumulative_df_S_W = data[2]
	accumulative_df_L_H_S_G = data[3]
	'''
	
	'''
	q1 = Queue()
	p1 =Process(target=seperate_process, args=(lig_water_dist_, [emp_df], file1, lig_wat_nam, score1, time_st, file_count1, accumulative_df_L_W, q1,))
	q2 = Queue()
	p2 =Process(target=seperate_process, args=(his_water_dist_, [emp_df], file1,	his_wat_nam, score1, time_st, file_count1, accumulative_df_H_W, q2,))
	q3 = Queue()
	p3 =Process(target=seperate_process, args=(ser_water_dist_, [emp_df], file1,	ser_wat_nam, score1, time_st, file_count1, accumulative_df_S_W, q3,))
	q4 = Queue()
	p4 =Process(target=seperate_process, args=(select_impt_dis_df_, [his_glu_dis_df_, his_ser_dis_df_], file1, lig_ser__his_glu_nam,\
			score1, time_st, file_count1, accumulative_df_L_H_S_G, q4,))
	
	#p1.Daemon = True
	
	p1.start()	
	accumulative_df_L_W = q1.get()
	p1.join()
	p2.start()	
	accumulative_df_H_W = q2.get()
	p2.join()
	p3.start()	
	accumulative_df_S_W = q3.get()
	p3.join()
	p4.start()	
	accumulative_df_L_H_S_G = q4.get()
	p4.join()
	
	'''
	
		
		
	accumulative_df_L_W = seperate_process(lig_water_dist_, [emp_df], file1,\
			lig_wat_nam, score1, time_st, file_count1, accumulative_df_L_W)
	accumulative_df_H_W = seperate_process(his_water_dist_, [emp_df], file1,\
			his_wat_nam, score1, time_st, file_count1, accumulative_df_H_W)
	accumulative_df_S_W = seperate_process(ser_water_dist_, [emp_df], file1,\
			ser_wat_nam, score1, time_st, file_count1, accumulative_df_S_W)
	accumulative_df_L_H_S_G = seperate_process(select_impt_dis_df_,\
			[his_glu_dis_df_, his_ser_dis_df_], file1, lig_ser__his_glu_nam,\
			score1, time_st, file_count1, accumulative_df_L_H_S_G)	
			
	return accumulative_df_L_W, accumulative_df_H_W, accumulative_df_S_W, accumulative_df_L_H_S_G, hevy_lig_dis_df_
	
	
	
###########################################################################################
def menux(ad1, ad2, ad3, ad4, score1, file1, prt_inf1, lig11, lig21, file_count1, time_st = '0'):

	dt_fms = impnt_dist_vs_score(prt_inf1,lig11, lig21)
	
	lig_water_dist_ = dt_fms[0]
	his_water_dist_ = dt_fms[1]
	ser_water_dist_ = dt_fms[2]
	his_glu_dis_df_ = dt_fms[3]
	his_ser_dis_df_ = dt_fms[4]
	select_impt_dis_df_ = dt_fms[5] ## Ligand to H_S_G
	
	lig_wat_nam = file1+'_lig_wt_impt_dis.csv'
	his_wat_nam = file1+'_his_wt_impt_dis.csv'
	ser_wat_nam = file1+'_ser_wt_impt_dis.csv'
	lig_ser__his_glu_nam = file1+'_lig_S_H_G_wt_impt_dis.csv'
	
	emp_df =pd.DataFrame()
	
	accumulative_df_L_W = ad1
	accumulative_df_H_W = ad2
	accumulative_df_S_W = ad3
	accumulative_df_L_H_S_G = ad4
	
			
	accumulative_df_L_W = seperate_process(lig_water_dist_, [emp_df], file1,\
			lig_wat_nam, score1, time_st, file_count1, accumulative_df_L_W)
	accumulative_df_H_W = seperate_process(his_water_dist_, [emp_df], file1,\
			his_wat_nam, score1, time_st, file_count1, accumulative_df_H_W)
	accumulative_df_S_W = seperate_process(ser_water_dist_, [emp_df], file1,\
			ser_wat_nam, score1, time_st, file_count1, accumulative_df_S_W)
	accumulative_df_L_H_S_G = seperate_process(select_impt_dis_df_,\
			[his_glu_dis_df_, his_ser_dis_df_], file1, lig_ser__his_glu_nam,\
			score1, time_st, file_count1, accumulative_df_L_H_S_G)	
			
	return accumulative_df_L_W, accumulative_df_H_W, accumulative_df_S_W, accumulative_df_L_H_S_G
	
###########################################################################################	
	
	
def seperate_process(main_df_get_apeended, df_list_to_append, file__, file_to_wirte_nam, score__, time_st__, file_count__, accumulative_df):
	
	transposed_df = df_append_transpose(main_df_get_apeended, df_list_to_append, file__, score__, time_st__, file_count__)
	accumulative_df = df_concatenate(accumulative_df, transposed_df)
	#accumulative_df.to_csv(file_to_wirte_nam, index=False)
	#q.put(accumulative_df)
	return accumulative_df
	#'''
	#print("Time Menu2: ", datetime.now()-t2)
	
###########################################################################################
	
def menu3_multiproc_pool(select_dist_df3_, one_pdb_,file_nm3_, score3_, file_count3_,time_stp_):
						################################# option 1 make new function here to multi process
	pdb_inf = get_pdb_info_into_df(one_pdb_, file_nm3_)
	
	prt_inf3_ = pdb_inf[0] ## only protein info ina df
	lig13_ = pdb_inf[1]['OP'] ## only the OP info
	lig23_ = pdb_inf[1]['WAT'] ## only the water info
	pro_lig_all3_ = pdb_inf[4] ## get all (protein and ligand info)
	
	select_dist_df3_ = menu_1(select_dist_df3_, score3_, file_nm3_, prt_inf3_, lig13_, lig23_, file_count3_,time_stp_ )
	return select_dist_df3_
					#####################################################################
###########################################################################################

#########    MAIN PROGRAM     #############################################################
## Take the directory path where the necessary files stored
## SELECT THE MENU TO IMPLEMENT
def main_finc():
	list_menu = menu_list()
	totstarttime = datetime.now()
	
	
	
	if list_menu[0] == True or list_menu[1] == True:
		path = input("Enter directory path to .pdbs: ")
		file_names = os.listdir(path)
	
		#ligand_id = input("Enter main ligand identifier: e.g. DCL: ")
	
		## "os.listdir(path)" gets the file names in the directory and stored in a list
	
		is_new_prot_run = True
		is_new_water_run = True
	
		select_dist_df = pd.DataFrame(columns=['Molecule_name','Score'])
		selected_atom_pair = []
		col_seq = np.ndarray(shape=1)
		file_count = 0
	
		## consider a file by file to do calculations
		for file in file_names:
			t2 = datetime.now()
	
	
			file_score = None
			if file.endswith(".pdb"): ## only .pdb files
				file_handle = open(os.path.join(path,file), "r") ## make connection to a file
				file_data = list(file_handle) ## store data of the into a list
	
				## get all the protein ligands info
				pdb_info_in_lists_1 = get_pdb_info_into_df(file_data, file)
	
				prt_inf = pdb_info_in_lists_1[0] ## only protein info ina df
				lig1 = pdb_info_in_lists_1[1]['OP'] ## only the OP info in df
				#if len(pdb_info_in_lists_1[1]) > 1:
	
				lig2 = pdb_info_in_lists_1[1]['WAT'] ## only the water info
				#print (lig2)
				#else:
				#lig2 = None
				score =pdb_info_in_lists_1[2]
				pro_lig_all = pdb_info_in_lists_1[4] ## get all (protein and ligand info)
	
				## get a center of mass or center of molecule roughly
				lig_centr_coord = get_aprox_centerof_molecule(lig1)
	
				## print important calculate and print important distances form ligand to protein and watres
				if list_menu[0]:
					menu_0(prt_inf, lig1, lig2, score, lig_centr_coord,is_new_prot_run, is_new_water_run, file)
					is_new_prot_run = False
					is_new_water_run = False
	
				if list_menu[1]:
					file_count += 1
					select_dist_df = menu_1(select_dist_df, score, file, prt_inf, lig1, lig2, file_count)
	
	
				file_handle.close()
	
	
	if list_menu[2]:
		menu_2()
	
	
	
	if list_menu[3] == True:
		print('All the concatenated pdb files and energy files should be under MAIN directory name')
		print('PDB files should have extention ".pdb" and energy files should be ".edat"')
		print('Make sure to rename the ENERGY file as with the name with corresponding PDB file')
		mainFolder = os.path.abspath(input('Enter the main directory name that has all the pdbs and edats sub folders: \n'))
		isEnergyfile_ = True
		pdb_names_list_ = sorted(read_sub_folders(mainFolder, '', 'pdb'))
		energy_names_list_ = sorted(read_sub_folders(mainFolder, '', 'edat'))
		
		processor_num = 5
		
		path_to_energ_file_ = ''
		
		if len(energy_names_list_) == 0:
			isEnergyfile_ = False
		#	print('Followings are the PDB names and EDAT names')
		argsmt = []
		for i in (range(len(pdb_names_list_))):
			pdb_path_with_file_ = pdb_names_list_[i]
			argsmt.append([pdb_path_with_file_, path_to_energ_file_, isEnergyfile_])
			#print(pdb_names_list_[i])
			#print(len(energy_names_list_))
		print('length argsmt ',len(argsmt))
		#pool_time = datetime.now()
		#pool_size = multiprocessing.cpu_count() * 2
		#pool = Pool(processes = processor_num)
		#pool = Pool(processes=pool_size, initializer=start_process, maxtasksperchild=2,)
		#pool.map(main_process, argsmt)
		#pool.close()
		#pool.join()
		
		#end_poolTime =  datetime.now() - pool_time
		
		#reg_time = datetime.now()
		#for arg in argsmt:
			
			#main_process(arg)
			
		#end_regTime =  datetime.now() - reg_time
		
		proc_time = datetime.now()
		p1 = Process(target=main_process, args=(argsmt[0],))
		p2 = Process(target=main_process, args=(argsmt[1],))
		p3 = Process(target=main_process, args=(argsmt[2],))
		p4 = Process(target=main_process, args=(argsmt[3],))
		#main_process(argsmt[2])
		p1.start() ## start the job
		p2.start() 
		p3.start() 
		p4.start()
		p1.join() ## wait until others are done
		p2.join() 
		p3.join()
		p4.join()
		proc_Time =  datetime.now() - proc_time
		
		#print('Pool tile', end_poolTime)
		#print('Regular time', end_regTime)
		print('Process time', proc_Time)
		print('Done ')
		
			#main_process(pdb_path_with_file_, path_to_energ_file_, isEnergyfile_)
		#		
		#		pdb_path_with_file = pdb_names_list_[i]
		#		pdb_head, pdb_tail = os.path.split(pdb_path_with_file)
		
		'''
		#args_for_processors = []
		#pdb_name = input('Enter concatenate PDB-file-name that derived from NAMD DCD-output-file:  '+ '\n')
		#pdb_name = 'AChE_DCV_08167_wb40A_neutral_pcm_chargesp1_selected_atoms_.pdb'
		energy_f_name = input('Enter file name of the energy (press enter if NO energy file available): ')
		
		curr_dir = os.getcwd()
		isEnergyfile = True
		
		if energy_f_name.strip() == '':
			print(energy_f_name)
			isEnergyfile = False
		file_names3 = os.listdir(curr_dir)
	
		if pdb_name in file_names3:
			pdb_path_with_file = os.path.join(curr_dir,pdb_name)
			path_to_energ_file = os.path.join(curr_dir,energy_f_name)
	
	
		else:
			print('Wanted PDB file or energy file is not found in current working directory'+ '\n')
			path_pdb = input('Enter path to derived PDB and energy file file (both should be in same directory):  '+ '\n')
			#path_pdb = "F:\\OneDrive\\MD_work\\MD_analysis\\get_data_using_CATDCD\\temp"
			print (path_pdb)
			pdb_path_with_file = os.path.join(path_pdb,pdb_name)
			
			if isEnergyfile:  ## in the cases of no energy file exists
				
				path_to_energ_file = os.path.join(path_pdb,energy_f_name)
	
			
			need_submenu = 'n' 
			'''
			#######################'input('''Sub menu for MD analysis:
					#Do you need 'important distances' as a seperate file along with important distance matrix?: (y/n)
					#''')
		## read energy file to panda
		#energy_df = pd.Dataframe()
		#main_process(pdb_path_with_file_, path_to_energ_file_, isEnergyfile_)	
		
		
##############################################################################
## this finction gives list of names based on given partial information
## e.g. partial file name and main directory	
## get a list of dcd names in subdirectories too
def read_sub_folders(base_path, base_nam, extention):
	names_list = []
	#print()
	#print('The '+extention+' file name/s that going to be process:')
	#print('--------------------------------------------------------')
	## get all the file names in subdirectories
	#path_with_file = os.path.join(base_path, base_nam+'*.'+extention)
	#for filename in glob.iglob(base_path+'/**/'+base_nam+'*.'+extention, recursive=True):
	for filename in glob.iglob(base_path+'/'+base_nam+'*.'+extention):
	
		names_list.append(filename)
		#log_head, log_tail = os.path.split(filename)
		#print(log_tail)
	#print('--------------------------------------------------------')
	#print()
	return names_list

##############################################################################
			
def main_process(args):
	pdb_path_with_file = args[0]
	path_to_energ_file = args[1]
	isEnergyfile = args[2]
	
	if True:

		## In this part, the contatenated pdb file that extracted from tne NAMD dcd file (used VMD for this),
		## is used. The concatenated pdb file break into individuals for each time step,
		## and further analysis carry out. Some of the part of here is repeting elsware (get_pdb_info_into_df), but
		## for the intergrity of the program, repetition is allowed.
		is_new_prot_run3 = True
		is_new_water_run3 = True
		#select_dist_df3 = pd.DataFrame(columns=['Molecule_name','Score'])
		
		accumulative_df_L_W_ = pd.DataFrame(columns=['Molecule_name','Score'])
		accumulative_df_H_W_ = pd.DataFrame(columns=['Molecule_name','Score'])
		accumulative_df_S_W_ = pd.DataFrame(columns=['Molecule_name','Score'])
		accumulative_df_L_H_S_G_ = pd.DataFrame(columns=['Molecule_name','Score']) 
		
		distance_count_dict = {}
		hevy_lig_dis_df_ = pd.DataFrame()
		file_count3 = 0
	
		one_pdb = list()
		time_stp = 0

	
		need_submenu = 'N'


		if isEnergyfile:
			energy_df = pd.read_csv(path_to_energ_file, delim_whitespace=True )
	
		#print (energy_df)
		#pdb_name = pdb_path_with_file.split('\\')[-1][:-4]  ## change this accordingly for linux
		#pdb_name = re.split('[/\\]', pdb_path_with_file)
		base_path, pdb_name = os.path.split(pdb_path_with_file)
		pdb_name =pdb_name[:-4]
		#print(pdb_name)
		#print(base_path)
		#print(os.path.join(base_path, pdb_name+'_lig_wt_impt_dis.csv'))
		
		#pdb_file_handle = open(path_with_file, 'r')
		with open(pdb_path_with_file, 'r') as pdb_file_handle:
			for i, linee in enumerate(pdb_file_handle):
	
				linee = linee.rstrip()
				if linee == '': continue
	
	
				if linee.startswith('ATOM') and not linee.startswith('END'):
					one_pdb.append(linee)
	
	
				elif linee.startswith('END'):
					
					if isEnergyfile:
						score3 = [energy_df.ix[time_stp, 'Nonbond'], 'NA'] ## here 'NA' is just for resolve 'interface_delta_enegy'
					else:
						score3 = ['NA', 'NA']
					
					## to resolve later error
					file_nm3 = ('MD_'+pdb_name+'_time_step_'+str(time_stp)+'.pdb')
					
					
					
					
					#####################################################################################

					## get all the protein ligands info
					pdb_inf = get_pdb_info_into_df(one_pdb, file_nm3)
	
					one_pdb = list()
					#print(pdb_file_handle.closed)
					#print(pdb_inf[0][pdb_inf[0]['res_nm'] == 'HIS' ])
					#break
	
	
	
					prt_inf3 = pdb_inf[0] ## only protein info ina df
					lig13 = pdb_inf[1]['OP'] ## only the OP info
					lig23 = pdb_inf[1]['WAT'] ## only the water info
					pro_lig_all3 = pdb_inf[4] ## get all (protein and ligand info)

					
					
					#if list_menu[1]:
					file_count3 += 1
					
					accumulative_df_L_W_, accumulative_df_H_W_, accumulative_df_S_W_, accumulative_df_L_H_S_G_, hevy_lig_dis_df__ =\
						menu_x(accumulative_df_L_W_, accumulative_df_H_W_,\
						accumulative_df_S_W_, accumulative_df_L_H_S_G_, score3,\
						file_nm3, prt_inf3, lig13, lig23, file_count3,time_stp )
					print('file count : ',file_count3)
					#print(len(accumulative_df_L_W_))
					
					hevy_lig_dis_df_temp = hevy_lig_dis_df__[['dist_name', 'dist']][hevy_lig_dis_df__['dist'] < 4.0].copy()
					hevy_lig_dis_df_ = (hevy_lig_dis_df_.append(hevy_lig_dis_df_temp, ignore_index=True))
					#distance_count_dict[]
					#write_to_csv(lig_water_dis_df, f_nam0[10:30]+'_lig_water.csv',f_nam0,is_new_water_run0,score0)
					
					time_stp += 1
					## get all the protein ligands info
					
					####################################################################################################
					
					
					
					'''############################################################################################################
					#for i in range(1, processor_num +1):
					#while count < processor_num:
					if count <= processor_num : 
						args_for_processors.append((select_dist_df3, one_pdb,file_nm3, score3, file_count3, time_stp))

						
						
						if count == processor_num: # or end of file:
							#print(len(args_for_processors))
							p = pool.starmap(menu3_multiproc_pool, args_for_processors)
							args_for_processors = []
							count = 0
							
							for item in p:
								#print(item['Time_step'])
								select_dist_df3 = select_dist_df3.append(item, ignore_index=True)
								#select_dist_df3.set_index('dist_name', inplace=True)	
								select_dist_df3.to_csv(file_nm3[10:30]+'_impotant_dist.csv', index=False)	
							print(len(select_dist_df3))
						count += 1
						#break					
					time_stp += 1
					file_count3 += 1
					#(select_dist_df3_, one_pdb_,file_nm3_, score3_, file_count3_,time_stp_)
					
						
					'''#################################################################################################
					#syntax for above: L = pool.starmap(func, [(1, 1), (2, 1), (3, 1)])
					
					'''################################# option 1 make new function here to multi process
					pdb_inf = get_pdb_info_into_df(one_pdb, file_nm3)
	
					one_pdb = list()
					#print(pdb_file_handle.closed)
					#print(pdb_inf[0][pdb_inf[0]['res_nm'] == 'HIS' ])
					#break
	
	
	
					prt_inf3 = pdb_inf[0] ## only protein info ina df
					lig13 = pdb_inf[1]['OP'] ## only the OP info
					lig23 = pdb_inf[1]['WAT'] ## only the water info
					pro_lig_all3 = pdb_inf[4] ## get all (protein and ligand info)
	
					
	

	
					#if list_menu[1]:
					
					select_dist_df3 = menu_1(select_dist_df3, score3, file_nm3, prt_inf3, lig13, lig23, file_count3,time_stp )
	
					'''#####################################################################

					
					
					
					## get a center of mass or center of molecule roughly
	
					
					if need_submenu == 'y':
						pdb_inf = get_pdb_info_into_df(one_pdb, file_nm3)
						prt_inf3 = pdb_inf[0] ## only protein info ina df
						lig13 = pdb_inf[1]['OP'] ## only the OP info
						lig23 = pdb_inf[1]['WAT'] ## only the water info
						pro_lig_all3 = pdb_inf[4] ## get all (protein and ligand info)
						lig_centr_coord3 = get_aprox_centerof_molecule(lig13)
						## print important calculate and print important distances form ligand to protein and watres
						### uncomment to get written important distances
						menu_0(prt_inf3, lig13, lig23, score3, lig_centr_coord3,is_new_prot_run3, is_new_water_run3,file_nm3)
						is_new_prot_run3 = False
						is_new_water_run3 = False
					
					
					one_pdb = list() ## empty one_pdb
		
		## counting the unique distances
		#dist_count = hevy_lig_dis_df_['dist_name'].groupby('dist_name').size()
		#print(dist_count.head(20))
		
		#dist_count = df[['col1', 'col2', 'col3', 'col4']].groupby(['col1', 'col2']).agg(['mean', 'count'])
		dist_count_1 = hevy_lig_dis_df_[['dist_name', 'dist']].groupby('dist_name').agg('count').reset_index()
		print(dist_count_1.head(20))
		
		
		
		lig_wat_nam = os.path.join(base_path, pdb_name+'_lig_wt_impt_dis.csv')
		his_wat_nam = os.path.join(base_path, pdb_name+'_his_wt_impt_dis.csv')
		ser_wat_nam = os.path.join(base_path, pdb_name+'_ser_wt_impt_dis.csv')
		lig_ser__his_glu_nam = os.path.join(base_path, pdb_name+'_lig_S_H_G_wt_impt_dis.csv')
		hevy_lig_dis_nam = os.path.join(base_path, pdb_name+'_close_contact_names.csv')
		hevy_lig_dis_count_nam = os.path.join(base_path, pdb_name+'_close_contact_nam_count.csv')
		#lig_wat_nam = pdb_name+'_lig_wt_impt_dis.csv'
		#his_wat_nam = pdb_name+'_his_wt_impt_dis.csv'
		#ser_wat_nam = pdb_name+'_ser_wt_impt_dis.csv'
		#lig_ser__his_glu_nam = pdb_name+'_lig_S_H_G_wt_impt_dis.csv'
		
				
		
		
		
		
					
		write_to_csvfile(accumulative_df_L_W_, lig_wat_nam)
		write_to_csvfile(accumulative_df_H_W_, his_wat_nam)
		write_to_csvfile(accumulative_df_S_W_, ser_wat_nam)
		write_to_csvfile(accumulative_df_L_H_S_G_, lig_ser__his_glu_nam)
		write_to_csvfile(hevy_lig_dis_df_, hevy_lig_dis_nam)
		write_to_csvfile(dist_count_1, hevy_lig_dis_count_nam)
		

	
	#print("Total_time: ", datetime.now()-totstarttime)
def write_to_csvfile(dff, fname):
	dff.to_csv(fname, index=False)
	

if __name__ == '__main__':
	main_finc()
	
	#pool = multiprocessing.Pool(processes = 3)
	#pool.map(main_finc())
				
				
				