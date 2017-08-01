# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 23:08:53 2017
This program get in selected time steps that have lowest RMS values with water,
and extract corresponding coordinates to the time steps in pdb fromat.
@author: Lasantha
"""
## 
#from sas7bdat import SAS7BDAT


#import subprocess
import os
import glob
import pandas as pd
import math
from datetime import datetime
#import platform
#import threading
#import multiprocessing
#subprocess.call('dir', shell=True)
## for linux
## subprocess.call('ls', shell=True)
#output = subprocess.check_output('dir', shell=True)
#print(output)
#subprocess.call('catdcd', shell=True)
#subprocess.Popen('cmd')
#subprocess.Popen(ROBCOPY)



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
	for filename in glob.iglob(base_path+'/**/'+base_nam+'*.'+extention, recursive=True):
		names_list.append(filename)
	#	log_head, log_tail = os.path.split(filename)
	#	print(log_tail)
	#print('--------------------------------------------------------')
	print()
	return names_list

##############################################################################

def read_SAS_dst(sas_dst):
	with SAS7BDAT(sas_dst) as sas:
		sas_df = sas.to_data_frame()
	return sas_df
		
##############################################################################		

def write_to_pdb(pdb_in_lst, pdb_nm_to_write):
	
	with open(pdb_nm_to_write, 'w') as pdb:
		for lin in pdb_in_lst:
			pdb.write(lin+'\n')
		
		pdb.write('END')
###########################################################################################
		
##############################################################################		

def write_to_csv(csv_in_lst, csv_nm_to_write):
	
	with open(csv_nm_to_write, 'w') as csv:
		line_num = 1
		for lin in csv_in_lst:
			if line_num == 1:
				csv.write(lin+'\n')
				line_num == 2
			else:               
				csv.write(lin)
              
		
		#csv.write('END')
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
	prot_found = False
	liga_found = False
	score_found = False
	other_lig_found = False

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
	
###########################################################################################
	
def cal_distance4(row,liginf ):
	
	t1 = datetime.now()

	if len(row) != 0 and len(liginf) != 0:
		row1 = (row.values)[0]
		lig = (liginf.values)[0]
		if len(liginf) > 2:
			lig = (liginf.values)
		
		## create a name to distance
		dis_nm = lig[8]+'/'+lig[4]+'/'+lig[5]+'/'+lig[6]+'/'+lig[7]+'_'+row1[8]+'/'+row1[4]+'/'+row1[5]+'/'+row1[6]+'/'+row1[7]
		#print (dis_nm)
		## assigning coordinates to x, y ,z
		x1 = float(row1[0])
		y1 = float(row1[1])
		z1 = float(row1[2])
		x2 = float(lig[0])
		y2 = float(lig[1])
		z2 = float(lig[2])
		#row[['chID','res_nm','resID','atm_type']]
	
		## calculate distance
		calcd_distns = round(math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2),3)
	
		#print("cal_distance3: pandas" , datetime.now()-t1)
		return [dis_nm, calcd_distns]
	else:
		return ['Epmty',100]

###########################################################################################
###########################################################################################


def give_reqrd_dist_and_pdbs(prt_inf_, lig1_, lig2_,atm_list_, pdb_count_, base_nm_, tme_stp):
	
	resi_list = ['437','199', '327']
	atm_lst = ['OG', 'OE1', 'OE2', 'ND1','ND2','NE1', 'NE2', 'OH'] ##***************** modify loc1
	hevy_atoms = ['P','S','O','C']

	only_impotant_atms = prt_inf_[prt_inf_['resID'].isin(resi_list) & prt_inf_['atm_type'].isin(atm_lst)]
	only_impotant_atms = only_impotant_atms[['x_cor','y_cor', 'z_cor','atmID','chID', 'res_nm', 'resID', 'atm_type', 'seg_id', 'atm_sym']]
	
	dist_dis = {'Time_step': tme_stp}
	
	
	lig1_ = lig1_[['x_cor','y_cor', 'z_cor','atmID','chID', 'res_nm', 'resID', 'atm_type', 'seg_id', 'atm_sym']]
	lig2_ = lig2_[['x_cor','y_cor', 'z_cor','atmID','chID', 'res_nm', 'resID', 'atm_type', 'seg_id', 'atm_sym']]
	
	
	## get important his atoms
	his_ND = only_impotant_atms[only_impotant_atms['resID'].isin(['437']) & prt_inf_['atm_type'].isin(['ND1', 'ND2'])]	
	his_NE = only_impotant_atms[only_impotant_atms['resID'].isin(['437']) & prt_inf_['atm_type'].isin(['NE1', 'NE2'])]
	## get important ser atoms
	wat_O = lig2_[lig2_['atm_sym'].isin(['O'])]

	ser_O =  only_impotant_atms[only_impotant_atms['resID']=='199']
	tyr_O =  only_impotant_atms[only_impotant_atms['resID']=='327']
	#lig2_only_water_O = lig2_[lig2_['atm_sym'].isin(['O'])]
	hist_list = str()
	wat_ser = str()
	wat_tyr = str()
	

	if lig1_.ix[lig1_.index[1],'res_nm'] in ['ACH']: ## WORK THIS LATER
		## selected ligand atoms for ACH
		#lig_inf3 = lig1_[lig1_['atm_type'].isin(['O1','O2','C6','N1','C1', 'C2','C4'])] 
		lig_P1 = lig1_[lig1_['atm_type'].isin(['C6'])] 
		#lig_inf_all = lig1_[lig1_['atm_sym'].isin(hevy_atoms)]
		lig_O1S1 = lig1_[lig1_['atm_type'].isin(['O1'])]  ## change to O2 if carbonyl O needed
		
	else:
		lig_P1 = lig1_[lig1_['atm_type'].isin(['P1'])]
		
		if lig1_.ix[lig1_.index[1],'res_nm'] in['TMT']:
			lig_O1S1 = lig1_[lig1_['atm_type'].isin(['S1'])]
		else:
			lig_O1S1 = lig1_[lig1_['atm_type'].isin(['O1'])]
	
	serO_OP_P1 = cal_distance4(ser_O, lig_P1)	
	hisND_OP_O1S1 = cal_distance4(his_ND, lig_O1S1)
	hisNE_OP_O1S1 = cal_distance4(his_NE, lig_O1S1)
	serO_HisND = cal_distance4(ser_O, his_ND)
	serO_HisNE = cal_distance4(ser_O, his_NE)
	tyrO_OP_P1 = cal_distance4(tyr_O, lig_P1)
		
	dist_dis = {serO_OP_P1[0] : serO_OP_P1[1], tyrO_OP_P1[0] : tyrO_OP_P1[1], hisND_OP_O1S1[0] : hisND_OP_O1S1[1], hisNE_OP_O1S1[0] : hisNE_OP_O1S1[1],\
				serO_HisND[0] : serO_HisND[1],serO_HisNE[0] : serO_HisNE[1], 'Time_step': tme_stp}
                				
	cutoff1 = 3.7
	cutoff2 = 3.7
	cutoff3 = 3.7
		
	if float(serO_OP_P1[1]) <= cutoff1:
		
		if float(serO_HisND[1])<= cutoff1 or float(serO_HisNE[1]) <= cutoff1:
			if float(hisND_OP_O1S1[1])<= cutoff1 or float(hisNE_OP_O1S1[1]) <= cutoff1:					
				#print('Hist', tme_stp)
				if tme_stp > 30:
					write_to_pdb(atm_list_, base_nm_[:-3]+'_hist_'+str(tme_stp)+'.pdb')
					hist_list = (str(tme_stp)+','+str(serO_OP_P1[1])+','+str(serO_HisND[1])+','+ str(serO_HisNE[1])+','+str(hisND_OP_O1S1[1])+','+ str(hisNE_OP_O1S1[1]))+ '\n'
					
					#print(serO_OP_P1[1],serO_HisND[1], serO_HisNE[1],hisND_OP_O1S1[1], hisNE_OP_O1S1[1] )
	for ligline in wat_O.index:
		wt_O_atm = wat_O.ix[ligline, ['x_cor','y_cor', 'z_cor','atmID','chID', 'res_nm', 'resID', 'atm_type', 'seg_id', 'atm_sym']]
			
		wat_OP_O = cal_distance4(lig_O1S1, wt_O_atm)
		wat_Ser_O = cal_distance4(ser_O, wt_O_atm)
		wat_Tyr_O = cal_distance4(tyr_O, wt_O_atm)
		
		dist_dis[wat_OP_O[0]] = wat_OP_O[1]
		dist_dis[wat_Ser_O[0]] = wat_Ser_O[1]
		dist_dis[wat_Tyr_O[0]] = wat_Tyr_O[1]
			
		#print(serO_OP_P1[1], wat_OP_O[1],wat_Ser_O[1] )
		if float(serO_OP_P1[1]) <= cutoff2 and float(wat_OP_O[1])<= cutoff2 and float(wat_Ser_O[1])<= cutoff2:
			#print('Water',tme_stp)
			if tme_stp > 30:
				write_to_pdb(atm_list_, base_nm_[:-3]+'_water_SER199_'+str(tme_stp)+'.pdb')
				wat_ser = wat_ser + (str(tme_stp)+','+str(serO_OP_P1[1])+','+ str(wat_OP_O[1])+','+str(wat_Ser_O[1])) + '\n'
				#print(float(serO_OP_P1[1]), float(wat_OP_O[1]),float(wat_Ser_O[1]))
				
		if float(tyrO_OP_P1[1]) <= cutoff3 and float(wat_OP_O[1])<= cutoff3 and float(wat_Tyr_O[1])<= cutoff3:
			#print('Water',tme_stp)
			if tme_stp > 30:
				write_to_pdb(atm_list_, base_nm_[:-3]+'_water_TYR327_'+str(tme_stp)+'.pdb')
				wat_tyr = wat_tyr + (str(tme_stp)+','+str(tyrO_OP_P1[1])+','+ str(wat_OP_O[1])+','+str(wat_Tyr_O[1])) + '\n'
				
				
				
				

	#data_frm = pd.DataFrame([dist_dis])

	return [hist_list,wat_ser, wat_tyr]



###########################################################################################

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

##############################################################################
## this function take the pdb file (set of time steps in one pdb),
## give out a LIST of LIST of pdb files
## main list contains the pdb info for  a time step and the sub list of the main list
## contain the ligand and protein+water(if any) seperated
def read_collective_pdb(colective_pdb_nm, pdb_bas_nm, lignd_id, new_dir):
	t1 = datetime.now()
	
	atm_list = []
	his_lst = ['tme_stp,serO_OP_P1,serO_HisND,serO_HisNE,hisND_OP_O1S1, hisNE_OP_O1S1']
	wat_ser_lst = ['tme_stp,serO_OP_P1, wat_OP_O,wat_Ser_O']
	wat_tyr_lst = ['tme_stp,tyrO_OP_P1, wat_OP_O,wat_tyr_O']
	pdb_count = 0
	distan_df = pd.DataFrame()
	base_pth, log_tail = os.path.split(colective_pdb_nm[0])
	#base_nm = os.path.join(base_pth, pdb_bas_nm)
	base_nm = os.path.join(new_dir, pdb_bas_nm)
	
	for pdb_nm in colective_pdb_nm:	
		with open(pdb_nm, 'r') as coll_pdb:
				
			for line in coll_pdb:
				line = line.rstrip()
				if line == '': continue
				
				if line.startswith('ATOM'):	
					atm_list.append(line)
	
				if line.startswith('END'):
					pdb_count += 1
					print(pdb_count)
					pdb_df = get_pdb_info_into_df(atm_list, 'fil_name', lignd_id)
					
					wat_df = pdb_df[1]['WAT']
					prt_df = pdb_df[0]
					lig_df = pdb_df[1]['OP']
					#pdb_nm_to_be = base_nm[:-3]+'_tracebk_ts_'+str(pdb_count)+'.pdb'
					
					one_line_lst = give_reqrd_dist_and_pdbs(prt_df, lig_df, wat_df,atm_list, pdb_count, base_nm, pdb_count)
					#print(distan_df)
					#print (len(), len(pdb_df[1]['OP']), len(pdb_df[1]['WAT']))
					#####distan_df = distan_df.append(one_line_df, ignore_index=True)
					if len(one_line_lst[0]) != 0:
						his_lst.append(one_line_lst[0])
					if len(one_line_lst[1]) != 0:
						wat_ser_lst.append(one_line_lst[1])
					if len(one_line_lst[2]) != 0:
						wat_tyr_lst.append(one_line_lst[2])
					
					#if pdb_count in need_pdb_num_lst:
					#print(base_nm[:-3]+'_tracebk_ts_'+str(pdb_count)+'.pdb')
					#write_to_pdb(atm_list, base_nm[:-3]+'_tracebk_ts_'+str(pdb_count)+'.pdb')
					atm_list = []
	
	
	write_to_csv(his_lst, base_nm[:-3]+'_hsit_short_dist_.csv')
	write_to_csv(wat_ser_lst, base_nm[:-3]+'_wat_ser199_short_dist_.csv')
	write_to_csv(wat_tyr_lst, base_nm[:-3]+'_wat_tyr327_short_dist_.csv')
	#distan_df.to_csv(base_nm+'_full_.csv', index=False)
	print ('Total time steps: ', pdb_count)
	print("traceback_pdb: ", datetime.now()-t1)
	
	
	
##############################################################################
def main_func():
	
	
	base_path_to_pdbs = input('Enter base path to the MD extracted pdb file location:\n ') 
	#base_path_to_pdbs ='E:\OneDrive\SAS_work\MD_SAS_analysis\ODM_analysis'
	pdb_base_nm = input('Enetr base name that common to all pdb files: \n') 
	#pdb_base_nm ='ache_odm_0153_wb40A_neutral'
	ligID = input('Enter ligand ID. e.g. OMT: \n')
	#ligID = 'ODM'
	#path_to_SAS = input('Enter base path to "'+ligID+'_tracebk_summary.sas7bdat" SAS data set: \n') 
	#'G:\OneDrive\SAS_practice\MD_SAS_analysis\OP_data\OMT_data\omt_tracebk_summary.sas7bdat'
	#path_n_name_to_SAS = os.path.join(path_to_SAS, ligID+'_tracebk_summary_2.sas7bdat')
	#sas_df = read_SAS_dst(path_n_name_to_SAS)
	
	directory = os.path.join(base_path_to_pdbs, 'react_lyk_pdbs')
	if not os.path.exists(directory):
		os.makedirs(directory)
	
	
	pdb_nm_lst = read_sub_folders(base_path_to_pdbs, pdb_base_nm, 'pdb')
	read_collective_pdb(pdb_nm_lst, pdb_base_nm, ligID, directory)
##############################################################################
	
	
if __name__ == '__main__':
	main_func()
	
	