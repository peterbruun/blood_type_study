# ----------------------------------------------- #
# Snakefile pipeline for the paper: 
# ABO/RhD Blood Groups and Phenome-Wide Disease Incidence in 481,298 Danish Patients
# The code is free to use but please considered citing the paper if you use the code

# -- HOW TO RUN -- #
# The pipeline is run in two steps
# First step estimates which phecodes to include in the analysis
# Second step runs the analysis on all included phecodes

# 1) First time the script is run the follwing snakemake command should be evoked
# run: snakemake -R --until filter_phecodes_and_blocklevels_level3_to_include
# And first run should be set to True
# first_run = True

# 2) At the second run the clean snakemake command should be evoked
# run: snakemake
# And the first run should be set to False
first_run = False

# ----------------------------------------------- #

# Define workdir
workdir: "/users/projects/bloodtype/"
# anaconda3/2021.5

import numpy as np
import pandas as pd
import sys, os, pickle, datetime, glob
import sqlite3

pd.set_option('max_info_columns', 10**10)
pd.set_option('display.max_columns', 1000)
pd.set_option('max_info_columns', 10**10)
pd.set_option('max_info_rows', 10**10)

os.chdir("/users/projects/bloodtype/")


if first_run:
	# Random list of phecodes for first run of script
	phecodes_to_include = ["286.3"]
else:
	# Second run: Can only be run after rule: filter_phecodes_and_blocklevels_level3_to_include
	phecodes_to_include = pd.read_csv("data/processed/included_phecodes.tsv",sep="\t",dtype=str)
	phecodes_to_include = list(phecodes_to_include.Phecodes.values)


# Define congenital/heridatary/paternal diseases and exclude pregnancy diagnosis from age at diagnosis analysis
# PheWAS
non_congenital_phecodes_to_include = phecodes_to_include.copy()
for code in phecodes_to_include:
	if (code[:3] in ["691","356","634","635","636","637","638","639","642","645","646","647","649","651","652","653","654","655","656","657","658","661","663","665","668","669","671","674","676","679","681", \
		"639","637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282"]) | (code[:5] in ["199.4","286.3","306.1","244.5","612.3","286.1","282.9","362.7","520.1","282.8","334.1","364.5"]):
		non_congenital_phecodes_to_include.remove(code)
# No cases for AB blood group thus exclude 728.1 from analysis of age at diagnosis
non_congenital_phecodes_to_include.remove("728.1")


rule all:
	input:
		# Final file
		"data/processed/bloodtype_final.pkl",
		# Table 1
		"results/table1.html",

		# Analysis input files for each code
		#expand("data/processed/PheWAS/{phecode}.tsv",phecode=phecodes_to_include),

		# -- all other vs One -- #
		expand("results/20220531/allVSOne/enter_registry/phewas/{blood_group}/estimates_{phecode}.tsv", blood_group =["A","B","AB","0"],phecode=phecodes_to_include),
		expand("results/20220531/allVSOne/enter_registry/phewas/RhD/RhD_estimates_{phecode}.tsv",phecode=phecodes_to_include),
		"results/20220531/allVSOne/enter_registry/phewas_estimates.tsv",
		
		# Age at diagnosis:
		# linear model
		"results/20220531/allVSOne/age_at_diagnosis/phewas_estimates.tsv",

		# Infiles
		bloodtypes = "/users/secureome/home/projects/bth/personal_folders/vicmuse/BloodType/bloodtypes_cleaned.tsv",
		t_adm = "/users/secureome/home/projects/registries/lpr2bth/2018/lpr/t_adm.tsv",
		t_diag = "/users/secureome/home/projects/registries/lpr2bth/2018/lpr/t_diag.tsv",
		ICD8_10_mapping = "/users/secureome/home/people/petras/base_data/ICD8_to_ICD10.tsv",
		t_person = "/users/secureome/home/projects/registries/lpr2bth/2018/cpr/t_person.tsv",
		phecodes = "/users/secureome/home/people/petras/base_data/Phecode_map_v1_2_icd10_beta.csv",
		phecodes_def = "/users/secureome/home/people/petras/base_data/phecode_def.csv",
		icd10_block_codes = "/users/secureome/home/classification/complete/icd10_eng_diag_chapters_all_update2016.tsv",
		lpr2bth_mapping = "data/raw/valid_lpr2bth.tsv",
		# Outfiles
		phecodes_to_include = "data/processed/included_phecodes.tsv",
		blocklevels_to_include = "data/processed/included_blocklevels.tsv",
		level3_to_include = "data/processed/included_level3.tsv",


rule preprocess_full_lpr_raw_ulrik: 
	input:
		t_adm = rules.all.input.t_adm,
		t_diag = rules.all.input.t_diag,
		bloodtypes = rules.all.input.bloodtypes,
	output:
		lpr = temp("data/processed/full_lpr_ulrik.pkl"),
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*200,
	threads: 1,
	script:
		"scripts/preprocess_lpr_raw.py"


rule convert_ICD8_toICD10:
	input:
		lpr = rules.preprocess_full_lpr_raw_ulrik.output.lpr,
		ICD8_10_mapping = rules.all.input.ICD8_10_mapping,
	output:
		lpr_mapped = temp("data/processed/lpr_ICD8_to_10.pkl"),
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,
	run:
		def convert_ICD8_to_10(df_master):
			# ICD10
			# Select only diagnosis with ICD10 codes
			df_icd10 = df_master[df_master.Diagnosis.str.contains("ICD10")].copy()
			# Strip ICD10 and D in front
			df_icd10["Diagnosis"] = df_icd10["Diagnosis"].str[7:]
			# Sort by cpr and INDDTO
			df_icd10.sort_values(by=['cpr_enc','D_INDDTO'],inplace=True,ascending=True)

			# ICD8
			# Select only diagnosis with ICD8 codes
			df_icd8 = df_master[df_master.Diagnosis.str.contains("ICD8")].copy()
			# Strip ICD8
			df_icd8["Diagnosis"] = df_icd8["Diagnosis"].str[5:]

			# Load mapping file
			icd8_map = pd.read_csv(input["ICD8_10_mapping"],sep="\t",dtype="str")
			# map to ICD8 to ICD10
			df_icd8.rename({"Diagnosis":"ICD-8"},axis=1,inplace=True)
			df_icd8 = pd.merge(df_icd8,icd8_map[["ICD-8","ICD-10 match","Score"]],on="ICD-8",how="left")
			df_icd8["Diagnosis"] = df_icd8["ICD-10 match"]
			# Drop non mapped icd8 codes
			df_icd8.dropna(subset=["Diagnosis"],inplace=True)
			# Drop mapping Scores of 5
			df_icd8 = df_icd8[df_icd8.Score != "5"].copy()
			df_icd8.drop(columns="Score",inplace=True)
			# concatenate icd10 and converted icd8
			df_icd8.drop(["ICD-8","ICD-10 match"],axis=1,inplace=True)
			df_diag = pd.concat([df_icd10,df_icd8],ignore_index=True,sort=False)

			df_diag.to_pickle(output["lpr_mapped"])

			return df_diag

		master = pd.read_pickle(input["lpr"])
		convert_ICD8_to_10(master)


rule phecode_map_to_lpr:
	input:
		lpr_mapped = rules.convert_ICD8_toICD10.output.lpr_mapped,
		phecodes = rules.all.input.phecodes,
	output:
		lpr_phecode = temp("data/processed/lpr_and_phecodes.pkl"),
		icd10_notmapped = temp("data/processed/icd10_notmapped_to_phecodes.tsv"),
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,
	run:
		phecodes = pd.read_csv(input.phecodes,sep=",",dtype="str")
		lpr = pd.read_pickle(input.lpr_mapped)

		# Make sure the length is only 4
		lpr["Diagnosis_phecode"] = lpr.Diagnosis.str[:4]
		# Remove commas to fit danish lpr format
		phecodes["ICD10"] = phecodes.ICD10.str.replace(".","")
		# Merge phecodes phenotype onto LPR
		lpr_phenotype = pd.merge(lpr,phecodes,left_on="Diagnosis_phecode",right_on="ICD10",how="left")
		
		# Map the ones which could be mapped to phecodes of length 4 to phecodes of length 3. eg. C809 -> C80
		mask = lpr_phenotype.PHECODE.isnull()
		lpr2_phenotype = lpr_phenotype.loc[mask,:].copy()
		lpr_phenotype = lpr_phenotype.loc[~mask,:].copy()
		lpr2_phenotype["Diagnosis_phecode"] = lpr2_phenotype.Diagnosis.str[:3]
		# Merge phecodes phenotype onto LPR
		lpr2_phenotype = pd.merge(lpr2_phenotype[["cpr_enc","D_INDDTO","Diagnosis","Diagnosis_phecode"]],phecodes,left_on="Diagnosis_phecode",right_on="ICD10",how="left")
		lpr2_phenotype.PHECODE.isnull().value_counts()

		# Concatenate df mapped with length 4 and length 3
		lpr_phenotype = pd.concat([lpr_phenotype,lpr2_phenotype],ignore_index=True,axis=0)
		lpr_phenotype.drop(["Diagnosis_phecode"],axis=1,inplace=True)
		del lpr, lpr2_phenotype, mask

		# Save the ones not matching
		not_found = pd.Series(lpr_phenotype[lpr_phenotype.PHECODE.isnull()].Diagnosis.unique())
		not_found.to_csv(output.icd10_notmapped,index=False,sep="\t")

		# Drop ICD10 column as it is redundant
		lpr_phenotype.drop(["ICD10"],axis=1,inplace=True)
		# Save to file
		lpr_phenotype.to_pickle(output["lpr_phecode"])


rule fix_wrongly_assigned_perinatal_and_pregnancy_diagnosis:
	input:
		lpr = rules.phecode_map_to_lpr.output.lpr_phecode,
		t_person = rules.all.input.t_person,
	output:
		lpr = temp("data/processed/lpr_and_phecodes_diagFixed.pkl"),
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,
	run:
		# Testing 
		#lpr = pd.read_pickle("data/processed/lpr_and_phecodes.pkl")
		#t_person = pd.read_csv("/users/secureome/home/projects/registries/lpr2bth/2018/cpr/t_person.tsv", dtype="str",sep="\t",usecols=["v_pnr_enc","C_KON","D_FODDATO","C_STATUS","D_STATUS_HEN_START","v_mor_pnr_enc"])
		# Snakemake
		lpr = pd.read_pickle(input.lpr)
		t_person = pd.read_csv(input.t_person, dtype="str",sep="\t",usecols=["v_pnr_enc","C_KON","D_FODDATO","C_STATUS","D_STATUS_HEN_START","v_mor_pnr_enc"])
		t_person["D_FODDATO"] = pd.to_datetime(t_person["D_FODDATO"],format="%Y-%m-%d",errors="coerce")
		t_person = t_person[t_person.v_pnr_enc != "MISSING_IN_BTH"].copy()
		t_person.drop_duplicates(inplace=True)

		# --- Correct diagnosis mistakenly given to newborn when should be assigned to mother and vice versa --- #

		# Subset lpr to does of O and P diagnosis
		lpr_subset = lpr.loc[(lpr.Diagnosis.str[0] == "O") | 
			(lpr.Diagnosis.str[0] == "P")].copy()

		lpr_subset = pd.merge(lpr_subset,t_person[["v_pnr_enc","D_FODDATO","v_mor_pnr_enc"]],left_on="cpr_enc",right_on = "v_pnr_enc",how="left")
		lpr_subset.drop("v_pnr_enc",axis=1,inplace=True)

		lpr_subset["D_age"] = lpr_subset["D_INDDTO"] - lpr_subset["D_FODDATO"]
		lpr_subset["D_age"] = lpr_subset["D_age"].dt.days


		# Contruct LPR for mother: N = 105
		mask_wrongly_assigned_newborn = (lpr_subset.Diagnosis.str[0] == "O") & (lpr_subset.D_age < 365.25*10)
		lpr_mother = lpr_subset.loc[mask_wrongly_assigned_newborn].copy()
		lpr_mother.drop("cpr_enc",axis=1,inplace=True)
		lpr_mother.rename({"v_mor_pnr_enc":"cpr_enc"},axis=1,inplace=True)
		lpr_mother = lpr_mother[list(lpr.columns)].copy()
		# Remove from full LPR
		lpr.drop(index=lpr_mother.index,inplace=True)

		# Wrongly assigned to mother: More complicated as mother can have more children
		# Only take perinatal diagnosis! The others are to unsure and doesnt fit for many of them
		# N = 13.236
		mask_wrongly_assigned_mother = ((lpr_subset.Diagnosis.str[0] == "P") & (lpr_subset.D_age > 365.25*10))
		# subset
		lpr_newborn = lpr_subset.loc[mask_wrongly_assigned_mother].copy()
		lpr_newborn = lpr_newborn[list(lpr.columns)]
		lpr_newborn.rename({"cpr_enc":"v_mor_pnr_enc2"},axis=1,inplace=True)

		# Remove from full LPR
		lpr.drop(index=lpr_newborn.index,inplace=True)
		# Get CPR of newborn to whom the diagnosis belong and add t_person to newborn with diagnosis
		lpr = pd.merge(lpr,t_person[["v_pnr_enc","D_FODDATO","v_mor_pnr_enc"]],left_on="cpr_enc",right_on = "v_pnr_enc",how="left")
		newborn_info = lpr[["cpr_enc","v_mor_pnr_enc","D_FODDATO"]].drop_duplicates().copy()
		lpr.drop(["v_pnr_enc","D_FODDATO","v_mor_pnr_enc"],axis=1,inplace=True)

		# Map diagnosis to newborn if it happend close to foddato
		# SQL to memory
		conn = sqlite3.connect(':memory:')

		#write the tables
		lpr_newborn.to_sql('lpr_newborn', conn, index=False,if_exists="replace")
		# Contruct mapping limits for diagnosis of the mother
		# 10 weeks after delivery
		newborn_info["D_FODDATO_last"] = newborn_info["D_FODDATO"] + pd.Timedelta("70 days")
		# 30 weeks before delivery
		newborn_info["D_FODDATO_first"] = newborn_info["D_FODDATO"] - pd.Timedelta("210 days")
		newborn_info.to_sql('newborn_info', conn, index=False,if_exists="replace")

		# Map mothers hemolytic reaction on foster to 30 weeks before and 4 weeks after delivery
		qry = '''
		SELECT * FROM newborn_info
		INNER JOIN lpr_newborn ON newborn_info.v_mor_pnr_enc = lpr_newborn.v_mor_pnr_enc2
		WHERE lpr_newborn.D_INDDTO >= newborn_info.D_FODDATO_first AND
		lpr_newborn.D_INDDTO <= newborn_info.D_FODDATO_last
		'''
		# Run query
		lpr_newborn = pd.read_sql_query(qry, conn)
		lpr_newborn.drop(["v_mor_pnr_enc","v_mor_pnr_enc2","D_FODDATO","D_FODDATO_first","D_FODDATO_last"],axis=1,inplace=True)
		cursor = conn.cursor()
		cursor.execute("DROP TABLE lpr_newborn")
		cursor.execute("DROP TABLE newborn_info")
		conn.close()

		# concatenate 2 reversed dataframes to add
		lpr_to_add = pd.concat([lpr_mother,lpr_newborn],ignore_index=True,axis=0)
		del t_person, lpr_mother, lpr_newborn

		# Add the right assigned diagnosis
		lpr = pd.concat([lpr,lpr_to_add], ignore_index=True, axis=0)
		lpr["D_INDDTO"] = pd.to_datetime(lpr["D_INDDTO"])
		lpr["D_UDDTO"] = pd.to_datetime(lpr["D_UDDTO"])
		del lpr_to_add

		lpr.to_pickle(output.lpr)


rule add_blocklevel_and_level3_label_to_lpr:
	input:
		icd10_block_codes = rules.all.input.icd10_block_codes,
		lpr_phecode = rules.fix_wrongly_assigned_perinatal_and_pregnancy_diagnosis.output.lpr,
	output:
		lpr_block = temp("data/processed/lpr_phecodes_block.pkl"),
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,

	run:
		''' testing
		lpr_phecode = pd.read_pickle("data/processed/lpr_ICD8_to_10.pkl")
		block_codes = pd.read_csv("/users/secureome/home/classification/complete/icd10_eng_diag_chapters_all_update2016.tsv",sep="\t",header=None,dtype=str)
		block_codes.rename({0:"Diagnosis_map",3:"ICDchapter",4:"block_name",5:"block_code"},axis=1,inplace=True)
		merge = pd.merge(lpr_phecode,block_codes.loc[:,["Diagnosis","block_name","block_code"]],on="Diagnosis",how="left")
		merge.to_pickle(output.lpr_block)
		'''
		lpr_phecode = pd.read_pickle(input.lpr_phecode)
		block_codes = pd.read_csv(input.icd10_block_codes,sep="\t",header=None,dtype=str)
		block_codes.rename({0:"Diagnosis_map",3:"ICDchapter",4:"block_name",5:"block_code"},axis=1,inplace=True)

		# Cut diagnosis code to max 3 length to match block level 3. 
		# As level this is enough to get block codes and chapterlevels
		lpr_phecode["Diagnosis_map"] = lpr_phecode.Diagnosis.str[:3]
		block_codes["Diagnosis_map"] = block_codes.Diagnosis_map.str[:3]
		block_codes.drop_duplicates("Diagnosis_map",inplace=True)

		# Merge and save 
		merge = pd.merge(lpr_phecode,block_codes.loc[:,["Diagnosis_map","block_name","block_code","ICDchapter"]],on="Diagnosis_map",how="left")

		merge.drop(["Diagnosis_map"],axis=1,inplace=True) # Keep full diagnosis code
		merge.to_pickle(output.lpr_block)


rule clean_bloodtype:
	input:
		bloodtypes = rules.all.input.bloodtypes,
		t_person = rules.all.input.t_person,
	output:
		outfile = "data/interim/bloodtypes_cleaned.pkl",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,
	run:
		#bloodtypes = pd.read_csv("/users/secureome/home/projects/bth/personal_folders/vicmuse/BloodType/bloodtypes_cleaned.tsv",sep="\t",names=["cpr_enc","Date_bt_meassured","Bloodtype","AB0","Rhesus"])
		bloodtypes = pd.read_csv(input["bloodtypes"],sep="\t",names=["cpr_enc","Date_bt_meassured","Bloodtype","AB0","Rhesus"])
		
		# First drop bloodtypes meassured out of BTH period
		start_of_BTH = pd.Timestamp(datetime.date(2006, 1, 1))
		end_of_BTH = pd.Timestamp(datetime.date(2018, 4, 10))
		bloodtypes["Date_bt_meassured"] = pd.to_datetime(bloodtypes["Date_bt_meassured"],format="%d-%m-%Y",errors="coerce")
		bloodtypes = bloodtypes[(bloodtypes.Date_bt_meassured>=start_of_BTH)&(bloodtypes.Date_bt_meassured<=end_of_BTH)].copy()

		# Remove cprs with changing bloodtypes
		# First drop duplicates with same meassured bloodtype
		bloodtypes.drop_duplicates(["cpr_enc","Bloodtype"],inplace=True,keep="first")
		# Then drop CPRs with changing bloodtypes
		bloodtypes.drop_duplicates(["cpr_enc"],inplace=True,keep=False)

		# Remove wierd bloodtypes
		mask = (bloodtypes.Bloodtype=="AB+DU")|(bloodtypes.Bloodtype=="0+DU")|(bloodtypes.Bloodtype=="B+DU")|(bloodtypes.Bloodtype=="A+DU")|(bloodtypes.Bloodtype=="XX")
		bloodtypes = bloodtypes.loc[~mask,:].copy()

		# Add birthdate and sex
		t_person = pd.read_csv(input["t_person"], dtype="str",sep="\t",usecols=["v_pnr_enc","C_KON","D_FODDATO","C_STATUS","D_STATUS_HEN_START"])
		#t_person = pd.read_csv("/users/secureome/home/projects/registries/lpr2bth/2018/cpr/t_person.tsv", dtype="str",sep="\t",usecols=["v_pnr_enc","C_KON","D_FODDATO","C_STATUS","D_STATUS_HEN_START"])
		
		# Only patient in BTH
		t_person = t_person[t_person.v_pnr_enc != "MISSING_IN_BTH"].copy()
		t_person.rename({"v_pnr_enc":"cpr_enc"},axis=1,inplace=True)
		t_person.drop_duplicates(inplace=True) # none duplicates!

		t_person["D_FODDATO"] = pd.to_datetime(t_person["D_FODDATO"],format="%Y-%m-%d",errors="coerce")
		t_person["D_STATUS_HEN_START"] = pd.to_datetime(t_person["D_STATUS_HEN_START"],format="%Y-%m-%d",errors="coerce")
		t_person["Birth_year"] = t_person["D_FODDATO"].dt.year

		# Drop patient dying or moving before 1977-1-1
		start_of_registry = pd.Timestamp(datetime.date(1977, 1, 1))
		t_person = t_person.loc[~((t_person.D_STATUS_HEN_START < start_of_registry)&(t_person.C_STATUS != 1))].copy()

		# Add death
		t_person["Dead"] = 0
		t_person.loc[t_person.C_STATUS == "90", "Dead"] = 1
		# Merge onto bloodtypes
		merge = pd.merge(bloodtypes,t_person,on="cpr_enc",how="inner")
		# save
		merge.to_pickle(output.outfile)


rule add_bloodtype_to_lpr:
	input:
		lpr = rules.add_blocklevel_and_level3_label_to_lpr.output.lpr_block,
		bloodtypes = rules.clean_bloodtype.output.outfile,
	output:
		lpr = temp("data/processed/bloodtype_and_lpr.pkl"),
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,
	run:
		# testing
		#lpr = pd.read_pickle("data/processed/lpr_phecodes_block.pkl")
		#bloodtypes = pd.read_pickle("data/interim/bloodtypes_cleaned.pkl")

		# Load lpr and bloodtypes
		lpr = pd.read_pickle(input.lpr)
		bloodtypes = pd.read_pickle(input.bloodtypes)

		# Include only A or B diagnosis
		lpr = lpr.loc[(lpr.DiagnosisType == "A")|(lpr.DiagnosisType == "B")]
		# Drop diagnosis before start of registry and end of registry
		start_of_registry = pd.Timestamp(datetime.date(1977, 1, 1))
		end_of_registry = pd.Timestamp(datetime.date(2018, 4, 10))
		lpr = lpr.loc[(lpr.D_INDDTO >= start_of_registry)&(lpr.D_INDDTO <= end_of_registry),:]

		# Merge bloodtypes and lpr
		merge = pd.merge(lpr,bloodtypes,on="cpr_enc",how="inner")
		del lpr, bloodtypes

		# Drop diagnosis given before date of birth
		# Non are dropped!
		merge = merge.loc[merge.D_INDDTO >= merge.D_FODDATO].copy()


		# Entry data
		merge["D_entry"] = merge["D_FODDATO"]
		merge.loc[merge.D_entry < start_of_registry,"D_entry"] = start_of_registry

		# Get days to death
		merge["Death_date"] = pd.NaT
		mask_dead = (merge.C_STATUS == "90")
		merge.loc[mask_dead, "Death_date"] = merge.loc[mask_dead, "D_STATUS_HEN_START"] 
		merge["Days_to_death"] = (merge.Death_date - merge.D_entry).dt.days

		# Get days to move
		merge["Move_date"] = pd.NaT
		mask_move = ((merge.C_STATUS != "90")&(merge.C_STATUS != "01"))
		merge.loc[mask_move, "Move_date"] = merge.loc[mask_move, "D_STATUS_HEN_START"] 
		merge["Days_to_move"] = (merge.Move_date - merge.D_entry).dt.days

		# Calculate days from registry entry to end of registry
		merge["Days_to_end_of_registry"] = (end_of_registry - merge.D_entry).dt.days

		# Days to diagnosis
		merge["Days_to_diagnosis"] = (merge.D_INDDTO - merge.D_entry).dt.days

		merge["age_at_entry"] = (start_of_registry - merge["D_FODDATO"]).dt.days
		merge["age_at_entry"] = merge["age_at_entry"] / 365.25
		merge.loc[merge.D_FODDATO > start_of_registry,"age_at_entry"] = 0

		# Year of entry
		merge["Year_of_entry"] = merge["D_FODDATO"].dt.year
		merge.loc[merge.D_FODDATO.dt.year < start_of_registry.year,"Year_of_entry"] = start_of_registry.year

		merge.to_pickle(output.lpr)


rule list_of_cprs_used_for_analysis:
	input:
		bloodtypes = rules.add_bloodtype_to_lpr.output.lpr,
	output:
		cpr_list = "data/processed/instudy_cpr_list.csv",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,
	run:
		# Load
		df = pd.read_pickle(input.bloodtypes)
		#df = pd.read_pickle("data/processed/bloodtype_and_lpr.pkl")
		# Get unique cprs
		df.drop_duplicates(subset=["cpr_enc"],inplace=True)
		# Save list for PJ to use
		df.cpr_enc.to_csv(output.cpr_list,sep="\t",index=False)
		#df.cpr_enc.to_csv("data/processed/instudy_cpr_list.csv",sep="\t",index=False)



rule fix_length3_phecodes:
	input:
		infile = rules.add_bloodtype_to_lpr.output.lpr,
		phecode_map = rules.all.input.phecodes,
	output:
		outfile = "data/processed/bloodtype_final.pkl",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,
	run:
		#data = pd.read_pickle("data/processed/bloodtype_and_lpr.pkl")
		#phecode_map = pd.read_csv("/users/secureome/home/people/petras/base_data/Phecode_map_v1_2_icd10_beta.csv",sep=",",dtype="str")

		data = pd.read_pickle(input.infile)
		phecode_map = pd.read_csv(input.phecode_map,sep=",",dtype="str")

		# Select length 3 phecodes
		phecode_map = phecode_map[phecode_map.PHECODE.str.len() == 3].copy()
		# Make subset of LPR of length 3 PHECODE to add to LPR
		sub_length3 = data.copy()
		sub_length3["PHECODE"] = sub_length3.PHECODE.str[:3]
		# Only include the length 3 phecodes from the phecode mapping file
		sub_length3[sub_length3.PHECODE.isin(phecode_map.PHECODE)].copy()
		# Drop duplicates
		sub_length3.sort_values(by=["cpr_enc","D_INDDTO"],inplace=True)
		sub_length3.drop_duplicates(["cpr_enc","D_INDDTO","PHECODE"],inplace=True)
		# Concatenate new set
		concat = pd.concat([data,sub_length3],axis=0,ignore_index=True)
		# Drop duplicates
		concat.sort_values(by=["cpr_enc","D_INDDTO"],inplace=True)
		concat.drop_duplicates(["cpr_enc","D_INDDTO","PHECODE"],inplace=True)

		concat.to_pickle(output.outfile)



rule create_table1:
	input:
		table1_file = rules.fix_length3_phecodes.output.outfile,
	output:
		table1 = "results/table1.html",
		#plot_age_at_entry = "results/Age_entry_count_plot.png",
		#plot_birth_year_count = "results/Birth_year_count_plot.png",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,

	script:
		"scripts/table1.py"



rule filter_phecodes_and_blocklevels_level3_to_include:
	input:
		lpr = rules.fix_length3_phecodes.output.outfile,
		phecodes = "/users/secureome/home/people/petras/base_data/phecode_def.csv",
	output:
		phecodes_to_include = "data/processed/included_phecodes.tsv",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*50,
	threads: 1,
	run:
		# Phecode inclusion list
		# Make inclusion list
		#data = pd.read_pickle("data/processed/bloodtype_and_lpr.pkl")
		#data = pd.read_pickle("data/processed/bloodtype_final.pkl")

		data = pd.read_pickle(input.lpr)
		# Load phecode definitions
		#phecode_def = pd.read_csv("/users/secureome/home/people/petras/base_data/phecode_def.csv",sep=",",dtype="str")
		phecode_def = pd.read_csv(input.phecodes,sep=",",dtype="str")
		# Map descriptive name to phecode
		phecode_def.rename({"phecode":"PHECODE"},inplace=True,axis=1)
		data = pd.merge(data,phecode_def,on="PHECODE",how="left")

		# Do not include injuries & poisonings (category 18) and phecodes => 1000
		# Specific for phecodes also not codes not mapped to category
		data = data.loc[(data.category.notnull())&(data.category != "injuries & poisonings")&(data.category != "symptoms")]
		#(data.PHECODE.astype(float) >= 1000).value_counts() # none

		# Sort by diagnosis date and keep first phecode
		data.sort_values(by=["cpr_enc","D_INDDTO"],inplace=True)
		data.dropna(subset=["PHECODE"],axis=0,inplace=True)
		data.drop_duplicates(subset=["cpr_enc","PHECODE"],inplace=True,keep="first")

		# Calculate prevalance of phecodes in study population
		top_phecodes = data.PHECODE.value_counts()/data.cpr_enc.nunique()
		# Select phecodes with >= 100
		mask_include_phecodes_count = data.PHECODE.value_counts() >= 100
		included_phecodes = top_phecodes[mask_include_phecodes_count].index
		pd.Series(data=included_phecodes,name="Phecodes",dtype=str).to_csv(output.phecodes_to_include,sep="\t",index=False,header=["Phecodes"])


## --- Log-linear Poisson Regression --- ##

rule prepare_PheWAS_input:
	input:
		lpr = rules.fix_length3_phecodes.output.outfile,
	output:
		phewas_data = "data/processed/PheWAS/{phecode}.tsv",
	params:
		phecode = "{phecode}",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*10,
	threads: 1,
	run:
		#data = pd.read_pickle("data/processed/bloodtype_final.pkl")
		data = pd.read_pickle(input.lpr)

		# Sort on PHECODE and D_INDDTO and remove dups ( keep first occurence)
		data.sort_values(by=["cpr_enc","D_INDDTO"],inplace=True)
		data.drop_duplicates(["cpr_enc","PHECODE"],inplace=True)


		# -- Make CASES and CONTROLS --- #
		# Make study data for phecode
		phecode = str(params.phecode)

		# Try to make a subset for phecode
		data["CASE"] = 0
		# Follow-up time is the minimum of death, migration and end of registry
		data["TIME"] = data.loc[:,["Days_to_death", "Days_to_move", "Days_to_end_of_registry"]].min(axis=1)
		# Set CASE
		data.loc[data.PHECODE == phecode,"CASE"] = 1
		# Set time to days to diagnosis
		data.loc[data.PHECODE == phecode,"TIME"] = data.loc[data.PHECODE == phecode,"Days_to_diagnosis"]

		# Set age at exit
		data["age_at_exit"] = data["age_at_entry"] + data["TIME"]/365.25

		# Define cases and controls
		cases = data.loc[data.CASE == 1].copy()
		cases.cpr_enc.nunique() # All are unique as expected
		# Not cases
		controls = data.loc[~(data.cpr_enc.isin(cases.cpr_enc))].copy()
		del data

		# Sort so the entry of control with lowest TIME is used
		controls.sort_values(by=["cpr_enc","TIME"],inplace=True)
		controls.drop_duplicates(subset=["cpr_enc"],keep="first",inplace=True)

		study_data = pd.concat([cases,controls],ignore_index=True,axis=0)
		del cases, controls

		# Drop dups -> none were present
		study_data.sort_values(by=["cpr_enc","TIME"],inplace=True)
		study_data.drop_duplicates(subset=["cpr_enc"],keep="first",inplace=True)

		study_data[["age_at_entry","age_at_exit","CASE","TIME","C_KON","AB0","Rhesus","Birth_year"]].to_csv(output.phewas_data,sep="\t",index=False,na_rep="NaN")
		del study_data	


## --- Run analysis --- ##

# 1. Analyses 

rule poisson_allVSOne_ABO:
	input:
		phewas_data = "data/processed/PheWAS/{phecode}.tsv",
	output:
		estimates = "results/20220531/allVSOne/enter_registry/phewas/{blood_group}/estimates_{phecode}.tsv",
	params:
		blood_group = "{blood_group}",
	resources:
		tmin = 60*2,
		mem_mb = 1024*6,
	threads: 1,
	run:
		# ---  RUN LOG-LINEAR POISSON REGRESSION --- #
		Poisson_CMB = [
			'Rscript scripts/poisson_allVSOne_ABO.R',
			'{input.phewas_data}',
			'{output.estimates}',
			'{params.blood_group}']
		shell(' '.join(Poisson_CMB))


rule poisson_allVSOne_RhD:
	input:
		phewas_data = "data/processed/PheWAS/{phecode}.tsv",
	output:
		estimates = "results/20220531/allVSOne/enter_registry/phewas/RhD/RhD_estimates_{phecode}.tsv",
	resources:
		tmin = 60*2,
		mem_mb = 1024*6,
	threads: 1,
	run:
		# ---  RUN LOG-LINEAR POISSON REGRESSION --- #
		Poisson_CMB = [
			'Rscript scripts/poisson_allVSOne_RhD.R',
			'{input.phewas_data}',
			'{output.estimates}']
		shell(' '.join(Poisson_CMB))


rule collect_allVSOne_phewas:
	input:
		AB0_phewas_estimates = expand("results/20220531/allVSOne/enter_registry/phewas/{blood_group}/estimates_{phecode}.tsv",phecode=phecodes_to_include,blood_group=["A","B","AB","0"]),
		RhD_phewas_estimates = expand("results/20220531/allVSOne/enter_registry/phewas/RhD/RhD_estimates_{phecode}.tsv",phecode=phecodes_to_include),
		phecodes = "/users/people/petras/base_data/phecode_def.csv",
	output:
		phewas = "results/20220531/allVSOne/enter_registry/phewas_estimates.tsv",
	resources:
		tmin = 60*2,
		mem_mb = 1024*30,
	threads: 1,
	script:
		"scripts/collect_results_allVSOne.py"



## --- Age of first diagnosis in hospital --- ##

rule prepare_phewas_age_at_diag:
	input:
		lpr = rules.fix_length3_phecodes.output.outfile,
	output:
		#phewas_data = temp("data/processed/age_at_diagnosis/phewas/{phecode}.tsv"),
		phewas_data = "data/processed/age_at_diagnosis/phewas/{phecode}.tsv",
	params:
		phecode = "{phecode}",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*20,
	threads: 1,
	run:
		data = pd.read_pickle(input.lpr)
		data = data[["cpr_enc","age_at_entry","D_INDDTO","D_FODDATO","PHECODE","Exl. Phecodes","Days_to_end_of_registry","Days_to_diagnosis","Days_to_move","Days_to_death","C_KON","AB0","Rhesus","Birth_year"]]

		# Find cases
		phecode = str(params.phecode)

		# Remove patient records without the phecode
		data = data.loc[data.PHECODE == phecode,:].copy()

		# Sort by diagnosis date and keep first phecode
		data.sort_values(by=["cpr_enc","D_INDDTO"],inplace=True)
		data.dropna(subset=["PHECODE"],axis=0,inplace=True)
		data.drop_duplicates(subset=["cpr_enc","PHECODE"],inplace=True,keep="first")

		# Estimate age at diagnosis for phecode
		data["age_at_diagnosis"] = 	(data.loc[data.PHECODE == phecode,"D_INDDTO"] - data.loc[data.PHECODE == phecode,"D_FODDATO"]).dt.days/365.25
		data["age_at_diagnosis"] = np.round(data["age_at_diagnosis"],2)
		# Remove patient getting the diagnosis on at the start of the registry, as we dont know when they actually got it
		start_of_registry = pd.Timestamp(datetime.date(1977, 1, 1))
		data = data.loc[data.D_INDDTO > start_of_registry,:].copy()

		# Save data to be used for linear regression
		data[["age_at_diagnosis","C_KON","AB0","Rhesus","Birth_year"]].to_csv(output.phewas_data,sep="\t",index=False,na_rep="NaN")


# Age of diagnosis

rule RUN_phewas_age_at_diag_allVSOne_ABO:
	input:
		phewas_data = "data/processed/age_at_diagnosis/phewas/{phecode}.tsv",
	output:
		estimates = "results/20220531/allVSOne/age_at_diagnosis/phewas/{blood_group}/estimates_{phecode}.tsv",
	params:
		blood_group = "{blood_group}",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*5,
	threads: 1,
	run:
		LinReg_CMB = [
			'Rscript scripts/linreg_allVSOne_ABO.R',
			'{input.phewas_data}',
			'{output.estimates}',
			'{params.blood_group}']
		shell(' '.join(LinReg_CMB))

rule RUN_phewas_age_at_diag_allVSOne_RhD:
	input:
		phewas_data = "data/processed/age_at_diagnosis/phewas/{phecode}.tsv",
	output:
		estimates = "results/20220531/allVSOne/age_at_diagnosis/phewas/RhD/RhD_estimates_{phecode}.tsv",
	resources:
		tmin = 60*24*1,
		mem_mb = 1024*5,
	threads: 1,
	run:
		LinReg_CMB = [
			'Rscript scripts/linreg_allVSOne_RhD.R',
			'{input.phewas_data}',
			'{output.estimates}']
		shell(' '.join(LinReg_CMB))


rule collect_allVSOne_phewas_linreg:
	input:
		AB0_phewas_estimates = expand("results/20220531/allVSOne/age_at_diagnosis/phewas/{blood_group}/estimates_{phecode}.tsv",phecode=non_congenital_phecodes_to_include,blood_group=["A","B","AB","0"]),
		RhD_phewas_estimates = expand("results/20220531/allVSOne/age_at_diagnosis/phewas/RhD/RhD_estimates_{phecode}.tsv",phecode=non_congenital_phecodes_to_include),
		phecodes = "/users/people/petras/base_data/phecode_def.csv",
	output:
		phewas = "results/20220531/allVSOne/age_at_diagnosis/phewas_estimates.tsv",
	resources:
		tmin = 60*2,
		mem_mb = 1024*30,
	threads: 1,
	script:
		"scripts/collect_results_allVSOne.py"





