

import numpy as np
import pandas as pd

from collections import defaultdict
from datetime import datetime

####################################################################################
## t_adm
# https://www.esundhed.dk/Dokumentation?rid=5&tid=7&vid=62

#d_inddto:
#	Variablebeskrivelse
# Kontaktens startdato defineret som indlaeggelsesdato for indlagte patienter, 
# dato for start paa deldoegnsbehandling, dato for start paa ambulant kontakt 
# for ambulante patienter og dato for ankomst for skadestuepatienter.

# d_uddto
#	Variablebeskrivelse
#Kontaktens afslutningsdato defineret som udskrivningsdato for indlagte patienter,
# dato for afslutning af deldoegnsbehandling og dato 
# for afslutning af ambulant forloeb for ambulante patienter. 
# Saafremt afslutningsdatoen for en skadestuekontakt ikke er indberettet, 
# saettes afslutningsdatoen lig med kontaktens startdato, 
# saaledes at skadestuekontakter altid vil have en afslutningsdato i registeret.
####################################################################################
	
#bloodflodet = "data/bloodflodet_cleaned.pkl"
#master_file = "/data/projects/deep_phenotyping/lpr/transfusion_master_table_enc_dbds.tsv"
#t_adm_file = "/data/projects/deep_phenotyping/lpr/transfusion_t_adm_enc_dbds.tsv"
#t_diag_file = "/data/projects/deep_phenotyping/lpr/transfusion_t_diag.asc"

def parse_hospital_registry(adm_file=snakemake.input["t_adm"], diag_file=snakemake.input["t_diag"], source='LPR'):
	# Modification of davids script where D_INDDTO is also included
	# https://github.com/DBDSDK/registry_data/blob/master/master_table.py
	source = "LPR"
	def prefix_icd(icd_code):
		if len(icd_code) == '':
			return icd_code
		elif icd_code.startswith('D'):
			return 'ICD10:%s' % icd_code
		else:
			return 'ICD8:%s' % icd_code

	COL_ORDER = ['cpr_enc', 'Source', 'D_INDDTO','D_UDDTO', 'DiagnosisType', 'Diagnosis']
	
	#adm_file = "/home/projects/registries/2018_newest/prepared_data/t_adm.tsv"
	#adm_file = "/home/projects/registries/lpr2bth/2018/classic_style_lpr/preprocessing_death_column/prepared_data/t_adm.tsv"
	t_adm = pd.read_csv(adm_file, sep='\t', dtype=str, usecols=['v_cpr_enc', 'k_recnum','d_inddto','d_uddto', 'c_adiag'],low_memory=False)
	t_adm.rename(columns={"c_adiag": "Diagnosis","v_cpr_enc":"cpr_enc","k_recnum":"K_RECNUM","d_inddto":"D_INDDTO","d_uddto":"D_UDDTO"}, inplace=True)
	t_adm.loc[t_adm.Diagnosis.isnull(), "Diagnosis"] = ''
	t_adm.Diagnosis = t_adm.Diagnosis.apply(prefix_icd)
	t_adm['DiagnosisType'] = 'A'
	t_adm['Source'] = source

	#combined_registry = [t_adm[COL_ORDER]]
	#del t_adm

	#diag_file = "/home/projects/registries/lpr2bth/2018/classic_style_lpr/preprocessing_death_column/prepared_data/t_diag.tsv"
	t_diag = pd.read_csv(diag_file, sep='\t', dtype=str, usecols=['V_RECNUM', 'C_DIAG', 'C_DIAGTYPE'],low_memory=False)
	t_diag.rename(columns={'C_DIAG': 'Diagnosis', 'C_DIAGTYPE': 'DiagnosisType'}, inplace=True)

	t_diag.Diagnosis = t_diag.Diagnosis.apply(prefix_icd)
	t_diag = t_diag.merge(t_adm[['cpr_enc', 'K_RECNUM', 'D_INDDTO','D_UDDTO']], left_on='V_RECNUM', right_on='K_RECNUM', how='left')
	t_diag['Source'] = source
	#combined_registry.append(t_diag[COL_ORDER])
	#del t_diag

	#x = pd.concat([t_adm[COL_ORDER],t_diag[COL_ORDER]]).sort_values(['cpr_enc', 'D_INDDTO']).drop_duplicates()

	return pd.concat([t_adm[COL_ORDER] , t_diag[COL_ORDER]]).sort_values(['cpr_enc' , 'D_INDDTO']).drop_duplicates()


# main
# Merge LPR
master = parse_hospital_registry()
master = master[master.cpr_enc.notnull()].copy()

# Remove whitesapce
master["DiagnosisType"] = master.DiagnosisType.str.strip()
# Remove H diagnosis
master = master[master.DiagnosisType != "H"].copy()
# convert to datetime
master["D_INDDTO"] = pd.to_datetime(master["D_INDDTO"],format="%Y-%m-%d",errors="coerce")
master["D_UDDTO"] = pd.to_datetime(master["D_UDDTO"],format="%Y-%m-%d",errors="coerce")
# Save LPR
master.to_pickle(snakemake.output["lpr"])



