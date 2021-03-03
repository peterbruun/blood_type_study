# collect analysis results of all diagnosis tested

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


# -- PHEWAS -- #
# Get list of all phewas results
estimates_list = snakemake.input["phewas_estimates"]
diagnosis_splitter = estimates_list[0].split("estimates_")[0] + "estimates_"

first_iter = True
for results in estimates_list:
	estimate = pd.read_csv(results,sep="\t")
	estimate["phecode"] = results.split(diagnosis_splitter)[1].strip(".tsv")
	if first_iter:
		all_estimates = estimate.copy()
		first_iter = False
	else:
		all_estimates = pd.concat([all_estimates,estimate],axis=0,ignore_index=True)

# FDR correction
all_estimates["p.value_FDR"] = np.nan
all_estimates["FDR_.05"] = np.nan
for term in all_estimates.term.unique():
	mask = (all_estimates.term==term)
	all_estimates.loc[mask,"p.value_FDR"] = fdrcorrection(all_estimates.loc[mask,"p.value"])[1]
	all_estimates.loc[mask,"FDR_.05"] = fdrcorrection(all_estimates.loc[mask,"p.value"])[0]


# Map descriptive name to phecode
phecodes = pd.read_csv(snakemake.input["phecodes"],sep=",",dtype="str")
all_estimates = pd.merge(all_estimates,phecodes,on="phecode",how="left")

# Save
all_estimates.to_csv(snakemake.output["phewas"],sep="\t",na_rep="NA",index=False)


# -- ICD LEVEL 3 -- #
# Get list of all phewas results
#estimates_list = [f for f in glob.glob("results/20210128/poisson/enter_registry/ICDlevel3/estimates_*")]
estimates_list = snakemake.input["level3_estimates"]
diagnosis_splitter = estimates_list[0].split("estimates_")[0] + "estimates_"


first_iter = True
for results in estimates_list:
	estimate = pd.read_csv(results,sep="\t")
	estimate["ICDlevel3"] = results.split(diagnosis_splitter)[1].strip(".tsv")
	if first_iter:
		all_estimates = estimate.copy()
		first_iter = False
	else:
		all_estimates = pd.concat([all_estimates,estimate],axis=0,ignore_index=True)

# FDR correction
all_estimates["p.value_FDR"] = np.nan
all_estimates["FDR_.05"] = np.nan
for term in all_estimates.term.unique():
	mask = (all_estimates.term==term)
	all_estimates.loc[mask,"p.value_FDR"] = fdrcorrection(all_estimates.loc[mask,"p.value"])[1]
	all_estimates.loc[mask,"FDR_.05"] = fdrcorrection(all_estimates.loc[mask,"p.value"])[0]


# Map descriptive name to level3
icdlevel3 = pd.read_csv(snakemake.input["icd10_block_codes"],sep="\t",dtype="str",header=None)
icdlevel3.rename({0:"Diagnosis",1:"Level3_def"},axis=1,inplace=True)
# Select only descriptive text from level 3
icdlevel3 = icdlevel3[icdlevel3.Diagnosis.str.len() == 3]
icdlevel3["ICDlevel3"] = icdlevel3.Diagnosis
# Merge info
all_estimates = pd.merge(all_estimates,icdlevel3[["ICDlevel3","Level3_def"]],on="ICDlevel3",how="left")

# Save
all_estimates.to_csv(snakemake.output["ICDlevel3"],sep="\t",na_rep="NA",index=False)



