# collect analysis results of all diagnosis tested

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


# -- PHEWAS -- #

# Collect ABO analyses
# Get list of all phewas results
estimates_AB0 = snakemake.input["AB0_phewas_estimates"]
estimates_RhD = snakemake.input["RhD_phewas_estimates"]

estimates_list = estimates_AB0 + estimates_RhD

# Testing
#import os
#os.chdir("/users/projects/bloodtype/")
#phecodes_to_include = pd.read_csv("data/processed/included_phecodes.tsv",sep="\t",dtype=str)
#phecodes_to_include = list(phecodes_to_include.Phecodes.values)
#estimates_list = ["results/20220531/allVSOne/enter_registry/phewas/A/estimates_"+phecode+".tsv" for phecode in phecodes_to_include] + \
#				 ["results/20220531/allVSOne/enter_registry/phewas/B/estimates_"+phecode+".tsv" for phecode in phecodes_to_include] + \
#				 ["results/20220531/allVSOne/enter_registry/phewas/AB/estimates_"+phecode+".tsv" for phecode in phecodes_to_include] + \
#				 ["results/20220531/allVSOne/enter_registry/phewas/0/estimates_"+phecode+".tsv" for phecode in phecodes_to_include] + \
#				 ["results/20220531/allVSOne/enter_registry/phewas/RhD/RhD_estimates_"+phecode+".tsv" for phecode in phecodes_to_include]
#diagnosis_splitter = estimates_list[0].split("estimates_")[0] + "estimates_"

first_iter = True
for results in estimates_list:
	estimate = pd.read_csv(results,sep="\t")
	estimate = estimate.loc[estimate.term == "factor(exposure)1",:]
	blood_group = results.split("/")[5]
	estimate["blood_group"] = blood_group
	if blood_group == "RhD":
		estimate["phecode"] = results.split("/")[6].strip(".tsv").strip("RhD_estimates_")
	else:
		estimate["phecode"] = results.split("/")[6].strip(".tsv").strip("estimates_")
	if first_iter:
		all_estimates = estimate.copy()
		first_iter = False
	else:
		all_estimates = pd.concat([all_estimates,estimate],axis=0,ignore_index=True)


# FDR adjusting for the TOTAL number of tests (5 blood groups * 1100 phecodes)
all_estimates["p.value_FDR"] = np.nan
all_estimates["FDR_.05"] = np.nan
all_estimates.loc[:,"p.value_FDR"] = fdrcorrection(all_estimates.loc[:,"p.value"])[1]
all_estimates.loc[:,"FDR_.05"] = fdrcorrection(all_estimates.loc[:,"p.value"])[0]


# Map descriptive name to phecode
phecodes = pd.read_csv(snakemake.input["phecodes"],sep=",",dtype="str")
all_estimates = pd.merge(all_estimates,phecodes,on="phecode",how="left")

# Save
all_estimates.to_csv(snakemake.output["phewas"],sep="\t",na_rep="NA",index=False)




