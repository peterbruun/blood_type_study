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

# FDR adjusting for the TOTAL number of tests (4 blood groups * 1100 phecodes)
all_estimates["p.value_FDR"] = np.nan
all_estimates["FDR_.05"] = np.nan
all_estimates.loc[:,"p.value_FDR"] = fdrcorrection(all_estimates.loc[:,"p.value"])[1]
all_estimates.loc[:,"FDR_.05"] = fdrcorrection(all_estimates.loc[:,"p.value"])[0]


# Map descriptive name to phecode
phecodes = pd.read_csv(snakemake.input["phecodes"],sep=",",dtype="str")
all_estimates = pd.merge(all_estimates,phecodes,on="phecode",how="left")

# Save
all_estimates.to_csv(snakemake.output["phewas"],sep="\t",na_rep="NA",index=False)

