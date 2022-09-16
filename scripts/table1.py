
# import numerical libraries
import pandas as pd
import numpy as np
from scipy import stats
import os, datetime

pd.set_option('max_info_columns', 10**10)
pd.set_option('display.max_columns', 1000)
pd.set_option('max_info_columns', 10**10)
pd.set_option('max_info_rows', 10**10)


# Plotting
os.environ['QT_QPA_PLATFORM']='offscreen'
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
plt.switch_backend('agg')
from pylab import savefig


# import tableone
# https://tableone.readthedocs.io/en/latest/bestpractice.html
from tableone import TableOne


## --------------------------------------------------------- ##

# Testing

os.chdir("/users/projects/bloodtype/")

# Load file
data = pd.read_pickle(snakemake.input["table1_file"])
#data = pd.read_pickle("data/processed/bloodtype_final.pkl")

data.drop_duplicates(subset=["cpr_enc"],inplace=True)

# Table 1 with stratification
data.Birth_year = data.Birth_year.astype("float")
bins = pd.IntervalIndex.from_tuples([(1900, 1905), (1905, 1910), (1910, 1915),(1915, 1920),(1920, 1925),
	(1925, 1930),(1930, 1935),(1935, 1940),(1940, 1945),(1945, 1950),(1950, 1955),(1955, 1960),(1960, 1965),
	(1965, 1970),(1970, 1975),(1975, 1980),(1980, 1985),(1985, 1990),(1990, 1995),(1995, 2000),(2000, 2005),(2005, 2010),(2010, 2015),(2015, 2020)])

data["Birth_year_bin"] = pd.cut(data.Birth_year,bins=bins)
data[data["Birth_year_bin"].isnull()]
data["Birth_year_bin"] = data["Birth_year_bin"].astype("str")

# Get follow up time without looking at PheCodes
data["Follow_up_time"] = data.loc[:,["Days_to_death", "Days_to_move", "Days_to_end_of_registry"]].min(axis=1)
data["Follow_up_time"] = data["Follow_up_time"]/365.25

# columns to summarize: this will be the order of the table
columns = ["AB0","Rhesus","C_KON","Birth_year","age_at_entry","Follow_up_time"]

# columns containing categorical variables
categorical = ["AB0","Rhesus","C_KON"]

# non-normal variables
#nonnormal continuous variables are summarized by 'median [Q1,Q3]' instead of mean (SD).
nonnormal = ["Birth_year",'age_at_entry',"Follow_up_time"]

# set decimal places for age to 0
decimals = {"age_at_entry": 0}

# optionally, a categorical variable for stratification
groupby = ["AB0","Rhesus"]



table1 = TableOne(data, columns=columns, categorical=categorical,
			nonnormal=nonnormal, label_suffix=True,
			decimals=decimals, smd=False,
			htest_name=True)

print(table1.tabulate(tablefmt = "fancy_grid"))

# Save to html
table1.to_html(snakemake.output["table1"])


## --------------------------------------------------------- ##
# Birth year distribution

plt.clf()
data.Birth_year = data.Birth_year.astype(int)
fig, ax = plt.subplots(figsize=(10, 7))
sns.set_style("whitegrid")
ax = sns.histplot(x="Birth_year", data=data, discrete =True)
ax.set(xlabel='Birth year', ylabel='Number of individuals')
#plt.xticks(rotation=0,size=10)
plt.tight_layout()
plt.savefig("results/Birth_year_count_plot.png",dpi=199)



# # Age at entry distribution
# plt.clf()
# data.age_at_entry = data.age_at_entry.astype(int)
# fig, ax = plt.subplots(figsize=(12, 8))
# #ax = sns.countplot(x="age_at_entry", data=data)
# ax = sns.distplot(data["age_at_entry"].round(0),kde=False)
# ax.set(xlabel="Age at entry in study", ylabel='Number of individuals')
# plt.xticks(rotation=0,size=10)
# #ax.set_xticklabels(labels, fontsize=14, rotation=30, ha= 'right')
# plt.tight_layout()
# plt.savefig("results/Age_entry_count_plot.png",dpi=199)



