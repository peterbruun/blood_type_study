import numpy as np
import pandas as pd
import re, os, pickle, datetime, sys
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


def time_to_event(df,days,target,duration_col):
	
	df[]

	days = int(days)
	# testing
	df[duration_col] = df[duration_col].astype("float")
	# Construct column
	df[target] = 0

	# Assign 1 if death occured for patient within N days
	df.loc[df[duration_col] <= days,target] = 1

	#df.loc[df[duration_col] > df["Days_to_next_trans"],target] = 0

	# Set people who doesnt die within N days to N days
	#df.loc[(df[duration_col]>days)|(df[duration_col].isnull()),duration_col] = days
	df.loc[(df[target]==0)|(df[duration_col].isnull()),duration_col] = days

	#mask_tran_before_death_or_followup = ((df["Days_to_next_trans"]<df[duration_col])&(df["Days_to_next_trans"]<days))
	#df.loc[mask_tran_before_death_or_followup,duration_col] = df.loc[mask_tran_before_death_or_followup,"Days_to_next_trans"]

	# Test if working
	#df[[duration_col,target,"Days_to_next_trans","death_org"]].head(20)
	#df.loc[df.death_org<28,[duration_col,target,"Days_to_next_trans","death_org"]].head(20)

	# End of data censoring
	# End date
	end_of_data = pd.Timestamp(datetime.date(2018, 4, 10))
	# Calculate days from transfusion to end of data
	df["Days_to_end_of_data"] = (end_of_data - df.transfusion_datetime).dt.days
	# If duration column is bigger than end of data
	mask_enddata_before_censoring = df["Days_to_end_of_data"] < df[duration_col]
	# Set potentiel cases to 0
	df.loc[mask_enddata_before_censoring,target] = 0
	# Set duration col to days to end of data
	df.loc[mask_enddata_before_censoring,duration_col] = df.loc[mask_enddata_before_censoring,"Days_to_end_of_data"]

	# Round if Days_to_next_trans to nearest whole day
	df[duration_col] = df[duration_col].round()
	# Round days to death
	df[duration_col] = df[duration_col].astype(int)

	# Old
	#mask_succeed_enddata = (df.transfusion_datetime + pd.Timedelta("{} days".format(str(days))) > end_of_data)
	#mask_not_dead = (df[target] == 0)
	#mask_enddate_smaller_than_next_trans = (df["Date_of_next_trans"] > end_of_data)
	#mask = ((mask_succeed_enddata)&(mask_not_dead)&mask_enddate_smaller_than_next_trans)
	# Set duration_col = end_of_data - transfusion_datetime
	# 3863 edited!!
	#df.loc[mask, duration_col] = (end_of_data - df.loc[mask, "transfusion_datetime"]).dt.days
	# Transfusion datetime not needed anymore
	df.drop(["Days_to_end_of_data"],axis=1,inplace=True)

	return df


data = pd.read_csv("data/processed/chapterlevel_lpr_bloodtypes.tsv",sep="\t",low_memory=False)
data["D_STATUS_HEN_START"] = pd.to_datetime(data["D_STATUS_HEN_START"],format="%Y-%m-%d",errors="coerce")
data["D_FODDATO"] = pd.to_datetime(data["D_FODDATO"],format="%Y-%m-%d",errors="coerce")

# Get days to death
data["Death_date"] = pd.NaT
mask_dead = (data.C_STATUS == 90)
data.loc[mask_dead, "Death_date"] = data.loc[mask_dead, "D_STATUS_HEN_START"] 
data["Days_to_death"] = (data.Death_date - data.D_FODDATO).dt.days

# Get days to move
data["Move_date"] = pd.NaT
mask_move = ((data.C_STATUS != 90)&(data.C_STATUS != 1))
data.loc[mask_move, "Move_date"] = data.loc[mask_move, "D_STATUS_HEN_START"] 
data["Days_to_move"] = (data.Move_date - data.D_FODDATO).dt.days

# Get days to end of registry
end_of_data = pd.Timestamp(datetime.date(2018, 4, 10))
# Calculate days from transfusion to end of data
data["Days_to_end_of_data"] = (end_of_data - data.D_FODDATO).dt.days

# Get days to diagnosis
chapters = ["Kap"+str(x) for x in range(1,22)]
for chapter in chapters:
	# To datetime
	data[chapter+"_D_INDDTO"] = pd.to_datetime(data[chapter+"_D_INDDTO"],format="%Y-%m-%d",errors="coerce")
	# Get days to chapter
	data["Days_to_"+chapter] = pd.NaT
	data["Days_to_"+chapter] = (data[chapter+"_D_INDDTO"] - data.D_FODDATO).dt.days

	# The ones without a event get their days_to_event set to end of data or date of death
	mask = (data[chapter] == 0)
	# First according to move
	data.loc[mask,"Days_to_"+chapter] = data.loc[mask,"Days_to_move"]
	# Then according to death
	data.loc[mask&data["Days_to_"+chapter].isnull(),"Days_to_"+chapter] = data.loc[mask,"Days_to_death"]
	# Last according to end of data
	data.loc[mask&data["Days_to_"+chapter].isnull(),"Days_to_"+chapter] = data.loc[mask,"Days_to_end_of_data"]


data.to_csv()



