args <- commandArgs(trailingOnly=TRUE)
# args[1] - model_input
# args[2] - estimates file

library(tidyverse)
library(broom)
library(survival)
library(rms)
#library(kableExtra)
#library(gtsummary)
#library(tab)


# Load data
data <- read_tsv(args[1])
setwd("/home/people/petras/bloodtype_diseases/")
#data = read_tsv("data/processed/PheWAS/174.tsv")
#data = read_tsv("data/processed/PheWAS/285.1.tsv")
#data = read_tsv("data/processed/ICD10level3/C51.tsv")
#data = read_tsv("data/processed/ICD10level3/N90.tsv")

diagnosis = str_remove(tail(str_split(args[1],"/")[[1]],1), ".tsv")

#diagnosis = "656.2"
#diagnosis = "619.3"
#diagnosis = "259.4"
#diagnosis = "637"

#diagnosis = "285.1"
#data = read_tsv(paste("data/processed/PheWAS/",diagnosis,".tsv",sep=""))

# load phecode definitions for sex specific diagnosis
phecodes_def = "/home/people/petras/base_data/phecode_def.csv"
phecodes_def <- read_csv(phecodes_def)
sex_phecodes <- phecodes_def %>% filter(sex == "Female"|sex == "Male") %>% select("phecode","sex")
# Replace with M and K
sex_phecodes[sex_phecodes$sex=="Female","sex"] = "K"
sex_phecodes[sex_phecodes$sex=="Male","sex"] = "M"


# Determine if its a single-sex analysis
if (diagnosis %in% pull(sex_phecodes["phecode"])) {
	single_sex = TRUE
	# PheWAS
	sex <- sex_phecodes %>% 
		filter(phecode == diagnosis) %>%
		select(sex) %>%
		pull()
	data <- data %>% filter(C_KON == sex)

} else if (substr(diagnosis, 1, 3) %in% c("N90","C51","C52","C53",
	"C54","C55","C56","C57","C58","D06","D25","D26","D27",
	"D28","D39","N70","N71","N72","N73","N74","N75","N76","N77","N80",
	"N81","N82","N83","N84","N85","N88","N89","N90","N92","N93","N94",
	"N95","N96","N97","N98","Q50","Q51","Q52") |
	substr(diagnosis, 1, 1) %in% c("O")) {
	single_sex = TRUE
	# Female ICD10
	sex <- "K"
	data <- data %>% filter(C_KON == sex)

} else if (substr(diagnosis, 1, 3) %in% c("C60","C61","C62","C63",
	"D29","D40","N40","N41","N42","N43","N44","N45","N46","N47",
	"N48","N49","N50","N51","Q53","Q54","Q55","Q64")) {
	single_sex = TRUE
	# Male ICD10
	sex <- "M"
	data <- data %>% filter(C_KON == sex)

} else if (data %>% filter(CASE==1) %>% distinct(C_KON) %>% count() == 1) {
	single_sex = TRUE
	# Subset cohort with only that sex
	sex <- unique(data[data$CASE==1,"C_KON"])[[1]][1]
	data <- data %>% filter(C_KON == sex)
} else {
	single_sex = FALSE
}

# Calculate knots using full data
#BIRTH_YEARS_KNOTS <- rcspline.eval(data$Birth_year, nk=5, knots.only=TRUE)

#Hereditary
#added = "282.9","356","362.7","520.1","665","282.8","282","334.1","364.5"  ,637 ,"658","661","657"
#added = "D55","D58","G60","K00","E80","G11","G60","D56","D57"

#PheWAS: 656 + 691 + 244.5 + 286.1 + 612.3 + 747-759.99
#ICD10: P + Q diagnosis
perinatal = FALSE
if (substr(diagnosis, 1, 3) %in% c("356","691","637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282") | 
	substr(diagnosis, 1, 5) %in% c("199.4","244.5","612.3","286.1","282.9","362.7","520.1","282.8","334.1","364.5") |  
	substr(diagnosis, 1, 1) %in% c("P","Q") | 
	substr(diagnosis, 1, 3) %in% c("E00","D66","D67","D68","D55","D58","G60","K00","E80","G11","G60","D56","D57")) {
	perinatal = TRUE
	# Dont filter if age at entry is smaller than age at exit
	# For this analysis we dont use time/person-years but population size
	data$pop = 1
} else {
	# Remove patients with no time in study / patient developing disease before bloodtype meassure
	data <- data[data$age_at_entry < data$age_at_exit,]
	# Split in 1-age intervals: 21 categories
	age_cut <- seq(1,115,1)
	data <- survSplit(Surv(age_at_entry,age_at_exit,CASE) ~., data= data, cut=age_cut, episode="age_group", id="id")
	# Calculate knots using full data
	AGE_KNOTS <- rcspline.eval(data$age_group, nk=5, knots.only=TRUE)
	# Make person-years column for each time window of patient
	data$pyrs = data$age_at_exit - data$age_at_entry
}

# Split birth year into groups of 5 year intervals
#data$Birth_year_bin <- cut(data$Birth_year,seq(1900,2020,5))

# Run analysis sex specific or not
if (perinatal){
	# Perinatal or congenital diagnosis
	if (single_sex) {
	# Only one sex
	# Group data into count table
	data <- data %>% group_by(Birth_year,AB0,Rhesus) %>% summarise(population = sum(pop), cases = sum(CASE))
	# Relative risk ratio
	fit_interact <- glm(cases ~ factor(Birth_year) + factor(AB0) + factor(Rhesus) + offset(log(population)),
		data,family = quasipoisson(link = "log"),maxit=200)

	} else {
	# Group data into count table
	data <- data %>% group_by(Birth_year,C_KON,AB0,Rhesus) %>% summarise(population = sum(pop), cases = sum(CASE))
	# Relative risk ratio
	fit_interact <- glm(cases ~ factor(Birth_year) + factor(C_KON) + factor(AB0) + factor(Rhesus) + offset(log(population)),
		data,family = quasipoisson(link = "log"),maxit=200)
	}
} else {
	if (single_sex) {
		# Only one sex
		# Group data into count table 
		data <- data %>% group_by(Birth_year,AB0,Rhesus,age_group) %>% summarise(pyrs = sum(pyrs), cases = sum(CASE))
		# Per 1000 person years
		data$pyrs <- data$pyrs/1000
		# Fit factor(Birth_year)
		fit_interact <- glm(cases ~ factor(Birth_year) + rcs(age_group,knots=AGE_KNOTS) + factor(AB0) + factor(Rhesus) + offset(log(pyrs)),
			data,family = quasipoisson(link = "log"),maxit=200)
	} else {
		# Both sex
		# Group data into count table 
		data <- data %>% group_by(Birth_year,AB0,Rhesus,C_KON,age_group) %>% summarise(pyrs = sum(pyrs), cases = sum(CASE))
		# Per 1000 person years
		data$pyrs <- data$pyrs/1000
		# Fit
		fit_interact <- glm(cases ~ factor(Birth_year) + rcs(age_group,knots=AGE_KNOTS) + factor(C_KON) + rcs(age_group,knots=AGE_KNOTS):factor(C_KON) + factor(AB0) + factor(Rhesus) + offset(log(pyrs)),
			data,family = quasipoisson(link = "log"),maxit=200)
	}
}


# Save estimates
tmp <- tidy(fit_interact)
tmp$estimate <- exp(tmp$estimate)
#tmp$conf.low = exp(tmp$conf.low)
#tmp$conf.high = exp(tmp$conf.high)
if (perinatal){
	tmp$pyrs = sum(data$population)
} else {
	tmp$pyrs = sum(data$pyrs) * 1000
}

tmp$cases = sum(data$cases)

# count for each blood group
cases <- data %>% group_by(AB0) %>% summarise(cases=sum(cases))

#O = as.character(unique(cases[cases$AB0 == "0",][2]))
A = as.character(unique(cases[cases$AB0 == "A",][2]))
B = as.character(unique(cases[cases$AB0 == "B",][2]))
AB = as.character(unique(cases[cases$AB0 == "AB",][2]))

cases <- data %>% group_by(Rhesus) %>% summarise(cases=sum(cases))
Rh_positiv = as.character(unique(cases[cases$Rhesus == "Positiv",][2]))

tmp$cases_A = A
tmp$cases_B = B
tmp$cases_AB = AB
tmp$cases_Rh = Rh_positiv

write_tsv(tmp,args[2])

# Testing


# # Group data into count table 
# data <- data %>% group_by(Birth_year,AB0,Rhesus,C_KON,age_group) %>% summarise(pyrs = sum(pyrs), cases = sum(CASE))
# # Per 1000 person years
# data$pyrs <- data$pyrs/1000
# # Fit
# fit_interact <- glm(cases ~ factor(Birth_year) + rcs(age_group,knots=AGE_KNOTS) + factor(C_KON) + rcs(age_group,knots=AGE_KNOTS):factor(C_KON) + factor(AB0) + factor(Rhesus) + offset(log(pyrs)),
# 	data,family = quasipoisson(link = "log"),maxit=200)

# fit_interact <- glm(cases ~ rcs(Birth_year,5) + rcs(age_group,knots=AGE_KNOTS) + factor(C_KON) + rcs(age_group,knots=AGE_KNOTS):factor(C_KON) + factor(AB0) + factor(Rhesus) + offset(log(pyrs)),
# 	data,family = quasipoisson(link = "log"),maxit=200)

# fit_interact <- glm(cases ~ rcs(Birth_year,knots=BIRTH_YEARS_KNOTS) + rcs(age_group,knots=AGE_KNOTS) + factor(C_KON) + rcs(age_group,knots=AGE_KNOTS):factor(C_KON) + factor(AB0) + factor(Rhesus) + offset(log(pyrs)),
# 	data,family = quasipoisson(link = "log"),maxit=200)



# tmp$estimate2 <- exp(tmp$estimate)

# tmp <- tmp %>% mutate(
#   CI_low = as.character(round(estimate2/exp(1.96*std.error),2)),
#   CI_high = as.character(round(estimate2*exp(1.96*std.error),2)),
#   IRR = as.character(round(estimate2,2)))

# tmp %>% filter(term %in% c("factor(C_KON)M","factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv"))

