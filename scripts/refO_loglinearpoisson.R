args <- commandArgs(trailingOnly=TRUE)
# args[1] - model_input
# args[2] - estimates file

library(tidyverse)
library(broom)
library(survival)
library(rms)


# Load data
setwd("/users/projects/bloodtype/")
data <- read_tsv(args[1])
#data = read_tsv("data/processed/PheWAS/174.tsv")

# Define exposure variable
#blood_group = as.character(args[3])
#data <- data %>% mutate(exposure = if_else(AB0 == blood_group, 1, 0))


diagnosis = str_remove(tail(str_split(args[1],"/")[[1]],1), ".tsv")
#diagnosis = "174"

# load phecode definitions for sex specific diagnosis
phecodes_def = "/users/people/petras/base_data/phecode_def.csv"
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

} else if (data %>% filter(CASE==1) %>% distinct(C_KON) %>% count() == 1) {
	single_sex = TRUE
	# Subset cohort with only that sex
	sex <- unique(data[data$CASE==1,"C_KON"])[[1]][1]
	data <- data %>% filter(C_KON == sex)
} else {
	single_sex = FALSE
}


perinatal = FALSE
if (substr(diagnosis, 1, 3) %in% c("356","691","637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282") | 
	substr(diagnosis, 1, 5) %in% c("199.4","244.5","612.3","286.1","282.9","362.7","520.1","282.8","334.1","364.5")) {
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
	# Make person-years column for each time window of patient
	data$pyrs = data$age_at_exit - data$age_at_entry
}

# Spline knots
BIRTH_YEARS_KNOTS <- c(1920,1940,1960,1980,2000)
# Calculate knots using full data
AGE_KNOTS <- c(20,40,60,70,90)


# Run analysis sex specific or not
if (perinatal){
	# Perinatal or congenital diagnosis
	if (single_sex) {
	# Only one sex
	# Group data into count table
	data <- data %>% group_by(Birth_year,AB0,Rhesus) %>% summarise(population = sum(pop), cases = sum(CASE))
	# Relative risk ratio
	fit_interact <- glm(cases ~ rcs(Birth_year,knots = BIRTH_YEARS_KNOTS) + factor(AB0) + factor(Rhesus) + offset(log(population)),
		data,family = quasipoisson(link = "log"),maxit=200)

	} else {
	# Group data into count table
	data <- data %>% group_by(Birth_year,C_KON,AB0,Rhesus) %>% summarise(population = sum(pop), cases = sum(CASE))
	# Relative risk ratio
	fit_interact <- glm(cases ~ rcs(Birth_year,knots = BIRTH_YEARS_KNOTS) + factor(C_KON) + factor(AB0) + factor(Rhesus) + offset(log(population)),
		data,family = quasipoisson(link = "log"),maxit=200)
	}
} else {
	if (single_sex) {
		# Only one sex
		# Group data into count table 
		data <- data %>% group_by(Birth_year,AB0,Rhesus,age_group) %>% summarise(pyrs = sum(pyrs), cases = sum(CASE))
		# Per 1000 person years
		data$pyrs <- data$pyrs/1000
		# Fit rcs(Birth_year,knots = BIRTH_YEARS_KNOTS)
		fit_interact <- glm(cases ~ rcs(Birth_year,knots = BIRTH_YEARS_KNOTS) + rcs(age_group,knots=AGE_KNOTS) + factor(AB0) + factor(Rhesus) + offset(log(pyrs)),
			data,family = quasipoisson(link = "log"),maxit=200)
	} else {
		# Both sex
		# Group data into count table 
		data <- data %>% group_by(Birth_year,AB0,Rhesus,C_KON,age_group) %>% summarise(pyrs = sum(pyrs), cases = sum(CASE))
		# Per 1000 person years
		data$pyrs <- data$pyrs/1000
		# Fit
		fit_interact <- glm(cases ~ rcs(Birth_year,knots = BIRTH_YEARS_KNOTS) + rcs(age_group,knots=AGE_KNOTS) + factor(C_KON) + rcs(age_group,knots=AGE_KNOTS):factor(C_KON) + factor(AB0) + factor(Rhesus) + offset(log(pyrs)),
			data,family = quasipoisson(link = "log"),maxit=200)
	}
}


# Save estimates
tmp <- tidy(fit_interact)
tmp <- tmp %>% 
	filter(term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv"))
tmp$estimate <- exp(tmp$estimate)
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
