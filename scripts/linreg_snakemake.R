args <- commandArgs(trailingOnly=TRUE)
# args[1] - model_input
# args[2] - estimates file

library(tidyverse)
library(broom)
#library(survival)
library(rms)
#library(kableExtra)

#library(gtsummary)
#library(tab)

# Load data
data <- read_tsv(args[1])
#setwd("/home/people/petras/bloodtype_diseases/")
#data = read_tsv("data/processed/age_at_diagnosis/phewas/285.1.tsv")

# load phecode definitions for sex specific diagnosis
phecodes_def = "/home/people/petras/base_data/phecode_def.csv"
phecodes_def <- read_csv(phecodes_def)
sex_phecodes <- phecodes_def %>% filter(sex == "Female"|sex == "Male") %>% select("phecode","sex")
# Replace with M and K
sex_phecodes[sex_phecodes$sex=="Female","sex"] = "K"
sex_phecodes[sex_phecodes$sex=="Male","sex"] = "M"

# Get diagnosis
diagnosis = str_remove(tail(str_split(args[1],"/")[[1]],1), ".tsv")
#diagnosis = "285.1"

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

} else if (length(unique(data$C_KON)) == 1) {
	single_sex = TRUE
	# Subset cohort with only that sex
	sex <- unique(data$C_KON)[[1]][1]
	data <- data %>% filter(C_KON == sex)
} else {
	single_sex = FALSE
}



# With splines
# If sex specific dont control for sex
if (single_sex) {
   fit <- glm(age_at_diagnosis ~ rcs(Birth_year,5) + factor(AB0) + factor(Rhesus),
	data,family = "gaussian",maxit=200)
} else {
   fit <- glm(age_at_diagnosis ~ rcs(Birth_year,5) + factor(C_KON) + factor(AB0) + factor(Rhesus),
	data,family = "gaussian",maxit=200)
}

# As category
# If sex specific dont control for sex
# if (single_sex) {
#    fit <- glm(age_at_diagnosis ~ factor(Birth_year) + factor(AB0) + factor(Rhesus),
# 	data,family = "gaussian",maxit=200)
# } else {
#    fit <- glm(age_at_diagnosis ~ factor(Birth_year) + factor(C_KON) + factor(AB0) + factor(Rhesus),
# 	data,family = "gaussian",maxit=200)
# }



#summary(fit_interact)
#exp(coef(fit_interact))


# Save estimates
tmp <- tidy(fit)
# Get size of population
tmp$N = dim(data)[1]
# Save
write_tsv(tmp,args[2])



