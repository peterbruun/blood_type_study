args <- commandArgs(trailingOnly=TRUE)
# args[1] - model_input
# args[2] - estimates file

library(tidyverse)
library(broom)
library(rms)


setwd("/users/projects/bloodtype/")

# Load data
data <- read_tsv(args[1])
#data = read_tsv("data/processed/age_at_diagnosis/phewas/285.1.tsv")
blood_group = as.character(args[3])
data <- data %>% mutate(exposure = if_else(AB0 == blood_group, 1, 0))

# load phecode definitions for sex specific diagnosis
phecodes_def = "/users/people/petras/base_data/phecode_def.csv"
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

} else if (length(unique(data$C_KON)) == 1) {
	single_sex = TRUE
	# Subset cohort with only that sex
	sex <- unique(data$C_KON)[[1]][1]
	data <- data %>% filter(C_KON == sex)
} else {
	single_sex = FALSE
}


BIRTH_YEARS_KNOTS <- c(1920,1940,1960,1980,2000)


# With splines
# If sex specific dont control for sex
if (single_sex) {
   fit <- glm(age_at_diagnosis ~ rcs(Birth_year,knots = BIRTH_YEARS_KNOTS) + factor(exposure) + factor(Rhesus),
	data,family = "gaussian",maxit=200)
} else {
   fit <- glm(age_at_diagnosis ~ rcs(Birth_year,knots = BIRTH_YEARS_KNOTS) + factor(C_KON) + factor(exposure) + factor(Rhesus),
	data,family = "gaussian",maxit=200)
}



#summary(fit_interact)
#exp(coef(fit_interact))


# Save estimates
tmp <- tidy(fit)
# Get size of population
tmp$N = dim(data)[1]
# Save
write_tsv(tmp,args[2])



