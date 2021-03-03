args <- commandArgs(trailingOnly=TRUE)
# args[1] - model_input
# args[2] - estimates file

library(tidyverse)
library(broom)
library(survival)
#library(kableExtra)
library(rms)
#library(gtsummary)
#library(tab)

# Load data
data <- read_tsv(args[1])

#setwd("/home/people/petras/bloodtype_diseases/")
#data = read_tsv("data/processed/allcause_mortality/allcause_mortality.tsv")

# Remove patients with no time in study / patient developing disease before bloodtype meassure
#data <- data[data$age_at_entry != data$age_at_exit,]
data <- data[data$age_at_entry < data$age_at_exit,]

#BIRTH_YEARS_KNOTS <- rcspline.eval(data$Birth_year, nk=5, knots.only=TRUE)

# Split in 1 age intervals
age_cut <- seq(1,110,1)
data <- survSplit(Surv(age_at_entry,age_at_exit,CASE) ~., data= data, cut=age_cut, episode="age_group", id="id")

AGE_KNOTS <- rcspline.eval(data$age_group, nk=5, knots.only=TRUE)

# Make person-years column for each time window of patient
data$pyrs = data$age_at_exit - data$age_at_entry

# Group data into count table 
data <- data %>% group_by(Birth_year,AB0,Rhesus,C_KON,age_group) %>% summarise(pyrs = sum(pyrs), cases = sum(CASE))
# Per 1000 person years
data$pyrs <- data$pyrs/1000

fit_interact <- glm(cases ~ factor(Birth_year) + rcs(age_group,knots=AGE_KNOTS) + factor(C_KON) + rcs(age_group,knots=AGE_KNOTS):factor(C_KON) + factor(AB0) + factor(AB0):rcs(age_group,knots=AGE_KNOTS) + factor(Rhesus) + factor(Rhesus):rcs(age_group,knots=AGE_KNOTS) + offset(log(pyrs)),
	data,family = quasipoisson(link = "log"))
#summary(fit_interact)
#exp(coef(fit_interact))

# Without interaction
fit <- glm(cases ~ factor(Birth_year) + rcs(age_group,knots=AGE_KNOTS) + factor(C_KON) + rcs(age_group,knots=AGE_KNOTS):factor(C_KON) + factor(AB0) + factor(Rhesus) + offset(log(pyrs)),
	data,family = quasipoisson(link = "log"))

# Test if interaction is significant
txt = anova(fit,fit_interact,test="Chisq") # its not significant

sink(args[3])
#sink("results/20210128/poisson/enter_registry/allcause_mortality/anova_test.txt")
print(txt)
sink()


# Save estimates
tmp <- tidy(fit)
tmp$estimate <- exp(tmp$estimate)
#tmp$conf.low = exp(tmp$conf.low)
#tmp$conf.high = exp(tmp$conf.high)
tmp$pyrs = sum(data$pyrs) * 1000
tmp$cases = sum(data$cases)
write_tsv(tmp,args[2])


# Save interaction model
tmp <- tidy(fit_interact)
tmp$estimate <- exp(tmp$estimate)
tmp$pyrs = sum(data$pyrs) * 1000
tmp$cases = sum(data$cases)
write_tsv(tmp,args[4])



