library(tidyverse)
library(kableExtra)
library(glue)
library(gt)

setwd("/users/projects/bloodtype/")
#setwd("/Users/pnr874/Desktop/PhD/Papers/Paper3 Bloodtypes/eLife/Data export")


## PheWAS ##

data = read_tsv("results/20220915/allVSOne/enter_registry/phewas_estimates.tsv")
data$phecode = as.character(data$phecode)
data <- data %>% mutate(phecode = str_trim(phecode))
data <- data %>% select(-c(term))


# Select subset with significant phecodes
significant_phecodes <- unique(data[data$blood_group %in% c("A","B","AB","0","RhD")&data$FDR_.05,]$phecode)
data <- data %>% filter(phecode %in% significant_phecodes & blood_group %in% c("A","B","AB","0","RhD"))


# Recalculate standard error adjusted for multiple testing
# According to: https://www.bmj.com/content/343/bmj.d2090
# 1. z-score
# Avoid exploding CIs by cutting high p-values
data <- data %>% mutate(p.value_FDR = ifelse(p.value_FDR > 0.97, 0.97,p.value_FDR))
data <- data %>% mutate(z = -0.862 + sqrt(0.743 - 2.404 * log(p.value_FDR)))
# 2. calculte standard error
# ignoring minus sign
data <- data %>% mutate(std.error = abs(log(estimate)/z))

# Recalculate confidence intervals
data <- data %>% mutate(
  CI_low = as.character(round(estimate/exp(qnorm(0.975)*std.error),2)),
  CI_high = as.character(round(estimate*exp(qnorm(0.975)*std.error),2)),
  IRR = as.character(round(estimate,2)))

# Make string of IRR and CI
data <- data %>% mutate(IRR_CI = paste(IRR," (",CI_low,", ",CI_high,")",sep=""))

data <- data %>% 
  mutate(category = str_to_title(category))


# Count unique phecodes
data %>% filter(FDR_.05) %>% distinct(phecode) %>% count()

# Count number of significant for each bloodtype
data %>% filter(FDR_.05) %>% count(blood_group)


# Median number of person years: Same for all studies as it is per phecode.
# Person-years
data %>% filter(blood_group == "A" & ! phecode %in% c("356","691","637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282","199.4","244.5","612.3","286.1","282.9","362.7","520.1","282.8","334.1","364.5")) %>%
          summarise(median = median(pyrs), Q1 = quantile(pyrs,0.25), Q3 = quantile(pyrs,0.75))


# Get number of positive and inverse associations for ABO blood groups
ABO <- data %>%
  filter(p.value_FDR < 0.05 & blood_group %in% c("A","B","AB","0")) %>%
  group_by(estimate > 1) %>%
  distinct(phecode) %>% 
  count()

# Get number of positive and inverse associations for Rhesus type
Rhesus <- data %>%
  filter(p.value_FDR < 0.05 & blood_group %in% c("RhD")) %>%
  group_by(estimate > 1) %>%
  distinct(phecode) %>% 
  count()

# Unique phecodes
unique_phecodes <- data %>%
  filter(p.value_FDR < 0.05 & blood_group %in% c("RhD","A","B","AB","0")) %>%
  #group_by(estimate > 1) %>%
  distinct(phecode) %>% 
  count()

unique_phecodes_by_group <- data %>%
  filter(p.value_FDR < 0.05 & blood_group %in% c("RhD","A","B","AB","0")) %>%
  group_by(blood_group) %>%
  distinct(phecode) %>% 
  count()

phecode_ABO <- data %>% 
  filter(p.value_FDR < 0.05 & blood_group %in% c("A","B","AB","0")) %>%
  distinct(phecode)
  # %>%
  #pull()

phecode_Rhesus <- data %>% 
  filter(p.value_FDR < 0.05 & blood_group ==  "RhD") %>%
  distinct(phecode) 

# Overlap between ABO and RhD
overlap <- phecode_ABO %>%
  filter(phecode %in% phecode_Rhesus$phecode) %>% count()

# Round person-years
data <- data %>% mutate(pyrs = round(pyrs,0))

# Make dataframe in order for kableextra with blood type as columns
data <- data %>% arrange(category_number, as.numeric(phecode)) %>%
  select(category,phecode,cases,pyrs,phenotype,blood_group,IRR_CI,p.value_FDR) %>% 
  #mutate(p.value_FDR = if_else(p.value_FDR<0.001,formatC(p.value_FDR,format="e",digits=2),as.character(round(p.value_FDR,3)))) %>%
  mutate(p.value_FDR = if_else(p.value_FDR<0.001,"<0.001",as.character(round(p.value_FDR,3)))) %>%
  pivot_wider(names_from = "blood_group",values_from = c("IRR_CI","p.value_FDR")) %>%
  select("category","phecode","phenotype","cases","pyrs",
    "IRR_CI_A","p.value_FDR_A",
    "IRR_CI_B","p.value_FDR_B",
    "IRR_CI_AB","p.value_FDR_AB",
    "IRR_CI_0","p.value_FDR_0",
    "IRR_CI_RhD","p.value_FDR_RhD")


# Put asterisk on prevalance ratio results
prevalence_phecodes3 = c("356","691","637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282")
prevalence_phecodes5 = c("199.4","244.5","612.3","286.1","282.9","362.7","520.1","282.8","334.1","364.5")
data <- data %>% 
  mutate(
    IRR_CI_A = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_A,"**",sep=""),IRR_CI_A),
    IRR_CI_B = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_B,"**",sep=""),IRR_CI_B),
    IRR_CI_AB = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_AB,"**",sep=""),IRR_CI_AB),
    IRR_CI_0 = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_0,"**",sep=""),IRR_CI_0),
    IRR_CI_RhD = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_RhD,"**",sep=""),IRR_CI_RhD),
    pyrs = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(pyrs,"*",sep=""),pyrs),
    )


data <- data %>%
  mutate(
    IRR_CI_A = cell_spec(IRR_CI_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    IRR_CI_B = cell_spec(IRR_CI_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    IRR_CI_AB = cell_spec(IRR_CI_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    IRR_CI_0 = cell_spec(IRR_CI_0,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_0,"<"))<0.05) | (p.value_FDR_0 =="<0.001"),T,F)),
    IRR_CI_RhD = cell_spec(IRR_CI_RhD,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_RhD,"<"))<0.05) | (p.value_FDR_RhD =="<0.001"),T,F)),
    p.value_FDR_A = cell_spec(p.value_FDR_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    p.value_FDR_B = cell_spec(p.value_FDR_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    p.value_FDR_AB = cell_spec(p.value_FDR_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    p.value_FDR_0 = cell_spec(p.value_FDR_0,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_0,"<"))<0.05) | (p.value_FDR_0 =="<0.001"),T,F)),
    p.value_FDR_RhD = cell_spec(p.value_FDR_RhD,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_RhD,"<"))<0.05) | (p.value_FDR_RhD =="<0.001"),T,F)))



# Print category order and names
for (category in unique(data$category)) {
   print(category)
   print(which(data$category %in% category))
}

# Remove category column
data <- data %>% select(-category)


brunak_palette = c(
"Undefined" = "#8dd3c7",
"infectious" = "#2FCFD3",
"neoplasms" = "#9E4F46",
"endocrine" = "#FF1AB9",
"hematopoietic" = "#F69EDC",
"mental" = "#2DC92D",
"neurological" = "#004DE6",
"sense" = "#FFC184",
"circulatory" = "#A00CC4",
"respiratory" = "#625D5D",
"digestive" = "#008495",
"genitourinary" = "#F9DF1E",
"pregnancy" = "#17CE95",
"dermatologic" = "#729656",
"musculoskeletal" = "#960866",
"congenital" = "#12601B",
"symptoms" = "#00007B",
"injuries" = "#FF0000")


# With diagnosis colors
kbl(data,"html",escape = FALSE,align=c("l","l","c","c","r","c","r","c","r","c","r","c","r","c"),
  col.names=c("Phecode","Phenotype","Cases","Person-years","IRR (95%CI)","P-value","IRR (95%CI)","P-value","IRR (95%CI)","P-value","IRR (95%CI)","P-value","IRR (95%CI)","P-value")) %>% 
  kable_classic("striped") %>%
  add_header_above(c(" " = 4, "Blood group A" = 2,"Blood group B" = 2,"Blood group AB" = 2,"Blood group 0" = 2,"Blood group RhD" = 2)) %>%
  pack_rows("Infectious Diseases", 1, 6, label_row_css = "background-color: #2FCFD3; color: #fff;",indent = FALSE) %>%
  pack_rows("Neoplasms", 7, 12, label_row_css = "background-color: #9E4F46; color: #fff;",indent = FALSE) %>%
  pack_rows("Endocrine/Metabolic", 13, 21, label_row_css = "background-color: #FF1AB9; color: #fff;",indent = FALSE) %>%
  pack_rows("Hematopoietic", 22, 27, label_row_css = "background-color: #F69EDC; color: #fff;",indent = FALSE) %>%
  pack_rows("Mental Disorders", 28, 28, label_row_css = "background-color: #2DC92D; color: #fff;",indent = FALSE) %>%
  pack_rows("Neurological", 29, 32, label_row_css = "background-color: #004DE6; color: #fff;",indent = FALSE) %>%
  pack_rows("Sense Organs", 33, 40, label_row_css = "background-color: #FFC184; color: #fff;",indent = FALSE) %>%
  pack_rows("Circulatory System", 41, 62, label_row_css = "background-color: #A00CC4; color: #fff;",indent = FALSE) %>%
  pack_rows("Respiratory", 63, 66, label_row_css = "background-color: #625D5D; color: #fff;",indent = FALSE) %>%
  pack_rows("Digestive", 67, 83, label_row_css = "background-color: #008495; color: #fff;",indent = FALSE) %>%
  pack_rows("Genitourinary", 84, 87, label_row_css = "background-color: #F9DF1E; color: #fff;",indent = FALSE) %>%
  pack_rows("Pregnancy Complications", 88, 95, label_row_css = "background-color: #17CE95; color: #fff;",indent = FALSE) %>%
  pack_rows("Dermatologic", 96, 98, label_row_css = "background-color: #729656; color: #fff;",indent = FALSE) %>%
  pack_rows("Musculoskeletal", 99, 108, label_row_css = "background-color: #960866; color: #fff;",indent = FALSE) %>%
  pack_rows("Congenital Anomalies", 109, 112, label_row_css = "background-color: #12601B; color: #fff;",indent = FALSE) %>%
  save_kable(file = "results/20220915/allVSOne/enter_registry/table2.html", self_contained = T)  





## Table of all phecode for suplementary

data = read_tsv("results/20220915/allVSOne/enter_registry/phewas_estimates.tsv")
data$phecode = as.character(data$phecode)

# Select only blood type info
data <- data %>% filter(blood_group %in% c("A","B","AB","0","RhD"))

# Recalculate standard error adjusted for multiple testing
# According to: https://www.bmj.com/content/343/bmj.d2090
# 1. z-score
data <- data %>% mutate(p.value_FDR = ifelse(p.value_FDR > 0.97, 0.97,p.value_FDR))
data <- data %>% mutate(z = -0.862 + sqrt(0.743 - 2.404 * log(p.value_FDR)))
# 2. calculte standard error
# ignoring minus sign
data <- data %>% mutate(std.error = abs(log(estimate)/z))


# Construct confidence intervals
data <- data %>% mutate(
  CI_low = as.character(round(estimate/exp(qnorm(0.975)*std.error),2)),
  CI_high = as.character(round(estimate*exp(qnorm(0.975)*std.error),2)),
  IRR = as.character(round(estimate,2)))

# Make string of IRR and CI
data <- data %>% mutate(IRR_CI = paste(IRR," (",CI_low,", ",CI_high,")",sep=""))

data[is.na(data$category),"category"] = "Undefined"
data <- data %>% 
  mutate(category = str_to_title(category))


# Round person-years
data <- data %>% mutate(pyrs = as.character(round(pyrs,0)))


# Make dataframe in order for kableextra with blood type as columns
data <- data %>% arrange(category_number, as.numeric(phecode)) %>%
  select(category,phecode,cases,pyrs,phenotype,blood_group,IRR_CI,p.value_FDR) %>% 
  #mutate(p.value_FDR = if_else(p.value_FDR<0.001,formatC(p.value_FDR,format="e",digits=2),as.character(round(p.value_FDR,3)))) %>%
  mutate(p.value_FDR = if_else(p.value_FDR<0.001,"<0.001",as.character(round(p.value_FDR,3)))) %>%
  pivot_wider(names_from = "blood_group",values_from = c("IRR_CI","p.value_FDR")) %>%
  select("category","phecode","phenotype","cases","pyrs",
    "IRR_CI_A","p.value_FDR_A",
    "IRR_CI_B","p.value_FDR_B",
    "IRR_CI_AB","p.value_FDR_AB",
    "IRR_CI_0","p.value_FDR_0",
    "IRR_CI_RhD","p.value_FDR_RhD")


prevalence_phecodes3 = c("356","691","637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282")
prevalence_phecodes5 = c("199.4","244.5","612.3","286.1","282.9","362.7","520.1","282.8","334.1","364.5")
data <- data %>% 
  mutate(
    IRR_CI_A = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_A,"**",sep=""),IRR_CI_A),
    IRR_CI_B = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_B,"**",sep=""),IRR_CI_B),
    IRR_CI_AB = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_AB,"**",sep=""),IRR_CI_AB),
    IRR_CI_0 = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_0,"**",sep=""),IRR_CI_0),
    IRR_CI_RhD = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_RhD,"**",sep=""),IRR_CI_RhD),
    pyrs = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(pyrs,"*",sep=""),pyrs),
    )


# Print category order and names
for (category in unique(data$category)) {
   print(category)
   print(which(data$category %in% category))
}

data <- data %>% select(-category)

data <- data %>%
  mutate(
    IRR_CI_A = cell_spec(IRR_CI_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    IRR_CI_B = cell_spec(IRR_CI_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    IRR_CI_AB = cell_spec(IRR_CI_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    IRR_CI_0 = cell_spec(IRR_CI_0,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_0,"<"))<0.05) | (p.value_FDR_0 =="<0.001"),T,F)),
    IRR_CI_RhD = cell_spec(IRR_CI_RhD,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_RhD,"<"))<0.05) | (p.value_FDR_RhD =="<0.001"),T,F)),
    p.value_FDR_A = cell_spec(p.value_FDR_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    p.value_FDR_B = cell_spec(p.value_FDR_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    p.value_FDR_AB = cell_spec(p.value_FDR_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    p.value_FDR_0 = cell_spec(p.value_FDR_0,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_0,"<"))<0.05) | (p.value_FDR_0 =="<0.001"),T,F)),
    p.value_FDR_RhD = cell_spec(p.value_FDR_RhD,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_RhD,"<"))<0.05) | (p.value_FDR_RhD =="<0.001"),T,F)))


kbl(data,"html",escape = F,col.names=c("Phecode","Phenotype","Cases","Person-years","IRR (95%CI)","P-value","IRR (95%CI)","P-value","IRR (95%CI)","P-value","IRR (95%CI)","P-value","IRR (95%CI)","P-value"),
    align=c("l","l","c","c","r","c","r","c","r","c","r","c","r","c")) %>% 
  kable_classic("striped") %>%
  add_header_above(c(" " = 4, "Blood group A" = 2,"Blood group B" = 2,"Blood group AB" = 2,"Blood group 0" = 2,"Blood group RhD" = 2)) %>%
  pack_rows("Infectious Diseases", 1, 54, label_row_css = "background-color: #2FCFD3; color: #fff;",indent = FALSE) %>%
  pack_rows("Neoplasms", 55, 170, label_row_css = "background-color: #9E4F46; color: #fff;",indent = FALSE) %>%
  pack_rows("Endocrine/Metabolic", 171, 287, label_row_css = "background-color: #FF1AB9; color: #fff;",indent = FALSE) %>%
  pack_rows("Hematopoietic", 288, 329, label_row_css = "background-color: #F69EDC; color: #fff;",indent = FALSE) %>%
  pack_rows("Mental Disorders", 330, 393, label_row_css = "background-color: #2DC92D; color: #fff;",indent = FALSE) %>%
  pack_rows("Neurological", 394, 461, label_row_css = "background-color: #004DE6; color: #fff;",indent = FALSE) %>%
  pack_rows("Sense Organs", 462, 564, label_row_css = "background-color: #FFC184; color: #fff;",indent = FALSE) %>%
  pack_rows("Circulatory System", 565, 701, label_row_css = "background-color: #A00CC4; color: #fff;",indent = FALSE) %>%
  pack_rows("Respiratory", 702, 769, label_row_css = "background-color: #625D5D; color: #fff;",indent = FALSE) %>%
  pack_rows("Digestive", 770, 903, label_row_css = "background-color: #008495; color: #fff;",indent = FALSE) %>%
  pack_rows("Genitourinary", 904, 1031, label_row_css = "background-color: #F9DF1E; color: #fff;",indent = FALSE) %>%
  pack_rows("Pregnancy Complications", 1032, 1084, label_row_css = "background-color: #17CE95; color: #fff;",indent = FALSE) %>%
  pack_rows("Dermatologic", 1085, 1162, label_row_css = "background-color: #729656; color: #fff;",indent = FALSE) %>%
  pack_rows("Musculoskeletal", 1163, 1260, label_row_css = "background-color: #960866; color: #fff;",indent = FALSE) %>%
  pack_rows("Congenital Anomalies", 1261, 1312, label_row_css = "background-color: #12601B; color: #fff;",indent = FALSE) %>%
  save_kable(file = "results/20220915/allVSOne/enter_registry/supplementary_table2.html", self_contained = T) 




# Age at first diagnosis #

data = read_tsv("results/20220915/allVSOne/age_at_diagnosis/phewas_estimates.tsv")
data$phecode = as.character(data$phecode)

# Select subset with significant phecodes
# Get significant phecodes
significant_phecodes <- unique(data[data$blood_group %in% c("A","B","AB","0","RhD")&round(data$p.value_FDR,3)<0.05,]$phecode)


# Recalculate standard error adjusted for multiple testing
# According to: https://www.bmj.com/content/343/bmj.d2090
# 1. z-score
data <- data %>% mutate(p.value_FDR = ifelse(p.value_FDR > 0.97, 0.97,p.value_FDR))
data <- data %>% mutate(z = -0.862 + sqrt(0.743 - 2.404 * log(p.value_FDR)))
# 2. calculte standard error
# ignoring minus sign
data <- data %>% mutate(std.error = abs(estimate/z))


# Construct confidence intervals
data <- data %>% mutate(
  CI_low = as.character(round(estimate-qnorm(0.975)*std.error,2)),
  CI_high = as.character(round(estimate+qnorm(0.975)*std.error,2)),
  Estimate = as.character(round(estimate,2)))

# Make string of IRR and CI
data <- data %>% mutate(Estimate_CI = paste(Estimate," (",CI_low,", ",CI_high,")",sep=""))

# Make dataframe in order for kableextra with blood type as columns
data <- data %>% arrange(as.numeric(phecode)) %>%
  select(phecode,phenotype,N,blood_group,Estimate_CI,p.value_FDR) %>% 
  #mutate(p.value_FDR = if_else(p.value_FDR<0.001,formatC(p.value_FDR,format="e",digits=2),as.character(round(p.value_FDR,3)))) %>%
  mutate(p.value_FDR = if_else(p.value_FDR<0.001,"<0.001",as.character(round(p.value_FDR,3)))) %>%
  pivot_wider(names_from = "blood_group",values_from = c("Estimate_CI","p.value_FDR")) %>%
  select("phecode","phenotype","N","Estimate_CI_A","p.value_FDR_A","Estimate_CI_B","p.value_FDR_B",
    "Estimate_CI_AB","p.value_FDR_AB","Estimate_CI_0","p.value_FDR_0","Estimate_CI_RhD","p.value_FDR_RhD")


# Color bold significant hits
data <- data %>%
  mutate(
    Estimate_CI_A = cell_spec(Estimate_CI_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    Estimate_CI_B = cell_spec(Estimate_CI_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    Estimate_CI_AB = cell_spec(Estimate_CI_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    Estimate_CI_0 = cell_spec(Estimate_CI_0,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_0,"<"))<0.05) | (p.value_FDR_0 =="<0.001"),T,F)),
    Estimate_CI_RhD = cell_spec(Estimate_CI_RhD,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_RhD,"<"))<0.05) | (p.value_FDR_RhD =="<0.001"),T,F)),
    p.value_FDR_A = cell_spec(p.value_FDR_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    p.value_FDR_B = cell_spec(p.value_FDR_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    p.value_FDR_AB = cell_spec(p.value_FDR_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    p.value_FDR_0 = cell_spec(p.value_FDR_0,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_0,"<"))<0.05) | (p.value_FDR_0 =="<0.001"),T,F)),
    p.value_FDR_RhD = cell_spec(p.value_FDR_RhD,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_RhD,"<"))<0.05) | (p.value_FDR_RhD =="<0.001"),T,F)))


# Save to file # Only works in R studio
#caption = "Associations between ABO blood group, Rhesus type and diagnosis rate ratio"
kbl(data,"html",escape = F,col.names=c("Phecode","Phenotype","N","Estimate (95%CI)","P-value","Estimate (95%CI)","P-value","Estimate (95%CI)","P-value","Estimate (95%CI)","P-value","Estimate (95%CI)","P-value"),
      align=c("l","l","c","c","c","c","c","c","c","c","c","c","c","c")) %>% 
  kable_classic("striped") %>%
  add_header_above(c(" " = 3, "Blood group A" = 2,"Blood group B" = 2,"Blood group AB" = 2,"Blood group 0" = 2,"Blood group RhD" = 2)) %>%
  save_kable(file = "results/20220915/allVSOne/age_at_diagnosis/suplementary_table3.html", self_contained = T)


# Only significante
data <- data %>% filter(phecode %in% significant_phecodes)
kbl(data,"html",escape = F,col.names=c("Phecode","Phenotype","N","Estimate (95%CI)","P-value","Estimate (95%CI)","P-value","Estimate (95%CI)","P-value","Estimate (95%CI)","P-value","Estimate (95%CI)","P-value"),
      align=c("l","l","c","c","c","c","c","c","c","c","c","c","c","c")) %>%  
  kable_classic("striped") %>%
  add_header_above(c(" " = 3, "Blood group A" = 2,"Blood group B" = 2,"Blood group AB" = 2,"Blood group 0" = 2,"Blood group RhD" = 2)) %>%
  save_kable(file = "results/20220915/allVSOne/age_at_diagnosis/table3.html", self_contained = T)


