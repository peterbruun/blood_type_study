library(tidyverse)
library(kableExtra)
library(glue)
library(gt)


setwd("/users/projects/bloodtype/")

## PheWAS ##

#data = read_tsv(snakemake@input[["phewas"]])

data = read_tsv("results/20220915/referenceO/enter_registry/phewas_estimates.tsv")
data$phecode = as.character(data$phecode)
data <- data %>% mutate(phecode = str_trim(phecode))

# Select subset with significant phecodes
significant_phecodes <- unique(data[data$term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv")&data$FDR_.05,]$phecode)
data <- data %>% filter(phecode %in% significant_phecodes & term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv"))


# Recalculate standard error adjusted for multiple testing
# According to: https://www.bmj.com/content/343/bmj.d2090

# Avoid exploding CIs by cutting high p-values
data <- data %>% mutate(p.value_FDR = ifelse(p.value_FDR > 0.97, 0.97,p.value_FDR))
# 1. adjusted  z-score
data <- data %>% mutate(z = -0.862 + sqrt(0.743 - 2.404 * log(p.value_FDR)))
# 2. calculte adjusted standard error
# ignoring minus sign
data <- data %>% mutate(std.error = abs(log(estimate)/z))

# Calculate adjusted confidence intervals
data <- data %>% mutate(
  CI_low = as.character(round(estimate/exp(qnorm(0.975)*std.error),2)),
  CI_high = as.character(round(estimate*exp(qnorm(0.975)*std.error),2)),
  IRR = as.character(round(estimate,2)))

# Make string of IRR and CI
data <- data %>% mutate(IRR_CI = paste(IRR," (",CI_low,", ",CI_high,")",sep=""))

#data[is.na(data$category),"category"] = "Other"
data <- data %>% 
  mutate(category = str_to_title(category))

# Count unique phecodes
data %>% filter(FDR_.05) %>% distinct(phecode) %>% count()

# Count number of significant for each bloodtype
data %>% filter(FDR_.05) %>% count(term)


# Median number of person years: Same for all studies as it is per phecode.
# Person-years
data %>% filter(term == "factor(AB0)A" & ! phecode %in% c("356","691","637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282","199.4","244.5","612.3","286.1","282.9","362.7","520.1","282.8","334.1","364.5")) %>%
          summarise(median = median(pyrs), Q1 = quantile(pyrs,0.25), Q3 = quantile(pyrs,0.75))


# Get number of positive and inverse associations for ABO blood groups
ABO <- data %>%
  filter(p.value_FDR < 0.05 & term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB")) %>%
  group_by(estimate > 1) %>%
  distinct(phecode) %>% 
  count()

# Get number of positive and inverse associations for Rhesus type
Rhesus <- data %>%
  filter(p.value_FDR < 0.05 & term %in% c("factor(Rhesus)Positiv")) %>%
  group_by(estimate > 1) %>%
  distinct(phecode) %>% 
  count()

# Unique phecodes
unique_phecodes <- data %>%
  filter(p.value_FDR < 0.05 & term %in% c("factor(Rhesus)Positiv","factor(AB0)A","factor(AB0)B","factor(AB0)AB")) %>%
  #group_by(estimate > 1) %>%
  distinct(phecode) %>% 
  count()

unique_phecodes_by_group <- data %>%
  filter(p.value_FDR < 0.05 & term %in% c("factor(Rhesus)Positiv","factor(AB0)A","factor(AB0)B","factor(AB0)AB")) %>%
  group_by(term) %>%
  distinct(phecode) %>% 
  count()

phecode_ABO <- data %>% 
  filter(p.value_FDR < 0.05 & term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB")) %>%
  distinct(phecode)
  # %>%
  #pull()

phecode_Rhesus <- data %>% 
  filter(p.value_FDR < 0.05 & term ==  "factor(Rhesus)Positiv") %>%
  distinct(phecode) 

# Overlap between ABO and RhD
overlap <- phecode_ABO %>%
  filter(phecode %in% phecode_Rhesus$phecode) %>% count()


# Make dataframe in order for kableextra with blood type as columns
data <- data %>% arrange(category_number, as.numeric(phecode)) %>%
  select(category,phecode,cases,cases_A,cases_B,cases_AB,cases_Rh,pyrs,phenotype,term,IRR_CI,p.value_FDR) %>% 
  #mutate(p.value_FDR = if_else(p.value_FDR<0.001,formatC(p.value_FDR,format="e",digits=2),as.character(round(p.value_FDR,3)))) %>%
  mutate(p.value_FDR = if_else(p.value_FDR<0.001,"<0.001",as.character(round(p.value_FDR,3)))) %>%
  pivot_wider(names_from = "term",values_from = c("IRR_CI","p.value_FDR")) %>%
  select("category","phecode","phenotype","cases","pyrs",
    "cases_A","IRR_CI_factor(AB0)A","p.value_FDR_factor(AB0)A",
    "cases_B","IRR_CI_factor(AB0)B","p.value_FDR_factor(AB0)B",
    "cases_AB","IRR_CI_factor(AB0)AB","p.value_FDR_factor(AB0)AB",
    "cases_Rh","IRR_CI_factor(Rhesus)Positiv","p.value_FDR_factor(Rhesus)Positiv")


# Rename columns not good with ( ) in column names
data <- data %>% rename(IRR_CI_A = "IRR_CI_factor(AB0)A",IRR_CI_B = "IRR_CI_factor(AB0)B",IRR_CI_AB = "IRR_CI_factor(AB0)AB",IRR_CI_Rhesus = "IRR_CI_factor(Rhesus)Positiv")
data <- data %>% rename(p.value_FDR_A = "p.value_FDR_factor(AB0)A",p.value_FDR_B = "p.value_FDR_factor(AB0)B",p.value_FDR_AB = "p.value_FDR_factor(AB0)AB",p.value_FDR_Rhesus = "p.value_FDR_factor(Rhesus)Positiv")

# Round person-years
data <- data %>% mutate(pyrs = round(pyrs,0))

# Put asterisk on prevalance ratio results
prevalence_phecodes3 = c("356","691","637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282")
prevalence_phecodes5 = c("199.4","244.5","612.3","286.1","282.9","362.7","520.1","282.8","334.1","364.5")
data <- data %>% 
  mutate(
    IRR_CI_A = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_A,"**",sep=""),IRR_CI_A),
    IRR_CI_B = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_B,"**",sep=""),IRR_CI_B),
    IRR_CI_AB = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_AB,"**",sep=""),IRR_CI_AB),
    IRR_CI_Rhesus = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_Rhesus,"**",sep=""),IRR_CI_Rhesus),
    pyrs = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(pyrs,"*",sep=""),pyrs),
    )


# Color bold significant hits
data <- data %>%
  mutate(
    IRR_CI_A = cell_spec(IRR_CI_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    IRR_CI_B = cell_spec(IRR_CI_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    IRR_CI_AB = cell_spec(IRR_CI_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    IRR_CI_Rhesus = cell_spec(IRR_CI_Rhesus,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_Rhesus,"<"))<0.05) | (p.value_FDR_Rhesus =="<0.001"),T,F)),
    p.value_FDR_A = cell_spec(p.value_FDR_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    p.value_FDR_B = cell_spec(p.value_FDR_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    p.value_FDR_AB = cell_spec(p.value_FDR_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    p.value_FDR_Rhesus = cell_spec(p.value_FDR_Rhesus,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_Rhesus,"<"))<0.05) | (p.value_FDR_Rhesus =="<0.001"),T,F)))



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
#"chapterColors8" = "#A7A7F9",
"circulatory" = "#A00CC4",
"respiratory" = "#625D5D",
"digestive" = "#008495",
"genitourinary" = "#F9DF1E",
"pregnancy" = "#17CE95",
"dermatologic" = "#729656",
"musculoskeletal" = "#960866",
"congenital" = "#12601B",
#"chapterColors17" = "#F68409",
"symptoms" = "#00007B",
"injuries" = "#FF0000")


# With new mapping
# With diagnosis colors
kbl(data,"html",escape = F,col.names=c("Phecode","Phenotype","Total events","Person-years","Events","IRR (95%CI)","P-value","Events","IRR (95%CI)","P-value","Events","IRR (95%CI)","P-value","Events","IRR (95%CI)","P-value")) %>% 
  kable_classic("striped") %>%
  add_header_above(c(" " = 4, "A" = 3,"B" = 3,"AB" = 3,"Rhesus" = 3)) %>%
  pack_rows("Infectious Diseases", 1, 5, label_row_css = "background-color: #2FCFD3; color: #fff;") %>%
  pack_rows("Neoplasms", 6, 13, label_row_css = "background-color: #9E4F46; color: #fff;") %>%
  pack_rows("Endocrine/Metabolic", 14, 22, label_row_css = "background-color: #FF1AB9; color: #fff;") %>%
  pack_rows("Hematopoietic", 23, 28, label_row_css = "background-color: #F69EDC; color: #fff;") %>%
  pack_rows("Mental Disorders", 29, 29, label_row_css = "background-color: #2DC92D; color: #fff;") %>%
  pack_rows("Neurological", 30, 33, label_row_css = "background-color: #004DE6; color: #fff;") %>%
  pack_rows("Sense Organs", 34, 41, label_row_css = "background-color: #FFC184; color: #fff;") %>%
  pack_rows("Circulatory System", 42, 64, label_row_css = "background-color: #A00CC4; color: #fff;") %>%
  pack_rows("Respiratory", 65, 67, label_row_css = "background-color: #625D5D; color: #fff;") %>%
  pack_rows("Digestive", 68, 82, label_row_css = "background-color: #008495; color: #fff;") %>%
  pack_rows("Genitourinary", 83, 87, label_row_css = "background-color: #F9DF1E; color: #fff;") %>%
  pack_rows("Pregnancy Complications", 88, 94, label_row_css = "background-color: #17CE95; color: #fff;") %>%
  pack_rows("Dermatologic", 95, 96, label_row_css = "background-color: #729656; color: #fff;") %>%
  pack_rows("Musculoskeletal", 97, 103, label_row_css = "background-color: #960866; color: #fff;") %>%
  pack_rows("Congenital Anomalies", 104, 106, label_row_css = "background-color: #12601B; color: #fff;") %>%
  #pack_rows("Symptoms", 159, 163, label_row_css = "background-color: #00007B; color: #fff;") %>%
  #pack_rows("Injuries & Poisonings", 72, 73, label_row_css = "background-color: #FF0000; color: #fff;") %>%
  #pack_rows("Other", 74, 75, label_row_css = "background-color: #8dd3c7; color: #fff;") %>%
  save_kable(file = "results/20220915/referenceO/enter_registry/table2.html", self_contained = T)  


## Table of all phecode for suplementary

data = read_tsv("results/20220915/referenceO/enter_registry/phewas_estimates.tsv")
data$phecode = as.character(data$phecode)

# Select only blood type info
data <- data %>% filter(term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv"))


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


# Make dataframe in order for kableextra with blood type as columns
data <- data %>% arrange(category_number, as.numeric(phecode)) %>%
  select(category,phecode,cases,cases_A,cases_B,cases_AB,cases_Rh,pyrs,phenotype,term,IRR_CI,p.value_FDR) %>% 
  mutate(p.value_FDR = if_else(p.value_FDR<0.001,"<0.001",as.character(round(p.value_FDR,3)))) %>%
  #mutate(p.value_FDR = if_else(p.value_FDR<0.001,formatC(p.value_FDR,format="e",digits=2),as.character(round(p.value_FDR,3)))) %>%
  pivot_wider(names_from = "term",values_from = c("IRR_CI","p.value_FDR")) %>%
  select("category","phecode","phenotype","cases","pyrs",
    "cases_A","IRR_CI_factor(AB0)A","p.value_FDR_factor(AB0)A",
    "cases_B","IRR_CI_factor(AB0)B","p.value_FDR_factor(AB0)B",
    "cases_AB","IRR_CI_factor(AB0)AB","p.value_FDR_factor(AB0)AB",
    "cases_Rh","IRR_CI_factor(Rhesus)Positiv","p.value_FDR_factor(Rhesus)Positiv")

# Rename columns not good with ( ) in column names
data <- data %>% rename(IRR_CI_A = "IRR_CI_factor(AB0)A",IRR_CI_B = "IRR_CI_factor(AB0)B",IRR_CI_AB = "IRR_CI_factor(AB0)AB",IRR_CI_Rhesus = "IRR_CI_factor(Rhesus)Positiv")
data <- data %>% rename(p.value_FDR_A = "p.value_FDR_factor(AB0)A",p.value_FDR_B = "p.value_FDR_factor(AB0)B",p.value_FDR_AB = "p.value_FDR_factor(AB0)AB",p.value_FDR_Rhesus = "p.value_FDR_factor(Rhesus)Positiv")


# Round person-years
data <- data %>% mutate(pyrs = as.character(round(pyrs,0)))

# # Put asterisk on prevalance ratio results
prevalence_phecodes3 = c("637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282")
prevalence_phecodes5 = c("244.5","612.3","286.1","282.9","356","362.7","520.1","282.8","334.1","364.5")

data <- data %>% 
  mutate(
    IRR_CI_A = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_A,"**",sep=""),IRR_CI_A),
    IRR_CI_B = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_B,"**",sep=""),IRR_CI_B),
    IRR_CI_AB = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_AB,"**",sep=""),IRR_CI_AB),
    IRR_CI_Rhesus = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(IRR_CI_Rhesus,"**",sep=""),IRR_CI_Rhesus),
    pyrs = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste(pyrs,"*",sep=""),pyrs),
    )

# Print category order and names
for (category in unique(data$category)) {
   print(category)
   print(which(data$category %in% category))
}

data <- data %>% select(-category)

# Color bold significant hits
data <- data %>%
  mutate(
    IRR_CI_A = cell_spec(IRR_CI_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    IRR_CI_B = cell_spec(IRR_CI_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    IRR_CI_AB = cell_spec(IRR_CI_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    IRR_CI_Rhesus = cell_spec(IRR_CI_Rhesus,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_Rhesus,"<"))<0.05) | (p.value_FDR_Rhesus =="<0.001"),T,F)),
    p.value_FDR_A = cell_spec(p.value_FDR_A,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_A,"<"))<0.05) | (p.value_FDR_A =="<0.001"),T,F)),
    p.value_FDR_B = cell_spec(p.value_FDR_B,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_B,"<"))<0.05) | (p.value_FDR_B =="<0.001"),T,F)),
    p.value_FDR_AB = cell_spec(p.value_FDR_AB,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_AB,"<"))<0.05) | (p.value_FDR_AB =="<0.001"),T,F)),
    p.value_FDR_Rhesus = cell_spec(p.value_FDR_Rhesus,"html",bold = ifelse((as.numeric(str_remove(p.value_FDR_Rhesus,"<"))<0.05) | (p.value_FDR_Rhesus =="<0.001"),T,F)))


kbl(data,"html",escape = F,col.names=c("Phecode","Phenotype","Total events","Person-years","Events","IRR (95%CI)","P-value","Events","IRR (95%CI)","P-value","Events","IRR (95%CI)","P-value","Events","IRR (95%CI)","P-value")) %>% 
  kable_classic("striped") %>%
  add_header_above(c(" " = 4, "A" = 3,"B" = 3,"AB" = 3,"Rhesus" = 3)) %>%
  pack_rows("Infectious Diseases", 1, 54, label_row_css = "background-color: #2FCFD3; color: #fff;") %>%
  pack_rows("Neoplasms", 55, 170, label_row_css = "background-color: #9E4F46; color: #fff;") %>%
  pack_rows("Endocrine/Metabolic", 171, 287, label_row_css = "background-color: #FF1AB9; color: #fff;") %>%
  pack_rows("Hematopoietic", 288, 329, label_row_css = "background-color: #F69EDC; color: #fff;") %>%
  pack_rows("Mental Disorders", 330, 393, label_row_css = "background-color: #2DC92D; color: #fff;") %>%
  pack_rows("Neurological", 394, 461, label_row_css = "background-color: #004DE6; color: #fff;") %>%
  pack_rows("Sense Organs", 462, 564, label_row_css = "background-color: #FFC184; color: #fff;") %>%
  pack_rows("Circulatory System", 565, 701, label_row_css = "background-color: #A00CC4; color: #fff;") %>%
  pack_rows("Respiratory", 702, 769, label_row_css = "background-color: #625D5D; color: #fff;") %>%
  pack_rows("Digestive", 770, 903, label_row_css = "background-color: #008495; color: #fff;") %>%
  pack_rows("Genitourinary", 904, 1031, label_row_css = "background-color: #F9DF1E; color: #fff;") %>%
  pack_rows("Pregnancy Complications", 1032, 1084, label_row_css = "background-color: #17CE95; color: #fff;") %>%
  pack_rows("Dermatologic", 1085, 1162, label_row_css = "background-color: #729656; color: #fff;") %>%
  pack_rows("Musculoskeletal", 1163, 1260, label_row_css = "background-color: #960866; color: #fff;") %>%
  pack_rows("Congenital Anomalies", 1261, 1312, label_row_css = "background-color: #12601B; color: #fff;") %>%
  #pack_rows("Symptoms", 1313, 1342, label_row_css = "background-color: #00007B; color: #fff;") %>%
  save_kable(file = "results/20220915/referenceO/enter_registry/suplement_table3.html", self_contained = T) 

