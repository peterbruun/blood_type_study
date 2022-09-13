# Manhattan, Volcano, and Forest Plots.
# Run in R terminal

library(tidyverse)
library(broom)
library(ggrepel)
library(PheWAS)
library(RColorBrewer)
library(ggpubr)

#setwd("/users/people/petras/bloodtype_diseasesV2/")
setwd("/users/secureome/home/people/petras/bloodtype_diseasesV2/")

# Disease color 
category = c(
	"Other",
	"infectious diseases",
	"neoplasms",
	"endocrine/metabolic",
	"hematopoietic",
	"mental disorders",
	"neurological",
	"sense organs",
	"circulatory system",
	"respiratory",
	"digestive",
	"genitourinary",
	"pregnancy complications",
	"dermatologic",
	"musculoskeletal",
	"congenital anomalies",
	"symptoms",
	"injuries & poisonings")

brunak_palette = c(
"Other" = "#8dd3c7",
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
"pregnancy" = "#17CE95",#12601B"
"dermatologic" = "#729656",
"musculoskeletal" = "#960866",
"congenital" = "#12601B",
#"chapterColors17" = "#F68409",
"symptoms" = "#00007B",
"injuries" = "#FF0000")
#"chapterColors20" = "#B7C94B",
#"chapterColors21" = "#000000")
	
mapping = data.frame(category,brunak_palette)


#data = read_tsv(snakemake@input[["phewas"]])
#data = read_tsv("results/20220531/quasi_spline_poisson/enter_registry/phewas_estimates.tsv")
data = read_tsv("results/20220531/quasi_spline_poisson/enter_registry/phewas_estimates.tsv")


# Map colour code onto
data[is.na(data$category),"category"] = "Other"
data <- left_join(data,mapping,by="category")

# Color all non-significant Phecodes grey
data <- data %>% mutate(brunak_palette = if_else(p.value_FDR>0.05,"#D3D3D3",as.character(brunak_palette)))

#all_estimates[all_estimates.term == "factor(AB0)A"]
data$phecode = as.character(data$phecode)

data %>% select(brunak_palette) %>% distinct() %>% count()

# testing
#trait= "factor(AB0)A"
#bloodtype = data[data$term == trait,]

# In GGPLOT
for (trait in c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv")) {
	bloodtype = data[data$term == trait,]

	# Volcano plot
	pdf(file=paste("results/20220531/quasi_spline_poisson/enter_registry/Volcano_",trait,".pdf",sep = ""), height = 8)
	print(bloodtype %>% ggplot(aes(x=log2(estimate),y=-log10(p.value_FDR))) + 
		geom_point(aes(color=FDR_.05)) +
		geom_text_repel(data = .%>% filter(FDR_.05==TRUE),aes(label=phecode)) +
		labs(title=paste("PheWAS: Bloodtype ",trait," (ref O)",sep="")) +
		xlab("log2(Hazard-ratio)") +
		ylab("-log10(FDR corrected p-value") +
		theme(text = element_text(size=10)))
	dev.off()

	#bloodtype = data[data$term == "factor(AB0)A",]

	# PheWAS plot
	bloodtype$description = bloodtype$phenotype
	bloodtype$phenotype = bloodtype$phecode
	bloodtype$p = bloodtype$p.value_FDR
	# Remove log 10
	#bloodtype$p = 10**(bloodtype$p.value_FDR)
	bloodtype$OR = bloodtype$estimate
	bloodtype$group = bloodtype$category
	#bloodtype$size = abs(bloodtype$estimate-1)
	bloodtype$size = bloodtype$estimate
	#bloodtype$direction = bloodtype$estimate
	bloodtype$color = as.character(bloodtype$brunak_palette)
	bloodtype$groupnum = bloodtype$category_number

	#palette = get_palette(palette = "default", 19)
	# Help found here:
	#https://github.com/PheWAS/PheWAS/blob/master/R/phenotypePlot.R
	# https://github.com/PheWAS/PheWAS/tree/master/R
	# https://github.com/PheWAS/PheWAS/blob/master/R/phenotypeManhattan.R
	#color.palette=brunak_palette

	bloodtype <- bloodtype %>% arrange(phecode)
	pdf(file=paste("results/20220531/quasi_spline_poisson/enter_registry/PheWAS_",trait,".pdf",sep = ""),width=13,height=7)
	if (trait == "factor(Rhesus)Positiv"){
		title = "Blood Group RhD (Reference: RhD Negative)"
	} else {
	title = paste("Blood Group ",strsplit(trait,")")[[1]][2]," (Reference: Blood Group O)",sep="")
	}
	print(phenotypeManhattan(bloodtype,title=title,
		significant.line=0.05,x.axis.label="PheWAS Diagnosis Code",
		annotate.phenotype.description=F,OR.direction=T,
		sort.by.category.value=T,use.color=T,
		annotate.phenotype=T,x.group.labels=T,
		sizes=F,OR.size=T,annotate.size=2.8,
		annotate.angle=0) + scale_size("Rate Ratio",range = c(2, 5),breaks=c(1,1.1,1.4,2,4)) +
		# Make y-axis log10-scale
		scale_y_log10(limits = c(1e-1,NA),name=expression(-log[10](italic(p)))) +
		#scale_y_continuous(trans='log2',name=expression(-log[10](italic(p)))) +
		annotation_logticks(sides="l") +
		#scale_x_continuous() +	
		#aes(reorder(group, p)) +
		# Remove direction legend
		guides(shape=FALSE))	
	dev.off()

}


# Forest Plot


## GGPLOT Forest plot ## 

# https://nightingalehealth.github.io/ggforestplot/articles/ggforestplot.html#installation
library(ggforestplot)
library(ggforce)


data = read_tsv("results/20220531/quasi_spline_poisson/enter_registry/phewas_estimates.tsv")
data$phecode = as.character(data$phecode)

# Select subset with significant phecodes

#significant_phecodes <- data %>% 
#	filter(term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv") & FDR_.05) %>% 
#	select(phecode) %>% unique() %>% pull()

significant_phecodes <- unique(data[data$term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv")&data$FDR_.05,]$phecode)
data <- data %>% filter(phecode %in% significant_phecodes & term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv"))


# Remove factor # Rename columns not good with ( ) in column names
data <- data %>%
  mutate(term = case_when(
    term == "factor(AB0)A"  ~ "Blood Group A",
    term == "factor(AB0)B"  ~ "Blood Group B",
    term == "factor(AB0)AB"  ~ "Blood Group AB",
    term == "factor(Rhesus)Positiv"  ~ "Blood Group RhD",
    TRUE ~ term
  ))

# Set order of term and phecodes
data[is.na(data$category),"category"] = "Other"
data <- data %>% 
  mutate(category = str_to_title(category))
group_order <- data %>% 
  arrange(as.numeric(phecode)) %>% 
  distinct(category) %>% 
  pull(category)
data <- data %>%
  arrange(as.numeric(phecode)) %>%
  mutate(Blood_type = factor(term,levels = c("Blood Group RhD","Blood Group AB","Blood Group B","Blood Group A")),category = factor(category, levels = group_order),name=phenotype)

      #
# Get estimate back to log-scale
data$estimate = log(data$estimate)
data[data$name=="Abnormality of organs and soft tissues of pelvis complicating pregnancy, childbirth, or the puerperium","name"] = "Abnormality of organs/soft tissues of pelvis complicating pregnancy"
# Every other
data$name = gsub("( \\S+) ", "\\1\n", data$name)

# Set * infront of phecodes that are prevalence ratios 

prevalence_phecodes3 = c("637","658","661","657","656","747","748","749","750","751","752","753","754","755","756","757","758","759","665","282")
prevalence_phecodes5 = c("244.5","612.3","286.1","282.9","356","362.7","520.1","282.8","334.1","364.5")

data <- data %>% 
  mutate(
     name = ifelse(substr(phecode, 1, 3) %in% prevalence_phecodes3 | substr(phecode, 1, 5) %in% prevalence_phecodes5,paste("*",name,sep=""),name))


# Save figure
pdf(file="results/20220531/quasi_spline_poisson/enter_registry/forestplot_phewas.pdf",width=8,height=60)#25
print(forestplot(
  df = data,
  estimate = estimate,
  se = std.error,
  pvalue = p.value_FDR,
  psignif = 0.05,
  xlab = "Incidence Rate Ratio (95% CI)",
  #title = "Association of AB0 blood group, Blood Group RhD type and\n diagnosis rate ratios",
  colour = Blood_type,
  logodds = TRUE) +
scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L),
    labels = c("Blood Group RhD","Blood Group AB","Blood Group B","Blood Group A")) +
facet_col(
    facets = ~category,
    scales = "free_y",
    space = "free") +
#theme_classic() +
theme(
  legend.position = "top",
  legend.title=element_blank(),
  axis.text=element_text(size=6),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background=element_blank()))
dev.off()

# #Plot figure without bone marrow, von will and isoimmuniation

# Remove bone marrow
data <- data %>% filter(phecode!="860")

pdf(file="results/20220531/quasi_spline_poisson/enter_registry/forestplot_phewas_not_bone.pdf",width=8,height=60)
print(forestplot(
  df = data,
  estimate = estimate,
  se = std.error,
  pvalue = p.value_FDR,
  psignif = 0.05,
  xlab = "Incidence Rate Ratio (95% CI)",
  #title = "Association of AB0 blood group, Blood Group RhD type and\n diagnosis rate ratios",
  colour = Blood_type,
  logodds = TRUE) +
scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L),
    labels = c("Blood Group RhD","Blood Group AB","Blood Group B","Blood Group A")) +
facet_col(
    facets = ~category,
    scales = "free_y",
    space = "free") +
#theme_classic() +
theme(
  legend.position = "top",
  legend.title=element_blank(),
  axis.text=element_text(size=6),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background=element_blank()))
dev.off()


## Plot significant hits with higheste IRRs
# https://datascienceplus.com/lattice-like-forest-plot-using-ggplot2-in-r/


data = read_tsv("results/20220531/quasi_spline_poisson/enter_registry/phewas_estimates.tsv")
data$phecode = as.character(data$phecode)

# Map colour code onto
data[is.na(data$category),"category"] = "Other"
data <- left_join(data,mapping,by="category")

# Select subset with significant phecodes

#significant_phecodes <- data %>% 
#	filter(term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv") & FDR_.05) %>% 
#	select(phecode) %>% unique() %>% pull()

data <- data %>% filter(data$FDR_.05 & term %in% c("factor(AB0)A","factor(AB0)B","factor(AB0)AB","factor(Rhesus)Positiv"))


# Remove factor # Rename columns not good with ( ) in column names
data <- data %>%
  mutate(term = case_when(
    term == "factor(AB0)A"  ~ "Blood Group A",
    term == "factor(AB0)B"  ~ "Blood Group B",
    term == "factor(AB0)AB"  ~ "Blood Group AB",
    term == "factor(Rhesus)Positiv"  ~ "Blood Group RhD",
    TRUE ~ term
  ))

# Set order of term and phecodes
data[is.na(data$category),"category"] = "Other"
data <- data %>% 
  mutate(category = str_to_title(category))
group_order <- data %>% 
  arrange(as.numeric(phecode)) %>% 
  distinct(category) %>% 
  pull(category)
data <- data %>%
  arrange(as.numeric(phecode)) %>%
  mutate(Blood_type = factor(term,levels = c("Blood Group RhD","Blood Group AB","Blood Group B","Blood Group A")),category = factor(category, levels = group_order),name=phenotype)

      #
# Get estimate back to log-scale
#data$estimate = log(data$estimate)
data[data$name=="Abnormality of organs and soft tissues of pelvis complicating pregnancy, childbirth, or the puerperium","name"] = "Abnormality of organs/soft tissues of pelvis complicating pregnancy"
# Every other
data$name = gsub("( \\S+) ", "\\1\n", data$name)

# Construct confidence intervals
data <- data %>% mutate(
	CI_low = estimate/exp(qnorm(0.975)*std.error),
 	CI_high = estimate*exp(qnorm(0.975)*std.error),
 	IRR = estimate)

# Sort on associations with biggest effect
inverse <- data %>% arrange(IRR)
positive <- data %>% arrange(desc(IRR))


p = ggplot(data=positive[1:20,],
    aes(x = reorder(name,IRR),y = IRR, ymin = CI_low, ymax = CI_high ))+
    #geom_point(aes(shape = term,col=category), size = 4) +
    geom_pointrange(aes(col=term))+
    geom_hline(aes(fill=term),yintercept =1, linetype=2)+
    xlab('term')+ ylab("Incidence Rate Ratio (95% Confidence Interval)") +
    labs(title = "Top 20 Postive Blood Group Associations") +
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high,col=term),width=0.5,cex=1)+ 
    #facet_wrap(~category,strip.position="left",nrow=9,scales = "free_y") +
    scale_colour_discrete("Blood Group") +
    theme(plot.title=element_text(size=16),#,face="bold"),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.title.y=element_blank(),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    coord_flip()
 p
ggsave("results/20220531/quasi_spline_poisson/enter_registry/top20_positive_IRRs.pdf",height=8)

p = ggplot(data=inverse[1:20,],
    aes(x = reorder(name,-IRR),y = IRR, ymin = CI_low, ymax = CI_high ))+
    #geom_point(aes(shape = term,col=category), size = 4) +
    geom_pointrange(aes(col=term))+
    geom_hline(aes(fill=term),yintercept =1, linetype=2)+
    xlab('term')+ ylab("Incidence Rate Ratio (95% Confidence Interval)") +
    labs(title = "Top 20 Inverse Blood Group associations") +
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high,col=term),width=0.5,cex=1)+ 
    #facet_wrap(~category,strip.position="left",nrow=9,scales = "free_y") +
    scale_colour_discrete("Blood Group") +
    theme(plot.title=element_text(size=16),#,face="bold"),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.title.y=element_blank(),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    coord_flip()
 p
ggsave("results/20220531/quasi_spline_poisson/enter_registry/top20_inverse_IRRs.pdf",height=8)

