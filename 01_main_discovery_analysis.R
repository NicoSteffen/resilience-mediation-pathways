# Load Packages -----------------------------------------------------------

library(foreign)
library(tidyverse)
library(easystats)
library(lme4)
library(lmerTest)
library(dendextend)
library(lavaan)
library(tidyr)
library(dplyr)
library(moments)
library(psych)
library(broom)

options(scipen = 999)

# Create custom functions -------------------------------------------------

roundallnumerics = function(df, digits){
  for(j in 1:ncol(df)){
    if(is.numeric(df[,j])){
      df[,j] = round(df[,j], digits)
    }
  }
  return(df)
}

p_labeller = function(vec){
  vec = as.numeric(vec)
  for(i in 1:length(vec)){
    if(is.na(vec[i]) == F & vec[i] < .001){
      vec[i] = "<.001***"
    }
    if(is.na(vec[i]) == F & vec[i] >= .001 & vec[i] < .01){
      vec[i] = paste0(vec[i], "**")
    }
    if(is.na(vec[i]) == F & vec[i] > .01 & vec[i] < .05){
      vec[i] = paste0(vec[i], "*")
    }
  }
  return(vec)
}

# Read data ---------------------------------------------------------------

# this dataset has multiple time points
swdata = foreign::read.spss()

# Clean data for patients -------------------------------------

container_TAQ_43 = swdata$TAQ_T0_43 # this item is merely a checkup item (not included in any subscales)

# Missings CDRISC entfernen

# Relevante Spalten extrahieren (alle CDRISC_T0 Items)
cdrisc_cols <- grep("CDRISC_T0", names(swdata), value = TRUE)

# Anzahl fehlender Werte pro Person berechnen
swdata$missing_CDRISC <- rowSums(is.na(swdata[, cdrisc_cols]))
table(swdata$missing_CDRISC)
# IDs mit mehr als 3 fehlenden Werten entfernen
swdata <- swdata[swdata$missing_CDRISC <= 3, ]


swdata = rename(swdata, c(ID = "Code",
                          Gender = "Geschlecht",
                          Age = "Alter",
                          Weight = "Gewicht_kg",
                          Height = "Größe_cm",
                          Comorbidities = "Anzahl_Diagnosen",
                          Age_onset = "Alter_bei_Ersterkrankung",
                          Hospitalizations = "Anzahl_Klinikaufenthalte",
                          Medication = "Aufnahne_Medikation_ja_nein",
                          Antidepressants = "Aufnahme_Antidepressiva",
                          Antipsychotics = "Aufnahme_Antipsychotika",
                          Benzodiazepines = "Aufnahme_Benzodiazepine"))

levels(swdata$Gender) = c("female", "male")
levels(swdata$Antidepressants) = c("yes", "no")
levels(swdata$Benzodiazepines) = c("yes", "no")
levels(swdata$Medication) = c("yes", "no")
levels(swdata$Antipsychotics) = c("yes", "no")

swdata$Antidepressants = relevel(swdata$Antidepressants, "no")
swdata$Benzodiazepines = relevel(swdata$Benzodiazepines, "no")
swdata$Medication = relevel(swdata$Medication, "no")
swdata$Antipsychotics = relevel(swdata$Antipsychotics, "no")

swdata = BBmisc::dropNamed(swdata, drop = c("filter_.", 
                                            "Uhrzeit", 
                                            "Bemerkungen",
                                            "Datum",
                                            "TAQ_T0_43"))

names(swdata) = gsub("TAQ_T0", "TAQ", names(swdata))


# Fetch relevant sample characteristics

characteristics = c("Gender",
                    "Age",
                    "Comorbidities",
                    "Age_onset",
                    "Hospitalizations",
                    "Medication",
                    "Antidepressants",
                    "Antipsychotics",
                    "Benzodiazepines")

# Select patient samples

table(swdata$Hauptdiagnose_ICD10_Kategorie)
table(swdata$Diagnose_ICD10_Welche)
table(swdata$Hauptdiagnose_Detail)

swdata$Diagnosis = NA
swdata$Diagnosis[swdata$Hauptdiagnose_Detail == "Bipolar"] = "Bipolar"
swdata$Diagnosis[swdata$Hauptdiagnose_ICD10_Kategorie == "Schizophrenie, schizotype und wahnhafte Störungen (F20-F29)"] = "Psychotic"
swdata$Diagnosis[swdata$Hauptdiagnose_Detail != "Bipolar" & swdata$Hauptdiagnose_ICD10_Kategorie == "Affektive Störungen (F30-F39)"] = "Depression"
swdata$Diagnosis[swdata$Hauptdiagnose_ICD10_Kategorie == "Persönlichkeits- und Verhaltensstörungen (F60-F69)"] = "PD"
swdata$Diagnosis[swdata$Hauptdiagnose_ICD10_Kategorie == "Neurotische, Belastungs- und somatoforme Störungen (F40-F49)"] = "Neurotic, stress-related and somatoform"

#remove unnecessary variables
patterns <- c("MSWS", "RS25", "SES", "SB", "BDI_T2", "BDI_T4", "BDI_T6", "BDI_T8")
swdata <- swdata[, !sapply(patterns, function(p) grepl(p, names(swdata))) %>% apply(1, any)]

# 2 TAQ Items "Family Secrets" wurden aus späteren Versionen ausgeschlossen
# https://complextrauma.org/wp-content/uploads/2019/03/TAQ-2019.pdf

# CAVE: In our (German) version there are only 3 age categories (7-12; 13-18; adult). The english version had 4 (0-6; 7-12; 13-18; adult)

swdata = swdata[, -c(which(grepl("TAQ_9", names(swdata))))]
swdata = swdata[, -c(which(grepl("TAQ_24", names(swdata))))]



# renumber TAQ items after 2-item exclusion

names(swdata)[grep("TAQ", names(swdata))] = paste0("TAQ_", rep(1:40, each = 3), rep(letters[1:3], 40))


# Explore data ------------------------------------------------------------

report_sample(swdata[, c(characteristics, "Diagnosis")], group_by = "Diagnosis")

swdata_f20 = swdata[swdata$Hauptdiagnose_ICD10_Kategorie == "Schizophrenie, schizotype und wahnhafte Störungen (F20-F29)",]
swdata_dep = swdata[swdata$Hauptdiagnose_ICD10_Kategorie == "Affektive Störungen (F30-F39)",]
swdata_dep = swdata_dep[swdata_dep$Hauptdiagnose_Detail != "Bipolar",]
swdata_dep = swdata_dep[swdata_dep$Hauptdiagnose_Detail != "Dystymie",]

swdata = swdata_dep 

levels(swdata$Gender) = c("Female", "Male")

flextable::save_as_docx(flextable::flextable(report_sample(swdata[, c(characteristics)])), n = T, path = "Table1.docx")

table(swdata$Benzodiazepines)
table(swdata$Hauptdiagnose_Detail)


# create data frame for alpha calculation (without recoding 1 to 0)
swdata2 = swdata

# convert item to numeric values
# recode 2nd TAQ item

for(j in which(grepl("TAQ_2", names(swdata)))){
  swdata[,j] = as.character(swdata[,j])
  swdata[,j][swdata[,j] == "Oft?/sehr stark"] = "0"
  swdata[,j][swdata[,j] == "Gelegentlich?/im Grunde schon"] = "0"
  swdata[,j][swdata[,j] == "Selten?/ein wenig"] = "2"
  swdata[,j][swdata[,j] == "Nie?/gar nicht"] = "3"
  swdata[,j][swdata[,j] == "Weiß nicht"] = "0"
  swdata[,j] = as.numeric(swdata[,j])
}


#Achtung: changed coding of TAQ according to scoring recommendations:


for(j in which(grepl("TAQ", names(swdata)))){
  swdata[,j] = as.character(swdata[,j])
  swdata[,j][swdata[,j] == "Oft?/sehr stark"] = "3"
  swdata[,j][swdata[,j] == "Gelegentlich?/im Grunde schon"] = "2"
  swdata[,j][swdata[,j] == "Selten?/ein wenig"] = "0"
  swdata[,j][swdata[,j] == "Nie?/gar nicht"] = "0"
  swdata[,j][swdata[,j] == "nein"] = "0"
  swdata[,j][swdata[,j] == "ja"] = "3"
  swdata[,j][swdata[,j] == "Weiß nicht"] = "0"
  swdata[,j] = as.numeric(swdata[,j])
}



for(j in which(grepl("CDRISC_T0", names(swdata)))){
  swdata[,j] = as.character(swdata[,j])
  swdata[,j][swdata[,j] == "überhaupt nicht wahr"] = "0"
  swdata[,j][swdata[,j] == "selten wahr"] = "1"
  swdata[,j][swdata[,j] == "manchmal wahr"] = "2"
  swdata[,j][swdata[,j] == "oft wahr"] = "3"
  swdata[,j][swdata[,j] == "fast immer wahr"] = "4"
  swdata[,j] = as.numeric(swdata[,j])
}


#Mean Imputation CDRISC
cdrisc_mat <- swdata[, cdrisc_cols]
row_means <- rowMeans(cdrisc_mat, na.rm = TRUE)

# Schleife oder apply zum Ersetzen der Missings
for (i in 1:nrow(cdrisc_mat)) {
  na_idx <- is.na(cdrisc_mat[i, ])
  cdrisc_mat[i, na_idx] <- row_means[i]
}

# Ersetzte Matrix wieder in den Datensatz einfügen
swdata[, cdrisc_cols] <- cdrisc_mat


swdata$CDRISC_25_T0_sum = rowSums(swdata[, grep("CDRISC_T0", names(swdata))])

# CDRISC 10 für Replication 

cdrisc10_items <- c(
  "CDRISC_T0_1", "CDRISC_T0_4", "CDRISC_T0_6", "CDRISC_T0_7", "CDRISC_T0_8",
  "CDRISC_T0_11", "CDRISC_T0_14", "CDRISC_T0_16", "CDRISC_T0_17", "CDRISC_T0_19"
)

swdata$CDRISC_10_sum <- rowSums(swdata[, cdrisc10_items], na.rm = TRUE)


# create TAQ subscales according to https://complextrauma.org/wp-content/uploads/2019/03/TAQ-Manual-Spinazzola-2019.pdf

# Hilfsfunktion zur Berechnung der Subskalenwerte
calc_subscale <- function(data, items) {
  apply(data[, items], 1, function(row) {
    if (all(is.na(row))) {
      return(NaN)
    } else if (all(is.na(row) | row == 0)) {
      return(0)
    } else {
      return(mean(ifelse(row == 0, NA, row), na.rm = TRUE))
    }
  })
}

# Liste mit Subskalen-Namen und den zugehörigen Items pro Altersgruppe
subscales <- list(
  TAQ_Competence = list(a = c("TAQ_3a", "TAQ_4a"), b = c("TAQ_3b", "TAQ_4b"), c = c("TAQ_3c", "TAQ_4c")),
  TAQ_Safety = list(a = c("TAQ_1a", "TAQ_5a", "TAQ_8a"), b = c("TAQ_1b", "TAQ_5b", "TAQ_8b"), c = c("TAQ_1c", "TAQ_5c", "TAQ_8c")),
  TAQ_Neglect = list(a = c("TAQ_2a", "TAQ_6a", "TAQ_7a", "TAQ_21a", "TAQ_27a"),
                     b = c("TAQ_2b", "TAQ_6b", "TAQ_7b", "TAQ_21b", "TAQ_27b"),
                     c = c("TAQ_2c", "TAQ_6c", "TAQ_7c", "TAQ_21c", "TAQ_27c")),
  TAQ_Separation = list(a = c("TAQ_10a", "TAQ_11a", "TAQ_12a", "TAQ_14a"),
                        b = c("TAQ_10b", "TAQ_11b", "TAQ_12b", "TAQ_14b"),
                        c = c("TAQ_10c", "TAQ_11c", "TAQ_12c", "TAQ_14c")),
  TAQ_Emotional_Abuse = list(a = c("TAQ_9a", "TAQ_16a", "TAQ_17a", "TAQ_18a", "TAQ_19a"),
                             b = c("TAQ_9b", "TAQ_16b", "TAQ_17b", "TAQ_18b", "TAQ_19b"),
                             c = c("TAQ_9c", "TAQ_16c", "TAQ_17c", "TAQ_18c", "TAQ_19c")),
  TAQ_Physical_Abuse = list(a = c("TAQ_28a", "TAQ_29a", "TAQ_30a"),
                            b = c("TAQ_28b", "TAQ_29b", "TAQ_30b"),
                            c = c("TAQ_28c", "TAQ_29c", "TAQ_30c")),
  TAQ_Sexual_Abuse = list(a = c("TAQ_35a", "TAQ_36a", "TAQ_37a", "TAQ_38a"),
                          b = c("TAQ_35b", "TAQ_36b", "TAQ_37b", "TAQ_38b"),
                          c = c("TAQ_35c", "TAQ_36c", "TAQ_37c", "TAQ_38c")),
  TAQ_Witnessing = list(a = c("TAQ_20a", "TAQ_22a", "TAQ_23a", "TAQ_24a", "TAQ_31a", "TAQ_34a"),
                        b = c("TAQ_20b", "TAQ_22b", "TAQ_23b", "TAQ_24b", "TAQ_31b", "TAQ_34b"),
                        c = c("TAQ_20c", "TAQ_22c", "TAQ_23c", "TAQ_24c", "TAQ_31c", "TAQ_34c")),
  TAQ_Impersonal_Trauma = list(a = c("TAQ_13a", "TAQ_15a", "TAQ_32a", "TAQ_33a", "TAQ_39a", "TAQ_40a"),
                               b = c("TAQ_13b", "TAQ_15b", "TAQ_32b", "TAQ_33b", "TAQ_39b", "TAQ_40b"),
                               c = c("TAQ_13c", "TAQ_15c", "TAQ_32c", "TAQ_33c", "TAQ_39c", "TAQ_40c")),
  TAQ_Alcohol_Drugs = list(a = c("TAQ_25a", "TAQ_26a"), b = c("TAQ_25b", "TAQ_26b"), c = c("TAQ_25c", "TAQ_26c"))
)

# Schleife zur Berechnung aller Subskalen
for (subscale in names(subscales)) {
  for (period in names(subscales[[subscale]])) {
    varname <- paste0(subscale, "_", period)
    swdata[[varname]] <- calc_subscale(swdata, subscales[[subscale]][[period]])
  }
}


swdata$TAQ_Competence_total = rowSums(swdata[, c("TAQ_Competence_a", "TAQ_Competence_b", "TAQ_Competence_c" )], na.rm = TRUE)
swdata$TAQ_Safety_total = rowSums(swdata[, c("TAQ_Safety_a", "TAQ_Safety_b", "TAQ_Safety_c")], na.rm = TRUE)
swdata$TAQ_Neglect_total = rowSums(swdata[, c("TAQ_Neglect_a", "TAQ_Neglect_b", "TAQ_Neglect_c")], na.rm = TRUE)
swdata$TAQ_Separation_total = rowSums(swdata[, c("TAQ_Separation_a", "TAQ_Separation_b", "TAQ_Separation_c")], na.rm = TRUE)
swdata$TAQ_Emotional_Abuse_total = rowSums(swdata[, c("TAQ_Emotional_Abuse_a", "TAQ_Emotional_Abuse_b", "TAQ_Emotional_Abuse_c")], na.rm = TRUE)
swdata$TAQ_Physical_Abuse_total = rowSums(swdata[, c("TAQ_Physical_Abuse_a", "TAQ_Physical_Abuse_b", "TAQ_Physical_Abuse_c")], na.rm = TRUE)
swdata$TAQ_Sexual_Abuse_total = rowSums(swdata[, c("TAQ_Sexual_Abuse_a", "TAQ_Sexual_Abuse_b", "TAQ_Sexual_Abuse_c")], na.rm = TRUE)
swdata$TAQ_Witnessing_total = rowSums(swdata[, c("TAQ_Witnessing_a", "TAQ_Witnessing_b", "TAQ_Witnessing_c")], na.rm = TRUE)
swdata$TAQ_Impersonal_Trauma_total = rowSums(swdata[, c("TAQ_Impersonal_Trauma_a", "TAQ_Impersonal_Trauma_b", "TAQ_Impersonal_Trauma_c")], na.rm = TRUE)
swdata$TAQ_Alcohol_Drugs_total = rowSums(swdata[, c("TAQ_Alcohol_Drugs_a", "TAQ_Alcohol_Drugs_b", "TAQ_Alcohol_Drugs_c")], na.rm = TRUE)


# Fetch Subscale names

TAQ_subscales = names(swdata[, grep("TAQ_[a-zA-Z]", names(swdata))])
TAQ_subscales = unique(gsub("TAQ_|_a|_b|_c", "", TAQ_subscales))


# Multiple Mediation with CDRisk ------------------------------------------

#create new TAQ_Pos + TAQ_Neg

swdata$TAQ_Pos = rowMeans(swdata[,c("TAQ_Safety_total","TAQ_Competence_total")])

#subscales pos
swdata$TAQ_Pos_a = rowMeans(swdata[,c("TAQ_Safety_a","TAQ_Competence_a")], na.rm = T)
swdata$TAQ_Pos_b = rowMeans(swdata[,c("TAQ_Safety_b","TAQ_Competence_b")], na.rm = T)
swdata$TAQ_Pos_c = rowMeans(swdata[,c("TAQ_Safety_c","TAQ_Competence_c")], na.rm = T)

swdata$TAQ_Neg = rowMeans(swdata[, c(
  "TAQ_Neglect_total",
  "TAQ_Separation_total",
  "TAQ_Emotional_Abuse_total",
  "TAQ_Physical_Abuse_total",
  "TAQ_Sexual_Abuse_total",
  "TAQ_Witnessing_total",
  "TAQ_Impersonal_Trauma_total",
  "TAQ_Alcohol_Drugs_total"
)])


#subscales neg

swdata$TAQ_Neg_a = rowMeans(swdata[, c(
  "TAQ_Neglect_a",
  "TAQ_Separation_a",
  "TAQ_Emotional_Abuse_a",
  "TAQ_Physical_Abuse_a",
  "TAQ_Sexual_Abuse_a",
  "TAQ_Witnessing_a",
  "TAQ_Impersonal_Trauma_a",
  "TAQ_Alcohol_Drugs_a"
)], na.rm = T)

swdata$TAQ_Neg_b = rowMeans(swdata[, c(
  "TAQ_Neglect_b",
  "TAQ_Separation_b",
  "TAQ_Emotional_Abuse_b",
  "TAQ_Physical_Abuse_b",
  "TAQ_Sexual_Abuse_b",
  "TAQ_Witnessing_b",
  "TAQ_Impersonal_Trauma_b",
  "TAQ_Alcohol_Drugs_b"
)], na.rm = T)

swdata$TAQ_Neg_c = rowMeans(swdata[, c(
  "TAQ_Neglect_c",
  "TAQ_Separation_c",
  "TAQ_Emotional_Abuse_c",
  "TAQ_Physical_Abuse_c",
  "TAQ_Sexual_Abuse_c",
  "TAQ_Witnessing_c",
  "TAQ_Impersonal_Trauma_c",
  "TAQ_Alcohol_Drugs_c"
)], na.rm = T)



# rep datensatz erstellen

vars <- c("TAQ_Pos", "TAQ_Neg", "CDRISC_10_sum", "BDI_T0_Gesamt")
main <- swdata[, names(swdata) %in% vars]


# Mediation total ---------------------------------------------------------



mediation.model <- "
CDRISC_25_T0_sum ~ a*TAQ_Pos + b*TAQ_Neg
BDI_T0_Gesamt ~ c*CDRISC_25_T0_sum
BDI_T0_Gesamt ~ d*TAQ_Pos + e*TAQ_Neg
ac := a*c
bc := b*c
total.pos := d + (a*c)
total.neg := e + (b*c)
"

fita <- sem(mediation.model, 
           data=swdata, 
           fixed.x = F, 
           std.ov = T, 
           missing = "fiml",
           se = "bootstrap", 
           bootstrap = 10000)

summary(fita)

#semTable::sem_tables(fita)

d = roundallnumerics(
  as.data.frame(
    parameterEstimates(fita)
  ), 3
)

term = paste(d$lhs, d$op, d$rhs)

term = gsub("BDI_T0_Gesamt", "BDI", term)
term = gsub("CDRISC_25_T0_sum", "CDRISC", term)
term = gsub("TAQ_Pos", "TAQ+", term)
term = gsub("TAQ_Neg", "TAQ-", term)
term = gsub("~~", "↔", term)
term = gsub("~", "←", term)
term = gsub("ac := a*c" , "TAQ+ indirect", term, fixed = T)
term = gsub("bc := b*c" , "TAQ- indirect", term, fixed = T)
term = gsub("total.pos := d+(a*c)" , "TAQ+ total", term, fixed = T)
term = gsub("total.neg := e+(b*c)" , "TAQ- total", term,  fixed = T)

d = BBmisc::dropNamed(d, drop = c("lhs", 
                                  "op",
                                  "rhs",
                                  "label"))

d = add_column(d, Term = term, .before = 1) 

d = d[-which(grepl("←1", d$Term)),]

d = d[-which(grepl("CDRISC ↔ CDRISC|BDI ↔ BDI", d$Term)),]

d = d[-c(6,8),]

names(d) = c("Term", "Point estimate", "Standard error", "z", "p-value", "CI lower", "CI upper")

d$`p-value` = p_labeller(d$`p-value`)

flextable::save_as_docx(flextable::flextable(d), path = "Table2.docx")
flextable::save_as_docx(flextable::flextable(as.data.frame(report_sample(swdata[, c(characteristics, "Diagnosis")], group_by = "Gender"))), path = "Table1.docx")







# Assemble table
table_summary <- swdata %>%
  reframe(
    Variable = c(
      "Major Depression (recurrent)",
      "Major Depression",
      #"Positive experience sum score [TAQ+]",
      #"Negative experience sum score [TAQ-]",
      #"Depressive symptoms sum score [BDI-II]",
      #"Resilience sum score [CD-RISC]",
      "Age of onset (years)",
      "Age (years)",
      "Comorbidities",
      "Hospitalizations",
      "Medication [yes], %",
      "Antidepressants [yes], %",
      "Antipsychotics [yes], %",
      "Benzodiazepines [yes], %"
    ),
    `All patients (n = 133)` = c(
      sum(Hauptdiagnose_Detail == "Rezidivierend Depressiv", na.rm = TRUE),
      sum(Hauptdiagnose_Detail == "Episodisch Depressiv (MD)", na.rm = TRUE),
      #mean(TAQ_Pos, na.rm = TRUE),
      #mean(TAQ_Neg, na.rm = TRUE),
      #mean(BDI_T0_Gesamt, na.rm = TRUE),
      #mean(CDRISC_25_T0_sum, na.rm = TRUE),
      mean(Age_onset, na.rm = TRUE),
      mean(Age, na.rm = TRUE),
      mean(Comorbidities, na.rm = TRUE),
      mean(Hospitalizations, na.rm = TRUE),
      mean(Medication == "yes", na.rm = TRUE) * 100,
      mean(Antidepressants == "yes", na.rm = TRUE) * 100,
      mean(Antipsychotics == "yes", na.rm = TRUE) * 100,
      mean(Benzodiazepines == "yes", na.rm = TRUE) * 100
    ),
    `Standard Deviation` = c(
      NA, NA,
      #sd(TAQ_Pos, na.rm = TRUE),
      #sd(TAQ_Neg, na.rm = TRUE),
      #sd(BDI_T0_Gesamt, na.rm = TRUE),
      #sd(CDRISC_25_T0_sum, na.rm = TRUE),
      sd(Age_onset, na.rm = TRUE),
      sd(Age, na.rm = TRUE),
      sd(Comorbidities, na.rm = TRUE),
      sd(Hospitalizations, na.rm = TRUE),
      rep(NA, 4)
    )#,
    #`Female (n = 75)` = c(
      #sum(Hauptdiagnose_Detail[Gender == "Female"] == "Rezidivierend Depressiv", na.rm = TRUE),
      #sum(Hauptdiagnose_Detail[Gender == "Female"] == "Episodisch Depressiv (MD)", na.rm = TRUE),
    #   mean(TAQ_Pos[Gender == "Female"], na.rm = TRUE),
    #   mean(TAQ_Neg[Gender == "Female"], na.rm = TRUE),
    #   mean(BDI_T0_Gesamt[Gender == "Female"], na.rm = TRUE),
    #   mean(CDRISC_25_T0_sum[Gender == "Female"], na.rm = TRUE),
    #   mean(Age_onset[Gender == "Female"], na.rm = TRUE),
    #   mean(Age[Gender == "Female"], na.rm = TRUE),
    #   mean(Comorbidities[Gender == "Female"], na.rm = TRUE),
    #   mean(Hospitalizations[Gender == "Female"], na.rm = TRUE),
    #   mean(Medication[Gender == "Female"] == "yes", na.rm = TRUE) * 100,
    #   mean(Antidepressants[Gender == "Female"] == "yes", na.rm = TRUE) * 100,
    #   mean(Antipsychotics[Gender == "Female"] == "yes", na.rm = TRUE) * 100,
    #   mean(Benzodiazepines[Gender == "Female"] == "yes", na.rm = TRUE) * 100
    # ),
    # `Female SD` = c(
    #   NA, NA,
    #   sd(TAQ_Pos[Gender == "Female"], na.rm = TRUE),
    #   sd(TAQ_Neg[Gender == "Female"], na.rm = TRUE),
    #   sd(BDI_T0_Gesamt[Gender == "Female"], na.rm = TRUE),
    #   sd(CDRISC_25_T0_sum[Gender == "Female"], na.rm = TRUE),
    #   sd(Age_onset[Gender == "Female"], na.rm = TRUE),
    #   sd(Age[Gender == "Female"], na.rm = TRUE),
    #   sd(Comorbidities[Gender == "Female"], na.rm = TRUE),
    #   sd(Hospitalizations[Gender == "Female"], na.rm = TRUE),
    #   rep(NA, 4)
    # ),
    # `Male (n = 58)` = c(
    #   sum(Hauptdiagnose_Detail[Gender == "Male"] == "Rezidivierend Depressiv", na.rm = TRUE),
    #   sum(Hauptdiagnose_Detail[Gender == "Male"] == "Episodisch Depressiv (MD)", na.rm = TRUE),
    #   mean(TAQ_Pos[Gender == "Male"], na.rm = TRUE),
    #   mean(TAQ_Neg[Gender == "Male"], na.rm = TRUE),
    #   mean(BDI_T0_Gesamt[Gender == "Male"], na.rm = TRUE),
    #   mean(CDRISC_25_T0_sum[Gender == "Male"], na.rm = TRUE),
    #   mean(Age_onset[Gender == "Male"], na.rm = TRUE),
    #   mean(Age[Gender == "Male"], na.rm = TRUE),
    #   mean(Comorbidities[Gender == "Male"], na.rm = TRUE),
    #   mean(Hospitalizations[Gender == "Male"], na.rm = TRUE),
    #   mean(Medication[Gender == "Male"] == "yes", na.rm = TRUE) * 100,
    #   mean(Antidepressants[Gender == "Male"] == "yes", na.rm = TRUE) * 100,
    #   mean(Antipsychotics[Gender == "Male"] == "yes", na.rm = TRUE) * 100,
    #   mean(Benzodiazepines[Gender == "Male"] == "yes", na.rm = TRUE) * 100
    # ),
    # `Male SD` = c(
    #   NA, NA,
    #   sd(TAQ_Pos[Gender == "Male"], na.rm = TRUE),
    #   sd(TAQ_Neg[Gender == "Male"], na.rm = TRUE),
    #   sd(BDI_T0_Gesamt[Gender == "Male"], na.rm = TRUE),
    #   sd(CDRISC_25_T0_sum[Gender == "Male"], na.rm = TRUE),
    #   sd(Age_onset[Gender == "Male"], na.rm = TRUE),
    #   sd(Age[Gender == "Male"], na.rm = TRUE),
    #   sd(Comorbidities[Gender == "Male"], na.rm = TRUE),
    #   sd(Hospitalizations[Gender == "Male"], na.rm = TRUE),
    #   rep(NA, 4)
    # ),
    # `p-value` = p_values
  )

# Optional: Round all numeric values
table_summary <- roundallnumerics(table_summary, 2)

# Print the final table
print(table_summary)



# Correlation analyses 

data <- swdata[, c("TAQ_Pos","TAQ_Pos_a","TAQ_Pos_b","TAQ_Pos_c", 
                   "TAQ_Neg", "TAQ_Neg_a", "TAQ_Neg_b", "TAQ_Neg_c", 
                   "CDRISC_25_T0_sum", "BDI_T0_Gesamt")]

descriptive_stats <- data.frame(
  Variable = c("TAQ_P","TAQ_PC", "TAQ_PY", "TAQ_PA", "TAQ_N", "TAQ_PC", "TAQ_PY", "TAQ_PA",
               "CDRISC", "BDI"),
  M = sapply(data, mean, na.rm = TRUE),
  SD = sapply(data, sd, na.rm = TRUE)#,
  #Skew = sapply(data, skewness, na.rm = TRUE),
  #Kurtosis = sapply(data, kurtosis, na.rm = TRUE)
)


# Datensatz mit nur transformierten TAQ Itemantworten erstellen
data_taq = swdata[, grep("^TAQ_[0-9]", colnames(swdata))]

#Erstelle Datensatz für alpha 
data_pos = swdata[,c("TAQ_3a", "TAQ_3b", "TAQ_3c", "TAQ_4a", "TAQ_4b", "TAQ_4c", 
                     "TAQ_1a", "TAQ_1b", "TAQ_1c", "TAQ_5a", "TAQ_5b", "TAQ_5c", 
                     "TAQ_8a", "TAQ_8b", "TAQ_8c")]

#Erstelle Datensatz für alpha 
data_neg <- data_taq[, !(colnames(data_taq) %in% colnames(data_pos))]


#Erstelle Datensatz für alpha 
data_cdrisc <- swdata[, grep("CDRISC_T0_", colnames(swdata))]

#Erstelle Datensatz für alpha 
data_bdi = swdata[, grep("BDI_T0", colnames(swdata))]
data_bdi = data_bdi %>% select(-BDI_T0_S_a, -BDI_T0_S_a_Aus, -BDI_T0_S_Aus, -BDI_T0_Gesamt)
data_bdi[] <- lapply(data_bdi, function(x) as.numeric(as.character(x)))


#berechne alpha
alpha(data_pos, use ="pairwise")$total$raw_alpha
alpha(data_neg, use ="pairwise")$total$raw_alpha
alpha(data_cdrisc, use ="pairwise")$total$raw_alpha
alpha(data_bdi, use ="pairwise")$total$raw_alpha


# Pearson correlations
correlation_matrix <- cor(data, use = "pairwise.complete.obs")
correlation_matrix <- as.data.frame(correlation_matrix)

# Combine all information into a final table
final_table <- descriptive_stats
final_table$TAQ_P <- correlation_matrix$TAQ_Pos
final_table$TAQ_PC <- correlation_matrix$TAQ_Pos_a
final_table$TAQ_PY <- correlation_matrix$TAQ_Pos_b
final_table$TAQ_PA <- correlation_matrix$TAQ_Pos_c
final_table$TAQ_N <- correlation_matrix$TAQ_Neg
final_table$TAQ_NC <- correlation_matrix$TAQ_Neg_a
final_table$TAQ_NY <- correlation_matrix$TAQ_Neg_b
final_table$TAQ_NA <- correlation_matrix$TAQ_Neg_c
final_table$CDRISC <- correlation_matrix$CDRISC_25_T0_sum
final_table$BDI <- correlation_matrix$BDI_T0_Gesamt

final_table = roundallnumerics(final_table,2)
flextable::save_as_docx(flextable::flextable(final_table), path = "Cor.Table.docx")

correlation_results <- Hmisc::rcorr(as.matrix(data))
p_values_matrix <- correlation_results$P
p_values_df <- as.data.frame(p_values_matrix)

#Benjamin Hochberg corection 
p_values_vector <- as.vector(p_values_matrix)
p_values_bh <- p.adjust(p_values_vector, method = "BH")

# Reshape the corrected p-values back into matrix form

p_values_bh_matrix <- matrix(p_values_bh, ncol = ncol(p_values_df), nrow = nrow(p_values_df))

# Convert corrected p-values matrices to data frames for easier viewing
p_values_bh_df <- as.data.frame(p_values_bh_matrix)

print(p_values_bh_df)

flextable::save_as_docx(flextable::flextable(p_values_bh_df), path = "p.values.Cor.Table.docx")



# Explorative -------------------------------------------------------------

# Regression analysis

model = lm(BDI_T0_Gesamt ~ Gender + Age + Age_onset + CDRISC_25_T0_sum, data = swdata)
summary(model)

# Regressions CDRISC and BDI ~ TAQ by age group


# Create a named list of formulas
predictors <- c("TAQ_Pos_a", "TAQ_Pos_b", "TAQ_Pos_c", 
                "TAQ_Neg_a", "TAQ_Neg_b", "TAQ_Neg_c")

# Function to run model and tidy result
get_regression_info <- function(dep_var, pred_var) {
  model <- lm(as.formula(paste(dep_var, "~", pred_var)), data = swdata)
  result <- tidy(model, conf.int = TRUE) %>% 
    filter(term == pred_var) %>%
    mutate(outcome = dep_var,
           predictor = pred_var) %>%
    select(outcome, predictor, estimate, std.error, statistic, p.value, conf.low, conf.high)
  return(result)
}

# Run regressions for CD-RISC and BDI
reg_results <- bind_rows(
  lapply(predictors, function(p) get_regression_info("CDRISC_25_T0_sum", p)),
  lapply(predictors, function(p) get_regression_info("BDI_T0_Gesamt", p))
)

# Rename for table display
reg_results <- reg_results %>%
  rename(`Point estimate` = estimate,
         `Standard error` = std.error,
         `t-value` = statistic,
         `p-value` = p.value,
         `CI lower` = conf.low,
         `CI upper` = conf.high) %>%
  mutate(`Regression` = paste(outcome, "←", predictor)) %>%
  select(`Regression`, `Point estimate`, `Standard error`, `t-value`, `p-value`, `CI lower`, `CI upper`)

# Round results
reg_results <- reg_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# View in console or export to Word/CSV
print(reg_results)

reg_results$`Standard error`




# Mediation models for subscales ------------------------------------------


# Model_a -----------------------------------------------------------------


mediation.modela <- "
CDRISC_25_T0_sum ~ a*TAQ_Pos_a + b*TAQ_Neg_a
BDI_T0_Gesamt ~ c*CDRISC_25_T0_sum
BDI_T0_Gesamt ~ d*TAQ_Pos_a + e*TAQ_Neg_a
ac := a*c
bc := b*c
total.pos := d + (a*c)
total.neg := e + (b*c)
"

fitb <- sem(mediation.modela, 
            data=swdata, 
            fixed.x = F, 
            std.ov = T, 
            missing = "fiml",
            se = "bootstrap", 
            bootstrap = 10000)

summary(fitb)

d_child <- roundallnumerics(
  as.data.frame(parameterEstimates(fitb)),
  3
)

term <- paste(d_child$lhs, d_child$op, d_child$rhs)

term <- gsub("BDI_T0_Gesamt", "BDI", term)
term <- gsub("CDRISC_25_T0_sum", "CDRISC", term)
term <- gsub("TAQ_Pos_a", "TAQ+ (child)", term)
term <- gsub("TAQ_Neg_a", "TAQ- (child)", term)
term <- gsub("~~", "↔", term)
term <- gsub("~", "←", term)
term <- gsub("ac := a*c", "TAQ+ indirect", term, fixed = TRUE)
term <- gsub("bc := b*c", "TAQ- indirect", term, fixed = TRUE)
term <- gsub("total.pos := d+(a*c)", "TAQ+ total", term, fixed = TRUE)
term <- gsub("total.neg := e+(b*c)", "TAQ- total", term, fixed = TRUE)

d_child <- BBmisc::dropNamed(d_child, drop = c("lhs", "op", "rhs", "label"))
d_child <- tibble::add_column(d_child, Term = term, .before = 1)
d_child <- d_child[!grepl("←1", d_child$Term), ]  # keep covariance lines!

names(d_child) <- c("Term", "Point estimate", "Standard error", "z", "p-value", "CI lower", "CI upper")
d_child$`p-value` <- p_labeller(d_child$`p-value`)

flextable::save_as_docx(flextable::flextable(d_child), path = "Mediation_Childhood.docx")



# Model_b -----------------------------------------------------------------


mediation.modelb <- "
CDRISC_25_T0_sum ~ a*TAQ_Pos_b + b*TAQ_Neg_b
BDI_T0_Gesamt ~ c*CDRISC_25_T0_sum
BDI_T0_Gesamt ~ d*TAQ_Pos_b + e*TAQ_Neg_b
ac := a*c
bc := b*c
total.pos := d + (a*c)
total.neg := e + (b*c)
"

fitc <- sem(mediation.modelb, 
            data=swdata, 
            fixed.x = F, 
            std.ov = T, 
            missing = "fiml",
            se = "bootstrap", 
            bootstrap = 10000)

summary(fitc)

d_youth <- roundallnumerics(
  as.data.frame(parameterEstimates(fitc)),
  3
)

term <- paste(d_youth$lhs, d_youth$op, d_youth$rhs)

term <- gsub("BDI_T0_Gesamt", "BDI", term)
term <- gsub("CDRISC_25_T0_sum", "CDRISC", term)
term <- gsub("TAQ_Pos_b", "TAQ+ (youth)", term)
term <- gsub("TAQ_Neg_b", "TAQ- (youth)", term)
term <- gsub("~~", "↔", term)
term <- gsub("~", "←", term)
term <- gsub("ac := a*c", "TAQ+ indirect", term, fixed = TRUE)
term <- gsub("bc := b*c", "TAQ- indirect", term, fixed = TRUE)
term <- gsub("total.pos := d+(a*c)", "TAQ+ total", term, fixed = TRUE)
term <- gsub("total.neg := e+(b*c)", "TAQ- total", term, fixed = TRUE)

d_youth <- BBmisc::dropNamed(d_youth, drop = c("lhs", "op", "rhs", "label"))
d_youth <- tibble::add_column(d_youth, Term = term, .before = 1)
d_youth <- d_youth[!grepl("←1", d_youth$Term), ]

names(d_youth) <- c("Term", "Point estimate", "Standard error", "z", "p-value", "CI lower", "CI upper")
d_youth$`p-value` <- p_labeller(d_youth$`p-value`)

flextable::save_as_docx(flextable::flextable(d_youth), path = "Mediation_Youth.docx")



# Model_c -----------------------------------------------------------------


mediation.modelc <- "
CDRISC_25_T0_sum ~ a*TAQ_Pos_c + b*TAQ_Neg_c
BDI_T0_Gesamt ~ c*CDRISC_25_T0_sum
BDI_T0_Gesamt ~ d*TAQ_Pos_c + e*TAQ_Neg_c
ac := a*c
bc := b*c
total.pos := d + (a*c)
total.neg := e + (b*c)
"

fitd <- sem(mediation.modelc, 
            data=swdata, 
            fixed.x = F, 
            std.ov = T, 
            missing = "fiml",
            se = "bootstrap", 
            bootstrap = 10000)

summary(fitd)

d_adult <- roundallnumerics(
  as.data.frame(parameterEstimates(fitd)),
  3
)

term <- paste(d_adult$lhs, d_adult$op, d_adult$rhs)

term <- gsub("BDI_T0_Gesamt", "BDI", term)
term <- gsub("CDRISC_25_T0_sum", "CDRISC", term)
term <- gsub("TAQ_Pos_c", "TAQ+ (adult)", term)
term <- gsub("TAQ_Neg_c", "TAQ- (adult)", term)
term <- gsub("~~", "↔", term)
term <- gsub("~", "←", term)
term <- gsub("ac := a*c", "TAQ+ indirect", term, fixed = TRUE)
term <- gsub("bc := b*c", "TAQ- indirect", term, fixed = TRUE)
term <- gsub("total.pos := d+(a*c)", "TAQ+ total", term, fixed = TRUE)
term <- gsub("total.neg := e+(b*c)", "TAQ- total", term, fixed = TRUE)

d_adult <- BBmisc::dropNamed(d_adult, drop = c("lhs", "op", "rhs", "label"))
d_adult <- tibble::add_column(d_adult, Term = term, .before = 1)
d_adult <- d_adult[!grepl("←1", d_adult$Term), ]

names(d_adult) <- c("Term", "Point estimate", "Standard error", "z", "p-value", "CI lower", "CI upper")
d_adult$`p-value` <- p_labeller(d_adult$`p-value`)

flextable::save_as_docx(flextable::flextable(d_adult), path = "Mediation_Adulthood.docx")



# Additional analyses -----------------------------------------------------

# TAQ subscales table  ----------------------------------------------------


# Create a summary table for TAQ subscales
taq_descriptives <- tibble(
  `TAQ developmental periods` = c(
    "TAQₚ", "TAQₚC", "TAQₚY", "TAQₚA",
    "TAQₙ", "TAQₙC", "TAQₙY", "TAQₙA"
  ),
  `Mean ± SD` = c(
    sprintf("%.2f ± %.2f", mean(swdata$TAQ_Pos, na.rm = TRUE), sd(swdata$TAQ_Pos, na.rm = TRUE)),
    sprintf("%.2f ± %.2f", mean(swdata$TAQ_Pos_a, na.rm = TRUE), sd(swdata$TAQ_Pos_a, na.rm = TRUE)),
    sprintf("%.2f ± %.2f", mean(swdata$TAQ_Pos_b, na.rm = TRUE), sd(swdata$TAQ_Pos_b, na.rm = TRUE)),
    sprintf("%.2f ± %.2f", mean(swdata$TAQ_Pos_c, na.rm = TRUE), sd(swdata$TAQ_Pos_c, na.rm = TRUE)),
    sprintf("%.2f ± %.2f", mean(swdata$TAQ_Neg, na.rm = TRUE), sd(swdata$TAQ_Neg, na.rm = TRUE)),
    sprintf("%.2f ± %.2f", mean(swdata$TAQ_Neg_a, na.rm = TRUE), sd(swdata$TAQ_Neg_a, na.rm = TRUE)),
    sprintf("%.2f ± %.2f", mean(swdata$TAQ_Neg_b, na.rm = TRUE), sd(swdata$TAQ_Neg_b, na.rm = TRUE)),
    sprintf("%.2f ± %.2f", mean(swdata$TAQ_Neg_c, na.rm = TRUE), sd(swdata$TAQ_Neg_c, na.rm = TRUE))
  )
)

# Print or export the table
print(taq_descriptives)

# Optional: Save as Word table
flextable::save_as_docx(flextable::flextable(taq_descriptives), path = "SI_TAQ_Descriptives.docx")


