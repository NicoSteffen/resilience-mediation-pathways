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



rep = foreign::read.spss("", to.data.frame = TRUE) 


#remove unnecessary Variables

prefixes <- c("KOMME", "NEO_FFI", "IPDE_", "MADRS_", "BSL_", 
              "IIPC_", "RSQ", "BIS_", "UCLA_", "SNI_", "BUL_", "STAI_", "VAR", "X.", "ES",
              "FNS", "Behav", "NTS", "pre_", "post_", "p2", "soe", "d_", "SICD")

cdrisk_vars <- paste0("CDRISK_", 1:25, "_alterFB")


remove_vars <- c(
  cdrisk_vars,
  unlist(lapply(prefixes, function(pref) grep(paste0("^", pref), names(rep), value = TRUE)))
)

rep <- rep[ , !(names(rep) %in% remove_vars)]

#remove TAQ Item 43
rep <- rep[, names(rep) != "TAQ_43"]

# remove missings CDRISC

cdrisc_cols <- grep("CDRISK", names(rep), value = TRUE)
rep$missing_CDRISC <- rowSums(is.na(rep[, cdrisc_cols])) # remove 7

# remove IDs with more than 3 missing
rep <- rep[rep$missing_CDRISC <= 3, ]

table(rep$group)


# Subset of BPD and CD
rep_subset <- rep[rep$group %in% c("BPD", "CD"), ]
rep_subset$group <- droplevels(rep_subset$group)

# meds
antidepressants <- c(
  "Agomelantin",         # Agomelatine
  "Amitriptylin",        # Tricyclic
  "Amioxid",             # Likely Amoxapin (tricyclic, antidepressant & antipsychotic profile)
  "Bupropion_Dosis",
  "Citalopram_Dosis",
  "Doxepin_Dosis",
  "Duloxetin",
  "Escitalopram_Dosis",
  "Fluoextin_Dosis",     # likely meant Fluoxetin
  "Mitrazapin_Dosis",    # likely Mirtazapin
  "Milnacipran",
  "Paroxetin_Dosis",
  "Sertralin",
  "Tianeurax",           # Tianeptin
  "Trazodon",
  "Trimipramin",
  "Tranylcypromin"
)


antipsychotics <- c(
  "Aripiprazol",
  "Melperon_Dosis",
  "Olanzapin",
  "Quetiapin_Dosis",
  "Risperidon",
  "Promethazin",
  "Prothipendyl"
)

mood_stabilizers_anticonvulsants <- c(
  "Ergenyl",        # Valproinsäure
  "Lamictal",       # Lamotrigin
  "Lithium",        # Klassiker unter den Stimmungsstabilisierern
  "Levetiracetam",  # Antikonvulsivum, off-label gelegentlich bei Aggression/Impulsivität
  "Oxcarbazepin",   # wie Carbamazepin – psychiatrisch relevant
  "Gabapentin",
  "Pregabalin"
)

other_psychopharm <- c(
  "adderall",
  "Methylphenidat",
  "Modafinil",
  "Clonidin",
  "Circadin",          # Melatonin – z.T. bei Schlafstörungen
  "Naltrexon"          # bei Sucht
)

psychopharmaka <- c(
  antidepressants,
  antipsychotics,
  mood_stabilizers_anticonvulsants,
  other_psychopharm
)




# new vars
rep$medication <- ifelse(
  rowSums(!is.na(rep[, psychopharmaka])) > 0, "yes", "no"
)

rep$antidepressants <- ifelse(
  rowSums(!is.na(rep[, antidepressants])) > 0, "yes", "no"
)

rep$antipsychotics <- ifelse(
  rowSums(!is.na(rep[, antipsychotics])) > 0, "yes", "no"
)

# in %
format_freq <- function(x) {
  tab <- table(x)
  total <- sum(tab)
  formatted <- paste0(tab, " (", round(100 * tab / total, 1), "%)")
  names(formatted) <- names(tab)
  return(formatted)
}

# Mean ± SD
format_mean_sd <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  sprintf("%.2f ± %.2f", m, s)
}

# Sex
geschlecht_tab <- format_freq(rep$Geschlecht)

# Group
group_tab <- format_freq(rep$group)

# Age (M ± SD)
alter_m_sd <- format_mean_sd(rep$Alter_echt)

# Age of onset (M ± SD)
alter_kb_m_sd <- format_mean_sd(rep$Alter_Krankheitsbeginn)

# Medication (yes/no)
medikation_tab <- format_freq(rep$medication)

# Antidepressiva (yes/no)
ad_tab <- format_freq(rep$antidepressants)

# Antipsychotika (yes/no)
ap_tab <- format_freq(rep$antipsychotics)


summary_table <- data.frame(
  Variable = c(
    paste("Geschlecht:", names(geschlecht_tab)),
    paste("Gruppe:", names(group_tab)),
    "Alter (echt)",
    "Alter bei Krankheitsbeginn",
    paste("Medikation:", names(medikation_tab)),
    paste("Antidepressiva:", names(ad_tab)),
    paste("Antipsychotika:", names(ap_tab))
  ),
  Wert = c(
    geschlecht_tab,
    group_tab,
    alter_m_sd,
    alter_kb_m_sd,
    medikation_tab,
    ad_tab,
    ap_tab
  ),
  row.names = NULL
)


print(summary_table, row.names = FALSE)


# Variablen
mean_sd <- function(x) {
  sprintf("%.2f ± %.2f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

# Patients (rep_bpd_cd)
taq_pos_pat <- mean_sd(rep_bpd_cd$TAQ_Pos)
taq_neg_pat <- mean_sd(rep_bpd_cd$TAQ_Neg)
bdi_pat      <- mean_sd(rep_bpd_cd$BDI_Score)
cdrisc_pat   <- mean_sd(rep_bpd_cd$CDRISC_10_sum)







# 2 TAQ Items "Family Secrets" wurden aus späteren Versionen ausgeschlossen
# https://complextrauma.org/wp-content/uploads/2019/03/TAQ-2019.pdf

# CAVE: In our (German) version there are only 3 age categories (7-12; 13-18; adult). The english version had 4 (0-6; 7-12; 13-18; adult)

rep = rep[, -c(which(grepl("TAQ_9", names(rep))))]
rep = rep[, -c(which(grepl("TAQ_24", names(rep))))]

# renumber TAQ items after 2-item exclusion

names(rep)[grep("TAQ", names(rep))] = paste0("TAQ_", rep(1:40, each = 3), rep(letters[1:3], 40))

# convert item to numeric values
# recode 2nd TAQ item

for(j in which(grepl("TAQ_2", names(rep)))){
  rep[,j] = as.character(rep[,j])
  rep[,j][rep[,j] == "oft/stark"] = "0"
  rep[,j][rep[,j] == "Gelegentlich/im Grunde schon"] = "0"
  rep[,j][rep[,j] == "selten/ein wenig"] = "2"
  rep[,j][rep[,j] == "nie/gar nicht"] = "3"
  rep[,j][rep[,j] == "weiß nicht"] = "0"
  rep[,j] = as.numeric(rep[,j])
}



for(j in which(grepl("TAQ", names(rep)))){
  rep[,j] = as.character(rep[,j])
  rep[,j][rep[,j] == "oft/stark"] = "3"
  rep[,j][rep[,j] == "Gelegentlich/im Grunde schon"] = "2"
  rep[,j][rep[,j] == "selten/ein wenig"] = "0"
  rep[,j][rep[,j] == "nie/gar nicht"] = "0"
  rep[,j][rep[,j] == "nein"] = "0"
  rep[,j][rep[,j] == "ja"] = "3"
  rep[,j][rep[,j] == "weiß nicht"] = "0"
  rep[,j] = as.numeric(rep[,j])
}

#View(rep[,grep("^TAQ_", names(rep))])

for(j in which(grepl("CDRISK", names(rep)))){
  rep[,j] = as.character(rep[,j])
  rep[,j][rep[,j] == "überhaupt nicht wahr"] = "0"
  rep[,j][rep[,j] == "selten wahr"] = "1"
  rep[,j][rep[,j] == "manchmal wahr"] = "2"
  rep[,j][rep[,j] == "oft wahr"] = "3"
  rep[,j][rep[,j] == "fast immer wahr"] = "4"
  rep[,j] = as.numeric(rep[,j])
}

#View(rep[,grep("^CDRISK_", names(rep))])

#Mean Imputation CDRISC
cdrisc_mat <- rep[, cdrisc_cols]
row_means <- rowMeans(cdrisc_mat, na.rm = TRUE)

# Schleife oder apply zum Ersetzen der Missings
for (i in 1:nrow(cdrisc_mat)) {
  na_idx <- is.na(cdrisc_mat[i, ])
  cdrisc_mat[i, na_idx] <- row_means[i]
}

# Ersetzte Matrix wieder in den Datensatz einfügen
rep[, cdrisc_cols] <- cdrisc_mat


rep$CDRISC_10_sum = rowSums(rep[, grep("CDRISK", names(rep))])


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
    rep[[varname]] <- calc_subscale(rep, subscales[[subscale]][[period]])
  }
}



rep$TAQ_Competence_total = rowSums(rep[, c("TAQ_Competence_a", "TAQ_Competence_b", "TAQ_Competence_c" )], na.rm = TRUE)
rep$TAQ_Safety_total = rowSums(rep[, c("TAQ_Safety_a", "TAQ_Safety_b", "TAQ_Safety_c")], na.rm = TRUE)
rep$TAQ_Neglect_total = rowSums(rep[, c("TAQ_Neglect_a", "TAQ_Neglect_b", "TAQ_Neglect_c")], na.rm = TRUE)
rep$TAQ_Separation_total = rowSums(rep[, c("TAQ_Separation_a", "TAQ_Separation_b", "TAQ_Separation_c")], na.rm = TRUE)
rep$TAQ_Emotional_Abuse_total = rowSums(rep[, c("TAQ_Emotional_Abuse_a", "TAQ_Emotional_Abuse_b", "TAQ_Emotional_Abuse_c")], na.rm = TRUE)
rep$TAQ_Physical_Abuse_total = rowSums(rep[, c("TAQ_Physical_Abuse_a", "TAQ_Physical_Abuse_b", "TAQ_Physical_Abuse_c")], na.rm = TRUE)
rep$TAQ_Sexual_Abuse_total = rowSums(rep[, c("TAQ_Sexual_Abuse_a", "TAQ_Sexual_Abuse_b", "TAQ_Sexual_Abuse_c")], na.rm = TRUE)
rep$TAQ_Witnessing_total = rowSums(rep[, c("TAQ_Witnessing_a", "TAQ_Witnessing_b", "TAQ_Witnessing_c")], na.rm = TRUE)
rep$TAQ_Impersonal_Trauma_total = rowSums(rep[, c("TAQ_Impersonal_Trauma_a", "TAQ_Impersonal_Trauma_b", "TAQ_Impersonal_Trauma_c")], na.rm = TRUE)
rep$TAQ_Alcohol_Drugs_total = rowSums(rep[, c("TAQ_Alcohol_Drugs_a", "TAQ_Alcohol_Drugs_b", "TAQ_Alcohol_Drugs_c")], na.rm = TRUE)


# Fetch Subscale names

TAQ_subscales = names(rep[, grep("TAQ_[a-zA-Z]", names(rep))])
TAQ_subscales = unique(gsub("TAQ_|_a|_b|_c", "", TAQ_subscales))




# Multiple Mediation with CDRisk ------------------------------------------

#create new TAQ_Pos + TAQ_Neg

rep$TAQ_Pos = rowMeans(rep[,c("TAQ_Safety_total","TAQ_Competence_total")])

#subscales pos
rep$TAQ_Pos_a = rowMeans(rep[,c("TAQ_Safety_a","TAQ_Competence_a")], na.rm = T)
rep$TAQ_Pos_b = rowMeans(rep[,c("TAQ_Safety_b","TAQ_Competence_b")], na.rm = T)
rep$TAQ_Pos_c = rowMeans(rep[,c("TAQ_Safety_c","TAQ_Competence_c")], na.rm = T)

rep$TAQ_Neg = rowMeans(rep[, c(
  "TAQ_Neglect_total",
  "TAQ_Separation_total",
  "TAQ_Emotional_Abuse_total",
  "TAQ_Physical_Abuse_total",
  "TAQ_Sexual_Abuse_total",
  "TAQ_Witnessing_total",
  "TAQ_Impersonal_Trauma_total",
  "TAQ_Alcohol_Drugs_total"
)])


#subscales pos

rep$TAQ_Neg_a = rowMeans(rep[, c(
  "TAQ_Neglect_a",
  "TAQ_Separation_a",
  "TAQ_Emotional_Abuse_a",
  "TAQ_Physical_Abuse_a",
  "TAQ_Sexual_Abuse_a",
  "TAQ_Witnessing_a",
  "TAQ_Impersonal_Trauma_a",
  "TAQ_Alcohol_Drugs_a"
)], na.rm = T)

rep$TAQ_Neg_b = rowMeans(rep[, c(
  "TAQ_Neglect_b",
  "TAQ_Separation_b",
  "TAQ_Emotional_Abuse_b",
  "TAQ_Physical_Abuse_b",
  "TAQ_Sexual_Abuse_b",
  "TAQ_Witnessing_b",
  "TAQ_Impersonal_Trauma_b",
  "TAQ_Alcohol_Drugs_b"
)], na.rm = T)

rep$TAQ_Neg_c = rowMeans(rep[, c(
  "TAQ_Neglect_c",
  "TAQ_Separation_c",
  "TAQ_Emotional_Abuse_c",
  "TAQ_Physical_Abuse_c",
  "TAQ_Sexual_Abuse_c",
  "TAQ_Witnessing_c",
  "TAQ_Impersonal_Trauma_c",
  "TAQ_Alcohol_Drugs_c"
)], na.rm = T)


# create subset for PDD and BPD

rep_bpd_cd <- subset(rep, group_1 %in% c("BPD", "CD"))



# Mediation total ---------------------------------------------------------

mediation.model <- "
CDRISC_10_sum ~ a*TAQ_Pos + b*TAQ_Neg
BDI_Score ~ c*CDRISC_10_sum
BDI_Score ~ d*TAQ_Pos + e*TAQ_Neg
ac := a*c
bc := b*c
total.pos := d + (a*c)
total.neg := e + (b*c)
"

set.seed(123)

fita <- sem(mediation.model, 
            data=rep_bpd_cd, 
            fixed.x = F, 
            std.ov = T, 
            missing = "fiml",
            se = "bootstrap", 
            bootstrap = 10000)

summary(fita)


semTable::sem_tables(fita)

d = roundallnumerics(
  as.data.frame(
    parameterEstimates(fita)
  ), 3
)

term = paste(d$lhs, d$op, d$rhs)

term = gsub("BDI_score", "BDI", term)
term = gsub("CDRISC_10_sum", "CDRISC", term)
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

flextable::save_as_docx(flextable::flextable(d), path = "Table_Replication_Clinical.docx")






