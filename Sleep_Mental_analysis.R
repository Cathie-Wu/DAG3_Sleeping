library(tidyverse)
library(vegan)
library(ape) 


# Load and clean ----------------------------------------------------------

load("metaphlan3_cleaned.RData")
metadata <- read.csv("DAG3_metadata_merged_ready_v27.csv")
sleep <- read_tsv("Sleep/LLs_sleep_data.tsv")

## Sleep phenotype ---------------
sleep_2 <- sleep %>% 
  select(DAG3_sampleID,psqi_globalscore_adu_c_1,psqi_qualityscore_adu_c_1) %>% 
  filter(DAG3_sampleID!="NA"
         & psqi_globalscore_adu_c_1!="NA"
         & psqi_globalscore_adu_c_1!= "$6") %>% 
  mutate(psqi = as.numeric(psqi_globalscore_adu_c_1), .keep= "unused") %>% 
  mutate(psqi_bi = as.numeric(psqi_qualityscore_adu_c_1), .keep= "unused") %>% 
  mutate(psqi_norm = qnorm((rank(psqi,na.last="keep")-0.5)/sum(!is.na(psqi)))) %>% 
  mutate(psqi_cat=cut(psqi, c(-1,5,10,20), right=T,labels=c("good","poor","very poor"))) %>% 
  mutate(psqi_cat_num=ifelse(psqi>5,ifelse(psqi>10,2,1),0))
         

## metadata -------------------
sleep_2 %>% inner_join(metadata,by="DAG3_sampleID")-> Sleep_meta

Sleep_meta$EXP.DIET.Probiotics <- factor(Sleep_meta$EXP.DIET.Probiotics,levels=c( "N", "Y.sometimes", "Y.often","Y.always"))
Sleep_meta$EXP.EARLYLIFE.LIVINGPLACE.child.1.4 <- factor(Sleep_meta$EXP.EARLYLIFE.LIVINGPLACE.child.1.4,levels=c("Other", "Village.rur", "Farmhouse", "Suburbs", "S.city.L.village", "City centre"))
Sleep_meta$MED.DISEASES.Neurological.Headaches.HowOften <- factor(Sleep_meta$MED.DISEASES.Neurological.Headaches.HowOften,levels=c("not at all", "several times a month",  "several times a week","every day"))
Sleep_meta$MED.HEALTH.RAND.Health.Change.1y <- factor(Sleep_meta$MED.HEALTH.RAND.Health.Change.1y,levels=c("much worse than a year ago","a little bit worse than a year ago", "about the same as a year ago", "a little bit better than a year ago", "much better than a year ago" ))
Sleep_meta$MED.HEALTH.RAND.Health.General <- factor(Sleep_meta$MED.HEALTH.RAND.Health.General,levels=c("poor", "mediocre", "good", "very good", "excellent"))
Sleep_meta$MED.INDICES.Fibrosis.Score.T1.lvl <- factor(Sleep_meta$MED.INDICES.Fibrosis.Score.T1.lvl,levels=c("Minimal fibrosis", "Intermediate fibrosis", "Severe fibrosis"))
Sleep_meta$MED.INDICES.FattyLiverIndex.T1.Class <- factor(Sleep_meta$MED.INDICES.FattyLiverIndex.T1.Class,levels=c("Low", "Intermediate", "High"))
Sleep_meta$MED.INDICES.NAFLD.T1.CLASS <- factor(Sleep_meta$MED.INDICES.NAFLD.T1.CLASS,levels=c("F0-F2", "Intermediate score", "F3-F4"))
Sleep_meta$SOCIOEC.INCOME.Income.Month <- factor(Sleep_meta$SOCIOEC.INCOME.Income.Month,levels=c("E.750.less", "E.750.to.1000", "E.1000.to.1500", "E.1500.to.2000", "E.2000.to.2500","E.2500.to.3000", "E.3000.to.3500", "E.3500.more" ))
Sleep_meta$META.POOP.COLLECTION_SEASON <- factor(Sleep_meta$META.POOP.COLLECTION_SEASON,levels=c("Spring", "Summer", "Fall", "Winter"))
Sleep_meta$MED.DISEASES.Rome3_IBS.factor <- factor(Sleep_meta$MED.DISEASES.Rome3_IBS.factor,levels=c("N", "IBS", "IBS.C", "IBS.D", "IBS.M"))

## mental -------------------
phenotypes <- read.csv("mental_health_sleep/01.11.2021_mental_health_Sio_subset.csv")
phenotypes$MDD_control <- factor(phenotypes$MDD_control, levels=c("control", "current_MDD"))
phenotypes$AnyDep_control <- factor(phenotypes$AnyDep_control, levels=c("control", "any_current_dep"))
phenotypes$LTMDDonly_control <- factor(phenotypes$LTMDDonly_control, levels=c("control", "lifetime_MDD_only"))
phenotypes$GAD_control <- factor(phenotypes$GAD_control, levels=c("control", "current_GAD"))
phenotypes$AnyAnx_control <- factor(phenotypes$AnyAnx_control, levels=c("control", "any_current_anx"))
phenotypes <- phenotypes %>% rename(LTGADonly_control=LTGADDonly_control)
phenotypes$LTGADonly_control <- factor(phenotypes$LTGADonly_control, levels=c("control", "lifetime_GAD_only"))
phenotypes$MED.MEDS.Antipsychotics_atypical_ATC_N05AH <- as.factor(phenotypes$MED.MEDS.Antipsychotics_atypical_ATC_N05AH)
phenotypes$MED.MEDS.Benzodiazepine_agonists_ATC_N05 <- as.factor(phenotypes$MED.MEDS.Benzodiazepine_agonists_ATC_N05)
phenotypes$MED.MEDS.Psychoanaleptics_ATC_N06 <- as.factor(phenotypes$MED.MEDS.Psychoanaleptics_ATC_N06)
phenotypes$MED.MEDS.Psycholeptics_ATC_N05 <- as.factor(phenotypes$MED.MEDS.Psycholeptics_ATC_N05)
phenotypes$MED.MEDS.Seratonine_uptake_inhibitors_selective_ATC_N06AB <- as.factor(phenotypes$MED.MEDS.Seratonine_uptake_inhibitors_selective_ATC_N06AB)
phenotypes$MED.MEDS.Tricyclic_antidepressants_ATC_N06AA <- as.factor(phenotypes$MED.MEDS.Tricyclic_antidepressants_ATC_N06AA)
phenotypes$CURRENT_TREATMENT.Mental_Health <- as.factor(phenotypes$CURRENT_TREATMENT.Mental_Health)
mental <- phenotypes %>% select(-X,-PSEUDOIDEXT,-PROJECT_PSEUDO_ID) %>% select(DAG3_sampleID:CURRENT_TREATMENT.Mental_Health)

Sleep_meta %>% select(-psqi,-psqi_bi) %>% inner_join(mental,by="DAG3_sampleID")-> Sleep_meta

## microbiome ------------------
colnames(taxa_matrix)[grepl("\\[", colnames(taxa_matrix))] <- "s__Collinsella_massiliensis"
colnames(df_species)[grepl("\\[", colnames(df_species))] <- "s__Collinsella_massiliensis"

df_species <-df_species %>% mutate(DAG3_sampleID=gsub("_metaphlan", "", df_species$ID)) %>% select(-1) %>% relocate(DAG3_sampleID)
  

taxa_matrix <-taxa_matrix %>% mutate(DAG3_sampleID=gsub("_metaphlan", "", taxa_matrix$ID)) %>% select(-1) %>% relocate(DAG3_sampleID)

save(df_species,taxa_matrix,Sleep_meta,file="Sleep_mental_clean.RData")


# alpha ----------------------------------------------------------

Sleep_meta %>% select(DAG3_sampleID,psqi_cat,psqi_norm,psqi_cat_num,ANTHRO.AGE,ANTHRO.BMI,ANTHRO.Sex,META.BATCH,META.DNA.conc.ng.ul,META.DNA.postclean.reads,META.POOP.BristolMean,META.POOP.Freq,AnyDep_control,AnyAnx_control,CURRENT_TREATMENT.Mental_Health)-> Sleep_meta_2

Sleep_meta_2[complete.cases(Sleep_meta_2),]-> Sleep_meta_2
print(dim(Sleep_meta_2))

Sleep_meta %>% count(AnyDep_control,AnyAnx_control)

alpha <- df_species %>% semi_join(Sleep_meta_2,by="DAG3_sampleID") 
alpha$alpha <- diversity(alpha[,2:dim(alpha)[2]], index = "shannon")

alpha %>% select("DAG3_sampleID","alpha") %>%  inner_join(Sleep_meta_2,by="DAG3_sampleID")-> alpha

model <- lm(alpha ~ psqi_cat_num + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq + AnyDep_control + CURRENT_TREATMENT.Mental_Health,data=alpha)

summary(model)
as.data.frame(summary(model)$coefficients) -> model_results

model_results <- model_results %>% rownames_to_column() %>% as_tibble()
write.csv(model_results,file = "output_mental/sleep_mental_alpha_model2.csv")

# Taxa --------------------------------------------------------------------

## CLR transformation ------------------------------------------------------
# df_species <-df_species %>% semi_join(Sleep_meta_2,by="DAG3_sampleID") 
# taxa_matrix <-taxa_matrix %>% semi_join(Sleep_meta_2,by="DAG3_sampleID")


print("filter out taxa present in more than 100 samples")

taxa_name <- c()
prevalence  <- c()
for (i in 1:dim(taxa_matrix)[2]){
  taxa_name <- c(taxa_name, colnames(taxa_matrix[i]))
  taxa <- taxa_matrix[[i]]
  prevalence  <- c(prevalence,length(taxa[!taxa == 0]))
}

taxa_prevalence  <- tibble(taxa=taxa_name,prevalence =prevalence)

filter <- taxa_prevalence %>% filter(prevalence>=100)%>% filter(taxa!="DAG3_sampleID") %>% pull(1)

print(paste(length(filter),"taxa selected"))

core_species <- df_species %>% select(-"DAG3_sampleID")
taxa_matrix2 <- taxa_matrix %>% select(-"DAG3_sampleID")

print("sudo count")
if(any(core_species==0)) core_species = core_species + min(core_species[core_species>0])/2
if(any(taxa_matrix2==0)) taxa_matrix2 = taxa_matrix2 + min(taxa_matrix2[taxa_matrix2>0])/2

print("calculate geometric mean for each sample")
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}
Gmean_core = apply(core_species, 1, gm_mean)
data_prepared = cbind(Gmean_core,taxa_matrix2)

print("CLR transformation")
data_transformed = t(apply(data_prepared,1,function(x){
  log(x / x[1])[-1]
}))

data_transformed <- as_tibble(data_transformed)
data_transformed <- data_transformed %>% select(filter)
data_transformed %>% mutate(DAG3_sampleID=taxa_matrix$DAG3_sampleID) -> data_transformed

final_data <- left_join(Sleep_meta_2,data_transformed)

## taxa ~ psqi_cat -------------------------------------------------------------
n <- c()
e <- c()
s <- c()
prev <- c()
r <- c()
p <- c()

for (col in filter){
  formula <- as.formula(paste(col,"~psqi_cat_num + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq + AnyDep_control + CURRENT_TREATMENT.Mental_Health", collapse= "" ))
  model <- lm(formula,data=final_data)
  result <- summary(model)
  n <- c(n,col)
  r <- c(r,result$adj.r.squared)
  p <- c(p,result$coefficients[2,4])
  e <- c(e,result$coefficients[2,1])
  s <- c(s,result$coefficients[2,2])
  prev_0 <- taxa_prevalence %>% filter(taxa==col) %>% pull(prevalence)
  prev <- c(prev,prev_0)
}
df_result <- tibble(Taxa=n,Prevalence=prev,Estimate=e,Std.Error=s,r=r,p_value=p,p_adjust = p.adjust(p, "fdr"))
df_result <- df_result %>% mutate(sig = ifelse(p_adjust< 0.05,T,F))

write.csv(df_result,file = "output_mental/sleep_mental_taxa_model_cat_num_adjusted.csv")

## taxa ~ psqi_norm -------------------------------------------------------------
n <- c()
e <- c()
s <- c()
prev <- c()
r <- c()
p <- c()

for (col in filter){
  formula <- as.formula(paste(col,"~psqi_norm + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq + AnyDep_control + CURRENT_TREATMENT.Mental_Health", collapse= "" ))
  model <- lm(formula,data=final_data)
  result <- summary(model)
  n <- c(n,col)
  r <- c(r,result$adj.r.squared)
  p <- c(p,result$coefficients[2,4])
  e <- c(e,result$coefficients[2,1])
  s <- c(s,result$coefficients[2,2])
  prev_0 <- taxa_prevalence %>% filter(taxa==col) %>% pull(prevalence)
  prev <- c(prev,prev_0)
}
df_result <- tibble(Taxa=n,Prevalence=prev,Estimate=e,Std.Error=s,r=r,p_value=p,p_adjust = p.adjust(p, "fdr"))
df_result <- df_result %>% mutate(sig = ifelse(p_adjust< 0.05,T,F))

write.csv(df_result,file = "output_mental/sleep_mental_taxa_model_norm_adjusted.csv")

## taxa ~ psqi_cat unadjusted -------------------------------------------------------------
n <- c()
e <- c()
s <- c()
prev <- c()
r <- c()
p <- c()

for (col in filter){
  formula <- as.formula(paste(col,"~psqi_cat_num + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq", collapse= "" ))
  model <- lm(formula,data=final_data)
  result <- summary(model)
  n <- c(n,col)
  r <- c(r,result$adj.r.squared)
  p <- c(p,result$coefficients[2,4])
  e <- c(e,result$coefficients[2,1])
  s <- c(s,result$coefficients[2,2])
  prev_0 <- taxa_prevalence %>% filter(taxa==col) %>% pull(prevalence)
  prev <- c(prev,prev_0)
}
df_result <- tibble(Taxa=n,Prevalence=prev,Estimate=e,Std.Error=s,r=r,p_value=p,p_adjust = p.adjust(p, "fdr"))
df_result <- df_result %>% mutate(sig = ifelse(p_adjust< 0.05,T,F))

write.csv(df_result,file = "output_mental/sleep_mental_taxa_model_cat_num_unadjusted.csv")

## taxa ~ psqi_norm unadjusted -------------------------------------------------------------
n <- c()
e <- c()
s <- c()
prev <- c()
r <- c()
p <- c()

for (col in filter){
  formula <- as.formula(paste(col,"~psqi_norm + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq", collapse= "" ))
  model <- lm(formula,data=final_data)
  result <- summary(model)
  n <- c(n,col)
  r <- c(r,result$adj.r.squared)
  p <- c(p,result$coefficients[2,4])
  e <- c(e,result$coefficients[2,1])
  s <- c(s,result$coefficients[2,2])
  prev_0 <- taxa_prevalence %>% filter(taxa==col) %>% pull(prevalence)
  prev <- c(prev,prev_0)
}
df_result <- tibble(Taxa=n,Prevalence=prev,Estimate=e,Std.Error=s,r=r,p_value=p,p_adjust = p.adjust(p, "fdr"))
df_result <- df_result %>% mutate(sig = ifelse(p_adjust< 0.05,T,F))

write.csv(df_result,file = "output_mental/sleep_mental_taxa_model_norm_unadjusted.csv")

# test ----------------------------------------------------------

Sleep_meta %>% select(DAG3_sampleID,psqi_cat,psqi_norm,psqi_cat_num,ANTHRO.AGE,ANTHRO.BMI,ANTHRO.Sex,META.BATCH,META.DNA.conc.ng.ul,META.DNA.postclean.reads,META.POOP.BristolMean,META.POOP.Freq)-> Sleep_meta_2

Sleep_meta_2[complete.cases(Sleep_meta_2),]-> Sleep_meta_2
print(dim(Sleep_meta_2))

final_data <- left_join(Sleep_meta_2,data_transformed)

## taxa ~ psqi_norm unadjusted -------------------------------------------------------------
n <- c()
e <- c()
s <- c()
prev <- c()
r <- c()
p <- c()

for (col in filter){
  formula <- as.formula(paste(col,"~psqi_norm + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq", collapse= "" ))
  model <- lm(formula,data=final_data)
  result <- summary(model)
  n <- c(n,col)
  r <- c(r,result$adj.r.squared)
  p <- c(p,result$coefficients[2,4])
  e <- c(e,result$coefficients[2,1])
  s <- c(s,result$coefficients[2,2])
  prev_0 <- taxa_prevalence %>% filter(taxa==col) %>% pull(prevalence)
  prev <- c(prev,prev_0)
}
df_result <- tibble(Taxa=n,Prevalence=prev,Estimate=e,Std.Error=s,r=r,p_value=p,p_adjust = p.adjust(p, "fdr"))
df_result <- df_result %>% mutate(sig = ifelse(p_adjust< 0.05,T,F))

write.csv(df_result,file = "output_mental/sleep_mental_taxa_model_norm_unadjusted_6563.csv")

## taxa ~ psqi_cat unadjusted -------------------------------------------------------------
n <- c()
e <- c()
s <- c()
prev <- c()
r <- c()
p <- c()

for (col in filter){
  formula <- as.formula(paste(col,"~psqi_cat_num + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq", collapse= "" ))
  model <- lm(formula,data=final_data)
  result <- summary(model)
  n <- c(n,col)
  r <- c(r,result$adj.r.squared)
  p <- c(p,result$coefficients[2,4])
  e <- c(e,result$coefficients[2,1])
  s <- c(s,result$coefficients[2,2])
  prev_0 <- taxa_prevalence %>% filter(taxa==col) %>% pull(prevalence)
  prev <- c(prev,prev_0)
}
df_result <- tibble(Taxa=n,Prevalence=prev,Estimate=e,Std.Error=s,r=r,p_value=p,p_adjust = p.adjust(p, "fdr"))
df_result <- df_result %>% mutate(sig = ifelse(p_adjust< 0.05,T,F))

write.csv(df_result,file = "output_mental/sleep_mental_taxa_model_cat_num_unadjusted_6563.csv")
