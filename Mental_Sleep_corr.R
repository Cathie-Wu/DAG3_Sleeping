library(tidyverse)

# Mental ------------------------------------------------------------------

phenotypes <- read.csv("mental_health_sleep/01.11.2021_mental_health_Sio_subset.csv")

#1. factorising variables

#a. case-controls
#NOTE: controls n = 5522

#MDD control
phenotypes$MDD_control <- factor(phenotypes$MDD_control, levels=c("control", "current_MDD"))
table(phenotypes$MDD_control, useNA = "ifany")
#MDD = 156

#AnyDep control
phenotypes$AnyDep_control <- factor(phenotypes$AnyDep_control, levels=c("control", "any_current_dep"))
table(phenotypes$AnyDep_control, useNA = "ifany")
#AnyDep = 226

#LTMDDonly_control
phenotypes$LTMDDonly_control <- factor(phenotypes$LTMDDonly_control, levels=c("control", "lifetime_MDD_only"))
table(phenotypes$LTMDDonly_control, useNA = "ifany")
#LTMDDonly = 749

#GAD control
phenotypes$GAD_control <- factor(phenotypes$GAD_control, levels=c("control", "current_GAD"))
table(phenotypes$GAD_control, useNA = "ifany")
#GAD = 339

#AnyAnx control
phenotypes$AnyAnx_control <- factor(phenotypes$AnyAnx_control, levels=c("control", "any_current_anx"))
table(phenotypes$AnyAnx_control, useNA = "ifany")
#AnyAnx = 385

#LTGADonly_control
# phenotypes <- phenotypes %>% mutate(LTGADonly_control=LTGADDonly_control,.keep="unused") %>% relocate(LTGADonly_control,.after = AnyAnx_control)
phenotypes <- phenotypes %>% rename(LTGADonly_control=LTGADDonly_control)
phenotypes$LTGADonly_control <- factor(phenotypes$LTGADonly_control, levels=c("control", "lifetime_GAD_only"))
table(phenotypes$LTGADonly_control, useNA = "ifany")
#LTGADonly = 243

#b. treatments
phenotypes$MED.MEDS.Antipsychotics_atypical_ATC_N05AH <- as.factor(phenotypes$MED.MEDS.Antipsychotics_atypical_ATC_N05AH)
table(phenotypes$MED.MEDS.Antipsychotics_atypical_ATC_N05AH, useNA = "ifany")
#N    Y
#7633   23

phenotypes$MED.MEDS.Benzodiazepine_agonists_ATC_N05 <- as.factor(phenotypes$MED.MEDS.Benzodiazepine_agonists_ATC_N05)
table(phenotypes$MED.MEDS.Benzodiazepine_agonists_ATC_N05, useNA = "ifany")
#N    Y
#7617   39

phenotypes$MED.MEDS.Psychoanaleptics_ATC_N06 <- as.factor(phenotypes$MED.MEDS.Psychoanaleptics_ATC_N06)
table(phenotypes$MED.MEDS.Psychoanaleptics_ATC_N06, useNA = "ifany")
#N    Y
#7623   33

phenotypes$MED.MEDS.Psycholeptics_ATC_N05 <- as.factor(phenotypes$MED.MEDS.Psycholeptics_ATC_N05)
table(phenotypes$MED.MEDS.Psycholeptics_ATC_N05, useNA = "ifany")
#N    Y
#7637   19

phenotypes$MED.MEDS.Seratonine_uptake_inhibitors_selective_ATC_N06AB <- as.factor(phenotypes$MED.MEDS.Seratonine_uptake_inhibitors_selective_ATC_N06AB)
table(phenotypes$MED.MEDS.Seratonine_uptake_inhibitors_selective_ATC_N06AB, useNA = "ifany")
#N    Y
#7422  234

phenotypes$MED.MEDS.Tricyclic_antidepressants_ATC_N06AA <- as.factor(phenotypes$MED.MEDS.Tricyclic_antidepressants_ATC_N06AA)
table(phenotypes$MED.MEDS.Tricyclic_antidepressants_ATC_N06AA, useNA = "ifany")
#N    Y
#7585   71

phenotypes$CURRENT_TREATMENT.Mental_Health <- as.factor(phenotypes$CURRENT_TREATMENT.Mental_Health)
table(phenotypes$CURRENT_TREATMENT.Mental_Health, useNA = "ifany")
#   N    Y
#7283  373

#c. covariates
phenotypes$ANTHRO.Sex <- factor(phenotypes$ANTHRO.Sex, levels = c("M", "F"))
str(phenotypes$ANTHRO.Sex)

phenotypes$META.BATCH <- factor(phenotypes$META.BATCH, levels = c ("dag3_batch1", "dag3_batch2", "dag3_batch3", "dag3_batch4", "dag3_batch5", "dag3_batch6", "dag3_batch7"))
str(phenotypes$META.BATCH)

mental <- phenotypes %>% select(-X,-PSEUDOIDEXT,-PROJECT_PSEUDO_ID) %>% select(DAG3_sampleID:CURRENT_TREATMENT.Mental_Health)


# count mental ------------------------------------------------------------
mental_count <- tibble(a=c(1:3))

for (i in 2:dim(mental)[2]){
  count <- mental[,i] %>% table(useNA = "ifany") %>% as.vector()
  if (length(count)==2){count <- c(count,NA)}
  mental_count <- mental_count %>% mutate(x=count)
  colnames(mental_count)[i] <- colnames(mental)[i]
}


mental_count1 <- mental_count %>% select(-a) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% as.tibble()
write.csv(mental_count1, "output_mental/mental_count1.csv")

# Sleep -------------------------------------------------------------------
sleep <- read_tsv("Sleep/LLs_sleep_data.tsv")

# Sleep phenotype
sleep_2 <- sleep %>% 
  select(DAG3_sampleID,psqi_globalscore_adu_c_1,psqi_qualityscore_adu_c_1) %>% 
  filter(DAG3_sampleID!="NA"
         & psqi_globalscore_adu_c_1!="NA"
         & psqi_globalscore_adu_c_1!= "$6") %>% 
  mutate(psqi = as.numeric(psqi_globalscore_adu_c_1), .keep= "unused") %>% 
  mutate(psqi_bi = as.numeric(psqi_qualityscore_adu_c_1), .keep= "unused") %>% 
  mutate(psqi_norm = qnorm((rank(psqi,na.last="keep")-0.5)/sum(!is.na(psqi)))) %>% 
  mutate(psqi_cat=cut(psqi, c(-1,5,10,20), right=T,labels=c("good","poor","very poor")))


sleep_2 %>% select(-psqi,-psqi_bi) %>% inner_join(mental,by="DAG3_sampleID")-> Sleep_mental

# count mental+sleep ------------------------
mental_count2 <- tibble(a=c(1:3))

for (i in 3:dim(Sleep_mental)[2]){
  count <- Sleep_mental[,i] %>% table(useNA = "ifany") %>% as.vector()
  if (length(count)==2){count <- c(count,NA)}
  mental_count2 <- mental_count2 %>% mutate(x=count)
  colnames(mental_count2)[i-1] <- colnames(Sleep_mental)[i]
}

mental_count2 <- mental_count2 %>% select(-a) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% as.tibble()
write.csv(mental_count2, "output_mental/mental_count2.csv")


# correlation cat------------------------
n <- c()
r <- c()
p <- c()

print("Start Loop")
for (column in colnames(Sleep_mental)){
  if (column %in% c("DAG3_sampleID","psqi_norm","psqi_cat")){ next }
  Sleep_mental %>% select(c(column, "psqi_cat")) %>% drop_na() -> Test_corr
  Test_corr <- lapply(Test_corr, function(x) as.numeric(x))
  cor.test(as_vector(Test_corr[1]),as_vector(Test_corr[2]), method = "spearman") -> test
  n <- c(n,column)
  r <- c(r,test$estimate)
  p <- c(p,test$p.value)
}

print("Combine Tibble")
sleep_mental_correlation <- tibble(mental=n,r=r,p_value=p)
sleep_mental_correlation <- sleep_mental_correlation %>% mutate(p_adjust = p.adjust(p_value, "fdr"), sig = ifelse(r>0.1 & p_adjust< 0.05,T,F))

write.csv(sleep_mental_correlation, "output_mental/sleep_mental_corr_cat.csv")

# correlation norm------------------------
n <- c()
r <- c()
p <- c()

print("Start Loop")
for (column in colnames(Sleep_mental)){
  if (column %in% c("DAG3_sampleID","psqi_norm","psqi_cat")){ next }
  Sleep_mental %>% select(c(column, "psqi_norm")) %>% drop_na() -> Test_corr
  Test_corr <- lapply(Test_corr, function(x) as.numeric(x))
  cor.test(as_vector(Test_corr[1]),as_vector(Test_corr[2]), method = "spearman") -> test
  n <- c(n,column)
  r <- c(r,test$estimate)
  p <- c(p,test$p.value)
}

print("Combine Tibble")
sleep_mental_correlation <- tibble(mental=n,r=r,p_value=p)
sleep_mental_correlation <- sleep_mental_correlation %>% mutate(p_adjust = p.adjust(p_value, "fdr"), sig = ifelse(r>0.1 & p_adjust< 0.05,T,F))

write.csv(sleep_mental_correlation, "output_mental/sleep_mental_corr_norm.csv")

# lm norm------------------------
n <- c()
p <- c()
r <- c()
e <- c()

print("PSQI - mental")
for (column in colnames(Sleep_mental)){
  if (column %in% c("DAG3_sampleID","psqi_norm","psqi_cat")){ next }
  Sleep_mental %>% select(c(column, "psqi_norm")) %>% drop_na() -> Test_lm
  model0 <- lm(psqi_norm~1,data=Test_lm)
  formula <- as.formula(paste(c( "psqi_norm~",column), collapse= "" ))
  model1<- lm(formula,data=Test_lm)
  result <- summary(model1)
  model <- anova(model0,model1)
  n <- c(n,column)
  r <- c(r,result$adj.r.squared)
  p <- c(p,model$`Pr(>F)`[2])
  e <- c(e,result$coefficients[2])
}

sleep_mental_lm <- tibble(mental=n,estimate=e,r=r,p_value=p)
sleep_mental_lm <- sleep_mental_lm %>% mutate(p_adjust = p.adjust(p_value, "fdr"), Sig = ifelse(r>0.01&p_adjust< 0.05,T,F))
write.csv(sleep_mental_lm, "output_mental/sleep_mental_lm_norm.csv")

# lm cat------------------------
Sleep_mental_num <- Sleep_mental
indx <- sapply(Sleep_mental_num, is.factor)
Sleep_mental_num[indx] <- lapply(Sleep_mental_num[indx], function(x) as.numeric(x))

n <- c()
p <- c()
r <- c()
e <- c()

print("PSQI - mental")
for (column in colnames(Sleep_mental_num)){
  if (column %in% c("DAG3_sampleID","psqi_norm","psqi_cat")){ next }
  Sleep_mental_num %>% select(c(column, "psqi_cat")) %>% drop_na() -> Test_lm
  model0 <- lm(psqi_cat~1,data=Test_lm)
  formula <- as.formula(paste(c( "psqi_cat~",column), collapse= "" ))
  model1<- lm(formula,data=Test_lm)
  result <- summary(model1)
  model <- anova(model0,model1)
  n <- c(n,column)
  r <- c(r,result$adj.r.squared)
  p <- c(p,model$`Pr(>F)`[2])
  e <- c(e,result$coefficients[2])
}

sleep_mental_lm <- tibble(mental=n,estimate=e,r=r,p_value=p)
sleep_mental_lm <- sleep_mental_lm %>% mutate(p_adjust = p.adjust(p_value, "fdr"), Sig = ifelse(r>0.01&p_adjust< 0.05,T,F))
write.csv(sleep_mental_lm, "output_mental/sleep_mental_lm_cat.csv")