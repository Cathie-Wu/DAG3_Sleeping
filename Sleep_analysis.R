library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(vegan)
library(ape) 


# Load and clean ----------------------------------------------------------

load("metaphlan3_cleaned.RData")
metadata <- read.csv("DAG3_metadata_merged_ready_v27.csv")
sleep <- read_tsv("Sleep/LLs_sleep_data.tsv")

# Sleep phenotype
sleep_2 <- sleep %>% 
  select(DAG3_sampleID,psqi_globalscore_adu_c_1,psqi_qualityscore_adu_c_1) %>% 
  filter(DAG3_sampleID!="NA"
         & psqi_globalscore_adu_c_1!="NA"
         & psqi_globalscore_adu_c_1!= "$6") %>% 
  mutate(psqi = as.numeric(psqi_globalscore_adu_c_1), .keep= "unused") %>% 
  mutate(psqi_bi = as.numeric(psqi_qualityscore_adu_c_1), .keep= "unused") %>% 
  mutate(psqi_norm = qnorm((rank(psqi,na.last="keep")-0.5)/sum(!is.na(psqi))))

# +metadata
sleep_2 %>% inner_join(metadata,by="DAG3_sampleID")-> Sleep_meta

# microbiome
colnames(taxa_matrix)[grepl("\\[", colnames(taxa_matrix))] <- "s__Collinsella_massiliensis"
colnames(df_species)[grepl("\\[", colnames(df_species))] <- "s__Collinsella_massiliensis"

df_species <-df_species %>% mutate(DAG3_sampleID=gsub("_metaphlan", "", df_species$ID)) %>% select(-1) %>% semi_join(Sleep_meta,by="DAG3_sampleID") 
df_species <-df_species %>% select("DAG3_sampleID",colnames(df_species))

taxa_matrix <-taxa_matrix %>% mutate(DAG3_sampleID=gsub("_metaphlan", "", taxa_matrix$ID)) %>% select(-1) %>% semi_join(Sleep_meta,by="DAG3_sampleID") 
taxa_matrix <-taxa_matrix %>% select("DAG3_sampleID",colnames(taxa_matrix))

save(df_species,taxa_matrix,Sleep_meta,file="Complete.RData")

# PSQI - phenotype ----------------------------------------------------------
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


n <- c()
p <- c()
r <- c()
e <- c()

print("PSQI - phenotype")
for (column in colnames(Sleep_meta)){
  if (column %in% c("DAG3_sampleID","psqi_norm","psqi")){ next }
  Sleep_meta %>% select(c(column, "psqi_norm")) %>% drop_na() -> Test_lm
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

sleep_norm_lm <- tibble(phenotype=n,estimate=e,r=r,p_value=p)
sleep_norm_lm <- sleep_norm_lm %>% mutate(p_adjust = p.adjust(p_value, "fdr"), use = ifelse(r>0.01&p_adjust< 0.05,T,F))
write_tsv(sleep_norm_lm,file = "output/sleep_norm_lm.tsv")

# Alpha  -----------------------------------------------------------------
alpha <- df_species
alpha$alpha <- diversity(alpha[,2:dim(alpha)[2]], index = "shannon")

Sleep_meta %>% select(DAG3_sampleID,psqi,psqi_bi,psqi_norm,ANTHRO.AGE,ANTHRO.BMI,ANTHRO.Sex,META.BATCH,META.DNA.conc.ng.ul,META.DNA.postclean.reads,META.POOP.BristolMean,META.POOP.Freq,META.POOP.COLLECTION_SEASON)-> Sleep_meta_2
Sleep_meta_2[complete.cases(Sleep_meta_2),]-> Sleep_meta_2

alpha %>% select("DAG3_sampleID","alpha") %>%  inner_join(Sleep_meta_2,by="DAG3_sampleID")-> alpha

model <- lm(alpha ~ psqi + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq + META.POOP.COLLECTION_SEASON,data=alpha)
summary(model)
as.data.frame(summary(model)$coefficients) -> model_results

model_results <- model_results %>% rownames_to_column() %>% as_tibble()
write_tsv(x = model_results,file = "output/psqi_alpha_model.tsv")

# Taxa --------------------------------------------------------------------

## CLR transformation ------------------------------------------------------

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
final_data <- left_join(Sleep_meta,data_transformed)


## taxa ~ psqi -------------------------------------------------------------
n <- c()
e <- c()
s <- c()
prev <- c()
r <- c()
p <- c()

for (col in filter){
  formula <- as.formula(paste(col,"~psqi_norm + ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.BATCH + META.DNA.conc.ng.ul + META.DNA.postclean.reads + META.POOP.BristolMean + META.POOP.Freq + META.POOP.COLLECTION_SEASON", collapse= "" ))
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
write_tsv(df_result,file = "output/taxa_psqi.tsv")
