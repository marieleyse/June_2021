setwd("Z:/AWAZANA/MELM/ALSPAC/data/")
ALSPAC = read.csv("Z:/AWAZANA/AWAZANA_LAB/ALSPAC/ALSPAC_all_combined_newPRS_Feb_2021Ruttercleanup.csv")
ALSPAC$Hyper_impulse_Rutter = -(ALSPAC$kj601 + ALSPAC$kj603)
ALSPAC$Inattention_Rutter = -(ALSPAC$kj618 + ALSPAC$kj629)
ALSPAC$Rutter = -(ALSPAC$kj601 + ALSPAC$kj603 + ALSPAC$kj618 + ALSPAC$kj629)

subsample_one_item_present1 <- ALSPAC[which((rowMeans(is.na(ALSPAC[c("kr850", "kr852", "kv8600", "jphyper", 
                                                                     "kqphyper", "n8365c", "tc4025c", "sa162b", "se162b",
                                                                     "Rutter",
                                                                     "f8bp036","f8ba036",
                                                                     "f8bp070",
                                                                     "f8ba070")]))) <= 13/14),]
dim(subsample_one_item_present1)

ADHDfactor.data1 <- subsample_one_item_present1[c("kr850", "kr852", "kv8600", "jphyper", 
                                                  "kqphyper", "n8365c", "tc4025c", "sa162b", "se162b",
                                                  "Rutter",
                                                  "f8bp036","f8ba036",
                                                  "f8bp070",
                                                  "f8ba070")]
str(ADHDfactor.data1)


# measurementInvariance(model= ADHD.fit2, 
#                       data = ADHDfactor.data2,
#                       estimator = "MLR", 
#                       missing = "FIML", 
#                       std.ov = TRUE, 
#                       std.lv = TRUE, 
#                       orthogonal = T, 
#                       verbose = T,  
#                       group = "gender_male")
# 
# capture.output(measurement_invariance, file = "model_differences_ADHD_model_fiml_*****.txt")
# o <- lavInspect(ADHD.fit2, what ='fit')
# capture.output(o, file = "fit_indices_ADHD_model_fiml_*****.txt")
# s <- summary(ADHD.fit2, standardized = T, fit.measures = T)
# capture.output(s, file = "model_output_ADHD_model_fiml_*****.txt")




sink('unifactor_adhd_model.txt')
ADHD.model1 <- '
ADHD =~ kr850 + kr852 + kv8600 + jphyper + kqphyper + n8365c + tc4025c + sa162b + se162b + 
Rutter +
f8bp036 + f8ba036 + f8bp070 + f8ba070'
ADHD.fit1 <- cfa(ADHD.model1, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)
o <- lavInspect(ADHD.fit1, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_onefactor.txt")
s <- summary(ADHD.fit1, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml_onefactor.txt")
sink()

sink("omega1.txt")
reliability(ADHD.fit1,
            omit.imps =c("no.conv", "no.se"))
sink()

ADHD_scores <- lavaan::lavPredict(ADHD.fit1)

### missForest imputation
# Impute
library(missForest)
set.seed(1234)
data_imp = missForest(ADHDfactor.data1, verbose = TRUE)$ximp

# Data transformations
data_imp = as.data.frame(apply(data_imp,2,as.numeric))

ALSPAC$kr850 = c(scale(ALSPAC$kr850))
ALSPAC$kr852 = c(scale(ALSPAC$kr852))
ALSPAC$kv8600 = c(scale(ALSPAC$kv8600))
ALSPAC$jphyper = c(scale(ALSPAC$jphyper))
ALSPAC$kqphyper = c(scale(ALSPAC$kqphyper))
ALSPAC$n8365c = c(scale(ALSPAC$n8365c))
ALSPAC$tc4025c = c(scale(ALSPAC$tc4025c))
ALSPAC$sa162b = c(scale(ALSPAC$sa162b))
ALSPAC$se162b = c(scale(ALSPAC$se162b))
ALSPAC$Rutter = c(scale(ALSPAC$Rutter))
ALSPAC$f8bp036 = c(scale(ALSPAC$f8bp036))
ALSPAC$f8ba036 = c(scale(ALSPAC$f8ba036))
ALSPAC$f8bp070 = c(scale(ALSPAC$f8bp070))
ALSPAC$f8ba070= c(scale(ALSPAC$f8ba070))

ADHDfactor.data1$kr850 = c(scale(ADHDfactor.data1$kr850))
ADHDfactor.data1$kr852 = c(scale(ADHDfactor.data1$kr852))
ADHDfactor.data1$kv8600 = c(scale(ADHDfactor.data1$kv8600))
ADHDfactor.data1$jphyper = c(scale(ADHDfactor.data1$jphyper))
ADHDfactor.data1$kqphyper = c(scale(ADHDfactor.data1$kqphyper))
ADHDfactor.data1$n8365c = c(scale(ADHDfactor.data1$n8365c))
ADHDfactor.data1$tc4025c = c(scale(ADHDfactor.data1$tc4025c))
ADHDfactor.data1$sa162b = c(scale(ADHDfactor.data1$sa162b))
ADHDfactor.data1$se162b = c(scale(ADHDfactor.data1$se162b))
ADHDfactor.data1$Rutter = c(scale(ADHDfactor.data1$Rutter))
ADHDfactor.data1$f8bp036 = c(scale(ADHDfactor.data1$f8bp036))
ADHDfactor.data1$f8ba036 = c(scale(ADHDfactor.data1$f8ba036))
ADHDfactor.data1$f8bp070 = c(scale(ADHDfactor.data1$f8bp070))
ADHDfactor.data1$f8ba070= c(scale(ADHDfactor.data1$f8ba070))

data_imp$kr850 = c(scale(data_imp$kr850))
data_imp$kr852 = c(scale(data_imp$kr852))
data_imp$kv8600 = c(scale(data_imp$kv8600))
data_imp$jphyper = c(scale(data_imp$jphyper))
data_imp$kqphyper = c(scale(data_imp$kqphyper))
data_imp$n8365c = c(scale(data_imp$n8365c))
data_imp$tc4025c = c(scale(data_imp$tc4025c))
data_imp$sa162b = c(scale(data_imp$sa162b))
data_imp$se162b = c(scale(data_imp$se162b))
data_imp$Rutter = c(scale(data_imp$Rutter))
data_imp$f8bp036 = c(scale(data_imp$f8bp036))
data_imp$f8ba036 = c(scale(data_imp$f8ba036))
data_imp$f8bp070 = c(scale(data_imp$f8bp070))
data_imp$f8ba070= c(scale(data_imp$f8ba070))

ALSPAC$ADHD_composite = ALSPAC$kr850 + ALSPAC$kr852 + ALSPAC$kv8600 + ALSPAC$jphyper + ALSPAC$kqphyper + ALSPAC$n8365c + ALSPAC$tc4025c + ALSPAC$sa162b + 
  ALSPAC$se162b + ALSPAC$Rutter + ALSPAC$f8bp036 + ALSPAC$f8ba036 + ALSPAC$f8bp070 + ALSPAC$f8ba070
  
ADHDfactor.data1$ADHD_composite2 = ADHDfactor.data1$kr850 + ADHDfactor.data1$kr852 + ADHDfactor.data1$kv8600 + 
  ADHDfactor.data1$jphyper + ADHDfactor.data1$kqphyper + ADHDfactor.data1$n8365c + ADHDfactor.data1$tc4025c + 
  ADHDfactor.data1$sa162b + ADHDfactor.data1$se162b + ADHDfactor.data1$Rutter + ADHDfactor.data1$f8bp036 + 
  ADHDfactor.data1$f8ba036 + ADHDfactor.data1$f8bp070 + ADHDfactor.data1$f8ba070
  
data_imp$ADHD_composite_imp = data_imp$kr850 + data_imp$kr852 + data_imp$kv8600 + data_imp$jphyper + data_imp$kqphyper + data_imp$n8365c + data_imp$tc4025c + data_imp$sa162b + 
  data_imp$se162b + data_imp$Rutter + data_imp$f8bp036 + data_imp$f8ba036 + data_imp$f8bp070 + data_imp$f8ba070

ADHDfactor.data1$ADHD_composite_wodawba =  
  ADHDfactor.data1$jphyper + ADHDfactor.data1$kqphyper + ADHDfactor.data1$n8365c + ADHDfactor.data1$tc4025c + 
  ADHDfactor.data1$sa162b + ADHDfactor.data1$se162b + ADHDfactor.data1$Rutter + ADHDfactor.data1$f8bp036 + 
  ADHDfactor.data1$f8ba036 + ADHDfactor.data1$f8bp070 + ADHDfactor.data1$f8ba070

data_imp$ADHD_composite_imp_wodawba = data_imp$jphyper + data_imp$kqphyper + data_imp$n8365c + data_imp$tc4025c + data_imp$sa162b + 
  data_imp$se162b + data_imp$Rutter + data_imp$f8bp036 + data_imp$f8ba036 + data_imp$f8bp070 + data_imp$f8ba070

NEW2 <- cbind(ADHDfactor.data1, data_imp)

NEW4 <- cbind(NEW2, ADHD_scores)

# NEW3 <- cbind(ALSPAC, NEW2)
# 
# NEW4 <- cbind(NEW3, ADHD_scores)

readr::write_csv(NEW4, file = "Z:/AWAZANA/MELM/ALSPAC/data/composite_test.csv")

NEW = read.csv("Z:/AWAZANA/MELM/ALSPAC/data/composite_test.csv")

png(filename = 'composite_test.png')

mydata.cor = stats::cor(NEW[, c("ADHD", "ADHD_composite2", "ADHD_composite_imp")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("ADHD", "ADHD_composite2", "ADHD_composite_imp")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.001, .01), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")

dev.off()

subsample_one_item_present2 <- ALSPAC[which((rowMeans(is.na(ALSPAC[c("kr850", "kr852", "kv8600", "jphyper", 
                                                                     "kqphyper", "n8365c", "tc4025c", "sa162b", "se162b",
                                                                     "Rutter",
                                                                     "f8bp036","f8ba036",
                                                                     "f8bp070",
                                                                     "f8ba070")]))) <= 13/14),]
dim(subsample_one_item_present2)

ADHDfactor.data2 <- subsample_one_item_present2[c("kr850", "kr852", "kv8600", "jphyper", 
                                                  "kqphyper", "n8365c", "tc4025c", "sa162b", "se162b",
                                                  "Rutter",
                                                  "f8bp036","f8ba036",
                                                  "f8bp070",
                                                  "f8ba070")]
str(ADHDfactor.data2)

sink('unifactor_adhd_model_wodawba.txt')
ADHD.model2 <- '
ADHD_wodawba =~ jphyper + kqphyper + n8365c + tc4025c + sa162b + se162b + 
Rutter +
f8bp036 + f8ba036 + f8bp070 + f8ba070'
ADHD.fit2 <- cfa(ADHD.model2, data = ADHDfactor.data2, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)
o <- lavInspect(ADHD.fit2, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_onefactor_wodawba.txt")
s <- summary(ADHD.fit2, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml_onefactor_wodawba.txt")
sink()

sink("omega1_wodawba.txt")
reliability(ADHD.fit2,
            omit.imps =c("no.conv", "no.se"))
sink()

sink('bifactor_adhd_model_adhd_informant_bi.txt')
ADHD.model_bi <- '
ADHD_bi =~ kr850 + kr852 + kv8600 + jphyper + kqphyper + n8365c + tc4025c + sa162b + se162b + 
Rutter +
f8bp036 + f8ba036 + f8bp070 + f8ba070
PARENT_bi =~ kr850 + kv8600 + jphyper + kqphyper + n8365c + tc4025c + Rutter 
TEACHER_bi =~ kr852 + sa162b + se162b
FIELD_bi =~ f8bp036 + f8ba036 + f8bp070 + f8ba070
'

ADHD.fit3 <- cfa(ADHD.model_bi, data = ADHDfactor.data2, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)
o <- lavInspect(ADHD.fit3, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_bifactor_adhd_informant_bi.txt")
s <- summary(ADHD.fit3, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml_bifactor_adhd_informant_bi.txt")
sink()

sink('omega_adhd_model_bi.txt')
semTools::reliability(ADHD.fit3,  
                      omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit3, "cov.lv")
sink()

semPaths(ADHD.fit3, "par", weighted = FALSE, nCharNodes = 7, shapeMan = "rectangle",
         sizeMan = 8, sizeMan2 = 5)
semPaths(ADHD.fit3, "std", weighted = FALSE, nCharNodes = 7, shapeMan = "rectangle",
         sizeMan = 8, sizeMan2 = 5)

#png(filename = 'adhd_latent_factor.png', width=1000, height=400)
pdf("adhd_latent_factor.pdf", width=1, height=1, compress = FALSE)
# cols <- wes_palette(name = "Zissou1", n = 5, type = "discrete")
# colorlist <- list(man = cols[2], lat = cols[5])
cols <- wes_palette(name = "FantasticFox1", n = 5, type = "discrete")
colorlist <- list(man = cols[3], lat = cols[4])
colorlist <- list(man = cols[3], lat = cols[4], man = cols[5])

#          color = colorlist, nCharNodes = 7, sizeMan = 7, sizeMan2 = 4,rotation=4,
#          shapeMan = "rectangle", intercepts =FALSE, residuals=FALSE, curve=0.5, fade = F, curvePivot = TRUE, arrows = 1,
# semPaths(ADHD.model_bi, what = "model", bifactor = "ADHD", whatLabels = "std", layout = "tree2", 
#          color = colorlist, nCharNodes = 0,rotation=2, shapeMan = "rectangle", sizeMan = 20, sizeMan2 = 2,
#          intercepts =FALSE, residuals=FALSE, curve=1, fade = F, curvePivot = TRUE, arrows = 1,
#          edge.label.cex=1, edge.label.position=0.65, mar = c(1,1,1,1), exoCov = FALSE, label.cex = 1, nodeLabels =c("Dominic ADHD 72m","CBCL Attention 48m", "    SDQ H 60m    ", " CBCL Attention 60m", "Conners ADHD 60m", "Conners ADHD 72m", "PAPA ADHD 72m",
#                                                                                                                     "    SDQ H 72m    ", "    SDQ H 60m    ", "Conners ADHD 60m  ", "Conners ADHD 72m", "    SDQ H 72m    ", "Conners ADHD 72m", "    SDQ H 72m    ", "ADHD", "Mother", "Father", "Teacher"))
semPaths(ADHD.fit3, what = "model", bifactor = "ADHD_bi", whatLabels = "std", layout = "tree2", 
         color = colorlist, nCharNodes = 0,rotation=2, shapeMan = "rectangle", sizeMan = 20, sizeMan2 = 2,
         intercepts =FALSE, residuals=FALSE, curve=1, fade = F, curvePivot = TRUE, arrows = 1,
         edge.label.cex=1, edge.label.position=0.65, mar = c(1,1,1,1), exoCov = FALSE, label.cex = 1)
dev.off()

pdf("adhd_latent_factor_numbered1.pdf", width=1, height=1, compress = FALSE)
cols <- wes_palette(name = "FantasticFox1", n = 5, type = "discrete")
colorlist <- list(man = cols[3], lat = cols[4])
colorlist <- list(man = cols[3], lat = cols[4], man = cols[5])
semPaths(ADHD.fit3, what = "model", bifactor = "ADHD_bi", whatLabels = "std", layout = "tree3", 
         color = colorlist, nCharNodes = 0,rotation=2, shapeMan = "rectangle", sizeMan = 20, sizeMan2 = 2,
         intercepts =FALSE, residuals=FALSE, curve=1, fade = F, curvePivot = TRUE, arrows = 1,
         edge.label.cex=1, edge.label.position=0.85, mar = c(1,1,1,1), exoCov = FALSE, label.cex = 1, nodeLabels =c("1","2", " 3 ", " 4", "5", "6", "7",
                       "8", " 9", "10", "11", "12", "13", "14 ", "15", " 16 ","17", "18"))
# 
dev.off()

pdf("adhd_latent_factor_labels.pdf", width=1, height=1, compress = FALSE)
cols <- wes_palette(name = "FantasticFox1", n = 5, type = "discrete")
colorlist <- list(man = cols[3], lat = cols[4])
colorlist <- list(man = cols[3], lat = cols[4], man = cols[5])
semPaths(ADHD.fit3, what = "model", bifactor = "ADHD_bi", whatLabels = "std", layout = "tree3", 
         color = colorlist, nCharNodes = 0,rotation=2, shapeMan = "rectangle", sizeMan = 26, sizeMan2 = 2,
         intercepts =FALSE, residuals=FALSE, curve=1, fade = F, curvePivot = TRUE, arrows = 1,
         edge.label.cex=1, edge.label.position=0.45, mar = c(1,1,1,1), exoCov = FALSE, label.cex = 1, nodeLabels =c("DAWBA parent 7yr", "DAWBA parent 10yr","SDQ 3yr 11m", "SDQ 6yr 9m", "SDQ 8yr 1m", "SDQ 9yr 7m", 
                                                                                                                    "Rutter 3yr 6m","DAWBA teacher 7yr",  "SDQ teacher 7yr", "SDQ teacher 10yr","Posting impulsivity distractibility 8yr 6m","Activities impulsivity distractibility 8yr 6m",
                                                                                                                    "Fidgety Posting 8yr 6m","Fidgety Activities 8yr 6m","ADHD", "Parent", "Teacher", "Field"))
dev.off()

# "kr850", "kv8600", "jphyper", "kqphyper", "n8365c", "tc4025c", "Rutter","kr852", "sa162b", "se162b", "f8bp036","f8ba036","f8bp070","f8ba070"
# 
# "DAWBA_parent_7", "DAWBA_parent_10","SDQ_47m", "SDQ_81m", "SDQ_97m", "SDQ_115m", 
# "Rutter","DAWBA_teacher_7",  "SDQ_teacher_7", "SDQ_teacher_10","Posting_impulsivity_distractibility","Activities_impulsivity_distractibility",
# "Fidgety_Posting","Fidgety_activities"

sink('bifactor_adhd_model_adhd_informant_hierarchy.txt')
ADHD.model_hierarchy <- '
ADHD_hierarchy =~ PARENT_hierarchy + TEACHER_hierarchy + FIELD_hierarchy
PARENT_hierarchy =~ kr850 + kv8600 + jphyper + kqphyper + n8365c + tc4025c + Rutter 
TEACHER_hierarchy =~ kr852 + sa162b + se162b
FIELD_hierarchy =~ f8bp036 + f8ba036 + f8bp070 + f8ba070
'

ADHD.fit4 <- cfa(ADHD.model_hierarchy, data = ADHDfactor.data2, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)
o <- lavInspect(ADHD.fit4, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_bifactor_adhd_informant_hierarchy.txt")
s <- summary(ADHD.fit4, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml_bifactor_adhd_informant_hierarchy.txt")
sink()

sink('bifactor_adhd_model_adhd_informant_wodawba.txt')
ADHD.model_bi_wodawba <- '
ADHD_bi_wodawba =~ kr850 + kr852 + kv8600 + jphyper + kqphyper + n8365c + tc4025c + sa162b + se162b + 
Rutter +
f8bp036 + f8ba036 + f8bp070 + f8ba070
PARENT_bi_wodawba =~ jphyper + kqphyper + n8365c + tc4025c + Rutter 
TEACHER_bi_wodawba =~ a*sa162b + a*se162b
FIELD_bi_wodawba =~ f8bp036 + f8ba036 + f8bp070 + f8ba070
'
ADHD.fit5 <- cfa(ADHD.model_bi_wodawba, data = ADHDfactor.data2, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)
o <- lavInspect(ADHD.fit5, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_bifactor_adhd_informant.txt")
s <- summary(ADHD.fit5, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml_bifactor_adhd_informant.txt")
sink()

sink('bifactor_adhd_model_adhd_informant_hierarchy_wodawba.txt')
ADHD.model_hierarchical_wodawba <- '
ADHD_hierarchical_wodawba =~ PARENT_hierarchical_wodawba + TEACHER_hierarchical_wodawba + FIELD_hierarchical_wodawba
PARENT_hierarchical_wodawba =~ jphyper + kqphyper + n8365c + tc4025c + Rutter 
TEACHER_hierarchical_wodawba =~ a*sa162b + a*se162b
FIELD_hierarchical_wodawba =~ f8bp036 + f8ba036 + f8bp070 + f8ba070
'
ADHD.fit6 <- cfa(ADHD.model_hierarchical_wodawba, data = ADHDfactor.data2, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)
o <- lavInspect(ADHD.fit6, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_bifactor_adhd_informant_hierarchy_hierarchical_wodawba.txt")
s <- summary(ADHD.fit6, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml_bifactor_adhd_informant_hierarchy_hierarchical_wodawba.txt")
sink()

ADHD_scores_2 <- lavaan::lavPredict(ADHD.fit2)

ADHD_scores_3 <- lavaan::lavPredict(ADHD.fit3)

ADHD_scores_4 <- lavaan::lavPredict(ADHD.fit4)

ADHD_scores_5 <- lavaan::lavPredict(ADHD.fit5)

ADHD_scores_6 <- lavaan::lavPredict(ADHD.fit6)

NEW2 <- cbind(NEW, ADHD_scores_2)

NEW3 <- cbind(NEW2, ADHD_scores_3)

NEW4 <- cbind(NEW3, ADHD_scores_4)

NEW5 <- cbind(NEW4, ADHD_scores_5)

NEW6 <- cbind(NEW5, ADHD_scores_6)

readr::write_csv(NEW6, file = "Z:/AWAZANA/MELM/ALSPAC/data/composite_test2.csv")

NEW = read.csv("Z:/AWAZANA/MELM/ALSPAC/data/composite_test2.csv")

# png(filename = 'composite_test2.png')
pdf("composite_test2.pdf", width=5, height=5, compress=FALSE)

# ADHD_wodawba ADHD_bi_wodawba PARENT_bi_wodawba TEACHER_bi_wodawba FIELD_bi_wodawba
# ADHD_bi PARENT_bi TEACHER_bi FIELD_bi
# ADHD_hierarchy PARENT_hierarchy TEACHER_hierarchy FIELD_hierarchy
# ADHD_hierarchical_wodawba PARENT_hierarchical_wodawba TEACHER_hierarchical_wodawba FIELD_hierarchical_wodawba

mydata.cor2 = stats::cor(NEW[, c("ADHD", "ADHD_composite2", "ADHD_composite_imp", "ADHD_composite_imp_wodawba", "ADHD_wodawba", "ADHD_bi_wodawba", "PARENT_bi_wodawba", "TEACHER_bi_wodawba", "FIELD_bi_wodawba", 
                                "ADHD_bi", "PARENT_bi", "TEACHER_bi", "FIELD_bi", 
                                "ADHD_hierarchy", "PARENT_hierarchy", "TEACHER_hierarchy", "FIELD_hierarchy", 
                                "ADHD_hierarchical_wodawba", "PARENT_hierarchical_wodawba", "TEACHER_hierarchical_wodawba", "FIELD_hierarchical_wodawba", 
                                "kr850", "kr852", "kv8600")],  method ="kendall",  use = "pairwise")

res2 <- corrplot::cor.mtest(NEW[, c("ADHD", "ADHD_composite2", "ADHD_composite_imp", "ADHD_composite_imp_wodawba", "ADHD_wodawba", "ADHD_bi_wodawba", "PARENT_bi_wodawba", "TEACHER_bi_wodawba", "FIELD_bi_wodawba", 
                                    "ADHD_bi", "PARENT_bi", "TEACHER_bi", "FIELD_bi", 
                                    "ADHD_hierarchy", "PARENT_hierarchy", "TEACHER_hierarchy", "FIELD_hierarchy", 
                                    "ADHD_hierarchical_wodawba", "PARENT_hierarchical_wodawba", "TEACHER_hierarchical_wodawba", "FIELD_hierarchical_wodawba",
                                    "ADHD_wodawba", "kr850", "kr852", "kv8600")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

corrplot::corrplot(mydata.cor2, p.mat = res2$p, method = "color", type = "upper",
                   sig.level = c(.001, .01), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")

dev.off()


subsample_one_item_present3 <- ALSPAC[which((rowMeans(is.na(ALSPAC[c("ID","kr850", "kr852", "kv8600", "jphyper", 
                                                                     "kqphyper", "n8365c", "tc4025c", "sa162b", "se162b",
                                                                     "Rutter",
                                                                     "f8bp036","f8ba036",
                                                                     "f8bp070",
                                                                     "f8ba070")]))) <= 13/15),]
dim(subsample_one_item_present3)

ADHDfactor.data3 <- subsample_one_item_present3[c("ID","kr850", "kr852", "kv8600", "jphyper", 
                                                  "kqphyper", "n8365c", "tc4025c", "sa162b", "se162b",
                                                  "Rutter",
                                                  "f8bp036","f8ba036",
                                                  "f8bp070",
                                                  "f8ba070")]
str(ADHDfactor.data3)

# "kwphyper", "tcphyper"
# "fjc1200", "fjc1201", 

#tb8600 6-band adhd at 13
#fh6870 6-band adhd at 15

subsample_one_item_present_old <- ALSPAC[which((rowMeans(is.na(ALSPAC[c("ID","fh6870","tb8600", "kj645")]))) <= 2/4),]
dim(subsample_one_item_present_old)

old <- subsample_one_item_present_old[c("ID","fh6870","tb8600", "kj645")]
str(old)

NEW2 <- cbind(NEW, ADHDfactor.data3)

NEW3 <- merge(NEW2, old, by = "ID", all=TRUE)

readr::write_csv(NEW3, file = "Z:/AWAZANA/MELM/ALSPAC/data/predict_older_dawba.csv")

NEW = read.csv("Z:/AWAZANA/MELM/ALSPAC/data/predict_older_dawba.csv")

pdf("correlation_olderbawba.pdf", width=5, height=5, compress=FALSE)

mydata.cor3 = stats::cor(NEW[, c("fh6870","tb8600", "ADHD", "ADHD_composite2", "ADHD_composite_imp", "ADHD_composite_imp_wodawba", "ADHD_wodawba", "ADHD_bi_wodawba", "PARENT_bi_wodawba", "TEACHER_bi_wodawba", "FIELD_bi_wodawba", 
                                "ADHD_bi", "PARENT_bi", "TEACHER_bi", "FIELD_bi", 
                                "ADHD_hierarchy", "PARENT_hierarchy", "TEACHER_hierarchy", "FIELD_hierarchy", 
                                "ADHD_hierarchical_wodawba", "PARENT_hierarchical_wodawba", "TEACHER_hierarchical_wodawba", "FIELD_hierarchical_wodawba", 
                                "kr850", "kr852", "kv8600")],  method ="kendall",  use = "pairwise")


res3 <- corrplot::cor.mtest(NEW[, c("fh6870","tb8600", 
                                    "ADHD", "ADHD_composite2", 
                                    "ADHD_composite_imp", "ADHD_composite_imp_wodawba", 
                                    "ADHD_wodawba", "ADHD_bi_wodawba", 
                                    "PARENT_bi_wodawba", "TEACHER_bi_wodawba", 
                                    "FIELD_bi_wodawba", "ADHD_bi", 
                                    "PARENT_bi", "TEACHER_bi", 
                                    "FIELD_bi", "ADHD_hierarchy", 
                                    "PARENT_hierarchy", "TEACHER_hierarchy", 
                                    "FIELD_hierarchy", "ADHD_hierarchical_wodawba", 
                                    "PARENT_hierarchical_wodawba", "TEACHER_hierarchical_wodawba", 
                                    "FIELD_hierarchical_wodawba", "kr850", 
                                    "kr852", "kv8600")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

colnames(mydata.cor3) <-c("DAWBA clinic 15","DAWBA parent 13", 
                         "ADHD", "ADHD composite score", 
                         "ADHD composite score imp", "ADHDcomposite imp wo DAWBA", 
                         "ADHD wo DAWBA", "ADHD bi wo DAWBA", 
                         "PARENT bi wo DAWBA", "TEACHER bi wo DAWBA", 
                         "FIELD bi wo DAWBA", "ADHD bi", 
                         "PARENT bi", "TEACHER bi", 
                         "FIELD bi", "ADHD hierarchy", 
                         "PARENT hierarchy", "TEACHER hierarchy", 
                         "FIELD hierarchy", "ADHD hierarchical wo DAWBA", 
                         "PARENT hierarchical wo DAWBA", "TEACHER hierarchical wo DAWBA", 
                         "FIELD hierarchical wo DAWBA", "DAWBA parent 7", 
                         "DAWBA teacher 7", "DAWBA parent 10")

rownames(mydata.cor3) <-c("DAWBA clinic 15","DAWBA parent 13", 
                         "ADHD", "ADHD composite score", 
                         "ADHD composite score imp", "ADHDcomposite imp wo DAWBA", 
                         "ADHD wo DAWBA", "ADHD bi wo DAWBA", 
                         "PARENT bi wo DAWBA", "TEACHER bi wo DAWBA", 
                         "FIELD bi wo DAWBA", "ADHD bi", 
                         "PARENT bi", "TEACHER bi", 
                         "FIELD bi", "ADHD hierarchy", 
                         "PARENT hierarchy", "TEACHER hierarchy", 
                         "FIELD hierarchy", "ADHD hierarchical wo DAWBA", 
                         "PARENT hierarchical wo DAWBA", "TEACHER hierarchical wo DAWBA", 
                         "FIELD hierarchical wo DAWBA", "DAWBA parent 7", 
                         "DAWBA teacher 7", "DAWBA parent 10")


# pdf("correlation_olderbawba.pdf", width=5, height=5, compress=FALSE)

corrplot::corrplot(mydata.cor3, p.mat = res3$p, method = "color", type = "upper",
                   sig.level = c(.001, .01), pch.cex = .6, tl.cex = .4,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")

dev.off()

pdf("correlation_olderbawba_fix.pdf", width=5, height=5, compress=FALSE)

corrplot::corrplot(mydata.cor3, p.mat = res3$p, method = "color", type = "full",
         sig.level = c(.001, .01), pch.cex = .6, tl.cex = .4,  tl.srt = 45,
         insig = "label_sig", pch.col = "black", tl.col = "black")

dev.off()