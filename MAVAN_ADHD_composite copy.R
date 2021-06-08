setwd("/Users/Marie-Elyse/Downloads")
# NEW = read.csv("MAVAN_48M_and_up_jun2020.csv")
# NEW = read.csv("MAVAN_48M_and_up_jun2020_new.csv")
# NEW = read.csv("NOV302020.csv", header=T)
NEW = read.csv("FEB2021.csv")

#Subset data according to who has at least one item present based on the updated (e.g. parecellated) number of items.
subsample_one_item_present1 <- NEW[which((rowMeans(is.na(NEW[c("CBCL48_sc6raw.x", "SDQ60_mother_hyperactivity","SDQ60_father_hyperactivity", "CBCL60_sc6raw.x","conners_mother_adhd_score.60m", "conners_father_adhd_score.60m", "conners_mother_adhd_score.72m", "conners_father_adhd_score.72m", "conners_teacher_adhd_score.72m", "PAPA_p4nadhd", "Dominic72_ADHD","SDQ72_mother_hyperactivity", "SDQ72_father_hyperactivity", "SDQ72_teacher_hyperactivity")]))) <= 13/14),]
dim(subsample_one_item_present1)

ADHDfactor.data1 <- subsample_one_item_present1[c("CBCL48_sc6raw.x", "SDQ60_mother_hyperactivity","SDQ60_father_hyperactivity", "CBCL60_sc6raw.x","conners_mother_adhd_score.60m", "conners_father_adhd_score.60m", "conners_mother_adhd_score.72m", "conners_father_adhd_score.72m", "conners_teacher_adhd_score.72m", "PAPA_p4nadhd", "Dominic72_ADHD","SDQ72_mother_hyperactivity", "SDQ72_father_hyperactivity", "SDQ72_teacher_hyperactivity")]
str(ADHDfactor.data1)

### missForest imputation
# Impute
library(missForest)
set.seed(1234)
data_imp = missForest(ADHDfactor.data1, verbose = TRUE)$ximp

# Data transformations
data_imp = as.data.frame(apply(data_imp,2,as.numeric))

NEW$std_CBCL48_sc6raw.x = c(scale(NEW$CBCL48_sc6raw.x))
NEW$std_SDQ60_mother_hyperactivity = c(scale(NEW$SDQ60_mother_hyperactivity))
NEW$std_SDQ60_father_hyperactivity = c(scale(NEW$SDQ60_father_hyperactivity))
NEW$std_CBCL60_sc6raw.x = c(scale(NEW$CBCL60_sc6raw.x))
NEW$std_conners_mother_adhd_score.60m = c(scale(NEW$conners_mother_adhd_score.60m))
NEW$std_conners_father_adhd_score.60m = c(scale(NEW$conners_father_adhd_score.60m))
NEW$std_conners_mother_adhd_score.72m = c(scale(NEW$conners_mother_adhd_score.72m))
NEW$std_conners_father_adhd_score.72m = c(scale(NEW$conners_father_adhd_score.72m))
NEW$std_conners_teacher_adhd_score.72m = c(scale(NEW$conners_teacher_adhd_score.72m))
NEW$std_PAPA_p4nadhd = c(scale(NEW$PAPA_p4nadhd))
NEW$std_Dominic72_ADHD = c(scale(NEW$Dominic72_ADHD))
NEW$std_SDQ72_mother_hyperactivity = c(scale(NEW$SDQ72_mother_hyperactivity))
NEW$std_SDQ72_father_hyperactivity = c(scale(NEW$SDQ72_father_hyperactivity))
NEW$std_SDQ72_teacher_hyperactivity = c(scale(NEW$SDQ72_teacher_hyperactivity))

ADHDfactor.data1$std_CBCL48_sc6raw.x = c(scale(ADHDfactor.data1$CBCL48_sc6raw.x))
ADHDfactor.data1$std_SDQ60_mother_hyperactivity = c(scale(data_imp$SDQ60_mother_hyperactivity))
ADHDfactor.data1$std_SDQ60_father_hyperactivity = c(scale(ADHDfactor.data1$SDQ60_father_hyperactivity))
ADHDfactor.data1$std_CBCL60_sc6raw.x = c(scale(ADHDfactor.data1$CBCL60_sc6raw.x))
ADHDfactor.data1$std_conners_mother_adhd_score.60m = c(scale(ADHDfactor.data1$conners_mother_adhd_score.60m))
ADHDfactor.data1$std_conners_father_adhd_score.60m = c(scale(ADHDfactor.data1$conners_father_adhd_score.60m))
ADHDfactor.data1$std_conners_mother_adhd_score.72m = c(scale(ADHDfactor.data1$conners_mother_adhd_score.72m))
ADHDfactor.data1$std_conners_father_adhd_score.72m = c(scale(ADHDfactor.data1$conners_father_adhd_score.72m))
ADHDfactor.data1$std_conners_teacher_adhd_score.72m = c(scale(ADHDfactor.data1$conners_teacher_adhd_score.72m))
ADHDfactor.data1$std_PAPA_p4nadhd = c(scale(ADHDfactor.data1$PAPA_p4nadhd))
ADHDfactor.data1$std_Dominic72_ADHD = c(scale(ADHDfactor.data1$Dominic72_ADHD))
ADHDfactor.data1$std_SDQ72_mother_hyperactivity = c(scale(ADHDfactor.data1$SDQ72_mother_hyperactivity))
ADHDfactor.data1$std_SDQ72_father_hyperactivity = c(scale(ADHDfactor.data1$SDQ72_father_hyperactivity))
ADHDfactor.data1$std_SDQ72_teacher_hyperactivity = c(scale(ADHDfactor.data1$SDQ72_teacher_hyperactivity))

data_imp$std_CBCL48_sc6raw.x = c(scale(data_imp$CBCL48_sc6raw.x))
data_imp$std_SDQ60_mother_hyperactivity = c(scale(data_imp$SDQ60_mother_hyperactivity))
data_imp$std_SDQ60_father_hyperactivity = c(scale(data_imp$SDQ60_father_hyperactivity))
data_imp$std_CBCL60_sc6raw.x = c(scale(data_imp$CBCL60_sc6raw.x))
data_imp$std_conners_mother_adhd_score.60m = c(scale(data_imp$conners_mother_adhd_score.60m))
data_imp$std_conners_father_adhd_score.60m = c(scale(data_imp$conners_father_adhd_score.60m))
data_imp$std_conners_mother_adhd_score.72m = c(scale(data_imp$conners_mother_adhd_score.72m))
data_imp$std_conners_father_adhd_score.72m = c(scale(data_imp$conners_father_adhd_score.72m))
data_imp$std_conners_teacher_adhd_score.72m = c(scale(data_imp$conners_teacher_adhd_score.72m))
data_imp$std_PAPA_p4nadhd = c(scale(data_imp$PAPA_p4nadhd))
data_imp$std_Dominic72_ADHD = c(scale(data_imp$Dominic72_ADHD))
data_imp$std_SDQ72_mother_hyperactivity = c(scale(data_imp$SDQ72_mother_hyperactivity))
data_imp$std_SDQ72_father_hyperactivity = c(scale(data_imp$SDQ72_father_hyperactivity))
data_imp$std_SDQ72_teacher_hyperactivity = c(scale(data_imp$SDQ72_teacher_hyperactivity))

NEW$ADHD_composite = NEW$CBCL48_sc6raw.x + NEW$SDQ60_mother_hyperactivity + NEW$SDQ60_father_hyperactivity + NEW$CBCL60_sc6raw.x + NEW$conners_mother_adhd_score.60m + NEW$conners_father_adhd_score.60m + NEW$conners_mother_adhd_score.72m + NEW$conners_father_adhd_score.72m + NEW$conners_teacher_adhd_score.72m + NEW$PAPA_p4nadhd + NEW$Dominic72_ADHD + NEW$SDQ72_mother_hyperactivity + NEW$SDQ72_father_hyperactivity + NEW$SDQ72_teacher_hyperactivity

ADHDfactor.data1$ADHD_composite2 = ADHDfactor.data1$CBCL48_sc6raw.x + ADHDfactor.data1$SDQ60_mother_hyperactivity + ADHDfactor.data1$SDQ60_father_hyperactivity + ADHDfactor.data1$CBCL60_sc6raw.x + ADHDfactor.data1$conners_mother_adhd_score.60m + ADHDfactor.data1$conners_father_adhd_score.60m + ADHDfactor.data1$conners_mother_adhd_score.72m + ADHDfactor.data1$conners_father_adhd_score.72m + ADHDfactor.data1$conners_teacher_adhd_score.72m + ADHDfactor.data1$PAPA_p4nadhd + ADHDfactor.data1$Dominic72_ADHD + ADHDfactor.data1$SDQ72_mother_hyperactivity + ADHDfactor.data1$SDQ72_father_hyperactivity + ADHDfactor.data1$SDQ72_teacher_hyperactivity

data_imp$ADHD_composite_imp = data_imp$std_CBCL48_sc6raw.x + data_imp$std_SDQ60_mother_hyperactivity + data_imp$std_SDQ60_father_hyperactivity + data_imp$std_CBCL60_sc6raw.x + data_imp$std_conners_mother_adhd_score.60m + data_imp$std_conners_father_adhd_score.60m + data_imp$std_conners_mother_adhd_score.72m + data_imp$std_conners_father_adhd_score.72m + data_imp$std_conners_teacher_adhd_score.72m + data_imp$std_PAPA_p4nadhd + data_imp$std_Dominic72_ADHD + data_imp$std_SDQ72_mother_hyperactivity + data_imp$std_SDQ72_father_hyperactivity + data_imp$std_SDQ72_teacher_hyperactivity

ADHD.model1 <- '
ADHD_factor =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + PAPA_p4nadhd + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity'

ADHD.fit1 <- lavaan::cfa(ADHD.model1, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)

ADHD_scores <- lavaan::lavPredict(ADHD.fit1)
# lavPredict(object, newdata = NULL, type = "lv", method = "EBM",
#            se = "none", acov = "none", label = TRUE, fsm = FALSE, 
#            append.data = FALSE, assemble = FALSE,
#            level = 1L, optim.method = "bfgs", ETA = NULL)

ADHD.model_woPAPA <- '
ADHD_factor_woPAPA =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity'

ADHD.fit_woPAPA <- lavaan::cfa(ADHD.model_woPAPA, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)

ADHD_scores_woPAPA <- lavaan::lavPredict(ADHD.fit_woPAPA)

NEW2 <- cbind(ADHDfactor.data1, data_imp)

NEW3 <- cbind(NEW, NEW2)

NEW4 <- cbind(NEW3, ADHD_scores)

NEW5 <- cbind(NEW4, ADHD_scores_woPAPA)

readr::write_csv(NEW5, file = "/Users/Marie-Elyse/Downloads/composite_test.csv")

NEW = read.csv("composite_test.csv")

pdf("composite_test.pdf", width=5, height=5, compress=FALSE)

mydata.cor = stats::cor(NEW[, c("ADHD",	"Mother",	"Father",	"Teacher", "ADHD_factor", "ADHD_factor_woPAPA",
                                "ADHD_composite",	"ADHD_composite2", "ADHD_composite_imp")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("ADHD",	"Mother",	"Father",	"Teacher", "ADHD_factor","ADHD_factor_woPAPA",
                                    "ADHD_composite",	"ADHD_composite2", "ADHD_composite_imp")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.01, .001), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

pdf("composite_test2.pdf", width=5, height=5, compress=FALSE)

mydata.cor = stats::cor(NEW[, c("PAPA_p4_adhd", "PAPA_p4nadhd", "ADHD",	"Mother",	"Father",	"Teacher", "ADHD_factor", "ADHD_factor_woPAPA",
                        "ADHD_composite2", "ADHD_composite_imp")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("PAPA_p4_adhd", "PAPA_p4nadhd", "ADHD",	"Mother",	"Father",	"Teacher", "ADHD_factor", "ADHD_factor_woPAPA",
                         "ADHD_composite2", "ADHD_composite_imp")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
         sig.level = c(.05, .01, .001), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
         insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

pdf("composite_test3.pdf", width=5, height=5, compress=FALSE)

mydata.cor = stats::cor(NEW[, c("PAPA_p4nadhd", "ADHD",	"Mother",	"Father",	"Teacher", "ADHD_factor", "ADHD_factor_woPAPA",
                                "ADHD_composite2", "ADHD_composite_imp")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("PAPA_p4nadhd", "ADHD",	"Mother",	"Father",	"Teacher", "ADHD_factor","ADHD_factor_woPAPA",
                                    "ADHD_composite2", "ADHD_composite_imp")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

colnames(mydata.cor) <-c("PAPA", "ADHD",	"Mother",	"Father",	"Teacher", "ADHD unifactor","ADHD factor wo PAPA",
                         "ADHD composite score", "ADHD composite score imp")

rownames(mydata.cor) <-c("PAPA", "ADHD",	"Mother",	"Father",	"Teacher", "ADHD unifactor","ADHD factor wo PAPA",
                         "ADHD composite score", "ADHD composite score imp")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.05, .01, .001), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

subsample_one_item_present1 <- NEW[which((rowMeans(is.na(NEW[c("CBCL48_sc6raw.x", "SDQ60_mother_hyperactivity","SDQ60_father_hyperactivity", "CBCL60_sc6raw.x","conners_mother_adhd_score.60m", "conners_father_adhd_score.60m", "conners_mother_adhd_score.72m", "conners_father_adhd_score.72m", "conners_teacher_adhd_score.72m", "PAPA_p4nadhd", "Dominic72_ADHD","SDQ72_mother_hyperactivity", "SDQ72_father_hyperactivity", "SDQ72_teacher_hyperactivity")]))) <= 13/14),]
dim(subsample_one_item_present1)

ADHDfactor.data1 <- subsample_one_item_present1[c("CBCL48_sc6raw.x", "SDQ60_mother_hyperactivity","SDQ60_father_hyperactivity", "CBCL60_sc6raw.x","conners_mother_adhd_score.60m", "conners_father_adhd_score.60m", "conners_mother_adhd_score.72m", "conners_father_adhd_score.72m", "conners_teacher_adhd_score.72m", "PAPA_p4nadhd", "Dominic72_ADHD","SDQ72_mother_hyperactivity", "SDQ72_father_hyperactivity", "SDQ72_teacher_hyperactivity")]
str(ADHDfactor.data1)

sink('bi_adhd_model.txt')
ADHD.model1 <- '
ADHD =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + PAPA_p4nadhd + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity
Mother =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_mother_adhd_score.72m + PAPA_p4nadhd + SDQ72_mother_hyperactivity
Father =~ SDQ60_father_hyperactivity + conners_father_adhd_score.60m + conners_father_adhd_score.72m + SDQ72_father_hyperactivity 
Teacher =~ conners_teacher_adhd_score.72m + SDQ72_teacher_hyperactivity'

ADHD.fit1 <- lavaan::cfa(ADHD.model1, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)

o <- lavInspect(ADHD.fit1, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_ADHD_rater.txt")
s <- summary(ADHD.fit1, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml__ADHD_rater.txt")
sink()

sink('omega_bi_adhd_model.txt')
semTools::reliability(ADHD.fit1, return.total = TRUE, dropSingle = TRUE,
                      omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit1, "cov.lv")
sink()

sink('bi_adhd_model_PAPAJUSTGEN.txt')
ADHD.model2 <- '
ADHD =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + PAPA_p4nadhd + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity
Mother =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_mother_adhd_score.72m + SDQ72_mother_hyperactivity
Father =~ SDQ60_father_hyperactivity + conners_father_adhd_score.60m + conners_father_adhd_score.72m + SDQ72_father_hyperactivity 
Teacher =~ conners_teacher_adhd_score.72m + SDQ72_teacher_hyperactivity'

ADHD.fit2 <- lavaan::cfa(ADHD.model2, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)

o <- lavInspect(ADHD.fit2, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_ADHD_rater_PAPAJUSTGEN.txt")
s <- summary(ADHD.fit2, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml__ADHD_rater_PAPAJUSTGEN.txt")
sink()

sink('omega_bi_adhd_model_PAPAJUSTGEN.txt')
semTools::reliability(ADHD.fit2, return.total = TRUE, dropSingle = TRUE,
                      omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit2, "cov.lv")
sink()

sink('bi_adhd_model_weigh.txt')
ADHD.model3 <- '
ADHD =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + PAPA_p4nadhd + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity
Mother =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_mother_adhd_score.72m + PAPA_p4nadhd + SDQ72_mother_hyperactivity
Father =~ SDQ60_father_hyperactivity + conners_father_adhd_score.60m + conners_father_adhd_score.72m + SDQ72_father_hyperactivity 
Teacher =~ a*conners_teacher_adhd_score.72m + a*SDQ72_teacher_hyperactivity'

ADHD.fit3 <- lavaan::cfa(ADHD.model3, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)

o <- lavInspect(ADHD.fit3, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_ADHD_rater_weigh.txt")
s <- summary(ADHD.fit3, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml__ADHD_rater_weigh.txt")
sink()

sink('omega_adhd_model_weigh.txt')
semTools::reliability(ADHD.fit3, return.total = TRUE, dropSingle = TRUE,
                      omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit3, "cov.lv")
sink()

sink('bi_adhd_model_PAPAJUSTGEN_weigh.txt')
ADHD.model4 <- '
ADHD =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + PAPA_p4nadhd + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity
Mother =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_mother_adhd_score.72m + SDQ72_mother_hyperactivity
Father =~ SDQ60_father_hyperactivity + conners_father_adhd_score.60m + conners_father_adhd_score.72m + SDQ72_father_hyperactivity 
Teacher =~ a*conners_teacher_adhd_score.72m + a*SDQ72_teacher_hyperactivity'

ADHD.fit4 <- lavaan::cfa(ADHD.model4, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)

o <- lavInspect(ADHD.fit4, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_ADHD_rater_PAPAJUSTGEN_weigh.txt")
s <- summary(ADHD.fit4, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml__ADHD_rater_PAPAJUSTGEN_weigh.txt")
sink()

sink('omega_adhd_model_PAPAJUSTGEN_weigh.txt')
semTools::reliability(ADHD.fit4, return.total = TRUE, dropSingle = TRUE,
                      omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit4, "cov.lv")
sink()

sink('bi_adhd_model_woPAPA_weigh.txt')
ADHD.model5 <- '
ADHD =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity
Mother =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_mother_adhd_score.72m + SDQ72_mother_hyperactivity
Father =~ SDQ60_father_hyperactivity + conners_father_adhd_score.60m + conners_father_adhd_score.72m + SDQ72_father_hyperactivity 
Teacher =~ a*conners_teacher_adhd_score.72m + a*SDQ72_teacher_hyperactivity'

ADHD.fit5 <- lavaan::cfa(ADHD.model4, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)

o <- lavInspect(ADHD.fit5, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_ADHD_rater_woPAPA_weigh.txt")
s <- summary(ADHD.fit5, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml__ADHD_rater_woPAPA_weigh.txt")
sink()

sink('omega_adhd_model_woPAPA_weigh.txt')
semTools::reliability(ADHD.fit5, return.total = TRUE, dropSingle = TRUE,
                      omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit5, "cov.lv")
sink()

ADHD_scores_bi <- lavaan::lavPredict(ADHD.fit1)
ADHD_scores_PAPAGEN <- lavaan::lavPredict(ADHD.fit2)
ADHD_scores_bi_w <- lavaan::lavPredict(ADHD.fit3)
ADHD_scores_PAPAGEN_w <- lavaan::lavPredict(ADHD.fit4)
ADHD_scores_woPAPA <- lavaan::lavPredict(ADHD.fit5)

#NEW2 <- cbind(NEW, ADHDfactor.data1)

NEW2 <- cbind(ADHDfactor.data1, ADHD_scores_bi)

NEW3 <- cbind(NEW2, ADHD_scores_PAPAGEN)

NEW4 <- cbind(NEW3, ADHD_scores_bi_w)

NEW5 <- cbind(NEW4, ADHD_scores_PAPAGEN_w)

NEW6 <- cbind(NEW5, ADHD_scores_woPAPA)

NEW = read.csv("FEB2021.csv")

NEW7 = read.csv("composite_test.csv")

NEW8 <- cbind(NEW6, NEW)

NEW9 <- cbind(NEW8, NEW7)

readr::write_csv(NEW9, file = "/Users/Marie-Elyse/Downloads/equal_weigh_test.csv")

NEW = read.csv("equal_weigh_test.csv")

pdf("equal_weigh_test.pdf", width=5, height=5, compress=FALSE)

# "PAPA_p4nadhd", "ADHD",	"Mother",	"Father",	"Teacher", "ADHD_factor", "ADHD_factor_woPAPA",
# "ADHD_composite2", "ADHD_composite_imp"

mydata.cor = stats::cor(NEW[, c("PAPA_p4nadhd", "ADHD_factor", "ADHD",	"Mother",	"Father",	"Teacher","ADHD.1",	"Mother.1",	"Father.1",	"Teacher.1","ADHD.2",	"Mother.2",	"Father.2",	"Teacher.2","ADHD.3",	"Mother.3",	"Father.3",	"Teacher.3","ADHD.4",	"Mother.4",	"Father.4",	"Teacher.4","ADHD_composite2", "ADHD_composite_imp")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("PAPA_p4nadhd", "ADHD_factor", "ADHD",	"Mother",	"Father",	"Teacher","ADHD.1",	"Mother.1",	"Father.1",	"Teacher.1","ADHD.2",	"Mother.2",	"Father.2",	"Teacher.2","ADHD.3",	"Mother.3",	"Father.3",	"Teacher.3","ADHD.4",	"Mother.4",	"Father.4",	"Teacher.4","ADHD_composite2", "ADHD_composite_imp")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

colnames(mydata.cor) <-c("PAPA", "ADHD unifactor", "ADHD",	"Mother",	"Father",	"Teacher","ADHD.1",	"Mother.1",	"Father.1",	"Teacher.1","ADHD.2",	"Mother.2",	"Father.2",	"Teacher.2","ADHD.3",	"Mother.3",	"Father.3",	"Teacher.3","ADHD.4",	"Mother.4",	"Father.4",	"Teacher.4","ADHD_composite2", "ADHD_composite_imp")

rownames(mydata.cor) <-c("")


corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.05, .01, .001), pch.cex = .4, tl.cex = .4,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

