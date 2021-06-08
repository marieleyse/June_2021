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

sink('uni_adhd_model.txt')
ADHD.model <- '
ADHD_factor =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + PAPA_p4nadhd + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity'

ADHD.fit <- lavaan::cfa(ADHD.model, data = ADHDfactor.data1, estimator = "MLR", missing = "FIML", std.ov = TRUE, std.lv = TRUE, orthogonal = T, verbose = T)

o <- lavInspect(ADHD.fit, what ='fit')
capture.output(o, file = "fit_indices_ADHD_model_fiml_onefactor.txt")
s <- summary(ADHD.fit, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml_onefactor.txt")
sink()

sink('omega_uni_adhd_model.txt')
semTools::reliability(ADHD.fit, return.total = TRUE, dropSingle = TRUE,
            omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit, "cov.lv")
sink()

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
s <- summary(ADHD.fit2, standardized = T, fit.measures = T)
capture.output(s, file = "model_output_ADHD_model_fiml__ADHD_rater_PAPAJUSTGEN_weigh.txt")
sink()

sink('omega_adhd_model_PAPAJUSTGEN_weigh.txt')
semTools::reliability(ADHD.fit4, return.total = TRUE, dropSingle = TRUE,
                      omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit4, "cov.lv")
sink()

ADHD_scores <- lavaan::lavPredict(ADHD.fit)
ADHD_scores_bi <- lavaan::lavPredict(ADHD.fit1)
ADHD_scores_PAPAGEN <- lavaan::lavPredict(ADHD.fit2)
ADHD_scores_bi_w <- lavaan::lavPredict(ADHD.fit3)
ADHD_scores_PAPAGEN_w <- lavaan::lavPredict(ADHD.fit4)

#NEW2 <- cbind(NEW, ADHDfactor.data1)

NEW2 <- cbind(ADHDfactor.data1, ADHD_scores)

NEW3 <- cbind(NEW2, ADHD_scores_bi)

NEW4 <- cbind(NEW3, ADHD_scores_PAPAGEN)

NEW5 <- cbind(NEW4, ADHD_scores_bi_w)

NEW6 <- cbind(NEW5, ADHD_scores_PAPAGEN_w)

readr::write_csv(NEW6, file = "/Users/Marie-Elyse/Downloads/equal_weigh_test.csv")

NEW = read.csv("equal_weigh_test.csv")

pdf("equal_weigh_test.pdf", width=5, height=5, compress=FALSE)

mydata.cor = stats::cor(NEW[, c("ADHD_factor", "ADHD",	"Mother",	"Father",	"Teacher","ADHD.1",	"Mother.1",	"Father.1",	"Teacher.1","ADHD.2",	"Mother.2",	"Father.2",	"Teacher.2","ADHD.3",	"Mother.3",	"Father.3",	"Teacher.3")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("ADHD_factor", "ADHD",	"Mother",	"Father",	"Teacher","ADHD.1",	"Mother.1",	"Father.1",	"Teacher.1","ADHD.2",	"Mother.2",	"Father.2",	"Teacher.2","ADHD.3",	"Mother.3",	"Father.3",	"Teacher.3")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.05, .01, .001), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

NEW = read.csv("equal_weigh_test.csv")

NEW2 = read.csv("FEB2021.csv")

NEW3 = cbind(NEW, NEW2)

readr::write_csv(NEW3, file = "/Users/Marie-Elyse/Downloads/equal_weigh_test2.csv")

NEW = read.csv("equal_weigh_test2.csv")

pdf("equal_weigh_test2.pdf", width=5, height=5, compress=FALSE)

mydata.cor = stats::cor(NEW[, c("ADHD_factor", "ADHD",	"Mother",	"Father",	"Teacher","ADHD.1",	"Mother.1",	"Father.1",	"Teacher.1","ADHD.2",	"Mother.2",	"Father.2",	"Teacher.2","ADHD.3",	"Mother.3",	"Father.3",	"Teacher.3","ADHD.4",	"Mother.4",	"Father.4",	"Teacher.4")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("ADHD_factor", "ADHD",	"Mother",	"Father",	"Teacher","ADHD.1",	"Mother.1",	"Father.1",	"Teacher.1","ADHD.2",	"Mother.2",	"Father.2",	"Teacher.2","ADHD.3",	"Mother.3",	"Father.3",	"Teacher.3","ADHD.4",	"Mother.4",	"Father.4",	"Teacher.4")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.05, .01, .001), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()


