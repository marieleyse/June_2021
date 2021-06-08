#LEGIT ADHD ~ G(PRS)xE(M+A) + covariates

setwd("/Users/Marie-Elyse/Downloads")
NEW = read.csv("FEB2021.csv")

#####
#####
#####
#####merge M and A

NEW2 = read.csv("A_M.csv")

#NEW3 = cbind(NEW, NEW2)

NEW3 =  merge(NEW, NEW2, by="PSCID")

readr::write_csv(NEW3, file = "/Users/Marie-Elyse/Downloads/APRIL2021.csv")

NEW = read.csv("APRIL2021.csv")

#####
#####
#####
#####make sure I use the appropriate ADHD factor (same structure as per ALSPAC)

#Subset data according to who has at least one item present based on the updated (e.g. parecellated) number of items.
subsample_one_item_present1 <- NEW[which((rowMeans(is.na(NEW[c("CBCL48_sc6raw.x", "SDQ60_mother_hyperactivity","SDQ60_father_hyperactivity", "CBCL60_sc6raw.x","conners_mother_adhd_score.60m", "conners_father_adhd_score.60m", "conners_mother_adhd_score.72m", "conners_father_adhd_score.72m", "conners_teacher_adhd_score.72m", "PAPA_p4nadhd", "Dominic72_ADHD","SDQ72_mother_hyperactivity", "SDQ72_father_hyperactivity", "SDQ72_teacher_hyperactivity")]))) <= 13/14),]
dim(subsample_one_item_present1)

ADHDfactor.data1 <- subsample_one_item_present1[c("CBCL48_sc6raw.x", "SDQ60_mother_hyperactivity","SDQ60_father_hyperactivity", "CBCL60_sc6raw.x","conners_mother_adhd_score.60m", "conners_father_adhd_score.60m", "conners_mother_adhd_score.72m", "conners_father_adhd_score.72m", "conners_teacher_adhd_score.72m", "PAPA_p4nadhd", "Dominic72_ADHD","SDQ72_mother_hyperactivity", "SDQ72_father_hyperactivity", "SDQ72_teacher_hyperactivity")]
str(ADHDfactor.data1)

sink('bi_adhd_model_weigh.txt')
ADHD.model3 <- '
ADHD.w =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + SDQ60_father_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_father_adhd_score.60m + conners_mother_adhd_score.72m + conners_father_adhd_score.72m + conners_teacher_adhd_score.72m + PAPA_p4nadhd + Dominic72_ADHD + SDQ72_mother_hyperactivity + SDQ72_father_hyperactivity + SDQ72_teacher_hyperactivity
Mother.w =~ CBCL48_sc6raw.x + SDQ60_mother_hyperactivity + CBCL60_sc6raw.x + conners_mother_adhd_score.60m + conners_mother_adhd_score.72m + PAPA_p4nadhd + SDQ72_mother_hyperactivity
Father.w =~ SDQ60_father_hyperactivity + conners_father_adhd_score.60m + conners_father_adhd_score.72m + SDQ72_father_hyperactivity 
Teacher.w =~ a*conners_teacher_adhd_score.72m + a*SDQ72_teacher_hyperactivity'

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

ADHD_scores_bi_w <- lavaan::lavPredict(ADHD.fit3)

NEW2 <- cbind(NEW, ADHD_scores_bi_w)

readr::write_csv(NEW2, file = "/Users/Marie-Elyse/Downloads/APRIL2021_BIFACTORw.csv")

NEW = read.csv("APRIL2021_BIFACTORw.csv")

plot_LEGIT = function(x, cov_values = NULL, gene_quant = c(.025,.50,.975), env_quant = c(.025,.50,.975), outcome_quant = c(.025,.50,.975), cols = c("#3288BD", "#CAB176", "#D53E4F"), ylab="Outcome", xlab="Environment", legtitle="Genetic score", leglab=NULL, xlim= NULL, ylim= NULL, x_at = NULL, y_at = NULL, cex.axis = 1.9, cex.lab=2, cex.main=2.2, cex.leg=2.2, legend="topleft", ...){
  
  # Better names (need to use x for S3 class consistency)
  object = x
  formula = object$formula
  
  # Checks
  
  # Quantile range
  gene_range = as.numeric(quantile(object$fit_main$data$G, gene_quant))
  env_range = as.numeric(quantile(object$fit_main$data$E, env_quant))
  formula_outcome = formula.tools::get.vars(formula)[1]
  outcome_range = as.numeric(quantile(object$fit_main$data[formula_outcome][,], outcome_quant))
  # Fix this attribute to prevent problems if trying to predict fit_main directly later on
  attr(object$fit_main$terms,"dataClasses")[attr(object$fit_main$terms,"dataClasses")=="nmatrix.1"] = "numeric"
  
  # Defaults
  if (is.null(ylim)){
    if (is.null(y_at)) ylim = c(outcome_range[1], outcome_range[length(outcome_range)])
    else ylim = c(y_at[1], y_at[length(y_at)])
  }
  else if (is.null(y_at)) y_at = seq(ylim[1], ylim[2], length.out=3) 
  if (is.null(xlim)){
    if (is.null(x_at)) xlim = c(env_range[1], env_range[length(env_range)])
    else xlim = c(x_at[1], x_at[length(x_at)])
  }
  else if (is.null(x_at)) x_at = seq(xlim[1], xlim[2], length.out=3) 
  
  # Plot
  op <- par(mfrow=c(1,1), mar=c(5.1, 5.1, 5.1, 4.1))
  graphics::plot(x=c(), y=c() ,ylab = ylab, xlab = xlab, ylim = ylim, xlim=xlim, cex.lab=cex.lab, cex.main=cex.main, pch=20, bty="n", xaxt="n", yaxt="n")
  if (is.null(y_at)) graphics::axis(2, at=outcome_range, labels=paste0(outcome_quant*100,"%"), cex.axis = cex.axis)
  else graphics::axis(2, at=y_at, cex.axis = cex.axis)
  if (is.null(x_at)) graphics::axis(1, at=env_range, labels=paste0(env_quant*100,"%"), cex.axis = cex.axis)
  else graphics::axis(1, at=x_at, cex.axis = cex.axis)
  E = seq(xlim[1],xlim[2], length.out=101)
  if (!is.null(object$crossover)) E_ = E - object$crossover
  else E_ = E
  # covariates
  covariates = formula.tools::get.vars(formula)[-1]
  covariates = covariates[covariates != "G" & covariates != "E"]
  if (!is.null(cov_values) && length(cov_values)!=length(covariates)) stop("cov_values doesn't have the correct number of covariates")
  
  # Predictions and plot
  for (j in 1:length(gene_range)){
    # making data for predictions
    G = gene_range[j]
    newdata_base = data.frame(E=E_, G=G)
    newdata = newdata_base
    for (i in 1:length(covariates)){
      if (identical(covariates, character(0))) newdata = newdata
      else if (is.null(cov_values)){
        newdata = cbind(newdata, rep(apply(object$fit_main$data[covariates[i]], 2, mean),101))
        colnames(newdata)[NCOL(newdata)] = covariates[i]
      }
      else{
        newdata = cbind(newdata, rep(cov_values[i], 101))
        colnames(newdata)[NCOL(newdata)] = names(cov_values)[i]
      }
    }
    # Prediction lines
    preds <- predict(object$fit_main, newdata = newdata, se.fit = TRUE, type = "link")
    ilink <- family(object$fit_main)$linkinv
    graphics::polygon(c(E,rev(E)),c(ilink(preds$fit-2*preds$se.fit),rev(ilink(preds$fit+2*preds$se.fit))),col=grDevices::adjustcolor(cols[j], alpha.f = 0.2),border = NA)
    graphics::lines(E,ilink(preds$fit), col=cols[j], lwd=2)
  }
  if (is.null(leglab)) leglab = paste0(gene_quant*100,"%")
  #legend(legend, legend=leglab,col = cols, lty=1, lwd=3, xpd = TRUE, cex = cex.leg, title=legtitle)
}

plot_LEGIT_3way = function(x, gender_range=c(0,1), gender_name="gender_male", gender_title=c("Female","Male"), cov_values = NULL, gene_quant = c(.1,.9), env_quant = c(.025,.50,.975), outcome_quant = c(.025,.50,.975), cols = c("blue", "deeppink"), ylab="ADHD", xlab="Maternal depression", legtitle="Genetic score", leglab=NULL, xlim= NULL, ylim= NULL, x_at = NULL, y_at = NULL, cex.axis = 0.8, cex.lab=0.8, cex.main=0.8, cex.leg=0.8, legend="topleft", ...){
  
  # Better names (need to use x for S3 class consistency)
  object = x
  formula = object$formula
  
  # Checks
  
  # Quantile range
  gene_range = as.numeric(quantile(object$fit_main$data$G, gene_quant))
  env_range = as.numeric(quantile(object$fit_main$data$E, env_quant))
  formula_outcome = formula.tools::get.vars(formula)[1]
  outcome_range = as.numeric(quantile(object$fit_main$data[formula_outcome][,], outcome_quant))
  # Fix this attribute to prevent problems if trying to predict fit_main directly later on
  attr(object$fit_main$terms,"dataClasses")[attr(object$fit_main$terms,"dataClasses")=="nmatrix.1"] = "numeric"
  
  # Defaults
  if (is.null(ylim)){
    if (is.null(y_at)) ylim = c(outcome_range[1], outcome_range[length(outcome_range)])
    else ylim = c(y_at[1], y_at[length(y_at)])
  }
  else if (is.null(y_at)) y_at = seq(ylim[1], ylim[2], length.out=3) 
  if (is.null(xlim)){
    if (is.null(x_at)) xlim = c(env_range[1], env_range[length(env_range)])
    else xlim = c(x_at[1], x_at[length(x_at)])
  }
  else if (is.null(x_at)) x_at = seq(xlim[1], xlim[2], length.out=3) 
  
  op <- par(mfrow=c(1,2), mar=c(5.1, 5.1, 5.1, 4.1))
  E = seq(xlim[1],xlim[2], length.out=101)
  if (!is.null(object$crossover)) E_ = E - object$crossover
  else E_ = E
  
  # covariates
  covariates = formula.tools::get.vars(formula)[-1]
  covariates = covariates[covariates != "G" & covariates != "E"]
  if (!is.null(cov_values) && length(cov_values)!=length(covariates)) stop("cov_values doesn't have the correct number of covariates")
  
  # Predictions and plot
  for (ff in 1:length(gender_range)){
    Gender = gender_range[ff]
    
    # Plot
    graphics::plot(x=c(), y=c() ,ylab = ylab, xlab = xlab, ylim = ylim, xlim=xlim, cex.lab=cex.lab, cex.main=cex.main, pch=10, bty="n", xaxt="n", yaxt="n", main=gender_title[ff], font.main = 1)
    if (is.null(y_at)) graphics::axis(2, at=outcome_range, labels=paste0(outcome_quant*100,"%"), cex.axis = cex.axis)
    else graphics::axis(2, at=y_at, cex.axis = cex.axis)
    if (is.null(x_at)) graphics::axis(1, at=env_range, labels=paste0(env_quant*100,"%"), cex.axis = cex.axis)
    else graphics::axis(1, at=x_at, cex.axis = cex.axis)
    
    for (j in 1:length(gene_range)){
      # making data for predictions
      G = gene_range[j]
      newdata = data.frame(E=E_, G=G, gender=Gender)
      colnames(newdata)[3] = gender_name
      for (i in 1:length(covariates)){
        if (identical(covariates, character(0))) newdata = newdata
        else if (is.null(cov_values)){
          newdata = cbind(newdata, rep(apply(object$fit_main$data[covariates[i]], 2, mean),101))
          colnames(newdata)[NCOL(newdata)] = covariates[i]
        }
        else{
          newdata = cbind(newdata, rep(cov_values[i], 101))
          colnames(newdata)[NCOL(newdata)] = names(cov_values)[i]
        }
      }
      # Prediction lines
      preds <- predict(object$fit_main, newdata = newdata, se.fit = TRUE, type = "link")
      ilink <- family(object$fit_main)$linkinv
      graphics::polygon(c(E,rev(E)),c(ilink(preds$fit-2*preds$se.fit),rev(ilink(preds$fit+2*preds$se.fit))),col=grDevices::adjustcolor(cols[j], alpha.f = 0.2),border = NA)
      graphics::lines(E,ilink(preds$fit), col=cols[j], lwd=2)
    }
    if (is.null(leglab)) leglab = paste0(gene_quant*100,"%")
    legend(legend, legend=leglab,col = cols, lty=1, lwd=3, xpd = TRUE, cex = cex.leg, title=legtitle)
  }
}


#####
#####
#####
#####
#####
#####
#####
#####
##### add M and A 

pdf("correlation_prs_adhd_m_a_covariates.pdf", width=5, height=5, compress=FALSE)
mydata.cor = stats::cor(NEW[, c("ADHD.w",	"Mother.w",	"Father.w",	"Teacher.w",
                                "PRS_0_001_adhd_child", "PC1", 
                                "PC2", "PC3", 
                                "gender_male", 
                                "auc_post_cesd", "AFactor_A",
                                "MFactor_GEN", "mom_age_birth",
                                "above_college","Hamilton")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("ADHD.w",	"Mother.w",	"Father.w",	"Teacher.w",
                                    "PRS_0_001_adhd_child", "PC1", 
                                    "PC2", "PC3", 
                                    "gender_male",	
                                    "auc_post_cesd", "AFactor_A",
                                    "MFactor_GEN", "mom_age_birth",
                                    "above_college","Hamilton")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

colnames(mydata.cor) <-c("ADHD",	"Mother",	"Father",	"Teacher",
                         "PRS_0_001_adhd_child", "PC1", 
                         "PC2", "PC3", 
                         "Child sex (male=1)",	 
                         "Postnatal depression", "A Factor",
                         "M Factor","Maternal age at birth",
                         "Education","Site")

rownames(mydata.cor) <-c("ADHD",	"Mother",	"Father",	"Teacher",
                         "PRS_0_001_adhd_child", "PC1", 
                         "PC2", "PC3", 
                         "Child sex (male=1)",	
                         "Postnatal depression", "A Factor",
                         "M Factor", "Maternal age at birth",
                         "Education","Site")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.001, .01), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

pdf("correlation_prs_adhd_m_a_covariates_p05.pdf", width=5, height=5, compress=FALSE)
corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.05), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

#########
#########z score

#####
#####
#####
#####
#####
#####
#####
#####

NEW$PRE_BY_POST = (NEW$auc_post_cesd*NEW$Pren_CESD)

NEW$std_ADHD = c(scale(NEW$ADHD))
NEW$std_ADHD.w = c(scale(NEW$ADHD.w))
NEW$std_PAPA_p4nadhd = c(scale(NEW$PAPA_p4nadhd))
NEW$std_conners_mother_hyperactivity_score.72m = c(scale(NEW$conners_mother_hyperactivity_score.72m))
NEW$std_PRS_0_001_adhd_child = c(scale(NEW$PRS_0_001_adhd_child))
NEW$std_PRS_0_3_ADHD_child = c(scale(NEW$PRS_0_3_ADHD_child))
NEW$std_auc_post_cesd = c(scale(NEW$auc_post_cesd))
NEW$std_Pren_CESD = c(scale(NEW$Pren_CESD))
NEW$std_PC1 = c(scale(NEW$PC1))
NEW$std_PC2 = c(scale(NEW$PC2))
NEW$std_PC3 = c(scale(NEW$PC3))
NEW$std_mom_age_birth = c(scale(NEW$mom_age_birth))
NEW$std_PRE_BY_POST = c(scale(NEW$PRE_BY_POST))
NEW$std_A = c(scale(NEW$AFactor_A))
NEW$std_M = c(scale(NEW$MFactor_GEN))

#####
#####
#####
#####
#####
#####
#####
#####
##### add M and A and remove site std_A std_M
##### option in LEGIT() and even glm() is family=quasibinomial(link = "logit")
##### It constrain to [0,1] so make sure that the outcome range is transformed to that range.
##### g+e+a+g*a+g*e+a*e

######
######
######2-way

G = NEW[,c("std_PRS_0_3_ADHD_child"), drop=FALSE]
E = NEW[,c("std_A", "std_M"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male + std_auc_post_cesd, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + postnatal depression + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()


comp = stats::complete.cases(NEW[c("std_ADHD.w", "std_A", "std_M",
                                   "gender_male",	"std_mom_age_birth", "std_auc_post_cesd",
                                   "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD.w - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD.w - mean(NEW_nomiss$std_ADHD.w))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male + std_auc_post_cesd, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_PRS_0_3_ADHD_child"), drop=FALSE]
E = NEW[,c("std_A", "std_M"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + std_auc_post_cesd, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + postnatal depression + std PC1-3'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_prs_a_m_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + postnatal depression + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD.w",
                                   "gender_male",	"std_A", "std_M", "std_mom_age_birth", "std_auc_post_cesd",
                                   "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD.w - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD.w - mean(NEW_nomiss$std_ADHD.w))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + std_auc_post_cesd, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

#######wo postnatal depression

G = NEW[,c("std_PRS_0_3_ADHD_child"), drop=FALSE]
E = NEW[,c("std_A", "std_M"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov_wopostnataldepression.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov_wopostnataldepression.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education  + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()


comp = stats::complete.cases(NEW[c("std_ADHD.w", "std_A", "std_M",
                                   "gender_male",	"std_mom_age_birth", 
                                   "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD.w - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD.w - mean(NEW_nomiss$std_ADHD.w))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_PRS_0_3_ADHD_child"), drop=FALSE]
E = NEW[,c("std_A", "std_M"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov_sex_wopostnataldepression.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov_sex_wopostnataldepression.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + std PC1-3'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_prs_a_m_cov_sex_sexfacet_wopostnataldepression.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD.w",
                                   "gender_male",	"std_A", "std_M", "std_mom_age_birth", 
                                   "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD.w - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD.w - mean(NEW_nomiss$std_ADHD.w))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

############
#############

G = NEW[,c("std_PRS_0_3_ADHD_child"), drop=FALSE]
E = NEW[,c("std_A", "std_M", "std_auc_post_cesd"),  drop=FALSE]

pdf("LEGIT_prs_a_m_post_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_post_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()


comp = stats::complete.cases(NEW[c("std_ADHD.w", "std_A", "std_M",
                                   "gender_male",	"std_mom_age_birth", "std_auc_post_cesd",
                                   "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD.w - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD.w - mean(NEW_nomiss$std_ADHD.w))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_PRS_0_3_ADHD_child"), drop=FALSE]
E = NEW[,c("std_A", "std_M", "std_auc_post_cesd"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov_post_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov_post_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + std PC1-3'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_prs_a_m_cov_post_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD.w",
                                   "gender_male",	"std_A", "std_M", "std_mom_age_birth", "std_auc_post_cesd",
                                   "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD.w - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD.w - mean(NEW_nomiss$std_ADHD.w))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD.w ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()



# ########adhd test
# 
# G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
# E = NEW[,c("std_A", "std_M"),  drop=FALSE]
# 
# pdf("LEGIT_prs_a_m_cov.pdf", width=5, height=5, compress=FALSE)
# sink(paste0('LEGIT_prs_a_m_cov.txt'))
# fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male + std_auc_post_cesd, genes=G,env=E)
# fit_genes = fit_legit_list$genes
# fit_env = fit_legit_list$env
# print(summary(fit_legit_list))
# print(fit_legit_list)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + postnatal depression + std PC1-3"), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# save.image(file='AUG_2020.RData')
# dev.off()
# 
# 
# comp = stats::complete.cases(NEW[c("std_ADHD", "std_A", "std_M",
#                                    "gender_male",	"std_mom_age_birth", "std_auc_post_cesd",
#                                    "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
# G_nomiss = G[comp, , drop = FALSE]
# E_nomiss = E[comp, , drop = FALSE]
# NEW_nomiss = NEW[comp, , drop = FALSE]
# 
# ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
# sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD.w))^2)
# R2 = 1 - ssres/sstotal
# 
# # # Cross-validation 5 times with 5 Folds
# # cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# # mean(cv_5folds$R2_cv)
# # print(mean(cv_5folds$R2_cv))
# # print(summary(cv_5folds))
# # print(cv_5folds)
# 
# # Leave-one-out cross-validation (Note: very slow)
# cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male + std_auc_post_cesd, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
# mean(cv_loo$R2_cv)
# print(mean(cv_loo$R2_cv))
# # print(summary(cv_loo))
# # print(cv_loo)
# 
# sink()
# 
# ######
# ######
# ######3-way
# 
# G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
# E = NEW[,c("std_A", "std_M"),  drop=FALSE]
# 
# pdf("LEGIT_prs_a_m_cov_sex.pdf", width=5, height=5, compress=FALSE)
# sink(paste0('LEGIT_prs_a_m_cov_sex.txt'))
# fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + std_auc_post_cesd, genes=G,env=E)
# fit_genes = fit_legit_list$genes
# fit_env = fit_legit_list$env
# print(summary(fit_legit_list))
# print(fit_legit_list)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + postnatal depression + std PC1-3'), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# dev.off()
# 
# pdf("LEGIT_prs_a_m_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
# plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
# mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + postnatal depression + std PC1-3",                   # Add main title
#       side = 3,
#       line = - 2,
#       outer = TRUE)
# dev.off()
# 
# comp = stats::complete.cases(NEW[c("std_ADHD",
#                                    "gender_male",	"std_A", "std_M", "std_mom_age_birth", "std_auc_post_cesd",
#                                    "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
# G_nomiss = G[comp, , drop = FALSE]
# E_nomiss = E[comp, , drop = FALSE]
# NEW_nomiss = NEW[comp, , drop = FALSE]
# 
# ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
# sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD.w))^2)
# R2 = 1 - ssres/sstotal
# 
# # # Cross-validation 5 times with 5 Folds
# # cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E, cv_iter=5, cv_folds=5)
# # mean(cv_5folds$R2_cv)
# # print(mean(cv_5folds$R2_cv))
# # print(summary(cv_5folds))
# # print(cv_5folds)
# 
# # Leave-one-out cross-validation (Note: very slow)
# cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + std_auc_post_cesd, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
# mean(cv_loo$R2_cv)
# print(mean(cv_loo$R2_cv))
# # print(summary(cv_loo))
# # print(cv_loo)
# 
# sink()
# 
# #######wo postnatal depression
# 
# G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
# E = NEW[,c("std_A", "std_M"),  drop=FALSE]
# 
# pdf("LEGIT_prs_a_m_cov_wopostnataldepression.pdf", width=5, height=5, compress=FALSE)
# sink(paste0('LEGIT_prs_a_m_cov_wopostnataldepression.txt'))
# fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E)
# fit_genes = fit_legit_list$genes
# fit_env = fit_legit_list$env
# print(summary(fit_legit_list))
# print(fit_legit_list)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education  + std PC1-3"), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# dev.off()
# 
# 
# comp = stats::complete.cases(NEW[c("std_ADHD", "std_A", "std_M",
#                                    "gender_male",	"std_mom_age_birth", 
#                                    "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
# G_nomiss = G[comp, , drop = FALSE]
# E_nomiss = E[comp, , drop = FALSE]
# NEW_nomiss = NEW[comp, , drop = FALSE]
# 
# ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
# sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
# R2 = 1 - ssres/sstotal
# 
# # # Cross-validation 5 times with 5 Folds
# # cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# # mean(cv_5folds$R2_cv)
# # print(mean(cv_5folds$R2_cv))
# # print(summary(cv_5folds))
# # print(cv_5folds)
# 
# # Leave-one-out cross-validation (Note: very slow)
# cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
# mean(cv_loo$R2_cv)
# print(mean(cv_loo$R2_cv))
# # print(summary(cv_loo))
# # print(cv_loo)
# 
# sink()
# 
# ######
# ######
# ######3-way
# 
# G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
# E = NEW[,c("std_A", "std_M"),  drop=FALSE]
# 
# pdf("LEGIT_prs_a_m_cov_sex_wopostnataldepression.pdf", width=5, height=5, compress=FALSE)
# sink(paste0('LEGIT_prs_a_m_cov_sex_wopostnataldepression.txt'))
# fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E)
# fit_genes = fit_legit_list$genes
# fit_env = fit_legit_list$env
# print(summary(fit_legit_list))
# print(fit_legit_list)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + std PC1-3'), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# dev.off()
# 
# pdf("LEGIT_prs_a_m_cov_sex_sexfacet_wopostnataldepression.pdf", width=10, height=5, compress=FALSE)
# plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
# mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + std PC1-3",                   # Add main title
#       side = 3,
#       line = - 2,
#       outer = TRUE)
# dev.off()
# 
# comp = stats::complete.cases(NEW[c("std_ADHD",
#                                    "gender_male",	"std_A", "std_M", "std_mom_age_birth", 
#                                    "above_college","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
# G_nomiss = G[comp, , drop = FALSE]
# E_nomiss = E[comp, , drop = FALSE]
# NEW_nomiss = NEW[comp, , drop = FALSE]
# 
# ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
# sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
# R2 = 1 - ssres/sstotal
# 
# # # Cross-validation 5 times with 5 Folds
# # cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E, cv_iter=5, cv_folds=5)
# # mean(cv_5folds$R2_cv)
# # print(mean(cv_5folds$R2_cv))
# # print(summary(cv_5folds))
# # print(cv_5folds)
# 
# # Leave-one-out cross-validation (Note: very slow)
# cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
# mean(cv_loo$R2_cv)
# print(mean(cv_loo$R2_cv))
# # print(summary(cv_loo))
# # print(cv_loo)
# 
# sink()
# 
sink('characteristic_description.txt')
print(summary(NEW[, c("ADHD",
                      "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                      "auc_post_cesd", "mom_age_birth",
                      "above_college","Hamilton","AFactor_A",	"MFactor_GEN",
                       "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                      "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                      "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                      "birth_size_percent2_x","AGE_BY_SITE_corr",
                      "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                      "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")]))

describeVEC <-Hmisc::describe(comp[, c("ADHD",
                                       "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                       "auc_post_cesd", "mom_age_birth",
                                       "above_college","Hamilton","AFactor_A",	"MFactor_GEN", "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                                       "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                       "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                       "birth_size_percent2_x","AGE_BY_SITE_corr",
                                       "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                       "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], transpose = TRUE)
print(describeVEC)

library(pastecs)
stat.describeVEC <-pastecs::stat.desc(NEW[, c("ADHD",
                                              "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                              "auc_post_cesd", "mom_age_birth",
                                              "above_college","Hamilton","AFactor_A",	"MFactor_GEN", "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                                              "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                              "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                              "birth_size_percent2_x","AGE_BY_SITE_corr",
                                              "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                              "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")])
print(stat.describeVEC)

###best
psych.describeVEC <-psych::describe(NEW[, c("ADHD",
                                            "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                            "auc_post_cesd", "mom_age_birth",
                                            "above_college","Hamilton","AFactor_A",	"MFactor_GEN", "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                                            "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                            "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                            "birth_size_percent2_x","AGE_BY_SITE_corr",
                                            "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                            "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")])
print(psych.describeVEC)
sink()

sink('characteristic_description_by_sex.txt')
psych.describeVEC.sex <-psych::describeBy(NEW[, c("ADHD",
                                                  "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                  "auc_post_cesd", "mom_age_birth",
                                                  "above_college","Hamilton","AFactor_A",	"MFactor_GEN", "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                                                  "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                  "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                  "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                  "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                  "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], NEW$gender_male)
print(psych.describeVEC.sex)

sink()

sink('characteristic_description_by_site.txt')
psych.describeVEC.site <-psych::describeBy(NEW[, c("ADHD",
                                                   "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                   "auc_post_cesd", "mom_age_birth",
                                                   "above_college","Hamilton", "AFactor_A",	"MFactor_GEN", "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                                                   "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                   "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                   "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                   "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                   "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], NEW$Hamilton)
print(psych.describeVEC.site)
sink()

sink('characteristic_description_by_adhd_dx.txt')
psych.describeVEC.ADHD.DX <-psych::describeBy(NEW[, c("ADHD",
                                                      "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                      "auc_post_cesd", "mom_age_birth",
                                                      "above_college","Hamilton","AFactor_A",	"MFactor_GEN", "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                                                      "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                      "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                      "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                      "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                      "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], NEW$PAPA_p4_adhd)
print(psych.describeVEC.ADHD.DX)
sink()

sink('characteristic_description_by_sibling.txt')
psych.describeVEC.Sibling <-psych::describeBy(NEW[, c("ADHD",
                                                      "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                      "auc_post_cesd", "mom_age_birth",
                                                      "above_college","Hamilton","AFactor_A",	"MFactor_GEN", "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                                                      "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                      "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                      "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                      "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                      "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], NEW$Sibling)
print(psych.describeVEC.Sibling)
sink()

sink('characteristic_description_by_adhd_median.txt')
NEW$ADHD_d<-sjmisc::dicho(
  NEW$ADHD,
  dich.by = "median",
  as.num = TRUE,
  var.label = NULL,
  val.labels = NULL,
  append = TRUE,
  suffix = "_d"
)

psych.describeVEC.ADHD <-psych::describeBy(NEW[, c("ADHD",
                                                   "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                   "auc_post_cesd", "mom_age_birth",
                                                   "above_college","Hamilton", "AFactor_A",	"MFactor_GEN", "PRS_0_001_adhd_child", "PRS_0_3_ADHD_child", "PC1", "PC2", "PC3",
                                                   "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                   "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                   "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                   "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                   "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], NEW$ADHD_d)
print(psych.describeVEC.ADHD)

sink()

sink('sex_differences_kk_adhd_demographics.txt')
kruskal.test(AFactor_A ~ gender_male, data=NEW)
kruskal.test(MFactor_GEN ~ gender_male, data=NEW)
kruskal.test(PAPA_p4nadhd ~ gender_male, data=NEW)
kruskal.test(conners_mother_hyperactivity_score.72m ~ PAPA_p4_adhd, data=NEW)
kruskal.test(ADHD ~ gender_male, data=NEW)
kruskal.test(Mother ~ gender_male, data=NEW)
kruskal.test(Father ~ gender_male, data=NEW)
kruskal.test(Teacher ~ gender_male, data=NEW)
kruskal.test(mom_age_birth ~ gender_male, data=NEW)
kruskal.test(Pren_income4 ~ gender_male, data=NEW)
kruskal.test(Pren_life_events ~ gender_male, data=NEW)
kruskal.test(QuintMat_w ~ gender_male, data=NEW)
kruskal.test(AGE_BY_SITE_corr~ gender_male, data=NEW)
kruskal.test(QuintSoc_w ~ gender_male, data=NEW)
kruskal.test(birth_size_percent2_x ~ gender_male, data=NEW)
kruskal.test(birth_wt_grams ~ gender_male, data=NEW)
kruskal.test(gestation_age_wks ~ gender_male, data=NEW)
kruskal.test(Pren_CESD ~ gender_male, data=NEW)
kruskal.test(HWB6_CESD ~ gender_male, data=NEW)
kruskal.test(HWB12_CESD ~ gender_male, data=NEW)
kruskal.test(HWB24_CESD ~ gender_male, data=NEW)
kruskal.test(HWB36_CESD ~ gender_male, data=NEW)
kruskal.test(HWB48_CESD ~ gender_male, data=NEW)
kruskal.test(HWB60_CESD ~ gender_male, data=NEW)
kruskal.test(HWB72_CESD ~ gender_male, data=NEW)
kruskal.test(auc_post_cesd ~ gender_male, data=NEW)

kruskal.test(B_DRD2_cat ~ gender_male, data=NEW)
kruskal.test(B_DRD3_rs6280_cat ~ gender_male, data=NEW)
kruskal.test(B_COMT_cat.y ~ gender_male, data=NEW)
kruskal.test(B_COMT_rs165599.y ~ gender_male, data=NEW)

kruskal.test(PRS_0_3_ADHD_child ~ gender_male, data=NEW)
kruskal.test(PRS_0_001_adhd_child ~ gender_male, data=NEW)
kruskal.test(PC1 ~ gender_male, data=NEW)
kruskal.test(PC2 ~ gender_male, data=NEW)
kruskal.test(PC3 ~ gender_male, data=NEW)
kruskal.test(birth_wt_grams ~ gender_male, data=NEW)
kruskal.test(gestation_age_wks ~ gender_male, data=NEW)
kruskal.test(Alcohol_During_Pregnancy ~ PAPA_p4_adhd, data=NEW)
kruskal.test(Smoking_During_Pregnancy ~ gender_male, data=NEW)

chisq.test(NEW$above_college, NEW$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$Sibling, NEW$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$PAPA_p4_adhd, NEW$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(NEW$B_DRD4, NEW$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DAT.y, NEW$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD2_rs1799978.y, NEW$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD1_hap.y, NEW$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

sink('ADHD_differences_kk_adhd_demographics.txt')
kruskal.test(AFactor_A ~ PAPA_p4_adhd, data=NEW)
kruskal.test(MFactor_GEN ~ PAPA_p4_adhd, data=NEW)
kruskal.test(PAPA_p4nadhd ~ PAPA_p4_adhd, data=NEW)
kruskal.test(ADHD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(Mother ~ PAPA_p4_adhd, data=NEW)
kruskal.test(Father ~ PAPA_p4_adhd, data=NEW)
kruskal.test(Teacher ~ PAPA_p4_adhd, data=NEW)
kruskal.test(mom_age_birth ~ PAPA_p4_adhd, data=NEW)
kruskal.test(Pren_income4 ~ PAPA_p4_adhd, data=NEW)
kruskal.test(Pren_life_events ~ PAPA_p4_adhd, data=NEW)
kruskal.test(QuintMat_w ~ PAPA_p4_adhd, data=NEW)
kruskal.test(AGE_BY_SITE_corr~ PAPA_p4_adhd, data=NEW)
kruskal.test(QuintSoc_w ~ PAPA_p4_adhd, data=NEW)
kruskal.test(birth_size_percent2_x ~ PAPA_p4_adhd, data=NEW)
kruskal.test(birth_wt_grams ~ PAPA_p4_adhd, data=NEW)
kruskal.test(gestation_age_wks ~ PAPA_p4_adhd, data=NEW)
kruskal.test(Pren_CESD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(HWB6_CESD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(HWB12_CESD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(HWB24_CESD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(HWB36_CESD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(HWB48_CESD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(HWB60_CESD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(HWB72_CESD ~ PAPA_p4_adhd, data=NEW)
kruskal.test(auc_post_cesd ~ PAPA_p4_adhd, data=NEW)

kruskal.test(B_DRD2_cat ~ PAPA_p4_adhd, data=NEW)
kruskal.test(B_DRD3_rs6280_cat ~ PAPA_p4_adhd, data=NEW)
kruskal.test(B_COMT_cat.y ~ PAPA_p4_adhd, data=NEW)
kruskal.test(B_COMT_rs165599.y ~ PAPA_p4_adhd, data=NEW)

kruskal.test(PRS_0_3_ADHD_child ~ PAPA_p4_adhd, data=NEW)
kruskal.test(PRS_0_001_adhd_child ~ PAPA_p4_adhd, data=NEW)
kruskal.test(PC1 ~ PAPA_p4_adhd, data=NEW)
kruskal.test(PC2 ~ PAPA_p4_adhd, data=NEW)
kruskal.test(PC3 ~ PAPA_p4_adhd, data=NEW)
kruskal.test(birth_wt_grams ~ PAPA_p4_adhd, data=NEW)
kruskal.test(gestation_age_wks ~ PAPA_p4_adhd, data=NEW)
kruskal.test(Smoking_During_Pregnancy ~ PAPA_p4_adhd, data=NEW)

chisq.test(NEW$above_college, NEW$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$gender_male, NEW$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$Sibling, NEW$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(NEW$B_DRD4, NEW$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DAT.y, NEW$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD2_rs1799978.y, NEW$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD1_hap.y, NEW$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

sink('Sibling_differences_kk_adhd_demographics.txt')
kruskal.test(AFactor_A ~ Sibling, data=NEW)
kruskal.test(MFactor_GEN ~ Sibling, data=NEW)
kruskal.test(PAPA_p4nadhd ~ Sibling, data=NEW)
kruskal.test(conners_mother_hyperactivity_score.72m ~ Sibling, data=NEW)
kruskal.test(ADHD ~ Sibling, data=NEW)
kruskal.test(Mother ~ Sibling, data=NEW)
kruskal.test(Father ~ Sibling, data=NEW)
kruskal.test(Teacher ~ Sibling, data=NEW)
kruskal.test(mom_age_birth ~ Sibling, data=NEW)
kruskal.test(Pren_income4 ~ Sibling, data=NEW)
kruskal.test(Pren_life_events ~ Sibling, data=NEW)
kruskal.test(QuintMat_w ~ Sibling, data=NEW)
kruskal.test(Smoking_During_Pregnancy ~ Sibling, data=NEW)
kruskal.test(AGE_BY_SITE_corr~ Sibling, data=NEW)
kruskal.test(QuintSoc_w ~ Sibling, data=NEW)
kruskal.test(birth_size_percent2_x ~ Sibling, data=NEW)
kruskal.test(birth_wt_grams ~ Sibling, data=NEW)
kruskal.test(gestation_age_wks ~ Sibling, data=NEW)
kruskal.test(Pren_CESD ~ Sibling, data=NEW)
kruskal.test(HWB6_CESD ~ Sibling, data=NEW)
kruskal.test(HWB12_CESD ~ Sibling, data=NEW)
kruskal.test(HWB24_CESD ~ Sibling, data=NEW)
kruskal.test(HWB36_CESD ~ Sibling, data=NEW)
kruskal.test(HWB48_CESD ~ Sibling, data=NEW)
kruskal.test(HWB60_CESD ~ Sibling, data=NEW)
kruskal.test(HWB72_CESD ~ Sibling, data=NEW)
kruskal.test(auc_post_cesd ~ Sibling, data=NEW)

kruskal.test(B_DRD2_cat ~ Sibling, data=NEW)
kruskal.test(B_DRD3_rs6280_cat ~ Sibling, data=NEW)
kruskal.test(B_COMT_cat.y ~ Sibling, data=NEW)
kruskal.test(B_COMT_rs165599.y ~ Sibling, data=NEW)

kruskal.test(PRS_0_3_ADHD_child ~ Sibling, data=NEW)
kruskal.test(PRS_0_001_adhd_child ~ Sibling, data=NEW)
kruskal.test(PC1 ~ Sibling, data=NEW)
kruskal.test(PC2 ~ Sibling, data=NEW)
kruskal.test(PC3 ~ Sibling, data=NEW)
kruskal.test(birth_wt_grams ~ Sibling, data=NEW)
kruskal.test(gestation_age_wks ~ Sibling, data=NEW)
kruskal.test(Alcohol_During_Pregnancy ~ Sibling, data=NEW)

chisq.test(NEW$above_college, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$gender_male, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$PAPA_p4_adhd, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(NEW$B_DRD4, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DAT.y, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD2_rs1799978.y, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD1_hap.y, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

sink('Site_differences_kk_adhd_demographics.txt')
kruskal.test(AFactor_A ~ Hamilton, data=NEW)
kruskal.test(MFactor_GEN ~ Hamilton, data=NEW)
kruskal.test(PAPA_p4nadhd ~ Hamilton, data=NEW)
kruskal.test(conners_mother_hyperactivity_score.72m ~ Hamilton, data=NEW)
kruskal.test(ADHD ~ Hamilton, data=NEW)
kruskal.test(Mother ~ Hamilton, data=NEW)
kruskal.test(Father ~ Hamilton, data=NEW)
kruskal.test(Teacher ~ Hamilton, data=NEW)
kruskal.test(mom_age_birth ~ Hamilton, data=NEW)
kruskal.test(Pren_income4 ~ Hamilton, data=NEW)
kruskal.test(Pren_life_events ~ Hamilton, data=NEW)
kruskal.test(QuintMat_w ~ Hamilton, data=NEW)
kruskal.test(AGE_BY_SITE_corr~ Hamilton, data=NEW)
kruskal.test(QuintSoc_w ~ Hamilton, data=NEW)
kruskal.test(birth_size_percent2_x ~ Hamilton, data=NEW)
kruskal.test(birth_wt_grams ~ Hamilton, data=NEW)
kruskal.test(gestation_age_wks ~ Hamilton, data=NEW)
kruskal.test(Pren_CESD ~ Hamilton, data=NEW)
kruskal.test(HWB6_CESD ~ Hamilton, data=NEW)
kruskal.test(HWB12_CESD ~ Hamilton, data=NEW)
kruskal.test(HWB24_CESD ~ Hamilton, data=NEW)
kruskal.test(HWB36_CESD ~ Hamilton, data=NEW)
kruskal.test(HWB48_CESD ~ Hamilton, data=NEW)
kruskal.test(HWB60_CESD ~ Hamilton, data=NEW)
kruskal.test(HWB72_CESD ~ Hamilton, data=NEW)
kruskal.test(auc_post_cesd ~ Hamilton, data=NEW)

kruskal.test(B_DRD2_cat ~ Hamilton, data=NEW)
kruskal.test(B_DRD3_rs6280_cat ~ Hamilton, data=NEW)
kruskal.test(B_COMT_cat.y ~ Hamilton, data=NEW)
kruskal.test(B_COMT_rs165599.y ~ Hamilton, data=NEW)

kruskal.test(PRS_0_3_ADHD_child ~ Hamilton, data=NEW)
kruskal.test(PRS_0_001_adhd_child ~ Hamilton, data=NEW)
kruskal.test(PC1 ~ Hamilton, data=NEW)
kruskal.test(PC2 ~ Hamilton, data=NEW)
kruskal.test(PC3 ~ Hamilton, data=NEW)
kruskal.test(birth_wt_grams ~ Hamilton, data=NEW)
kruskal.test(gestation_age_wks ~ Hamilton, data=NEW)
kruskal.test(Alcohol_During_Pregnancy ~ Hamilton, data=NEW)
kruskal.test(Smoking_During_Pregnancy ~ Hamilton, data=NEW)

chisq.test(NEW$above_college, NEW$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$gender_male, NEW$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$PAPA_p4_adhd, NEW$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(NEW$B_DRD4, NEW$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DAT.y, NEW$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD2_rs1799978.y, NEW$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD1_hap.y, NEW$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

NEW$ADHD_d<-dicho(
  NEW$ADHD,
  dich.by = "median",
  as.num = TRUE,
  var.label = NULL,
  val.labels = NULL,
  append = TRUE,
  suffix = "_d"
)
head(NEW$ADHD_d)

sink('ADHD_factor_differences_kk_adhd_demographics.txt')
kruskal.test(AFactor_A ~ ADHD_d, data=NEW)
kruskal.test(MFactor_GEN ~ ADHD_d, data=NEW)
kruskal.test(PAPA_p4nadhd ~ ADHD_d, data=NEW)
kruskal.test(conners_mother_hyperactivity_score.72m ~ ADHD_d, data=NEW)
kruskal.test(ADHD ~ ADHD_d, data=NEW)
kruskal.test(Mother ~ ADHD_d, data=NEW)
kruskal.test(Father ~ ADHD_d, data=NEW)
kruskal.test(Teacher ~ ADHD_d, data=NEW)
kruskal.test(mom_age_birth ~ ADHD_d, data=NEW)
kruskal.test(Pren_income4 ~ ADHD_d, data=NEW)
kruskal.test(Pren_life_events ~ ADHD_d, data=NEW)
kruskal.test(QuintMat_w ~ ADHD_d, data=NEW)
kruskal.test(AGE_BY_SITE_corr~ ADHD_d, data=NEW)
kruskal.test(QuintSoc_w ~ ADHD_d, data=NEW)
kruskal.test(birth_size_percent2_x ~ ADHD_d, data=NEW)
kruskal.test(birth_wt_grams ~ ADHD_d, data=NEW)
kruskal.test(gestation_age_wks ~ ADHD_d, data=NEW)
kruskal.test(Pren_CESD ~ ADHD_d, data=NEW)
kruskal.test(HWB6_CESD ~ ADHD_d, data=NEW)
kruskal.test(HWB12_CESD ~ ADHD_d, data=NEW)
kruskal.test(HWB24_CESD ~ ADHD_d, data=NEW)
kruskal.test(HWB36_CESD ~ ADHD_d, data=NEW)
kruskal.test(HWB48_CESD ~ ADHD_d, data=NEW)
kruskal.test(HWB60_CESD ~ ADHD_d, data=NEW)
kruskal.test(HWB72_CESD ~ ADHD_d, data=NEW)
kruskal.test(auc_post_cesd ~ ADHD_d, data=NEW)

kruskal.test(B_DRD2_cat ~ ADHD_d, data=NEW)
kruskal.test(B_DRD3_rs6280_cat ~ ADHD_d, data=NEW)
kruskal.test(B_COMT_cat.y ~ ADHD_d, data=NEW)
kruskal.test(B_COMT_rs165599.y ~ ADHD_d, data=NEW)

kruskal.test(PRS_0_3_ADHD_child ~ ADHD_d, data=NEW)
kruskal.test(PRS_0_001_adhd_child ~ ADHD_d, data=NEW)
kruskal.test(PC1 ~ ADHD_d, data=NEW)
kruskal.test(PC2 ~ ADHD_d, data=NEW)
kruskal.test(PC3 ~ ADHD_d, data=NEW)
kruskal.test(birth_wt_grams ~ ADHD_d, data=NEW)
kruskal.test(gestation_age_wks ~ ADHD_d, data=NEW)
kruskal.test(Alcohol_During_Pregnancy ~ ADHD_d, data=NEW)
kruskal.test(Smoking_During_Pregnancy ~ ADHD_d, data=NEW)

chisq.test(NEW$gender_male, NEW$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$above_college, NEW$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$Sibling, NEW$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$PAPA_p4_adhd, NEW$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(NEW$B_DRD4, NEW$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DAT.y, NEW$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD2_rs1799978.y, NEW$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD1_hap.y, NEW$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

