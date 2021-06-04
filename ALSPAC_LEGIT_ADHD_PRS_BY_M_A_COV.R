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

sink('omega_adhd_model_weigh.txt')
semTools::reliability(ADHD.fit3, return.total = TRUE, dropSingle = TRUE,
                      omit.imps = c("no.conv", "no.se"))
lavaan::lavInspect(ADHD.fit3, "cov.lv")
sink()

########predict factor
ADHD_scores <- lavaan::lavPredict(ADHD.fit3)

########merge items with factor
NEW2 <- cbind(ADHDfactor.data2, ADHD_scores)

########merge items with id
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

########merge items with factor
NEW3 <- cbind(NEW2, ADHDfactor.data3["ID"])

########save data
readr::write_csv(NEW3, file = "Z:/AWAZANA/MELM/ALSPAC/data/questionnaires_bifactor.csv")


############need to merge this data to ALSPAC

########adding dawba

NEW = read.csv("Z:/AWAZANA/MELM/ALSPAC/data/questionnaires_bifactor.csv")

NEW2 <- cbind(NEW, ADHDfactor.data2)

ALSPAC = read.csv("Z:/AWAZANA/AWAZANA_LAB/ALSPAC/ALSPAC_all_combined_newPRS_Feb_2021Ruttercleanup.csv")

# subsample_one_item_present_old <- ALSPAC[which((rowMeans(is.na(ALSPAC[c("ID","fh6870","tb8600", "kj645", "FJCI200", "FJCI201", "f179", "f180", "e375", "e376a", "g274", "g275", "h184a","h184b")]))) <= 2/4),]
# dim(subsample_one_item_present_old)

old <- ALSPAC[c("ID","fh6870","tb8600", "kj645", "FJCI200", "FJCI201", "f179", "f180", "e375", "e376a", "g274", "g275", "h184a","h184b", "e391","f201")]
str(old)

NEW3 <- merge(NEW2, old, by = "ID", all=TRUE)

NEW3$e391[NEW3$e391 == -1] = NA 
NEW3$f201[NEW3$f201 == -1] = NA

######average depression based on epds at 8 weeks and 8 months
NEW3$avg_matdep_post = rowMeans(NEW3[, c("e391", "f201")], na.rm=TRUE)

######calculate AUC for depression on CCEI during first 3 years postbirth
t1 = c(2,8,21,33)

# initialize new AUC variable
NEW3$ccei_post_auc= rep(0,NROW(NEW3))

for (i in 1:NROW(NEW3)){
  y = NEW3[i,c("f180", "e376a","g275","h184b")]
  print(y)
  whichismissing = is.na(y)
  print(whichismissing)
  print(sum(whichismissing))
  print(!whichismissing)
  print(t1[!whichismissing])
  print(y[!whichismissing])
  print('----------')
  if(sum(whichismissing) < 4) NEW3$ccei_post_auc[i] = DescTools::AUC(t1[!whichismissing],y[!whichismissing])
  else NEW3$ccei_post_auc[i] = NA
}

NEW3$ccei_post_auc2= rep(0,NROW(NEW3))

for (i in 1:NROW(NEW3)){
  y = NEW3[i,c("f180", "e376a","g275","h184b")]
  print(y)
  whichismissing = is.na(y)
  print(whichismissing)
  print(sum(whichismissing))
  print(!whichismissing)
  print(t1[!whichismissing])
  print(y[!whichismissing])
  print('----------')
  if(sum(whichismissing) < 3) NEW3$ccei_post_auc2[i] = DescTools::AUC(t1[!whichismissing],y[!whichismissing])
  else NEW3$ccei_post_auc2[i] = NA
}

NEW3$ccei_post_auc3= rep(0,NROW(NEW3))

for (i in 1:NROW(NEW3)){
  y = NEW3[i,c("f180", "e376a","g275","h184b")]
  print(y)
  whichismissing = is.na(y)
  print(whichismissing)
  print(sum(whichismissing))
  print(!whichismissing)
  print(t1[!whichismissing])
  print(y[!whichismissing])
  print('----------')
  if(sum(whichismissing) < 2) NEW3$ccei_post_auc3[i] = DescTools::AUC(t1[!whichismissing],y[!whichismissing])
  else NEW3$ccei_post_auc3[i] = NA
}

sink('Z:/AWAZANA/MELM/ALSPAC/data/auc_info.txt')
print(summary(NEW3$ccei_post_auc))
print(summary(NEW3$ccei_post_auc2))
print(summary(NEW3$ccei_post_auc3))

str(NEW3$ccei_post_auc)
str(NEW3$ccei_post_auc2)
str(NEW3$ccei_post_auc3)
sink()

readr::write_csv(NEW3, file = "Z:/AWAZANA/MELM/ALSPAC/data/predict_older_dawba.csv")

NEW = read.csv("Z:/AWAZANA/MELM/ALSPAC/data/predict_older_dawba.csv")

# ## Variables you wanted (Marie-Elyse)
# your.data <- alspac_new[c("ID", prs_names_new, "A","M","GPF","sex","avg_age_z","pc1","pc2","pc3","pc4","pat_ed","mat_ed","mat_age")]
# 
# readr::write_csv(your.data, file = "Z:/AWAZANA/MELM/ALSPAC/data/prs_m_a_covariates.csv")

NEW2 = read.csv("Z:/AWAZANA/MELM/ALSPAC/data/prs_covariates.csv")

#merge data questionnaires and prs covariates

NEW3 <- merge(NEW, NEW2, by = "ID", all=TRUE)

NEW4 = read.csv("Z:/AWAZANA/AWAZANA_LAB/ALSPAC/data_alspac_alex_no_imp_with_matpsy.csv")

NEW4$A <- NEW4$fscore_gena 
NEW4$M <- NEW4$MatPsy 

NEW5 <- merge(NEW3, NEW4[c("ID","A","M")], by = "ID", all=TRUE)

readr::write_csv(NEW5, file = "Z:/AWAZANA/MELM/ALSPAC/data/ALSPAC_LEGIT_ADHD_PRS_BY_M_A.csv")

######identify prs threshold
#0.3 from Alex meta analysis paper


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

plot_LEGIT_3way = function(x, gender_range=c(0,1), gender_name="sex", gender_title=c("Male","Female"), cov_values = NULL, gene_quant = c(.1,.9), env_quant = c(.025,.50,.975), outcome_quant = c(.025,.50,.975), cols = c("blue", "deeppink"), ylab="ADHD", xlab="Maternal depression", legtitle="Genetic score", leglab=NULL, xlim= NULL, ylim= NULL, x_at = NULL, y_at = NULL, cex.axis = 0.8, cex.lab=0.8, cex.main=0.8, cex.leg=0.8, legend="topleft", ...){
  # # #changed the gender_name and the order of title since boys are 0 in ALSPAC and girls are 1
  
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


NEW = read.csv("Z:/AWAZANA/MELM/ALSPAC/data/ALSPAC_LEGIT_ADHD_PRS_BY_M_A.csv")

#####
#####
#####adhd_child "5e08","1e08","1e07","1e06","1e05","1e04","1e03","1e02","5e02","1e01","2e01","3e01","4e01","5e01","1"
#####

NEW$std_ADHD_bi = c(scale(NEW$ADHD_bi))
NEW$std_adhd_child_1e05 = c(scale(NEW$adhd_child_1e05))
NEW$std_adhd_child_00001 = c(scale(NEW$adhd_child_00001))
NEW$std_adhd_child_0001 = c(scale(NEW$adhd_child_0001))
NEW$std_adhd_child_001 = c(scale(NEW$adhd_child_001))
NEW$std_adhd_child_005 = c(scale(NEW$adhd_child_005))
NEW$std_adhd_child_01 = c(scale(NEW$adhd_child_01))
NEW$std_adhd_child_02 = c(scale(NEW$adhd_child_02))
NEW$std_adhd_child_03 = c(scale(NEW$adhd_child_03))
NEW$std_adhd_child_04 = c(scale(NEW$adhd_child_04))
NEW$std_adhd_child_05 = c(scale(NEW$adhd_child_05))
NEW$std_adhd_child_1 = c(scale(NEW$adhd_child_1))
NEW$std_A = c(scale(NEW$A))
NEW$std_M = c(scale(NEW$M))
NEW$std_avg_age_z = c(scale(NEW$avg_age_z))
NEW$std_pc1 = c(scale(NEW$pc1))
NEW$std_pc2 = c(scale(NEW$pc2))
NEW$std_pc3 = c(scale(NEW$pc3))
NEW$std_pc4 = c(scale(NEW$pc4))
NEW$std_pat_ed = c(scale(NEW$pat_ed))
NEW$std_mat_ed = c(scale(NEW$mat_ed))
NEW$std_mat_age = c(scale(NEW$mat_age))
NEW$std_fh6870 = c(scale(NEW$fh6870))
NEW$std_tb8600 = c(scale(NEW$tb8600))
NEW$std_kj645 = c(scale(NEW$kj645))
NEW$std_ccei_post_auc = c(scale(NEW$ccei_post_auc))
NEW$std_ccei_post_auc2 = c(scale(NEW$ccei_post_auc2))
NEW$std_ccei_post_auc3 = c(scale(NEW$ccei_post_auc3))

sink('ADHD_DX_characteristic_summary.txt')
SUMMARYVEC <- summary(NEW[ , c("std_ADHD_bi",
                                  "sex",	"std_A", "std_M", 
                                  "std_mat_ed", "std_mat_age", "std_ccei_post_auc", "std_ccei_post_auc2", "std_ccei_post_auc3",
                                  "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")],na.rm = TRUE) 

print(SUMMARYVEC)
sink()

sink('ADHD_DX_characteristic_summary2.txt')

nVEC <- (NEW[ , c("std_ADHD_bi",
                     "sex",	"std_A", "std_M", 
                     "std_mat_ed", "std_mat_age", "std_ccei_post_auc", "std_ccei_post_auc2", "std_ccei_post_auc3",
                     "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")])  
print(nVEC)
sink()

sink('ADHD_DX_characteristic_description3.txt')

describeVEC <- describe(NEW[ , c("std_ADHD_bi",
                                    "sex",	"std_A", "std_M", 
                                    "std_mat_ed", "std_mat_age", "std_ccei_post_auc", "std_ccei_post_auc2", "std_ccei_post_auc3",
                                    "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")]) 

print(describeVEC)
sink()


pdf("correlation_prs_adhd_m_a_covariates.pdf", width=5, height=5, compress=FALSE)
mydata.cor = stats::cor(NEW[, c("std_ADHD_bi",
                                "sex",	"std_A", "std_M", 
                                "std_mat_ed", "std_mat_age", "std_ccei_post_auc", "std_ccei_post_auc2", "std_ccei_post_auc3",
                                "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("std_ADHD_bi",
                                    "sex",	"std_A", "std_M", 
                                    "std_mat_ed", "std_mat_age", "std_ccei_post_auc", "std_ccei_post_auc2", "std_ccei_post_auc3",
                                    "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

# colnames(mydata.cor) <-c("ADHD",	"Mother",	"Father",	"Teacher",
#                          "PRS_0_001_adhd_child", "PC1", 
#                          "PC2", "PC3", 
#                          "Child sex (male=1)",	"Prenatal depression", "Prenatal by Postnatal depression", 
#                          "Postnatal depression", "Maternal age at birth",
#                          "Education","Site")
# 
# rownames(mydata.cor) <-c("ADHD",	"Mother",	"Father",	"Teacher",
#                          "PRS_0_001_adhd_child", "PC1", 
#                          "PC2", "PC3", 
#                          "Child sex (male=1)",	"Prenatal depression", "Prenatal by Postnatal depression", 
#                          "Postnatal depression", "Maternal age at birth",
#                          "Education","Site")
# 
corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.001, .01), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

pdf("correlation_prs_adhd_m_a_covariates_p05.pdf", width=5, height=5, compress=FALSE)
corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.05), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()



######2-way

G = NEW[,c("std_adhd_child_03"), drop=FALSE]
E = NEW[,c("std_A", "std_M"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex + std_ccei_post_auc, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + std postnatal depression + std PC1-4"), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
dev.off()


comp = stats::complete.cases(NEW[c("std_ADHD_bi",
                                   "sex",	"std_A", "std_M", 
                                   "std_mat_ed", "std_mat_age", "std_ccei_post_auc",
                                   "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD_bi - mean(NEW_nomiss$std_ADHD_bi))^2)
R2 = 1 - ssres/sstotal
print(R2)

std_noise = sd(NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))
print(std_noise)


# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex + std_ccei_post_auc, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex + std_ccei_post_auc, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_adhd_child_03"), drop=FALSE]
E = NEW[,c("std_A", "std_M"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + std_ccei_post_auc, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + std postnatal depression + std PC1-4'), cex.main = 0.8)
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
plot_LEGIT_3way(fit_legit_list, gender_name="sex", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + std postnatal depression + std PC1-4",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD_bi",
                                   "sex",	"std_A", "std_M", 
                                   "std_mat_ed", "std_mat_age", "std_ccei_post_auc",
                                   "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD_bi - mean(NEW_nomiss$std_ADHD_bi))^2)
R2 = 1 - ssres/sstotal
print(R2)

std_noise = sd(NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))
print(std_noise)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + std_ccei_post_auc, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + std_ccei_post_auc, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()


######2-way

G = NEW[,c("std_adhd_child_03"), drop=FALSE]
E = NEW[,c("std_A", "std_M", "std_ccei_post_auc"),  drop=FALSE]

pdf("LEGIT_prs_a_m_postnatal_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_postnatal_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + std PC1-4"), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
dev.off()


comp = stats::complete.cases(NEW[c("std_ADHD_bi",
                                   "sex",	"std_A", "std_M", 
                                   "std_mat_ed", "std_mat_age", "std_ccei_post_auc",
                                   "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD_bi - mean(NEW_nomiss$std_ADHD_bi))^2)
R2 = 1 - ssres/sstotal
print(R2)

std_noise = sd(NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))
print(std_noise)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_adhd_child_03"), drop=FALSE]
E = NEW[,c("std_A", "std_M", "std_ccei_post_auc"),  drop=FALSE]

pdf("LEGIT_prs_a_m_postnatal_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_postnatal_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + std postnatal depression + std PC1-4'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_prs_a_m_postnatal_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="sex", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + std PC1-4",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD_bi",
                                   "sex",	"std_A", "std_M", 
                                   "std_mat_ed", "std_mat_age", "std_ccei_post_auc",
                                   "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD_bi - mean(NEW_nomiss$std_ADHD_bi))^2)
R2 = 1 - ssres/sstotal
print(R2)

std_noise = sd(NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))
print(std_noise)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

############ auc made up of 3 time points
######2-way

G = NEW[,c("std_adhd_child_03"), drop=FALSE]
E = NEW[,c("std_A", "std_M"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov_post3.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov_post3.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex + std_ccei_post_auc3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + std postnatal depression + std PC1-4"), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
dev.off()


comp = stats::complete.cases(NEW[c("std_ADHD_bi",
                                   "sex",	"std_A", "std_M", 
                                   "std_mat_ed", "std_mat_age", "std_ccei_post_auc3",
                                   "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD_bi - mean(NEW_nomiss$std_ADHD_bi))^2)
R2 = 1 - ssres/sstotal
print(R2)

std_noise = sd(NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))
print(std_noise)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex + std_ccei_post_auc, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex + std_ccei_post_auc3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_adhd_child_03"), drop=FALSE]
E = NEW[,c("std_A", "std_M"),  drop=FALSE]

pdf("LEGIT_prs_a_m_cov_post3_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_cov_post3_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + std_ccei_post_auc3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + std postnatal depression + std PC1-4'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_prs_a_m_cov_post3_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="sex", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + std postnatal depression + std PC1-4",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD_bi",
                                   "sex",	"std_A", "std_M", 
                                   "std_mat_ed", "std_mat_age", "std_ccei_post_auc3",
                                   "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD_bi - mean(NEW_nomiss$std_ADHD_bi))^2)
R2 = 1 - ssres/sstotal
print(R2)

std_noise = sd(NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))
print(std_noise)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + std_ccei_post_auc3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + std_ccei_post_auc, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

######
######post auc made up of 3 time points
######2-way

G = NEW[,c("std_adhd_child_03"), drop=FALSE]
E = NEW[,c("std_A", "std_M", "std_ccei_post_auc3"),  drop=FALSE]

pdf("LEGIT_prs_a_m_postnatal3_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_postnatal3_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + std PC1-4"), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
dev.off()


comp = stats::complete.cases(NEW[c("std_ADHD_bi",
                                   "sex",	"std_A", "std_M", 
                                   "std_mat_ed", "std_mat_age", "std_ccei_post_auc3",
                                   "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD_bi - mean(NEW_nomiss$std_ADHD_bi))^2)
R2 = 1 - ssres/sstotal
print(R2)

std_noise = sd(NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))
print(std_noise)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4 + sex, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_adhd_child_03"), drop=FALSE]
E = NEW[,c("std_A", "std_M", "std_ccei_post_auc3"),  drop=FALSE]

pdf("LEGIT_prs_a_m_postnatal3_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_a_m_postnatal3_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + std postnatal depression + std PC1-4'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_prs_a_m_postnatal3_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="sex", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + std PC1-4",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD_bi",
                                   "sex",	"std_A", "std_M", 
                                   "std_mat_ed", "std_mat_age", "std_ccei_post_auc3",
                                   "std_adhd_child_03", "std_pc1", "std_pc2", "std_pc3", "std_pc4")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD_bi - mean(NEW_nomiss$std_ADHD_bi))^2)
R2 = 1 - ssres/sstotal
print(R2)

std_noise = sd(NEW_nomiss$std_ADHD_bi - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))
print(std_noise)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD_bi ~ G*E*sex + std_mat_age + std_mat_ed + std_pc1 + std_pc2 + std_pc3 + std_pc4, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
# print(cv_loo)

sink()

