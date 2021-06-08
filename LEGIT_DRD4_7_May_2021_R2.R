setwd("/Users/Marie-Elyse/Downloads")
# NEW = read.csv("MAVAN_48M_and_up_jun2020.csv")
# NEW = read.csv("MAVAN_48M_and_up_jun2020_new.csv")
NEW = read.csv("FEB2021.csv")

#Need to redo with "B_DRD4" 

NEW$PRE_BY_POST = (NEW$auc_post_cesd*NEW$Pren_CESD)

#need to recode
NEW$B_DRD3_rs6280_cat[NEW$baby_DRD3_rs6280_=="C C"] <- "0"
NEW$B_DRD3_rs6280_cat[NEW$baby_DRD3_rs6280_=="T C"] <- "0.5"
NEW$B_DRD3_rs6280_cat[NEW$baby_DRD3_rs6280_=="T T"] <- "1"

NEW$B_DRD3_rs6280_cat <- as.numeric(as.character(NEW$B_DRD3_rs6280_cat))

# B_DRD4 = .;
# if baby_DRD4_exonIII_ NE '2 7' then B_DRD4 = 0;
# if baby_DRD4_exonIII_ NE '3 7' then B_DRD4 = 0;
# if baby_DRD4_exonIII_ NE '4 7' then B_DRD4 = 0;
# if baby_DRD4_exonIII_ NE '5 7' then B_DRD4 = 0;
# if baby_DRD4_exonIII_ NE '6 7' then B_DRD4 = 0;
# if baby_DRD4_exonIII_ NE '7 7' then B_DRD4 = 0;
# if baby_DRD4_exonIII_ NE '7 8' then B_DRD4 = 0;
# if baby_DRD4_exonIII_ = ' ' then B_DRD4 = .;

# 0 needs to 1 and 1 needs to be 0
# seems like it was fine afterall

# oldvalues <- c("0", "1")
# newvalues <- c("1","0")  
# 
# NEW$B_DRD4 <- newvalues[ match(NEW$B_DRD4, oldvalues) ]
# NEW$B_DRD4 <- as.numeric(as.character(NEW$B_DRD4))

# NEW$B_COMT_cat

# if baby_COMT_rs4680_ = 'G G' then B_COMT_cat = 0; Val/Val
# else if baby_COMT_rs4680_ = 'G A' then B_COMT_cat = .50;
# else if baby_COMT_rs4680_ = 'A A' then B_COMT_cat = 1; Met/Met
# else B_COMT_cat = .;

# if baby_COMT_rs4680_ = 'G G' then B_COMT_cat = 1; Val/Val
# else if baby_COMT_rs4680_ = 'G A' then B_COMT_cat = .50;
# else if baby_COMT_rs4680_ = 'A A' then B_COMT_cat = 0; Met/Met
# else B_COMT_cat = .;

oldvalues <- c("0", "0.5", "1")
newvalues <- c("1","0.5", "0")

NEW$B_COMT_cat.y <- newvalues[ match(NEW$B_COMT_cat.y, oldvalues) ]
NEW$B_COMT_cat.y <- as.numeric(as.character(NEW$B_COMT_cat.y))

count(NEW$B_COMT_cat.y)

WOSIBS = read.csv("JUNE2021.csv")

WOSIBS$PRE_BY_POST = (WOSIBS$auc_post_cesd*WOSIBS$Pren_CESD)

#need to recode
WOSIBS$B_DRD3_rs6280_cat[WOSIBS$baby_DRD3_rs6280_=="C C"] <- "0"
WOSIBS$B_DRD3_rs6280_cat[WOSIBS$baby_DRD3_rs6280_=="T C"] <- "0.5"
WOSIBS$B_DRD3_rs6280_cat[WOSIBS$baby_DRD3_rs6280_=="T T"] <- "1"

WOSIBS$B_DRD3_rs6280_cat <- as.numeric(as.character(WOSIBS$B_DRD3_rs6280_cat))

oldvalues <- c("0", "0.5", "1")
newvalues <- c("1","0.5", "0")

WOSIBS$B_COMT_cat.y <- newvalues[ match(NEW$B_COMT_cat.y, oldvalues) ]
WOSIBS$B_COMT_cat.y <- as.numeric(as.character(NEW$B_COMT_cat.y))

count(WOSIBS$B_COMT_cat.y)


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

## Example with 3 way GxExgender
#library(LEGIT)
#data = example_2way(100, sigma = 1, logit = FALSE, seed = NULL)
#gender_male=data$G[,1]
#mydata = cbind(data$data, gender_male)
#genes=data$G[,-1]
#env=data$E
#fit = LEGIT(mydata, genes, env, y ~ G*E*gender_male)
#plot_LEGIT_3way(fit, gender_name="gender_male", cex.leg=1.5)

# plot_LEGIT_3way = function(x, gender_range=c(0,1), gender_name="gender_male", gender_title=c("Female","Male"), cov_values = NULL, gene_quant = c(.025,.50,.975), env_quant = c(.025,.50,.975), outcome_quant = c(.025,.50,.975), cols = c("#3288BD", "#CAB176", "#D53E4F"), ylab="Outcome", xlab="Environment", legtitle="Genetic score", leglab=NULL, xlim= NULL, ylim= NULL, x_at = NULL, y_at = NULL, cex.axis = 0.8, cex.lab=0.8, cex.main=0.8, cex.leg=0.8, legend="topleft", ...){

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

# pdf("LEGIT_dopamine_comt_pre_post_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
# plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
# # title(main = "std ADHD ~ GxExSEX \n + std Mom age at birth + education + site")
# mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
#       side = 3,
#       line = - 2,
#       outer = TRUE)
# dev.off()

#undefined collumn need to fix
# pdf("correlation_factor_gene.pdf", width=5, height=5, compress=FALSE)
# mydata.cor = stats::cor(NEW[, c("ADHD",	"Mother",	"Father",	"Teacher",
#                          "gender_male", "B_HTTLPR_2",	"B_OXT_pep1",	"B_TPH2",	"B_HTR1A",
#                          "B_HTR1B_best",	"B_HTR2A_alt_r",	"B_5HTR1B_rs130058",	
#                          "B_DRD4", "B_DRD4_78",	"B_DAT",
#                          "B_DRD2",	"B_DRD1_hap",
#                          "B_DRD2_rs1799978", "B_DRD3_rs6280", "B_GR_rs10052957", "B_COMT", "B_COMT_rs165599",
#                          "B_BDNF_r", "B_DRD4",	"B_DAT.y",
#                          "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#                          "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#                          "B_COMT_cat.y", "B_COMT_rs165599.y", 
#                          "B_BDNF_cat","PC1",  "PC2", "PC3", "gender_male",	"Pren_CESD", "PRE_BY_POST", 
#                          "auc_post_cesd", "mom_age_birth",
#                          "above_college","Hamilton")],  method ="kendall",  use = "pairwise")
# 
# res1 <- corrplot::cor.mtest(NEW[, c("ADHD",	"Mother",	"Father",	"Teacher",
#                           "gender_male", "B_HTTLPR_2",	"B_OXT_pep1",	"B_TPH2",	"B_HTR1A",
#                           "B_HTR1B_best",	"B_HTR2A_alt_r",	"B_5HTR1B_rs130058",	
#                           "B_DRD4", "B_DRD4_78",	"B_DAT",
#                           "B_DRD2",	"B_DRD1_hap",
#                           "B_DRD2_rs1799978", "B_DRD3_rs6280", "B_GR_rs10052957", "B_COMT", "B_COMT_rs165599",
#                           "B_BDNF_r", "B_DRD4",	"B_DAT.y",
#                           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#                           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#                           "B_COMT_cat.y", "B_COMT_rs165599.y", 
#                           "B_BDNF_cat","PC1", 
#                           "PC2", "PC3", 
#                           "gender_male",	"Pren_CESD", "PRE_BY_POST", 
#                           "auc_post_cesd", "mom_age_birth",
#                           "above_college","Hamilton")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")
# 
# corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
#          sig.level = c(.001, .01), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
#          insig = "label_sig", pch.col = "black", tl.col = "black")
# dev.off()

# DRD1 haplotype (rs686, rs4532 and rs265981), DRD2 (rs1800497), DRD2 (rs1799978), DRD3 (rs6280), DRD4 VNTR, DAT1 VNTR, COMT (rs4680), COMT (rs165599)

pdf("correlation_factor_gene.pdf", width=5, height=5, compress=FALSE)
mydata.cor = stats::cor(NEW[, c("ADHD",	"Mother",	"Father",	"Teacher",
                                "B_DRD4",	"B_DAT.y",
                                "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                "B_COMT_cat.y", "B_COMT_rs165599.y", 
                                "PRS_0_001_adhd_child", "PC1", 
                                "PC2", "PC3", 
                                "gender_male","Pren_CESD", "PRE_BY_POST", 
                                "auc_post_cesd", "mom_age_birth",
                                "above_college","Hamilton")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(NEW[, c("ADHD",	"Mother",	"Father",	"Teacher",
                                    "B_DRD4",	"B_DAT.y",
                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                    "B_COMT_cat.y", "B_COMT_rs165599.y", 
                                    "PRS_0_001_adhd_child", "PC1", 
                                    "PC2", "PC3", 
                                    "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                    "auc_post_cesd", "mom_age_birth",
                                    "above_college","Hamilton")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

colnames(mydata.cor) <-c("ADHD",	"Mother",	"Father",	"Teacher",
                         "DRD4 7R",	"DAT1 10R",
                         "DRD2 (rs1800497)",	"DRD2 (rs1799978)",  
                         "DRD1 haplotype (rs686, rs4532 and rs265981)", "DRD3 (rs6280)", 
                         "COMT (rs4680)", "COMT (rs165599)",
                         "PRS_0_001_adhd_child", "PC1", 
                         "PC2", "PC3", 
                         "Child sex (male=1)",	"Prenatal depression", "Prenatal by Postnatal depression", 
                         "Postnatal depression", "Maternal age at birth",
                         "Education","Site")

rownames(mydata.cor) <-c("ADHD",	"Mother",	"Father",	"Teacher",
                         "DRD4 7R",	"DAT1 10R",
                         "DRD2 (rs1800497)",	"DRD2 (rs1799978)",  
                         "DRD1 haplotype (rs686, rs4532 and rs265981)", "DRD3 (rs6280)", 
                         "COMT (rs4680)", "COMT (rs165599)",
                         "PRS_0_001_adhd_child", "PC1", 
                         "PC2", "PC3", 
                         "Child sex (male=1)",	"Prenatal depression", "Prenatal by Postnatal depression", 
                         "Postnatal depression", "Maternal age at birth",
                         "Education","Site")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.001, .01), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

pdf("correlation_factor_gene_p05.pdf", width=5, height=5, compress=FALSE)
corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.05), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

###WOSIBS

pdf("correlation_factor_gene_WOSIBS.pdf", width=5, height=5, compress=FALSE)
mydata.cor = stats::cor(WOSIBS[, c("ADHD",	"Mother",	"Father",	"Teacher",
                                "B_DRD4",	"B_DAT.y",
                                "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                "B_COMT_cat.y", "B_COMT_rs165599.y", 
                                "PRS_0_001_adhd_child", "PC1", 
                                "PC2", "PC3", 
                                "gender_male","Pren_CESD", "PRE_BY_POST", 
                                "auc_post_cesd", "mom_age_birth",
                                "above_college","Hamilton")],  method ="kendall",  use = "pairwise")

res1 <- corrplot::cor.mtest(WOSIBS[, c("ADHD",	"Mother",	"Father",	"Teacher",
                                    "B_DRD4",	"B_DAT.y",
                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                    "B_COMT_cat.y", "B_COMT_rs165599.y", 
                                    "PRS_0_001_adhd_child", "PC1", 
                                    "PC2", "PC3", 
                                    "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                    "auc_post_cesd", "mom_age_birth",
                                    "above_college","Hamilton")], conf.level = .95, method ="kendall",  alternative = "two.sided", exact=FALSE, use = "pairwise")

colnames(mydata.cor) <-c("ADHD",	"Mother",	"Father",	"Teacher",
                         "DRD4 7R",	"DAT1 10R",
                         "DRD2 (rs1800497)",	"DRD2 (rs1799978)",  
                         "DRD1 haplotype (rs686, rs4532 and rs265981)", "DRD3 (rs6280)", 
                         "COMT (rs4680)", "COMT (rs165599)",
                         "PRS_0_001_adhd_child", "PC1", 
                         "PC2", "PC3", 
                         "Child sex (male=1)",	"Prenatal depression", "Prenatal by Postnatal depression", 
                         "Postnatal depression", "Maternal age at birth",
                         "Education","Site")

rownames(mydata.cor) <-c("ADHD",	"Mother",	"Father",	"Teacher",
                         "DRD4 7R",	"DAT1 10R",
                         "DRD2 (rs1800497)",	"DRD2 (rs1799978)",  
                         "DRD1 haplotype (rs686, rs4532 and rs265981)", "DRD3 (rs6280)", 
                         "COMT (rs4680)", "COMT (rs165599)",
                         "PRS_0_001_adhd_child", "PC1", 
                         "PC2", "PC3", 
                         "Child sex (male=1)",	"Prenatal depression", "Prenatal by Postnatal depression", 
                         "Postnatal depression", "Maternal age at birth",
                         "Education","Site")

corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.001, .01), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()

pdf("correlation_factor_gene_p05_WOSIBS.pdf", width=5, height=5, compress=FALSE)
corrplot::corrplot(mydata.cor, p.mat = res1$p, method = "color", type = "upper",
                   sig.level = c(.05), pch.cex = .8, tl.cex = .6,  tl.srt = 45,
                   insig = "label_sig", pch.col = "black", tl.col = "black")
dev.off()


##BDNF was captured using SNP rs6265 

# baby_BDNF_rs56164415_#too many missing
# B_BDNF_r

#B_DRD3_rs6280_cat	B_BDNF_cat	B_DRD2_cat	B_DRD1_hap.y	B_DRD4	B_DAT.y	B_DRD2_rs1799978.y	B_COMT_cat.y	B_COMT_rs165599.y

#########
#########
#########
#########
#########
#########
#########
######### data selection

all_item_present <- NEW[which((rowMeans(is.na(NEW[c("PSCID","B_DRD4",	"B_DAT.y",
                                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                    "B_COMT_cat.y", "B_COMT_rs165599.y", 
                                                    "PRS_0_001_adhd_child", "PC1", 
                                                    "PC2", "PC3")]))) <= 0/13),]
dim(all_item_present)
#154

gene <- all_item_present[c("PSCID","B_DRD4",	"B_DAT.y",
                           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                           "B_COMT_cat.y", "B_COMT_rs165599.y", 
                           "PRS_0_001_adhd_child", "PC1", 
                           "PC2", "PC3")]
str(gene)

all_item_present <- NEW[which((rowMeans(is.na(NEW[c("PSCID","B_DRD4",	"B_DAT.y",
                                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                    "B_COMT_cat.y", "B_COMT_rs165599.y")]))) <= 0/9),]
dim(all_item_present)
#185

geneb <- all_item_present[c("PSCID","B_DRD4",	"B_DAT.y",
                           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                           "B_COMT_cat.y", "B_COMT_rs165599.y")]
str(geneb)

#WOSIBS

all_item_present <- WOSIBS[which((rowMeans(is.na(WOSIBS[c("PSCID","B_DRD4",	"B_DAT.y",
                                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                    "B_COMT_cat.y", "B_COMT_rs165599.y")]))) <= 0/9),]
dim(all_item_present)
#185

gene_WOSIBS <- all_item_present[c("PSCID","B_DRD4",	"B_DAT.y",
                            "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                            "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                            "B_COMT_cat.y", "B_COMT_rs165599.y")]
str(gene_WOSIBS)

all_item_present <- NEW[which((rowMeans(is.na(NEW[c("PSCID","PRS_0_001_adhd_child", "PC1", 
                                                    "PC2", "PC3")]))) <= 0/5),]
dim(all_item_present)
#236

genec <- all_item_present[c("PSCID","PRS_0_001_adhd_child", "PC1", 
                           "PC2", "PC3")]
str(genec)

some_item_present <- NEW[which((rowMeans(is.na(NEW[c("PSCID","ADHD","PAPA_p4nadhd","conners_mother_hyperactivity_score.72m",
                                                     "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                     "auc_post_cesd", "mom_age_birth",
                                                     "above_college","Hamilton")]))) <= 9/10),]
dim(some_item_present)
#394

environment_covariates <- some_item_present[c("PSCID","ADHD","PAPA_p4nadhd","conners_mother_hyperactivity_score.72m","gender_male", "Pren_CESD", "PRE_BY_POST", 
                                              "auc_post_cesd", "mom_age_birth",
                                              "above_college","Hamilton")]
str(environment_covariates)

#########
#########
#########
#########
#########
#########
#########
######### merge data for imputation

data_miss_gene=merge(environment_covariates, geneb, by=c("PSCID"))
data_miss=merge(environment_covariates, gene, by=c("PSCID"))
data_miss_prs=merge(environment_covariates, genec, by=c("PSCID"))
data_miss_gene_WOSIBS=merge(environment_covariates, gene_WOSIBS, by=c("PSCID"))

#########
#########
#########
#########
#########
#########
#########imputation

# Setting up the variables classes properly
data_miss_gene$gender_male = factor(data_miss_gene$gender_male)
data_miss_gene$Hamilton = factor(data_miss_gene$Hamilton)
data_miss_gene$above_college = factor(data_miss_gene$above_college)
data_miss_gene$B_DRD4 = factor(data_miss_gene$B_DRD4)
data_miss_gene$B_DAT = factor(data_miss_gene$B_DAT.y)
data_miss_gene$B_DRD2_cat = factor(data_miss_gene$B_DRD2_cat)
data_miss_gene$B_COMT_cat = factor(data_miss_gene$B_COMT_cat.y)
data_miss_gene$B_DRD1_hap.y = factor(data_miss_gene$B_DRD1_hap.y)
data_miss_gene$B_DRD2_rs1799978.y = factor(data_miss_gene$B_DRD2_rs1799978.y)
data_miss_gene$B_DRD3_rs6280_cat = factor(data_miss_gene$B_DRD3_rs6280_cat)
data_miss_gene$B_COMT_rs165599.y = factor(data_miss_gene$B_COMT_rs165599.y)

# Setting up the variables classes properly

data_miss_gene_WOSIBS$gender_male = factor(data_miss_gene_WOSIBS$gender_male)
data_miss_gene_WOSIBS$Hamilton = factor(data_miss_gene_WOSIBS$Hamilton)
data_miss_gene_WOSIBS$above_college = factor(data_miss_gene_WOSIBS$above_college)
data_miss_gene_WOSIBS$B_DRD4 = factor(data_miss_gene_WOSIBS$B_DRD4)
data_miss_gene_WOSIBS$B_DAT = factor(data_miss_gene_WOSIBS$B_DAT.y)
data_miss_gene_WOSIBS$B_DRD2_cat = factor(data_miss_gene_WOSIBS$B_DRD2_cat)
data_miss_gene_WOSIBS$B_COMT_cat = factor(data_miss_gene_WOSIBS$B_COMT_cat.y)
data_miss_gene_WOSIBS$B_DRD1_hap.y = factor(data_miss_gene_WOSIBS$B_DRD1_hap.y)
data_miss_gene_WOSIBS$B_DRD2_rs1799978.y = factor(data_miss_gene_WOSIBS$B_DRD2_rs1799978.y)
data_miss_gene_WOSIBS$B_DRD3_rs6280_cat = factor(data_miss_gene_WOSIBS$B_DRD3_rs6280_cat)
data_miss_gene_WOSIBS$B_COMT_rs165599.y = factor(data_miss_gene_WOSIBS$B_COMT_rs165599.y)

# Setting up the variables classes properly
data_miss_prs$gender_male = factor(data_miss_prs$gender_male)
data_miss_prs$Hamilton = factor(data_miss_prs$Hamilton)
data_miss_prs$above_college = factor(data_miss_prs$above_college)

# Setting up the variables classes properly
data_miss$gender_male = factor(data_miss$gender_male)
data_miss$Hamilton = factor(data_miss$Hamilton)
data_miss$above_college = factor(data_miss$above_college)
data_miss$B_DRD4 = factor(data_miss$B_DRD4)
data_miss$B_DAT = factor(data_miss$B_DAT.y)
data_miss$B_DRD2_cat = factor(data_miss$B_DRD2_cat)
data_miss$B_COMT_cat = factor(data_miss$B_COMT_cat.y)
data_miss$B_DRD1_hap.y = factor(data_miss$B_DRD1_hap.y)
data_miss$B_DRD2_rs1799978.y = factor(data_miss$B_DRD2_rs1799978.y)
data_miss$B_DRD3_rs6280_cat = factor(data_miss$B_DRD3_rs6280_cat)
data_miss$B_COMT_rs165599.y = factor(data_miss$B_COMT_rs165599.y)

library(missForest)
set.seed(1234)
data_miss_ = data_miss
data_miss_$PSCID = NULL
data_imp2_ = missForest(data_miss_, verbose = TRUE)$ximp

set.seed(1234)
data_miss_gene_ = data_miss_gene
data_miss_gene_$PSCID = NULL
data_imp2_gene_ = missForest(data_miss_gene_, verbose = TRUE)$ximp

set.seed(1234)
data_miss_gene_WOSIBS_ = data_miss_gene_WOSIBS
data_miss_gene_WOSIBS_$PSCID = NULL
data_imp2_gene_WOSIBS_ = missForest(data_miss_gene_WOSIBS_, verbose = TRUE)$ximp

set.seed(1234)
data_miss_prs_ = data_miss_prs
data_miss_prs_$PSCID = NULL
data_imp2_prs_ = missForest(data_miss_prs_, verbose = TRUE)$ximp

# Data transformations
data_imp2 = as.data.frame(apply(data_imp2_,2,as.numeric))
data_imp2$PSCID = data_miss$PSCID

data_imp2_prs = as.data.frame(apply(data_imp2_prs_,2,as.numeric))
data_imp2_prs$PSCID = data_miss_prs$PSCID

data_imp2_gene = as.data.frame(apply(data_imp2_gene_,2,as.numeric))
data_imp2_gene$PSCID = data_miss_gene$PSCID

data_imp2_gene_WOSIBS = as.data.frame(apply(data_imp2_gene_WOSIBS_,2,as.numeric))
data_imp2_gene_WOSIBS$PSCID = data_imp2_gene_WOSIBS$PSCID

# data_imp2_gene$gender_male = factor(data_imp2_gene$gender_male)
# data_imp2_gene$Hamilton = factor(data_imp2_gene$Hamilton)
# data_imp2_gene$above_college = factor(data_imp2_gene$above_college)
# data_imp2_gene$B_DRD4_78 = factor(data_imp2_gene$B_DRD4)
# data_imp2_gene$B_DAT = factor(data_imp2_gene$B_DAT.y)
# data_imp2_gene$B_DRD2_cat = factor(data_imp2_gene$B_DRD2_cat)
# data_imp2_gene$B_COMT_cat = factor(data_imp2_gene$B_COMT_cat.y)
# data_imp2_gene$B_DRD1_hap.y = factor(data_imp2_gene$B_DRD1_hap.y)
# data_imp2_gene$B_DRD2_rs1799978.y = factor(data_imp2_gene$B_DRD2_rs1799978.y)
# data_imp2_gene$B_DRD3_rs6280_cat = factor(data_imp2_gene$B_DRD3_rs6280_cat)
# data_imp2_gene$B_COMT_rs165599.y = factor(data_imp2_gene$B_COMT_rs165599.y)

# data_imp2_gene$gender_male = (as.integer(data_imp2_gene$gender_male)
# print(class(data_imp2_gene$gender_male))
data_imp2_gene$gender_male = as.integer(data_imp2_gene$gender_male)
data_imp2_gene$Hamilton = as.integer(data_imp2_gene$Hamilton)
data_imp2_gene$above_college = as.integer(data_imp2_gene$above_college)

data_imp2_gene_WOSIBS$gender_male = as.integer(data_imp2_gene_WOSIBS$gender_male)
data_imp2_gene_WOSIBS$Hamilton = as.integer(data_imp2_gene_WOSIBS$Hamilton)
data_imp2_gene_WOSIBS$above_college = as.integer(data_imp2_gene_WOSIBS$above_college)

#########
#########
#########
#########
#########
#########
#########z score

NEW$std_ADHD = c(scale(NEW$ADHD))
NEW$std_PAPA_p4nadhd = c(scale(NEW$PAPA_p4nadhd))
NEW$std_conners_mother_hyperactivity_score.72m = c(scale(NEW$conners_mother_hyperactivity_score.72m))
NEW$std_PRS_0_001_adhd_child = c(scale(NEW$PRS_0_001_adhd_child))
NEW$std_auc_post_cesd = c(scale(NEW$auc_post_cesd))
NEW$std_Pren_CESD = c(scale(NEW$Pren_CESD))
NEW$std_PC1 = c(scale(NEW$PC1))
NEW$std_PC2 = c(scale(NEW$PC2))
NEW$std_PC3 = c(scale(NEW$PC3))
NEW$std_mom_age_birth = c(scale(NEW$mom_age_birth))
NEW$std_PRE_BY_POST = c(scale(NEW$PRE_BY_POST))

WOSIBS$std_ADHD = c(scale(WOSIBS$ADHD))
WOSIBS$std_PAPA_p4nadhd = c(scale(WOSIBS$PAPA_p4nadhd))
WOSIBS$std_conners_mother_hyperactivity_score.72m = c(scale(WOSIBS$conners_mother_hyperactivity_score.72m))
WOSIBS$std_PRS_0_001_adhd_child = c(scale(WOSIBS$PRS_0_001_adhd_child))
WOSIBS$std_auc_post_cesd = c(scale(WOSIBS$auc_post_cesd))
WOSIBS$std_Pren_CESD = c(scale(WOSIBS$Pren_CESD))
WOSIBS$std_PC1 = c(scale(WOSIBS$PC1))
WOSIBS$std_PC2 = c(scale(WOSIBS$PC2))
WOSIBS$std_PC3 = c(scale(WOSIBS$PC3))
WOSIBS$std_mom_age_birth = c(scale(WOSIBS$mom_age_birth))
WOSIBS$std_PRE_BY_POST = c(scale(WOSIBS$PRE_BY_POST))

data_imp2$std_ADHD = c(scale(data_imp2$ADHD))
data_imp2$std_PAPA_p4nadhd = c(scale(data_imp2$PAPA_p4nadhd))
data_imp2$std_conners_mother_hyperactivity_score.72m = c(scale(data_imp2$conners_mother_hyperactivity_score.72m))
data_imp2$std_PRS_0_001_adhd_child = c(scale(data_imp2$PRS_0_001_adhd_child))
data_imp2$std_auc_post_cesd = c(scale(data_imp2$auc_post_cesd))
data_imp2$std_Pren_CESD = c(scale(data_imp2$Pren_CESD))
data_imp2$std_PC1 = c(scale(data_imp2$PC1))
data_imp2$std_PC2 = c(scale(data_imp2$PC2)) 
data_imp2$std_PC3 = c(scale(data_imp2$PC3)) 
data_imp2$std_mom_age_birth = c(scale(data_imp2$mom_age_birth))
data_imp2$std_PRE_BY_POST = c(scale(data_imp2$PRE_BY_POST))

data_imp2_gene$std_ADHD = c(scale(data_imp2_gene$ADHD))
data_imp2_gene$std_PAPA_p4nadhd = c(scale(data_imp2_gene$PAPA_p4nadhd))
data_imp2_gene$std_conners_mother_hyperactivity_score.72m = c(scale(data_imp2_gene$conners_mother_hyperactivity_score.72m))
data_imp2_gene$std_auc_post_cesd = c(scale(data_imp2_gene$auc_post_cesd))
data_imp2_gene$std_Pren_CESD = c(scale(data_imp2_gene$Pren_CESD))
data_imp2_gene$std_mom_age_birth = c(scale(data_imp2_gene$mom_age_birth))
data_imp2_gene$std_PRE_BY_POST = c(scale(data_imp2_gene$PRE_BY_POST))

data_imp2_gene_WOSIBS$std_ADHD = c(scale(data_imp2_gene_WOSIBS$ADHD))
data_imp2_gene_WOSIBS$std_PAPA_p4nadhd = c(scale(data_imp2_gene_WOSIBS$PAPA_p4nadhd))
data_imp2_gene_WOSIBS$std_conners_mother_hyperactivity_score.72m = c(scale(data_imp2_gene_WOSIBS$conners_mother_hyperactivity_score.72m))
data_imp2_gene_WOSIBS$std_auc_post_cesd = c(scale(data_imp2_gene_WOSIBS$auc_post_cesd))
data_imp2_gene_WOSIBS$std_Pren_CESD = c(scale(data_imp2_gene_WOSIBS$Pren_CESD))
data_imp2_gene_WOSIBS$std_mom_age_birth = c(scale(data_imp2_gene_WOSIBS$mom_age_birth))
data_imp2_gene_WOSIBS$std_PRE_BY_POST = c(scale(data_imp2_gene_WOSIBS$PRE_BY_POST))

data_imp2_prs$std_ADHD = c(scale(data_imp2_prs$ADHD))
data_imp2_prs$std_PAPA_p4nadhd = c(scale(data_imp2_prs$PAPA_p4nadhd))
data_imp2_prs$std_conners_mother_hyperactivity_score.72m = c(scale(data_imp2_prs$conners_mother_hyperactivity_score.72m))
data_imp2_prs$std_PRS_0_001_adhd_child = c(scale(data_imp2_prs$PRS_0_001_adhd_child))
data_imp2_prs$std_auc_post_cesd = c(scale(data_imp2_prs$auc_post_cesd))
data_imp2_prs$std_Pren_CESD = c(scale(data_imp2_prs$Pren_CESD))
data_imp2_prs$std_PC1 = c(scale(data_imp2_prs$PC1))
data_imp2_prs$std_PC2 = c(scale(data_imp2_prs$PC2)) 
data_imp2_prs$std_PC3 = c(scale(data_imp2_prs$PC3)) 
data_imp2_prs$std_mom_age_birth = c(scale(data_imp2_prs$mom_age_birth))
data_imp2_prs$std_PRE_BY_POST = c(scale(data_imp2_prs$PRE_BY_POST))

######## demographics

sink('characteristic_description.txt')
print(summary(NEW[, c("ADHD",
                       "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                       "auc_post_cesd", "mom_age_birth",
                       "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                       "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                       "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                       "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                       "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                       "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                       "birth_size_percent2_x","AGE_BY_SITE_corr",
                       "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                       "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")]))

describeVEC <-Hmisc::describe(NEW[, c("ADHD",
                                       "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                       "auc_post_cesd", "mom_age_birth",
                                       "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                       "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                       "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                       "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
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
                                               "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                               "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                               "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                               "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
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
                                             "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                             "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                             "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                             "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
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
                                                  "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                  "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                  "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                  "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
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
                                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
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
                                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                      "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
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
                                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                      "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
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
                                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                                                   "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                   "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                   "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                   "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                   "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], NEW$ADHD_d)
print(psych.describeVEC.ADHD)

sink()

sink('characteristic_description_WOSIBS.txt')
print(summary(WOSIBS[, c("ADHD",
                      "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                      "auc_post_cesd", "mom_age_birth",
                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                      "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                      "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                      "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                      "birth_size_percent2_x","AGE_BY_SITE_corr",
                      "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                      "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")]))

describeVEC_WOSIBS <-Hmisc::describe(WOSIBS[, c("ADHD",
                                       "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                       "auc_post_cesd", "mom_age_birth",
                                       "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                       "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                       "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                       "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                                       "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                       "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                       "birth_size_percent2_x","AGE_BY_SITE_corr",
                                       "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                       "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], transpose = TRUE)
print(describeVEC_WOSIBS)

library(pastecs)
stat.describeVEC_WOSIBS <-pastecs::stat.desc(WOSIBS[, c("ADHD",
                                              "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                              "auc_post_cesd", "mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                                              "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                              "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                              "birth_size_percent2_x","AGE_BY_SITE_corr",
                                              "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                              "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")])
print(stat.describeVEC_WOSIBS)

###best
psych.describeVEC_WOSIBS <-psych::describe(WOSIBS[, c("ADHD",
                                            "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                            "auc_post_cesd", "mom_age_birth",
                                            "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                            "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                            "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                            "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                                            "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                            "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                            "birth_size_percent2_x","AGE_BY_SITE_corr",
                                            "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                            "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")])
print(psych.describeVEC_WOSIBS)
sink()


sink('characteristic_description_by_sex_WOSIBS.txt')
psych.describeVEC.sex_WOSIBS <-psych::describeBy(WOSIBS[, c("ADHD",
                                                  "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                  "auc_post_cesd", "mom_age_birth",
                                                  "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                  "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                  "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                  "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                                                  "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                  "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                  "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                  "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                  "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], WOSIBS$gender_male)
print(psych.describeVEC.sex_WOSIBS)

sink()

sink('characteristic_description_by_site_WOSIBS.txt')
psych.describeVEC.site_WOSIBS <-psych::describeBy(WOSIBS[, c("ADHD",
                                                   "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                   "auc_post_cesd", "mom_age_birth",
                                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                                                   "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                   "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                   "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                   "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                   "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], WOSIBS$Hamilton)
print(psych.describeVEC.site_WOSIBS_WOSIBS)
sink()

sink('characteristic_description_by_adhd_dx_WOSIBS.txt')
psych.describeVEC.ADHD.DX_WOSIBS <-psych::describeBy(WOSIBS[, c("ADHD",
                                                      "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                      "auc_post_cesd", "mom_age_birth",
                                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                      "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                                                      "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                      "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                      "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                      "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                      "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], WOSIBS$PAPA_p4_adhd)
print(psych.describeVEC.ADHD.DX_WOSIBS)
sink()

# sink('characteristic_description_by_sibling_WOSIBS.txt')
# psych.describeVEC.Sibling_WOSIBS <-psych::describeBy(WOSIBS[, c("ADHD",
#                                                       "gender_male",	"Pren_CESD", "PRE_BY_POST", 
#                                                       "auc_post_cesd", "mom_age_birth",
#                                                       "above_college","Hamilton","B_DRD4",	"B_DAT.y",
#                                                       "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#                                                       "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#                                                       "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
#                                                       "Pren_income4", "QuintMat_w", "QuintSoc_w", 
#                                                       "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
#                                                       "birth_size_percent2_x","AGE_BY_SITE_corr",
#                                                       "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
#                                                       "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], WOSIBS$Sibling)
# print(psych.describeVEC.Sibling_WOSIBS)
# sink()

sink('characteristic_description_by_adhd_median_WOSIBS.txt')
WOSIBS$ADHD_d<-sjmisc::dicho(
  WOSIBS$ADHD,
  dich.by = "median",
  as.num = TRUE,
  var.label = NULL,
  val.labels = NULL,
  append = TRUE,
  suffix = "_d"
)

psych.describeVEC.ADHD_WOSIBS <-psych::describeBy(WOSIBS[, c("ADHD",
                                                   "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                                   "auc_post_cesd", "mom_age_birth",
                                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "PRS_0_001_adhd_child", "PC1", "PC2", "PC3",
                                                   "Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                   "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                   "birth_size_percent2_x","AGE_BY_SITE_corr",
                                                   "Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                   "PAPA_p4_adhd", "PAPA_p4nadhd","conners_mother_hyperactivity_score.72m")], WOSIBS$ADHD_d)
print(psych.describeVEC.ADHD_WOSIBS)
sink()

sink('sex_differences_kk_adhd_demographics.txt')
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

#kruskal.test(PRS_0_001_adhd_child ~ Sibling, data=NEW)
#kruskal.test(PC1 ~ Sibling, data=NEW)
#kruskal.test(PC2 ~ Sibling, data=NEW)
#kruskal.test(PC3 ~ Sibling, data=NEW)
kruskal.test(birth_wt_grams ~ Sibling, data=NEW)
kruskal.test(gestation_age_wks ~ Sibling, data=NEW)
kruskal.test(Alcohol_During_Pregnancy ~ Sibling, data=NEW)

chisq.test(NEW$above_college, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$gender_male, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$PAPA_p4_adhd, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(NEW$B_DRD4, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DAT.y, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
#chisq.test(NEW$B_DRD2_rs1799978.y, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(NEW$B_DRD1_hap.y, NEW$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

sink('Site_differences_kk_adhd_demographics.txt')
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

##############################

sink('sex_differences_kk_adhd_demographics_WOSIBS.txt')
kruskal.test(PAPA_p4nadhd ~ gender_male, data=WOSIBS)
kruskal.test(conners_mother_hyperactivity_score.72m ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(ADHD ~ gender_male, data=WOSIBS)
kruskal.test(Mother ~ gender_male, data=WOSIBS)
kruskal.test(Father ~ gender_male, data=WOSIBS)
kruskal.test(Teacher ~ gender_male, data=WOSIBS)
kruskal.test(mom_age_birth ~ gender_male, data=WOSIBS)
kruskal.test(Pren_income4 ~ gender_male, data=WOSIBS)
kruskal.test(Pren_life_events ~ gender_male, data=WOSIBS)
kruskal.test(QuintMat_w ~ gender_male, data=WOSIBS)
kruskal.test(AGE_BY_SITE_corr~ gender_male, data=WOSIBS)
kruskal.test(QuintSoc_w ~ gender_male, data=WOSIBS)
kruskal.test(birth_size_percent2_x ~ gender_male, data=WOSIBS)
kruskal.test(birth_wt_grams ~ gender_male, data=WOSIBS)
kruskal.test(gestation_age_wks ~ gender_male, data=WOSIBS)
kruskal.test(Pren_CESD ~ gender_male, data=WOSIBS)
kruskal.test(HWB6_CESD ~ gender_male, data=WOSIBS)
kruskal.test(HWB12_CESD ~ gender_male, data=WOSIBS)
kruskal.test(HWB24_CESD ~ gender_male, data=WOSIBS)
kruskal.test(HWB36_CESD ~ gender_male, data=WOSIBS)
kruskal.test(HWB48_CESD ~ gender_male, data=WOSIBS)
kruskal.test(HWB60_CESD ~ gender_male, data=WOSIBS)
kruskal.test(HWB72_CESD ~ gender_male, data=WOSIBS)
kruskal.test(auc_post_cesd ~ gender_male, data=WOSIBS)

kruskal.test(B_DRD2_cat ~ gender_male, data=WOSIBS)
kruskal.test(B_DRD3_rs6280_cat ~ gender_male, data=WOSIBS)
kruskal.test(B_COMT_cat.y ~ gender_male, data=WOSIBS)
kruskal.test(B_COMT_rs165599.y ~ gender_male, data=WOSIBS)

kruskal.test(PRS_0_001_adhd_child ~ gender_male, data=WOSIBS)
kruskal.test(PC1 ~ gender_male, data=WOSIBS)
kruskal.test(PC2 ~ gender_male, data=WOSIBS)
kruskal.test(PC3 ~ gender_male, data=WOSIBS)
kruskal.test(birth_wt_grams ~ gender_male, data=WOSIBS)
kruskal.test(gestation_age_wks ~ gender_male, data=WOSIBS)
kruskal.test(Alcohol_During_Pregnancy ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(Smoking_During_Pregnancy ~ gender_male, data=WOSIBS)

chisq.test(WOSIBS$above_college, WOSIBS$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
#chisq.test(WOSIBS$Sibling, WOSIBS$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$PAPA_p4_adhd, WOSIBS$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(WOSIBS$B_DRD4, WOSIBS$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DAT.y, WOSIBS$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DRD2_rs1799978.y, WOSIBS$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DRD1_hap.y, WOSIBS$gender_male, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

sink('ADHD_differences_kk_adhd_demographics_WOSIBS.txt')
kruskal.test(PAPA_p4nadhd ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(ADHD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(Mother ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(Father ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(Teacher ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(mom_age_birth ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(Pren_income4 ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(Pren_life_events ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(QuintMat_w ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(AGE_BY_SITE_corr~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(QuintSoc_w ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(birth_size_percent2_x ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(birth_wt_grams ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(gestation_age_wks ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(Pren_CESD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(HWB6_CESD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(HWB12_CESD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(HWB24_CESD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(HWB36_CESD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(HWB48_CESD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(HWB60_CESD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(HWB72_CESD ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(auc_post_cesd ~ PAPA_p4_adhd, data=WOSIBS)

kruskal.test(B_DRD2_cat ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(B_DRD3_rs6280_cat ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(B_COMT_cat.y ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(B_COMT_rs165599.y ~ PAPA_p4_adhd, data=WOSIBS)

kruskal.test(PRS_0_001_adhd_child ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(PC1 ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(PC2 ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(PC3 ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(birth_wt_grams ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(gestation_age_wks ~ PAPA_p4_adhd, data=WOSIBS)
kruskal.test(Smoking_During_Pregnancy ~ PAPA_p4_adhd, data=WOSIBS)

chisq.test(WOSIBS$above_college, WOSIBS$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$gender_male, WOSIBS$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
#chisq.test(WOSIBS$Sibling, WOSIBS$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(WOSIBS$B_DRD4, WOSIBS$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DAT.y, WOSIBS$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DRD2_rs1799978.y, WOSIBS$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DRD1_hap.y, WOSIBS$PAPA_p4_adhd, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()
# 
# sink('Sibling_differences_kk_adhd_demographics_WOSIBS.txt')
# kruskal.test(PAPA_p4nadhd ~ Sibling, data=WOSIBS)
# kruskal.test(conners_mother_hyperactivity_score.72m ~ Sibling, data=WOSIBS)
# kruskal.test(ADHD ~ Sibling, data=WOSIBS)
# kruskal.test(Mother ~ Sibling, data=WOSIBS)
# kruskal.test(Father ~ Sibling, data=WOSIBS)
# kruskal.test(Teacher ~ Sibling, data=WOSIBS)
# kruskal.test(mom_age_birth ~ Sibling, data=WOSIBS)
# kruskal.test(Pren_income4 ~ Sibling, data=WOSIBS)
# kruskal.test(Pren_life_events ~ Sibling, data=WOSIBS)
# kruskal.test(QuintMat_w ~ Sibling, data=WOSIBS)
# kruskal.test(Smoking_During_Pregnancy ~ Sibling, data=WOSIBS)
# kruskal.test(AGE_BY_SITE_corr~ Sibling, data=WOSIBS)
# kruskal.test(QuintSoc_w ~ Sibling, data=WOSIBS)
# kruskal.test(birth_size_percent2_x ~ Sibling, data=WOSIBS)
# kruskal.test(birth_wt_grams ~ Sibling, data=WOSIBS)
# kruskal.test(gestation_age_wks ~ Sibling, data=WOSIBS)
# kruskal.test(Pren_CESD ~ Sibling, data=WOSIBS)
# kruskal.test(HWB6_CESD ~ Sibling, data=WOSIBS)
# kruskal.test(HWB12_CESD ~ Sibling, data=WOSIBS)
# kruskal.test(HWB24_CESD ~ Sibling, data=WOSIBS)
# kruskal.test(HWB36_CESD ~ Sibling, data=WOSIBS)
# kruskal.test(HWB48_CESD ~ Sibling, data=WOSIBS)
# kruskal.test(HWB60_CESD ~ Sibling, data=WOSIBS)
# kruskal.test(HWB72_CESD ~ Sibling, data=WOSIBS)
# kruskal.test(auc_post_cesd ~ Sibling, data=WOSIBS)
# 
# kruskal.test(B_DRD2_cat ~ Sibling, data=WOSIBS)
# kruskal.test(B_DRD3_rs6280_cat ~ Sibling, data=WOSIBS)
# kruskal.test(B_COMT_cat.y ~ Sibling, data=WOSIBS)
# kruskal.test(B_COMT_rs165599.y ~ Sibling, data=WOSIBS)
# 
# kruskal.test(birth_wt_grams ~ Sibling, data=WOSIBS)
# kruskal.test(gestation_age_wks ~ Sibling, data=WOSIBS)
# kruskal.test(Alcohol_During_Pregnancy ~ Sibling, data=WOSIBS)
# 
# chisq.test(WOSIBS$above_college, WOSIBS$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
# chisq.test(WOSIBS$gender_male, WOSIBS$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
# chisq.test(WOSIBS$PAPA_p4_adhd, WOSIBS$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
# 
# chisq.test(WOSIBS$B_DRD4, WOSIBS$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
# chisq.test(WOSIBS$B_DAT.y, WOSIBS$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
# #chisq.test(WOSIBS$B_DRD2_rs1799978.y, WOSIBS$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
# chisq.test(WOSIBS$B_DRD1_hap.y, WOSIBS$Sibling, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
# 
# sink()

sink('Site_differences_kk_adhd_demographics_WOSIBS.txt')
kruskal.test(PAPA_p4nadhd ~ Hamilton, data=WOSIBS)
kruskal.test(conners_mother_hyperactivity_score.72m ~ Hamilton, data=WOSIBS)
kruskal.test(ADHD ~ Hamilton, data=WOSIBS)
kruskal.test(Mother ~ Hamilton, data=WOSIBS)
kruskal.test(Father ~ Hamilton, data=WOSIBS)
kruskal.test(Teacher ~ Hamilton, data=WOSIBS)
kruskal.test(mom_age_birth ~ Hamilton, data=WOSIBS)
kruskal.test(Pren_income4 ~ Hamilton, data=WOSIBS)
kruskal.test(Pren_life_events ~ Hamilton, data=WOSIBS)
kruskal.test(QuintMat_w ~ Hamilton, data=WOSIBS)
kruskal.test(AGE_BY_SITE_corr~ Hamilton, data=WOSIBS)
kruskal.test(QuintSoc_w ~ Hamilton, data=WOSIBS)
kruskal.test(birth_size_percent2_x ~ Hamilton, data=WOSIBS)
kruskal.test(birth_wt_grams ~ Hamilton, data=WOSIBS)
kruskal.test(gestation_age_wks ~ Hamilton, data=WOSIBS)
kruskal.test(Pren_CESD ~ Hamilton, data=WOSIBS)
kruskal.test(HWB6_CESD ~ Hamilton, data=WOSIBS)
kruskal.test(HWB12_CESD ~ Hamilton, data=WOSIBS)
kruskal.test(HWB24_CESD ~ Hamilton, data=WOSIBS)
kruskal.test(HWB36_CESD ~ Hamilton, data=WOSIBS)
kruskal.test(HWB48_CESD ~ Hamilton, data=WOSIBS)
kruskal.test(HWB60_CESD ~ Hamilton, data=WOSIBS)
kruskal.test(HWB72_CESD ~ Hamilton, data=WOSIBS)
kruskal.test(auc_post_cesd ~ Hamilton, data=WOSIBS)

kruskal.test(B_DRD2_cat ~ Hamilton, data=WOSIBS)
kruskal.test(B_DRD3_rs6280_cat ~ Hamilton, data=WOSIBS)
kruskal.test(B_COMT_cat.y ~ Hamilton, data=WOSIBS)
kruskal.test(B_COMT_rs165599.y ~ Hamilton, data=WOSIBS)

kruskal.test(PRS_0_001_adhd_child ~ Hamilton, data=WOSIBS)
kruskal.test(PC1 ~ Hamilton, data=WOSIBS)
kruskal.test(PC2 ~ Hamilton, data=WOSIBS)
kruskal.test(PC3 ~ Hamilton, data=WOSIBS)
kruskal.test(birth_wt_grams ~ Hamilton, data=WOSIBS)
kruskal.test(gestation_age_wks ~ Hamilton, data=WOSIBS)
kruskal.test(Alcohol_During_Pregnancy ~ Hamilton, data=WOSIBS)
kruskal.test(Smoking_During_Pregnancy ~ Hamilton, data=WOSIBS)

chisq.test(WOSIBS$above_college, WOSIBS$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$gender_male, WOSIBS$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$PAPA_p4_adhd, WOSIBS$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(WOSIBS$B_DRD4, WOSIBS$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DAT.y, WOSIBS$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DRD2_rs1799978.y, WOSIBS$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DRD1_hap.y, WOSIBS$Hamilton, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

WOSIBS$ADHD_d<-dicho(
  WOSIBS$ADHD,
  dich.by = "median",
  as.num = TRUE,
  var.label = NULL,
  val.labels = NULL,
  append = TRUE,
  suffix = "_d"
)
head(WOSIBS$ADHD_d)

sink('ADHD_factor_differences_kk_adhd_demographics_WOSIBS.txt')
kruskal.test(PAPA_p4nadhd ~ ADHD_d, data=WOSIBS)
kruskal.test(conners_mother_hyperactivity_score.72m ~ ADHD_d, data=WOSIBS)
kruskal.test(ADHD ~ ADHD_d, data=WOSIBS)
kruskal.test(Mother ~ ADHD_d, data=WOSIBS)
kruskal.test(Father ~ ADHD_d, data=WOSIBS)
kruskal.test(Teacher ~ ADHD_d, data=WOSIBS)
kruskal.test(mom_age_birth ~ ADHD_d, data=WOSIBS)
kruskal.test(Pren_income4 ~ ADHD_d, data=WOSIBS)
kruskal.test(Pren_life_events ~ ADHD_d, data=WOSIBS)
kruskal.test(QuintMat_w ~ ADHD_d, data=WOSIBS)
kruskal.test(AGE_BY_SITE_corr~ ADHD_d, data=WOSIBS)
kruskal.test(QuintSoc_w ~ ADHD_d, data=WOSIBS)
kruskal.test(birth_size_percent2_x ~ ADHD_d, data=WOSIBS)
kruskal.test(birth_wt_grams ~ ADHD_d, data=WOSIBS)
kruskal.test(gestation_age_wks ~ ADHD_d, data=WOSIBS)
kruskal.test(Pren_CESD ~ ADHD_d, data=WOSIBS)
kruskal.test(HWB6_CESD ~ ADHD_d, data=WOSIBS)
kruskal.test(HWB12_CESD ~ ADHD_d, data=WOSIBS)
kruskal.test(HWB24_CESD ~ ADHD_d, data=WOSIBS)
kruskal.test(HWB36_CESD ~ ADHD_d, data=WOSIBS)
kruskal.test(HWB48_CESD ~ ADHD_d, data=WOSIBS)
kruskal.test(HWB60_CESD ~ ADHD_d, data=WOSIBS)
kruskal.test(HWB72_CESD ~ ADHD_d, data=WOSIBS)
kruskal.test(auc_post_cesd ~ ADHD_d, data=WOSIBS)

kruskal.test(B_DRD2_cat ~ ADHD_d, data=WOSIBS)
kruskal.test(B_DRD3_rs6280_cat ~ ADHD_d, data=WOSIBS)
kruskal.test(B_COMT_cat.y ~ ADHD_d, data=WOSIBS)
kruskal.test(B_COMT_rs165599.y ~ ADHD_d, data=WOSIBS)

kruskal.test(PRS_0_001_adhd_child ~ ADHD_d, data=WOSIBS)
kruskal.test(PC1 ~ ADHD_d, data=WOSIBS)
kruskal.test(PC2 ~ ADHD_d, data=WOSIBS)
kruskal.test(PC3 ~ ADHD_d, data=WOSIBS)
kruskal.test(birth_wt_grams ~ ADHD_d, data=WOSIBS)
kruskal.test(gestation_age_wks ~ ADHD_d, data=WOSIBS)
kruskal.test(Alcohol_During_Pregnancy ~ ADHD_d, data=WOSIBS)
kruskal.test(Smoking_During_Pregnancy ~ ADHD_d, data=WOSIBS)

chisq.test(WOSIBS$gender_male, WOSIBS$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$above_college, WOSIBS$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
#chisq.test(WOSIBS$Sibling, WOSIBS$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$PAPA_p4_adhd, WOSIBS$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

chisq.test(WOSIBS$B_DRD4, WOSIBS$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DAT.y, WOSIBS$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DRD2_rs1799978.y, WOSIBS$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant
chisq.test(WOSIBS$B_DRD1_hap.y, WOSIBS$ADHD_d, correct=FALSE) # if p value above 0.05 variables are independent if below they are dependant

sink()

###########demographics for complete PRS
###########demographics for complete Dopamine


some_item_present <- NEW[which((rowMeans(is.na(NEW[c("PSCID","Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                                     "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                                     "birth_size_percent2_x","AGE_BY_SITE_corr","Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                                     "PAPA_p4_adhd")]))) <= 11/12),]
dim(some_item_present)

additional_variables <- some_item_present[c("PSCID","Pren_income4", "QuintMat_w", "QuintSoc_w", 
                                              "birth_wt_grams", "gestation_age_wks", "Pren_life_events", 
                                              "birth_size_percent2_x","AGE_BY_SITE_corr","Smoking_During_Pregnancy","Alcohol_During_Pregnancy", 
                                              "PAPA_p4_adhd")]

demo_gene=merge(additional_variables, data_miss_gene, by=c("PSCID"))
demo=merge(additional_variables, data_miss, by=c("PSCID"))
demo_prs=merge(additional_variables, data_miss_prs, by=c("PSCID"))

comp = stats::complete.cases(NEW[c("PSCID","ADHD",
                                   "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                   "auc_post_cesd", "mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

sink('characteristic_description_withoutPCs.txt')
print(summary(comp[, c("ADHD",
                       "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                       "auc_post_cesd", "mom_age_birth",
                       "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                       "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                       "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                       "B_COMT_cat.y", "B_COMT_rs165599.y")]))

describeVEC <-Hmisc::describe(comp[, c("ADHD",
                                      "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                      "auc_post_cesd", "mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], transpose = TRUE)
print(describeVEC)

library(pastecs)
stat.describeVEC <-pastecs::stat.desc(comp[, c("ADHD",
                                              "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                              "auc_post_cesd", "mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")])
print(stat.describeVEC)

###best
psych.describeVEC <-psych::describe(comp[, c("ADHD",
                                            "gender_male",	"Pren_CESD", "PRE_BY_POST", 
                                            "auc_post_cesd", "mom_age_birth",
                                            "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                            "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                            "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                            "B_COMT_cat.y", "B_COMT_rs165599.y")])
print(psych.describeVEC)
sink()

#"Pren_income4", "QuintMat_w", "QuintSoc_w", "birth_wt_grams", "gestation_age_wks", "Pren_life_events", "birth_size_percent2_x","AGE_BY_SITE_corr","Smoking_During_Pregnancy","Alcohol_During_Pregnancy",
#"PAPA_p4nadhd", "PAPA_p4_adhd", "conners_mother_hyperactivity_score.72m", "conners_father_hyperactivity_score.72m", "conners_teacher_hyperactivity_score.72m", "conners_mother_hyperactivity_score.60m", "conners_father_hyperactivity_score.60m",

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)


comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

comp = stats::complete.cases(data_imp2_gene[c("std_ADHD",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)

comp = stats::complete.cases(data_imp2_prs[c("std_ADHD",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)

######
######
######
######
######2-way
######
######
######
######
######

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

sink(paste0('LEGIT_dopamine_comt_pre_post_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

# # R2
# ssres = sum((NEW$std_ADHD - predict(fit_legit_list, data, G,E))^2)
# sstotal = sum((NEW$std_ADHD - mean(NEW$std_ADHD))^2)
# R2 = 1 - ssres/sstotal

pdf("LEGIT_dopamine_comt_pre_post_cov.pdf", width=5, height=5, compress=FALSE)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Example with 3 way GxExgender
#library(LEGIT)
#data = example_2way(100, sigma = 1, logit = FALSE, seed = NULL)
#gender_male=data$G[,1]
#mydata = cbind(data$data, gender_male)
#genes=data$G[,-1]
#env=data$E
#fit = LEGIT(mydata, genes, env, y ~ G*E*gender_male)
#plot_LEGIT_3way(fit, gender_name="gender_male", cex.leg=1.5)

#cv_folds = 10

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
# print(summary(cv_loo))
# print(cv_loo)

sink()

# pdf("LEGIT_dopamine_comt_pre_post_cov_facet_test.pdf", width=5, height=5, compress=FALSE)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character("std ADHD ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
# legend("topleft",
#        title = "Genetic makeup",
#        c("gene 0", "gene 1"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# save.image(file='AUG_2020.RData')
# dev.off()
# dev.off()
# sink()
# dev.off()

######
###### gene + pc

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_pc.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ G*E + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_pc.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()


######3-way

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_sex.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way + pc

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_sex_pc.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_sex_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_sex_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_pc_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
# title(main = "std ADHD ~ GxExSEX \n + std Mom age at birth + education + site")
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST pc

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_pc.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_SEX.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_SEX.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_SEX.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST pc

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_pc.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_pc_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######2-way

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_prs_pre_post_cov.png'))
pdf("LEGIT_prs_pre_post_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_pre_post_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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


comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_prs_pre_post_cov_sex.png'))
pdf("LEGIT_prs_pre_post_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_pre_post_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3'), cex.main = 0.8)
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

pdf("LEGIT_prs_pre_post_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######2-way PREXPOST

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_prs_pre_post_INT.png'))
pdf("LEGIT_prs_pre_post_INT.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_pre_post_INT.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######3-way/4-way PREXPOST

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_prs_pre_post_INT_sex.png'))
pdf("LEGIT_prs_pre_post_INT_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_pre_post_INT_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

pdf("LEGIT_prs_pre_post_INT_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

##############
##############
##############
##############
##############
##############
##############
##############
# imputation
######
######### gene imp

G = data_imp2_gene[,c("B_DRD4",	"B_DAT.y",
                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                      "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_imp_gene.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_imp_gene.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_imp_gene.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
#Error in glm.fit(x = numeric(0), y = numeric(0), weights = NULL, start = NULL,  :
#object 'fit' not found
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_imp_gene.RData')
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_ADHD",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_ADHD - mean(data_imp2_gene_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_ADHD ~ std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######### all imp

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_imp_all.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_ADHD - mean(data_imp2_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_imp_all.RData')
dev.off()

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
###### gene + pc imp

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_pc_imp.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_pc_imp.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_ADHD - mean(data_imp2_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way imp

G = data_imp2_gene[,c("B_DRD4",	"B_DAT.y",
                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                      "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_sex_imp.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_ADHD",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_ADHD - mean(data_imp2_gene_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way imp all

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_imp_all_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_ADHD - mean(data_imp2_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()


######3-way imp all + pc

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.png'))
pdf("LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_ADHD - mean(data_imp2_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp

G = data_imp2_gene[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_imp.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_ADHD",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_ADHD - mean(data_imp2_gene_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp all

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_imp_all.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_imp_all_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_ADHD - mean(data_imp2_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp all + pc

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_pc_imp.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_ADHD - mean(data_imp2_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp

G = data_imp2_gene[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_imp.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_ADHD",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_ADHD - mean(data_imp2_gene_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp all

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_ADHD - mean(data_imp2_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp all pc

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.png'))
pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_ADHD",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_ADHD - mean(data_imp2_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT_cv(data=data_imp2, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######2-way imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_prs_pre_post_cov_imp.png'))
pdf("LEGIT_prs_pre_post_cov_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_pre_post_cov_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

comp = stats::complete.cases(data_imp2_prs[c("std_ADHD",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_ADHD - mean(data_imp2_prs_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######3-way imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('LEGIT_prs_pre_post_cov_sex_imp.png'))
pdf("LEGIT_prs_pre_post_cov_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_pre_post_cov_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3'), cex.main = 0.8)
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

pdf("LEGIT_prs_pre_post_cov_sex_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_prs[c("std_ADHD",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_ADHD - mean(data_imp2_prs_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_prs_pre_post_INT_imp.png'))
pdf("LEGIT_prs_pre_post_INT_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_pre_post_INT_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

comp = stats::complete.cases(data_imp2_prs[c("std_ADHD",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_ADHD - mean(data_imp2_prs_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######3-way/4-way PREXPOST imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('LEGIT_prs_pre_post_INT_sex_imp.png'))
pdf("LEGIT_prs_pre_post_INT_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_prs_pre_post_INT_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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


pdf("LEGIT_prs_pre_post_INT_sex_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_prs[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_ADHD - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_ADHD - mean(data_imp2_prs_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_ADHD ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######
######
######PAPA outcome
######
######
######
######
######2-way

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
###### gene + pc

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_pc.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ G*E + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_pc.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way + pc

G = NEW[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_pc.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_pc_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST pc

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_pc.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)
# 
# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST pc

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_pc.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######2-way

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_prs_pre_post_cov.png'))
pdf("PAPA_LEGIT_prs_pre_post_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_prs_pre_post_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_prs_pre_post_cov_sex.png'))
pdf("PAPA_LEGIT_prs_pre_post_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_prs_pre_post_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std PAPA ADHD ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######2-way PREXPOST

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_prs_pre_post_INT.png'))
pdf("PAPA_LEGIT_prs_pre_post_INT.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_prs_pre_post_INT.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST) \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######3-way/4-way PREXPOST

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_prs_pre_post_INT_sex.png'))
pdf("PAPA_LEGIT_prs_pre_post_INT_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_prs_pre_post_INT_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_prs_pre_post_INT_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_PAPA_p4nadhd",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

########################

G = data_imp2_gene[,c("B_DRD4",	"B_DAT.y",
                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                      "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_imp_gene.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_imp_gene.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_imp_gene.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, maxiter=1000, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
#Error in glm.fit(x = numeric(0), y = numeric(0), weights = NULL, start = NULL,  :
#object 'fit' not found
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ G*E + \n std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_imp_gene.RData')
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_PAPA_p4nadhd",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_PAPA_p4nadhd - mean(data_imp2_gene_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######### all imp

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_imp_all.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_imp_all.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_PAPA_p4nadhd",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_PAPA_p4nadhd - mean(data_imp2_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
###### gene + pc imp

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_pc_imp.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_pc_imp.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_PAPA_p4nadhd",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_PAPA_p4nadhd - mean(data_imp2_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way imp

G = data_imp2_gene[,c("B_DRD4",	"B_DAT.y",
                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                      "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_imp.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_PAPA_p4nadhd",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_PAPA_p4nadhd - mean(data_imp2_gene_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way imp all

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_imp_all_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_PAPA_p4nadhd",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_PAPA_p4nadhd - mean(data_imp2_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way imp all + pc

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_PAPA_p4nadhd",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_PAPA_p4nadhd - mean(data_imp2_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp

G = data_imp2_gene[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_imp.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_PAPA_p4nadhd",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_PAPA_p4nadhd - mean(data_imp2_gene_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp all

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_imp_all.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD")
title(as.character("std PAPA ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_PAPA_p4nadhd",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_PAPA_p4nadhd - mean(data_imp2_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp all + pc

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_pc_imp.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_PAPA_p4nadhd",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_PAPA_p4nadhd - mean(data_imp2_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp

G = data_imp2_gene[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_imp.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_PAPA_p4nadhd",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_PAPA_p4nadhd - mean(data_imp2_gene_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp all

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_PAPA_p4nadhd",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_PAPA_p4nadhd - mean(data_imp2_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp all pc

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.png'))
pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()


comp = stats::complete.cases(data_imp2[c("std_PAPA_p4nadhd",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_PAPA_p4nadhd - mean(data_imp2_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######2-way imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_prs_pre_post_cov_imp.png'))
pdf("PAPA_LEGIT_prs_pre_post_cov_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_prs_pre_post_cov_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

comp = stats::complete.cases(data_imp2_prs[c("std_PAPA_p4nadhd",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_PAPA_p4nadhd - mean(data_imp2_prs_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######3-way imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_prs_pre_post_cov_sex_imp.png'))
pdf("PAPA_LEGIT_prs_pre_post_cov_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_prs_pre_post_cov_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std PAPA ADHD ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(data_imp2_prs[c("std_PAPA_p4nadhd",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_PAPA_p4nadhd - mean(data_imp2_prs_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######2-way PREXPOST imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_prs_pre_post_INT_imp.png'))
pdf("PAPA_LEGIT_prs_pre_post_INT_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_prs_pre_post_INT_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

pdf("PAPA_LEGIT_prs_pre_post_INT_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_prs[c("std_PAPA_p4nadhd",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_PAPA_p4nadhd - mean(data_imp2_prs_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######3-way/4-way PREXPOST imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('PAPA_LEGIT_prs_pre_post_INT_sex_imp.png'))
pdf("PAPA_LEGIT_prs_pre_post_INT_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_LEGIT_prs_pre_post_INT_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("PAPA_LEGIT_prs_pre_post_INT_sex_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ADHD ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_prs[c("std_PAPA_p4nadhd",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_PAPA_p4nadhd - mean(data_imp2_prs_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######
######
######Conners mother rated hyperactivity at 72m outcome
######
######
######
######
######2-way

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
# dev.off()
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
###### gene + pc

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_pc.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ G*E + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_pc.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way + pc

G = NEW[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_pc.png'))
pdf("Conners hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_pc_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)
# 
# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST pc

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_pc.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)
# 
# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST
G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST pc

G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_pc.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_pc.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_pc.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######2-way

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_prs_pre_post_cov.png'))
pdf("Conners_hyperactivity_LEGIT_prs_pre_post_cov.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_prs_pre_post_cov.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######3-way

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_prs_pre_post_cov_sex.png'))
pdf("Conners_hyperactivity_LEGIT_prs_pre_post_cov_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_prs_pre_post_cov_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std Conners hyperactivity ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######2-way PREXPOST

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_prs_pre_post_INT.png'))
pdf("Conners_hyperactivity_LEGIT_prs_pre_post_INT.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_prs_pre_post_INT.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

pdf("Conners_hyperactivity_LEGIT_prs_pre_post_INT_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######
######3-way/4-way PREXPOST

G = NEW[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_prs_pre_post_INT_sex.png'))
pdf("Conners_hyperactivity_LEGIT_prs_pre_post_INT_sex.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_prs_pre_post_INT_sex.txt'))
fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_prs_pre_post_INT_sex_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(NEW[c("std_conners_mother_hyperactivity_score.72m",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = NEW[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

#############

G = data_imp2_gene[,c("B_DRD4",	"B_DAT.y",
                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                      "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_imp_gene.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_imp_gene.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_imp_gene.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, maxiter=1000, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
#Error in glm.fit(x = numeric(0), y = numeric(0), weights = NULL, start = NULL,  :
#object 'fit' not found
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ G*E + \n std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_imp_gene.RData')
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_conners_mother_hyperactivity_score.72m",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######### all imp

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_imp_all.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_cov_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_imp_all.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_conners_mother_hyperactivity_score.72m",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
###### gene + pc imp

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat",
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_pc_imp.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020_pc_imp.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_conners_mother_hyperactivity_score.72m",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way imp

G = data_imp2_gene[,c("B_DRD4",	"B_DAT.y",
                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                      "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_imp.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_conners_mother_hyperactivity_score.72m",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way imp all

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_imp_all_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_conners_mother_hyperactivity_score.72m",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way imp all + pc

G = data_imp2[,c("B_DRD4",	"B_DAT.y",
                 "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                 "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                 "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_cov_sex_pc_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_conners_mother_hyperactivity_score.72m",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp

G = data_imp2_gene[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_imp.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_conners_mother_hyperactivity_score.72m",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp all

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_imp_all.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_conners_mother_hyperactivity_score.72m",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST imp all + pc

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_pc_imp.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0", "gene 1"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(data_imp2[c("std_conners_mother_hyperactivity_score.72m",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp

G = data_imp2_gene[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_imp.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene[c("std_conners_mother_hyperactivity_score.72m",
                                              "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                              "std_auc_post_cesd", "std_mom_age_birth",
                                              "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                              "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_gene_nomiss = data_imp2_gene[comp, , drop = FALSE]

ssres = sum((data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_gene_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_gene_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_gene))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp all

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_imp_all_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_conners_mother_hyperactivity_score.72m",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp all pc

G = data_imp2[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.png'))
pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_dopamine_comt_pre_post_INT_SEX_pc_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2[c("std_conners_mother_hyperactivity_score.72m",
                                         "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                         "std_auc_post_cesd", "std_mom_age_birth",
                                         "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                         "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                         "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                         "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)

G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_nomiss = data_imp2[comp, , drop = FALSE]

ssres = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######2-way imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_prs_pre_post_cov_imp.png'))
pdf("Conners_hyperactivity_LEGIT_prs_pre_post_cov_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_prs_pre_post_cov_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

comp = stats::complete.cases(data_imp2_prs[c("std_conners_mother_hyperactivity_score.72m",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)
# 
# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######3-way imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_prs_pre_post_cov_sex_imp.png'))
pdf("Conners_hyperactivity_LEGIT_prs_pre_post_cov_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_prs_pre_post_cov_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character('std Conners hyperactivity ~ GxE + sex \n + std Mom age at birth + education + site + std PC1-3'), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

comp = stats::complete.cases(data_imp2_prs[c("std_conners_mother_hyperactivity_score.72m",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######
######2-way PREXPOST imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_prs_pre_post_INT_imp.png'))
pdf("Conners_hyperactivity_LEGIT_prs_pre_post_INT_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_prs_pre_post_INT_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
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

pdf("Conners_hyperactivity_LEGIT_prs_pre_post_INT_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_prs[c("std_conners_mother_hyperactivity_score.72m",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3 + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST imp

G = data_imp2_prs[,c("std_PRS_0_001_adhd_child"), drop=FALSE]
E = data_imp2_prs[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

# png(filename = paste0('Conners_hyperactivity_LEGIT_prs_pre_post_INT_sex_imp.png'))
pdf("Conners_hyperactivity_LEGIT_prs_pre_post_INT_sex_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_LEGIT_prs_pre_post_INT_sex_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
save.image(file='AUG_2020.RData')
dev.off()

pdf("Conners_hyperactivity_LEGIT_prs_pre_post_INT_sex_imp_sexfacet.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n +std Mom age at birth + education + site + std PC1-3",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_prs[c("std_conners_mother_hyperactivity_score.72m",
                                             "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                             "std_auc_post_cesd", "std_mom_age_birth",
                                             "above_college","Hamilton","std_PRS_0_001_adhd_child", "std_PC1", "std_PC2", "std_PC3")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
data_imp2_prs_nomiss = data_imp2_prs[comp, , drop = FALSE]

ssres = sum((data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, data_imp2_prs_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m - mean(data_imp2_prs_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# # Cross-validation 5 times with 5 Folds
# cv_5folds = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=5, cv_folds=5)
# mean(cv_5folds$R2_cv)
# print(mean(cv_5folds$R2_cv))
# print(summary(cv_5folds))
# print(cv_5folds)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_prs, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_PC1 + std_PC2 + std_PC3 + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(data_imp2_prs))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

###################
###################
###################
###################
###################WOSIBS
###################

G = WOSIBS[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

sink(paste0('LEGIT_dopamine_comt_pre_post_cov_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

pdf("LEGIT_dopamine_comt_pre_post_cov_WOSIBS.pdf", width=5, height=5, compress=FALSE)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way

G = WOSIBS[,c("B_DRD4",	"B_DAT.y",
           "B_DRD2_cat",	"B_DRD2_rs1799978.y",
           "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
           "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_sex_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_sexfacet_WOSIBS.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("LEGIT_dopamine_comt_pre_post_INT_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST

G = WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_sexfacet_WOSIBS.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_ADHD",
                                   "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                   "std_auc_post_cesd", "std_mom_age_birth",
                                   "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                   "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                   "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                   "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

###################
###################
###################
###################
###################
################### dopamine only from only one family with imp

G = data_imp2_gene_WOSIBS[,c("B_DRD4",	"B_DAT.y",
              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
              "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

sink(paste0('LEGIT_dopamine_comt_pre_post_cov_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

pdf("LEGIT_dopamine_comt_pre_post_cov_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_ADHD",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way

G = data_imp2_gene_WOSIBS[,c("B_DRD4",	"B_DAT.y",
              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
              "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_cov_sex_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_cov_sex_sexfacet_WOSIBS_imp.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_ADHD",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = data_imp2_gene_WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("LEGIT_dopamine_comt_pre_post_INT_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_ADHD",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST

G = data_imp2_gene_WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_sexfacet_WOSIBS_imp.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_ADHD",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

###################
###################
###################
################### PAPA
################### dopamine only from only one family
###################

G = WOSIBS[,c("B_DRD4",	"B_DAT.y",
              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
              "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

sink(paste0('PAPA_dopamine_comt_pre_post_cov_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

pdf("PAPA_dopamine_comt_pre_post_cov_WOSIBS.pdf", width=5, height=5, compress=FALSE)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_PAPA_p4nadhd",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way

G = WOSIBS[,c("B_DRD4",	"B_DAT.y",
              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
              "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

pdf("PAPA_dopamine_comt_pre_post_cov_sex_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_dopamine_comt_pre_post_cov_sex_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("PAPA_dopamine_comt_pre_post_cov_sex_sexfacet_WOSIBS.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_PAPA_p4nadhd",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("PAPA_dopamine_comt_pre_post_INT_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_dopamine_comt_pre_post_INT_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_PAPA_p4nadhd",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST

G = WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("PAPA_dopamine_comt_pre_post_INT_SEX_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_dopamine_comt_pre_post_INT_SEX_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("PAPA_dopamine_comt_pre_post_INT_SEX_sexfacet_WOSIBS.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_PAPA_p4nadhd",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd- predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

###################
###################
###################
###################
################### PAPA
################### dopamine only from only one family with imputation

G = data_imp2_gene_WOSIBS[,c("B_DRD4",	"B_DAT.y",
                             "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                             "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                             "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

sink(paste0('PAPA_dopamine_comt_pre_post_cov_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

pdf("PAPA_dopamine_comt_pre_post_cov_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_PAPA_p4nadhd",
                                                     "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                                     "std_auc_post_cesd", "std_mom_age_birth",
                                                     "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                     "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                     "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                     "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way

G = data_imp2_gene_WOSIBS[,c("B_DRD4",	"B_DAT.y",
                             "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                             "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                             "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

pdf("PAPA_dopamine_comt_pre_post_cov_sex_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_dopamine_comt_pre_post_cov_sex_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("PAPA_dopamine_comt_pre_post_cov_sex_sexfacet_WOSIBS_imp.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_PAPA_p4nadhd",
                                                     "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                                     "std_auc_post_cesd", "std_mom_age_birth",
                                                     "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                     "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                     "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                     "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = data_imp2_gene_WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("PAPA_dopamine_comt_pre_post_INT_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_dopamine_comt_pre_post_INT_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_PAPA_p4nadhd",
                                                     "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                                     "std_auc_post_cesd", "std_mom_age_birth",
                                                     "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                     "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                     "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                     "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_PAPA_p4nadhd ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST

G = data_imp2_gene_WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("PAPA_dopamine_comt_pre_post_INT_SEX_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('PAPA_dopamine_comt_pre_post_INT_SEX_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std PAPA", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std PAPA ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("PAPA_dopamine_comt_pre_post_INT_SEX_sexfacet_WOSIBS_imp.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std PAPA ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_PAPA_p4nadhd",
                                                     "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                                     "std_auc_post_cesd", "std_mom_age_birth",
                                                     "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                     "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                     "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                     "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_PAPA_p4nadhd - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_PAPA_p4nadhd - mean(NEW_nomiss$std_PAPA_p4nadhd))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_PAPA_p4nadhd ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

###################
###################
###################
################### Conners hyperactivity
################### dopamine only from only one family
###################

G = WOSIBS[,c("B_DRD4",	"B_DAT.y",
              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
              "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_cov_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

pdf("Conners_hyperactivity_dopamine_comt_pre_post_cov_WOSIBS.pdf", width=5, height=5, compress=FALSE)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_conners_mother_hyperactivity_score.72m",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way

G = WOSIBS[,c("B_DRD4",	"B_DAT.y",
              "B_DRD2_cat",	"B_DRD2_rs1799978.y",
              "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
              "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

pdf("Conners_activity_dopamine_comt_pre_post_cov_sex_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_cov_sex_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("Conners_hyperactivity_dopamine_comt_pre_post_cov_sex_sexfacet_WOSIBS.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_conners_mother_hyperactivity_score.72m",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("Conners_hyperactivity_dopamine_comt_pre_post_INT_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_INT_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_conners_mother_hyperactivity_score.72m",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST

G = WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("Conners_hyperactivity_dopamine_comt_pre_post_INT_SEX_WOSIBS.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_INT_SEX_WOSIBS.txt'))
fit_legit_list = LEGIT::LEGIT(data=WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("Conners_hyperactivity_dopamine_comt_pre_post_INT_SEX_sexfacet_WOSIBS.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(WOSIBS[c("std_conners_mother_hyperactivity_score.72m",
                                      "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                      "std_auc_post_cesd", "std_mom_age_birth",
                                      "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                      "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                      "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                      "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

###################
###################
###################
###################
################### Conners
################### dopamine only from only one family with imputation

G = data_imp2_gene_WOSIBS[,c("B_DRD4",	"B_DAT.y",
                             "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                             "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                             "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_cov_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)

pdf("Conners_hyperactivity_dopamine_comt_pre_post_cov_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ G*E + std Mom age at birth + education + site + sex"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_conners_mother_hyperactivity_score.72m",
                                                     "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                                     "std_auc_post_cesd", "std_mom_age_birth",
                                                     "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                     "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                     "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                     "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way

G = data_imp2_gene_WOSIBS[,c("B_DRD4",	"B_DAT.y",
                             "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                             "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                             "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]

pdf("Conners_hyperactivity_dopamine_comt_pre_post_cov_sex_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_cov_sex_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college  + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("Conners_hyperactivity_dopamine_comt_pre_post_cov_sex_sexfacet_WOSIBS_imp.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxExSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_conners_mother_hyperactivity_score.72m",
                                                     "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                                     "std_auc_post_cesd", "std_mom_age_birth",
                                                     "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                     "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                     "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                     "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######2-way PREXPOST

G = data_imp2_gene_WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("Conners_hyperactivity_dopamine_comt_pre_post_INT_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_INT_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topleft",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_conners_mother_hyperactivity_score.72m",
                                                     "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                                     "std_auc_post_cesd", "std_mom_age_birth",
                                                     "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                     "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                     "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                     "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()

######3-way/4-way PREXPOST

G = data_imp2_gene_WOSIBS[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
E = data_imp2_gene_WOSIBS[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]

pdf("Conners_hyperactivity_dopamine_comt_pre_post_INT_SEX_WOSIBS_imp.pdf", width=5, height=5, compress=FALSE)
sink(paste0('Conners_hyperactivity_dopamine_comt_pre_post_INT_SEX_WOSIBS_imp.txt'))
fit_legit_list = LEGIT::LEGIT(data=data_imp2_gene_WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E)
fit_genes = fit_legit_list$genes
fit_env = fit_legit_list$env
print(summary(fit_legit_list))
print(fit_legit_list)
plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std Conners hyperactivity", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
title(as.character("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site"), cex.main = 0.8)
legend("topright",
       title = "Genetic makeup",
       c("gene 0.1", "gene 0.9"),
       col = c("blue", "deeppink"),
       cex = 0.8,
       lwd = 1, 
       lty = 1,
       bty = "n")
dev.off()

pdf("Conners_hyperactivity_dopamine_comt_pre_post_INT_SEX_sexfacet_WOSIBS_imp.pdf", width=10, height=5, compress=FALSE)
plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
mtext("std Conners hyperactivity ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site",                   # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

comp = stats::complete.cases(data_imp2_gene_WOSIBS[c("std_conners_mother_hyperactivity_score.72m",
                                                     "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
                                                     "std_auc_post_cesd", "std_mom_age_birth",
                                                     "above_college","Hamilton","B_DRD4",	"B_DAT.y",
                                                     "B_DRD2_cat",	"B_DRD2_rs1799978.y",
                                                     "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
                                                     "B_COMT_cat.y", "B_COMT_rs165599.y")], G, E)
G_nomiss = G[comp, , drop = FALSE]
E_nomiss = E[comp, , drop = FALSE]
NEW_nomiss = data_imp2_gene_WOSIBS[comp, , drop = FALSE]

ssres = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
sstotal = sum((NEW_nomiss$std_conners_mother_hyperactivity_score.72m - mean(NEW_nomiss$std_conners_mother_hyperactivity_score.72m))^2)
R2 = 1 - ssres/sstotal
print(R2)

# Leave-one-out cross-validation (Note: very slow)
cv_loo = LEGIT::LEGIT_cv(data=data_imp2_gene_WOSIBS, formula = std_conners_mother_hyperactivity_score.72m ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
mean(cv_loo$R2_cv)
print(mean(cv_loo$R2_cv))
print(summary(cv_loo))
print(cv_loo)

sink()


# ######
# ###### gene + pc
# 
# G = NEW[,c("B_DRD4",	"B_DAT.y",
#            "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#            "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#            "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
# E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]
# 
# png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_pc.png'))
# pdf("LEGIT_dopamine_comt_pre_post_cov_pc.pdf", width=5, height=5, compress=FALSE)
# sink(paste0('LEGIT_dopamine_comt_pre_post_cov_pc.txt'))
# fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E)
# fit_genes = fit_legit_list$genes
# fit_env = fit_legit_list$env
# print(summary(fit_legit_list))
# print(fit_legit_list)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character("std ADHD ~ G*E + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
# legend("topleft",
#        title = "Genetic makeup",
#        c("gene 0.1", "gene 0.9"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# save.image(file='AUG_2020_pc.RData')
# dev.off()
# 
# comp = stats::complete.cases(NEW[c("std_ADHD",
#                                    "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
#                                    "std_auc_post_cesd", "std_mom_age_birth",
#                                    "above_college","Hamilton","B_DRD4",	"B_DAT.y",
#                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#                                    "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
# G_nomiss = G[comp, , drop = FALSE]
# E_nomiss = E[comp, , drop = FALSE]
# NEW_nomiss = NEW[comp, , drop = FALSE]
# 
# ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
# sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
# R2 = 1 - ssres/sstotal
# print(R2)
# 
# # Leave-one-out cross-validation (Note: very slow)
# cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + std_PC1 + std_PC2 + std_PC3  + Hamilton + gender_male, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
# mean(cv_loo$R2_cv)
# print(mean(cv_loo$R2_cv))
# print(summary(cv_loo))
# print(cv_loo)
# 
# sink()
# 


# ######3-way + pc
# 
# G = NEW[,c("B_DRD4",	"B_DAT.y",
#            "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#            "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#            "B_COMT_cat.y", "B_COMT_rs165599.y"), drop=FALSE]
# E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd"),  drop=FALSE]
# 
# # png(filename = paste0('LEGIT_dopamine_comt_pre_post_cov_sex_pc.png'))
# pdf("LEGIT_dopamine_comt_pre_post_cov_sex_pc.pdf", width=5, height=5, compress=FALSE)
# sink(paste0('LEGIT_dopamine_comt_pre_post_cov_sex_pc.txt'))
# fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
# fit_genes = fit_legit_list$genes
# fit_env = fit_legit_list$env
# print(summary(fit_legit_list))
# print(fit_legit_list)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0.1", "gene 0.9"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# save.image(file='AUG_2020.RData')
# dev.off()
# 
# pdf("LEGIT_dopamine_comt_pre_post_cov_sex_pc_sexfacet.pdf", width=10, height=5, compress=FALSE)
# plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
# # title(main = "std ADHD ~ GxExSEX \n + std Mom age at birth + education + site")
# mtext("std ADHD ~ GxExSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
#       side = 3,
#       line = - 2,
#       outer = TRUE)
# dev.off()
# 
# comp = stats::complete.cases(NEW[c("std_ADHD",
#                                    "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
#                                    "std_auc_post_cesd", "std_mom_age_birth",
#                                    "above_college","Hamilton","B_DRD4",	"B_DAT.y",
#                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#                                    "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
# G_nomiss = G[comp, , drop = FALSE]
# E_nomiss = E[comp, , drop = FALSE]
# NEW_nomiss = NEW[comp, , drop = FALSE]
# 
# ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
# sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
# R2 = 1 - ssres/sstotal
# print(R2)
# 
# # Leave-one-out cross-validation (Note: very slow)
# cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
# mean(cv_loo$R2_cv)
# print(mean(cv_loo$R2_cv))
# print(summary(cv_loo))
# print(cv_loo)
# 
# sink()

# ######2-way PREXPOST pc
# 
# G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
# E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]
# 
# # png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_pc.png'))
# pdf("LEGIT_dopamine_comt_pre_post_INT_pc.pdf", width=5, height=5, compress=FALSE)
# sink(paste0('LEGIT_dopamine_comt_pre_post_INT_pc.txt'))
# fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
# fit_genes = fit_legit_list$genes
# fit_env = fit_legit_list$env
# print(summary(fit_legit_list))
# print(fit_legit_list)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(0,1), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character("std ADHD ~ GxE(PRExPOST) + sex \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
# legend("topleft",
#        title = "Genetic makeup",
#        c("gene 0.1", "gene 0.9"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# save.image(file='AUG_2020.RData')
# dev.off()
# 
# comp = stats::complete.cases(NEW[c("std_ADHD",
#                                    "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
#                                    "std_auc_post_cesd", "std_mom_age_birth",
#                                    "above_college","Hamilton","B_DRD4",	"B_DAT.y",
#                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#                                    "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
# G_nomiss = G[comp, , drop = FALSE]
# E_nomiss = E[comp, , drop = FALSE]
# NEW_nomiss = NEW[comp, , drop = FALSE]
# 
# ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
# sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
# R2 = 1 - ssres/sstotal
# print(R2)
# 
# # Leave-one-out cross-validation (Note: very slow)
# cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E + std_mom_age_birth + above_college + Hamilton + gender_male + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
# mean(cv_loo$R2_cv)
# print(mean(cv_loo$R2_cv))
# print(summary(cv_loo))
# print(cv_loo)
# 
# sink()

# ######3-way/4-way PREXPOST pc
# 
# G = NEW[,c("B_DRD4","B_DAT.y", "B_DRD2_cat","B_DRD2_rs1799978.y", "B_DRD1_hap.y","B_DRD3_rs6280_cat", "B_COMT_cat.y","B_COMT_rs165599.y"),drop=FALSE]
# E = NEW[,c("std_Pren_CESD", "std_auc_post_cesd", "std_PRE_BY_POST"),  drop=FALSE]
# 
# # png(filename = paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_pc.png'))
# pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_pc.pdf", width=5, height=5, compress=FALSE)
# sink(paste0('LEGIT_dopamine_comt_pre_post_INT_SEX_pc.txt'))
# fit_legit_list = LEGIT::LEGIT(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E)
# fit_genes = fit_legit_list$genes
# fit_env = fit_legit_list$env
# print(summary(fit_legit_list))
# print(fit_legit_list)
# plot_LEGIT(fit_legit_list, env_quant=c(.1,.9), gene_quant=c(.1,.9), col = c("blue", "deeppink"), family="Arial", ylab="std ADHD", cex.axis = 1, cex.lab=1, cex.main=0.8, cex.leg=0.8)
# title(as.character("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3"), cex.main = 0.8)
# legend("topright",
#        title = "Genetic makeup",
#        c("gene 0.1", "gene 0.9"),
#        col = c("blue", "deeppink"),
#        cex = 0.8,
#        lwd = 1, 
#        lty = 1,
#        bty = "n")
# save.image(file='AUG_2020.RData')
# dev.off()
# 
# pdf("LEGIT_dopamine_comt_pre_post_INT_SEX_pc_sexfacet.pdf", width=10, height=5, compress=FALSE)
# plot_LEGIT_3way(fit_legit_list, gender_name="gender_male", cex.leg=0.8)
# mtext("std ADHD ~ GxE(PRExPOST)xSEX \n + std Mom age at birth + education + site + std PC1-3",                   # Add main title
#       side = 3,
#       line = - 2,
#       outer = TRUE)
# dev.off()
# 
# comp = stats::complete.cases(NEW[c("std_ADHD",
#                                    "gender_male",	"std_Pren_CESD", "std_PRE_BY_POST", 
#                                    "std_auc_post_cesd", "std_mom_age_birth",
#                                    "above_college","Hamilton","B_DRD4",	"B_DAT.y",
#                                    "B_DRD2_cat",	"B_DRD2_rs1799978.y",
#                                    "B_DRD1_hap.y", "B_DRD3_rs6280_cat", 
#                                    "B_COMT_cat.y", "B_COMT_rs165599.y", "std_PC1", "std_PC2", "std_PC3")], G, E)
# G_nomiss = G[comp, , drop = FALSE]
# E_nomiss = E[comp, , drop = FALSE]
# NEW_nomiss = NEW[comp, , drop = FALSE]
# 
# ssres = sum((NEW_nomiss$std_ADHD - predict(fit_legit_list, NEW_nomiss, G_nomiss, E_nomiss))^2)
# sstotal = sum((NEW_nomiss$std_ADHD - mean(NEW_nomiss$std_ADHD))^2)
# R2 = 1 - ssres/sstotal
# print(R2)
# 
# 
# # Leave-one-out cross-validation (Note: very slow)
# cv_loo = LEGIT::LEGIT_cv(data=NEW, formula = std_ADHD ~ G*E*gender_male + std_mom_age_birth + above_college + Hamilton + std_PC1 + std_PC2 + std_PC3, genes=G,env=E, cv_iter=1, cv_folds=NROW(NEW))
# mean(cv_loo$R2_cv)
# print(mean(cv_loo$R2_cv))
# print(summary(cv_loo))
# print(cv_loo)
# 
# sink()

