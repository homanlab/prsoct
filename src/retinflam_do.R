# created on Wed Jun 14 12:03:36 2023
# Finn Rabe, <finn dot rabe at bli dot uzh dot ch>
#-----------------------------------------------------------------------

## Load data and create timestamp for makefile
# setwd("retinflam/src")
source("retinflam_load.R")
file.create("../output/R/retinflam_do.Rout")

# run python script that create filtered dataframes and figures
# use_python("/usr/local/bin/python")
source_python("./retinflam_do.py")s

df_norm <- read.csv("../output/data/df_normtest_macula.csv")
df_norm <- subset(df_norm, select = c("Phenotype", "Statistic", "p"))
df_norm$p <- format_p(df_norm$p, stars = )


# function to round p vals
p_round <- function(x, n = 2) {
    max(abs(round(x, n)), abs(signif(x, 1))) * sign(x)
}

# function converts p-values to asterisk
signif.num <- function(x) {
    symnum(x,
        corr = FALSE, na = FALSE, legend = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols = c("***", "**", "*", ".", "")
    )
}

# load filtered overall macula data frame and compute rlm
df_macula <- read.csv("../output/data/df_mlm_macula.csv")
paired_test_result <- t.test(df_macula$Macula_right, df_macula$Macula_left, paired = TRUE)
ttest_macula_stat <- round(paired_test_result$statistic[[1]], digits = 2)
ttest_macula_pval <- paired_test_result$p.value

## Associations between overall macular thickness and PRSSZ while controlling for confounds
prsvar <- "PRSSZ"
cov <- "Age + Age_squared + Sex + Smoking_status + Alcohol_drinker_status + BMI + Townsend_index +Hypertension + Diabetes_mellitus+Image_quality_mean+Genotype_array + Genetic.PC1 + Genetic.PC2 + Genetic.PC3 + Genetic.PC4 + Genetic.PC5 + Genetic.PC6 + Genetic.PC7 + Genetic.PC8 + Genetic.PC9 + Genetic.PC10"
mean_macula <- rlm(as.formula(paste("Macula_mean ~ ", prsvar, "+", cov)),
    data = df_macula, psi = psi.huber
)
n_model_left <- nobs(mean_macula)
pc_dd <- data.frame(summary(mean_macula)$coefficients)
mean_macula_b <- round(pc_dd$Value[2], digits = 2)
mean_macula_se <- round(pc_dd$Std..Error[2], digits = 5)
mean_macula_conflow <- round(confint.default(object = mean_macula, parm = prsvar, level = 0.95)[1], digits = 2)
mean_macula_confhigh <- round(confint.default(object = mean_macula, parm = prsvar, level = 0.95)[2], digits = 2)
fpc <- f.robftest(mean_macula, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
mean_macula_F <- fpcFval[[1]]
mean_macula_p <- round(fpc$p.value, digits = 5)
pc1bp <- bptest(mean_macula)
mean_macula_pcbppval <- format_p(map(pc1bp$p.value, 1)[[1]], stars = FALSE)
mean_macula_pcbpstat <- round(pc1bp$statistic[[1]], digits = 2)
meanmacularfwer <- p.adjust(c(mean_macula_p), method = "holm")

## Same association as above, just with imputed dataset
df_macula <- read.csv("../output/data/df_mlm_macula_miss.csv")
vars_to_test <- c("Macula_mean", "PRSSZ", "Age", "Age_squared", "Sex", "Smoking_status", "Alcohol_drinker_status", "BMI", "Townsend_index", "Hypertension", "Diabetes_mellitus", "Image_quality_mean", "Genotype_array", "Genetic.PC1", "Genetic.PC2", "Genetic.PC3", "Genetic.PC4", "Genetic.PC5", "Genetic.PC6", "Genetic.PC7", "Genetic.PC8", "Genetic.PC9", "Genetic.PC10")
df_miss <- df_macula[, vars_to_test]
imp <- mice(df_miss, m = 5, maxit = 50, method = "pmm", seed = 500)
completed_data <- complete(imp)
mean_macula_imp <- rlm(as.formula(paste("Macula_mean ~ ", prsvar, "+", cov)),
    data = completed_data, psi = psi.huber
)
pc_dd <- data.frame(summary(mean_macula_imp)$coefficients)
mean_macula_imp_b <- round(pc_dd$Value[2], digits = 2)
mean_macula_imp_se <- round(pc_dd$Std..Error[2], digits = 5)
mean_macula_imp_conflow <- round(confint.default(object = mean_macula_imp, parm = prsvar, level = 0.95)[1], digits = 2)
mean_macula_imp_confhigh <- round(confint.default(object = mean_macula_imp, parm = prsvar, level = 0.95)[2], digits = 2)
fpc <- f.robftest(mean_macula_imp, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
mean_macula_imp_F <- fpcFval[[1]]
mean_macula_imp_p <- round(fpc$p.value, digits = 5)

## Associations between left/right macular thickness and PRSSZ while controlling for confounds
cov_left <- "Age + Age_squared + Sex + Smoking_status + Alcohol_drinker_status + BMI + Townsend_index+ Hypertension + Diabetes_mellitus + OCT_quality_left +Genotype_array +Genetic.PC1+Genetic.PC2+ Genetic.PC3+ Genetic.PC4+ Genetic.PC5+Genetic.PC6+ Genetic.PC7+ Genetic.PC8+ Genetic.PC9+ Genetic.PC10"
m_macula_left <- rlm(as.formula(paste("Macula_left ~ ", prsvar, "+", cov_left)),
    data = df_macula, psi = psi.huber
)
pc_dd <- data.frame(summary(m_macula_left)$coefficients)
lmacula_b <- round(pc_dd$Value[2], digits = 2)
lmacula_se <- round(pc_dd$Std..Error[2], digits = 5)
lmacula_conflow <- round(confint.default(object = m_macula_left, parm = prsvar, level = 0.95)[1], digits = 2)
lmacula_confhigh <- round(confint.default(object = m_macula_left, parm = prsvar, level = 0.95)[2], digits = 2)
fpc <- f.robftest(m_macula_left, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
# lmacula_F <- map(fpcFval, 1)[[1]]
lmacula_F <- fpcFval[[1]]
lmacula_p <- round(fpc$p.value, digits = 5)

cov_right <- "Age + Age_squared + Sex + Smoking_status + Alcohol_drinker_status + BMI + Townsend_index+ Hypertension + Diabetes_mellitus + OCT_quality_left +Genotype_array +Genetic.PC1+Genetic.PC2+ Genetic.PC3+ Genetic.PC4+ Genetic.PC5+Genetic.PC6+ Genetic.PC7+ Genetic.PC8+ Genetic.PC9+ Genetic.PC10"
m_macula_right <- rlm(as.formula(paste("Macula_right ~ ", prsvar, "+", cov_right)),
    data = df_macula, psi = psi.huber
)
pc_dd <- data.frame(summary(m_macula_right)$coefficients)
rmacula_b <- round(pc_dd$Value[2], digits = 2)
rmacula_se <- round(pc_dd$Std..Error[2], digits = 5)
rmacula_conflow <- round(confint.default(object = m_macula_right, parm = prsvar, level = 0.95)[1], digits = 2)
rmacula_confhigh <- round(confint.default(object = m_macula_right, parm = prsvar, level = 0.95)[2], digits = 2)
fpc <- f.robftest(m_macula_right, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
# rmacula_F <- map(fpcFval, 1)[[1]]
rmacula_F <- fpcFval[[1]]
rmacula_p <- round(fpc$p.value, digits = 5)
tab_model(m_macula_left, m_macula_right, file = "../output/figures/TableS1_PCRLM_Macula.html", show.fstat = TRUE)
macularfwer <- p.adjust(c(lmacula_p, rmacula_p), method = "holm")

## Associations between outer retina (INL-RPE) mean thickness and PRSSZ while controlling for confounds
df_outerret <- read.csv("../output/data/df_mlm_inner_outer.csv")
prsvar <- "PRSSZ"
cov <- "Age + Age_squared + Sex + Smoking_status + Alcohol_drinker_status + BMI + Townsend_index+ Hypertension + Diabetes_mellitus+Image_quality_mean+Genotype_array + Genetic.PC1 + Genetic.PC2 + Genetic.PC3 + Genetic.PC4 + Genetic.PC5 + Genetic.PC6 + Genetic.PC7 + Genetic.PC8 + Genetic.PC9 + Genetic.PC10"
mean_outerret <- rlm(as.formula(paste("INL_RPE_mean ~ ", prsvar, "+", cov)),
    data = df_outerret, psi = psi.huber
)
pc_dd <- data.frame(summary(mean_outerret)$coefficients)
mean_outerret_b <- round(pc_dd$Value[2], digits = 2)
mean_outerret_se <- round(pc_dd$Std..Error[2], digits = 5)
mean_outerret_conflow <- round(confint.default(object = mean_outerret, parm = prsvar, level = 0.95)[1], digits = 2)
mean_outerret_confhigh <- round(confint.default(object = mean_outerret, parm = prsvar, level = 0.95)[2], digits = 2)
fpc <- f.robftest(mean_outerret, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
mean_outerret_F <- fpcFval[[1]]
mean_outerret_p <- round(fpc$p.value, digits = 5)
adjust_macularfwer <- p.adjust(c(mean_macula_p, mean_outerret_p), method = "holm")


# Regression results for all inner/outer phenotypes
df_filt <- read.csv("../output/data/df_mlm_inner_outer.csv")
ret_list <- c("RNFL_mean", "GC_IPL_mean", "INL_mean", "INL_RPE_mean")
prsvar <- "PRSSZ"
df_phenotypes <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), c("Phenotype", "Coef", "CI(low)", "CI(high)", "SE", "BPstat", "BPp", "F", "p"))
cov <- "Age + Age_squared + Sex + Smoking_status + Alcohol_drinker_status + BMI + Townsend_index+ Hypertension + Diabetes_mellitus+Image_quality_mean+Genotype_array + Genetic.PC1 + Genetic.PC2 + Genetic.PC3 + Genetic.PC4 + Genetic.PC5 + Genetic.PC6 + Genetic.PC7 + Genetic.PC8 + Genetic.PC9 + Genetic.PC10"
pvals_pc <- c()
i <- 0
for (var_name in ret_list) {
    # Create the formula string
    pc_str <- paste(var_name, "~", prsvar, "+", cov)
    # Convert the string to a formula
    pc_obj <- as.formula(pc_str)
    # Run the regression
    m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)

    # pc_intercept<- round(map(m_pc$coefficients["(Intercept)"], 1)[[1]],digits = 2)
    pc_intercept <- round(m_pc$coefficients["(Intercept)"][[1]], digits = 2)
    pc1bp <- bptest(m_pc)
    tab1_pcbppval <- format_p(map(pc1bp$p.value, 1)[[1]], stars = FALSE)
    # tab1_pcbpstat <- round(map(pc1bp$statistic, 1)[[1]],digits = 2)
    tab1_pcbpstat <- round(pc1bp$statistic[[1]], digits = 2)

    pc_dd <- data.frame(summary(m_pc)$coefficients)
    # pc_b <- round(pc_dd$Value[2],digits = 2)
    pc_b <- round(pc_dd$Value[2], digits = 3)
    pc_se <- round(pc_dd$Std..Error[2], digits = 5)
    lgc_conflow <- round(confint.default(object = m_pc, parm = prsvar, level = 0.95)[1], digits = 3)
    lgc_confhigh <- round(confint.default(object = m_pc, parm = prsvar, level = 0.95)[2], digits = 3)

    fpc <- f.robftest(m_pc, var = prsvar)
    fpcFval <- round(fpc$statistic, digits = 2)
    fpcFval <- map(fpcFval, 1)[[1]]
    fpcpval <- round(fpc$p.value, digits = 5)

    pvals_pc <- append(pvals_pc, fpcpval)
    data <- unlist(list(var_name, pc_b, lgc_conflow, lgc_confhigh, pc_se, tab1_pcbpstat, tab1_pcbppval, fpcFval, fpcpval), recursive = FALSE)
    df_it <- as.data.frame(t(data))
    new <- rep(i, ncol(df_it))
    df_phenotypes[nrow(df_it) + i, ] <- df_it
    i <- i + 1
}

## Pathway-specific associations
df_filt <- read.csv("../output/data/df_mlm_macula.csv")
df_prs <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("PRS", "Phenotype", "Coef", "CI(low)", "CI(high)", "SE", "F", "p"))
prs_list <- c("PRSSZ", "PRSSZ_0_001", "PRSSZ_0_05", "PRSSZ_0_1", "PRSSZ_0_2", "PRSSZ_0_3", "PRSSZ_0_4", "PRSSZ_0_5", "PRSSZ_1") # ,"PRSSZ_bestfit")
prs_list_mod <- c("PRSSZ", "PRSSZ 0.001", "PRSSZ 0.05", "PRSSZ 0.1", "PRSSZ 0.2", "PRSSZ 0.3", "PRSSZ 0.4", "PRSSZ 0.5", "PRSSZ 1") # ,"PRSSZ_bestfit")
# prs_list_mod <- rep(prs_list_mod, each = 2)
pc_list <- "Macula_mean"
cov <- "Age + Age_squared + Sex + Smoking_status + Alcohol_drinker_status + BMI + Townsend_index+ Hypertension + Diabetes_mellitus+Image_quality_mean+Genotype_array + Genetic.PC1 + Genetic.PC2 + Genetic.PC3 + Genetic.PC4 + Genetic.PC5 + Genetic.PC6 + Genetic.PC7 + Genetic.PC8 + Genetic.PC9 + Genetic.PC10"
# Loop over the new variable names
i <- 0
pfdr_cols <- c()
for (var_name in prs_list) {
    pvals_pc <- c()
    for (pc in pc_list) {
        # Create the formula string
        pc_str <- paste(pc, "~", var_name, "+", cov)
        # Convert the string to a formula
        pc_obj <- as.formula(pc_str)
        # Run the regression
        m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)

        # pc_intercept<- round(map(m_pc$coefficients["(Intercept)"], 1)[[1]],digits = 2)
        pc_intercept <- round(m_pc$coefficients["(Intercept)"][[1]], digits = 2)
        pc1bp <- bptest(m_pc)
        tab1_pcbppval <- round(map(pc1bp$p.value, 1)[[1]], digits = 2)
        # tab1_pcbpstat <- round(map(pc1bp$statistic, 1)[[1]],digits = 2)
        tab1_pcbpstat <- round(pc1bp$statistic[[1]], digits = 2)

        pc_dd <- data.frame(summary(m_pc)$coefficients)
        pc_b <- round(pc_dd$Value[2], digits = 2)
        pc_se <- round(pc_dd$Std..Error[2], digits = 5)
        lgc_conflow <- round(confint.default(object = m_pc, parm = var_name, level = 0.95)[1], digits = 3)
        lgc_confhigh <- round(confint.default(object = m_pc, parm = var_name, level = 0.95)[2], digits = 3)

        fpc <- f.robftest(m_pc, var = var_name)
        fpcFval <- round(fpc$statistic, digits = 2)
        # fpcFval <- map(fpcFval, 1)[[1]]
        fpcFval <- fpcFval[[1]]
        fpcpval <- round(fpc$p.value, digits = 4)

        pvals_pc <- append(pvals_pc, fpcpval)
        # data <- unlist(list(var_name,pc,pc_b,lgc_conflow,lgc_confhigh,pc_se,tab1_pcbpstat,tab1_pcbppval,fpcFval,fpcpval),recursive = FALSE)
        data <- unlist(list(var_name, pc, pc_b, lgc_conflow, lgc_confhigh, pc_se, fpcFval, fpcpval), recursive = FALSE)
        df_it <- as.data.frame(t(data))
        new <- rep(i, ncol(df_it))
        df_prs[nrow(df_it) + i, ] <- df_it
        i <- i + 1
    }
    pc_fdr <- p.adjust(pvals_pc, method = "holm")
    pfdr_cols <- append(pfdr_cols, pc_fdr)
}
# p_ast <- signif.num(pfdr_cols)
# pfdr_cols_ast <- paste(pfdr_cols,p_ast,sep="")
df_prs$PRS <- prs_list_mod
df_prs$Coef <- as.numeric(df_prs$Coef)
df_prs$F <- as.numeric(df_prs$F)
df_prs["pFWER"] <- pfdr_cols # pfdr_cols_ast


df_filt <- read.csv("../output/data/df_mlm_inner_outer.csv")
df_pathprs <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("GeseaCode", "Pathway", "Phenotype", "Coef", "CI(low)", "CI(high)", "SE", "BPstat", "BPp", "F", "p"))
cpathprs_code <- c(
    "NEUROINFLAM_PRS", "ACUTEINFLAM_PRS", "CHROINFLAM_PRS", "TGFB_PRS", "WNT_PRS", "CATENIN_PRS", "DOPPOSREG_PRS",
    "ABNOVAS_PRS", "CORART_PRS"
)
pc_list <- c("RNFL_mean", "GC_IPL_mean", "INL_mean", "INL_RPE_mean")
pathprs_list <- c("M24927", "M6557", "M15140", "M18933", "M25305", "M17761", "M24111", "M43559", "M36658")
cpathprs_code <- rep(cpathprs_code, each = length(pc_list))
# Loop over the new variable names
i <- 0
pfdr_cols <- c()
for (var_name in pathprs_list) {
    pvals_pc <- c()
    for (pc in pc_list) {
        # Create the formula string
        pc_str <- paste(pc, "~", var_name, "+", cov)
        # Convert the string to a formula
        pc_obj <- as.formula(pc_str)
        # Run the regression
        m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)

        # pc_intercept<- round(map(m_pc$coefficients["(Intercept)"], 1)[[1]],digits = 2)
        pc_intercept <- round(m_pc$coefficients["(Intercept)"][[1]], digits = 2)
        pc1bp <- bptest(m_pc)
        tab1_pcbppval <- round(map(pc1bp$p.value, 1)[[1]], digits = 2)
        # tab1_pcbpstat <- round(map(pc1bp$statistic, 1)[[1]],digits = 2)
        tab1_pcbpstat <- round(pc1bp$statistic[[1]], digits = 2)

        pc_dd <- data.frame(summary(m_pc)$coefficients)
        pc_b <- round(pc_dd$Value[2], digits = 2)
        pc_se <- round(pc_dd$Std..Error[2], digits = 5)
        lgc_conflow <- round(confint.default(object = m_pc, parm = var_name, level = 0.95)[1], digits = 3)
        lgc_confhigh <- round(confint.default(object = m_pc, parm = var_name, level = 0.95)[2], digits = 3)

        fpc <- f.robftest(m_pc, var = var_name)
        fpcFval <- round(fpc$statistic, digits = 2)
        # fpcFval <- map(fpcFval, 1)[[1]]
        fpcFval <- fpcFval[[1]]
        fpcpvalnum <- as.numeric(fpc$p.value)
        fpcpval <- round(fpcpvalnum, digits = 5)

        pvals_pc <- append(pvals_pc, fpcpval)
        data <- unlist(list(var_name, cpathprs_code[i + 1], pc, pc_b, lgc_conflow, lgc_confhigh, pc_se, tab1_pcbpstat, tab1_pcbppval, fpcFval, fpcpval), recursive = FALSE)
        df_it <- as.data.frame(t(data))
        new <- rep(i, ncol(df_it))
        df_pathprs[nrow(df_it) + i, ] <- df_it
        i <- i + 1
    }
    pc_fdr <- p.adjust(pvals_pc, method = "holm")
    pfdr_cols <- append(pfdr_cols, pc_fdr)
}
# p_ast <- signif.num(pfdr_cols)
# pfdr_cols_ast <- paste(pfdr_cols,p_ast,sep="")
df_pathprs$Coef <- as.numeric(df_pathprs$Coef)
df_pathprs$p <- as.numeric(df_pathprs$p)
df_pathprs$F <- as.numeric(df_pathprs$F)
df_pathprs["pFWER"] <- pfdr_cols # pfdr_cols_ast
df_pathprs["Groups"] <- c(
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Developmental gene expression", "Developmental gene expression", "Developmental gene expression",
    "Developmental gene expression", "Developmental gene expression", "Developmental gene expression",
    "Microvasculature", "Microvasculature",
    "Microvasculature", "Microvasculature"
)
df_pathprs_ordered <- df_pathprs[order(df_pathprs$Pathway, df_pathprs$Phenotype), ]

## add competitive p vals
df_pthcomp <- read.csv("../output/data/pathway_comp.csv")
# df_pathcomp_reord <- df_pthcomp %>% arrange(factor(Pathway, levels = unique(df_pathprs$Pathway)))
df_pathcomp_reord <- df_pthcomp[order(df_pthcomp$Pathway), ]
df_pathcomp_reord["GeseaPathwayCode"] <- df_pathprs_ordered$GeseaCode
df_pathcomp_reord["Groups"] <- c(
    "Microvasculature", "Microvasculature", "Microvasculature", "Microvasculature",
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Developmental", "Developmental", "Developmental", "Developmental",
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Microvasculature", "Microvasculature", "Microvasculature", "Microvasculature",
    "Developmental", "Developmental", "Developmental", "Developmental",
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Developmental", "Developmental", "Developmental", "Developmental"
)
df_pathcomp_reord <- subset(df_pathcomp_reord, select = c("Pathway", "GeseaPathwayCode", "Num_SNP", "Phenotype", "Estimate", "SE", "selfcontained.p", "competitive.p", "Groups"))
# correct for multiple comparisons
# df_pathcomp_reord <- df_pathcomp_reord %>%
#     group_by(Pathway) %>%
#     mutate(selfcon.pFWER = if_else(PC %in% c("PC1", "PC2"), p.adjust(selfcontained.p, method = "holm"), NA_real_)) %>%
#     mutate(comp.pFWER = if_else(PC %in% c("PC1", "PC2"), p.adjust(competitive.p, method = "holm"), NA_real_)) %>%
#     ungroup()
df_pathcomp_reord <- as.data.frame(df_pathcomp_reord)
df_pathcomp_reord["CIhigh"] <- round(df_pathcomp_reord$Estimate + (1.96 * df_pathcomp_reord$SE), digits = 3)
df_pathcomp_reord["CIlow"] <- round(df_pathcomp_reord$Estimate - (1.96 * df_pathcomp_reord$SE), digits = 3)
df_pathcomp_ord <- df_pathcomp_reord[order(df_pathcomp_reord$Pathway, df_pathcomp_reord$Phenotype), ]
df_pathcomp_ord$Pathway <- gsub("_", " ", df_pathcomp_ord$Pathway)
write.csv(df_pathcomp_ord, file = paste0("../output/data/pathwayreg_results.csv"), row.names = TRUE)


## Mediation analysis for GCIPL and INL-RPE and Neuroinflammatory pathway PRSSZ
df_filt <- read.csv("../output/data/df_mlm_inner_outer.csv")

prsvar <- "M24927" # Neuroinflamm pathway
# y_name <- c("GC_IPL_mean", "INL_RPE_mean")
y_name <- c("GC_IPL_mean")
# inflam_markers <- c("Monocytes", "Neutrophils", "SII", "NLR", "PLR", "LMR", "CRP")
inflam_markers <- c("SII", "NLR", "PLR", "LMR", "CRP")

cov <- "Age + Age_squared + Sex + Smoking_status + Alcohol_drinker_status + BMI + Townsend_index+ Hypertension + Diabetes_mellitus + Fasting_time+Image_quality_mean+Genotype_array + Genetic.PC1 + Genetic.PC2 + Genetic.PC3 + Genetic.PC4 + Genetic.PC5 + Genetic.PC6 + Genetic.PC7 + Genetic.PC8 + Genetic.PC9 + Genetic.PC10"
i <- 0
for (y in y_name) {
    results_list <- list()
    for (inflam_marker in inflam_markers) {
        model_total <- rlm(as.formula(paste(y, "~", prsvar, "+", cov)), data = df_filt, psi = psi.huber)
        model_direct <- rlm(as.formula(paste(y, "~", prsvar, "+", inflam_marker, "+", cov)), data = df_filt, psi = psi.huber)
        model_mediator <- rlm(as.formula(paste(inflam_marker, "~", prsvar, "+", cov)), data = df_filt, psi = psi.huber)
        model_outcome <- rlm(as.formula(paste(y, "~", inflam_marker, "+", cov)), data = df_filt, psi = psi.huber)

        a_coefficient <- coef(model_mediator)[prsvar]
        b_coefficient <- coef(model_outcome)[inflam_marker]
        indirect_effect <- a_coefficient * b_coefficient
        total_effect_coef <- coef(model_total)[prsvar]
        prop_mediated <- round((abs(indirect_effect) / abs(total_effect_coef)) * 100, 2)

        boot_indirect <- function(data, indices) {
            d <- data[indices, ]
            med_model <- rlm(as.formula(paste(inflam_marker, "~", prsvar, "+", cov)), data = d, psi = psi.huber)
            out_model <- rlm(as.formula(paste(y, "~", prsvar, "+", inflam_marker, "+", cov)), data = d, psi = psi.huber)
            a <- coef(med_model)[prsvar]
            b <- coef(out_model)[inflam_marker]
            return(a * b)
        }

        boot_results <- boot(data = df_filt, statistic = boot_indirect, R = 1000)
        boot_ci <- boot.ci(boot_results, type = "perc", conf = 0.95)
        medCI_low <- boot_ci$percent[4]
        medCI_high <- boot_ci$percent[5]
        boot_p_value <- boot.pval(boot_res = boot_results, type = "perc", theta_null = 0)

        add_asterisks <- function(p_value) {
            if (p_value < 0.001) {
                return("***")
            } else if (p_value < 0.01) {
                return("**")
            } else if (p_value < 0.05) {
                return("*")
            } else {
                return("")
            }
        }

        results_list[[inflam_marker]] <- list(
            total_effect = round(c(coef(model_total)[prsvar], confint.default(model_total, parm = prsvar)[1, ], f.robftest(model_total, var = prsvar)$p.value), digits = 2),
            direct_effect = round(c(coef(model_direct)[prsvar], confint.default(model_direct, parm = prsvar)[1, ], f.robftest(model_direct, var = prsvar)$p.value), digits = 2),
            path_a = round(c(coef(model_mediator)[prsvar], confint.default(model_mediator, parm = prsvar)[1, ], f.robftest(model_mediator, var = prsvar)$p.value), digits = 2),
            path_b = round(c(coef(model_outcome)[inflam_marker], confint.default(model_outcome, parm = inflam_marker)[1, ], f.robftest(model_outcome, var = inflam_marker)$p.value), digits = 2),
            indirect_effect = round(c(indirect_effect, boot_ci$percent[4:5], boot_p_value), digits = 3)
        )
    }

    format_result <- function(result, digits = 2) {
        estimate <- format(result[1], digits = digits)
        ci_lower <- format(result[2], digits = digits)
        ci_upper <- format(result[3], digits = digits)
        p_value <- result[4]
        asterisks <- add_asterisks(p_value)

        paste0(estimate, asterisks, " [", ci_lower, ", ", ci_upper, "]")
    }

    result_matrix <- matrix(nrow = 5, ncol = length(inflam_markers))
    rownames(result_matrix) <- c("Path A", "Path B", "Total Effect", "Direct Effect", "Indirect Effect")
    colnames(result_matrix) <- inflam_markers

    for (i in seq_along(inflam_markers)) {
        marker <- inflam_markers[i]
        result_matrix[, i] <- c(
            format_result(results_list[[marker]]$path_a),
            format_result(results_list[[marker]]$path_b),
            format_result(results_list[[marker]]$total_effect),
            format_result(results_list[[marker]]$direct_effect),
            format_result(results_list[[marker]]$indirect_effect)
        )
    }
    write.csv(result_matrix, file = paste0("../output/data/mediation_analysis_", y, ".csv"), row.names = TRUE)
    i <- i + 1
}

# Create mediation table as pdf
df_gcipl <- read.csv("../output/data/mediation_analysis_GC_IPL_mean.csv", row.names = NULL)
colnames(df_gcipl)[1] <- ""
med_table <- tableGrob(df_gcipl, rows = NULL)

# Save the table as a PDF
pdf("../output/figures/mediation_tab_gcipl.pdf", width = 15, height = 2)
grid.draw(med_table)
dev.off()


## Associations between individuals subfield's thickness and PRSSZ while controlling for confounds
df_filt <- read.csv("../output/data/df_mlm.csv")
df_ret <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c("Macular subfield", "hm", "Coef", "CI(low)", "CI(high)", "SE", "BPstat", "BPp", "F", "p"))
# ret_list <- c("GC_IPL_left","GC_IPL_right","Macula_left","Macula_right", "RNFL_left","RNFL_right")
ret_list <- c(
    "Inner_Inferior_left", "Inner_Inferior_right",
    "Outer_Inferior_left", "Outer_Inferior_right",
    "Inner_Nasal_left", "Inner_Nasal_right",
    "Outer_Nasal_left", "Outer_Nasal_right",
    "Inner_Superior_left", "Inner_Superior_right",
    "Outer_Superior_left", "Outer_Superior_right",
    "Inner_Temporal_left", "Inner_Temporal_right",
    "Outer_Temporal_left", "Outer_Temporal_right",
    "Central_left", "Central_right"
)
ret_list_mod <- gsub("_", " ", ret_list)
cov <- "Age + Age_squared + Sex + Smoking_status + Alcohol_drinker_status + BMI + Townsend_index+ Hypertension +Diabetes_mellitus+Genotype_array +Genetic.PC1+Genetic.PC2+ Genetic.PC3+ Genetic.PC4+ Genetic.PC5+Genetic.PC6+ Genetic.PC7+ Genetic.PC8+ Genetic.PC9+ Genetic.PC10"
# Loop over the new variable names
i <- 0
prsvar <- "PRSSZ"
pvals_pc <- c()
for (var_name in ret_list) {
    # Get eye label
    lat_lst <- as.list(strsplit(var_name, "_")[[1]])[-1]
    lat <- lat_lst[length(lat_lst)]
    octqc <- paste0("OCT_quality_", lat)

    # Create the formula string
    pc_str <- paste(var_name, "~", prsvar, "+", cov, "+", octqc)
    # Convert the string to a formula
    pc_obj <- as.formula(pc_str)
    # Run the regression
    m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)

    # pc_intercept<- round(map(m_pc$coefficients["(Intercept)"], 1)[[1]],digits = 2)
    pc_intercept <- round(m_pc$coefficients["(Intercept)"][[1]], digits = 2)
    pc1bp <- bptest(m_pc)
    tab1_pcbppval <- format_p(map(pc1bp$p.value, 1)[[1]], stars = FALSE)
    # tab1_pcbpstat <- round(map(pc1bp$statistic, 1)[[1]],digits = 2)
    tab1_pcbpstat <- round(pc1bp$statistic[[1]], digits = 2)

    pc_dd <- data.frame(summary(m_pc)$coefficients)
    # pc_b <- round(pc_dd$Value[2],digits = 2)
    pc_b <- pc_dd$Value[2]
    pc_se <- round(pc_dd$Std..Error[2], digits = 5)
    lgc_conflow <- round(confint.default(object = m_pc, parm = prsvar, level = 0.95)[1], digits = 3)
    lgc_confhigh <- round(confint.default(object = m_pc, parm = prsvar, level = 0.95)[2], digits = 3)

    fpc <- f.robftest(m_pc, var = prsvar)
    fpcFval <- round(fpc$statistic, digits = 2)
    fpcFval <- map(fpcFval, 1)[[1]]
    fpcpval <- round(fpc$p.value, digits = 5)

    pvals_pc <- append(pvals_pc, fpcpval)
    data <- unlist(list(var_name, lat, pc_b, lgc_conflow, lgc_confhigh, pc_se, tab1_pcbpstat, tab1_pcbppval, fpcFval, fpcpval), recursive = FALSE)
    df_it <- as.data.frame(t(data))
    new <- rep(i, ncol(df_it))
    df_ret[nrow(df_it) + i, ] <- df_it
    i <- i + 1
}
pc_fdr <- p.adjust(df_ret$p, method = "holm")
# pfdr_cols_ast <- paste(pc_fdr,p_ast,sep="")
df_ret$Coef <- as.numeric(df_ret$Coef)
df_ret$F <- as.numeric(df_ret$F)
df_ret$p <- as.numeric(df_ret$p)
df_ret["Macular subfield"] <- ret_list_mod
df_ret["pFWER"] <- pc_fdr # pfdr_cols_ast
df_ret["Groups"] <- c(
    "Inner Inferior", "Inner Inferior",
    "Outer Inferior", "Outer Inferior",
    "Inner Nasal", "Inner Nasal",
    "Outer Nasal", "Outer Nasal",
    "Inner Superior", "Inner Superior",
    "Outer Superior", "Outer Superior",
    "Inner Temporal", "Inner Temporal",
    "Outer Temporal", "Outer Temporal",
    "Central", "Central"
)

# test if F-stat is significantly different between left and right eye
lh_coef <- df_ret[df_ret$hm == "left", "Coef"]
rh_coef <- df_ret[df_ret$hm == "right", "Coef"]
normt_lh <- shapiro.test(lh_coef)
normt_rh <- shapiro.test(rh_coef)
paired_test_result <- t.test(lh_coef, rh_coef, paired = TRUE)
ttest_subf_stat <- round(paired_test_result$statistic[[1]], digits = 2)
ttest_subf_pval <- paired_test_result$p.value
ttest_subf_df <- paired_test_result$df

allfvals <- df_ret$Coef
color_range <- colorRampPalette(c("#065535", "#d3ffce"))(length(allfvals))[rank(allfvals)]
df_ret_col <- df_ret
df_ret_col["Fcolor"] <- color_range

subf_lbl <- c("Inner Superior", "Outer Superior", "Inner Temporal", "Outer Temporal", "Inner Inferior", "Outer Inferior", "Inner Nasal", "Outer Nasal")
# subf_num <- c('2','6','5','9','4','8','3','7')
subf_num <- c("IS", "OS", "IT", "OT", "II", "OI", "IN", "ON")
pie_data <- data.frame(
    subfield = subf_lbl,
    position = c(2020, 2021, 2020, 2021, 2020, 2021, 2020, 2021),
    thickness = c(10, 10, 10, 10, 10, 10, 10, 10)
)

# Run for each eye seperately
Fleft <- df_ret_col[df_ret_col$hm == "left", ]
Fright <- df_ret_col[df_ret_col$hm == "right", ]

df_subf_right <- Fright %>% arrange(factor(Groups, levels = subf_lbl))
df_subf_right$Coef <- round(df_subf_right$Coef, digits = 3)
# subfvals <- df_subf_sort$mean_F
color_code <- df_subf_right$Fcolor
# color_code = colorRampPalette(c('#d3ffce', '#065535'))(length(subfvals))[rank(subfvals)]
# color_code = c('black','black','grey','grey','white', 'white','white','white','white')


## Subfield associations results fundus plots (right eye)
subfheatmap_right <- pie_data %>%
    ggplot(aes(x = position, y = thickness, fill = subfield)) +
    geom_col(
        position = "fill", width = 1,
        color = "white", show.legend = FALSE
    ) +
    coord_polar(theta = "y", start = 225 * pi / 180) +
    lims(x = c(2019, 2022)) +
    scale_fill_manual(
        values = color_code[-length(color_code)],
        breaks = subf_lbl
    ) +
    geom_point(aes(x = 2019, y = 0),
        size = 49,
        shape = 21, fill = color_code[9], colour = "white"
    ) +
    geom_text(
        label = "CS", x = 2019, y = 0, size = 12, fontface = 2,
        show.legend = FALSE
    ) +
    geom_textpath(
        position = position_fill(vjust = .5),
        angle = 90, alpha = 1,
        aes(color = "black", label = subf_num), color = "black",
        size = 12, fontface = 2, show.legend = FALSE
    ) +
    theme_void()

ggsave("../output/figures/ETDRS_heatmap_right.png",
    subfheatmap_right,
    bg = "transparent",
    width = 3042, height = 3042, dpi = 300, units = "px"
)

img_inset <- image_read("../output/figures/ETDRS_heatmap_right.png")
img_inset <- image_scale(img_inset, "70%x")
img <- image_read("../data/retina_bg/macula_2D_right.png")

img_with_inset <- img %>% image_composite(
    img_inset,
    operator = "Atop",
    offset = "-300-200",
    gravity = "Center"
)
image_write(img_with_inset, "../output/figures/msubf_heatmap_right.png")


## Subfield associations results fundus plots (left eye)
df_subf_left <- Fleft %>% arrange(factor(Groups, levels = subf_lbl))
df_subf_left$Coef <- round(df_subf_left$Coef, digits = 3)
color_code <- df_subf_left$Fcolor
# color_code = colorRampPalette(c('#d3ffce', '#065535'))(length(subfvals))[rank(subfvals)]

subfheatmap_left <- pie_data %>% ggplot(aes(x = position, y = thickness, fill = subfield)) +
    geom_col(position = "fill", width = 1, color = "white", show.legend = FALSE) +
    coord_polar(theta = "y", start = 225 * pi / 180) +
    lims(x = c(2019, 2022)) +
    scale_fill_manual(values = color_code[-length(color_code)], breaks = subf_lbl) +
    geom_point(aes(x = 2019, y = 0), size = 49, shape = 21, fill = color_code[9], colour = "white") +
    geom_text(label = "CS", x = 2019, y = 0, size = 12, fontface = 2, show.legend = FALSE) +
    geom_textpath(
        position = position_fill(vjust = .5), angle = 90, alpha = 1,
        aes(color = "black", label = subf_num), color = "black", size = 12, fontface = 2, show.legend = FALSE
    ) +
    theme_void()

ggsave("../output/figures/ETDRS_heatmap_left.png", subfheatmap_left, bg = "transparent", width = 3042, height = 3042, dpi = 300, units = "px")

img_inset <- image_read("../output/figures/ETDRS_heatmap_left.png")
img_inset <- image_scale(img_inset, "70%x")
img <- image_read("../data/retina_bg/macula_2D_left.png")

img_with_inset <- img %>% image_composite(
    img_inset,
    operator = "Atop",
    offset = "+300-200",
    gravity = "Center"
)
image_write(img_with_inset, "../output/figures/msubf_heatmap_left.png")


## Diagram of exclusion/inclusion numbers
results <- read.csv("../output/data/retinflam_results_m.csv")
miss_vals <- read.csv("../output/data/df_missing_vals.csv")
# create participant diagramme
graph <- grViz("digraph {
graph [layout = dot, rankdir = TB]
 node [shape = rectangle, fontname = Helvetica]
rec1 [label = '@@1\n @@2']
rec2 [label = '@@3\n-@@4\n-@@5\n-@@6\n-@@7\n-@@8']
rec3 [label = '@@9\n @@10']
rec4 [label = '@@11\n @@12']
rec5 [label = '@@13\n @@14']
rec6 [label = '@@15\n-@@16\n-@@17\n @@18\n @@19\n @@20\n @@21']
rec7 [label = '@@22\n @@23']
 # edge definitions with the node IDs
rec1 -> rec2 -> rec3 -> rec4 -> rec6 -> rec7
rec3 -> rec5
}
[1]:  paste0('Individuals with macular measurements')
[2]:  paste0('(n = ', results$n_total[1], ')')
[3]:  paste0('Excluded (n = ', results$n_ethn[1]+results$n_qen_qc[1]+results$n_eyedis[1]+results$n_my_hyperopic[1]+results$n_octqc[1], ')')
[4]:  paste0('SNP QC/Sample QC (n = ', results$n_qen_qc[1], ')')
[5]:  paste0('Non British/Irish ancestry (n = ', results$n_ethn[1], ')')
[6]:  paste0('Eye diseases/disorders (n = ', results$n_eyedis[1], ')')
[7]:  paste0('Highly myopic/hyperopic eyes (n = ', results$n_my_hyperopic[1], ')')
[8]:  paste0('OCT image quality QC (n = ', results$n_octqc[1], ')')
[9]:  paste0('Individuals')
[10]: paste0('(n = ', results$n_total[1]-results$n_ethn[1]-results$n_qen_qc[1]-results$n_eyedis[1]-results$n_my_hyperopic[1]-results$n_octqc[1], ')')
[11]: paste0('Not Diagnosed')
[12]: paste0('(n = ', results$n_total[1]-results$n_ethn[1]-results$n_qen_qc[1]-results$n_eyedis[1]-results$n_my_hyperopic[1]-results$n_octqc[1]-results$n_diag[1], ')')
[13]: paste0('Diagnosed ICD-10, F20-29')
[14]: paste0('(n = ', results$n_diag[1], ')')
[15]: paste0('Excluded (n = ', results$n_nondiag_nan+results$n_antipsy[1], ')')
[16]: paste0('Using antipsychotics (n = ', results$n_antipsy[1], ')')
[17]: paste0('Incomplete data (n = ', results$n_nondiag_nan, '):')
[18]: paste0('(', miss_vals[[1]][1], ' n = ', miss_vals$Missing.values[1], ')')
[19]: paste0('(', miss_vals[[1]][2], ' n = ', miss_vals$Missing.values[2], ')')
[20]: paste0('(', miss_vals[[1]][3], ' n = ', miss_vals$Missing.values[3], ')')
[21]: paste0('(', miss_vals[[1]][4], ' n = ', miss_vals$Missing.values[4], ')')
[22]: paste0('Individuals')
[23]: paste0('(n = ', n_model_left, ')')
") # exchange line 16 with sample size used by macula thickness measurements
graph %>%
    export_svg() %>%
    charToRaw() %>%
    rsvg_pdf("../output/figures/Diagramme1.pdf")


## Define all variables for markdown
diag1 <- "../output/figures/Diagramme1.pdf"
mediat_tab <- "../output/figures/mediation_tab_gcipl.pdf"

fig1 <- "../data/retina_bg/Fig1_eye_anat.png"
fig2 <- "../output/figures/Fig1a_overall_retinap_partialreg.png"
fig2b <- "../output/figures/Fig1c_inner_outer_retinap_partialreg.png"
fig3a <- "../output/figures/msubf_heatmap_left.png"
fig3b <- "../output/figures/msubf_heatmap_right.png"

figs1a <- "../output/figures/Appx_FigS1_overall_retinap_partialreg.png"
figs2 <- "../output/figures/Appx_FigS2_inner_outer_regmatrix.png"
figs3 <- "../output/figures/Appx_FigureS1_CRP_logtransform.png"

p1c <- readPNG(fig3a)
p1d <- readPNG(fig3b)
p1cp <- rasterGrob(p1c, interpolate = TRUE)
p1dp <- rasterGrob(p1d, interpolate = TRUE)
figs1 <- cowplot::plot_grid(p1cp, p1dp, nrow = 1, label_size = 25)
figs1new <- ggsave(filename = "../output/figures/retinflam_figs1_subfields.png", figs1)

# collect sample stats and transform into table 1
columns <- c("Characteristic", "N", "Mean (SD)")
df_nsex <- data.frame(matrix(ncol = length(columns)))
colnames(df_nsex) <- columns
df_nsex[nrow(df_nsex) + 1, ] <- list("Female", results$female, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Male", results$male, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Diabetes Mellitus", results$n_hypertension, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Hypertension", results$n_diabetestwo, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Current smoker", results$curr_smoker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Previous smoker", results$prev_smoker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Non smoker", results$non_smoker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Current alcohol drinker", results$curr_alcohol_drinker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Previous alcohol drinker", results$prev_alcohol_drinker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Non alcohol drinker", results$non_alcohol_drinker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Age", n_model_left, results$mean_age)
df_nsex[nrow(df_nsex) + 1, ] <- list("BMI", n_model_left, results$mean_bmi)
df_nsex[nrow(df_nsex) + 1, ] <- list("Townsend_index", n_model_left, results$mean_townsend_index)
# append retinal phenotypes measures
print("here")
df_samp_retphen <- read_csv("../output/data/retphen_samp.csv",
    col_names = TRUE,
    show_col_types = FALSE
)
for (i in 1:2) {
    row <- df_samp_retphen[i, ]
    df_nsex[nrow(df_nsex) + 1, ] <- list(row$retinal_phenotype, n_model_left, row$thickness)
}
df_nsex <- df_nsex[-(1), ]
