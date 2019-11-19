library(tidyverse)
library(plyr)
library(doParallel)

cl = makeCluster(8, 'PSOCK')
registerDoParallel(cl)


setwd('/Users/alicia/daly_lab/GINGER/prs/ukbb/')

pheno_code <- read.delim('pheno_code_ukb_code.txt', sep='\t')

covariates = c("isFemale", "age", "age_squared", "age_isFemale", "age_squared_isFemale", paste0("PC", seq(20)))

# get population labels from here
date(); ukbb_pca <- read.table('/Users/alicia/daly_lab/ukbb_diverse_pops/pca/globalref_ukbb_pca_pops_rf_50.txt.gz', header=T) %>%
  select(s, pop); date() 
# get PCs and unrelated indicator from here
date(); pca_unrel <- read.table('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis/ukbb_covariates_all.txt.bgz', header=T) %>%
  filter(used_in_pca_calculation=='true') %>%
  left_join(ukbb_pca, by='s'); date()

afr_unrel <- subset(pca_unrel, pop=='AFR') %>%
  select(s)

afr_labels <- read_delim(gzfile('afrref_ukbb_pca_pops_rf_30_6PCs.txt.gz'), delim='\t') %>%
  right_join(afr_unrel)

table(afr_labels$pop)

# Bootstrap 95% CI around r^2
bootstrap_CIs <- function(my_subset) {
  sample_subset <- sample_n(my_subset, nrow(my_subset), replace=TRUE)
  model1 <- lm(as.formula(sprintf("pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=sample_subset)
  model0 <- lm(as.formula(sprintf("pheno ~ %s", paste(covariates, collapse = " + "))), data=sample_subset)
  sample_r2 = summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  return(sample_r2)
}

# Get R^2 and p-value for nested model
compute_r2 <- function(dataset, phenotype, p) {
  print(paste(phenotype, p))
  my_subset <- subset(dataset, p_threshold==p & pheno==phenotype)
  model1 <- lm(as.formula(sprintf("pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=my_subset)
  model0 <- lm(as.formula(sprintf("pheno ~ %s", paste(covariates, collapse = " + "))), data=my_subset)
  r2 <- summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  pval <- anova(model1, model0)$'Pr(>F'[2]
  date(); bs_means <- replicate(100, bootstrap_CIs(my_subset)); date()
  delta <- bs_means - r2
  d = quantile(delta, c(0.025, 0.975))
  ci = r2 - c(d[1], d[2])
  r2_vec <- c(r2, rev(ci), pval, my_subset$pheno[1], my_subset$p_threshold[1])
  names(r2_vec) <- c('r2', '2.5%', '97.5%', 'p', 'pheno', 'p_threshold')
  return(r2_vec)
}

# Read in and wrangle PRS + phenotype + covariate data
read_phenos <- function(pheno_code) {
  print(pheno_code)
  prs <- read_delim(gzfile(paste0('ukb31063.gwas_holdout.', pheno_code, '_PRS.txt.bgz')), delim='\t') %>%
    left_join(pca_unrel) %>%
    select(-sex)
  
  prs <- prs %>%
    gather('p_threshold', 'PRS', num_range('s', 1:10))
  prs$pheno <- pheno_code
  return(prs)
}

# "30750", 
pheno_codes <- c("4080", "4079", "21001", "21002", "50", "49", "48", "30620", "30600", "30610", "30650", "30840", "30690", "30730", "30760", "30780", "30870", "30000", "30010", "30020", "30030", "30040", "30050", "30060", "30070", "30080", "30100", "30200", "30180", "30190", "30210", "30220", "30120", "30140", "30130", "30150", "30160")
p_thresholds <- paste0('s', 1:10)

all_prs <- map_df(pheno_codes, read_phenos)

pheno_test <- c("4080", "4079")

date(); test <- compute_r2(all_prs, '4080', 's1'); date()

# need to split up by population
# run this for all phenos, takes ~20 minutes for all phenotypes
date(); r2_combined <- ldply(pheno_test, function(x) { ldply(p_thresholds, function(y) { compute_r2(all_prs, x, y) } ) }, 
                             .parallel=T, .paropts=list(.export = c("compute_r2", "bootstrap_CIs", "p_thresholds", "pheno_test", "covariates"),
                                                        .packages = c('dplyr', 'plyr'))); date()

