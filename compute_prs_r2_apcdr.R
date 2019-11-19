library(tidyverse)
library(cowplot)

setwd('/Users/alicia/daly_lab/GINGER/prs/apcdr/prs')

true_phenos <- read.delim('../gwas_phenotypes_28Oct14.txt', header=T, sep=' ') %>% 
  gather('pheno_name', 'true_pheno', SBP:BASO)
pheno_codes <- read.csv('../apcdr_ukb_pheno_info.csv', header=T)
pca_covar <- read.table('../apcdr2_pca.eigenvec', header=T)

true_phenos_covar <- true_phenos %>%
  right_join(pca_covar, by=c('GWAS_ID'='IID'))

p_thresholds <- c('p5e8', 'pe6', 'pe4', 'pe3', 'pe2', 'p05', 'p1', 'p2', 'p5', 'all')
#covariates = c("isFemale", "age", "age_squared", "age_isFemale", "age_squared_isFemale", paste0("PC", seq(20)))
covariates = c("sex_surv", "age_cal_surv", paste0("PC", seq(10)))

# Bootstrap 95% CI around r^2
bootstrap_CIs <- function(my_subset) {
  sample_subset <- sample_n(my_subset, nrow(my_subset), replace=TRUE)
  model1 <- lm(as.formula(sprintf("true_pheno ~ SCORE + %s", paste(covariates, collapse = " + "))), data=sample_subset)
  model0 <- lm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=sample_subset)
  sample_r2 = summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  return(sample_r2)
}

# Get R^2 and p-value for nested model
compute_r2 <- function(phenotype, p, dataset=prs_agg) {
  print(paste(phenotype, p))
  my_subset <- subset(dataset, p_threshold==p & prs_pheno==phenotype)
  model1 <- lm(as.formula(sprintf("true_pheno ~ SCORE + %s", paste(covariates, collapse = " + "))), data=my_subset)
  model0 <- lm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=my_subset)
  r2 <- summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  pval <- anova(model1, model0)$'Pr(>F'[2]
  date(); bs_means <- replicate(100, bootstrap_CIs(my_subset)); date()
  delta <- bs_means - r2
  d = quantile(delta, c(0.025, 0.975))
  ci = r2 - c(d[1], d[2])
  r2_vec <- c(r2, rev(ci), pval, my_subset$prs_pheno[1], my_subset$p_threshold[1])
  names(r2_vec) <- c('r2', '2.5%', '97.5%', 'p', 'pheno', 'p_threshold')
  return(data.frame(t(r2_vec)))
}

# read in and wrangle PRS + phenotype + covariate data
read_prs_merge_phenos <- function(pheno_name, p) {
  prs_pheno_code <- pheno_codes$pheno_code[which(pheno_codes$ukb_prefix == pheno_name)]
  #print(paste(pheno_name, pheno_code, p))
  prs <- read.table(paste0('rsid_', pheno_name, '.gwas.imputed_v3.both_sexes.', p, '.profile.gz'), header=T)
  prs$prs_pheno <- prs_pheno_code
  prs$p_threshold <- p
  prs_true <- prs %>%
    left_join(true_phenos_covar, by=c('IID'='GWAS_ID', 'prs_pheno'='pheno_name'))
  prs_true$true_pheno <- as.numeric(prs_true$true_pheno)
  return(prs_true)
}

argList <- list(x = pheno_codes$ukb_prefix, y = p_thresholds)
crossArg <- cross_df(argList)
prs_agg <- map2_dfr(crossArg$x, crossArg$y, read_prs_merge_phenos)

argList <- list(x = pheno_codes$pheno_code, y = p_thresholds)
crossArg <- cross_df(argList)
#argList <- list(prs_agg, crossArg$x, crossArg$y)
#r2_combined <- pmap_dfr(argList, compute_r2)
date(); r2_combined <- map2_dfr(crossArg$x, crossArg$y, compute_r2); date()
# Error: Argument 1 must have names
# error because return type not right
# this takes [1] "Tue Oct 22 22:21:46 2019"
# to 

write.table(r2_combined, 'apcdr_from_ukbb.r2.txt', quote = F, row.names = F, sep='\t')

max_r2 <- r2_combined %>%
  group_by(pheno) %>%
  arrange(desc(r2)) %>%
  slice(1)

max_r2$r2 <- as.numeric(max_r2$r2)
max_r2$CI95_low <- as.numeric(max_r2$X2.5.)
max_r2$CI95_high <- as.numeric(max_r2$X97.5.)
max_r2 <- max_r2 %>% arrange(desc(r2))
max_r2$pheno <- factor(max_r2$pheno, levels=max_r2$pheno)

pd <- position_dodge(0.5)
p1 <- ggplot(max_r2, aes(x=pheno, y=r2)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=0) +
  labs(x='Phenotype', y=bquote('Max'~R^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=10),
        text = element_text(size=12, color='black'),
        axis.text.y = element_text(size=10))

save_plot('apcdr_from_ukbb_r2_top.pdf', p1, base_height=6, base_width=10)

r2_combined$rand <- rnorm(nrow(r2_combined), 0, 1)
r2_combined[sample(1:nrow(r2_combined), round(nrow(r2_combined) * 0.99, 0)), 'rand'] <- 0
