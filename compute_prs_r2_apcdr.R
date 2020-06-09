library(tidyverse)
library(cowplot)
library(tictoc)

setwd('/Users/alicia/daly_lab/GINGER/prs/apcdr/prs')

# Bootstrap 95% CI around r^2
bootstrap_CIs <- function(my_subset) {
  sample_subset <- sample_n(my_subset, nrow(my_subset), replace=TRUE)
  model1 <- lm(as.formula(sprintf("pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=sample_subset)
  model0 <- lm(as.formula(sprintf("pheno ~ %s", paste(covariates, collapse = " + "))), data=sample_subset)
  sample_r2 = summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  return(sample_r2)
}

# Get R^2 and p-value for nested model
compute_r2 <- function(phenotype, p, dataset=prs_agg) {
  print(paste(phenotype, p)) # check this
  my_subset <- subset(dataset, p_threshold==p & pheno_code==phenotype)
  model1 <- lm(as.formula(sprintf("pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=my_subset)
  model0 <- lm(as.formula(sprintf("pheno ~ %s", paste(covariates, collapse = " + "))), data=my_subset)
  r2 <- summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  pval <- anova(model1, model0)$'Pr(>F'[2]
  date(); bs_means <- replicate(100, bootstrap_CIs(my_subset)); date()
  delta <- bs_means - r2
  d = quantile(delta, c(0.025, 0.975))
  ci = r2 - c(d[1], d[2])
  r2_vec <- c(r2, rev(ci), pval, my_subset$pheno_code[1], my_subset$p_threshold[1], my_subset$analysis[1])
  names(r2_vec) <- c('r2', 'CI_2.5', 'CI_97.5%', 'p', 'pheno', 'p_threshold', 'analysis')
  return(data.frame(t(r2_vec)))
}

read_phenos <- function(pheno_code, analysis, dirname) {
  print(pheno_code)
  prs_file <- paste0(dirname, pheno_code, '_apcdr_PRS.txt.bgz')
  if (file.exists(prs_file)) {
    prs <- read_delim(gzfile(prs_file), delim='\t', col_types=cols(pheno=col_double())) %>%
      inner_join(pca_covar, by=c('s'='IID'))
    
    prs <- prs %>%
      gather('p_threshold', 'PRS', num_range('s', 1:10)) %>%
      mutate(analysis=analysis)
    prs$pheno_code <- as.character(pheno_code)
    tbl_df(prs)
    return(prs)
  }
}

# true_phenos <- read.delim('../gwas_phenotypes_28Oct14.txt', header=T, sep=' ') %>% 
#   gather('pheno_name', 'true_pheno', SBP:BASO)
pheno_codes <- read.csv('../apcdr_ukb_pheno_info.csv', header=T)
pca_covar <- read.table('../apcdr2_pca.eigenvec', header=T) %>% select(-FID)

ukb_prs <- map_df(pheno_codes$pheno_code, function(x) read_phenos(x, analysis='UKB', dirname = 'ukb_10k_eur_holdout/ukb31063.gwas_holdout.'))
bbj_ukb_prs <- map_df(pheno_codes$pheno_code, function(x) read_phenos(x, analysis='BBJ+UKB', dirname = 'bbj_ukb_10k_eur_holdout_meta/'))
bbj_page_ukb_prs <- map_df(pheno_codes$pheno_code, function(x) read_phenos(x, analysis='BBJ+PAGE+UKB', dirname = 'bbj_page_ukb_10k_eur_holdout_meta/'))

p_thresholds <- paste0('s', 1:10)
#covariates = c("isFemale", "age", "age_squared", "age_isFemale", "age_squared_isFemale", paste0("PC", seq(20)))
covariates = c("age", "sex", paste0("PC", seq(10)))

# argList <- list(x = pheno_codes$ukb_prefix, y = p_thresholds)
# crossArg <- cross_df(argList)
# prs_agg <- map2_dfr(crossArg$x, crossArg$y, read_prs_merge_phenos)

argList <- list(x = pheno_codes$pheno_code, y = p_thresholds)
crossArg <- cross_df(argList)

run_r2 <- function(crossArg, dataset=ukb_prs) {
  argList <- list(x = unique(dataset$pheno_code), y = unique(dataset$p_threshold))
  crossArg <- cross_df(argList)
  tic(); r2_combined <- pmap(list(crossArg$x, crossArg$y), function(x, y) { compute_r2(x, y, dataset) } ); toc()
  r2_combined2 <- do.call(rbind.data.frame, r2_combined) %>%
    mutate_at(vars(r2:p), function(x) as.numeric(as.character(x)))
  r2_combined2$r2 <- as.numeric(as.character(r2_combined2$r2))
  return(r2_combined2)
}

r2_comb_ukb <- run_r2(crossArg, dataset=ukb_prs)
write.table(r2_comb_ukb, 'apcdr.ukb_prs.r2.txt', quote = F, row.names = F, sep='\t'); 
r2_comb_bbj_ukb <- run_r2(crossArg, dataset=bbj_ukb_prs)
write.table(r2_comb_bbj_ukb, 'apcdr.bbj_ukb_prs.r2.txt', quote = F, row.names = F, sep='\t'); 
r2_comb_bbj_page_ukb <- run_r2(crossArg, dataset=bbj_page_ukb_prs)
write.table(r2_comb_bbj_page_ukb, 'apcdr.bbj_page_ukb_prs.r2.txt', quote = F, row.names = F, sep='\t'); 
r2_combined <- bind_rows(r2_comb_ukb, r2_comb_bbj_ukb, r2_comb_bbj_page_ukb)

write.table(r2_combined, 'apcdr_from_ukbb.r2.txt', quote = F, row.names = F, sep='\t')


# Plot R^2 results --------------------------------------------------------

r2_combined <- read.table('apcdr_from_ukbb.r2.txt', header=T)
phenos_in_all <- as.character(unique(subset(r2_combined, analysis=='BBJ+PAGE+UKB')$pheno))
max_r2 <- r2_combined %>%
  group_by(pheno, analysis) %>%
  arrange(desc(r2)) %>%
  filter(pheno %in% phenos_in_all) %>%
  slice(1)

max_r2$r2 <- as.numeric(max_r2$r2)
max_r2$CI95_low <- as.numeric(max_r2$CI_2.5)
max_r2$CI95_high <- as.numeric(max_r2$CI_97.5)
max_r2 <- max_r2 %>% arrange(desc(r2))
max_r2$pheno <- factor(max_r2$pheno, levels=max_r2$pheno)

pd <- position_dodge(0.5)
p1 <- ggplot(max_r2, aes(x=pheno, y=r2, color=analysis)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=0, position=pd) +
  labs(x='Phenotype', y=bquote('Max'~R^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=10),
        text = element_text(size=12, color='black'),
        axis.text.y = element_text(size=10))

save_plot('apcdr_from_ukbb_r2_top.pdf', p1, base_height=6, base_width=10)

r2_combined$rand <- rnorm(nrow(r2_combined), 0, 1)
r2_combined[sample(1:nrow(r2_combined), round(nrow(r2_combined) * 0.99, 0)), 'rand'] <- 0


# Split some Uganda PRS by phenos -----------------------------------------

# low accuracy in WC PRS specific to APCDR due to malnutrition?
apcdr_phenos <- read.delim('../gwas_phenotypes_28Oct14.txt', header=T, sep=' ') %>%
  mutate(weight_range=case_when(BMI < 18.5 ~ 'Underweight\n(BMI<18.5)',
                                BMI >= 18.5 & BMI < 25 ~ 'Normal weight\n(18.5<=BMI<25)',
                                BMI >= 25 & BMI < 30 ~ 'Overweight\n(25<=BMI<30)',
                                BMI >= 30 ~ 'Obese\n(BMI>=30)'))

p_weight <- ggplot(subset(apcdr_phenos, weight_range!='NA'), aes(x=as.numeric(as.character(WC)), fill=weight_range)) +
  geom_density(alpha=0.8) +
  scale_fill_brewer(palette='Set1') +
  labs(x='Waist circumference', y='Density', fill='Weight range') +
  theme_classic() +
  theme(axis.text = element_text(size = 14, color='black'),
        text = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = c(0.8, 0.8))

ggsave('apcdr_bmi_wc.pdf', height=5, width=5)  


apcdr_phenos <- apcdr_phenos %>%
  select(GWAS_ID, weight_range)
prs_agg_wc <- subset(prs_agg, prs_pheno=='WC') %>%
  left_join(apcdr_phenos, by=c('IID'='GWAS_ID')) %>%
  mutate(prs_pheno=weight_range)
prs_agg_wc$sex_surv <- as.numeric(prs_agg_wc$sex_surv)
prs_agg_wc$weight_range <- factor(prs_agg_wc$weight_range, levels=c('Underweight\n(BMI<18.5)',
                                                                    'Normal weight\n(18.5<=BMI<25)',
                                                                    'Overweight\n(25<=BMI<30)',
                                                                    'Obese\n(BMI>=30)'))



date(); r2_comb_wc <- map_dfr(levels(prs_agg_wc$weight_range), function(x) { 
  map_dfr(p_thresholds, function(y) { 
    compute_r2(x, y, prs_agg_wc)}) } ); date()

max_r2_wc <- r2_comb_wc %>%
  group_by(pheno) %>%
  arrange(desc(r2)) %>%
  slice(1)
max_r2_wc$r2 <- as.numeric(max_r2_wc$r2)
max_r2_wc$CI95_low <- as.numeric(max_r2_wc$X2.5.)
max_r2_wc$CI95_high <- as.numeric(max_r2_wc$X97.5.)
max_r2_wc <- max_r2_wc %>% arrange(desc(r2))
max_r2_wc$pheno <- factor(max_r2_wc$pheno, levels=max_r2_wc$pheno)

p_wc <- ggplot(max_r2_wc, aes(x=pheno, y=r2)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=0) +
  labs(x='Phenotype', y=bquote('Max'~R^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=10),
        text = element_text(size=12, color='black'),
        axis.text.y = element_text(size=10))
