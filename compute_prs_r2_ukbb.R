library(tidyverse)
library(furrr)
library(tictoc)
library(cowplot)
# library(plyr)
# library(doParallel)
# 
# cl = makeCluster(8, 'PSOCK')
# registerDoParallel(cl)

setwd('/Users/alicia/daly_lab/GINGER/prs/ukbb/')

pheno_code <- read.delim('pheno_code_ukb_code.txt', sep='\t')

covariates = c("isFemale", "age", "age_squared", "age_isFemale", "age_squared_isFemale", paste0("PC", seq(10)))

# get population labels from here
date(); ukbb_pca <- read.table('/Users/alicia/daly_lab/ukbb_diverse_pops/pca/globalref_ukbb_pca_pops_rf_50.txt.gz', header=T) %>%
  select(s, pop); date() 
# get GWAS holdout inds
eur_gwas <- read.table('ukb31063.gwas_samples.gwas_vs_holdout.txt', header=T)
# get PCs and unrelated indicator from here
date(); pca_unrel <- read.table('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis/ukbb_covariates_all.txt.bgz', header=T) %>%
  filter(used_in_pca_calculation=='true') %>%
  left_join(ukbb_pca, by='s') %>%
  left_join(eur_gwas, by='s')%>%
  subset(pop=='EUR'&in_gwas==FALSE | pop!='EUR'); date()

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
compute_r2 <- function(phenotype, population, p, dataset=all_prs) { 
  print(paste(phenotype, p, population, date()))
  my_subset <- subset(dataset, p_threshold==p & pheno_code==phenotype & pop==population)
  model1 <- lm(as.formula(sprintf("pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=my_subset)
  model0 <- lm(as.formula(sprintf("pheno ~ %s", paste(covariates, collapse = " + "))), data=my_subset)
  r2 <- summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  pval <- anova(model1, model0)$'Pr(>F'[2]
  date(); bs_means <- replicate(100, bootstrap_CIs(my_subset)); date()
  delta <- bs_means - r2
  d = quantile(delta, c(0.025, 0.975), na.rm=T)
  ci = r2 - c(d[1], d[2])
  r2_vec <- c(r2, rev(ci), pval, phenotype, p, population)
  names(r2_vec) <- c('r2', 'CI_2.5', 'CI_97.5', 'p', 'pheno', 'p_threshold', 'pop')
  return(data.frame(t(r2_vec)))
}

# Read in and wrangle PRS + phenotype + covariate data
read_phenos <- function(pheno_code, AFR=FALSE) {
  print(pheno_code)
  prs <- read_delim(gzfile(paste0('ukb31063.gwas_holdout.', pheno_code, '_PRS.txt.bgz')), delim='\t') %>%
    inner_join(pca_unrel) %>%
    select(-sex)
  
  if(AFR) {
    prs <- prs %>% 
      select(-(PC1:PC20), -pop) %>%
      inner_join(afr_labels)
  }
  
  prs <- prs %>%
    gather('p_threshold', 'PRS', num_range('s', 1:10))
  prs$pheno_code <- pheno_code
  return(prs)
}

pheno_codes <- c("4080", "4079", "21001", "21002", "50", "49", "48", "30620", "30600", "30610", "30650", "30840", "30690", "30730", "30760", "30780", "30870", "30000", "30010", "30020", "30030", "30040", "30050", "30060", "30070", "30080", "30100", "30200", "30180", "30190", "30210", "30220", "30120", "30140", "30130", "30150", "30160")
p_thresholds <- paste0('s', 1:10)

all_prs <- map_df(pheno_codes, read_phenos)
afr_prs <- map_df(pheno_codes, function(x) read_phenos(x, TRUE))

pops <- as.character(unique(all_prs$pop))
pops <- pops[!pops %in% c('OCE', 'oth')]
afr_pops <- as.character(unique(afr_prs$pop))

test <- compute_r2('50', 'EUR', 's5')

run_r2 <- function(which_pop, phenotypes, dataset=all_prs) {
  argList <- list(x = phenotypes, y = which_pop, z = p_thresholds)
  crossArg <- cross_df(argList)
  tic(); r2_combined <- pmap(list(crossArg$x, crossArg$y, crossArg$z), function(x, y, z) { compute_r2(x, y, z, dataset) } ); toc()
  r2_combined2 <- do.call(rbind.data.frame, r2_combined) %>%
    mutate_at(vars(r2:p), function(x) as.numeric(as.character(x)))
  r2_combined2$r2 <- as.numeric(as.character(r2_combined2$r2))
  return(r2_combined2)
}

r2_combined <- run_r2(pops, pheno_codes)
r2_combined_afr <- run_r2(afr_pops, pheno_codes, afr_prs)

write.table(r2_combined, 'ukbb_holdout.cont.r2.txt', quote = F, row.names = F, sep='\t')
write.table(r2_combined_afr, 'ukbb_holdout.afr.r2.txt', quote = F, row.names = F, sep='\t')


# read R2 values and plot -------------------------------------------------

r2_combined <- read.table('ukbb_holdout.cont.r2.txt', header=T)
r2_combined_afr <- read.table('ukbb_holdout.afr.r2.txt', header=T)
r2_combined_apcdr <- read.table('../apcdr/prs/apcdr_from_ukbb.r2.txt', header=T) %>%
  mutate(pop='Uganda', analysis='APCDR Uganda', pheno_code=pheno) %>%
  group_by(pheno, pop) %>%
  arrange(desc(r2)) %>%
  slice(1) %>%
  ungroup() %>%
  select(-pheno)
colnames(r2_combined_apcdr)[2:3] <- c('CI_2.5', 'CI_97.5')

summarize_r2 <- function(r2) {
  max_r2 <- r2 %>%
    group_by(pheno, pop) %>%
    arrange(desc(r2)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate_at(vars(pheno), function(x) as.numeric(as.character(x)))
  max_r2 <- max_r2 %>% 
    left_join(pheno_code, by=c('pheno'='ukb_code'))
  pop_order <- max_r2 %>% group_by(pop) %>% summarize(r2=mean(r2)) %>% arrange(desc(r2))
  max_r2$pop <- factor(max_r2$pop, levels=as.character(pop_order$pop))
  pheno_order <- max_r2 %>% arrange(desc(r2)) %>% ungroup() %>% select(pheno_code) %>% unique()
  max_r2$pheno_code <- factor(max_r2$pheno_code, levels=as.character(pheno_order$pheno_code))
  return(max_r2)
}

max_r2 <- summarize_r2(r2_combined) %>% mutate(analysis='UKB continental ancestry')
max_r2_afr <- summarize_r2(r2_combined_afr) %>% mutate(analysis='UKB African ancestry')
max_r2_comb <- max_r2 %>% bind_rows(max_r2_afr) %>% bind_rows(r2_combined_apcdr)
max_r2_comb$analysis <- factor(max_r2_comb$analysis, levels=c('UKB continental ancestry', 'UKB African ancestry', 'APCDR Uganda'))

r2_summarize <- max_r2_comb %>%
  ungroup() %>%
  group_by(pheno_code) %>%
  dplyr::mutate(best_eur_pred=max(r2 * (pop=='EUR'))) %>%
  group_by(pop, pheno_code) %>%
  dplyr::mutate(rel_eur = r2 / best_eur_pred) %>%
  subset(CI_97.5-CI_2.5 < 0.1)
mean_r2 <- r2_summarize %>%
  group_by(pop) %>%
  summarise(mean_r2=mean(r2), mean_rel_eur = mean(rel_eur), sem_rel_eur=sd(rel_eur)/sqrt(n()),
            median_r2=median(r2), median_rel_eur=median(rel_eur), mad_rel_eur=mad(rel_eur)/sqrt(n())) %>%
  arrange(desc(mean_rel_eur))
r2_summarize <- r2_summarize %>%
  left_join(mean_r2, by='pop')

r2_summarize$pop <- factor(r2_summarize$pop, levels = mean_r2$pop)

style_sheet <- read.csv('/Users/alicia/daly_lab/grants/neurogap/U01 ancestral pops/pca/tgp_color_style_sheet.csv', header=T)
color_vec <- as.character(style_sheet$Color)
names(color_vec) <- style_sheet$Population

p1 <- ggplot(r2_summarize, aes(x=pop, y=rel_eur, fill=pop)) +
  facet_wrap(~analysis, scales='free_x', shrink=T) +
  geom_violin() +
  scale_fill_manual(values=color_vec) +
  geom_point(aes(x=pop, y=median_rel_eur), color='black') +
  geom_errorbar(aes(ymin=median_rel_eur - mad_rel_eur, ymax=median_rel_eur + mad_rel_eur),color='black', width=0.1) +
  labs(x='Population', y='Prediction accuracy\n(Relative to Europeans)') +
  theme_bw() +
  guides(fill=F) +
  theme(axis.text = element_text(color='black', angle=45, hjust=1),
        text = element_text(size=16))

  #r2_summarize %>% ungroup %>% count(pop, analysis) %>% arrange((analysis))

#ggsave('ukb_continents_within_africa_rel_r2.pdf', p1, width=8, height=6)
ggsave('ukb_continent_africa_apcdr_rel_r2.pdf', p1, width=10, height=6)



# Compute prediction accuracy from meta-analysis PRS ----------------------

# Read in and wrangle PRS + phenotype + covariate data
read_phenos_meta <- function(pheno_code, AFR=FALSE) {
  print(pheno_code)
  prs <- read_delim(gzfile(paste0('apcdr_ukb_10k_eur_holdout_meta/', pheno_code, '_PRS.txt.bgz')), delim='\t') %>%
    inner_join(pca_unrel) %>%
    select(-sex)
  
  if(AFR) {
    prs <- prs %>% 
      select(-(PC1:PC20), -pop) %>%
      inner_join(afr_labels)
  }
  
  prs <- prs %>%
    gather('p_threshold', 'PRS', num_range('s', 1:10))
  prs$pheno_code <- pheno_code
  return(prs)
}

phenos <- c('ALP', 'ALT', 'AST', 'Albumin', 'BASO', 'BMI', 'Bilirubin', 'Cholesterol', 'DBP', 'EOS', 'GGT', 'HCT', 'HC', 'HDL', 'HT', 'Hb', 'LDL', 'LYMPH', 'MCHC', 'MCH', 'MCV', 'MONO', 'NEU', 'PLT', 'RBC', 'RDW', 'SBP', 'TG', 'WBC', 'WC', 'WT')
all_prs_meta <- map_df(phenos, read_phenos_meta)
all_prs_meta$pop <- as.character(all_prs_meta$pop)

r2_meta <- run_r2(pops, phenos, dataset=all_prs_meta)
write.table(r2_meta, 'apcdr_ukb_10k_eur_holdout_meta/ukbb_holdout.ukb_ugr.cont.r2.txt', quote = F, row.names = F, sep='\t')

# read R2 values from meta-analysis and plot -------------------------------------------------

r2_meta <- read.table('apcdr_ukb_10k_eur_holdout_meta/ukbb_holdout.ukb_ugr.cont.r2.txt', header=T)
# r2_combined_afr <- read.table('ukbb_holdout.afr.r2.txt', header=T)
# r2_combined_apcdr <- read.table('../apcdr/prs/apcdr_from_ukbb.r2.txt', header=T) %>%
#   mutate(pop='Uganda', analysis='APCDR Uganda', pheno_code=pheno) %>%
#   group_by(pheno, pop) %>%
#   arrange(desc(r2)) %>%
#   slice(1) %>%
#   ungroup() %>%
#   select(-pheno)
# colnames(r2_combined_apcdr)[2:3] <- c('CI_2.5', 'CI_97.5')

summarize_r2_meta <- function(r2) {
  max_r2 <- r2 %>%
    group_by(pheno, pop) %>%
    arrange(desc(r2)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(pheno_code = pheno)
  print(head(max_r2))
  pop_order <- max_r2 %>% group_by(pop) %>% summarize(r2=mean(r2)) %>% arrange(desc(r2))
  max_r2$pop <- factor(max_r2$pop, levels=as.character(pop_order$pop))
  pheno_order <- max_r2 %>% arrange(desc(r2)) %>% ungroup() %>% select(pheno_code) %>% unique()
  max_r2$pheno_code <- factor(max_r2$pheno_code, levels=as.character(pheno_order$pheno_code))
  return(max_r2)
}

max_r2_meta <- summarize_r2_meta(r2_meta) %>% mutate(analysis='UKB+UGR meta')
# #max_r2_afr <- summarize_r2(r2_combined_afr) %>% mutate(analysis='UKB African ancestry')
# max_r2_comb <- max_r2 %>% bind_rows(max_r2_afr) %>% bind_rows(r2_combined_apcdr)
# max_r2_comb$analysis <- factor(max_r2_comb$analysis, levels=c('UKB continental ancestry', 'UKB African ancestry', 'APCDR Uganda'))

r2_summarize_meta <- max_r2_meta %>%
  ungroup() %>%
  group_by(pheno_code) %>%
  dplyr::mutate(best_eur_pred=max(r2 * (pop=='EUR'))) %>%
  group_by(pop, pheno_code) %>%
  dplyr::mutate(rel_eur = r2 / best_eur_pred) %>%
  subset(CI_97.5-CI_2.5 < 0.1)
mean_r2_meta <- r2_summarize_meta %>%
  group_by(pop) %>%
  summarise(mean_r2=mean(r2), mean_rel_eur = mean(rel_eur), sem_rel_eur=sd(rel_eur)/sqrt(n()),
            median_r2=median(r2), median_rel_eur=median(rel_eur), mad_rel_eur=mad(rel_eur)/sqrt(n())) %>%
  arrange(desc(mean_rel_eur))
r2_summarize_meta <- r2_summarize_meta %>%
  left_join(mean_r2, by='pop')

r2_summarize_meta$pop <- factor(r2_summarize_meta$pop, levels = mean_r2$pop)

style_sheet <- read.csv('/Users/alicia/daly_lab/grants/neurogap/U01 ancestral pops/pca/tgp_color_style_sheet.csv', header=T)
color_vec <- as.character(style_sheet$Color)
names(color_vec) <- style_sheet$Population

p2 <- ggplot(r2_summarize_meta, aes(x=pop, y=rel_eur, fill=pop)) +
  #facet_wrap(~analysis, scales='free_x', shrink=T) +
  geom_violin() +
  scale_fill_manual(values=color_vec) +
  geom_point(aes(x=pop, y=median_rel_eur), color='black') +
  geom_errorbar(aes(ymin=median_rel_eur - mad_rel_eur, ymax=median_rel_eur + mad_rel_eur),color='black', width=0.1) +
  labs(x='Population', y='Prediction accuracy\n(Relative to Europeans)') +
  theme_bw() +
  guides(fill=F) +
  theme(axis.text = element_text(color='black', angle=45, hjust=1),
        text = element_text(size=16))

ggsave('apcdr_ukb_10k_eur_holdout_meta/ukb_meta_rel_r2.pdf', p1, width=10, height=6)

r2_summarize_comb <- bind_rows(r2_summarize %>% select(-pheno), r2_summarize_meta %>% select(-pheno))
r2_summarize_comb2 <- subset(r2_summarize_comb, analysis %in% c('UKB+UGR meta', 'UKB continental ancestry')) %>%
  mutate(analysis_order = case_when(analysis == 'UKB continental ancestry' & pop == 'EUR' ~ 1,
                                    analysis == 'UKB+UGR meta' & pop == 'EUR' ~ 2,
                                    analysis == 'UKB continental ancestry' & pop == 'AMR' ~ 5,
                                    analysis == 'UKB+UGR meta' & pop == 'AMR' ~ 6,
                                    analysis == 'UKB continental ancestry' & pop == 'MID' ~ 9,
                                    analysis == 'UKB+UGR meta' & pop == 'MID' ~ 10,
                                    analysis == 'UKB continental ancestry' & pop == 'CSA' ~ 13,
                                    analysis == 'UKB+UGR meta' & pop == 'CSA' ~ 14,
                                    analysis == 'UKB continental ancestry' & pop == 'EAS' ~ 17,
                                    analysis == 'UKB+UGR meta' & pop == 'EAS' ~ 18,
                                    analysis == 'UKB continental ancestry' & pop == 'AFR' ~ 21,
                                    analysis == 'UKB+UGR meta' & pop == 'AFR' ~ 22
         ),
         analysis_order_jitter = jitter(analysis_order, amount=0.09))
r2_summarize_comb2 <- r2_summarize_comb2 %>%
  mutate(ids = paste(pop, pheno_code, sep='_'))

xbreaks = c(seq(1,21,by=4), seq(2,22, by=4))
xbreaks = xbreaks[order(xbreaks)]
library(gghalves)

p3 <- ggplot(r2_summarize_comb2, aes(y=rel_eur)) +
  geom_point(aes(x=analysis_order_jitter, color=pop), size=1.5, alpha=0.6) +
  geom_line(aes(x=analysis_order_jitter, id=ids), color='lightgrey') +
  
  geom_half_violin(data = r2_summarize_comb2 %>% filter(analysis_order %% 2 == 1), 
                   aes(x=analysis_order, y=rel_eur, group=analysis_order, fill=pop), 
                   position = position_nudge(x=-0.3),
                   draw_quantiles = c(0.25,.5,.75)) +
  geom_half_violin(data = r2_summarize_comb2 %>% filter(analysis_order %% 2 == 0), 
                   aes(x=analysis_order, y=rel_eur, group=analysis_order, fill=pop), 
                   position = position_nudge(x=0.3),
                   draw_quantiles = c(0.25,.5,.75), side='r') +
  # Define additional plot settings
  theme_classic() +
  scale_fill_manual(values=color_vec, name='Population') +
  scale_color_manual(values=color_vec, name='Population') +
  scale_x_continuous(breaks=xbreaks, labels=rep(c('UKB EUR', 'UKB EUR+UGR'), 6)) +
  xlab('Discovery cohort') + ylab('Prediction accuracy\n(Relative to Europeans)') +
  theme(axis.text.x=element_text(angle=45, vjust = 1, hjust = 1),
        axis.text = element_text(color='black', size=10))
  
ggsave('apcdr_ukb_10k_eur_holdout_meta/ukb_only_ugr_meta_rel_r2.pdf', p3, width=10, height=6)


#r2_summarize %>% ungroup %>% count(pop, analysis) %>% arrange((analysis))

#ggsave('ukb_continents_within_africa_rel_r2.pdf', p1, width=8, height=6)



# Plot all traits for all ancestries and AFR ancestries (max R2) -----------------------------

pd <- position_dodge(0.4) # move them .05 to the left and right
ylims = c(min(max_r2_afr$CI_2.5, max_r2$CI_2.5), max(max_r2_afr$CI_97.5, max_r2$CI_97.5))
p1 <- ggplot(max_r2, aes(x=pheno_code, y=r2, color=pop)) +
  geom_point(position=pd) +
  #facet_grid(~sumstats) +
  geom_errorbar(aes(ymin=CI_2.5, ymax=CI_97.5), width=.2, position=pd) +
  scale_color_brewer(palette='Dark2', name='Populations') +
  labs(x='Phenotype', y=bquote(Maximum~R^2), title='Continental ancestries') +
  ylim(ylims) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
        text = element_text(size=16, color='black'),
        legend.justification = c(1, 1), legend.position = c(1, 1))

p2 <- ggplot(max_r2_afr, aes(x=pheno_code, y=r2, color=pop)) + #, shape=pop
  geom_point(position=pd) +
  #facet_grid(~sumstats) +
  geom_errorbar(aes(ymin=CI_2.5, ymax=CI_97.5), width=.5, position=pd) +
  scale_color_brewer(palette='Dark2', name='Populations') +
  #scale_shape_manual(name='Population') +
  labs(x='Phenotype', y=bquote(Maximum~R^2), title='African ancestries') +
  ylim(ylims) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
        text = element_text(size=16, color='black'),
        legend.justification = c(1, 1), legend.position = c(1, 1))

p <- plot_grid(p1, p2, labels = c('a', 'b'))
ggsave('ukbb_diverse_pops_afr_r2.pdf', p, width=12, height=6)


# Compare East Africans across cohorts ------------------------------------

apcdr_ukb_east_r2 <- bind_rows(subset(max_r2_afr, pop=='East'), r2_combined_apcdr)
apcdr_ukb_east_r2$pop <- recode_factor(apcdr_ukb_east_r2$pop, East='UKB East Africans')
pheno_order <- apcdr_ukb_east_r2 %>% arrange(desc(r2)) %>% ungroup() %>% select(pheno_code) %>% unique()
apcdr_ukb_east_r2$pheno_code <- factor(apcdr_ukb_east_r2$pheno_code, levels=as.character(pheno_order$pheno_code))

apcdr_r2 <- subset(apcdr_ukb_east_r2, pop=='UKB East Africans') %>% arrange(pheno_code)
east_r2 <- subset(apcdr_ukb_east_r2, pop=='Uganda') %>% subset(pheno_code !='HbA1c') %>% arrange(pheno_code)
dist_diffs <- data.frame(pheno=apcdr_r2$pheno_code,
           apcdr_higher=apcdr_r2$CI_2.5 - east_r2$CI_97.5,
           east_higher=east_r2$CI_2.5 - apcdr_r2$CI_97.5) %>%
  rowwise() %>%
  mutate(max_diff = max(apcdr_higher, east_higher))
dist_diffs <- dist_diffs %>%
  arrange(desc(max_diff))
apcdr_ukb_east_r2$pheno_code <- factor(apcdr_ukb_east_r2$pheno_code, levels=dist_diffs$pheno)
apcdr_ukb_east_r2 <- apcdr_ukb_east_r2 %>%
  group_by(pop) %>%
  mutate(median=median(r2), mad=mad(r2))

p_uganda_east <- ggplot(apcdr_ukb_east_r2, aes(x=pop, y=r2, fill=pop)) +
  geom_violin() +
  geom_point(aes(x=pop, y=median), color='black') +
  geom_errorbar(aes(ymin=median - mad, ymax=median + mad),color='black', width=0.1) +
  labs(x='Population', y=bquote(Maximum~trait~R^2)) +
  theme_bw() +
  guides(fill=F) +
  theme(axis.text = element_text(color='black', angle=45, hjust=1),
        text = element_text(size=16))
ggsave('apcdr_ukb_east_overall.pdf', p_uganda_east, height=5, width=4)

p3 <- ggplot(apcdr_ukb_east_r2, aes(x=pheno_code, y=r2, color=pop)) + #, shape=pop
  geom_point(position=pd) +
  #facet_grid(~sumstats) +
  geom_errorbar(aes(ymin=CI_2.5, ymax=CI_97.5), width=.5, position=pd) +
  scale_color_brewer(palette='Dark2', name='Dataset') +
  #scale_shape_manual(name='Population') +
  labs(x='Phenotype', y=bquote(Maximum~R^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10),
        axis.text = element_text(color='black'),
        text = element_text(size=16, color='black'),
        legend.justification = c(1, 1), legend.position = c(1, 1))
ggsave('apcdr_ukbb_east_afr_r2.pdf', p3, width=8, height=6)

date(); test <- compute_r2(all_prs, '50', 's5'); date()




# need to split up by population
# run this for all phenos, takes ~20 minutes for all phenotypes
date(); r2_combined <- ldply(pheno_test, function(x) { ldply(p_thresholds, function(y) { compute_r2(all_prs, x, y) } ) }, 
                             .parallel=T, .paropts=list(.export = c("compute_r2", "bootstrap_CIs", "p_thresholds", "pheno_test", "covariates"),
                                                        .packages = c('dplyr', 'plyr'))); date()

