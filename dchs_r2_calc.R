library(tidyverse)

setwd('/Users/alicia/daly_lab/GINGER/prs/drakenstein')

# Read in phenos, PCA, and PRS
phenos <- read.delim('Ginger_project.data_19AUG2020.txt', header=T)
pca_file <- read.table('pca/DCHS_OUTCOME_PCA20.eigenvec', header=F) #Need to rename "V2" to "IID" in PCA file to match PRS & phenos files
colnames(pca_file) <- c("V1","IID", paste0('PC', 1:20))
ethnicity <- read.delim('DCHS_ethnicity_15July2020.txt', header=T)

# get a vector of all phenotypes
pheno_names <- phenos %>% dplyr::select(-BARCODE) %>% colnames(.)
p <- rev(c('all', 'p5', 'p2', 'p1', 'p05', 'pe2', 'pe3', 'pe4', 'pe6', 'p5e8'))
covariates <- paste0('PC', 1:20)
prs_prefixes <- c('depression_MTAG_DEP', 'depression_ex23andMe', 'depression_to10K', 'education_excl23andMe', 
            'education_to10K', 'neuroticism_GWAS', 'neuroticism_MTAG', 'ptsd_all_freeze2', 
            'substances_DrinksPerWeek', 'substances_SmokingInitiation', 'well_being_MTAG_SWB', 'all_height')
ethnicities <- as.character(unique(ethnicity$Ethnicity))

# read in how to pair up phenos and PRS (i.e. GWAS discovery analyses)
pheno_prs_pairs <- read.delim('pheno_prs_pairs.txt', header=T)

read_prs <- function(pheno, p_threshold) {
  filename = paste0('prs/', pheno, '_prs_out.', p_threshold, '.profile')
  print(filename)
  if(file.exists(filename)){
    prs_depression <- read.table(filename, header=T) %>%
      mutate(p_threshold=p_threshold, prs_prefix=pheno)
  }
}

# read in all PRS for all phenos
argList = list(prs_prefixes = prs_prefixes, p = p)
crossArg = cross_df(argList) # need to cross args to get all combos read
prs_phenos <- pmap(list(crossArg$prs_prefixes, crossArg$p), read_prs)
prs_all_phenos <- do.call(rbind.data.frame, prs_phenos) # bind together data frame (not list)
# prs_depression_to10K <- subset(prs_all_phenos, pheno=='depression_to10K') # get a test pheno

# combine measured phenos with PCA and ethnicity data. remove individuals geno'd twice: 101239, 300134
pheno_pca_ethnicity <- phenos %>%
  mutate('IID'=BARCODE) %>%
  left_join(pca_file, by='IID') %>%
  left_join(ethnicity, by='IID') %>%
  subset(!BARCODE %in% c('101239', '300134')) %>% # these individuals show up twice in PCA and have two different PCA positions. 
# Different secondary IDs (labeled V1). Hmm...
  gather('pheno_name', 'true_pheno', maternal_height_anc:Prenatal_smoking_cotinine2)
  # sex is all female (moms)
# they should have age. 

prs_pheno_pca_ethnicity <- pheno_pca_ethnicity %>%
  left_join(prs_all_phenos, by='IID')

# Bootstrap 95% CI around r^2
bootstrap_CIs <- function(my_subset) {
  sample_subset <- sample_n(my_subset, nrow(my_subset), replace=TRUE)
  model1 <- lm(as.formula(sprintf("true_pheno ~ SCORE + %s", paste(covariates, collapse = " + "))), data=sample_subset)
  model0 <- lm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=sample_subset)
  sample_r2 = summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  return(sample_r2)
}

# Get R^2 and p-value for nested model
compute_r2 <- function(p, prs_prefix, pheno, ethnicity, dataset=prs_pheno_pca_ethnicity) { #need R dataframe for each phenotype
  # input requirements: dataset is a data frame with the phenotype of interest named true_pheno
  # all PRS tested (incl column labeled p_threshold are included in data frame
  print(paste(p, prs_prefix, pheno, ethnicity))
  my_subset <- subset(dataset, p_threshold==p & prs_prefix==prs_prefix & pheno_name==pheno & Ethnicity==ethnicity)
  model1 <- lm(as.formula(sprintf("true_pheno ~ SCORE + %s", paste(covariates, collapse = " + "))), data=my_subset)
  model0 <- lm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=my_subset)
  r2 <- summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  pval <- anova(model1, model0)$'Pr(>F'[2]
  date(); bs_means <- replicate(100, bootstrap_CIs(my_subset)); date()
  delta <- bs_means - r2
  d = quantile(delta, c(0.025, 0.975), na.rm=T)
  ci = r2 - c(d[1], d[2])
  r2_vec <- c(r2, rev(ci), pval, p) #phenotype, 
  names(r2_vec) <- c('r2', 'CI_2.5', 'CI_97.5', 'p', 'p_threshold') #'pheno', 
  return(data.frame(t(r2_vec)))
  return(r2)
}

# make a 2nd argList for computing r2
argList = list(p_thresholds = p, prs = prs_prefixes, pheno = pheno_names, ethnicities = ethnicities)
# possibly need to match up list of prs prefixes and phenos
crossArg = cross_df(argList) %>% # need to cross args to get all combos read
  inner_join(pheno_prs_pairs, by=c('prs'='prs_prefix', 'pheno'='pheno'))
date(); prs_r2 <- pmap(list(crossArg$p_thresholds, crossArg$prs, crossArg$pheno, crossArg$ethnicities), compute_r2); date()
prs_r2_all <- do.call(rbind.data.frame, prs_r2) # bind together data frame (not list)
prs_r2_all <- as.data.frame(crossArg) %>% dplyr::select(-p_thresholds) %>% bind_cols(prs_r2_all)
write.table(prs_r2_all, 'dchs_r2.txt', row.names=F, quote=F, sep='\t')


# Plot R^2 for all phenos -------------------------------------------------

prs_r2_all <- read.delim('dchs_r2.txt', header=T) %>%
  left_join(pheno_prs_pairs, by=c('prs'='prs_prefix', 'pheno'='pheno')) %>%
  mutate(Ethnicity=case_when(ethnicities=='African ancestry' ~ 'Black/African',
                            ethnicities=='Mixed ancestry' ~ 'Mixed')) %>%
  filter(pheno_rename != 'Education')

prs_r2_all$pheno_rename <- factor(prs_r2_all$pheno_rename, levels=unique(prs_r2_all$pheno_rename)[c(5,1,2,3,4)])

prs_r2_top <- prs_r2_all %>%
  group_by(pheno_rename, Ethnicity) %>%
  arrange(p) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(r2)

p1 <- ggplot(prs_r2_top, aes(x=Ethnicity, y=r2, color=Ethnicity)) +
  facet_wrap(~pheno_rename, nrow=2) +
  geom_errorbar(aes(ymin=CI_2.5, ymax=CI_97.5), width=0.1) +
  geom_point() +
  labs(x='Phenotypes', y=bquote(Max~R^2)) +
  scale_color_brewer(palette='Set1') +
  theme_classic() +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(color='black'))

ggsave('dchs_multipheno_r2_all2.pdf', p1, height=6, width=10)
ggsave('dchs_multipheno_r2_all2.png', p1, height=6, width=10)


# Height by tertiles ------------------------------------------------------

# split up Mixed ancestry by PC1
mixed_height <- prs_pheno_pca_ethnicity %>%
  filter(Ethnicity=='Mixed ancestry' & pheno_name=='maternal_height_anc' & 
           p_threshold=='p05' & prs_prefix=='all_height') %>%
  mutate(Ethnicity=ntile(PC1, 2))
  
# make a 2nd argList for computing r2
argList = list(p_thresholds = 'p05', prs = 'all_height', pheno = 'maternal_height_anc', ethnicities = 1:2)
# possibly need to match up list of prs prefixes and phenos
crossArg = cross_df(argList) %>% # need to cross args to get all combos read
  inner_join(pheno_prs_pairs, by=c('prs'='prs_prefix', 'pheno'='pheno'))
date(); prs_r2 <- pmap(list(crossArg$p_thresholds, crossArg$prs, crossArg$pheno, crossArg$ethnicities), dataset=mixed_height, compute_r2); date()
prs_r2_all <- do.call(rbind.data.frame, prs_r2) # bind together data frame (not list)
prs_r2_all <- as.data.frame(crossArg) %>% dplyr::select(-p_thresholds) %>% bind_cols(prs_r2_all)


# Plot PCA ethnicity data -------------------------------------------------

pca_ethnicity <- pheno_pca_ethnicity %>% 
  select(-c(pheno_name, true_pheno)) %>%
  mutate(Ethnicity=case_when(Ethnicity=='Black/African ancestry' ~ 'Black/African',
                             Ethnicity=='Mixed ancestry' ~ 'Mixed'))
  unique()

p2 <- ggplot(pca_ethnicity, aes(x=PC1, y=PC2, color=Ethnicity)) +
  geom_point() +
  scale_color_brewer(palette='Set1') +
  theme_classic() +
  theme(text = element_text(size=14))

ggsave('dchs_pheno_ethnicity_pc1_2.pdf', p2, height=6, width=6)
