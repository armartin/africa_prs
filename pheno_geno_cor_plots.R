library(tidyverse)
library(corrplot)
library(ade4)

setwd('/Users/alicia/daly_lab/GINGER/prs/apcdr')


# Phenotype correlations (UKB) --------------------------------------------

pheno_codes <- read.csv('apcdr_ukb_pheno_info.csv', header=T) %>%
  select(pheno_code, ukb_code)
ukb_pheno_cor <- read.table(gzfile('../ukbb/pairwise_correlations.txt.bgz'), header=T)
ukb_pheno_cor <- merge(ukb_pheno_cor, pheno_codes, by.x='pheno_i', by.y='ukb_code')
colnames(ukb_pheno_cor)[ncol(ukb_pheno_cor)] <- 'pheno_i_name'
ukb_pheno_cor <- merge(ukb_pheno_cor, pheno_codes, by.x='pheno_j', by.y='ukb_code')
colnames(ukb_pheno_cor)[ncol(ukb_pheno_cor)] <- 'pheno_j_name'

ukb_cor <- ukb_pheno_cor %>% 
  filter(pheno_i_name %in% pheno_codes$pheno_code & pheno_j_name %in% pheno_codes$pheno_code) %>%
  filter(!pheno_i_name %in% c('LYMPHPr', 'NEUPr', 'BASOPr', 'EOSPr', 'MONOPr')) %>%
  filter(!pheno_j_name %in% c('LYMPHPr', 'NEUPr', 'BASOPr', 'EOSPr', 'MONOPr')) %>%
  select(c(entry, pheno_i_name, pheno_j_name)) %>%
  spread(key=pheno_i_name, value=entry) %>%
  select(-pheno_j_name)
rownames(ukb_cor) <- colnames(ukb_cor)
ukb_cor <- as.matrix(ukb_cor)
#ukb_cor <- as.matrix(ukb_cor^2)

pdf('pheno_corr_ukbb.pdf')
corrplot(ukb_cor, method='circle', order='hclust', tl.col='black')
dev.off()

ukb_clust_order <- colnames(ukb_cor)[corrMatOrder(ukb_cor, order="hclust")]
write.table(ukb_clust_order, 'ukb_clust_order.txt', row.names = F, col.names = F, quote = F)
ukb_clust_order <- as.character(read.table('ukb_clust_order.txt')$V1)

# Phenotype correlations (GPC) --------------------------------------------

phenos <- read.delim('gwas_phenotypes_28Oct14.txt', header=T, sep=' ')
phenos_only <- phenos %>% 
  select(SBP:BASO)
phenos_only <- phenos_only %>%
  select(noquote(order(colnames(phenos_only)))) %>%
  select(-c(LYMPHPr, NEUPr, BASOPr, EOSPr, MONOPr, WHR)) %>%
  mutate_if(is.factor, as.numeric) %>%
  select(ukb_clust_order)

M <- cor(phenos_only, use='pairwise.complete.obs')
colnames(M)[corrMatOrder(M, order="hclust")]

pdf('pheno_corr_apcdr.pdf')
corrplot(M, method='circle', tl.col='black')
dev.off()

mantel.rtest(dist(ukb_cor), dist(M), nrepet=9999)


# Heritability comparisons ------------------------------------------------

h2_apcdr <- read.table('h2_apcdr.txt', header=T) %>%
  mutate(study='Uganda GPC', se_low=h2-SE_h2, se_high=h2+SE_h2) %>%
  select(phenotype, h2, SE_h2, study, se_low, se_high) %>%
  filter(!(phenotype %in% c('HbA1c', 'WHR')))
h2_ukb <- read.table('h2_ukb.txt', header=T, sep='\t') %>%
  mutate(study='UK Biobank', phenotype=pheno_code, SE_h2=SE,
         se_low=h2-SE, se_high=h2+SE) %>%
  select(phenotype, h2, SE_h2, study, se_low, se_high) %>%
  filter(!(phenotype %in% c('LYMPHPr', 'NEUPr', 'BASOPr', 'EOSPr', 'MONOPr')))
h2_comb <- h2_apcdr %>%
  bind_rows(h2_ukb) %>%
  arrange(desc(h2))

h2_comb$phenotype <- factor(h2_comb$phenotype, levels=(subset(h2_comb, study=='Uganda GPC') %>% arrange(desc(h2)))$phenotype)
pd <- position_dodge(0.5)
p_h2 <- ggplot(h2_comb, aes(x=phenotype, y=h2, color=study)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=se_low, ymax=se_high), width=0, position=pd) +
  theme_classic() +
  scale_color_brewer(palette='Set1') +
  labs(x='Phenotype', y=expression(h^2)) +
  guides(color=guide_legend(title="Study")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10),
        axis.text = element_text(color='black'),
        text = element_text(size=16, color='black'),
        legend.justification = c(1, 1), legend.position = c(1, 1))

ggsave('h2_apcdr_ukb.pdf', p_h2, height=7, width=7)

# Genetic correlations (UKB) ----------------------------------------------

rg <- read.table('ldsc/ukb31063.gwas_holdout.ldsc.txt', header=T) %>%
  select(p1:rg)
rg2 <- rg %>%
  mutate(r1=p2, r2=p1) %>%
  mutate(p1=r1, p2=r2) %>%
  select(-c(r1,r2))
rg_comb <- rg %>%
  bind_rows(rg2)

rg_comb <- merge(rg_comb, pheno_codes, by.x='p1', by.y='ukb_code')
colnames(rg_comb)[ncol(rg_comb)] <- 'p1_name'
rg_comb <- merge(rg_comb, pheno_codes, by.x='p2', by.y='ukb_code')
colnames(rg_comb)[ncol(rg_comb)] <- 'p2_name'
rg1 <- data.frame(rg=1,p1_name=unique(rg_comb$p1_name), p2_name=unique(rg_comb$p1_name))
rg_comb <- rg_comb %>%
  select(-c(p1, p2)) %>%
  bind_rows(rg1) %>%
  filter(!p1_name %in% c('LYMPHPr', 'NEUPr', 'BASOPr', 'EOSPr', 'MONOPr')) %>%
  filter(!p2_name %in% c('LYMPHPr', 'NEUPr', 'BASOPr', 'EOSPr', 'MONOPr'))
  
rg_comb$p1_name <- factor(rg_comb$p1_name, levels=ukb_clust_order)
rg_comb$p2_name <- factor(rg_comb$p2_name, levels=ukb_clust_order)
# select(ukb_clust_order)
rg_spread <- rg_comb %>%
  spread(key=p1_name, value=rg) %>%
  select(-p2_name)
rownames(rg_spread) <- colnames(rg_spread)
rg_spread <- as.matrix(rg_spread)
#ukb_cor <- as.matrix(ukb_cor^2)

pdf('genetic_corr_ukbb.pdf')
corrplot(rg_spread, method='circle', tl.col='black')
dev.off()
