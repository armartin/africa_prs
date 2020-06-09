library(tidyverse)
library(ggrepel)

setwd('/Users/alicia/daly_lab/GINGER/prs/ukbb/comparison_results')
bmi <- read.table(gzfile('/Users/alicia/daly_lab/GINGER/prs/ukbb/comparison_results/BMI.bgz'), header=T, sep='\t')

gwis <- bmi %>% 
  select(chrom:nearest_genes, p_UKB, p_BBJ_PAGE_UKB) %>% #pull out UKB only versus BBJ+PAGE+UKB meta
  filter(p_UKB < 5e-8 & p_BBJ_PAGE_UKB < 5e-8) %>% #keep only genome-wide significant associations in both datasets
  mutate(log10p_UKB=-log10(p_UKB), log10p_BBJ_PAGE_UKB=-log10(p_BBJ_PAGE_UKB)) %>% #convert to -log10(p)
  arrange(desc(log10p_BBJ_PAGE_UKB)) %>% #sort table from most to least significant for exploratory data analysis
  group_by(nearest_genes) %>%
  # create 1 label for each nearest gene cluster, color gene clusters the same
  mutate(label=ifelse(log10p_UKB > 100 & log10p_BBJ_PAGE_UKB==max(log10p_BBJ_PAGE_UKB), 
                      as.character(nearest_genes), NA),
         annot=ifelse(log10p_UKB > 100, as.character(nearest_genes), NA))


# make a QQ-like plot
p_cohorts <- ggplot(gwis, aes(x=log10p_UKB, y=log10p_BBJ_PAGE_UKB, color=annot, label=label)) + 
  geom_point() +
  geom_abline(slope=1) + 
  geom_text_repel() +
  #scale_color_brewer(palette='Set1') +
  theme_classic() +
  labs(x=expression(-log[10]*"(p)"~UKB),
       y=expression(-log[10]*"(p)"~BBJ+PAGE+UKB)) +
  guides(color=F) +
  theme(text = element_text(size=14))

ggsave('UKB_BBJ-PAGE-UKB_BMI_annot.png', p_cohorts)
