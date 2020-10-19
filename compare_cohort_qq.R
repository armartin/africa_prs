#!/usr/bin/env Rscript

#library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(argparse)
library(RColorBrewer)

## get arguments
parser <- ArgumentParser()
parser$add_argument('--pheno', type='character')
parser$add_argument('--cohort1', type='character')
parser$add_argument('--cohort2', type='character')
args <- parser$parse_args()

print(paste(args$pheno, args$cohort1, args$cohort2, sep=', '))

#setwd('/Users/alicia/daly_lab/GINGER/prs/apcdr/prs/sumstats')

date(); comb <- read.delim(gzfile(paste0(args$pheno, '.bgz')), header=T); date()

# choose a downsampling bin size for plotting purposes, do filtering and downsampling
n_bins = 10
comb_filt <- comb %>%
  filter(!is.na(p_UKB) & !is.na(p_BBJ_PAGE_UKB)) %>%
  filter(p_UKB < 0.5 & p_BBJ_PAGE_UKB < 0.5 & p_UKB > 0 & p_BBJ_PAGE_UKB != 0) %>%
  mutate(log10p_UKB=-log10(p_UKB), log10p_BBJ_PAGE_UKB=-log10(p_BBJ_PAGE_UKB)) %>%
  mutate(diff_p = log10p_BBJ_PAGE_UKB - log10p_UKB) %>% #
  #arrange(desc(log10p_BBJ_PAGE_UKB)) %>%
  arrange(desc(diff_p)) %>%
  mutate(x_grouping=ntile(log10p_UKB, n_bins),
         y_grouping=ntile(log10p_BBJ_PAGE_UKB, n_bins))

# find genome-wide significant associations
# fit model, identify associations above prediction interval
gwis <- comb_filt %>%
  filter(p_BBJ_PAGE_UKB < 5e-8 & p_UKB < 5e-8) %>%
  #arrange(desc(log10p_BBJ_PAGE_UKB))
  arrange(desc(diff_p))
gwis.lm <- lm(log10p_BBJ_PAGE_UKB ~ log10p_UKB, data=gwis)
pred <- as.data.frame(predict(gwis.lm,se.fit=TRUE, level=0.99,interval="predict"))
above <- data.frame(gwis,
                   pred_above = gwis$log10p_BBJ_PAGE_UKB > pred$fit.upr) %>%
  filter(pred_above)

# label top variant per top 10 genes with name
above$nearest_genes <- factor(above$nearest_genes, levels=as.character(unique(above$nearest_genes)))
n_genes = 10
top_variant_labs <- above %>%
  group_by(nearest_genes) %>%
  slice(1) %>%
  head(n_genes)

sample_size = 500
to_downsample <- comb_filt %>%
  count(x_grouping, y_grouping) %>%
  mutate(downsample = case_when(n > sample_size & (x_grouping < n_bins | y_grouping < n_bins) ~ TRUE,
                   TRUE ~ FALSE)) %>%
  ungroup()
  
date(); comb_filt_downsample <- comb_filt %>% 
  left_join(to_downsample) %>%
  filter(downsample) %>%
  group_by(x_grouping, y_grouping) %>%
  sample_n(sample_size) %>%
  ungroup() %>%
  union_all(comb_filt %>%
              left_join(to_downsample) %>%
              filter(!downsample)) %>%
  mutate(lab_gene_name = case_when(chrom %in% top_variant_labs$chrom & pos %in% top_variant_labs$pos &
           ref %in% top_variant_labs$ref & alt %in% top_variant_labs$alt ~ nearest_genes),
         color_variants = case_when(chrom %in% above$chrom & pos %in% above$pos &
           ref %in% above$ref & alt %in% above$alt ~ nearest_genes)); date()
comb_filt_downsample$color_variants <- as.character(comb_filt_downsample$color_variants)
comb_filt_downsample$color_variants[which(is.na(comb_filt_downsample$color_variants))] <- 'NA'

color_vec <- rep('grey', length(unique(comb_filt_downsample$color_variants)))
names(color_vec) <- as.character(unique(comb_filt_downsample$color_variants))
color_labels <- as.character(unique(subset(comb_filt_downsample, !is.na(lab_gene_name))$lab_gene_name))
color_vec[which(names(color_vec) %in% color_labels)] <- colorRampPalette(brewer.pal(7, 'Set1'))(length(color_labels))

p1 <- ggplot(comb_filt_downsample, aes(x=log10p_UKB, y=log10p_BBJ_PAGE_UKB, color=color_variants, label=lab_gene_name)) + 
  geom_point()  +
  geom_abline(intercept=0) +
  geom_text_repel(nudge_y=10) +
  scale_color_manual(values=color_vec) +
  labs(x=expression(-log[10]*"(p)"~UKB),
       y=expression(-log[10]*"(p)"~UKB+BBJ+PAGE)) +
  theme_classic() +
  guides(color=F) +
  theme(text = element_text(size=16),
        axis.text = element_text(color='black'))

ggsave(paste0(args$pheno, '_', args$cohort1, '_vs_', args$cohort2, '.png'), p1, height=7, width=7)
