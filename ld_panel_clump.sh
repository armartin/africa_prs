#!/bin/bash
#$ -b y
#$ -wd /home/unix/armartin/atgu/armartin/ginger/apcdr/sumstats
#$ -w e
#$ -e /humgen/atgu1/fs03/armartin/ginger/apcdr/sumstats/clump_panel/logs
#$ -o /humgen/atgu1/fs03/armartin/ginger/apcdr/sumstats/clump_panel/logs
#$ -l h_vmem=8G
#$ -N clump
#$ -t 1-2

PHENO_LINE=`sed '1d' meta_analyze_ukb_bbj_ugr.txt | sed "${SGE_TASK_ID}q;d"`
PHENO=`echo $PHENO_LINE | awk '{ print $1 }'`
ANALYSIS_NAME=$1
CLUMP_PREFIX=clump_panel/${ANALYSIS_NAME}_${PHENO}

# make LD panel in plink
echo "Making LD panel"
echo ${CLUMP_PREFIX}

plink --bfile /humgen/atgu1/fs03/shared_resources/1kG/integrated/20130502/ALL.1KG_phase3.20130502.genotypes.maf005 \
--keep ${CLUMP_PREFIX}.keep \
--make-bed \
--out ${CLUMP_PREFIX}

# run clumping
echo "Running clumping"
plink --bfile ${CLUMP_PREFIX} \
--clump ${ANALYSIS_NAME}/${PHENO}.meta \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.1 \
--clump-kb 500 \
--out ${CLUMP_PREFIX}

# modify clumped files
echo "Modifying clumped files"
sed -i '/^\s*$/d' ${CLUMP_PREFIX}.clumped
bgzip -f ${CLUMP_PREFIX}.clumped
mv ${CLUMP_PREFIX}.clumped.gz ${CLUMP_PREFIX}.clumped.bgz

# clean up intermediate LD panel
echo "Cleaning up intermediate LD panel files"
rm ${CLUMP_PREFIX}.bed
rm ${CLUMP_PREFIX}.bim
rm ${CLUMP_PREFIX}.fam