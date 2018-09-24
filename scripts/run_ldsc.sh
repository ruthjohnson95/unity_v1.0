#!/bin/sh

LDSC_DIR=$1
HM3_PATH=$2
REF_LD_PATH=$3
TRAIT1=$4
TRAIT2=$5
file1=$6
file2=$7
N1=$8
N2=$9
OUTDIR=${10}

TRAIT1=$(basename $file1 | cut -d. -f1)
TRAIT2=$(basename $file2 | cut -d. -f1)


python $LDSC_DIR/munge_sumstats.py --sumstats $file1  --N $N1 --out $OUTDIR/'ldsc'.$TRAIT1 --merge-allele $HM3_PATH --ignore BETA,OR,SE,BETA_STD

python $LDSC_DIR/munge_sumstats.py --sumstats $file2  --N $N2 --out $OUTDIR/'ldsc'.$TRAIT2 --merge-allele $HM3_PATH --ignore BETA,OR,SE,BETA_STD

python $LDSC_DIR/ldsc.py --rg $OUTDIR/'ldsc'.$TRAIT1.sumstats.gz,$OUTDIR/'ldsc'.$TRAIT2.sumstats.gz --ref-ld-chr $REF_LD_PATH --w-ld-chr $REF_LD_PATH  --out $OUTDIR/'ldsc'.$TRAIT1.$TRAIT2