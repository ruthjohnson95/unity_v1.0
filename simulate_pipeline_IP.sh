#!/usr/bin/env bash

#$ -cwd
#$ -j y
#$ -l h_data=5G,h_rt=2:00:00
#$ -o unity_v1.log

source /u/local/Modules/default/init/modules.sh
module load python/2.7

sim_yaml=$1

if [ -z "$sim_yaml" ]
then
    sim_yaml="sim.yaml"
fi 

# parse yaml file 
SIM_NAME=$(cat $sim_yaml | grep "SIM NAME" | cut -d':' -f2 | tr -d " \t\n\r" )
H_SIM=$(cat $sim_yaml | grep "H SIM" | cut -d':' -f2 | tr -d " \t\n\r" )
P_SIM=$(cat $sim_yaml | grep "P SIM" | cut -d':' -f2 | tr -d " \t\n\r" )
N=$(cat $sim_yaml | grep "SAMPLE SIZE" | cut -d':' -f2 | tr -d " \t\n\r" )
OUT_DIR=$(cat $sim_yaml | grep "OUTDIR" | cut -d':' -f2 | tr -d " \t\n\r" )
LD_LIST=$(cat $sim_yaml | grep "LD LIST" | cut -d':' -f2 | tr -d " \t\n\r" )
RSID_LIST=$(cat $sim_yaml | grep "RSID LIST" | cut -d':' -f2 | tr -d " \t\n\r" )
LD_DIR=$(cat $sim_yaml | grep "LD DIR" | cut -d':' -f2 | tr -d " \t\n\r" )
RSID_DIR=$(cat $sim_yaml | grep "RSID DIR" | cut -d':' -f2 | tr -d " \t\n\r" )
M=$(cat $sim_yaml | grep "NUM SNPS" | cut -d':' -f2 | tr -d " \t\n\r" )
SEED=$(cat $sim_yaml | grep "SEED" | cut -d':' -f2 | tr -d " \t\n\r" )
PLINK=$(cat $sim_yaml | grep "PLINK" | cut -d':' -f2 | tr -d " \t\n\r" )
REF_PATH=$(cat $sim_yaml | grep "REF BIM" | cut -d':' -f2 | tr -d " \t\n\r")

# make simulation output dir 
SIMDIR=$OUT_DIR'/'$SIM_NAME
mkdir -p $SIMDIR

# make header 
echo "- - - - - - - - - - UNITY v1.0 simulation - - - - - - - - -"
echo "SIM NAME: "$SIM_NAME
echo "H SIM: "$H_SIM
echo "P SIM: "$P_SIM
echo "N: "$N

if [ "$LD_LIST" = "NA" ] && [ "$RSID_LIST" = "NA" ]
then # no LD 
    echo "NUM SNPS: "$M
    echo "Not simulating with LD"
else
    echo "LD LIST: "$LD_LIST
    echo "RSID LIST: "$RSID_LIST
fi

echo "Outputing simulated gwas to directory: "$SIMDIR

echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"


# checks if LD_LIST and RSID_LIST is not NA 
if [ "$LD_LIST" = "NA" ] && [ "$RSID_LIST" = "NA" ]
then
    # simulate without LD 
    echo "no LD"
else
    # get chr/BP from plink
    for rsid_file in `cat $RSID_LIST`
    do
	DATE=`date '+%Y-%m-%d %H:%M:%S'`
	echo $DATE" Extracting CHR/BP information for snps: "$rsid_file
	snp_list=$RSID_DIR/$rsid_file
	file_prefix="${rsid_file%.*}"
	CHR=$(echo $file_prefix | cut -d'.' -f1 | sed 's/[^0-9]*//g' )
	$PLINK --bfile $REF_PATH.$CHR --extract $snp_list --make-just-bim --biallelic-only --out $SIMDIR/$file_prefix

    done
fi

# simulate gwas effects
python scripts/simulate.py --sim_name $SIM_NAME --h_sim $H_SIM --p_sim $P_SIM --N $N --rsid_list $RSID_LIST --ld_list $LD_LIST --outdir $SIMDIR --bim_dir $SIMDIR --ld_dir $LD_DIR --seed $SEED --M $M

# if LD, combine all files into one and then prune
if [ "$LD_LIST" != "NA" ] && [ "$RSID_LIST" != "NA" ]
then
    
    # contacenate all files into one 
    GWAS_FILE=$SIMDIR/$SIM_NAME'.gwas'

    echo $GWAS_FILE
    COUNTER=0
    for file in `ls $SIMDIR/chr*.gwas`
    do
	if [ "$COUNTER" -lt 1 ]
	then
	    # get the file header
	    cat $file > $GWAS_FILE
	else
	    # skip the file header 
	    tail -n +2 $file >> $GWAS_FILE
	fi

	COUNTER=$((COUNTER+1))

    done
fi 

# perform LD pruning 
WINDOW=5000
BP_HEAD="BP"
#python scripts/prune_ld.py --gwas_file $GWAS_FILE --window $WINDOW --bp_head $BP_HEAD

# run inference 
ITS=1000
GWAS_FILE=$SIMDIR/$SIM_NAME'.gwas'
python src/main.py  --s $SEED --H $H_SIM --N $N --id $SIM_NAME --its $ITS --gwas_file $GWAS_FILE --outdir $OUT_DIR
