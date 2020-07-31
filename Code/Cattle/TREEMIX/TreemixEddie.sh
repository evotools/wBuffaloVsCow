#!/bin/bash
#$ -cwd
#$ -N TreeMix
#$ -l h_vmem=128G
#$ -pe sharedmem 1
#$ -R y
#$ -l h_rt=400:00:00

. /etc/profile.d/modules.sh
module load igmm/compilers/gcc/5.5.0
module load python/2.7.10
module load R
module load igmm/apps/plink/1.90b4


while getopts ":i:n:r:l:s:w:v" opt; do
  case $opt in
    i) indata=${OPTARG};;
    n) nmig=${OPTARG};;
    r) rscript=${OPTARG};;
    l) accountLD=${OPTARG};;
    s) size=${OPTARG};;
    w) winsizebp=${OPTARG};;
    v) version=${OPTARG};;
  esac
done


if [ $nmig == 0 ]; then
    opts=""
else
    opts="-m $SGE_TASK_ID"
fi

k="-k 1"
if [ $accountLD == 'y' ]; then
    if [ -e $indata.bim ]; then
		win=`wc -l $indata.bim | awk -v sz=$size -v len=$winsizebp '{print int(len/(sz/$1))}'`
    elif [ -e $indata.map ]; then
		win=`wc -l $indata.map | awk -v sz=$size -v len=$winsizebp '{print int(len/(sz/$1))}'`
    fi
    opts=$opts" -k "$win
    k="-k "$win
fi

mkdir './N_Migration_'${SGE_TASK_ID}_${accountLD}_${winsizebp}
echo "Active options: -i ./treemix_input.gz -o './N_Migration_'${SGE_TASK_ID}_${accountLD}_${winsizebp}/treemix_output_$SGE_TASK_ID $opts"

treemix -i ./treemix_input".gz" $opts -o  './N_Migration_'${SGE_TASK_ID}_${accountLD}_${winsizebp}/treemix_output_$SGE_TASK_ID
Rscript $rscript/Treemix_plots.R './N_Migration_'${SGE_TASK_ID}_${accountLD}_${winsizebp}/"treemix_output_"$SGE_TASK_ID poporder "N_Migration_"$SGE_TASK_ID_${accountLD}_${winsizebp}/PLOT_N_$SGE_TASK_ID $rscript

if [ $nmig == 0 ]; then
	threepop -i ./treemix_input.gz $k > threepopoutput.txt
	fourpop -i ./treemix_input.gz $k > fourpopoutput.txt
fi
