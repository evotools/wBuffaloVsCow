#!/bin/bash
#$ -N phylip
#$ -cwd                  
#$ -l h_rt=11:59:59 
#$ -pe sharedmem 1
#$ -l h_vmem=32G
#$ -R y
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -hold_jid prepibs,ibsdist

# Initialise the environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/plink/1.90b4
module load python/2.7.10


while getopts ":o:n" opt; do
  case $opt in
    o) outgr=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done

export PATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/phylip/exe:$PATH

bootstr=`ls infile_* | wc -l | awk '{print $1}'`
cat infile_* > infile
zip -9 AllDistMatrices.zip infile_*
rm infile_*
if [ -z "$outgr" ]; then
        echo "Estimate unrooted tree."
        echo -e "M\n$bootstr\n135\ny\n" | neighbor
        mkdir Neighbor
        mv out* Neighbor
        cp Neighbor/outtree ./intree
        echo -e "y\n" | consense
else
        echo "Estimate rooted tree on: "$outgr
        ogr_n=`head -1 Outgroups_indexes.txt`
        echo -e "O\n$ogr_n\nM\n$bootstr\n135\ny\n" | neighbor
        mkdir Neighbor
        mv out* Neighbor
        cp Neighbor/outtree ./intree
        echo -e "R\ny\n" | consense
fi
echo "Done."
