#!/bin/bash
# All SGE arguments
#$ -N flo
#$ -cwd
#$ -pe sharedmem 12
#$ -l h_vmem=8G
#$ -l h_rt=160:00:00
#$ -R y
#$ -o ./flo.out
#$ -e ./flo.err
#$ -P roslin_ctlgh

sourcefa=$1
tgtfa=$2

src=`realpath sourcefa`
tgt=`realpath tgtfa`


sed "s/THREADS/$NSLOTS/g" flo_opts.template.yaml > flo_opts.tmp1.yaml
sed "s#SRCFA#$src#g" flo_opts.tmp1.yaml > flo_opts.tmp2.yaml
sed "s#TGTFA#$tgt#g" flo_opts.tmp2.yaml > flo_opts.yaml && rm flo_opts.tmp1.yaml flo_opts.tmp2.yaml
rake -f /PATH/TO/flo/Rakefile`