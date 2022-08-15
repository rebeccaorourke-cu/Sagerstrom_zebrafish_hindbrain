#!/bin/bash
#BSUB -J cellranger
#BSUB -o logs/cellranger.%J.out
#BSUB -e logs/cellranger.%J.err
#BSUB -n 12
#BSUB -R "select[mem>47] rusage[mem=47]"
#BSUB -R "span[hosts=1]" 

mkdir -p logs

. /usr/share/Modules/init/bash
module load modules modules-init
module load cellranger_arc/1.0.1 

cellranger-arc --version
REF=cellranger_arc_genomes/zebrafishPlusGFP/GRCz11

SAMPLES=(HB13hpf)
SAMPLE=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
ID=${SAMPLE}

cellranger-arc count --id=${ID} --reference=${REF} --libraries=library.csv --localmem 47 --localcores 12

