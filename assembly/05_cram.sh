#!/bin/bash 

mkdir -p cram


R1=hic_R1.fastq.gz  
R2=hic_R2.fastq.gz

sample_name=santarcticus
prefix=hic

samtools import -@8 -r ID:${prefix} -r CN:OmniC -r PU:${prefix} -r SM:${sample_name} ${R1} ${R2} -o cram/${prefix}.cram --write-index
