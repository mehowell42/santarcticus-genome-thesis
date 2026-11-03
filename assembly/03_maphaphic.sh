awk '/^S/{print ">"$2;print $3}' sa_hifism_omnic/santarcticus.asm.hic.hap1.p_ctg.gfa > sa_hap1.hic.p_ctg.fa

awk '/^S/{print ">"$2;print $3}' sa_hifism_omnic/santarcticus.asm.hic.hap2.p_ctg.gfa > sa_hap2.hic.p_ctg.fa

cat sa_hap1.hic.p_ctg.fa sa_hap2.hic.p_ctg.fa > sa_hap12.hic.p_ctg.fa

bwa index sa_hap12.hic.p_ctg.fa

hicdir="AGRF_CAGRF25040459-2_2353N2LT3/"

cat AGRF_STAGSCAGRF26020170-2_23JGHTLT3/PC00173-A6_23JGHTLT3_CTGATCGTAT-ACATATGCGC_L003_R1.fastq.gz > hic_R1.fastq.gz
cat AGRF_STAGSCAGRF26020170-2_23JGHTLT3/PC00173-A6_23JGHTLT3_CTGATCGTAT-ACATATGCGC_L003_R2.fastq.gz > hic_R2.fastq.gz


bwa mem -5SP -t 28 sa_hap12.hic.p_ctg.fa hic_R1.fastq.gz hic_R2.fastq.gz | \
	samblaster | samtools view - -@ 28 -S -h -b -F 3340 -o HiC.bam
