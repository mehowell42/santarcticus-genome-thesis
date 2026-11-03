SB=/fast/jd139945/santarcticus/analysis/earlgrey/earlgrey_sandbox
LIB=consensi.terrier_merged.fa
GENOME=GCA_035222895.1_eeStrDroe1_p1.0_genomic.fna

singularity exec -C --bind $(pwd):/data --bind $(pwd)/rm_tmp:/tmp --bind $(pwd)/rm_tmp:/var/tmp   --home /data --pwd /data --env TMPDIR=/tmp $SB   RepeatMasker -a -lib /data/$LIB -pa 8 -xsmall -gff -dir /data/rm_sneu /data/$GENOME 2>&1 | tee rm_sneu.log

RMDIR=$(singularity exec -C $SB bash -c 'dirname $(readlink -f $(which RepeatMasker))')

singularity exec -C --bind $(pwd):/data --home /data --pwd /data $SB a.align

singularity exec -C --bind $(pwd):/data --home /data --pwd /data $SB   perl $RMDIR/util/calcDivergenceFromAlign.pl -s /data/rm_sneu/sneu.divsum /data/rm_sneu/GCA_035222895.1_eeStrDroe1_p1.0_genomic.fna.align

singularity exec -C --bind $(pwd):/data --home /data --pwd /data $SB   perl $RMDIR/util/createRepeatLandscape.pl -div /data/rm_sneu/sneu.divsum -g 885292883 > rm_sneu/sneu_landscape.html
