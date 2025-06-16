BASE=`basename $1 _1.fq.gz`
REF=$3
ANNOT=$4

STAR --genomeDir $3/genome_mm10/ \
--sjdbGTFfile $4/mm10/mm10.ncbiRefSeq.gtf \
--runThreadN 8 \
--outBAMsortingThreadN 8 \
--readFilesCommand zcat \
--readFilesIn $1 $2 \
--outFileNamePrefix ./${BASE}_mm10_ \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--outFilterMultimapNmax 20 \
--outFilterType BySJout \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--limitBAMsortRAM 40000000000

samtools index ./${BASE}_mm10_Aligned.sortedByCoord.out.bam
