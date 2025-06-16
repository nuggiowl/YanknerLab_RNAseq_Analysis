input_bam=$1
ANNOT=$2

htseq-count -f bam -r pos -s reverse $input_bam $ANNOT/mm10/mm10.ncbiRefSeq.gtf
