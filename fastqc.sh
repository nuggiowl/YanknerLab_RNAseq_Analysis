mkdir fastqc
cd fastqc
fastqc -o . --extract -f fastq -t 6 ../*.fastq.gz
multiqc .
