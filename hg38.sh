cd ~/microbiome/hg38/
bwa mem -t 8 ./index/hg38.fa ../data/SRR8314071_1.fastq ../data/SRR8314071_2.fastq > ./sam/SRR8314071.sam
samtools view -@ 8 -Sb ./sam/SRR8314071.sam | samtools sort -@ 8 -o ./bam/SRR8314071_sorted.bam
samtools index -@ 8 ./bam/SRR8314071_sorted.bam ./bam/SRR8314071_sorted.bam.bai
samtools view -b -f 4 ./bam/SRR8314071_sorted.bam > ./bam/SRR8314071_unmapped.bam
bedtools bamtofastq -i ./bam/SRR8314071_unmapped.bam -fq ./processed_data/SRR8314071_unmapped_1.fastq -fq2 ./processed_data/SRR8314071_unmapped_2.fastq