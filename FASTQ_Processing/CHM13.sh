cd ~/microbiome/CHM13/
bowtie2 -x ./index/chm13v2.0 -1 ../hg38/processed_data/SRR8314071_unmapped_1.fastq -2 ../hg38/processed_data/SRR8314071_unmapped_2.fastq -S ./sam/SRR8314071_CHM13.sam -p 8
samtools view -Sb ./sam/SRR8314071_CHM13.sam | samtools sort -o ./bam/SRR8314071_CHM13_sorted.bam
samtools index ./bam/SRR8314071_CHM13_sorted.bam ./bam/SRR8314071_CHM13_sorted.bam.bai
samtools view -b -f 4 ./bam/SRR8314071_CHM13_sorted.bam > ./bam/SRR8314071_CHM13_unmapped.bam
bedtools bamtofastq -i ./bam/SRR8314071_CHM13_unmapped.bam -fq ./processed_data/SRR8314071_final_unmapped_1.fastq -fq2 ./processed_data/SRR8314071_final_unmapped_2.fastq
