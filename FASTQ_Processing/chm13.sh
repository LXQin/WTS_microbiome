# Define variables
base_dir="home/nfs/pengy3/microbiome/CHM13"
chm13_index="home/nfs/pengy3/microbiome/CHM13/index/chm13v2.0"
input_dir="home/nfs/pengy3/microbiome/hg38/unmapped_files"
bam_dir="home/nfs/pengy3/microbiome/CHM13/bam"
unmapped_dir="home/nfs/pengy3/microbiome/CHM13/unmapped_files"
sample_id="SRR8314071"

# Alignment with Bowtie2 and sorting with Samtools
bowtie2 -x "$chm13_index" -1 "${input_dir}/${sample_id}_unmapped_1.fastq" -2 "${input_dir}/${sample_id}_unmapped_2.fastq" -p 8 | \
samtools sort -O bam -o "${bam_dir}/${sample_id}_CHM13_sorted.bam"

# Index the sorted BAM file
samtools index "${bam_dir}/${sample_id}_CHM13_sorted.bam"

# Extract unmapped reads
samtools view -b -f 4 "${bam_dir}/${sample_id}_CHM13_sorted.bam" > "${bam_dir}/${sample_id}_CHM13_unmapped.bam"

# Convert BAM to FASTQ
bedtools bamtofastq -i "${bam_dir}/${sample_id}_CHM13_unmapped.bam" \
-fq "${unmapped_dir}/${sample_id}_final_unmapped_1.fastq" \
-fq2 "${unmapped_dir}/${sample_id}_final_unmapped_2.fastq"
