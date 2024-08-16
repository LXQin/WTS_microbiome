# Define variables
base_dir=~/microbiome/CHM13
index_dir="$base_dir/index/chm13v2.0"
input_dir=~/microbiome/hg38/unmapped_files
bam_dir="$base_dir/bam"
unmapped_dir="$base_dir/unmapped_files"
sample_id="SRR8314071"

# Change to base directory
cd "$base_dir"

# Alignment with Bowtie2 and sorting with Samtools
bowtie2 -x "$index_dir" -1 "${input_dir}/${sample_id}_unmapped_1.fastq" -2 "${input_dir}/${sample_id}_unmapped_2.fastq" -p 8 | \
samtools sort -O bam -o "${bam_dir}/${sample_id}_CHM13_sorted.bam"

# Index the sorted BAM file
samtools index "${bam_dir}/${sample_id}_CHM13_sorted.bam"

# Extract unmapped reads
samtools view -b -f 4 "${bam_dir}/${sample_id}_CHM13_sorted.bam" > "${bam_dir}/${sample_id}_CHM13_unmapped.bam"

# Convert BAM to FASTQ
bedtools bamtofastq -i "${bam_dir}/${sample_id}_CHM13_unmapped.bam" \
-fq "${unmapped_dir}/${sample_id}_final_unmapped_1.fastq" \
-fq2 "${unmapped_dir}/${sample_id}
