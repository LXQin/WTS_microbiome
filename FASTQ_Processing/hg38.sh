# Define variables
input_dir="/home/nfs/pengy3/microbiome/data"
trimmed_dir="/home/nfs/pengy3/microbiome/trimmed_files"
hg38_index="/home/nfs/pengy3/microbiome/hg38/index/hg38.fa"
bam_dir="/home/nfs/pengy3/microbiome/hg38/bam"
unmapped_dir="/home/nfs/pengy3/microbiome/hg38/unmapped_files"
sample_id="SRR8314071"

# Trim Galore for paired-end reads
trim_galore --paired "${input_dir}/${sample_id}_1.fastq" "${input_dir}/${sample_id}_2.fastq" -o "$trimmed_dir"

# Alignment with BWA and sorting with samtools
bwa mem -t 8 "$hg38_index" "${trimmed_dir}/${sample_id}_1_val_1.fq" "${trimmed_dir}/${sample_id}_2_val_2.fq" | \
samtools sort -@ 8 -O bam -o "${bam_dir}/${sample_id}_sorted.bam"

# Index the sorted BAM file
samtools index -@ 8 "${bam_dir}/${sample_id}_sorted.bam"

# Filter out mitochondrial reads (chrM) and save the result
samtools view -h "${bam_dir}/${sample_id}_sorted.bam" | grep -v "chrM" | \
samtools view -b -o "${bam_dir}/${sample_id}_no_mit.bam"

# Extract unmapped reads
samtools view -b -f 4 "${bam_dir}/${sample_id}_no_mit.bam" > "${bam_dir}/${sample_id}_unmapped.bam"

# Convert BAM to FASTQ
bedtools bamtofastq -i "${bam_dir}/${sample_id}_unmapped.bam" \
-fq "${unmapped_dir}/${sample_id}_unmapped_1.fastq" \
-fq2 "${unmapped_dir}/${sample_id}_unmapped_2.fastq"
