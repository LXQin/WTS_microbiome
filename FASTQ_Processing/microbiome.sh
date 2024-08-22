#!/bin/bash

# Set base and fastq directory and reference genome database
fastq_dir="/home/nfs/pengy3/microbiome/data"
base_dir="/home/nfs/pengy3/microbiome_test"
hg38_ref="/home/nfs/pengy3/microbiome/hg38/index/hg38.fa"
CHM13_ref="/home/nfs/pengy3/microbiome/CHM13/index/chm13v2.0"
kraken2_ref="/home/nfs/pengy3/microbiome/kraken2/16GB_standard_db"

echo "Base Directory: $base_dir"
echo "FASTQ Directory: $fastq_dir"
echo "Reference Directory:"
echo "hg38: $hg38_ref"
echo "CHM13: $CHM13_ref"
echo "kraken2 standard database: $kraken2_ref"

for file in "$fastq_dir"/*_1.fastq; do
  # Extract sample id
  sample_id=$(basename "$file" _1.fastq)
  if [ -f "$fastq_dir/${sample_id}_2.fastq" ]; then
    # 0. Output paired sample id
    echo "Processing sample with files:"
    echo "$fastq_dir/${sample_id}_1.fastq"
    echo "$fastq_dir/${sample_id}_2.fastq"
    
    # 1. hg38 & trimmed
    ## Set directory
    input_dir="$fastq_dir"
    output_dir="$base_dir/hg38"
    trimmed_dir="$base_dir/trimmed_files"

    ## Create necessary directories
    if [ ! -d "$output_dir" ]; then
      echo "Creating directories needed for hg38."
      mkdir -p "$output_dir/bam"
      mkdir -p "$output_dir/unmapped_files"
      mkdir -p "$trimmed_dir"
    fi
    
    ## trimmed
    echo "Starting trimming process"
    trim_galore --paired "$input_dir/${sample_id}_1.fastq" "$input_dir/${sample_id}_2.fastq" -o $trimmed_dir
    trimmed_r1="${trimmed_dir}/${sample_id}_1_val_1.fq"
    trimmed_r2="${trimmed_dir}/${sample_id}_2_val_2.fq"
    
    ## bwa
    echo "Starting bwa"
    bwa mem -t 8 "$hg38_ref" "$trimmed_r1" "$trimmed_r2" | samtools sort -@ 8 -O bam -o "$output_dir/bam/${sample_id}_sorted.bam"
    samtools index -@ 8 "$output_dir/bam/${sample_id}_sorted.bam"
    samtools view -h "$output_dir/bam/${sample_id}_sorted.bam" | grep -v "chrM" | samtools view -b -o "$output_dir/bam/${sample_id}_no_mit.bam" -
    samtools view -b -f 4 -f 8 "$output_dir/bam/${sample_id}_no_mit.bam" > "$output_dir/bam/${sample_id}_unmapped.bam"
    bedtools bamtofastq -i "$output_dir/bam/${sample_id}_unmapped.bam" -fq "$output_dir/unmapped_files/${sample_id}_unmapped_1.fastq" -fq2 "$output_dir/unmapped_files/${sample_id}_unmapped_2.fastq"

    # 2. CHM13
    ## Set directory
    input_dir="$base_dir/hg38/unmapped_files"
    output_dir="$base_dir/CHM13"
    
    ## Create necessary directories
    if [ ! -d "$output_dir" ]; then
      echo "Creating directories needed for CHM13."
      mkdir -p "$output_dir/bam"
      mkdir -p "$output_dir/unmapped_files"
    fi
    
    ## CHM13
    echo "Starting CHM13 alignment"
    bowtie2 -x "$CHM13_ref" -1 "$input_dir/${sample_id}_unmapped_1.fastq" -2 "$input_dir/${sample_id}_unmapped_2.fastq" -p 8 | samtools sort -O bam -o "$output_dir/bam/${sample_id}_CHM13_sorted.bam"
    samtools index "$output_dir/bam/${sample_id}_CHM13_sorted.bam"
    samtools view -b -f 4 -f 8 "$output_dir/bam/${sample_id}_CHM13_sorted.bam" > "$output_dir/bam/${sample_id}_CHM13_unmapped.bam"
    bedtools bamtofastq -i "$output_dir/bam/${sample_id}_CHM13_unmapped.bam" \
      -fq "$output_dir/unmapped_files/${sample_id}_final_unmapped_1.fastq" -fq2 "$output_dir/unmapped_files/${sample_id}_final_unmapped_2.fastq"

    # 3. kraken2
    ## Set directory
    input_dir="$base_dir/CHM13/unmapped_files"
    output_dir="$base_dir/kraken2/reports"
    
    ## Create necessary directories
    if [ ! -d "$output_dir" ]; then
      echo "Creating directories needed for Kraken2."
      mkdir -p "$output_dir"
    fi
    
    ## 16GB standard database of kraken2
    echo "Starting kraken2 classification"
    kraken2 --db "$kraken2_ref" --paired "$input_dir/${sample_id}_final_unmapped_1.fastq" "$input_dir/${sample_id}_final_unmapped_2.fastq" --report "$output_dir/${sample_id}_report.txt"
    
    # 4. Delete unused files to minimize memory utilizity
    find "$trimmed_dir" -type f -delete
    find "$base_dir/hg38" -type f -delete
    find "$base_dir/CHM13/unmapped_files" -type f -delete
  fi
done
