#!/bin/bash

# Set base and fastq directory and reference genome database
## fastq_dir is the input file path (path to FASTQ), please move all the FASTQ files into this path (without any subfolders)
fastq_dir="/home/nfs/pengy3/microbiome/data"
## base_dir is the output file path, all the bam and report files generated will be stored here
base_dir="/home/nfs/pengy3/microbiome_test"
## Path to reference genome database
hg38_ref="/home/nfs/pengy3/microbiome/hg38/index/hg38.fa"
CHM13_ref="/home/nfs/pengy3/microbiome/CHM13/index/chm13v2.0"
kraken2_ref_16="/home/nfs/pengy3/microbiome/kraken2/16GB_standard_db"
kraken2_ref_8="/home/nfs/pengy3/microbiome/kraken2/8GB_standard_db"
## Common sample ID list
samples_dir="/home/nfs/pengy3/microbiome/sh_files/filter_ID.csv"

echo "Base Directory: $base_dir"
echo "FASTQ Directory: $fastq_dir"
echo "Reference Directory:"
echo "hg38: $hg38_ref"
echo "CHM13: $CHM13_ref"
echo "kraken2 16GB standard database: $kraken2_ref_16"
echo "kraken2 8GB standard database: $kraken2_ref_8"

# Read common sample ID list
sample_IDs=$(tail -n +2 "$samples_dir")
for file in "$fastq_dir"/*_1.fastq; do
  # Extract sample id
  sample_id=$(basename "$file" _1.fastq)
  
  # Check samples in the sample list have paired FASTQ files
  if [ -f "$fastq_dir/${sample_id}_2.fastq" ] && echo "$sample_IDs" | grep -qw "$sample_id"; then
    # 0. Output paired sample id
    echo "Processing paired sample $sample_id"
    
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
    
    # exclude mitochondria genes
    samtools view -h "$output_dir/bam/${sample_id}_sorted.bam" | grep -v "chrM" | samtools view -b -o "$output_dir/bam/${sample_id}_no_mit.bam" -
    samtools view -b -f 4 -f 8 "$output_dir/bam/${sample_id}_no_mit.bam" > "$output_dir/bam/${sample_id}_unmapped.bam"
    # conserve mitochondria genes
    # samtools view -b -f 4 -f 8 "$output_dir/bam/${sample_id}_sorted.bam" > "$output_dir/bam/${sample_id}_unmapped.bam"
    
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
    output_dir="$base_dir/kraken2"
    
    ## Create necessary directories
    if [ ! -d "$output_dir" ]; then
      echo "Creating directories needed for Kraken2."
      mkdir -p "$output_dir/reports"
      mkdir -p "$output_dir/counts"
    fi
    
    ## 16GB standard database of kraken2
    echo "Starting kraken2 classification"
    kraken2 --db "$kraken2_ref_16" --paired "$input_dir/${sample_id}_final_unmapped_1.fastq" "$input_dir/${sample_id}_final_unmapped_2.fastq" --report "$output_dir/reports/${sample_id}_report_16.txt"
    kraken2 --db "$kraken2_ref_8" --paired "$input_dir/${sample_id}_final_unmapped_1.fastq" "$input_dir/${sample_id}_final_unmapped_2.fastq" --report "$output_dir/reports/${sample_id}_report_8.txt"
    
    # 4. Delete unused files to minimize memory utility
    echo "Start to delete unnecessary intermediate files"
    find "$trimmed_dir" -type f -delete
    find "$base_dir/hg38" -type f -delete
    find "$base_dir/CHM13/unmapped_files" -type f -delete

    # 5. Bracken
    echo "From report to counts using bracken"
    bracken -d "$kraken2_ref_16" -i "$output_dir/reports/${sample_id}_report_16.txt" -o "$output_dir/counts/${sample_id}_counts_16.txt" -l S
    bracken -d "$kraken2_ref_8" -i "$output_dir/reports/${sample_id}_report_8.txt" -o "$output_dir/counts/${sample_id}_counts_8.txt" -l S
  fi
done
