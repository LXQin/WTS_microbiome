base_dir="/home/nfs/pengy3/microbiome"
input_dir="/home/nfs/pengy3/microbiome/CHM13/unmapped_files"
minikraken2_db="/home/nfs/pengy3/microbiome/kraken2/minikraken2"
sample_id="SRR8314071"
kraken2 --db "$minikraken2_db" --paired "${input_dir}/${sample_id}_final_unmapped_1.fastq" "${input_dir}/${sample_id}_final_unmapped_2.fastq" --report minikraken2_report.txt
bracken -d "$minikraken2_db" -i minikraken2_report.txt -o bracken_species_counts.txt -l S
