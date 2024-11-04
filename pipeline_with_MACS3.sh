#!/bin/bash
#SBATCH --job-name=fastqc_job
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --partition=open

source path/to/conda.sh
conda activate processing_atac-seq

# Define the input directory and output directory
input_dir="path/to/reads"
fastqc_output_dir="path/to/fastqc_output"
skewer_output_dir="path/to/skewer_output"
bwa_output_dir="path/to/bwa_output"
post_align_output_dir="path/to/post_align_output"
reference_genome="path/to/hg38.fa"
macs3_output_dir="path/to/macs3_output"   


post_alignment_script="path/to/post_alignment.py"

## FASTQC 

# Create the output directory if it doesn't exist
mkdir -p "$fastqc_output_dir"

# Loop through each R1 file in the input directory
for read1 in "$input_dir"/*_R1.fastq.gz; do
    # Get the base name of the file (without path and "_R1.fastq.gz")
    base_name=$(basename "$read1" _R1.fastq.gz)
    
    # Define the corresponding R2 file (should have the same base name)
    read2="$input_dir"/"${base_name}"_R2.fastq.gz
    
    # Check if both R1 and R2 files exist
    if [[ -f "$read2" ]]; then
        echo "Running FastQC for $read1 and $read2"

        # Run FastQC on both R1 and R2 and place the output in the output directory
        fastqc "$read1" "$read2" -o "$fastqc_output_dir"
    else
        echo "Missing R2 file for $read1. Skipping..."
    fi
done

echo "FastQC analysis complete. Results stored in $fastqc_output_dir."

## SKEWER 

mkdir -p "$skewer_output_dir"

for read1 in "$input_dir"/*_R1.fastq.gz; do
    # Get the base name of the file (without path and "_R1.fastq.gz")
    base_name=$(basename "$read1" _R1.fastq.gz)

    # Define the corresponding R2 file (should have the same base name)
    read2="$input_dir"/"${base_name}"_R2.fastq.gz

    # Check if both R1 and R2 files exist
    if [[ -f "$read2" ]]; then
        echo "Running Skewer for $read1 and $read2"

        # Run Skewer on the paired-end reads
        skewer -m pe -f sanger -o "$skewer_output_dir"/"${base_name}" "$read1" "$read2"

        echo "Skewer completed for $base_name"
    else
        echo "Missing R2 file for $read1. Skipping..."
    fi
done

echo "Skewer processing complete. Results stored in $skewer_output_dir."

## BWA

mkdir -p "$bwa_output_dir"

# Loop through the trimmed Skewer output files and run BWA MEM
for trimmed_read1 in "$skewer_output_dir"/*-trimmed-pair1.fastq; do
    # Get the base name of the file (without path and "-trimmed-pair1.fastq")
    base_name=$(basename "$trimmed_read1" -trimmed-pair1.fastq)
    
    # Define the corresponding R2 trimmed file
    trimmed_read2="$skewer_output_dir/${base_name}-trimmed-pair2.fastq"
    
    # Check if both trimmed paired-end files exist
    if [[ -f "$trimmed_read2" ]]; then
        echo "Running BWA MEM for $trimmed_read1 and $trimmed_read2"

        sam_output="$bwa_output_dir/${base_name}.sam"
        bam_output="$bwa_output_dir/${base_name}.bam"

        # Run BWA MEM on the trimmed reads and store the output in SAM format
        bwa mem -t 4 "$reference_genome" "$trimmed_read1" "$trimmed_read2" > "$bwa_output_dir/${base_name}.sam"

        echo "Converting SAM to BAM for $base_name"
        samtools view -Sb "$sam_output" | samtools sort -o "$bam_output"

        # Optional: Index the BAM file (useful for downstream analysis)
        echo "Indexing BAM file for $base_name"
        samtools index "$bam_output"

        echo "BWA MEM alignment and BAM conversion completed for $base_name"
    else
        echo "Missing trimmed R2 file for $trimmed_read1. Skipping..."
    fi
done

echo "BWA MEM and BAM processing complete. Results stored in $bwa_output_dir."

## Post Alignment Adustments 

mkdir -p "$post_align_output_dir"

for BAM_FILE in "$bwa_output_dir"/*.bam; do
    base_name=$(basename "$BAM_FILE" .bam)
    echo "Running post-alignment adjustments for $base_name"
    
    # Use the local variable to call post_alignment.py
    python "$post_alignment_script" -b "$BAM_FILE" -s "$base_name" -o "$post_align_output_dir"
done

echo "Post-alignment processing complete. Results stored in $post_align_output_dir."

## MACS3
mkdir -p "$macs3_output_dir"

for BAM_FILE in "$post_align_output_dir"/*.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    SAMPLE_OUTPUT_DIR="$macs3_output_dir/$SAMPLE_NAME"
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    echo "Processing sample: $SAMPLE_NAME with MACS3"

    # Run MACS3 callpeak
    macs3 callpeak -t "$BAM_FILE" -f BAMPE -g hs -n "$SAMPLE_NAME" --outdir "$SAMPLE_OUTPUT_DIR" --nomodel --shift -100 --extsize 200 -B$

    echo "MACS3 peak calling complete for $SAMPLE_NAME"
done

echo "MACS3 processing complete. Results stored in $macs3_output_dir."
