#!/bin/bash

#Author: Ajinkya Khilari
# Ensure correct number of arguments
if [ "$#" -ne 10 ]; then
  echo "Usage: $0 <path_to_sample_directories> <path_to_reference_genome> <path_to_primer_bed_file> <threads> <min_length> <max_length> <clair3_env> <clair3_model> <quality_score> <output_directory>"
  exit 1
fi

# Assign command-line arguments to variables
SAMPLE_DIR_PATH="$1"
REF_GENOME="$2"
PRIMER_BED="$3"
THREADS="$4"
MIN_LENGTH="$5"
MAX_LENGTH="$6"
CLAIR3_ENV_NAME="$7"    # Clair3 environment name
CLAIR3_MODEL="$8"
QUALITY_SCORE="$9"
OUTPUT_DIR="${10}"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Process each sample
echo "Processing samples from ${SAMPLE_DIR_PATH}" 
for sample_dir in "${SAMPLE_DIR_PATH}"/*/; do
  barcode=$(basename "$sample_dir")
  echo "Processing $barcode..."

  barcode_output_dir="${OUTPUT_DIR}/${barcode}"
  mkdir -p "${barcode_output_dir}"

  # Merge FASTQ files
  cat "${sample_dir}"*.fastq > "${barcode_output_dir}/${barcode}.fastq"

  # Quality filtering with fastp (using fastp_env)
  conda run -n fastp_env fastp -i "${barcode_output_dir}/${barcode}.fastq" \
        -o "${barcode_output_dir}/${barcode}_filt.fastq" \
        -q "${QUALITY_SCORE}" -l "${MIN_LENGTH}" --max_len1 "${MAX_LENGTH}" --thread "${THREADS}"
  if [[ $? -ne 0 ]]; then
    echo "Error: fastp failed for $barcode" >> "${barcode_output_dir}/error.log"
    continue
  fi

  # Index reference genome (using samtools_env)
  conda run -n samtools_env samtools faidx "${REF_GENOME}"

  # Minimap2 alignment (using minimap2_env and samtools_env)
  conda run -n minimap2_env minimap2 -ax map-ont -t "${THREADS}" "${REF_GENOME}" "${barcode_output_dir}/${barcode}_filt.fastq" | \
  conda run -n samtools_env samtools view -bS -F 4 - | conda run -n samtools_env samtools sort -@ "${THREADS}" -o "${barcode_output_dir}/${barcode}.sorted.bam"
  if [[ $? -ne 0 ]]; then
    echo "Error: minimap2 alignment failed for $barcode" >> "${barcode_output_dir}/error.log"
    continue
  fi

  # Index BAM file (using samtools_env)
  conda run -n samtools_env samtools index "${barcode_output_dir}/${barcode}.sorted.bam"

  # Primer trimming with ivar (using ivar_env)
  conda run -n ivar_env ivar trim -i "${barcode_output_dir}/${barcode}.sorted.bam" -b "${PRIMER_BED}" -e -m 1 -s 4 -q 0 -p "${barcode_output_dir}/${barcode}.primerclipped"
  if [[ $? -ne 0 ]]; then
    echo "Error: ivar trim failed for $barcode" >> "${barcode_output_dir}/error.log"
    continue
  fi

  # Sort & index primer-trimmed BAM (using samtools_env)
  conda run -n samtools_env samtools sort -@ "${THREADS}" "${barcode_output_dir}/${barcode}.primerclipped.bam" -o "${barcode_output_dir}/${barcode}.primerclipped.sorted.bam"
  conda run -n samtools_env samtools index "${barcode_output_dir}/${barcode}.primerclipped.sorted.bam"

  # Clair3 Variant Calling (using clair3_env)
  conda run -n "${CLAIR3_ENV_NAME}" run_clair3.sh --bam_fn "${barcode_output_dir}/${barcode}.primerclipped.sorted.bam" \
                --ref_fn "${REF_GENOME}" \
                --model_path "${CLAIR3_MODEL}" \
                --threads "${THREADS}" \
                --platform ont \
                --output "${barcode_output_dir}/clair3_output" \
                --include_all_ctgs --haploid_sensitive
  if [[ $? -ne 0 ]]; then
    echo "Error: Clair3 variant calling failed for $barcode" >> "${barcode_output_dir}/error.log"
    continue
  fi

  # Index Clair3 VCF (using tabix_env)
  if [[ -f "${barcode_output_dir}/clair3_output/merge_output.vcf.gz" ]]; then
    conda run -n tabix_env tabix -p vcf "${barcode_output_dir}/clair3_output/merge_output.vcf.gz" -f
    if [[ $? -ne 0 ]]; then
      echo "Error: Tabix VCF indexing failed for $barcode" >> "${barcode_output_dir}/error.log"
      continue
    fi
  else
    echo "Error: Clair3 VCF file not found for $barcode" >> "${barcode_output_dir}/error.log"
    continue
  fi

  # Identify Low-Quality Variant Positions & Create Mask File (using bcftools_env)
  conda run -n bcftools_env bcftools query -f '%CHROM\t%POS0\t%POS\n' \
                 -i 'QUAL<10 || FORMAT/GQ<10 || FORMAT/AF<0.5' \
                 "${barcode_output_dir}/clair3_output/merge_output.vcf.gz" > \
                 "${barcode_output_dir}/${barcode}.variant_mask.bed"

  # Identify Low-Coverage Regions (<20X) & Add to Mask (using samtools_env)
  conda run -n samtools_env samtools depth -a "${barcode_output_dir}/${barcode}.primerclipped.sorted.bam" | \
  awk '$3 < 20 {print $1, $2-1, $2}' OFS="\t" > "${barcode_output_dir}/${barcode}.coverage_mask.bed"

  # Combine Mask Files
  cat "${barcode_output_dir}/${barcode}.variant_mask.bed" "${barcode_output_dir}/${barcode}.coverage_mask.bed" | \
  sort -k1,1 -k2,2n | uniq > "${barcode_output_dir}/${barcode}.mask.bed"

  # Generate Consensus with Masking (using bcftools_env)
  conda run -n bcftools_env bcftools consensus -f "${REF_GENOME}" \
                     -m "${barcode_output_dir}/${barcode}.mask.bed" \
                     -o "${barcode_output_dir}/${barcode}.final.consensus.fasta" \
                     "${barcode_output_dir}/clair3_output/merge_output.vcf.gz"
  if [[ $? -ne 0 ]]; then
    echo "Error: bcftools consensus failed for $barcode" >> "${barcode_output_dir}/error.log"
    continue
  fi

  echo "Processing for $barcode completed."
done

echo "All samples processed successfully."
