#!/bin/bash                              
set -e

ml Singularity/3.5.3                     
export PATH=$SINGULARITYROOT/bin/:$PATH

BASE_BUCKET="/fh/scratch/delete90/nelson_p/james/"

# Load the module                                                                                                                                 
ml nextflow

NXF_VER=20.01.0 nextflow \
    -c ../configs/nextflow.config \
    run \
    -resume \
    main.nf \
    -profile hpc \
    -work-dir $BASE_BUCKET/pdx/matched/work3/ \
    --pdx true \
    --matched true \
    --cnvnator false \
    --gatk_single_var false \
    --input_csv ../manifests/lucap_manifest.txt \
    --input_vcf ./vcfs.txt \
    --output_folder /fh/scratch/delete90/nelson_p/james/funcotating/cnvkit2/\
    --pdx_reference $BASE_BUCKET/references/mm10/GRCm38.primary_assembly.genome.fa \
    --reference $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa \
    --reference_index $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa.fai \
    --reference_dict $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.dict \
    --rear $BASE_BUCKET/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --rear_index $BASE_BUCKET/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi \
    --indels $BASE_BUCKET/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --indels_index $BASE_BUCKET/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
    --common_variants $BASE_BUCKET/references/hg38/af-only-gnomad.hg38.vcf.gz \
    --common_variants_index $BASE_BUCKET/references/hg38/af-only-gnomad.hg38.vcf.gz.tbi \
    --ref_name hg38 \
    --contig_dict $BASE_BUCKET/references/gatk/Homo_sapiens_assembly38.dict \
    --input_beds local_capture.csv \
    --clinvar /fh/scratch/delete90/nelson_p/james/references/hg38/clinvar_20200520_chr.vcf \
    --refFlat /fh/scratch/delete90/nelson_p/james/references/hg38/refFlat.txt \ 
    --coords /home/jehopkin/nf-fh-pcp-wes/coords.txt \


