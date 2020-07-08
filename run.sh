#!/bin/bash                              
set -e

ml Singularity/3.5.3                     
export PATH=$SINGULARITYROOT/bin/:$PATH

BASE_BUCKET="/fh/scratch/delete90/nelson_p/james/"

# Load the module                                                                                                                                 
ml nextflow

NXF_VER=20.01.0 nextflow \
    run \
    -resume \
    main.nf \
    -profile hpc \
    -work-dir $BASE_BUCKET/pdx/matched/work3/ \
    --input_csv ../manifests/lucap_manifest.csv \
    --output_folder /fh/scratch/delete90/nelson_p/james/cnvkit/\
    --reference $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa \
    --reference_index $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa.fai \
    --reference_dict $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.dict \
    --input_beds local_capture.csv \
    --refFlat /fh/scratch/delete90/nelson_p/james/references/hg38/refFlat.txt \ 


