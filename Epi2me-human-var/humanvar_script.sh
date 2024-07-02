#basecaller might give some problems, use tools to check ur merged bam file to see what basecaller was used, if it has khz prefix like mine did just skip the prefix. At some point i needed to fix chromsome naming for seomthing here by adding chr (dont remember what sorry...). The people in the github for the pipeline are super helpfull dont be afraid to ask for help. i also forgot to add CNV flag, probably good information. (worked with v2.2.4) (Nextflow/23.10.0)

# Run the workflow 
nextflow run epi2me-labs/wf-human-variation \
    --bam '/path/to/merged.bam' \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0' \
    --mod True \
    --ref '/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa' \
    --sample_name 'sample_name' \
    --snp True \
    --sv True \
    --phased True \
    --sex XX \
    -profile singularity \
    --outdir '/path/to/output_directory' \
    --tr_bed '/path/to/human_GRCh38_no_alt_analysis_set.trf.bed' 
