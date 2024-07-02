#super simple, swap out the aligners according to the nf-core/rnaseq github. Really stupid that bootstrapping for aligners is not automaticly included i sadly made the mistake of not bootstrapping and didnt realise.... check the out extra flags that can be added to aligners to bootstrap (this is needed to do analysis with sleuth)


mkdir rnaseq_salmon_out
ml use /path/to/module/directory/
ml load Nextflow/23.10.0

nextflow run nf-core/rnaseq \
  --input samplesheet.csv \
  -profile singularity \
  --pseudo_aligner salmon \
  --fasta /path/to/genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
  --gtf /path/to/genome/Homo_sapiens.GRCh38.112.gtf.gz \
  --save_reference \
  --outdir rnaseq_salmon_out
