#This creates a merged BAM file for human variation pipeline and a BED file for methylation analysis. Possible errors here are headers and chromosome labeling easy to fix. This also creates other files usefull for various different downstream applications which can be deleted if not required (takes allot of space)

# Log file
log_file="script_output.log"

# Redirect all output to the log file
exec > >(tee -a $log_file) 2>&1

# Increase the limit of open files
ulimit -n 5000

# Directories containing BAM files
input_dir="path/to/input_directory"

# Output directory
output_dir="path/to/output_directory"

# Reference genome
reference_genome="path/to/reference_genome.fa"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Merge BAM files
merged_bam="$output_dir/merged.bam"
echo "Merging BAM files from $input_dir into $merged_bam"
samtools merge -n -c -@ 32 -o $merged_bam $(ls $input_dir/*.bam | grep -v ".bam.bai")

# Convert merged BAM to FASTQ
fastq_file="$output_dir/output.fq"
echo "Converting $merged_bam to FASTQ $fastq_file"
samtools fastq -TMM,ML --threads 32 $merged_bam > $fastq_file

# Align FASTQ to reference genome
aligned_bam="$output_dir/output_aligned.bam"
echo "Aligning $fastq_file to reference genome $reference_genome"
minimap2 -k17 -ax map-ont --secondary=yes -t 32 -y $reference_genome $fastq_file | samtools view -Sb - --threads 32 | samtools sort - --threads 32 > $aligned_bam

# Filtering data for mapping quality
filtered_bam="$output_dir/output_aligned_and_MQ_filtered.bam"
echo "Filtering $aligned_bam for mapping quality"
samtools view -h -b -q 30 --threads 32 $aligned_bam > $filtered_bam

# Filtering data to get rid of signal noise
filtered_bam_2="$output_dir/output_aligned_and_filtered_2.bam"
echo "Filtering $filtered_bam to get rid of signal noise"
samtools view -@ 32 -h $filtered_bam | awk 'length($10) >= 400 || $1 ~ /^@/' | samtools view -@ 32 -b -o $filtered_bam_2 -

# Index the filtered BAM file
echo "Indexing $filtered_bam"
samtools index $filtered_bam

# Convert BAM to BED
bed_file="$output_dir/$(basename ${output_dir}).bed"
echo "Converting $filtered_bam to BED $bed_file"
modbam2bed -e -t 32 $reference_genome $filtered_bam > $bed_file
