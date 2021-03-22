BASE="${HOME}/ont-case-study"

# Set up input data
INPUT_DIR="${BASE}/input/data"
REF="GRCh38_no_alt.chr20.fa"
BAM="HG002_ONT_50x_2_GRCh38.chr20.bam"

# Set the number of CPUs to use
THREADS="4"

# Set up output directory
OUTPUT_DIR="${BASE}/output"
OUTPUT_PREFIX="HG002_ONT_50x_2_GRCh38_PEPPER_Margin_DeepVariant.chr20"
OUTPUT_VCF="HG002_ONT_50x_2_GRCh38_PEPPER_Margin_DeepVariant.chr20.vcf.gz"

## Create local directory structure
#mkdir -p "${OUTPUT_DIR}"
#mkdir -p "${INPUT_DIR}"

# Download the data to input directory
#wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_ONT_50x_2_GRCh38.chr20.bam
#wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_ONT_50x_2_GRCh38.chr20.bam.bai
#wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr20.fa
#wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr20.fa.fai
## Pull the docker image to sigularity, this is a 6.6GB download
singularity pull docker://kishwars/pepper_deepvariant:r0.4

# The pull command creates pepper_deepvariant_r0.4.sif file locally

# Run PEPPER-Margin-DeepVariant
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t ${THREADS} \
--ont
# Optional parameters:
# -s HG002 # optional: Sets Sample Name
# --gvcf # optional: Produces gVCF output
# --phased_output # optional: Produces phased output
