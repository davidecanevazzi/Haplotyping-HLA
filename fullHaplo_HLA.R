### Haplotyping pipeline for HLA given minION fastq reads ###



# Set up input data
REFERENCE= "chr6.fa"
prefix = "name_of_the _sample"
FASTQ = "name_of_fastq_file"
SAM = paste0(prefix,".sam")
BAM = paste0(prefix,".bam")
SORTED_BAM = paste0(prefix,".sorted.bam")
INPUT_DIR = "/input/data"
OUTPUT_DIR = "/output"
OUTPUT_VCF = paste0(prefix,".vcf")


# Set the number of CPUs to use
THREADS="4"


#mapping against the reference
system(paste0("minimap2 -a -z 600,200 -x map-ont ", REFERENCE , " ", FASTQ, " >",SAM)) 

#data conversion and indexinx with samtools
system(paste0("samtools view -bS ",SAM, " > ",BAM))  #convert .sam>.bam
system(paste0("rm ",SAM))
system(paste0("samtools sort ",BAM," -o ", SORTED_BAM)) #sort the .bam file
system(paste0("samtools index ",SORTED_BAM)) #index the sorted .bam file

#from now pepper-margin-deepvariant



# The pull command creates pepper_deepvariant_r0.4.sif file locally

system(paste0("singularity pull docker://kishwars/pepper_deepvariant:r0.4"))

system(paste0("singularity exec --bind /usr/lib/locale/ \
              pepper_deepvariant_r0.4.sif \
              run_pepper_margin_deepvariant call_variant \
              -b ",BAM, " \
              -f ", REFERENCE," \
              -o ", OUTPUT_DIR, " \
              -p ", prefix, " \
              -t ${THREADS} \
              --ont"))



