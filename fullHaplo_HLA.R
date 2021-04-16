### Haplotyping pipeline for HLA given minION fastq reads ###

arg=commandArgs(trailingOnly = TRUE)
  
  FASTQ=arg[1]

  prefix=sub("\\..*", "", FASTQ)

  #you have to discover how to move around 

  # Set up input data
  REFERENCE= "chr6.fa"
  SAM = paste0(prefix,".sam")
  BAM = paste0(prefix,".bam")
  SORTED_BAM = paste0(prefix,".sorted.bam")
  
  if (is.na(arg[2])) {
    arg[2]="pmd_output"
  }
  OUTPUT_DIR = arg[2]
  OUTPUT_VCF = paste0(prefix,".vcf")
  
  
  # Set the number of CPUs to use
  THREADS="4"

  
  FASTQ
  prefix
  REFERENCE
  SAM
  BAM
  SORTED_BAM
  OUTPUT_DIR
  OUTPUT_VCF
  
  system(paste0("echo ---------------- Mapping with minimap2 [1/4] ---------------- "))
  

if(FALSE){
  #mapping against the reference
  
  system(paste0("echo ---------------- Mapping with minimap2 [1/4] ---------------- "))
  
  system(paste0("minimap2 -a -z 600,200 -x map-ont ", REFERENCE , " ", FASTQ, " >",SAM)) 
  
  system(paste0("echo ---------------- Samtools indexing and sorting [2/4] ---------------- "))
  #data conversion and indexinx with samtools
  system(paste0("samtools view -bS ",SAM, " > ",BAM))  #convert .sam>.bam
  system(paste0("rm ",SAM))
  system(paste0("samtools sort ",BAM," -o ", SORTED_BAM)) #sort the .bam file
  system(paste0("samtools index ",SORTED_BAM)) #index the sorted .bam file
  
  #from now pepper-margin-deepvariant
  
  HLA.bam=paste0(prefix,"HLA.bam")
  
  # The pull command creates pepper_deepvariant_r0.4.sif file locally
  
  #system(paste0("singularity pull docker://kishwars/pepper_deepvariant:r0.4"))
  
  system(paste0("samtools view",SORTED_BAM,"chr6:29940532-31355875 >",HLA.bam)) #select only the HLA genes
  
  system(paste0("echo ---------------- Executing Pepper-Margin-Deepvariant [3/4] ---------------- "))
  
  system(paste0("singularity exec --bind /usr/lib/locale/ \
              pepper_deepvariant_r0.4.sif \
              run_pepper_margin_deepvariant call_variant \
              -b ",HLA.bam , " \
              --phased_output \
              -f ", REFERENCE," \
              -o ", OUTPUT_DIR, " \
              -p ", prefix, " \
              -t ${THREADS} \
              --ont"))
  
  #From here haplotyping with de-novo assembly with flye
  
  HAPLOTAGGED.bam=list.files(paste0(OUTPUT_DIR, "/intermediate_files"), pattern="*MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam", full.names=TRUE)[1]
  
  HLA_A.bam=paste0(prefix,"HLA_A.bam")
  HLA_B.bam=paste0(prefix,"HLA_B.bam")
  HLA_C.bam=paste0(prefix,"HLA_C.bam")
  
  
  
  #here I create subset of the haplotagged bam for each gene
  system(paste0("samtools view", HAPLOTAGGED.bam,"chr6:29941532-29946870 >",HLA_A.bam))
  system(paste0("samtools view",HAPLOTAGGED.bam,"chr6:31352875-31358179 >",HLA_B.bam))
  system(paste0("samtools view",HAPLOTAGGED.bam,"chr6:31267749-31273092  >",HLA_C.bam))
  
  
  
  system(paste0("echo ---------------- Executing Flye [4/4] ---------------- "))
  
  
  
  
  
  
}
  
system(paste0("echo ---------------- Finished ---------------- "))

