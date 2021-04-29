### Haplotyping pipeline for HLA given minION fastq reads ###


#  ------------------------------------ MODULE ACTIVATION -------------------------
#module load samtools
#module load minimap2
#module load flye
#----------------------------------------------------------------------------------


minimap=FALSE
samtools=FALSE
pmd=FALSE
flye=TRUE


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
  arg[2]="output"
}

OUTPUT_DIR = arg[2]
VCF = paste0(prefix,".vcf")


# Set the number of CPUs to use
THREADS="12"


paste0('fastq file --> ', FASTQ)
paste0('prefix --> ', prefix)
paste0('SAM file --> ', SAM)
paste0('BAM file --> ', BAM)
paste0('output dir --> ', OUTPUT_DIR)
paste0('VCF file --> ',VCF)

if(minimap){
  #mapping against the reference
  
  system(paste0("echo ---------------- Mapping with minimap2 [1/4] ---------------- "))
  
  system(paste0("minimap2 -a -z 600,200 -x map-ont ", REFERENCE , " ", FASTQ, " >",SAM)) 
}

if(samtools){
  system(paste0("echo ---------------- Samtools indexing and sorting [2/4] ---------------- "))
  #data conversion and indexinx with samtools
  system(paste0("samtools view -bS ",SAM, " > ",BAM))  #convert .sam>.bam
  system(paste0("rm ",SAM))
  system(paste0("samtools sort ",BAM," -o ", SORTED_BAM)) #sort the .bam file
  system(paste0("samtools index ",SORTED_BAM)) #index the sorted .bam file
}

if(pmd){
  #from now pepper-margin-deepvariant
  prefix=paste0(prefix,'_')
  HLA.bam=paste0(prefix,"HLA.bam")
  
  # The pull command creates pepper_deepvariant_r0.4.sif file locally
  
  #system(paste0("singularity pull docker://kishwars/pepper_deepvariant:r0.4"))
  
  system(paste0("samtools view -bS",SORTED_BAM," chr6:29940532-31355875 >",HLA.bam)) #select only the HLA genes
  
  system(paste0("echo ---------------- Executing Pepper-Margin-Deepvariant [3/4] ---------------- "))
  
  system(paste0("singularity exec --bind /usr/lib/locale/ \
              pepper_deepvariant_r0.4.sif \
              run_pepper_margin_deepvariant call_variant \
              -b ",HLA.bam , " \
              --phased_output \
              -f ", REFERENCE," \
              -o ", OUTPUT_DIR, " \
              -p ", prefix, " \
              -t ${",THREADS,"} \
              --ont"))

}
  #From here haplotyping with de-novo assembly with flye

if(flye){  
  HAPLOTAGGED.bam=list.files(paste0(OUTPUT_DIR, "/intermediate_files"), pattern="*MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam", full.names=TRUE)[1]

  HLA_A.bam=paste0(prefix,"_HLA_A.bam")
  HLA_B.bam=paste0(prefix,"_HLA_B.bam")
  HLA_C.bam=paste0(prefix,"_HLA_C.bam")
  
  #here I create subset of the haplotagged bam for each gene
  system(paste0("samtools view -bS ", HAPLOTAGGED.bam," chr6:29941532-29946870 >",HLA_A.bam))
  system(paste0("samtools view -bS ", HAPLOTAGGED.bam," chr6:31352875-31358179 >",HLA_B.bam))
  system(paste0("samtools view -bS ", HAPLOTAGGED.bam," chr6:31267749-31273092 >",HLA_C.bam))

  #31353872-31357188 
  #mica
  
  
  #then I execute flye for each gene
  system(paste0("echo ---------------- Executing Flye [4/4] ---------------- "))
  
  #split the haplotypes
  system(paste0("bamtools split -in ", HLA_A.bam, " -tag HP"))
  system(paste0("bamtools split -in ", HLA_B.bam, " -tag HP"))
  system(paste0("bamtools split -in ", HLA_C.bam, " -tag HP"))
  
  #extract prefixes
  flye_prefixA=sub("\\..*", "", HLA_A.bam)
  flye_prefixB=sub("\\..*", "", HLA_B.bam)
  flye_prefixC=sub("\\..*", "", HLA_C.bam)
  
  #names of .fa files
  A1=paste0(prefix,'A1.fa')
  A2=paste0(prefix,'A2.fa')
  B1=paste0(prefix,'B1.fa')
  B2=paste0(prefix,'B2.fa')
  C1=paste0(prefix,'C1.fa')
  C2=paste0(prefix,'C2.fa')
  
  #convert each haplotype from .bam to .fa
  system(paste0("samtools bam2fq ", flye_prefixA, ".TAG_HP_1.bam | seqtk seq -A > ", A1))
  system(paste0("samtools bam2fq ", flye_prefixA, ".TAG_HP_2.bam | seqtk seq -A > ", A2))
  system(paste0("samtools bam2fq ", flye_prefixB, ".TAG_HP_1.bam | seqtk seq -A > ", B1))
  system(paste0("samtools bam2fq ", flye_prefixB, ".TAG_HP_2.bam | seqtk seq -A > ", B2))
  system(paste0("samtools bam2fq ", flye_prefixC, ".TAG_HP_1.bam | seqtk seq -A > ", C1))
  system(paste0("samtools bam2fq ", flye_prefixC, ".TAG_HP_2.bam | seqtk seq -A > ", C2))
}
  
if(FALSE){  
    
  #out dirs
  oA1=paste0(OUTPUT_DIR,"/flyeA1/")
  oA2=paste0(OUTPUT_DIR,"/flyeA2/")
  oB1=paste0(OUTPUT_DIR,"/flyeB1/")
  oB2=paste0(OUTPUT_DIR,"/flyeB2/")
  oC1=paste0(OUTPUT_DIR,"/flyeC1/")
  oC2=paste0(OUTPUT_DIR,"/flyeC2/")
  
  #execute de-novo assembly with flye
  system(paste0("flye --nano-raw ", A1, " --out-dir ", oA1, "--threads 4 -m 1000"))
  system(paste0("flye --nano-raw ", A2, " --out-dir ", oA2, "--threads 4 -m 1000"))
  system(paste0("flye --nano-raw ", B1, " --out-dir ", oB1, "--threads 4 -m 1000"))
  system(paste0("flye --nano-raw ", B2, " --out-dir ", oB2, "--threads 4 -m 1000"))
  system(paste0("flye --nano-raw ", C1, " --out-dir ", oC1, "--threads 4 -m 1000"))
  system(paste0("flye --nano-raw ", C2, " --out-dir ", oC2, "--threads 4 -m 1000"))
  
}

system(paste0("echo ---------------- Finished ---------------- "))


