### Haplotyping pipeline for HLA given minION fastq reads ###


#  ------------------------------------ MODULE ACTIVATION -------------------------
#module load samtools
#module load minimap2
#module load flye
#----------------------------------------------------------------------------------


minimap=FALSE
samtools=FALSE
pmd=TRUE
flye=FALSE


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
  prefix=paste0(prefix,'_')
  HLA.bam=paste0(prefix,"HLA.bam")
  system(paste0("samtools view -bS",SORTED_BAM," chr6:29940532-33099696 >",HLA.bam)) #select only the HLA genes
}

if(pmd){
  #from now pepper-margin-deepvariant
 
  # The pull command creates pepper_deepvariant_r0.4.sif file locally
  
  #system(paste0("singularity pull docker://kishwars/pepper_deepvariant:r0.4"))
  

  
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
  
  #class1
    HLA_A.bam=paste0(prefix,"HLA_A.bam")
    HLA_B.bam=paste0(prefix,"HLA_B.bam")
    HLA_C.bam=paste0(prefix,"HLA_C.bam")
  #class2
    HLA_DPA1.bam=paste0(prefix,"HLA_DPA1.bam")
    HLA_DPB1.bam=paste0(prefix,"HLA_DPB1.bam")
    HLA_DQA1.bam=paste0(prefix,"HLA_DQA1.bam")
    HLA_DQB1.bam=paste0(prefix,"HLA_DQB1.bam")
    HLA_DRB1.bam=paste0(prefix,"HLA_DRB1.bam")
    
  #here I create subset of the haplotagged bam for each gene
  
  #class1
    system(paste0("samtools view -bS", HAPLOTAGGED.bam," chr6:29941532-29946870 >",HLA_A.bam))
    system(paste0("samtools view -bS", HAPLOTAGGED.bam," chr6:31352875-31358179 >",HLA_B.bam))
    system(paste0("samtools view -bS", HAPLOTAGGED.bam," chr6:31267749-31273092 >",HLA_C.bam))
    
  #class2
    system(paste0("samtools view -bS", HAPLOTAGGED.bam," chr6:33064569-33080748 >",HLA_DPA1.bam))
    system(paste0("samtools view -bS", HAPLOTAGGED.bam," chr6:33075990-33089696 >",HLA_DPB1.bam))
    system(paste0("samtools view -bS", HAPLOTAGGED.bam," chr6:32637406-32654846 >",HLA_DQA1.bam))
    system(paste0("samtools view -bS", HAPLOTAGGED.bam," chr6:32659467-32666657 >",HLA_DQB1.bam))
    system(paste0("samtools view -bS", HAPLOTAGGED.bam," chr6:32578775-32589848 >",HLA_DRB1.bam))
  
  #then I execute flye for each gene
  system(paste0("echo ---------------- Executing Flye [4/4] ---------------- "))
  
  #split the haplotypes
  
  #class1
    system(paste0("bamtools split -in ", HLA_A.bam, " -tag HP"))
    system(paste0("bamtools split -in ", HLA_B.bam, " -tag HP"))
    system(paste0("bamtools split -in ", HLA_C.bam, " -tag HP"))
    
  #class2
    system(paste0("bamtools split -in ", HLA_DPA1.bam, " -tag HP"))
    system(paste0("bamtools split -in ", HLA_DPB1.bam, " -tag HP"))
    system(paste0("bamtools split -in ", HLA_DQA1.bam, " -tag HP"))
    system(paste0("bamtools split -in ", HLA_DQB1.bam, " -tag HP"))
    system(paste0("bamtools split -in ", HLA_DRB1.bam, " -tag HP"))
  
  #extract prefixes
  
  #class1
    flye_prefixA=sub("\\..*", "", HLA_A.bam)
    flye_prefixB=sub("\\..*", "", HLA_B.bam)
    flye_prefixC=sub("\\..*", "", HLA_C.bam)
  #class2
    flye_prefixDPA1=sub("\\..*", "", HLA_DPA1.bam)
    flye_prefixDPB1=sub("\\..*", "", HLA_DPB1.bam)
    flye_prefixDQA1=sub("\\..*", "", HLA_DQA1.bam)
    flye_prefixDQB1=sub("\\..*", "", HLA_DQB1.bam)
    flye_prefixDRB1=sub("\\..*", "", HLA_DRB1.bam)
  
  #names of .fa files
  
  #class1
    A1=paste0(prefix,'A1.fa')
    A2=paste0(prefix,'A2.fa')
    B1=paste0(prefix,'B1.fa')
    B2=paste0(prefix,'B2.fa')
    C1=paste0(prefix,'C1.fa')
    C2=paste0(prefix,'C2.fa')
  #class2
    DPA1_1=paste0(prefix,'DPA1_1.fa')
    DPA1_2=paste0(prefix,'DPA1_2.fa')
    DPB1_1=paste0(prefix,'DPB1_1.fa')
    DPB2_2=paste0(prefix,'DPB1_2.fa')
    DQA1_1=paste0(prefix,'DQA1_1.fa')
    DQA1_2=paste0(prefix,'DQA1_2.fa')
    DQB1_1=paste0(prefix,'DQB1_1.fa')
    DQB1_2=paste0(prefix,'DQB1_2.fa')
    DRB1_1=paste0(prefix,'DRB1_1.fa')
    DRB1_2=paste0(prefix,'DRB1_2.fa')
  
  #convert each haplotype from .bam to .fa
    
  #class1
    system(paste0("samtools bam2fq ", flye_prefixA, "TAG_HP_1.bam | seqtk seq -A > ", A1))
    system(paste0("samtools bam2fq ", flye_prefixA, "TAG_HP_2.bam | seqtk seq -A > ", A2))
    system(paste0("samtools bam2fq ", flye_prefixB, "TAG_HP_1.bam | seqtk seq -A > ", B1))
    system(paste0("samtools bam2fq ", flye_prefixB, "TAG_HP_2.bam | seqtk seq -A > ", B2))
    system(paste0("samtools bam2fq ", flye_prefixC, "TAG_HP_1.bam | seqtk seq -A > ", C1))
    system(paste0("samtools bam2fq ", flye_prefixC, "TAG_HP_2.bam | seqtk seq -A > ", C2))
  #class2
    system(paste0("samtools bam2fq ", flye_prefixDPA1_1, "TAG_HP_1.bam | seqtk seq -A > ", DPA1_1))
    system(paste0("samtools bam2fq ", flye_prefixDPA1_2, "TAG_HP_2.bam | seqtk seq -A > ", DPA1_2))
    system(paste0("samtools bam2fq ", flye_prefixDPB1_1, "TAG_HP_1.bam | seqtk seq -A > ", DPB1_1))
    system(paste0("samtools bam2fq ", flye_prefixDPB1_2, "TAG_HP_2.bam | seqtk seq -A > ", DPB1_2))
    system(paste0("samtools bam2fq ", flye_prefixDQA1_1, "TAG_HP_1.bam | seqtk seq -A > ", DQA1_1))
    system(paste0("samtools bam2fq ", flye_prefixDQA1_2, "TAG_HP_2.bam | seqtk seq -A > ", DQA1_2))
    system(paste0("samtools bam2fq ", flye_prefixDQB1_1, "TAG_HP_1.bam | seqtk seq -A > ", DQB1_1))
    system(paste0("samtools bam2fq ", flye_prefixDQB1_2, "TAG_HP_2.bam | seqtk seq -A > ", DQB1_2))
    system(paste0("samtools bam2fq ", flye_prefixDRB1_1, "TAG_HP_1.bam | seqtk seq -A > ", DRB1_1))
    system(paste0("samtools bam2fq ", flye_prefixDRB1_2, "TAG_HP_2.bam | seqtk seq -A > ", DRB1_2))
}

if(FALSE){  
  
  #out dirs
  
  #class1
    oA1=paste0(OUTPUT_DIR,"/flyeA1/")
    oA2=paste0(OUTPUT_DIR,"/flyeA2/")
    oB1=paste0(OUTPUT_DIR,"/flyeB1/")
    oB2=paste0(OUTPUT_DIR,"/flyeB2/")
    oC1=paste0(OUTPUT_DIR,"/flyeC1/")
    oC2=paste0(OUTPUT_DIR,"/flyeC2/")
  #class2
    oDPA1_1=paste0(OUTPUT_DIR,"/flyeDPA1_1/")
    oDPA1_2=paste0(OUTPUT_DIR,"/flyeDPA1_2/")
    oDPB1_1=paste0(OUTPUT_DIR,"/flyeDPB1_1/")
    oDPB1_2=paste0(OUTPUT_DIR,"/flyeDPB1_2/")
    oDQA1_1=paste0(OUTPUT_DIR,"/flyeDQA1_1/")
    oDQA1_2=paste0(OUTPUT_DIR,"/flyeDQA1_2/")
    oDQB1_1=paste0(OUTPUT_DIR,"/flyeDQB1_1/")
    oDQB1_2=paste0(OUTPUT_DIR,"/flyeDQB1_2/")
    oDRB1_1=paste0(OUTPUT_DIR,"/flyeDRB1_1/")
    oDRB1_2=paste0(OUTPUT_DIR,"/flyeDRB1_2/")
  
  #execute de-novo assembly with flye
  
  #class1
    system(paste0("flye --nano-raw ", A1, " --out-dir ", oA1, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", A2, " --out-dir ", oA2, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", B1, " --out-dir ", oB1, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", B2, " --out-dir ", oB2, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", C1, " --out-dir ", oC1, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", C2, " --out-dir ", oC2, "--threads 4 -m 1000 -i 2"))
  #class2
    system(paste0("flye --nano-raw ", DPA1_1, " --out-dir ", oDPA1_1, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DPA1_2, " --out-dir ", oDPA1_2, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DPB1_1, " --out-dir ", oDPB1_1, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DPB1_2, " --out-dir ", oDPB1_2, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DQA1_1, " --out-dir ", oDQA1_1, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DQA1_2, " --out-dir ", oDQA1_2, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DQB1_1, " --out-dir ", oDQB1_1, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DQB1_2, " --out-dir ", oDQB1_2, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DRB1_1, " --out-dir ", oDRB1_1, "--threads 4 -m 1000 -i 2"))
    system(paste0("flye --nano-raw ", DRB1_2, " --out-dir ", oDRB1_2, "--threads 4 -m 1000 -i 2"))
}

system(paste0("echo ---------------- Finished ---------------- "))
