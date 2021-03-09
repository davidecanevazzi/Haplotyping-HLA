
#This is the code for the aplotyping pipeline

reference="GRCh38_latest_genomic.fna.gz"
#system(paste0("cd ", reference))

aln.sam="aln2.sam"   #name of the .sam file
aln.bam="aln2.bam"   #name of the .bam file
aln.sorted.bam="aln2.sorted.bam" #name of the sorted bam file"
HLA="HLA_A_B_ashke.bam" #name of the file in the region of interest
region="NC_000006.12:29942532-31353875" #region of interest

system(paste0("samtools view -bS ",aln.sam, " > ",aln.bam))  #convert .sam>.bam
system(paste0("sort ",aln.bam," -o ", aln.sorted.bam)) #sort the .bam file
system(paste0("samtools index ",aln.sorted.bam)) #index the sorted .bam file
system(paste0("samtools view -b ", aln.sorted.bam, " ", region, " > ", HLA )) #selects only the region of interest
system(paste0("samtools index ",HLA)) #index of the selected reads for IGV visualization


#samtools view -bS aln1.sam > aln1.bam
#samtools sort aln1.bam -o aln1.sorted.bam
#samtools index sample.sorted.bam

#samtools view aln1.sorted.bam "NT_167249.2:1,181,883-31,357,583" > HLA_A_B.bam
#samtools view -b aln1.sorted.bam "NT_167249.2:1,181,883-31,357,583" > HLA_A_B1.bam

#samtools view -b aln1.sorted.bam "chr6:29,908,309-31,326,956" > HLA_A_B.bam
#samtools index HLA_A_B.bam 


#system(paste0("cp ", home_dir, "/", sample_name, ".fasta ", sample_dir, "/", sample_name, "_decont.fasta"))

