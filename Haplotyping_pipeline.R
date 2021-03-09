
#This is the code for the aplotyping pipeline

reference="GRCh38_latest_genomic.fna.gz"
#system(paste0("cd ", reference))

aln.sam="aln2.sam"
aln.bam="aln2.bam"
aln.sorted.bam="aln2.sorted.bam"
HLA="HLA_A_B_ashke.bam"
region="NC_000006.12:29942532-31353875"

system(paste0("samtools view -bS ",aln.sam, " > ",aln.bam))
system(paste0("sort ",aln.bam," -o ", aln.sorted.bam))
system(paste0("samtools index ",aln.sorted.bam))
system(paste0("samtools view -b ", aln.sorted.bam, " ", region, " > ", HLA ))
system(paste0("samtools index ",HLA))


#samtools view -bS aln1.sam > aln1.bam
#samtools sort aln1.bam -o aln1.sorted.bam
#samtools index sample.sorted.bam

#samtools view aln1.sorted.bam "NT_167249.2:1,181,883-31,357,583" > HLA_A_B.bam
#samtools view -b aln1.sorted.bam "NT_167249.2:1,181,883-31,357,583" > HLA_A_B1.bam

#samtools view -b aln1.sorted.bam "chr6:29,908,309-31,326,956" > HLA_A_B.bam
#samtools index HLA_A_B.bam 


#system(paste0("cp ", home_dir, "/", sample_name, ".fasta ", sample_dir, "/", sample_name, "_decont.fasta"))
