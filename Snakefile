#how it works:
#module load samtools minimap2 PEPPER DEEPVARIANT MARGIN bamtools tabix
#module load gcc/6.3.0 mkl/11.3.3 PYTHON/3.6.0 GLIBC/2.19 CONDA/4.5.11_PYTHON3
#source activate /home/devel/talioto/miniconda3/envs/flye-v2.8.3/


#mapping and indexing



rule minimap:
    input:
        "data/chr6.fa",
        "data/{sample}.fastq"
    output:
        "data/{sample}.sam"
    shell:
        "minimap2 -a -t 24 -z 600,200 -x map-ont {input} > {output}"

rule samtools_view:
    input:
        "data/{sample}.sam"
    output:
        "data/{sample}.bam"
    shell:
        "samtools view -@12 -bS {input} > {output}"

rule samtools_sort:
    input:
        "data/{sample}.bam"
    output:
        "data/{sample}.sorted.bam"
    shell:
        "santools sort -@ 12 {input} -o {output}"

rule samtools_index:
    input:
        "data/{sample}.bam"
    output:
        "data/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule HLA_selection:
    input:
        bam= "data/{sample}.bam",
        bai= "data/{sample}.bam.bai"
    output:
        "data/{sample}_HLA.bam"
    shell:
        "samtools view -bS {input.bam} chr6:29940532-33099696 > {output}"

#export HDF5_USE_FILE_LOCKING='FALSE'

#run pepper margin

rule pepper:
    input:
        bam="data/{sample}_HLA.bam" ,
        ref="data/chr6.fa"
    output:
        "output/logs/{sample}_pepper_snp.log"
    shell:
        "singularity exec /apps/PEPPER/0.4/pepper pepper_snp  call_variant -b  {input.bam} -f  {input.ref} -t 4 -m ../../../../../scratch/production/DAT/apps/PEPPER_MARGIN_DEEPVARIANT/0.4/pepper_models/PEPPER_SNP_R941_ONT_V4.pkl -o output/pepper_snp -r chr6 -s Sample -w 4 -bs 64 --ont 2>&1|tee {output}"

rule gzip:
    input:
        "output/logs/{sample}_pepper_snp.log",
    output:
        "output/pepper_snp/{sample}_log_gzip.txt"
    shell:
        "gzip output/pepper_snp/PEPPER_SNP_OUTPUT.vcf|tee {output}"

rule mv:
    input:
        "output/pepper_snp/{sample}_log_gzip.txt"
    output:
        "output/pepper_snp/{sample}_PEPPER_SNP_OUTPUT.vcf.gz"
    shell:
        "mv output/pepper_snp/PEPPER_SNP_OUTPUT.vcf.gz {output}"

rule tabix:
    input:
        "output/pepper_snp/{sample}_PEPPER_SNP_OUTPUT.vcf.gz"
    output:
        "output/pepper_snp/{sample}_PEPPER_SNP_OUTPUT.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

rule margin:
    input:
        bam="data/{sample}_HLA.bam",
        vcf="output/pepper_snp/{sample}_PEPPER_SNP_OUTPUT.vcf.gz",
        vcfi="output/pepper_snp/{sample}_PEPPER_SNP_OUTPUT.vcf.gz.tbi",
        ref="data/chr6.fa"
    output:
        "output/logs/{sample}_margin_haplotag.log"
    shell:
        "time margin phase {input.bam}  {input.ref}  {input.vcf} ../../../../../apps/MARGIN/2.2.2/params/misc/allParams.ont_haplotag.json -t 4 -r chr6 -V -o output/MARGIN_PHASED.PEPPER_SNP_MARGIN 2>&1 | tee {output}"


rule samtools_index2:
    input:
        "output/logs/{sample}_margin_haplotag.log"
    output:
        "output/logs/{sample}_index2.log"
    shell:
        "samtools index output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam| tee {output}"


#division into HLA-genes

rule select_A:
    input:
        bam="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        bai="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        lo="output/logs/{sample}_index2.log"
    output:
        "output/{sample}_HLA_A.bam"
    shell:
        "samtools view -bS {input.bam} chr6:29941532-29946870 > {output}"

rule select_B:
    input:
        bam="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        bai="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        lo="output/logs/{sample}_index2.log"
    output:
        "output/{sample}_HLA_B.bam"
    shell:
        "samtools view -bS {input.bam} chr6:31352875-31358179 > {output}"

rule select_C:
    input:
        bam="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        bai="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        lo="output/logs/{sample}_index2.log"
    output:
        "output/{sample}_HLA_C.bam"
    shell:
        "samtools view -bS {input.bam} chr6:31267749-31273092  > {output}"

rule select_DPA1:
    input:
        bam="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        bai="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        lo="output/logs/{sample}_index2.log"
    output:
        "output/{sample}_HLA_DPA1.bam"
    shell:
        "samtools view -bS {input.bam} chr6:33064569-33080748 > {output}"

rule select_DPB1:
    input:
        bam="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        bai="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        lo="output/logs/{sample}_index2.log"
    output:
        "output/{sample}_HLA_DPB1.bam"
    shell:
        "samtools view -bS {input.bam} chr6:33075990-33089696 > {output}"

rule select_DQA1:
    input:
        bam="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        bai="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        lo="output/logs/{sample}_index2.log"
    output:
        "output/{sample}_HLA_DQA1.bam"
    shell:
        "samtools view -bS {input.bam} chr6:32637406-32654846 > {output}"

rule select_DQB1:
    input:
        bam="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        bai="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        lo="output/logs/{sample}_index2.log"
    output:
        "output/{sample}_HLA_DQB1.bam"
    shell:
        "samtools view -bS {input.bam} chr6:32659467-32666657> {output}"

rule select_DRB1:
    input:
        bam="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        bai="output/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        lo="output/logs/{sample}_index2.log"
    output:
        "output/{sample}_HLA_DRB1.bam"
    shell:
        "samtools view -bS {input.bam} chr6:32578775-32589848 > {output}"


rule split_A:
    input:
        "output/{sample}_HLA_A.bam"
    output:
        "output/{sample}_HLA_A.TAG_HP_1.bam",
        "output/{sample}_HLA_A.TAG_HP_2.bam"
    shell:
        "bamtools split -in {input} -tag HP"

rule split_B:
    input:
        "output/{sample}_HLA_B.bam"
    output:
        "output/{sample}_HLA_B.TAG_HP_1.bam",
        "output/{sample}_HLA_B.TAG_HP_2.bam"
    shell:
        "bamtools split -in {input} -tag HP"

rule split_C:
    input:
        "output/{sample}_HLA_C.bam"
    output:
        "output/{sample}_HLA_C.TAG_HP_1.bam",
        "output/{sample}_HLA_C.TAG_HP_2.bam"
    shell:
        "bamtools split -in {input} -tag HP"

rule split_DPA1:
    input:
        "output/{sample}_HLA_DPA1.bam"
    output:
        "output/{sample}_HLA_DPA1.TAG_HP_1.bam",
        "output/{sample}_HLA_DPA1.TAG_HP_2.bam"
    shell:
        "bamtools split -in {input} -tag HP"

rule split_DPB1:
    input:
        "output/{sample}_HLA_DPB1.bam"
    output:
        "output/{sample}_HLA_DPB1.TAG_HP_1.bam",
        "output/{sample}_HLA_DPB1.TAG_HP_2.bam"
    shell:
        "bamtools split -in {input} -tag HP"

rule split_DQA1:
    input:
        "output/{sample}_HLA_DQA1.bam"
    output:
        "output/{sample}_HLA_DQA1.TAG_HP_1.bam",
        "output/{sample}_HLA_DQA1.TAG_HP_2.bam"
    shell:
        "bamtools split -in {input} -tag HP"

rule split_DQB1:
    input:
        "output/{sample}_HLA_DQB1.bam"
    output:
        "output/{sample}_HLA_DQB1.TAG_HP_1.bam",
        "output/{sample}_HLA_DQB1.TAG_HP_2.bam"
    shell:
        "bamtools split -in {input} -tag HP"

rule split_DRB1:
    input:
        "output/{sample}_HLA_DRB1.bam"
    output:
        "output/{sample}_HLA_DRB1.TAG_HP_1.bam",
        "output/{sample}_HLA_DRB1.TAG_HP_2.bam"
    shell:
        "bamtools split -in {input} -tag HP"


rule to_fasta_A1:
    input:
        "output/{sample}_HLA_A.TAG_HP_1.bam"
    output:
        "output/{sample}_HLA_A_1.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"


rule to_fasta_A2:
    input:
        "output/{sample}_HLA_A.TAG_HP_2.bam"
    output:
        "output/{sample}_HLA_A_2.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_B1:
    input:
        "output/{sample}_HLA_B.TAG_HP_1.bam"
    output:
        "output/{sample}_HLA_B_1.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_B2:
    input:
        "output/{sample}_HLA_B.TAG_HP_2.bam"
    output:
        "output/{sample}_HLA_B_2.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_C1:
    input:
        "output/{sample}_HLA_C.TAG_HP_1.bam"
    output:
        "output/{sample}_HLA_C_1.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_C2:
    input:
        "output/{sample}_HLA_C.TAG_HP_2.bam"
    output:
        "output/{sample}_HLA_C_2.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DPA1_1:
    input:
        "output/{sample}_HLA_DPA1.TAG_HP_1.bam"
    output:
        "output/{sample}_HLA_DPA1_1.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DPA1_2:
    input:
        "output/{sample}_HLA_DPA1.TAG_HP_2.bam"
    output:
        "output/{sample}_HLA_DPA1_2.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DPB1_1:
    input:
        "output/{sample}_HLA_DPB1.TAG_HP_1.bam"
    output:
        "output/{sample}_HLA_DPB1_1.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DPB1_2:
    input:
        "output/{sample}_HLA_DPB1.TAG_HP_2.bam"
    output:
        "output/{sample}_HLA_DPB1_2.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DQA1_1:
    input:
        "output/{sample}_HLA_DQA1.TAG_HP_1.bam"
    output:
        "output/{sample}_HLA_DQA1_1.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DQA1_2:
    input:
        "output/{sample}_HLA_DQA1.TAG_HP_2.bam"
    output:
        "output/{sample}_HLA_DQA1_2.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DQB1_1:
    input:
        "output/{sample}_HLA_DQB1.TAG_HP_1.bam"
    output:
        "output/{sample}_HLA_DQB1_1.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DQB1_2:
    input:
        "output/{sample}_HLA_DQB1.TAG_HP_2.bam"
    output:
        "output/{sample}_HLA_DQB1_2.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DRB1_1:
    input:
        "output/{sample}_HLA_DRB1.TAG_HP_1.bam"
    output:
        "output/{sample}_HLA_DRB1_1.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"

rule to_fasta_DRB1_2:
    input:
        "output/{sample}_HLA_DRB1.TAG_HP_2.bam"
    output:
        "output/{sample}_HLA_DRB1_2.fa"
    shell:
        "samtools fasta  {input}  -0 {output}"


##### assembly with flye


rule flye_A_1:
    input:
        "output/{sample}_HLA_A_1.fa"
    output:
        "output/flye_A_1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_A_1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_A_2:
    input:
        "output/{sample}_HLA_A_2.fa"
    output:
        "output/flye_A_2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_A_2/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_B_1:
    input:
        "output/{sample}_HLA_B_1.fa"
    output:
        "output/flye_B_1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_B_1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_B_2:
    input:
        "output/{sample}_HLA_B_2.fa"
    output:
        "output/flye_B_2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_B_2/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_C_1:
    input:
        "output/{sample}_HLA_C_1.fa"
    output:
        "output/flye_C_1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_C_1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_C_2:
    input:
        "output/{sample}_HLA_C_2.fa"
    output:
        "output/flye_C_2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_C_2/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DPA_1:
    input:
        "output/{sample}_HLA_DPA1_1.fa"
    output:
        "output/flye_DPA1_1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DPA1_1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DPA1_2:
    input:
        "output/{sample}_HLA_DPA1_2.fa"
    output:
        "output/flye_DPA1_2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DPA1_2/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DPB1_1:
    input:
        "output/{sample}_HLA_DPB1_1.fa"
    output:
        "output/flye_DPB1_1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DPB1_1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DPB1_2:
    input:
        "output/{sample}_HLA_DPB1_2.fa"
    output:
        "output/flye_DPB1_2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DPB1_2/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DQA1_1:
    input:
        "output/{sample}_HLA_DQA1_1.fa"
    output:
        "output/flye_DQA1_1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DQA1_1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DQA1_2:
    input:
        "output/{sample}_HLA_DQA1_2.fa"
    output:
        "output/flye_DQA1_2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DQA1_2/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DQB1_1:
    input:
        "output/{sample}_HLA_DQB1_1.fa"
    output:
        "output/flye_DQB1_1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DQB1_1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DQB1_2:
    input:
        "output/{sample}_HLA_DQB1_2.fa"
    output:
        "output/flye_DQB1_2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DQB1_2/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DRB1_1:
    input:
        "output/{sample}_HLA_DRB1_1.fa"
    output:
        "output/flye_DRB1_1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DRB1_1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_DRB1_2:
    input:
        "output/{sample}_HLA_DRB1_2.fa"
    output:
        "output/flye_DRB1_2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir output/flye_DRB1_2/ --threads 4 -m 1000 -i 2|tee {output}"
