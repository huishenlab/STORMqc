#!/bin/bash
#SBATCH -J gather_rna_resources
#SBATCH -o gather_rna_resources.log
#SBATCH --ntasks 16
#SBATCH --time 24:00:00
#SBATCH --mem=128G

cd ${SLURM_SUBMIT_DIR}

# Load needed tools
module load bbc2/fastx_toolkit/fastx_toolkit-0.0.13
module load bbc2/STAR/STAR-2.7.11a

# ERCC spike-ins
#-----------------------------------------------------------------------------------------------------------------------
# Download FASTA and GTF for ERCC spike-ins
wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip

# ERCC GTF needs to be slightly modified to be able to append to GENCODE GTF
python patch_ercc_gtf.py

# Clean up
rm ERCC92.zip
#-----------------------------------------------------------------------------------------------------------------------

# Human reference
#-----------------------------------------------------------------------------------------------------------------------
# Download GTF
wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz | \
gunzip > gencode.v45.primary_assembly.annotation.gtf

# Download FASTA
wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz | \
gunzip > GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
#-----------------------------------------------------------------------------------------------------------------------

# Human enhancers and repeats
#-----------------------------------------------------------------------------------------------------------------------
# Enhancers
wget -O- https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz | \
zcat | \
awk 'BEGIN{ OFS="\t" }{ print $1, $2, $3 }' | \
sort -k1,1 -k2,2n -k3,3n | \
gzip > F5.hg38.enhancers.sorted.bed.gz

# Repeats
wget -O- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz | \
zcat | \
awk 'BEGIN{ OFS="\t" }{ print $6, $7, $8 }' | \
sort -k1,1 -k2,2n -k3,3n | \
gzip > rmsk.hg38.sorted.bed.gz
#-----------------------------------------------------------------------------------------------------------------------

# Mouse reference
#-----------------------------------------------------------------------------------------------------------------------
# Download GTF
wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.primary_assembly.annotation.gtf.gz | \
gunzip > gencode.vM34.primary_assembly.annotation.gtf

# Download FASTA
wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/GRCm39.primary_assembly.genome.fa.gz | \
gunzip > GRCm39.primary_assembly.genome.fa
samtools faidx GRCm39.primary_assembly.genome.fa
#-----------------------------------------------------------------------------------------------------------------------

# Mouse enhancers and repeats
#-----------------------------------------------------------------------------------------------------------------------
# Enhancers
wget -O- https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.bed.gz | \
zcat | \
awk 'BEGIN{ OFS="\t" }{ print $1, $2, $3 }' | \
sort -k1,1 -k2,2n -k3,3n | \
gzip > F5.mm10.enhancers.sorted.bed.gz

# Repeats
wget -O- https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz | \
zcat | \
awk 'BEGIN{ OFS="\t" }{ print $6, $7, $8 }' | \
sort -k1,1 -k2,2n -k3,3n | \
gzip > rmsk.mm10.sorted.bed.gz
#-----------------------------------------------------------------------------------------------------------------------

# Merge files
#-----------------------------------------------------------------------------------------------------------------------
# Human
cat GRCh38.primary_assembly.genome.fa ERCC92.fa > GRCh38_ERCC.merged.fa
cat gencode.v45.primary_assembly.annotation.gtf ERCC92.genes.patched.gtf > GRCh38_ERCC.merged.gtf

awk 'BEGIN{ OFS="\t" }{ if ($3 == "gene") { print $1, $4, $5 } }' GRCh38_ERCC.merged.gtf | \
sort -k1,1 -k2,2n -k3,3n | \
gzip > GRCh38_ERCC.sorted.bed.gz

zcat GRCh38_ERCC.sorted.bed.gz F5.hg38.enhancers.sorted.bed.gz rmsk.hg38.sorted.bed.gz | \
sort -k1,1 -k2,2n -k3,3n | \
bedtools merge -i - | \
gzip > GRCh38_ERCC_enhancers_rmsk.merged.sorted.bed.gz

# Mouse
cat GRCm39.primary_assembly.genome.fa ERCC92.fa > GRCm39_ERCC.merged.fa
cat gencode.vM34.primary_assembly.annotation.gtf ERCC92.genes.patched.gtf > GRCm39_ERCC.merged.gtf

awk 'BEGIN{ OFS="\t" }{ if ($3 == "gene") { print $1, $4, $5 } }' GRCm39_ERCC.merged.gtf | \
sort -k1,1 -k2,2n -k3,3n | \
gzip > GRCm39_ERCC.sorted.bed.gz

zcat GRCm39_ERCC.sorted.bed.gz F5.mm10.enhancers.sorted.bed.gz rmsk.mm10.sorted.bed.gz | \
sort -k1,1 -k2,2n -k3,3n | \
bedtools merge -i - | \
gzip > GRCm39_ERCC_enhancers_rmsk.merged.sorted.bed.gz
#-----------------------------------------------------------------------------------------------------------------------

# Create index
#-----------------------------------------------------------------------------------------------------------------------
# Human
mkdir star_index_human
STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ${SLURM_SUBMIT_DIR}/star_index_human \
    --genomeFastaFiles GRCh38_ERCC.merged.fa \
    --sjdbGTFfile GRCh38_ERCC.merged.gtf

# Mouse
mkdir star_index_mouse
STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ${SLURM_SUBMIT_DIR}/star_index_mouse \
    --genomeFastaFiles GRCm39_ERCC.merged.fa \
    --sjdbGTFfile GRCm39_ERCC.merged.gtf
#-----------------------------------------------------------------------------------------------------------------------

# Create list of gene IDs
#-----------------------------------------------------------------------------------------------------------------------
# Human
cat GRCh38_ERCC.merged.gtf | \
awk -v rrna="gene_ids.hg38.rrna.txt" -v mito="gene_ids.hg38.mito.txt" -v ercc="gene_ids.hg38.ercc.txt" 'BEGIN{ FS="\t" }{
    if ($3 == "gene") {
        split($9,a,"; ")
        for (i in a) {
            split(a[i],b," ")
            gsub(/"/, "", b[2])
            d[b[1]] = b[2]
        }
        if (match(d["gene_type"], "rRNA")) {
            print d["gene_id"], d["gene_name"] > rrna
        }
        if (match(d["gene_type"], "Mt_tRNA")) {
            print d["gene_id"], d["gene_name"] > mito
        }
        if (match(d["gene_type"], "ercc_control")) {
            print d["gene_id"], d["gene_name"] > ercc
        }
    }
}'

# Mouse
cat GRCm39_ERCC.merged.gtf | \
awk -v rrna="gene_ids.mm10.rrna.txt" -v mito="gene_ids.mm10.mito.txt" -v ercc="gene_ids.mm10.ercc.txt" 'BEGIN{ FS="\t" }{
    if ($3 == "gene") {
        split($9,a,"; ")
        for (i in a) {
            split(a[i],b," ")
            gsub(/"/, "", b[2])
            d[b[1]] = b[2]
        }
        if (match(d["gene_type"], "rRNA")) {
            print d["gene_id"], d["gene_name"] > rrna
        }
        if (match(d["gene_type"], "Mt_tRNA")) {
            print d["gene_id"], d["gene_name"] > mito
        }
        if (match(d["gene_type"], "ercc_control")) {
            print d["gene_id"], d["gene_name"] > ercc
        }
    }
}'
#-----------------------------------------------------------------------------------------------------------------------
