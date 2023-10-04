# Evaluating Intragenomic Conflict in the Altruistic, Pheromone-Mediated Retinue Behavior in Honey Bees: allele-specific read counts

Generate read counts in F2 workers at F1 SNPs within Amel_HAv3.1 RefSeq genes

## Contents

-   [Requirements](#requirements)
-   [Process F1 WGS reads](#process-f1-WGS-reads)
    -   Retrieve data from SRA with [`sra-tools`]
    -   Trim adapters with [`fastp`]
    -   Align with [`BWA-MEM`]
    -   Identify SNPs with [`freebayes`]
    -   Generate F1 reference genomes with [`gatk`]
-   [Process F2 RNA-seq reads](#process-f2-rna-seq-libraries)
    -   Retrieve data from SRA with [`sra-tools`]
    -   Trim adapters with [`fastp`]
    -   Align to parent reference genomes with [`tophat2`]
-   [Calculate F2 read coverage at F1 SNPs](#calculate-f2-read-coverage-at-f1-snps)
    -   Generate BED file of SNP positions within genes
    -   Calculate read coverage at each SNP

# Requirements

##### miniconda

<https://docs.conda.io/en/latest/miniconda.html>

##### sra-tools

```
conda create --name sra-tools
conda install -c bioconda -n sra-tools sra-tools
```

##### fastp

```
conda create --name fastp
conda install -c bioconda -n fastp fastp
```

##### bwa

```
conda create --name bwa
conda install -c bioconda -n bwa bwa
```

##### tophat2

```
conda create --name tophat2
conda install -c bioconda -n tophat2 tophat
```

##### freebayes

```
conda create --name freebayers
conda install -c bioconda -n freebayes freebayes
```

##### samtools

```
conda create --name samtools
conda install -c bioconda -n samtools samtools
```

##### bcftools

```
conda create --name bcftools
conda install -c bioconda -n bcftools bcftools
```

##### gatk

```
conda create --name gatk
conda install -c bioconda -n gatk gatk
```

##### bedtools

```
conda create --name bedtools
conda install -c bioconda -n bedtools bedtools
```

# Process F1 WGS reads

### Sample metadata

| SRA         |   ID  |  Block  |
|-------------|-------|---------|
| SRR24049731 | Y12D  | Y12xO20 |
| SRR24049730 | Y12Q  | Y12xO20 |
| SRR24049747 | O20D  | Y12xO20 |
| SRR24049736 | O20Q  | Y12xO20 |
| SRR24049759 | B4D   | B4xW36  |
| SRR24049758 | B4Q   | B4xW36  |
| SRR24049733 | W36D  | B4xW36  |
| SRR24049732 | W36Q  | B4xW36  |
| SRR14654187 | W4D   | LB11xW4 |
| SRR14654186 | W4Q   | LB11xW4 |
| SRR14654189 | LB11D | LB11xW4 |
| SRR14654188 | LB11Q | LB11xW4 |

## Define variables and make directories

```
DIR_WD="set your working directory here"
THREADS="set number of threads for multithreaded operations here"

DIR_RAW=${DIR_WD}/raw
DIR_TRIM=${DIR_WD}/trim
DIR_INDEX=${DIR_WD}/index
DIR_BWA=${DIR_WD}/BWA
DIR_VARIANTS=${DIR_WD}/variants
DIR_ARG=${DIR_WD}/ARGs
mkdir ${DIR_RAW} ${DIR_TRIM} ${DIR_INDEX} ${DIR_BWA} ${DIR_VARIANTS} ${DIR_ARG}

FILES=("SRR24049731" "SRR24049730" "SRR24049747" "SRR24049736" \
       "SRR24049758" "SRR24049759" "SRR24049733" "SRR24049732" \
       "SRR14654186" "SRR14654187" "SRR14654188" "SRR14654189")

### Pick one F2 WGS library to serve as the reference for formatting headers
### If only running a subset of samples, must be one of the libraries that are used
PREF=("SRR24049731")
```

## Retrieve F1 WGS reads

```
conda activate sra-tools

for FILE in "${FILES[@]}"
do
  prefetch -O ${DIR_RAW} ${FILE}
  fasterq-dump -O ${DIR_RAW} ${DIR_RAW}/${FILE}.sra
  rm {FILE}.sra
done

conda deactivate
```

## Trim adapters from F1 WGS reads

```
conda activate fastp

for FILE in "${FILES[@]}"
do
  fastp -w ${THREADS} -i ${FILE}_1.fastq -I ${FILE}_2.fastq \
  -o ${DIR_TRIM}/${FILE}_1.fastq -O ${DIR_TRIM}/${FILE}_2.fastq
done

conda deactivate
```

## Generate BWA alignment index for the Amel_HAv3.1 reference genome

```
cd ${DIR_INDEX}

wget -O Amel_HAv3.1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
gunzip Amel_HAv3.1.fna.gz
mv Amel_HAv3.1.fna Amel_HAv3.1.fasta

wget -O Amel_HAv3.1.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz
gunzip Amel_HAv3.1.gff.gz

conda activate bwa
bwa index Amel_HAv3.1.fasta
conda deactivate
```

## Align F1 WGS reads to the Amel_HAv3.1 reference genome

1. Align with [`bwa mem`]
2. Sort and output bam file with [`samtools`]

```
cd ${DIR_TRIM}

conda activate bwa

for FILE in "${FILES[@]}"
do
  bwa mem -t ${THREADS} ${DIR_INDEX}/Amel_HAv3.1.fasta ${FILE}_1.fastq ${FILE}_2.fastq > ${DIR_BWA}/${FILE}.sam
done

conda deactivate

conda activate samtools

cd ${DIR_BWA}

for FILE in "${FILES[@]}"
do
  samtools view -@ ${THREADS} -O BAM ${FILE}.sam | samtools sort -@ ${THREADS} -O BAM > ${FILE}.bam
done

conda deactivate
```

## Identify SNPs

Using freebayes to account for differences in ploidy between `DIPLOID` and `HAPLOID` samples.

```
DIPLOID=("SRR24049730" "SRR24049736" "SRR24049759" "SRR24049732" "SRR14654186" "SRR14654188")
HAPLOID=("SRR24049731" "SRR24049747" "SRR24049758" "SRR24049733" "SRR14654187" "SRR14654189")

cd ${DIR_INDEX}

conda activate samtools

samtools faidx Amel_HAv3.1.fasta

conda deactivate

cd ${DIR_BWA}

conda activate freebayes

for FILE in "${DIPLOID[@]}"
do
  freebayes -f ${DIR_INDEX}/Amel_HAv3.1.fasta ${FILE}.bam > ${DIR_VARIANTS}/${FILE}.vcf
done

for FILE in "${HAPLOID[@]}"
do
  freebayes -f ${DIR_INDEX}/Amel_HAv3.1.fasta ${FILE}.bam -p 1 > ${DIR_VARIANTS}/${FILE}.vcf
done

conda deactivate
```

## Generate F1 reference genomes

1. Remove heterozygous variants and indels with [`bcftools filter`]
2. Integrate homozygous variants into Amel_HAv3.1 for each sample, separately, using [`gatk FastaAlternateReferenceMaker`]
3. Replace headers in fasta files created by gatk with headers from original reference genome fasta file
4. Create [`tophat2`] alignment index for each F1 reference genome

```
INDEX="${DIR_INDEX}/Amel_HAv3.1.fasta"

cd ${DIR_INDEX}

conda activate gatk

gatk CreateSequenceDictionary -R Amel_HAv3.1.fasta -O Amel_HAv3.1.dict

cd ${DIR_VARIANTS}

for FILE in "${FILES[@]}"
do
   conda activate bcftools
   ### '-e' = "excluxe"; '-i' = "include"; 'GT' = "genotype"; 'TYPE' = "variant type"
   bcftools filter -e 'GT="het"' ${FILE}.vcf \
   | bcftools filter -i 'TYPE="snp"' - > ${FILE}_homozygous_snps.vcf
   bgzip -c ${FILE}_homozygous_snps.vcf > ${FILE}_homozygous_snps.vcf.gz
   tabix -p vcf ${FILE}_homozygous_snps.vcf.gz
   conda deactivate
   conda activate gatk
   gatk FastaAlternateReferenceMaker -R ${INDEX} -O ${DIR_ARG}/${FILE}.fasta -V ${FILE}_homozygous_snps.vcf.gz
   conda deactivate
done

conda deactivate

cd ${DIR_ARG}

### Get lines in PREF beginning with ">" | strip leading ">"
grep ">" ${PREF}.fasta | sed 's/>//g' > ${DIR_INDEX}/bad_headers.txt
### Get lines beginning with ">" | strip text at first space | strip leading ">"
grep ">" ${DIR_INDEX}/Amel_HAv3.1.fasta | sed 's/\s.*$//' | sed 's/>//g' > ${DIR_INDEX}/good_headers.txt
### Combine "bad_headers.txt" and "good_headers.txt" as columns in a tab-separated file
paste -d"\t" ${DIR_INDEX}/bad_headers.txt ${DIR_INDEX}/good_headers.txt > ${DIR_INDEX}/replace_headers.tsv

for FILE in "${FILES[@]}"
do
  ### For each line beginning with ">" in "${FILE}.fasta"
  #### replace the header with the corresponding header from "replace_headers.tsv"
  awk 'FNR==NR{  a[">"$1]=$2;next}$1 in a{  sub(/>/,">"a[$1]"|",$1)}1' \
  ${DIR_INDEX}/replace_headers.tsv ${FILE}.fasta | sed 's/:.*//' > ${FILE}.fa
  rm ${FILE}.fasta
done

conda activate tophat2

for FILE in "${FILES[@]}"
do
  bowtie2-build ${FILE}.fa ${FILE}
done

conda deactivate
```

# Process F2 RNA-seq libraries

### Sample metadata

| SRA         | Queen |   Block   |    Phenotype   |
|-------------|-------|-----------|----------------|
| SRR24049740 | Y12   |  Y12xO20  |   Responsive   |
| SRR24049739 | Y12   |  Y12xO20  |   Responsive   |
| SRR24049738 | Y12   |  Y12xO20  |   Responsive   |
| SRR24049737 | Y12   |  Y12xO20  |  Unresponsive  |
| SRR24049735 | Y12   |  Y12xO20  |  Unresponsive  |
| SRR24049734 | Y12   |  Y12xO20  |  Unresponsive  |
| SRR24049753 | O20   |  Y12xO20  |   Responsive   |
| SRR24049752 | O20   |  Y12xO20  |   Responsive   |
| SRR24049751 | O20   |  Y12xO20  |   Responsive   |
| SRR24049750 | O20   |  Y12xO20  |  Unresponsive  |
| SRR24049749 | O20   |  Y12xO20  |  Unresponsive  |
| SRR24049748 | O20   |  Y12xO20  |  Unresponsive  |
| SRR24049729 | B4    |  B4xW36   |   Responsive   |
| SRR24049728 | B4    |  B4xW36   |   Responsive   |
| SRR24049757 | B4    |  B4xW36   |   Responsive   |
| SRR24049756 | B4    |  B4xW36   |  Unresponsive  |
| SRR24049755 | B4    |  B4xW36   |  Unresponsive  |
| SRR24049754 | B4    |  B4xW36   |  Unresponsive  |
| SRR24049746 | W36   |  B4xW36   |   Responsive   |
| SRR24049745 | W36   |  B4xW36   |   Responsive   |
| SRR24049744 | W36   |  B4xW36   |   Responsive   |
| SRR24049743 | W36   |  B4xW36   |  Unresponsive  |
| SRR24049742 | W36   |  B4xW36   |  Unresponsive  |
| SRR24049741 | W36   |  B4xW36   |  Unresponsive  |
| SRR14654177 |  W4   |  LB11xB4  |   Responsive   |
| SRR14654175 |  W4   |  LB11xB4  |   Responsive   |
| SRR14654183 |  W4   |  LB11xB4  |   Responsive   |
| SRR14654176 |  W4   |  LB11xB4  |  Unresponsive  |
| SRR14654174 |  W4   |  LB11xB4  |  Unresponsive  |
| SRR14654182 |  W4   |  LB11xB4  |  Unresponsive  |
| SRR14654185 | LB11  |  LB11xB4  |   Responsive   |
| SRR14654181 | LB11  |  LB11xB4  |   Responsive   |
| SRR14654179 | LB11  |  LB11xB4  |   Responsive   |
| SRR14654184 | LB11  |  LB11xB4  |  Unresponsive  |
| SRR14654180 | LB11  |  LB11xB4  |  Unresponsive  |
| SRR14654178 | LB11  |  LB11xB4  |  Unresponsive  |



## Define variables and make directories

```
DIR_TOPHAT2=${DIR_WD}/tophat2
DIR_SORT=${DIR_WD}/sort_RNA
DIR_COUNTS=${DIR_WD}/counts
mkdir ${DIR_TOPHAT2} ${DIR_SORT} ${DIR_COUNTS}

FILES=("SRR14654174" "SRR14654175" "SRR14654176" "SRR14654177" "SRR14654178" \
       "SRR14654179" "SRR14654180" "SRR14654181" "SRR14654182" "SRR14654183" \
       "SRR14654184" "SRR14654185" "SRR24049740" "SRR24049739" "SRR24049738" \
       "SRR24049737" "SRR24049735" "SRR24049734" "SRR24049753" "SRR24049752" \
       "SRR24049751" "SRR24049750" "SRR24049749" "SRR24049748" "SRR24049729" \
       "SRR24049728" "SRR24049757" "SRR24049756" "SRR24049755" "SRR24049754" \
       "SRR24049746" "SRR24049745" "SRR24049744" "SRR24049743" "SRR24049742" "SRR24049741")
       
Y12=("SRR24049740" "SRR24049739" "SRR24049738" "SRR24049737" "SRR24049735" "SRR24049734")
O20=("SRR24049753" "SRR24049752" "SRR24049751" "SRR24049750" "SRR24049749" "SRR24049748")
B4=("SRR24049729" "SRR24049728" "SRR24049757" "SRR24049756" "SRR24049755" "SRR24049754")
W36=("SRR24049746" "SRR24049745" "SRR24049744" "SRR24049743" "SRR24049742" "SRR24049741")
LB11=("SRR14654185" "SRR14654184" "SRR14654181" "SRR14654180" "SRR14654179" "SRR14654178")
W4=("SRR14654177" "SRR14654176" "SRR14654175" "SRR14654174" "SRR14654183" "SRR14654182")
```

## Retrieve F2 RNA-seq libraries

```
conda activate sra-tools

cd ${DIR_RAW}

for FILE in "${FILES[@]}"
do
  prefetch ${FILE}
  fasterq-dump -O ${DIR_RAW} ${FILE}
done

conda deactivate
```

## Trim adapters from F2 RNA-seq libraries

```
conda activate fastp

for FILE in "${FILES[@]}"
do
  fastp -w ${THREADS} -i ${FILE}_1.fastq -I ${FILE}_2.fastq \
  -o ${DIR_TRIM}/${FILE}_1.fastq -O ${DIR_TRIM}/${FILE}_2.fastq
done

conda deactivate
```

## Align F2 RNA-seq libraries to respective F1 genomes

```
cd ${DIR_TRIM}

conda activate tophat2

# Y12
for FILE in "${Y12[@]}"
do
  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/Y12Q_${FILE} \
  ${DIR_ARG}/SRR24049730 ${FILE}.fastq

  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/Y12D_${FILE} \
  ${DIR_ARG}/SRR24049731 ${FILE}.fastq
done

# O20
for FILE in "${O20[@]}"
do
  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/O20Q_${FILE} \
  ${DIR_ARG}/SRR24049736 ${FILE}.fastq

  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/O20D_${FILE} \
  ${DIR_ARG}/SRR24049747 ${FILE}.fastq
done

# B4
for FILE in "${B4[@]}"
do
  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/B4Q_${FILE} \
  ${DIR_ARG}/SRR24049758 ${FILE}.fastq

  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/B4D_${FILE} \
  ${DIR_ARG}/SRR24049759 ${FILE}.fastq
done

# W36
for FILE in "${W36[@]}"
do
  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/W36Q_${FILE} \
  ${DIR_ARG}/SRR24049732 ${FILE}.fastq

  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/W36D_${FILE} \
  ${DIR_ARG}/SRR24049733 ${FILE}.fastq
done

# LB11
for FILE in "${LB11[@]}"
do
  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/LB11Q_${FILE} \
  ${DIR_ARG}/SRR14654188 ${FILE}.fastq

  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/LB11D_${FILE} \
  ${DIR_ARG}/SRR14654189 ${FILE}.fastq
done

# W4
for FILE in "${W4[@]}"
do
  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/W4Q_${FILE} \
  ${DIR_ARG}/SRR14654186 ${FILE}.fastq

  tophat2 --library-type fr-firststrand --b2-very-sensitive -p ${THREADS} -o ${DIR_TOPHAT2}/W4D_${FILE} \
  ${DIR_ARG}/SRR14654187 ${FILE}.fastq
done

conda deactivate
```

# Count F2 reads at F1 SNPs

## Filter SNPs

### Remove SNPs that are the same between parents

```
cd ${DIR_VARIANTS}

conda activate bedtools

### bedtools arguments: '-header' = keep header; '-v' = get SNPs in '-a' that are not in 'b'

# Y12
bedtools intersect -header -v -a SRR24049730_homozygous_snps.vcf.gz \
-b SRR24049731_homozygous_snps.vcf.gz > SRR24049730_outer.vcf

bedtools intersect -header -v -a SRR24049731_homozygous_snps.vcf.gz \
-b SRR24049730_homozygous_snps.vcf.gz > SRR24049731_outer.vcf

# O20
bedtools intersect -header -v -a SRR24049736_homozygous_snps.vcf.gz \
-b SRR24049747_homozygous_snps.vcf.gz > SRR24049736_outer.vcf

bedtools intersect -header -v -a SRR24049747_homozygous_snps.vcf.gz \
-b SRR24049736_homozygous_snps.vcf.gz > SRR24049747_outer.vcf

# B4
bedtools intersect -header -v -a SRR24049759_homozygous_snps.vcf.gz \
-b SRR24049758_homozygous_snps.vcf.gz > SRR24049759_outer.vcf

bedtools intersect -header -v -a SRR24049758_homozygous_snps.vcf.gz \
-b SRR24049759_homozygous_snps.vcf.gz > SRR24049758_outer.vcf

# W36
bedtools intersect -header -v -a SRR24049732_homozygous_snps.vcf.gz \
-b SRR24049733_homozygous_snps.vcf.gz > SRR24049732_outer.vcf

bedtools intersect -header -v -a SRR24049733_homozygous_snps.vcf.gz \
-b SRR24049732_homozygous_snps.vcf.gz > SRR24049733_outer.vcf

# LB11
bedtools intersect -header -v -a SRR14654188_homozygous_snps.vcf.gz \
-b SRR14654189_homozygous_snps.vcf.gz > SRR14654188_outer.vcf

bedtools intersect -header -v -a SRR14654189_homozygous_snps.vcf.gz \
-b SRR14654188_homozygous_snps.vcf.gz > SRR14654189_outer.vcf

# W4
bedtools intersect -header -v -a SRR14654186_homozygous_snps.vcf.gz \
-b SRR14654187_homozygous_snps.vcf.gz > SRR14654186_outer.vcf

bedtools intersect -header -v -a SRR14654187_homozygous_snps.vcf.gz \
-b SRR14654186_homozygous_snps.vcf.gz > SRR14654187_outer.vcf

conda deactivate

conda activate bcftools

### grep arguments: "-v '^#'" = get lines that do not begin with '#"
### cat = combine files

# Y12
grep -v '^#' SRR24049730_outer.vcf | cat SRR24049731_outer.vcf - > Y12.vcf
bgzip -c Y12.vcf > Y12.vcf.gz

# O20
grep -v '^#' SRR24049736_outer.vcf | cat SRR24049747_outer.vcf - > O20.vcf
bgzip -c O20.vcf > O20.vcf.gz

# B4
grep -v '^#' SRR24049759_outer.vcf | cat SRR24049758_outer.vcf - > B4.vcf
bgzip -c B4.vcf > B4.vcf.gz

# W36
grep -v '^#' SRR24049732_outer.vcf | cat SRR24049733_outer.vcf - > W36.vcf
bgzip -c W36.vcf > W36.vcf.gz

# LB11
grep -v '^#' SRR14654188_outer.vcf | cat SRR14654189_outer.vcf - > LB11.vcf
bgzip -c LB11.vcf > LB11.vcf.gz

# W4
grep -v '^#' SRR14654186_outer.vcf | cat SRR14654187_outer.vcf - > W4.vcf
bgzip -c W4.vcf > W4.vcf.gz

conda deactivate
```

### Keep variants that are concordant between crosses (consensus variants)

```
conda activate bedtools

### bedtools arguments: '-u' = get SNPs in '-a' that are in '-b'

# Y12xO20
bedtools intersect -header -u -a Y12.vcf.gz -b O20.vcf.gz > Y12xO20_aSet.vcf

# B4xW36
bedtools intersect -header -u -a B4.vcf.gz -b W36.vcf.gz > B4xW36_aSet.vcf

# LB11xW4
bedtools intersect -header -u -a LB11.vcf.gz -b W4.vcf.gz > LB11xW4_aSet.vcf

### get lines that do not begin with '#" | 
#### keep only columns 1 and 2, repeating column 2 as column 3 |
#### populate column 4 with sequential numbers (1-n) |
#### in column 4, paste "snp_" to beginning of each line
 
# Y12xO20
grep -v '^#' Y12xO20_aSet.vcf | awk -v OFS="\t" '{print $1, $2, $2}' \
 | awk -v OFS="\t" '$4=(FNR FS $4)' \
 | awk -v OFS="\t" '{print $1, $2, $3, "snp_"$4}' > Y12xO20_aSet.bed
 
# B4xW36
grep -v '^#' B4xW36_aSet.vcf | awk -v OFS="\t" '{print $1, $2, $2}' \
 | awk -v OFS="\t" '$4=(FNR FS $4)' \
 | awk -v OFS="\t" '{print $1, $2, $3, "snp_"$4}' > B4xW36_aSet.bed

# LB11xW4
grep -v '^#' LB11xW4_aSet.vcf | awk -v OFS="\t" '{print $1, $2, $2}' \
 | awk -v OFS="\t" '$4=(FNR FS $4)' \
 | awk -v OFS="\t" '{print $1, $2, $3, "snp_"$4}' > LB11xW4_aSet.bed

conda deactivate
```

## Generate BED file of SNP positions within genes

BED files generated from this procedure for each block, containing SNPs within Amel_HAv3.1 RefSeq genes, can be found in the [SNP_bed_files.zip](https://github.com/sbresnahan/IGC-retinue/blob/main/SNP_bed_files.zip) archive in this repository.

```
conda activate bedtools

cd ${DIR_INDEX}

### Get lines with "gene" in column 3 | keep columns 1-8
awk '$3 == "gene" { print $0 }' Amel_HAv3.1.gff | awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8}' > Amel_HAv3.1_genes.txt
### Get lines with "gene" in column 3 | keep column 9 | get substring beginning with 'GeneID' |
#### cut text before and including first appearance of ":" | cut text after and including first appearance of ";" |
#### cut text after and including first appearance of "," | paste "LOC" to beginning of each line
awk '$3 == "gene" { print $0 }' Amel_HAv3.1.gff | awk '{print $9}' | grep -Po 'GeneID[^\s]*' \
| cut -d':' -f2 | cut -d';' -f1 | cut -d',' -f1 | sed -e 's/^/LOC/' > Amel_HAv3.1_geneIDs.txt
#### Combine "Amel_HAv3.1_genes.txt" and "Amel_HAv3.1_geneIDs.txt" as columns in a tab-separated file
paste -d'\t' Amel_HAv3.1_genes.txt Amel_HAv3.1_geneIDs.txt > Amel_HAv3.1_genes.gff3

### bedtools intersect arguments: '-wb' = get genes in '-a' that intersect with SNPs in -'b', output coordinates in '-b' |
#### keep columns 10-12, 13 joined with column 9 by ":", & columns 6-7 |
#### remove rows containing "NC_001566.1" (corresponding to mitochondrial sequences) |
#### sort by chromosome and position

# Y12xO20
bedtools intersect -wb -a Amel_HAv3.1_genes.gff3 -b ${DIR_VARIANTS}/Y12xO20_aSet.bed \
| awk -v OFS="\t" '{print $10, $11, $12, $13 ":" $9, $6, $7}' \
| grep -v '^NC_001566.1' > ${DIR_VARIANTS}/Y12xO20_SNPs_for_analysis.bed

# B4xW36
bedtools intersect -wb -a Amel_HAv3.1_genes.gff3 -b ${DIR_VARIANTS}/B4xW36_aSet.bed \
| awk -v OFS="\t" '{print $10, $11, $12, $13 ":" $9, $6, $7}' \
| grep -v '^NC_001566.1' > ${DIR_VARIANTS}/B4xW36_SNPs_for_analysis.bed

# LB11xW4
bedtools intersect -wb -a Amel_HAv3.1_genes.gff3 -b ${DIR_VARIANTS}/LB11xW4_aSet.bed \
| awk -v OFS="\t" '{print $10, $11, $12, $13 ":" $9, $6, $7}' \
| grep -v '^NC_001566.1' > ${DIR_VARIANTS}/LB11xW4_SNPs_for_analysis.bed

conda deactivate

cd ${DIR_VARIANTS}

# Y12xO20
sort --parallel=${THREADS} -k1,1 -k2,2n Y12xO20_SNPs_for_analysis.bed > Y12xO20_SNPs_for_analysis_sorted.bed

# B4xW36
sort --parallel=${THREADS} -k1,1 -k2,2n B4xW36_SNPs_for_analysis.bed > B4xW36_SNPs_for_analysis_sorted.bed

# LB11xW4
sort --parallel=${THREADS} -k1,1 -k2,2n LB11xW4_SNPs_for_analysis.bed > LB11xW4_SNPs_for_analysis_sorted.bed
```

## Calculate F2 read coverage at F1 SNPs

1. Filter tophat2 alignment files to remove primary alignments with mismatches (`-tag XM:0`) and secondary alignments (`-isPrimaryAlignment true`) with [`bamtools filter`]
2. Convert bam to coordinate-sorted bed file with [`bedtools bamtobed`]
3. Intersect F2 read alignment bed file with corresponding F1 SNP:gene files, accounting for strandedness (`-S`), and count read coverage at SNP:genes with [`bedtools intersect`]

```
# Y12xO20

## Y12
for FILE in "${Y12[@]}"
do
  conda activate bamtools
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/Y12Q_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/Y12Q_${FILE}/nomm_hits.bam
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/Y12D_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/Y12D_${FILE}/nomm_hits.bam

  conda deactivate
  conda activate bedtools

  bedtools bamtobed -i ${DIR_ALIGN}/Y12Q_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/Y12Q_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/Y12Q_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/Y12Q_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/Y12xO20_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/Y12Q_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/Y12Q_${FILE}.txt

  bedtools bamtobed -i ${DIR_ALIGN}/Y12D_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/Y12D_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/Y12D_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/Y12D_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/Y12xO20_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/Y12D_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/Y12D_${FILE}.txt

  conda deactivate
done

## O20
for FILE in "${O20[@]}"
do
  conda activate bamtools
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/O20Q_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/O20Q_${FILE}/nomm_hits.bam
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/O20D_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/O20D_${FILE}/nomm_hits.bam

  conda deactivate
  conda activate bedtools

  bedtools bamtobed -i ${DIR_ALIGN}/O20Q_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/O20Q_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/O20Q_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/O20Q_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/Y12xO20_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/O20Q_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/O20Q_${FILE}.txt

  bedtools bamtobed -i ${DIR_ALIGN}/O20D_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/O20D_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/O20D_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/O20D_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/Y12xO20_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/O20D_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/O20D_${FILE}.txt

  conda deactivate
done


# B4xW36

## B4
for FILE in "${B4[@]}"
do
  conda activate bamtools
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/B4Q_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/B4Q_${FILE}/nomm_hits.bam
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/B4D_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/B4D_${FILE}/nomm_hits.bam

  conda deactivate
  conda activate bedtools

  bedtools bamtobed -i ${DIR_ALIGN}/B4Q_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/B4Q_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/B4Q_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/B4Q_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/B4xW36_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/B4Q_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/B4Q_${FILE}.txt

  bedtools bamtobed -i ${DIR_ALIGN}/B4D_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/B4D_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/B4D_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/B4D_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/B4xW36_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/B4D_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/B4D_${FILE}.txt

  conda deactivate
done

## W36
for FILE in "${W36[@]}"
do
  conda activate bamtools
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/W36Q_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/W36Q_${FILE}/nomm_hits.bam
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/W36D_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/W36D_${FILE}/nomm_hits.bam

  conda deactivate
  conda activate bedtools

  bedtools bamtobed -i ${DIR_ALIGN}/W36Q_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/W36Q_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/W36Q_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/W36Q_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/B4xW36_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/W36Q_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/W36Q_${FILE}.txt

  bedtools bamtobed -i ${DIR_ALIGN}/W36D_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/W36D_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/W36D_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/W36D_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/B4xW36_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/W36D_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/W36D_${FILE}.txt

  conda deactivate
done

# LB11xW4

## LB11
for FILE in "${LB11[@]}"
do
  conda activate bamtools
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/LB11Q_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/LB11Q_${FILE}/nomm_hits.bam
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/LB11D_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/LB11D_${FILE}/nomm_hits.bam

  conda deactivate
  conda activate bedtools

  bedtools bamtobed -i ${DIR_ALIGN}/LB11Q_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/LB11Q_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/LB11Q_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/LB11Q_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/LB11xW4_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/LB11Q_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/LB11Q_${FILE}.txt

  bedtools bamtobed -i ${DIR_ALIGN}/LB11D_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/LB11D_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/LB11D_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/LB11D_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/LB11xW4_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/LB11D_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/LB11D_${FILE}.txt

  conda deactivate
done

## W4
for FILE in "${W4[@]}"
do
  conda activate bamtools
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/W4Q_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/W4Q_${FILE}/nomm_hits.bam
  bamtools filter -tag XM:0 -isPrimaryAlignment true -in ${DIR_ALIGN}/W4D_${FILE}/accepted_hits.bam \
  -out ${DIR_ALIGN}/W4D_${FILE}/nomm_hits.bam

  conda deactivate
  conda activate bedtools

  bedtools bamtobed -i ${DIR_ALIGN}/W4Q_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/W4Q_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/W4Q_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/W4Q_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/LB11xW4_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/W4Q_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/W4Q_${FILE}.txt

  bedtools bamtobed -i ${DIR_ALIGN}/W4D_${FILE}/nomm_hits.bam \
  > ${DIR_ALIGN}/W4D_${FILE}/nomm_hits.bed

  sort --parallel=${THREADS} -k1,1 -k2,2n ${DIR_ALIGN}/W4D_${FILE}/nomm_hits.bed \
  > ${DIR_SORT}/W4D_${FILE}_nomm_hits_sorted.bed

  bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/LB11xW4_SNPs_for_analysis_sorted.bed \
  -b ${DIR_SORT}/W4D_${FILE}_nomm_hits_sorted.bed \
  > ${DIR_COUNTS}/W4D_${FILE}.txt

  conda deactivate
done
```
