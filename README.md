# HLA genotyping of KOLF2.1

```bash
# on helix:
# retrieve mapped reads from google cloud bucket
module load google-cloud-sdk
gcloud init     # log in using account with download permissions to gs://singlecellindi
gsutil -m cp gs://singlecellindi/WGS/Jax/IlluminaWGS/hg38/bams/KOLF2-ARID2-A2.bam .
gsutil -m cp gs://singlecellindi/WGS/Jax/IlluminaWGS/hg38/bams/KOLF2-ARID2-A2.bai .


# on biowulf:
# request interactive session
sinteractive --mem=40g --cpus-per-task=4 --gres=lscratch:100


# on compute node
module load xHLA/2018-04-04
module load bedtools/2.30.0
module load samtools/1.16.1
module load edirect/17.1
module load bwa/0.7.17


# set directory variables
dir=$(pwd)
bamfilename='KOLF2-ARID2-A2.bam'
bamfile="${dir}/${bamfilename}"
inbam='KOLF.bam'
threads=4


# work in lscratch
TMPDIR="/lscratch/${SLURM_JOB_ID}"
cd "${TMPDIR}"


# retrieve sambama/0.8.2 binary
wget -O sambamba.gz \
    'https://github.com/biod/sambamba/releases/download/v0.8.2/sambamba-0.8.2-linux-amd64-static.gz'
gunzip sambamba.gz
chmod +x sambamba


# retrieve xHLA git repo
git clone https://github.com/humanlongevity/HLA.git && cd HLA


# link original files to current directory
ln -s "${bamfile}" "${inbam}"
ln -s "${bamfile%.bam}.bai" "${inbam}.bai"


# Pull out unmapped, chr6:29844528-33100696, HLA contig, or alt contig reads
../sambamba view \
    -f "bam" -h -p -l 0 -t ${threads} \
    -F "unmapped or mate_is_unmapped or (ref_name == 'chr6' and (position > 29844528 and position < 33100696)) or ref_name =~ /^HLA|chr6    *alt/" \
    -o HLA.bam \
    ${inbam}


# convert from bam to paired fastq
../sambamba sort -p -n -t ${threads} -o - HLA.bam | \
    bamToFastq -i /dev/stdin -fq "extracted.1.fq" -fq2 "extracted.2.fq"
# some number of reads could not find mate in BAM and were skipped,
# leaving ~1.5M reads mapped to chr6


# retrieve chr6 assembly for mapping (not included in xHLA repo)
esearch -db nucleotide -query 'NC_000006.12' | efetch -format fasta > data/chr6/hg38.chr6.fna


# modify header to match xHLA expectation
sed -i 's/^>.*$/>chr6/' data/chr6/hg38.chr6.fna && \
bwa index data/chr6/hg38.chr6.fna


# re-map reads to chr6 reference for xHLA genotyping
bwa mem -t ${threads} \
'data/chr6/hg38.chr6.fna' \
'extracted.1.fq' \
'extracted.2.fq' | \
samtools view -b - | ../sambamba sort -t ${threads} -o KOLF2-hg38.chr6.bam /dev/stdin && \
../sambamba index KOLF2-hg38.chr6.bam

# sanity check of mapped reads
#   $ samtools idxstats full.bam
#   chr6    170805979       1551266 217344
#   *       0       0       1361782


# extract chr6:29844528-33100696
samtools view -o KOLF2-HLA.bam -b KOLF2-hg38.chr6.bam chr6:29844528-33100696
../sambamba index KOLF2-HLA.bam

# sanity check of mapped reads
#   $ samtools idxstats KOLF-hla.bam
#   chr6    170805979       928971  6382
#   *       0       0       0


xhla --sample_id KOLF2 --input_bam_path KOLF2-HLA.bam --output_path KOLF_HLA

samtools depth KOLF2-HLA.bam > depth.txt
```

## xHLA solution:
```
A*01:01
A*11:01
B*08:01
B*51:01
C*07:01
C*15:02
DPB1*09:01
DPB1*14:01
DQB1*03:02
DQB1*05:01
DRB1*01:01
DRB1*04:04
```

## Visualizing depth
```R
# before running R, on biowulf shell:
# module load R/3.6.3

library(data.table)
library(ggplot2)

dat <- fread('depth.txt')
setnames(dat, c('contig','position','depth'))

dat.ag <- dat[, .N, by=depth]
g <- ggplot(dat.ag[depth<100], aes(x=depth, y=N)) + geom_bar(stat='identity', position='dodge')

ggsave(g, file='depth.png')
```

![read-depth](/depth.png)
