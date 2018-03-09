# Fetch raw data

## Download VCF calls for autosomes from 1000 Genomes

```bash
mkdir -p 1000_genomes
cd 1000_genomes
PREFIX='ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502'
export SUFFIX='phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'
rm -f 1000_genome_links.txt
for CHROM in {1..22}; do
    for EXT in "vcf.gz" "vcf.gz.tbi"; do
        echo $PREFIX/ALL.chr$CHROM.$SUFFIX.$EXT >> 1000_genome_links.txt
    done
done
aria2c --input-file 1000_genome_links.txt \
    --split 5 --max-concurrent-downloads 10 --max-connection-per-server 16
touch *.tbi
cd ..
```

# Prepare filtering

## Install dependencies (bedops)

```bash
wget -q https://github.com/bedops/bedops/releases/download/v2.4.30/bedops_linux_x86_64-v2.4.30.tar.bz2
tar --extract --file=bedops_linux_x86_64-v2.4.30.tar.bz2
mv bin/* /usr/bin/
```

## Prepare exome mask

```bash
wget -q 'ftp://ftp.ensembl.org/pub/grch37/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'
gzip -d Homo_sapiens.GRCh37.87.gtf.gz
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' \
    Homo_sapiens.GRCh37.87.gtf | gtf2bed - > annotation.bed
grep 'gene_biotype "protein_coding"' annotation.bed \
    | awk '{print "chr" $0}' | sed -e 's|chrMT|chrM|g' | bedtools sort > coding.bed
awk '{if ($8 == "exon") print $0}' coding.bed > exons.bed
sed -e 's|^chr||' exons.bed | cut -f-3 > exons.filter.bed
bedtools sort -i exons.filter.bed | bedtools merge -i - > exons.merged.bed
```

## Split mask by chromosome
```bash
mkdir -p exome_mask
awk -F$'\t' '{print $0 > "exome_mask/" $1 ".bed"}' exons.merged.bed
```

# Filter VCF by exome

```bash
mkdir -p 1000_exome
filter_vcf () {
    INFILE=$1
    MASK=$2
    OUTFILE=$3
    zcat $INFILE | bedtools intersect -sorted -a - -b $MASK -header | bgzip -c > $OUTFILE
    tabix $OUTFILE
}
export -f filter_vcf
parallel filter_vcf \
    1000_genomes/ALL.chr{}.$SUFFIX.vcf.gz exome_mask/{}.bed 1000_exome/{}.vcf.gz ::: {1..22}
```

# Downsample genotypes

* `0|0` → 0
* `0|1`, `1|0` → 1
* `1|1` → 2

```bash
mkdir -p 1000_downsampled
downsample_vcf () {
    INFILE=$1
    OUTFILE=$2
    bcftools query -f '[%GT ]\n' $INFILE \
        | sed -e 's.0|0.0.g' -e 's.0|1.1.g' -e 's.1|0.1.g' -e 's.1|1.2.g' \
        | tr -d ' ' | gzip -c > $OUTFILE
}
export -f downsample_vcf
parallel downsample_vcf \
     1000_exome/{}.vcf.gz 1000_downsampled/{}.txt.gz ::: {1..22}
```
