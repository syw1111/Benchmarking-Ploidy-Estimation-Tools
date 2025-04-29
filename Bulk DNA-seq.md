## TITAN
```
source activate titan
snakemake -s TitanCNA.snakefile --cores 24
```
## Sequenza
```
source activate Sequenza
sequenza-utils bam2seqz -gc home/WuLab/Benchmark/hg38.gc50Base.txt.gz -n home/Benchmark/Alignment/Data1/Preprocess_Data1/Normal_1.hg38.rmdup.filter.sorted.bam -t home/Benchmark/Alignment/Data1/Preprocess_Data1/Tumor_1.hg38.rmdup.filter.sorted.bam -F home/WuLab/Benchmark/human/genome.fa -o Tumor_rep1.seqz.gz -f illumina
sequenza-utils seqz_binning -w 50 -s Tumor_rep1.seqz.gz -o Tumor_rep1.small.seqz.gz
zcat Tumor_rep1.small.seqz.gz|grep -v '_random'|grep -v 'chrM' | bgzip > Tumor_rep1.nochrM.small.seqz.gz
tabix -f -s 1 -b 2 -e 2 -S 1 Tumor_rep1.nochrM.small.seqz.gz
mkdir Tumor_rep1
Rscript Tumor-1.R
```
**Tumor-1.R**
```
library(sequenza)
test <- sequenza.extract("Tumor_rep1.nochrM.small.seqz.gz", verbose = FALSE)
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test,
    cp.table = CP, sample.id = "Tumor_rep1",
            out.dir="Tumor_rep1")
```
## Sclust
```
Sclust bamprocess -t Tumor_rep1/Tumor_rep1_T.bam -n Tumor_rep1/Tumor_rep1_N.bam -o Tumor_rep1/Tumor_rep1 -build hg38 -part 2 -r chr1
...
Sclust bamprocess -t Tumor_rep1/Tumor_rep1_T.bam -n Tumor_rep1/Tumor_rep1_N.bam -o Tumor_rep1/Tumor_rep1 -build hg38 -part 2 -r chrX
Sclust bamprocess -i Tumor_rep1/Tumor_rep1 -o Tumor_rep1/Tumor_rep1
input_vcf="home/Benchmark/Alignment/Data1/Preprocess_Data1/GATK_SNV/Tumor_1_mutation/Tumor_1.hg38.filtered.vcf.gz"
output_vcf="home/WuLab/Benchmark/Sclust/Purity/Tumor_rep1/Tumor_rep1_mutations.vcf"
zcat $input_vcf | grep -v '^##' > tmp_data.txt
zcat $input_vcf | grep '^##' > $output_vcf
echo "##INFO=<ID=DP,Number=1,Type=Integer,Description='Read Depth Tumor'>" >> $output_vcf
echo "##INFO=<ID=DP_N,Number=1,Type=Integer,Description='Read Depth Normal'>" >> $output_vcf
echo "##INFO=<ID=AF,Number=A,Type=Float,Description='Allelic Frequency Tumor'>" >> $output_vcf
echo "##INFO=<ID=AF_N,Number=A,Type=Float,Description='Allelic Frequency Normal'>" >> $output_vcf
echo "##INFO=<ID=FR,Number=1,Type=Float,Description='Forward-Reverse Score'>" >> $output_vcf
echo "##INFO=<ID=TG,Number=1,Type=String,Description='Target Name (Genome Partition)'>" >> $output_vcf
echo "##INFO=<ID=DB,Number=0,Type=Flag,Description='dbSNP Membership'>" >> $output_vcf
grep '^#CHROM' tmp_data.txt >> $output_vcf
awk -v OFS="\t" '
BEGIN {
        header = "DP=.;DP_N=.;AF=.;AF_N=.;FR=.;TG=.;DB=."
}
{
        if ($1 ~ /^#/) next;
        chrom=$1; pos=$2; id=$3; ref=$4; alt=$5; qual=$6; filter=$7; info=$8;
        split($9, format_tags, ":");
        split($10, tumor_data, ":");
        split($11, normal_data, ":");
        for (i in format_tags) {
                if (format_tags[i] == "DP") {
                        dp_tumor = tumor_data[i];
                        dp_normal = normal_data[i];
                }
                if (format_tags[i] == "AF") {
                        af_tumor = tumor_data[i];
                        af_normal = normal_data[i];
                }
        }
        new_info = "DP=" dp_tumor ";DP_N=" dp_normal ";AF=" af_tumor ";AF_N=" af_normal ";FR=.;TG=.;DB=."
        print chrom, pos, id, ref, alt, qual, filter, new_info
}' tmp_data.txt >> $output_vcf
rm tmp_data.txt
Sclust cn -rc Tumor_rep1/Tumor_rep1_rcount.txt -snp Tumor_rep1/Tumor_rep1_snps.txt -vcf Tumor_rep1/Tumor_rep1_mutations.vcf -o Tumor_rep1/Tumor_rep1
```
## PURPLE
```
mkdir Tumor_rep1
java -jar -Xmx8G home/software/WGD_software/cobalt-1.14.1.jar -reference Normal_rep1 -reference_bam home/Benchmark/Alignment/Data1/Preprocess_Data1/Normal_1.hg38.rmdup.filter.sorted.bam -tumor Tumor_rep1 -tumor_bam home/Benchmark/Alignment/Data1/Preprocess_Data1/Tumor_1.hg38.rmdup.filter.sorted.bam -gc_profile home/ref/Human/hg38/HMFTools/v6_0/ref/38/copy_number/GC_profile.1000bp.38.cnp -output_dir Tumor_rep1 -threads 8
java -jar -Xmx8G home/software/WGD_software/amber-3.9.jar -reference Normal_rep1 -reference_bam home/Benchmark/Alignment/Data1/Preprocess_Data1/Normal_1.hg38.rmdup.filter.sorted.bam -tumor Tumor_rep1 -tumor_bam home/Benchmark/Alignment/Data1/Preprocess_Data1/Tumor_1.hg38.rmdup.filter.sorted.bam -output_dir Tumor_rep1 -loci home/ref/Human/hg38/HMFTools/v6_0/ref/38/GermlineHetPon.38.vcf.gz -threads 6 -ref_genome_version V38
java -jar -Xmx8G home/software/WGD_software/hmftools-purple-v3.8.4/purple_v3.8.4.jar -reference Normal_rep1 -tumor Tumor_rep1 -amber Tumor_rep1 -cobalt Tumor_rep1 -gc_profile home/ref/Human/hg38/HMFTools/v6_0/ref/38/copy_number/GC_profile.1000bp.38.cnp -ref_genome_version 38 -ref_genome home/WuLab/Benchmark/human/genome.fa -ensembl_data_dir home/ref/Human/hg38/HMFTools/v6_0/ref/38/common/ensembl_data -output_dir Tumor_rep1
```
## Accucopy
```
HOST_DIR=home/WuLab/Benchmark/Accucopy/
CONTAINER_DIR=/mnt
docker run --rm -v ${HOST_DIR}:${CONTAINER_DIR} -w /usr/local/Accucopy imagehub.cstcloud.cn/polyactis/accucopy python main.py -c ${CONTAINER_DIR}/configure -t ${CONTAINER_DIR}/Tumor_1.hg38.rmdup.filter.sorted.bam -n ${CONTAINER_DIR}/Normal_1.hg38.rmdup.filter.sorted.bam -o ${CONTAINER_DIR}/Tumor_rep1/output --snp_output_dir ${CONTAINER_DIR}/Tumor_rep1/snps
```


