## ASCAT
```
source activate r4.3.1
# For NGS in Tumor-normal paired mode
Rscript Run_ASCAT.r ${tumor} ${normal}
# For NGS in Tumor-only mode
Rscript Run_ASCAT_TumorOnly.r ${tumor}
# For TGS in Tumor-normal paired mode
Rscript Run_ASCAT_LongRead.r ${tumor} ${normal}
```
**Run_ASCAT.r**
```
args = (commandArgs(TRUE))
library(ASCAT)
Tumour.bam <- paste(args[[1]], sep="")
Normal.bam <- paste(args[[2]], sep="")
ascat.prepareHTS(
  tumourseqfile = Tumour.bam,
  normalseqfile = Normal.bam,
  tumourname = "Tumour_mix",
  normalname = "HCC1395BL",
  allelecounter_exe = "alleleCounter",
  alleles.prefix = "G1000_alleles_hg38_chr",
  loci.prefix = "G1000_loci_hg38_chr",
  gender = "XX",
  genomeVersion = "hg38",
  nthreads = 8,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt")
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = "XX", genomeVersion = "hg38")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_G1000_hg38.txt", replictimingfile = "RT_G1000_hg38.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = TRUE)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
output_file <- "purity_ploidy_results.txt"
file_conn <- file(output_file, "w")
writeLines(paste("Purity: ", ascat.output$purity), file_conn)
writeLines(paste("Ploidy: ", ascat.output$ploidy), file_conn)
close(file_conn)
```
**Run_ASCAT_TumorOnly.r**
```
args = (commandArgs(TRUE))
library(ASCAT)
Tumour.bam <- paste(args[[1]], sep="")
ascat.prepareHTS(
  tumourseqfile = Tumour.bam,
  tumourname = "Tumour_mix",
  allelecounter_exe = "alleleCounter",
  alleles.prefix = "G1000_alleles_hg38_chr",
  loci.prefix = "G1000_loci_hg38_chr",
  gender = "XX",
  genomeVersion = "hg38",
  nthreads = 16,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt")
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", gender = "XX", genomeVersion = "hg38")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_G1000_hg38.txt", replictimingfile = "RT_G1000_hg38.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
gg = ascat.predictGermlineGenotypes(ascat.bc, platform = "WGS_hg38_50X")
ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=gg)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, write_segments = TRUE)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
output_file <- "purity_ploidy_results.txt"
file_conn <- file(output_file, "w")
writeLines(paste("Purity: ", ascat.output$purity), file_conn)
writeLines(paste("Ploidy: ", ascat.output$ploidy), file_conn)
close(file_conn)
```
**Run_ASCAT_LongRead.r**
```
args = (commandArgs(TRUE))
library(ASCAT)
Tumour.bam <- paste(args[[1]], sep="")
Normal.bam <- paste(args[[2]], sep="")
ascat.prepareHTS(
  tumourseqfile = Tumour.bam,
  normalseqfile = Normal.bam,
  tumourname = "Tumour_mix",
  normalname = "HCC1395BL",
  allelecounter_exe = "alleleCounter",
  alleles.prefix = "G1000_alleles_hg38_chr",
  loci.prefix = "G1000_loci_hg38_chr",
  gender = "XX",
  genomeVersion = "hg38",
  nthreads = 8,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt",
  loci_binsize = 500,
  min_base_qual= 10,
  additional_allelecounter_flags="-f 0")
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = "XX", genomeVersion = "hg38")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_G1000_hg38.txt", replictimingfile = "RT_G1000_hg38.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = TRUE)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
output_file <- "purity_ploidy_results.txt"
file_conn <- file(output_file, "w")
writeLines(paste("Purity: ", ascat.output$purity), file_conn)
writeLines(paste("Ploidy: ", ascat.output$ploidy), file_conn)
close(file_conn)
```

## CNAnorm
```
gcbase=CNAnorm/gc1000base_38.txt
Tumor=$1
Normal=$2
prefix=$3
mkdir ${prefix} && cd ${prefix}
mkdir tmp
#Get windows
perl cox_ctdna_cnv/bam2windows.pl --samtools-path samtools -ts -cs -d ./tmp -r 1000 -gc ${gcbase} ${Tumor} ${Normal} > ${prefix}.tab.tmp
less ${prefix}.tab.tmp|grep -v '_' > ${prefix}.tab
#Run CNAnorm
source activate r4.3.1
Rscript Run_CANnorm.r ${prefix}
```
**Run_CANnorm.r**
```
args = (commandArgs(TRUE))
library(CNAnorm)
input.tab <- paste(args[[1]],".tab", sep="",collapse="")
da <- read.table(input.tab, sep = '\t', header = T)
CN <- dataFrame2object(da)
toSkip <- c("chrY","chrM")
CN <- gcNorm(CN, exclude = toSkip)
CN <- addSmooth(CN, lambda = 7)
CN <- peakPloidy(CN, exclude = toSkip)
pdf(file="PeakPloidy.pdf")
p1 <-plotPeaks(CN, special1 = 'chrX', special2 ='chrY')
print(p1)
dev.off()
CN <- addDNACopy(CN)
CN <- validation(CN)
CN <- discreteNorm(CN)
pdf(file="CNVprofile.pdf")
p2 <-plotGenome(CN, superimpose = 'DNACopy')
print(p2)
dev.off()
exportTable(CN, file = "CNAnorm_table.tab", show = 'ploidy')
file_conn <- file("tumContent.txt", "w")
writeLines(paste("Purity: ", CN@Res@suggested.tumContent), file_conn)
close(file_conn)
```

## ABSOLUTE
```
seg=$1
maf=$2
sample_name=$3
r-4.1.3/bin/Rscript RunAbsolute.r ${seg} ${sample_name}_out ${sample_name} ${maf}
```
**RunAbsolute.r**
```
args = (commandArgs(TRUE))
seg.dat <- paste(args[[1]], sep="")
out.path <- paste(args[[2]], sep="")
sp.name <- paste(args[[3]], sep="")
maf <- paste(args[[4]], sep="")
library(ABSOLUTE)
RunAbsolute(
  seg.dat.fn = seg.dat,
  maf.fn = maf,
  min.mut.af = 0,
  min.ploidy = 0.5,
  max.ploidy = 4,
  max.sigma.h = 0.2,
  platform = "Illumina_WES",
  copy_num_type = "total",
  sigma.p = 0,
  results.dir = out.path,
  primary.disease = "BLCA",
  sample.name = sp.name,
  max.as.seg.count = 1500,
  max.non.clonal = 1,
  max.neg.genome = 0.005
)
```

## absCNseq
```
source activate r4.3.1
Rscript run_absCNseq.r ${CN} ${SNP} ${Tumor_name} ${run_path}
```
**run_absCNseq.r**
```
args = (commandArgs(TRUE))
library(absCNseq)
CN.file <- paste(args[[1]], sep="")
snv.file <- paste(args[[2]], sep="")
sample.name <- paste(args[[3]], sep="")
result.dir <- paste(args[[4]], sep="")
abscnseq.output <- run.absCNSeq(CN.file,snv.fn = snv.file, result.dir,sample.name, seq.type = "WGS")
output_file <- "purity_ploidy_results.txt"
file_conn <- file(output_file, "w")
writeLines(paste("Purity: ", abscnseq.output$searchRes[1,"alpha"]), file_conn)
writeLines(paste("Ploidy: ", abscnseq.output$searchRes[1,"tau"]), file_conn)
close(file_conn)
```
## PyLOH
```
Tumor=$1
Normal=$2
Tumor_name=$3
segment=$4
mkdir ${Tumor_name} && cd ${Tumor_name}
run_path=$(pwd)
python seg2bed.py ${segment} segment.bed --seg_length 1000000
python PyLOH.py preprocess genome.fa ${Normal} ${Tumor} ${Tumor_name} --min_depth 20 --min_base_qual 10 --min_map_qual 10 --process_num 36 --segments_bed segment.bed
python PyLOH.py run_model ${Tumor_name} --allelenumber_max 2 --max_iters 100 --stop_value 1e-7
```

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
## FACETS
```
source activate cnv_facets
tumor=$1
normal=$2
prefix=$3
mkdir ${prefix} && cd ${prefix}
Rscript cnv_facets.R -t ${tumor} -n ${normal} -vcf hg38_addchr.vcf.gz -o ${prefix} --gbuild hg38
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
**Tumor only mode**
```
mkdir Tumor_rep1
java -jar -Xmx8G home/software/WGD_software/cobalt-1.14.1.jar -tumor_only_diploid_bed home/ref/Human/hg38/HMFTools/v6_0/ref/38/copy_number/DiploidRegions.38.bed.gz -tumor Tumor_rep1 -tumor_bam home/Benchmark/Alignment/Data1/Preprocess_Data1/Tumor_1.hg38.rmdup.filter.sorted.bam -gc_profile home/ref/Human/hg38/HMFTools/v6_0/ref/38/copy_number/GC_profile.1000bp.38.cnp -output_dir Tumor_rep1 -threads 24
java -jar -Xmx8G home/software/WGD_software/amber-3.9.jar -tumor Tumor_rep1 -tumor_bam home/Benchmark/Alignment/Data1/Preprocess_Data1/Tumor_1.hg38.rmdup.filter.sorted.bam -output_dir Tumor_rep1 -loci home/ref/Human/hg38/HMFTools/v6_0/ref/38/GermlineHetPon.38.vcf.gz -threads 24 -ref_genome_version V38
java -jar -Xmx8G home/software/WGD_software/hmftools-purple-v3.8.4/purple_v3.8.4.jar -tumor Tumor_rep1 -amber Tumor_rep1 -cobalt Tumor_rep1 -gc_profile home/ref/Human/hg38/HMFTools/v6_0/ref/38/copy_number/GC_profile.1000bp.38.cnp -ref_genome_version 38 -ref_genome home/WuLab/Benchmark/human/genome.fa -ensembl_data_dir home/ref/Human/hg38/HMFTools/v6_0/ref/38/common/ensembl_data -output_dir Tumor_rep1 -threads 24
```
## Accucopy
```
HOST_DIR=home/WuLab/Benchmark/Accucopy/
CONTAINER_DIR=/mnt
docker run --rm -v ${HOST_DIR}:${CONTAINER_DIR} -w /usr/local/Accucopy imagehub.cstcloud.cn/polyactis/accucopy python main.py -c ${CONTAINER_DIR}/configure -t ${CONTAINER_DIR}/Tumor_1.hg38.rmdup.filter.sorted.bam -n ${CONTAINER_DIR}/Normal_1.hg38.rmdup.filter.sorted.bam -o ${CONTAINER_DIR}/Tumor_rep1/output --snp_output_dir ${CONTAINER_DIR}/Tumor_rep1/snps
```


