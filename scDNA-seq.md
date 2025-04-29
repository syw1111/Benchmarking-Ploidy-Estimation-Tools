## HMMcopy
```
mkdir single_cell_pipeline
cd  single_cell_pipeline
docker run -w $PWD -v home/WuLab/Benchmark/HMMcopy_final/single_cell_pipeline/:home/WuLab/Benchmark/HMMcopy_final/single_cell_pipeline/ quay.io/singlecellpipeline/single_cell_pipeline_hmmcopy:v0.8.14 sh TN1.sh
```
**TN1.sh**
```
single_cell hmmcopy \
        --input_yaml home/WuLab/Benchmark/HMMcopy_final/single_cell_pipeline/UID-NNA-TN1.yaml \
        --library_id NAVINACT --maxjobs 8 --nocleanup \
        --sentinel_only --loglevel DEBUG  --submit local \
        --tmpdir tmp --pipelinedir pipeline --out_dir output_TN1/ \
        --config_override '{"refdir": "home/WuLab/Benchmark/HMMcopy_final/single_cell_pipeline/mondrian-ref-20-22/",  "hmmcopy": {"gc_wig_file": "home/WuLab/Benchmark/HMMcopy_final/single_cell_pipeline/mondrian-ref-20-22/human/GRCh37-lite.gc.ws_500000.wig", "map_wig_file": "home/WuLab/Benchmark/HMMcopy_final/single_cell_pipeline/mondrian-ref-20-22/human/GRCh37-lite.map.ws_125_to_500000.wig", "ref_genome": "home/WuLab/Benchmark/HMMcopy_final/single_cell_pipeline/mondrian-ref-20-22/human/hg19.fa"}}'
```
**UID-NNA-TN1.yaml**
```
UID-NNA-TN1_SLX-NAVINACT_000001_TN1-S3-C298-SRR13981519:
  bam: home/WuLab/Benchmark/HMMcopy_final/single_cell_pipeline/data/TN1/UID-NNA-TN1_SLX-NAVINACT_000001_TN1-S3-C298-SRR13981519.bam
  is_control: no
  condition: B
  pick_met: C1
  column: 1
  row: 1
  img_col: 1
  index_i5: i5-20
  index_i7: i7-28
  primer_i5: GTATAG
  primer_i7: CTATCT
  sample_id: UID-NNA-TN1
  library_id: SLX-NAVINACT
```
## Ginkgo
```
source activate ginkgo
sample_name=$1
bed_path=$2
mkdir ${sample_name}_output && cd ${sample_name}_output
run_path=$(pwd)
mkdir bed_data && cd bed_data
ln -s ${bed_path}/*.gz ./
cd ${run_path}
ls ${run_path}/bed_data | grep .bed.gz$ > ${run_path}/bed_data/list
ulimit -n 4096
bash ~/analyze.sh ${run_path}/bed_data config.example
```
## AneuFinder
```
bams_path=$1
project_name=$2
source activate R-4.4.0
mkdir ${project_name} && cd ${project_name}
run_path=$(pwd)
mkdir BAM_PATH && cd BAM_PATH && ln -s ${bams_path}/*.bam ./ && cd ${run_path}
Rscript run_aneufinder.R -b ${run_path}/BAM_PATH  -d ${run_path} -w 100e3
python get_aneufinderCN.py ${run_path}/BROWSERFILES/method-edivisive/binsize_1e+05_stepsize_1e+05_CNV.bed.gz aneufinder_cnv.tsv 100000
sample_num=$(less aneufinder_cnv.tsv|awk '{print NF}'|head -n1)
for i in $(seq 4 $sample_num)
do
sample_name=$(less aneufinder_cnv.tsv | grep 'CHROM' | awk -v col=$i '{print $col}')
ploidy=$(less aneufinder_cnv.tsv | grep -v 'CHROM' | awk -v col=$i '{sum += $col} END {print sum/NR}')
echo -e "${sample_name}\t${ploidy}" >> result.txt
done
```


## SCOPE
```
source activate R-4.4.0
Rscript Run.R
```
**Run.R**
```
library(SCOPE)
library(WGSmapp)
bamfolder <- "home/WuLab/Benchmark/SRA018951/BAM"
bamFile <- list.files(bamfolder, pattern = "*.hg19.sorted.rg.dedup.bam$", full.names = TRUE)
sampname_raw <- sapply(basename(bamFile), function(x) {
                  sub("\\.hg19\\.sorted\\.rg\\.dedup\\.bam$", "", x)  
})
sampname <- unname(sapply(strsplit(sampname_raw, ".", fixed = TRUE), "[", 1))
bambedObj <- get_bam_bed(bamdir = bamFile, sampname = sampname)
ref_raw <- bambedObj$ref
mapp <- get_mapp(ref_raw)
gc <- get_gc(ref_raw)
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 20, seq = 'single-end')
Y_raw <- coverageObj$Y
QCmetric_raw <- get_samp_QC(bambedObj)
qcObj <- perform_qc(Y_raw = Y_raw, sampname_raw = sampname, ref_raw = ref_raw, QCmetric_raw = QCmetric_raw)
Y <- qcObj$Y
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric
Gini <- get_gini(Y)
normObj.sim <- normalize_codex2_ns_noK(Y_qc = Y,gc_qc = ref$gc,norm_index = which(Gini<=0.12))
pdf("SRA018951_plot.pdf")
ploidy.sim <- initialize_ploidy(Y = Y,Yhat = normObj.sim$Yhat, ref = ref, SoS.plot=TRUE)
dev.off()
ploidy.sim
```
## SeCNV
```
sample_name=$1
bam_path=$2
source activate python3.9
mkdir ${sample_name}_output
cd ${sample_name}_output
run_path=$(pwd)
ref_file=hg19.fa
bam_pattern="*dedup.bam"
bin_size=500000
ref_type=hg19
min_ploidy=1.5
max_ploidy=5.0
topK="auto_set"
sigma="auto_set"
normal_cell="None"
mkdir input_tmp
cp ${bam_path}/*.bam* input_tmp/
i=1
while [ $(ls ${run_path}/input_tmp/* | wc -l) -gt 0 ]
do
mkdir ${sample_name}_input_${i} && mkdir ${sample_name}_output_${i}
ls ${run_path}/input_tmp/*.bam|head -n100|while read line
do
mv ${line}* ${sample_name}_input_${i}/
done
output_path_sub=${run_path}/${sample_name}_output_${i}
bam_path_sub=${run_path}/${sample_name}_input_${i}
#Step 1
python preprocess.py ${output_path_sub} ${ref_file} ${bam_path_sub} ${bam_pattern} ${bin_size} ${ref_type}
#Step 2
python call_cn.py ${output_path_sub} ${min_ploidy} ${max_ploidy} ${topK} ${sigma} ${normal_cell}
i=$((i+1))
done
```
## rcCAE
```
source activate rccae
mappability_file=hg19_50.bw
hg19_ref=hg19.fa
bams_dir=$1
sample_name=$2
mkdir ${sample_name}_output && cd ${sample_name}_output
run_path=$(pwd)
ls ${bams_dir}/*.dedup.bam|awk -F '/' '{print $NF}'|sed 's/.bam//g' > bams_name.txt
# Step1
bams_name=$(ls ${run_path}/bams_name.txt)
/sibcb1/wuweilab1/songyawei/software/WGD_scDNAseq_software/rccae/prep/bin/prepInput --bam ${bams_dir} --ref ${hg19_ref} --map ${mappability_file} --barcode ${bams_name} --output prep.txt --chrlist chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22
# Step2
python rccae/cae/train.py --input prep.txt --output train_result
#Step3
/sibcb1/wuweilab1/songyawei/software/WGD_scDNAseq_software/rccae/hmm/run_SCHMM.sh MCR_R2016b_glnxa64_installer/MUTsigcv/v91 train_result/lrc.txt result 10
```
**mappability_file generation**
```
gemlibrary=GEM-binaries-Linux-x86_64-core_2-20130406-045632
bigwig=UCSC
hg19_ref=hg19.fa
hg19_chromsize=hg19.chrom.sizes.txt
chmod +x $gemlibrary/bin/gem* $bigwig/wigToBigWig
export PATH=$PATH:$gemlibrary/bin:$bigwig
gem-indexer -T 10 -c dna -i ${hg19_ref} -o hg19_index
gem-mappability -T 10 -I hg19_index.gem -l 50 -o hg19_50
gem-2-wig -I hg19_index.gem -i hg19_50.mappability -o hg19_50
wigToBigWig hg19_50.wig ${hg19_chromsize} hg19_50.bw
mappability_file=$(pwd)/hg19_50.bw
```


## CNVeil
```
source activate CNVeil
python3 home/software/WGD_scDNAseq_software/CNVeil-main/Python/CNVeil.py -bam home/WuLab/Benchmark/SRA018951/BAM --reference home/Reference/Human/hg19/hg19.fa --reftype hg19 --seq_type single-end --output_dir output -t 10 -crm q -pm cl
```
## scAbsolute
```
source activate Copynumber
snakemake --snakefile workflow/Snakefile_absolute --use-conda
```
