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
