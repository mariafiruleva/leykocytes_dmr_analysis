source("https://bioconductor.org/biocLite.R")
if (!require("BiocManager")) {
  install.packages("BiocManager") }
BiocManager::install("DMRcate")
BiocManager::install("rtracklayer")

library(DMRcate)
library(bsseq)

bismarkBSseq_1 <-BiocGenerics::combine(
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_106lym_106lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_116lym_116lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_17lym_17lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_196lym_196lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_215lym_215lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_234lym_234lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_264lym_264lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_266lym_266lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_279lym_279lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_285lym_285lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_287lym_287lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_308lym_308lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_330lym_330lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_352lym_352lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_387lym_387lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_390lym_390lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_398lym_398lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_423lym_423lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_426lym_426lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_44lym_44lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_457lym_457lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_465lym_465lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_469lym_469lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_488lym_488lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_490lym_490lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_497lym_497lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_66lym_66lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_72lym_72lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_lymphocytes_74lym_74lym.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_second_batch_lymphocytes_243_lym_243_lym_dedup.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_second_batch_lymphocytes_396_lym_396_lym_dedup.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_second_batch_lymphocytes_132_lym_132_lym_dedup.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_second_batch_lymphocytes_226_lym_226_lym_dedup.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_second_batch_lymphocytes_463_limf_463_limf_dedup.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_second_batch_lymphocytes_478_lym_478_lym_dedup.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE),
  bsseq::read.bismark(files = "results_rrbs_second_batch_lymphocytes_492_lym_492_lym_dedup.bedgraph.gz.bismark.cov.gz",
                      rmZeroCov = TRUE,
                      strandCollapse = FALSE,
                      verbose = TRUE)
)

covered_loci <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq_1, type="Cov")==0) == 0)
bismarkBSseq_1 <- bismarkBSseq_1[covered_loci, ]

covariates <- read.table("smoking6_binary.tsv", sep = "\t", header = T) 

#smoking6_binary
model = model.matrix(~smoking6_binary,data=covariates) 
design = as.data.frame(covariates$smoking6_binary)
colnames(design) <- "smoking6_binary"

DMLfit = DMLfit.multiFactor(bismarkBSseq_1, design=design, formula=~smoking6_binary)
colnames(DMLfit$X)
DMLtest.smoking6_binary = DMLtest.multiFactor(DMLfit, coef=2)

DMLtest.smoking6_binary <- DMLtest.smoking6_binary[,c(3,1,2,4,5)]
colnames(DMLtest.smoking6_binary) <- c("stat", "chr", "pos", "diff", "fdr")

wgbsannot <- cpg.annotate("sequencing", DMLtest.smoking6_binary, design = model)
wgbs.DMRs<- dmrcate(wgbsannot, lambda = 1000, C = 50, pcutoff = 0.05, mc.cores = 1)
wgbs.ranges <- extractRanges(wgbs.DMRs, genome = "hg38")

nrow(as.data.frame(wgbs.ranges)) #how many sig DMRs are in this model
write.table(wgbs.ranges, "smoking6_binary_DMRcate.tsv", sep = "\t", row.names = F)

##DSS instead of DMRcate for the DMRs calling
#smoking6_binary
model = model.matrix(~smoking6_binary,data=covariates) 
design = as.data.frame(covariates$smoking6_binary)
colnames(design) <- "smoking6_binary"

DMLfit = DMLfit.multiFactor(bismarkBSseq_1, design=design, formula=~smoking6_binary)
colnames(DMLfit$X)
DMLtest.smoking6_binary = DMLtest.multiFactor(DMLfit, coef=2)

DSS_dmrs <- callDMR(DMLtest.smoking6_binary, p.threshold=0.05)
write.table(DSS_dmrs, "smoking_last6months_DSS.tsv", sep = "\t", row.names = F)
