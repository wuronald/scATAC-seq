# module load R/4.1.0

library(DiffBind) # DiffBind_3.2.1

# loading previous DiffBind analysis
load("2023-03-27_all_normoxia_vs_hypoxia_DiffBind.RData")

normread <- dba.peakset(res.cnt, bRetrieve = TRUE, DataType = DBA_DATA_FRAME, writeFile = 'normread.txt')
datac    <- dba.count(res.cnt, peaks=NULL, score=DBA_SCORE_READS)
rawread  <- dba.peakset(datac, bRetrieve = TRUE, DataType = DBA_DATA_FRAME, writeFile = 'rawread.txt' )

# extract colnames
sample_cols <- paste(res.cnt$samples$Tissue, res.cnt$samples$Treatment, res.cnt$samples$Replicate, sep = "_")

# adding appropriate colnames to raw read table
colnames(rawread)[(1:3)] <- tolower(colnames(rawread)[(1:3)]) # lower case chr, start, end
colnames(rawread)[-(1:3)] <- sample_cols # add sample names to colnames

# Use GenomicRanges to convert the peak names to GRanges object
library(GenomicRanges) #  GenomicRanges_1.44.0 

df <- rawread   
bulk_gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end)) # DiffBind uses 1-based coordinates
df <- df[, which(!(colnames(df) %in% c("chr", "start", "end")))]
rownames(df) <- paste(bulk_gr)

# convert coordinates to hg38
library(liftOver)
library(rtracklayer)

# download chain file from UCSC:
# save to /cluster/home/rwu/R/x86_64-pc-linux-gnu-library/4.1/liftOver/extdata
# wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
# gunzip hg19ToHg38.over.chain.gz
# import chain file provided by UCSC:

path1 <- system.file(package="liftOver","extdata","hg19ToHg38.over.chain")
ch = import.chain(path1)

# perform liftOver (convert hg19 to hg38)
bulk_gr_hg38 <- liftOver(bulk_gr, ch)
bulk_gr_hg38 <- unlist(bulk_gr_hg38) # convert CompressedGRangesList to GRanges
genome(bulk_gr_hg38) <- "hg38"

# IMPORTANT: Some regions may not lift over successfully
# Keep only the rows that successfully lifted over to exactly one region (0 or >1 are discarded)
keep_idx <- which(lengths(liftOver(bulk_gr, ch)) == 1)
df <- df[keep_idx, ]
bulk_gr_hg38 <- bulk_gr_hg38[keep_idx]

# subset df to remove coordinate columns and add hg38-based rownames
df <- df[, which(!(colnames(df) %in% c("chr", "start", "end")))]
rownames(df) <- paste(bulk_gr_hg38)

# Create SummarizedExperiment Object
seBulk <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(df)), rowRanges = bulk_gr_hg38)
seBulk

# save SummarizedExperiment object
saveRDS(seBulk, file = "2023-03-27_seBulk_hg38_all_normoxia_vs_hypoxia_DiffBind.rds")
