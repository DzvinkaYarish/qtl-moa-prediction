l = '~/local/R_libs/'

library("GenomicRanges", lib=l)
library("GenomicFeatures", lib=l)
library("dplyr", lib=l)
library("tidyr", lib=l)


args <- commandArgs(trailingOnly = TRUE)
fpath <- args[1]
outpath <- args[2]

vars <- read.csv(fpath, sep='\t')

vars <- separate(data = vars, col = variant, into = c("chr", "pos", "ref", "alt"), sep = "_")
vars <- vars %>%
  mutate(chr = sub(paste0("^", "chr"), "", chr)) %>%
  mutate(pos = as.integer(pos))

vars_ranges <- GRanges(seqnames = vars$chr, strand = c("*"),
              ranges = IRanges(start = vars$pos, width = 1))

txdb <- makeTxDbFromGFF("~/Desktop/Thesis/data/gencode.v39.annotation.nochr.gtf", format="gtf")


# txdb <- makeTxDbFromGFF("~/Desktop/Thesis/data/gencode.v44.basic.annotation.gtf", format="gtf")

# txdb <- makeTxDbFromGFF("/home/dzvinka/Desktop/Thesis/data/Homo_sapiens.GRCh38.105.gtf", format="gtf")

exons_gencode <- exons(txdb)
genes_gencode <- genes(txdb)

nearest_exons <- nearest(vars_ranges, exons_gencode)

min_distances <- c()

for (row in 1:nrow(vars)) {
    var_indx <- rownames(vars[row,])
    ex_indx <- nearest_exons[row]

    s <- abs(start(ranges(exons_gencode[ex_indx, ])) - vars[row, "pos"])
    e <- abs(end(ranges(exons_gencode[ex_indx, ])) - vars[row, "pos"])
    min_distances <- append(min_distances, min(s, e))

}

vars$dist_closest_junction <- min_distances

overlaps <- findOverlaps(vars_ranges, genes_gencode)

vars$gene_id <- NaN

vars[queryHits(overlaps),]$gene_id <- genes_gencode[subjectHits(overlaps),]$gene_id

vars$is_in_gene <- FALSE
vars[unique(queryHits(overlaps)),]$is_in_gene <- TRUE

write.table(vars, outpath, row.names=FALSE, sep="\t", quote = FALSE)
