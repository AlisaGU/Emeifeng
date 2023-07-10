#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library("seqinr")
library("data.table")
library("dplyr")
library("Biostrings")
library("parallel")
# 2. functions ------------------------------------------------------------ TODO:
read.gff <- function(gff_name = NULL) {
    data <- fread(cmd = paste0("awk -F\"\\t\" '$3==\"CDS\"' ", gff_name), header = F, stringsAsFactors = F)
    colnames(data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    data$attribute <- sapply(data$attribute, function(x) {
        gsub("Parent=", "", regmatches(x, regexpr("Parent=[^;]*", x)))
    })
    return(data)
}

get_sequence <- function(data = NULL) {
    cds <- paste(CDS_seq[rev(data$cdseqname)], collapse = "")
    pep <- c2s(seqinr::translate(s2c(cds)))
    result <- list(cds = cds, pep = pep)
}

wrapperForWrite.fasta <- function(seqs, outfile) {
    if (file.exists(outfile)) {
        file.remove(outfile)
    }
    namesOfSeq <- names(seqs)
    for (i in 1:length(seqs)) {
        seq <- seqs[i]
        write.fasta(s2c(seq), file.out = outfile, open = "a", names = names(seq))
    }
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
gff_name <- Args[6]
CDS_seqName<- Args[7]
cds_for_each_mRNA_outfilename <- Args[8]
pep_for_each_mRNA_outfilename <- Args[9]

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
CDS_seq <- read.fasta(
    file = CDS_seqName,
    as.string = TRUE, seqtype = "DNA", set.attributes = F
)

gff <- read.gff(gff_name = gff_name)
gff$cdseqname <- paste(paste(gff$seqname, paste(gff$start, gff$end, sep = "-"), sep = "_"), gff$strand, sep = ":")
gff_split <- split(gff, f = gff$attribute)

seq1 <- mclapply(gff_split, function(data) {
    get_sequence(data = data)
}, mc.cores = 8)

cdsSeqs <- sapply(seq1, function(x) {
    x[[1]]
})
pepSeqs <- sapply(seq1, function(x) {
    x[[2]]
})
wrapperForWrite.fasta(cdsSeqs, cds_for_each_mRNA_outfilename)
wrapperForWrite.fasta(pepSeqs, pep_for_each_mRNA_outfilename)
