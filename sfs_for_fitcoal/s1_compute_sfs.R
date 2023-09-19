#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(if (!require("data.table")) install.packages("data.table")))
suppressWarnings(suppressMessages(if (!require("stringr")) install.packages("stringr")))
suppressWarnings(suppressMessages(if (!require("parallel")) install.packages("parallel")))
suppressWarnings(suppressMessages(if (!require("ggplot2")) install.packages("ggplot2")))
# 2. functions ------------------------------------------------------------ TODO:
convert <- function(x) {
    one_index <- which(x == "1")
    zero_index <- which(x == "0")
    x[one_index] <- 0
    x[zero_index] <- 1
    return(x)
}

add_missing_sample_count <- function(vcf = NULL) {
    missing_count <- mclapply(vcf[, 3:ncol(vcf)], function(x) {
        str_count(x, fixed("."))
    }, mc.cores = 5)
    missing_count <- do.call(cbind, missing_count)
    vcf$missing_chro_count <- rowSums(missing_count)
    return(vcf)
}

summary_missing <- function(vcf_missing_site_count = NULL) {
    missing_count <- vcf_missing_site_count$missing_chro_count
    fre <- table(missing_count)
    ratio <- fre / sum(fre)
    cum_ratio <- cumsum(ratio)
    result <- rbind(fre, ratio, cum_ratio)
    return(result)
}

compute_SFS_under_specified_threshold <- function(vcf_missing_site_count = NULL, max_missing_site = NULL) {
    index_of_specified_threshold <- which(vcf_missing_site_count$missing_chro_count == max_missing_site)
    genotypes_under_specified_threshold <- vcf_missing_site_count[index_of_specified_threshold, -c(1, 2, ncol(vcf_missing_site_count)), with = FALSE]
    genotypes_under_specified_threshold <- sapply(genotypes_under_specified_threshold, as.numeric)
    # derived_allele_count <- mclapply(genotypes_under_specified_threshold, function(x) {
    #     str_count(x, "1")
    # }, mc.cores = 5)
    # derived_allele_count <- do.call(cbind, derived_allele_count)
    # derived_allele_total_count <- rowSums(derived_allele_count)
    derived_allele_total_count <- rowSums(genotypes_under_specified_threshold)
    result <- table(derived_allele_total_count)
    return(result)
}

convert_SFS_to_FitCoal_no_missing_mode <- function(SFS = NULL, outfile = NULL) {
    SFS <- unlist(SFS)
    cat(SFS, "\n", file = outfile, sep = "\t", append = T)
}

convert_SFS_to_FitCoal_missing_mode <- function(SFS = NULL, outfile = NULL) {
    sapply(SFS, function(x) {
        result <- c("SampleSize", length(x) + 1, x)
        cat(result, "\n", file = outfile, sep = "\t", append = T)
    })
}

transfer_SFS_to_ggplot_format <- function(SFS = NULL) {
    data <- data.frame(1:length(SFS), SFS)
    colnames(data) <- c("epsilon", "value")
    return(data)
}

plot_SFS <- function(SFS = NULL, title = NULL) {
    data <- transfer_SFS_to_ggplot_format(SFS = SFS)
    p <- ggplot(data = data) +
        geom_bar(aes(x = epsilon, y = value), stat = "identity", fill = "#76c2da", color = "#407bc9") +
        labs(title = title) +
        # geom_vline(xintercept = 36, color = "red", linetype = "dotted") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "#f3852b", face = "bold",
                size = 18, vjust = 1, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            strip.background = element_rect(fill = "#83a7e9d0"),
            strip.text = element_text(colour = "black", size = 15)
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:

Args <- commandArgs()
data_filename <- Args[6]

setwd(dirname(data_filename))
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:

data_anc <- fread(data_filename, header = F, stringsAsFactors = F)

vcf <- data_anc[, -c(3, 4)]

colnames(vcf) <- c("chrom", "pos", paste("s", 1:(ncol(vcf) - 2), sep = "_"))

vcf_missing_site_count <- add_missing_sample_count(vcf = vcf)
missing_summary <- summary_missing(vcf_missing_site_count = vcf_missing_site_count)
write.table(missing_summary, "missing_summary.txt", quote = F, sep = "\t")


max_missing_sites <- as.numeric(colnames(missing_summary))
SFS_under_specified_threshold <- lapply(max_missing_sites, function(max_missing_site) {
    SFS <- compute_SFS_under_specified_threshold(vcf_missing_site_count = vcf_missing_site_count, max_missing_site = max_missing_site)
    index_removed <- unlist(sapply(c(0, ncol(vcf) - 2 - max_missing_site), function(x) {
        which(x == names(SFS))
    }))
    SFS <- SFS[-index_removed]

    result <- rep(0, ncol(vcf) - 2 - max_missing_site - 1)
    names(result) <- 1:(ncol(vcf) - 2 - max_missing_site - 1)
    result[match(names(SFS), names(result))] <- SFS
    return(result)
})
names(SFS_under_specified_threshold) <- max_missing_sites
pdf(file = paste0(paste(max_missing_sites, collapse = "_"), ".pdf"), width = 7, height = 7)
for (ith in max_missing_sites) {
    SFS.i <- SFS_under_specified_threshold[[as.character(ith)]]
    p <- plot_SFS(SFS = SFS.i, title = as.character(ith))
    plot(p)
}
garbage <- dev.off()

outfile <- "unfolded_sfs.txt"
if (file.exists(outfile)) {
    file.remove(outfile)
}
if (length(max_missing_sites) == 1) {
    convert_SFS_to_FitCoal_no_missing_mode(SFS = SFS_under_specified_threshold, outfile = outfile)
} else {
    convert_SFS_to_FitCoal_missing_mode(SFS = SFS_under_specified_threshold, outfile = outfile)
}
