library(StrataFinder)
library(ggplot2)
library(patchwork)
library(scales)

data_file_path <- "data/fungus_gnats/Bcop.Holo2.XpX.pool.Illumina.indHetdens.100kbwins.csv" # Where my fungus is

fit_type <- "exponential"
min_len <- 5
min_bps <- 0 # Minimum number of breakpoints to check for
max_bps <- 10 # Maximum number of breakpoints to check for (currently this exits early if the p-value < 0.05)
n_bps <- 7 # If positive, will estimate this number, if 0 or less, will find optimal number

figure_save_path <- paste("figures/", "fungus_gnats/", fit_type, "/", sep="")
if (n_bps > 0) {
    result_save_file <- paste("results/", "fungus_gnats/fungus_gnats_", fit_type, "_", n_bps, ".csv", sep="")
} else {
    result_save_file <- paste("results/", "fungus_gnats/fungus_gnats_", fit_type, ".csv", sep="")
}

df_all <- read.table(data_file_path, sep=",", header=TRUE)
df_all <- df_all[!is.na(df_all$het_H2_XpX_pool_reseq.Bcop_v2),]
df_all$counts <- df_all$het_H2_XpX_pool_reseq.Bcop_v2 * 10000 # We must have whole numbers
df_all$start <- (df_all$start-1)

sf <- new("StrataFinder")
sf$fit_type <- fit_type
sf$min_len <- min_len

result_df <- data.frame()

# for (s in unique(df_all$scaffold) ){
for (s in c("X")){
    print(s)

    df <- df_all[df_all$scaffold == s,] 
    print(head(df))
    if (sf$fit_type == "quadratic") {
        sf$x <- as.integer(df$start)
    } else {
        sf$x <- as.integer(df$start)
    }
    sf$y <- as.numeric(df$counts)

    sf$precomputeSSRMat()
    if (n_bps > 0) {
        curr_results <- sf$findNBreakpoints(n_bps)
        curr_results$n_bps <- n_bps
    } else {
        curr_results <- sf$findBreakpoints(min_bps, max_bps, early_stop=TRUE, quiet=FALSE)
    }
    
    print(curr_results)
    curr_df <- data.frame("scaffold"=s, "n_bps"=curr_results$n_bps, "bps"=paste(unlist(curr_results$bps), collapse=", "))
    result_df <- rbind(result_df, curr_df)

    if (!is.null(figure_save_path)) {
        p <- sf$plot(show_error=FALSE)
        p <- p + 
                    scale_x_continuous(labels=comma, breaks=sf$x[sf$x %% 1e7 == 0]) + #xlab("Genomic position (kb)") +
                    ylim(0, 230) +
                    scale_y_continuous(breaks=c(0, 50, 100, 150, 200), labels=c("0", "0.005", "0.010", "0.015", "0.020")) #ylab("Density of heterozygous sites") +
        p <- p + theme(text=element_text(family="TT Arial")) + ylab("")
        if (n_bps > 0) {
            ggsave(paste(figure_save_path, s, "_", n_bps, ".svg", sep=""), width=7, height=2, p)
        } else {
            ggsave(paste(figure_save_path, s, ".svg", sep=""), width=7, height=3, p)
        }
    }
}
if (!is.null(result_save_file)) {
    write.csv(result_df, result_save_file, row.names=FALSE)
}