library(StrataFinder)
library(ggplot2)

##########
# Estimates breakpoint locations for all data files in a folder.
# Optionally also generates figures for all files.
##########

# # Load the data
species_group <- "outgroups"
kb <- 100
lg <- 14
fit_type <- "exponential"
data_file_path <- paste("data/", species_group, "/", kb, "kb/LG", lg, "/", sep="")

show_error=FALSE

min_len <- 5
min_bps <- 0
max_bps <- 10
n_bps <- 0

if (n_bps > 0) {
    result_save_file <- paste("results/", species_group, "/", species_group, "_", fit_type, "_LG_", lg, "_", kb, "kb_", n_bps, "bps.csv", sep="")
    figure_save_path <- paste("figures/", species_group, "_", n_bps, "/LG", lg, "/", fit_type, "/", kb, "kb/", sep="")
} else {
    result_save_file <- paste("results/", species_group, "/", species_group, "_", fit_type, "_LG_", lg, "_", kb, "kb.csv", sep="")
    figure_save_path <- paste("figures/", species_group, "/LG", lg, "/", fit_type, "/", kb, "kb/", sep="")
}

sf <- new("StrataFinder")
sf$fit_type <- fit_type
sf$min_len <- min_len

data_files <- list.files(path=data_file_path, full.names=TRUE)

result_df <- data.frame()

for (f in data_files) {
    species <- substr(f, nchar(data_file_path)+2, nchar(f)-4)
    print(species)

    df <- read.table(f, sep="\t", header=FALSE)

    if (sf$fit_type == "quadratic") {
        sf$x <- df$V2 / 1000
    }
    else {
        sf$x <- df$V2
    }
    sf$y <- df$V4

    sf$precomputeSSRMat()

    if (n_bps > 0) {
        curr_results <- sf$findNBreakpoints(n_bps)
        curr_results$n_bps <- n_bps
    } else {
        curr_results <- sf$findBreakpoints(min_bps, max_bps, early_stop=TRUE, quiet=FALSE, test="t")
    }
    curr_df <- data.frame("species"=species, 
                          "n_bps"=curr_results$n_bps, 
                          "bps_idxs"=paste(unlist(curr_results$bps_idxs), sep=", ", collapse=", "),
                          "bps"=paste(unlist(curr_results$bps), sep=", ", collapse=", "))
    result_df <- rbind(result_df, curr_df)

    if (!is.null(figure_save_path)) {
        p <- sf$plot(show_error=show_error)
        # p <- p + ggtitle(species)
        p <- p + scale_y_continuous(breaks=c(0, 200, 400, 600, 800, 1000, 1200))
        ggsave(paste(figure_save_path, species, ".svg", sep=""), width=7, height=3, p)
    }
}

print(result_df)
write.csv(result_df, result_save_file, row.names=FALSE)