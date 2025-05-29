library(StrataFinder)
library(ggplot2)

data_file_path <- "data/sticklebacks/ENSEMBL_X_to_Y_v2.csv" # Stickelback

fit_type <- "exponential"
min_len <- 5

single <- TRUE
min_bps <- 0 # Minimum number of breakpoints to check for
max_bps <- 10 # Maximum number of breakpoints to check for (currently this exits early if the p-value < 0.05)
n_bps <- 2

figure_save_file <- paste0("figures/sticklebacks/", "sticklebacks_", fit_type, ".svg")
result_save_file <- paste0("results/sticklebacks/", "sticklebacks_", fit_type, ".csv")

df <- read.table(data_file_path, sep=",", header=TRUE)
df <- df[df$ds < 0.5,]
df <- df[order(df$HiC_XPosition1),]
if (fit_type == "quadratic") {
    x <- as.integer(df$HiC_XPosition1/1000)
} else {
    x <- as.integer(df$HiC_XPosition1)
}
y <- as.numeric(df$ds*10000)

sf <- new("StrataFinder")
sf$fit_type <- fit_type
sf$min_len <- min_len
sf$x <- x
sf$y <- y

sf$precomputeSSRMat()

if (single) {
    results <- sf$findNBreakpoints(n_bps)
    results <- data.frame("bps"=gsub("^c\\(|\\)$", "", gsub("integer\\(0\\)", "", paste(results$bps, sep=","))))
    figure_save_file <- paste0(substring(figure_save_file, 1, nchar(figure_save_file)-4), "_", n_bps, ".svg")
    result_save_file <- paste0(substring(result_save_file, 1, nchar(result_save_file)-4), "_", n_bps, ".csv")
} else {
    results <- sf$findBreakpoints(min_bps, max_bps, early_stop=TRUE, quiet=FALSE) # Change early_stop to TRUE if you want to evaluate all bp number possibilities
    results <- data.frame("n_bps"=results$n_bps, "bps"=gsub("^c\\(|\\)$", "", gsub("integer\\(0\\)", "", paste(results$bps, sep=","))))
}
print(results)

if (!is.null(figure_save_file)) {
    p <- sf$plot(show_error=FALSE)
    title <- paste("Sticklebacks, ", sf$fit_type, " fits", sep="")
    if (single) {
        title <- paste(title, ", ", n_bps, " bps", sep="")
    }
    if (fit_type == "quadratic") {
        p <- p + coord_cartesian(xlim=c(0,max(df$HiC_XPosition1/1000))) +
                 scale_x_continuous(breaks=c(5000, 10000, 15000, 20000), labels=c("5,000,000", "10,000,000", "15,000,000", "20,000,000")) 
    } else {
        p <- p + coord_cartesian(xlim=c(0,max(df$HiC_XPosition1))) +
                 scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000), labels=c("5,000,000", "10,000,000", "15,000,000", "20,000,000")) 

    }
    p <- p + ggtitle(title)
    p <- p + theme(text=element_text(family="TT Arial")) + ylab("")
    ggsave(figure_save_file, width=7, height=3, p)
}

if (!is.null(result_save_file)) {
    write.csv(results, result_save_file, row.names=FALSE)
}