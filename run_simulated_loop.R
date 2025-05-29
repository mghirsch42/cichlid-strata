library(StrataFinder)
library(ggplot2)
library(optparse)

# Generate command line flags
option_list <- list(
    make_option(c("-d", "--data_path"), action="store", type="character", default=""),
    make_option(c("-s", "--save_file"), action="store", type="character", default=""),
    make_option(c("-w", "--figure_save_path"), action="store", type="character", default=""),
    make_option(c("-f", "--fit_type"), action="store", type="character", default="linear"),
    make_option(c("-u", "--min_bps"), action="store", type="integer", default=0),
    make_option(c("-v", "--max_bps"), action="store", type="integer", default=5),
    make_option(c("-m", "--min_len"), action="store", type="integer", default=5),
    make_option(c("-e", "--early_stop"), action="store_true", default=FALSE),
    make_option(c("-q", "--quiet"), action="store_true", default=FALSE),
    make_option(c("-t", "--test"), action="store", type="character", default="t"),
    make_option("--show_error", action="store_true", default=FALSE)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Error checking and implementing defaults.
if (trimws(opt$data_path) == "") {
    print("No data path input. Please input a valid data path.")
    quit()
} else if (!dir.exists(opt$data_path)) {
    print("Data path does not exist. Please create the directory.")
}
save_path <- strsplit(opt$save_file, "/")[[1]]
save_path <- save_path[1:(length(save_path)-1)]
save_path <- paste0(save_path, "/", collapse="/")
if (trimws(opt$save_file) == "") {
    print("No save file input. Will not save results.")
    opt$save_file <- NULL
} else if (!dir.exists(save_path)) {
    print("Save path does not exist. Please create the directory.")
    quit()
}
if (trimws(opt$figure_save_path) == "") {
    print("No figure save path input. Will not save figures")
    opt$figure_save_path <- NULL
} else if (!file.exists(opt$figure_save_path)) {
    print("Figure save path does not exist. Please create the directory.")
    quit()
}
if (opt$min_bps < 0 | opt$max_bps < 0) {
    print("Number of breakpoints cannot be negative.")
    quit()
}
if (opt$test!="t" & opt$test!="f" & opt$test!="either" & opt$test!="both") {
    print("Test type must be 't', 'f', 'either', or 'both'.")
    quit()
}

writeLines(paste0("Running StrataFinder over simulated data with:\n",
                "\tData path: ", opt$data_path, "\n",
                "\tSave file: ", opt$save_file, "\n",
                "\tFigure save path: ", opt$figure_save_path, "\n",
                "\tFit type: ", opt$fit_type, "\n",
                "\tMinimum number of breakpoints: ", opt$min_bp, "\n",
                "\tMaximum number of breakpoints: ", opt$max_bp, "\n",
                "\tMinimum length of segment between breakpoints: ", opt$min_len, "\n",
                "\tTest type: ", opt$test))

sf <- new("StrataFinder")
sf$fit_type <- opt$fit_type
sf$min_len <- opt$min_len

data_files <- list.files(path=opt$data_path, full.names=TRUE)

result_df <- data.frame()


for (f in data_files) {

    sim_idx <- substr(f, nchar(opt$data_path)+2, nchar(f)-4)
    if (sim_idx == "bps") { next }
    if (as.integer(sim_idx) %% 10 == 0) {
        print(sim_idx)
    }

    df <- read.table(f, sep=",", header=TRUE)

    if (sf$fit_type == "quadratic") {
        sf$x <- df$x / 1000
    }
    else {
        sf$x <- df$x
    }
    sf$y <- df$y
    sf$precomputeSSRMat()
    curr_results <- sf$findStrata(opt$min_bps, opt$max_bps, test=opt$test, early_stop=opt$early_stop, quiet=opt$quiet)
    curr_df <- data.frame("sim_idx"=sim_idx, 
                          "n_bps"=curr_results$n_bps, 
                          "bps_idxs"=paste(unlist(curr_results$bps_idxs), collapse=","),
                          "bps"=paste(unlist(curr_results$bps), collapse=","))
    result_df <- rbind(result_df, curr_df)
    if (!is.null(opt$figure_save_path)) {
        p <- sf$plot(show_error=FALSE)
        ggsave(paste0(opt$figure_save_path, sim_idx, ".png"), width=7, height=3, p)
    }
}

print(result_df)

if (!is.null(opt$save_file)) {
    write.csv(result_df, opt$save_file, row.names=FALSE)
}