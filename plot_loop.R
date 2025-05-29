library(StrataFinder)
library(ggplot2)
library(scales)


# Load the data
species_group <- "tropheini"
kb <- 100
lg <- 19
n_bps <- 3
fit_type <- "quadratic"
data_file_path <- paste("data/", species_group, "/", kb, "kb/LG", lg, "/", sep="")
# result_save_file <- paste("results/", species_group, "/", species_group, "_", fit_type, "_LG_", lg, "_", kb, "kb_", n_bps, "bps.csv", sep="")
# figure_save_path <- paste0("figures/", species_group, "_", n_bps, "/LG", lg, "/", fit_type, "/", kb, "kb/")
figure_save_path <- paste0("figures/", species_group, "/LG", lg, "/", fit_type, "/", kb, "kb/")
min_len <- 5

data_files <- list.files(path=data_file_path, full.names=TRUE)

for (f in data_files) {
    species <- substr(f, nchar(data_file_path)+2, nchar(f)-4)
    print(species)

    data_file <- paste(data_file_path, species, ".txt", sep="")
    data_df <- read.table(data_file, sep="\t", header=FALSE)

    # result_save_file <- paste("results/", species_group, "/", species_group, "_", fit_type, "_LG_", lg, "_", kb, "kb_", n_bps, "bps.csv", sep="")
    result_save_file <- paste("results/", species_group, "/", species_group, "_", fit_type, "_LG_", lg, "_", kb, "kb.csv", sep="")
    results <- read.table(result_save_file, sep=",", header=TRUE)
    bps <- results[results$species == species,"bps"]
    if(grepl(",", bps, fixed=TRUE))
        bps <- unlist(strsplit(bps, ", "))
    bps <- as.integer(bps)

    bps_idxs <- c()
    if (sum(is.na(bps)) == 0) {
        for (i in 1:length(bps)) {
            if (fit_type == "quadratic") {
                bps_idxs <- c(bps_idxs, rownames(data_df[(data_df$V2)/1000 == bps[i],])) 
            } else {
                bps_idxs <- c(bps_idxs, rownames(data_df[data_df$V2 == bps[i],])) 
            }
        }
        bps_idxs <- as.integer(bps_idxs)
    }
    sf <- new("StrataFinder")
    sf$fit_type <- fit_type
    sf$min_len <- min_len
    if (fit_type == "quadratic") {
        sf$x <- (data_df$V2 / 1000)
    } else {
        sf$x <- data_df$V2
    }
    sf$y <- data_df$V4
    if (sum(is.na(bps)) == 0) {
        sf$bps <- bps
        sf$bps_idxs <- as.numeric(bps_idxs)
    }
    res <- sf$piecewiseRegression()
    p <- sf$plot()
    if (fit_type == "quadratic") {
        p <- p + coord_cartesian(xlim=c(0,max(data_df$V2/1000))) +
                 scale_x_continuous(breaks=c(0, 10000, 20000, 30000), labels=c("0", "10,000,000", "20,000,000", "30,000,000")) 
    } else {
    p <- p + coord_cartesian(xlim=c(0,max(data_df$V2))) +
             scale_x_continuous(breaks=c(0, 10000000, 20000000, 30000000), labels=c("0", "10,000,000", "20,000,000", "30,000,000")) 

    }
    p <- p + theme(text=element_text(family="TT Arial")) + ylab("")
            #  ggtitle(species)
    # Size for the plots to put on trees: 7 x 2.5
    ggsave(paste(figure_save_path, species, ".svg", sep=""), width=7, height=2.5, p)
}