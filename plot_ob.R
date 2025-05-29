library(StrataFinder)
library(ggplot2)
library(scales)
library(gridExtra)
library(patchwork)

fit_type <- "linear"
figure_save_path <- "figures/paper_figures/ob/"

run <- function(sf, data_df, result_file, species, fit_type) {
    results <- read.table(result_file, sep=",", header=TRUE)
    bps <- get_bps(results, species)
    bps_idxs <- get_bps_idxs(data_df, bps)
    if (fit_type == "quadratic") { sf$x <- (data_df$V2 / 1000)
    } else { sf$x <- data_df$V2 }
    sf$y <- data_df$V4
    sf$bps <- bps
    sf$bps_idxs <- as.numeric(bps_idxs)
    res <- sf$piecewiseRegression()
    return(sf)
}

get_bps <- function(results, species) {
    bps <- results[results$species == species,"bps"]
    if(grepl(",", bps, fixed=TRUE))
        bps <- unlist(strsplit(bps, ","))
    bps <- as.integer(bps)
    return(bps)
}
get_bps_idxs <- function(data_df, bps) {
    bps_idxs <- c()
    for (i in 1:length(bps)) {
        if (fit_type == "quadratic") { bps_idxs <- c(bps_idxs, rownames(data_df[(data_df$V2)/1000 == bps[i],])) 
        } else { bps_idxs <- c(bps_idxs, rownames(data_df[data_df$V2 == bps[i],]))}
    }
    bps_idxs <- as.integer(bps_idxs)
    return(bps_idxs)
}

plotter <- function(sf, data_df, fit_type) {
    p <- sf$plot(show_error=TRUE)
    if (fit_type == "quadratic") {
        p <- p + coord_cartesian(xlim=c(0,max(data_df$V2/1000))) +
                    scale_x_continuous(breaks=c(0, 10000, 20000, 30000, 40000), labels=c("0", "10,000,000", "20,000,000", "30,000,000", "40,000,000")) 
    } else {
        p <- p + coord_cartesian(xlim=c(0,max(data_df$V2))) +
                    scale_x_continuous(breaks=c(0, 10000000, 20000000, 30000000, 40000000), labels=c("0", "10,000,000", "20,000,000", "30,000,000", "40,000,000")) 
    }
    # p <- p + theme(text=element_text(family="TT Arial")) + ylab("")
    return(p)
}

# Labeotropheus_Ltrewavasae_Maison_BBM_OBF
sf <- new("StrataFinder")
sf$fit_type <- fit_type

species <- "Labeotropheus_Ltrewavasae_Maison_BBM_OBF"
data_file_path <- "data/ob/100kb/LG5/Labeotropheus_Ltrewavasae_Maison_BBM_OBF.txt"
data_llm <- read.table(data_file_path, sep="\t", header=FALSE)

sf <- run(sf, data_llm, paste0("results/ob/ob_", fit_type, "_LG_5_100kb_2bps.csv"), "Labeotropheus_Ltrewavasae_Maison_BBM_OBF", fit_type)
p_llm <- plotter(sf, data_llm, fit_type)
p_llm <- p_llm + plot_annotation(title="Labeotropheus Ltrewavasae Maison Reef BBM x OBF")
ggsave(paste0(figure_save_path, "Labeotropheus_Ltrewavasae_Maison_BBM_OBF.svg"), width=5, height=3.5, p_llm)


# Metriaclima Mzebra NkhataBay M OBF
sf <- new("StrataFinder")
sf$fit_type <- fit_type

species <- "Metriaclima_Mzebra_NkhataBay_M_OBF"
data_file_path <- "data/ob/100kb/LG5/Metriaclima_Mzebra_NkhataBay_M_OBF.txt"
data_mz <- read.table(data_file_path, sep="\t", header=FALSE)

sf <- run(sf, data_mz, paste0("results/ob/ob_", fit_type, "_LG_5_100kb_2bps.csv"), "Metriaclima_Mzebra_NkhataBay_M_OBF", fit_type)
p_mz <- plotter(sf, data_mz, fit_type)
p_mz <- p_mz + plot_annotation(title="Metriaclima Zebra Nkhata Bay BBM x OBF")
ggsave(paste0(figure_save_path, "Metriaclima_Mzebra_NkhataBay_M_OBF.svg"), width=5, height=3.5, p_mz)
