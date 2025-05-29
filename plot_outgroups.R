library(StrataFinder)
library(ggplot2)
library(scales)
library(gridExtra)
library(patchwork)


fit_type <- "exponential"
figure_save_path <- "figures/paper_figures/outgroups/"

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
    p <- sf$plot()
    if (fit_type == "quadratic") {
        p <- p + coord_cartesian(xlim=c(0,max(data_df$V2/1000))) +
                    scale_x_continuous(breaks=c(0, 10000, 20000, 30000), labels=c("0", "10,000,000", "20,000,000", "30,000,000")) 
    } else {
        p <- p + coord_cartesian(xlim=c(0,max(data_df$V2))) +
                    scale_x_continuous(breaks=c(0, 10000000, 20000000, 30000000), labels=c("0", "10,000,000", "20,000,000", "30,000,000")) 
    }
    p <- p + theme(text=element_text(family="TT Arial", size=6)) + ylab("")
    return(p)
}

# Chromidotilapia guntheri LG 19
sf <- new("StrataFinder")
sf$fit_type <- fit_type

species <- "Chromidotilapia guntheri"
data_file_path <- "data/outgroups/100kb/LG19/Chromidotilapia_guntheri.txt"
data_cg <- read.table(data_file_path, sep="\t", header=FALSE)

sf <- run(sf, data_cg, paste0("results/outgroups/outgroups_", fit_type, "_LG_19_100kb_1bps.csv"), "Chromidotilapia_guntheri", fit_type)
p_cg_1 <- plotter(sf, data_cg, fit_type) + ylim(-2, max(data_cg$V4+5))
ggsave(paste0(figure_save_path, "Chromidotilapia_guntheri_1bps.svg"), width=5, height=3, p_cg_1)

sf <- run(sf, data_cg, paste0("results/outgroups/outgroups_", fit_type, "_LG_19_100kb_2bps.csv"), "Chromidotilapia_guntheri", fit_type)
p_cg_2 <- plotter(sf, data_cg, fit_type) + ylim(-2, max(data_cg$V4+5))
ggsave(paste0(figure_save_path, "Chromidotilapia_guntheri_2bps.svg"), width=5, height=3, p_cg_2)

sf <- run(sf, data_cg, paste0("results/outgroups/outgroups_", fit_type, "_LG_19_100kb_3bps.csv"), "Chromidotilapia_guntheri", fit_type)
p_cg_3 <- plotter(sf, data_cg, fit_type) + ylim(-2, max(data_cg$V4+5))
ggsave(paste0(figure_save_path, "Chromidotilapia_guntheri_3bps.svg"), width=5, height=3, p_cg_3)

sf <- run(sf, data_cg, paste0("results/outgroups/outgroups_", fit_type, "_LG_19_100kb_4bps.csv"), "Chromidotilapia_guntheri", fit_type)
p_cg_4 <- plotter(sf, data_cg, fit_type) + ylim(-2, max(data_cg$V4+5))
ggsave(paste0(figure_save_path, "Chromidotilapia_guntheri_4bps.svg"), width=5, height=3, p_cg_4)

p_cg <- p_cg_1 / p_cg_2 / p_cg_3 / p_cg_4
p_cg <- p_cg + plot_annotation(title="Chromidotilapia guntheri")


# Oreochromis niloticus Amherst LG 1
sf <- new("StrataFinder")
sf$fit_type <- fit_type

species <- "Oreochromis_niloticus_Amherst"
data_file_path <- "data/outgroups/100kb/LG1/Oreochromis_niloticus_Amherst.txt"
data_on <- read.table(data_file_path, sep="\t", header=FALSE)

sf <- run(sf, data_on, paste0("results/outgroups/outgroups_", fit_type, "_LG_1_100kb_1bps.csv"), "Oreochromis_niloticus_Amherst", fit_type)
p_on_1 <- plotter(sf, data_on, fit_type) + ylim(-2, max(data_on$V4+5))
ggsave(paste0(figure_save_path, "Oreochromis_niloticus_Amherst_1bps.svg"), width=5, height=3, p_on_1)

sf <- run(sf, data_on, paste0("results/outgroups/outgroups_", fit_type, "_LG_1_100kb_2bps.csv"), "Oreochromis_niloticus_Amherst", fit_type)
p_on_2 <- plotter(sf, data_on, fit_type) + ylim(-2, max(data_on$V4+5))
ggsave(paste0(figure_save_path, "Oreochromis_niloticus_Amherst_2bps.svg"), width=5, height=3, p_on_2)

sf <- run(sf, data_on, paste0("results/outgroups/outgroups_", fit_type, "_LG_1_100kb_3bps.csv"), "Oreochromis_niloticus_Amherst", fit_type)
p_on_3 <- plotter(sf, data_on, fit_type) + ylim(-2, max(data_on$V4+5))
ggsave(paste0(figure_save_path, "Oreochromis_niloticus_Amherst_3bps.svg"), width=5, height=3, p_on_3)

sf <- run(sf, data_on, paste0("results/outgroups/outgroups_", fit_type, "_LG_1_100kb_4bps.csv"), "Oreochromis_niloticus_Amherst", fit_type)
p_on_4 <- plotter(sf, data_on, fit_type) + ylim(-2, max(data_on$V4+5))
ggsave(paste0(figure_save_path, "Oreochromis_niloticus_Amherst_4bps.svg"), width=5, height=3, p_on_3)

p_on <- p_on_1 / p_on_2 / p_on_3 / p_on_4
p_on <- p_on + plot_annotation(title="Oreochromis niloticus Amherst")


# Oreochromis mossambicus LG 14
sf <- new("StrataFinder")
sf$fit_type <- fit_type

species <- "Oreochromis_mossambicus"
data_file_path <- "data/outgroups/100kb/LG14/Oreochromis_mossambicus.txt"
data_on <- read.table(data_file_path, sep="\t", header=FALSE)

sf <- run(sf, data_on, paste0("results/outgroups/outgroups_", fit_type, "_LG_14_100kb_1bps.csv"), "Oreochromis_mossambicus", fit_type)
p_om_1 <- plotter(sf, data_on, fit_type) + ylim(-2, max(data_on$V4+5))
ggsave(paste0(figure_save_path, "Oreochromis_mossambicus_1bps.svg"), width=5, height=3, p_om_1)

sf <- run(sf, data_on, paste0("results/outgroups/outgroups_", fit_type, "_LG_14_100kb_2bps.csv"), "Oreochromis_mossambicus", fit_type)
p_om_2 <- plotter(sf, data_on, fit_type) + ylim(-2, max(data_on$V4+5))
ggsave(paste0(figure_save_path, "Oreochromis_mossambicus_2bps.svg"), width=5, height=3, p_om_2)

sf <- run(sf, data_on, paste0("results/outgroups/outgroups_", fit_type, "_LG_14_100kb_3bps.csv"), "Oreochromis_mossambicus", fit_type)
p_om_3 <- plotter(sf, data_on, fit_type) + ylim(-2, max(data_on$V4+5))
ggsave(paste0(figure_save_path, "Oreochromis_mossambicus_3bps.svg"), width=5, height=3, p_om_3)

sf <- run(sf, data_on, paste0("results/outgroups/outgroups_", fit_type, "_LG_14_100kb_4bps.csv"), "Oreochromis_mossambicus", fit_type)
p_om_4 <- plotter(sf, data_on, fit_type) + ylim(-2, max(data_on$V4+5))
ggsave(paste0(figure_save_path, "Oreochromis_mossambicus_4bps.svg"), width=5, height=3, p_om_3)

p_om <- p_om_1 / p_om_2 / p_om_3 / p_om_4
p_om <- p_om + plot_annotation(title="Oreochromis mossambicus")

p <- p_cg | p_om | p_on
# p <- grid.arrange(p_cg_1, p_cg_2, nrow=2)
ggsave(paste0(figure_save_path, "outgroups_", fit_type, "_raw.svg"), width=6.5, height=4, p)
