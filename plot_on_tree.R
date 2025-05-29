library(ggtree)
library(ggplot2)
library(tidytree)
library(ggimage)

lg <- 19
fit_type <- "linear"
kb <- 100
n_bps <- 3

# Load the data
tree_file <- "data/tree2.txt"
dotplt_path <- paste0("figures/tropheini", "/LG", lg, "/", fit_type, "/", kb, "kb/")
fig_save_path <- paste0("figures/paper_figures/tropheini_tree_lg", lg, "_fit/", fit_type, ".svg")
tree <- read.tree(tree_file)

files <- list.files(dotplt_path)
files <- files[grepl(".svg", files, fixed = TRUE)]
species <- gsub(".svg", "", files)
species <- species[species %in% c("Tropheus_sp_black", "Tropheus_brichardi", "Tropheus_polli") == FALSE]

tree_df <- as_tibble(tree)
to_drop <- tree_df$label[!(tree_df$label %in% species) & (!is.na(tree_df$label))]
tree <- drop.tip(tree, to_drop)

p <- ggtree(tree) + geom_tiplab(size=3) +
        geom_tiplab(aes(image=paste0(dotplt_path, label, ".svg")), geom="image", offset=30, size=0.18) +
        xlim(0, 60)
p <- rotate(p, rootnode(tree))
ggsave(fig_save_path, width=6.5, height=7.75, p)