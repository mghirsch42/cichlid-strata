library(StrataFinder)
library(ggplot2)
library(patchwork)
library(scales)


data_file_path <- "data/fungus_gnats/Bcop.Holo2.XpX.pool.Illumina.indHetdens.100kbwins.csv" # Where my fungus is

previous_bps <- data.frame("bps"=c(4050001, 13050001, 13550001, 37350001, 38350001, 40000001, 40100001, 49500001, 49600001, 52150001, 52550001, 54050001, 55950001, 57250001, 58450001, 61950001))

figure_save_path <- "figures/paper_figures/fungus_gnats_previous.svg"

df_all <- read.table(data_file_path, sep=",", header=TRUE)
df_all <- df_all[!is.na(df_all$het_H2_XpX_pool_reseq.Bcop_v2),]
df_all$counts <- df_all$het_H2_XpX_pool_reseq.Bcop_v2

df <- df_all[df_all$scaffold == "X",] 
x <- as.integer(df$start)
y <- as.numeric(df$counts)
print(y)
p <- ggplot() + geom_point(aes(x=x, y=y), size=0.5) +
                xlim(0, max(x)) + xlab("") +
                geom_vline(xintercept=previous_bps$bps, color="red") +
                ylim(0, .023) +
                scale_y_continuous(breaks=c(0, .0050, .010, .015, .020), labels=c("0", "0.005", "0.010", "0.015", "0.020")) +
                scale_x_continuous(labels=comma, breaks=c(0, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7))
p <- p + theme(text=element_text(family="TT Arial")) + ylab("")
ggsave(paste(figure_save_path, ".svg", sep=""), width=7, height=2, p)
