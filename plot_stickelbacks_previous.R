library(StrataFinder)
library(ggplot2)


data_file_path <- "data/sticklebacks/ENSEMBL_X_to_Y_v2.csv" # Stickelback

previous_bps <- data.frame("bps"=c(6890000, 12500000))
figure_save_file <- "figures/paper_figures/sticklebacks_previous.svg"

df <- read.table(data_file_path, sep=",", header=TRUE)
df <- df[df$ds < 0.5,]
df <- df[order(df$HiC_XPosition1),]
x <- as.integer(df$HiC_XPosition1)
y <- as.numeric(df$ds)

p <- ggplot() + geom_point(aes(x=x, y=y), size=0.5) +
                xlim(0, max(x)) + xlab("") + ylab("sex-specific snp density") +
                geom_vline(xintercept=previous_bps$bps, color="blue") +
                coord_cartesian(xlim=c(0,max(df$HiC_XPosition1))) +
                scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000), labels=c("5,000,000", "10,000,000", "15,000,000", "20,000,000")) 

p <- p + theme(text=element_text(family="TT Arial")) + ylab("")
ggsave(figure_save_file, width=7, height=3, p)
