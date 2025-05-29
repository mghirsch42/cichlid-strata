library(ggplot2)
library(ggpubr)
library(patchwork)

data_fit_type <- "exponential"
data_bps <- c(0, 2, 4)
data_noise <- c(0.01, 0.05)
est_fit_type <- c("linear", "exponential", "quadratic")
tests <- c("t", "f", "either", "both")

results_df <- data.frame()

for (eft in est_fit_type) {
    for (dbp in data_bps) {
        for (dn in data_noise) {
            for (test in tests) {
                results_file <- paste0("results/simulation/data_", data_fit_type, "/est_", eft, "/nbps_", dbp, "_noise_", dn, "_test_", test, ".csv")
                temp <- read.csv(results_file)
                temp$est_fit_type <- eft
                temp$true_n_bps <- dbp
                temp$true_n_bps <- cut(temp$true_n_bps, breaks=c(-Inf, data_bps), labels=data_bps)
                temp$noise <- dn
                temp$noise <- cut(temp$noise, breaks=c(-Inf, data_noise), labels=data_noise)
                temp$test <- test
                true_results_file <- paste0("data/simulation/", data_fit_type, "/nbps_", dbp, "/noise_", dn, "/bps.csv")
                temp2 <- read.csv(true_results_file)
                temp2 <- subset(temp2, select=c("sim_idx", "bps"))
                colnames(temp2)[colnames(temp2) == "bps"] <- "true_bps_idxs"
                temp3 <- merge(temp, temp2, by="sim_idx")
                results_df <- rbind(results_df, temp3)
            }
        }
    }
}
results_df[results_df$test=="t", "test"] <- "t-test"
results_df[results_df$test=="f", "test"] <- "f-test"
results_df$test <- factor(results_df$test, levels = c("t-test", "f-test", "either", "both"), ordered=TRUE)
# print(head(results_df))
# hist(results_df$n_bps)
# quit()
# bp_comp <- subset(results_df[results_df$n_bps == results_df$true_n_bps], select=c("sim_idx", "bps_idxs", "true_bps_idxs"))
bp_comp <- results_df[results_df$n_bps == results_df$true_n_bps,]
print(head(bp_comp))

compare_bps <- function(row) {
    true_bps <- lapply(strsplit(row["true_bps_idxs"][[1]], ","), as.integer)[[1]]
    est_bps <- lapply(strsplit(row["bps_idxs"][[1]], ","), as.integer)[[1]]
    a <- abs(true_bps-est_bps)
    return (sum(a))
}

# k <- apply(bp_comp, 1, compare_bps)
# print(paste("total runs", nrow(results_df)))
# print(paste("n correct n bps", length(k)))
# print(paste("n exact locations", sum(k==0)))
# print(paste("median", median(k)))
# print(paste("mean", mean(k)))
# print(paste("sd", sd(k)))
# print(paste("max", max(k)))

# hist(k)
# quit()
main_palette <- c("#FFFFFF", "#FFFFFF", "#FFFFFF")
noise_palette <- c("#E69F00", "#56B4E9")
# test_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
test_palette <- c("#00d073", "#F0E442", "#3880d0", "#e95e46", "#CC79A7")

df0 <- results_df[results_df$true_n_bps==0,]
df2 <- results_df[results_df$true_n_bps==2,]
df4 <- results_df[results_df$true_n_bps==4,]

p0_all <- ggplot(data=df0, aes(x=est_fit_type, y=n_bps, fill=est_fit_type)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=main_palette) + scale_y_continuous(breaks=seq(0, max(df0$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")
p0_noise <- ggplot(data=df0, aes(x=est_fit_type, y=n_bps, fill=noise)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=noise_palette) + scale_y_continuous(breaks=seq(0, max(df0$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")
p0_test <- ggplot(data=df0, aes(x=est_fit_type, y=n_bps, fill=test)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=test_palette) + scale_y_continuous(breaks=seq(0, max(df0$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")

p2_all <- ggplot(data=df2, aes(x=est_fit_type, y=n_bps, fill=est_fit_type)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=main_palette) + scale_y_continuous(breaks=seq(0, max(df2$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")
p2_noise <- ggplot(data=df2, aes(x=est_fit_type, y=n_bps, fill=noise)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=noise_palette) + scale_y_continuous(breaks=seq(0, max(df2$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")
p2_test <- ggplot(data=df2, aes(x=est_fit_type, y=n_bps, fill=test)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=test_palette) + scale_y_continuous(breaks=seq(0, max(df2$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")

p4_all <- ggplot(data=df4, aes(x=est_fit_type, y=n_bps, fill=est_fit_type)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=main_palette) + scale_y_continuous(breaks=seq(0, max(df4$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")
p4_noise <- ggplot(data=df4, aes(x=est_fit_type, y=n_bps, fill=noise)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=noise_palette) + scale_y_continuous(breaks=seq(0, max(df4$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")
p4_test <- ggplot(data=df4, aes(x=est_fit_type, y=n_bps, fill=test)) + geom_boxplot(outlier.size=.5) + scale_x_discrete(limits=c("linear", "exponential", "quadratic")) + scale_fill_manual(values=test_palette) + scale_y_continuous(breaks=seq(0, max(df4$n_bps), by=2)) + theme(text=element_text(family="TT Arial", size=8), legend.position="bottom") + ylab("estimated breakpoints") + xlab("estimation fit type")

plots <- list(p0_all, p0_noise, p0_test, p2_all, p2_noise, p2_test, p4_all, p4_noise, p4_test)
plt <- ggarrange(plotlist=plots, ncol=3, nrow=3)

ggsave(paste0("figures/paper_figures/simulation/", data_fit_type, "_raw.svg"), plt, width=6, height=5.75)

# quit()
# plots <- list()

# for (i in 1:3) {
#     eft <- est_fit_type[i]
#     curr_df <- results_df[results_df$est_fit_type == eft,]
#     # Plot just breakpoints
#     plt_all <- ggplot(data=curr_df, aes(x=true_n_bps, y=n_bps, group=true_n_bps)) +
#             geom_boxplot() +
#             ggtitle(paste(data_fit_type, "data,", eft, "estimation"))
#     plots[[3*(i-1)+1]] <- plt_all

#     # Plot by noise
#     plt_noise <- ggplot(data=curr_df, aes(x=true_n_bps, y=n_bps, fill=noise)) +
#             geom_boxplot() +
#             ggtitle(paste(data_fit_type, "data,", eft, "estimation"))
#     plots[[3*(i-1)+2]] <- plt_noise

#     # Plot by test
#     plt_test <- ggplot(data=curr_df, aes(x=true_n_bps, y=n_bps, fill=test)) +
#             geom_boxplot() +
#             ggtitle(paste(data_fit_type, "data,", eft, "estimation"))
#     plots[[3*(i-1)+3]] <- plt_test
    
# }

# plots <- as.list(plots)
# print(names(plots[2])) 
# quit()

# plt <- ggarrange(plotlist=plots, ncol=3, nrow=3)
# ggsave(paste0("figures/paper_figures/simulation/", data_fit_type, "_data.svg"), plt, width=10, height=7)

