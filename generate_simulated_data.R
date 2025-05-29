library(StrataFinder)
library(ggplot2)
library(purrr)
library(mcp)
library(optparse)
library(stringr)

# Generate command line flags
option_list <- list(
    make_option(c("-d", "--data_file"), action="store", type="character", default=""),
    make_option(c("-s", "--save_path"), action="store", type="character", default=""),
    make_option(c("-f", "--fit_type"), action="store", type="character", default="linear"),
    make_option(c("-b", "--n_bps"), action="store", type="integer", default=1),
    make_option(c("-m", "--min_len"), action="store", type="integer", default=5),
    make_option(c("-n", "--noise"), action="store", type="numeric", default=.5),
    make_option(c("-q", "--n_samples"), action="store", type="integer", default=10)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Error checking and implementing defaults.
if (opt$data_file == "") {
    print("No data file input. Tying to use default data file: data/tropheini/100kb/LG5/Pseudosimochromis_marginatus.txt")
    if (file.exists("data/tropheini/100kb/LG5/Pseudosimochromis_marginatus.txt")) {
        print("\tDefault file exists and will be used.")
    } else {
        print("\tDefault file does not exist. Please provide an appropriate data file.")
        quit()
    }
    opt$data_file <- "data/tropheini/100kb/LG5/Pseudosimochromis_marginatus.txt"
}
if (opt$save_path == "") {
    print("No save path input. Will not save results.")
    opt$save_path <- NULL
} else if (!file.exists(opt$save_path)) {
    print("Save path does not exist. Please create the directory.")
    quit()
}
if (opt$n_bps < 0) {
    print("Number of breakpoints cannot be negative.")
    quit()
}

writeLines(paste0("Generating simulated data with:\n",
                "\tData file for x points: ", opt$data_file, "\n",
                "\tSave path: ", opt$save_path, "\n",
                "\tFit type: ", opt$fit_type, "\n",
                "\tNumber of breakpoints: ", opt$n_bps, "\n",
                "\tMinimum length of segment between breakpoints: ", opt$min_len, "\n",
                "\tNoise (standard deviation of normal distribution): ", opt$noise, "\n",
                "\tNumber of samples: ", opt$n_samples))

# Functions for the different fit types
linear <- function(x, m, b, s) {
    data <- m*x + b + rnorm(n=1, sd=s)
    data[data<0] <- 0
    return(data) 
}
exponential <- function(x, m, b, s) {
    data <- exp(m*x+b) + rnorm(n=1, sd=s)
    data[data<0] <- 0
    return (data)
}
quadratic <- function(x, mx, mx2, b, s) {
    data <- b + mx*x + mx2*x*x + rnorm(n=1, sd=s)
    data[data<0] <- 0
    return (data)
}

next_segment <- function(i, bps, n_bps, df, fit_type, min_len, noise) {
    flag <- TRUE
    final_idx <- 0
    while (flag) {
        # Sample coefficients and intercept mean and standard deviation from tropheini estimates
        if (fit_type == "linear") {
            # 100kb
            # mx <- 4.42e-7 + rnorm(n=1, sd=1.6e-05)
            # b <- 56 + rnorm(n=1, sd=272)
            # mx2 <- 0

            # 10 kb
            mx <- 8.24e-7 + rnorm(n=1, sd=4.45e-6)
            b <- -16.15 + rnorm(n=1, sd=115)
            mx2 <- 0
        }
        if (fit_type == "exponential") {
            # 100 kb
            # mx <- 2.79e-9 + rnorm(1, sd=1.02e-07)
            # b <- 3.17 + rnorm(n=1, sd=2.56)
            # mx2 <- 0

            # 10 kb
            mx <- 6.52e-7 + rnorm(1, sd=4.58e-6)
            b <- 1.17 + rnorm(n=1, sd=3.85)
            mx2 <- 0
        }
        if (fit_type == "quadratic") {
            # 100 kb
            # mx <- 0.10 + rnorm(n=1, sd=0.76)
            # mx2 <- -2.44e-06 + rnorm(n=1, sd=1.31e-5)
            # b <- -1379.73 + rnorm(n=1, sd=11189.2)
            
            # 10 kb
            mx <- 0.021 + rnorm(n=1, sd=0.10)
            mx2 <- -3.83e-7 + rnorm(n=1, sd=1.76e-6)
            b <- -312 + rnorm(n=1, sd=1539)
        }
        if (i == 0 & i == n_bps) {
            xi <- df$V2[1:nrow(df)]
        }
        else if (i == n_bps) {
            xi <- df$V2[bps[n_bps]:nrow(df)]

        } else if (i == 0) {
            xi <- df$V2[1:bps[i+1]]

        } else {
            xi <- df$V2[(bps[i]+1):bps[i+1]]
        }
        if (fit_type == "quadratic") {
            yi <- generate_y(xi, fit_type, mx=mx, mx2=mx2, b=b, noise=0)
        
        } else {
            yi <- generate_y(xi, fit_type, mx=mx, b=b, noise=0)
        }
        # If we have sufficient non zero data points and appropriate maximum value, exit the loop
        if ((max(yi) < 1000) & (sum(yi==0) < 0.01*length(yi))) {
            flag = FALSE
        } 
    }
    return (list("x"=xi, "y"=yi, "mx"=mx, mx2=mx2, "b"=b))
}

# Function to generate y data based on the fit type
generate_y <- function(x, fit_type, mx=0, mx2=0, b=0, noise=0) {
    if (fit_type == "linear") {
        y <- unlist(map(x, linear, mx, b, noise))
    } else if (fit_type == "exponential") {
        y <- unlist(map(x, exponential, mx, b, noise))
    } else if (fit_type == "quadratic") {
        y <- unlist(map(x, quadratic, mx, mx2, b, noise))
    }
    return (y)
}


df <- read.table(opt$data_file, sep="\t")
if (opt$fit_type == "quadratic") {
    df$V2 <- df$V2 / 1000
}
results_bps <- data.frame()

for (s in 1:opt$n_samples) {
    if (opt$n_bps > 1) {
        flag <- TRUE
        while (flag) {
            flag <- FALSE
            bps <- sample(x=seq.int(opt$min_len, nrow(df)-(opt$min_len)), size=opt$n_bps)
            bps <- bps[order(bps)]
            for (i in 1:(length(bps)-1)) {
                if (bps[i+1] < (bps[i]+10)) {
                    flag <- TRUE
                }
            }
        }
    } else {
        bps <- sample(x=seq.int(opt$min_len, nrow(df)-(opt$min_len)), size=opt$n_bps)
        bps <- bps[order(bps)]
    }
    x <- c()
    y <- c()
    mxs <- c()
    mx2s <- c()
    bs <- c()
    for (i in 0:opt$n_bps) {
        seg <- next_segment(i, bps, opt$n_bps, df, opt$fit_type, opt$min_len, opt$noise)
        x <- c(x, seg$x)
        y <- c(y, seg$y)
        mxs <- c(mxs, seg$mx)
        mx2s <- c(mx2s, seg$mx2)
        bs <- c(bs, seg$b)
    }
    max_y <- max(y)
    x <- c()
    y <- c()
    for (i in 0:opt$n_bps) {
        if (i == 0 & i == opt$n_bps) {
            xi <- df$V2[1:nrow(df)]
        } else if (i == opt$n_bps) {
            xi <- df$V2[bps[opt$n_bps]:nrow(df)]
        } else if (i == 0) {
            xi <- df$V2[1:bps[i+1]]
        } else {
            xi <- df$V2[(bps[i]+1):bps[i+1]]
        }
        if (opt$fit_type == "quadratic") {
            yi <- generate_y(xi, opt$fit_type, mx=mxs[i+1], mx2=mx2s[i+1], b=bs[i+1], noise=opt$noise*max_y)
        } else {
        yi <- generate_y(xi, opt$fit_type, mx=mxs[i+1], b=bs[i+1], noise=opt$noise*max_y)

        }
        x <- c(x, xi)
        y <- c(y, yi)
    }
    if (opt$fit_type == "quadratic") {
        x <- x*1000
    }
    df2 <- data.frame("x"=x, "y"=y)
    curr_df <- data.frame(
                "sim_idx"=s,
                "bps"=paste(bps, collapse=","), 
                "intercepts"=paste(bs, collapse=","), 
                "coefs"=paste(mxs, collapse=","))
    if (opt$fit_type == "quadratic") {
        curr_df$coefs2 <- paste(mx2s, collapse=",")
    }
    results_bps <- rbind(results_bps, curr_df)
    if (!is.null(opt$save_path)) {
        write.csv(df2, paste0(opt$save_path, s, ".csv"), row.names=FALSE)
    }
}
print(results_bps)
if (!is.null(opt$save_path)) {
    write.csv(results_bps, paste0(opt$save_path, "bps.csv"), row.names=FALSE)
}