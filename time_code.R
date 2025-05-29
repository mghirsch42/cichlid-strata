library(StrataFinder)
library(ggplot2)
library(purrr)

kb <- 10
data_file <- paste0("data/tropheini/", kb, "kb/LG5/Pseudosimochromis_marginatus.txt")
save_path <- "time_results/"
fit_type <- "exponential"
n_samples <- 100
noise <- 0.01
min_len <- 5

save_file <- paste0(save_path, kb, "kb_", fit_type, ".txt")
print(save_file)

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
            mx <- 4.42e-7 + rnorm(n=1, sd=1.6e-05)
            b <- 56 + rnorm(n=1, sd=272)
            mx2 <- 0
        }
        if (fit_type == "exponential") {
            mx <-2.79e-9 + rnorm(1, sd=1.02e-07)
            b <- 3.17 + rnorm(n=1, sd=2.56)
            mx2 <- 0
        }
        if (fit_type == "quadratic") {
            mx <- 0.10 + rnorm(n=1, sd=0.76)
            mx2 <- -2.44e-06 + rnorm(n=1, sd=1.31e-5)
            b <- -1379.73 + rnorm(n=1, sd=11189.2)
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


df <- read.table(data_file, sep="\t")
if (fit_type == "quadratic") {
    df$V2 <- df$V2 / 1000
}

times <- c()
n_bps <- 16 # number of breakpoints to simulate
for (s in 1:n_samples) {
    if (n_bps > 1) {
        flag <- TRUE
        while (flag) {
            flag <- FALSE
            bps <- sample(x=seq.int(min_len, nrow(df)-(min_len)), size=n_bps)
            bps <- bps[order(bps)]
            for (i in 1:(length(bps)-1)) {
                if (bps[i+1] < (bps[i]+10)) {
                    flag <- TRUE
                }
            }
        }
    } else {
        bps <- sample(x=seq.int(min_len, nrow(df)-(min_len)), size=n_bps)
        bps <- bps[order(bps)]
    }
    x <- c()
    y <- c()
    mxs <- c()
    mx2s <- c()
    bs <- c()
    for (i in 0:n_bps) {
        seg <- next_segment(i, bps, n_bps, df, fit_type, min_len, noise)
        x <- c(x, seg$x)
        y <- c(y, seg$y)
        mxs <- c(mxs, seg$mx)
        mx2s <- c(mx2s, seg$mx2)
        bs <- c(bs, seg$b)
    }
    max_y <- max(y)
    x <- c()
    y <- c()
    for (i in 0:n_bps) {
        if (i == 0 & i == n_bps) {
            xi <- df$V2[1:nrow(df)]
        } else if (i == n_bps) {
            xi <- df$V2[bps[n_bps]:nrow(df)]
        } else if (i == 0) {
            xi <- df$V2[1:bps[i+1]]
        } else {
            xi <- df$V2[(bps[i]+1):bps[i+1]]
        }
        if (fit_type == "quadratic") {
            yi <- generate_y(xi, fit_type, mx=mxs[i+1], mx2=mx2s[i+1], b=bs[i+1], noise=noise*max_y)
        } else {
        yi <- generate_y(xi, fit_type, mx=mxs[i+1], b=bs[i+1], noise=noise*max_y)

        }
        x <- c(x, xi)
        y <- c(y, yi)
    }
        
    sf <- new("StrataFinder")
    sf$fit_type <- fit_type
    sf$min_len <- min_len
    sf$x <- df$V2
    sf$y <- df$V4

    t1 <- Sys.time()
    sf$precomputeSSRMat()
    t2 <- Sys.time()
    curr_results <- sf$findNStrata(2, quiet=TRUE)
    t3 <- Sys.time()
    curr_results <- sf$findNStrata(4, quiet=TRUE)
    t4 <- Sys.time()
    curr_results <- sf$findNStrata(8, quiet=TRUE)
    t5 <- Sys.time()
    curr_results <- sf$findNStrata(16, quiet=TRUE)
    t6 <- Sys.time()
    curr_results <- sf$findNStrata(32, quiet=TRUE)
    t7 <- Sys.time()

    ssr_time <- difftime(t2, t1, units="secs")
    run2 <- difftime(t3, t2, units="secs")
    run4 <- difftime(t4, t3, units="secs")
    run8 <- difftime(t5, t4, units="secs")
    run16 <- difftime(t6, t5, units="secs")
    run32 <- difftime(t7, t6, units="secs")

    times <- rbind(times, data.frame("ssr_time"=ssr_time, "R2"=run2, "R4"=run4, "R8"=run8, "R16"=run16, "R32"=run32))

    print(paste(s, ssr_time, run2, run4, run8, run16, run32))
}

print(times)
write.table(times, save_file, sep=",")

