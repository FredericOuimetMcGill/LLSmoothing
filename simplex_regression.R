###########################################################################################
## Local linear smoothing for regression surfaces on the simplex using Dirichlet kernels ##
###########################################################################################

## Written by Frederic Ouimet (August 2024)

require("LaplacesDemon") # for the Dirichlet distribution
require("ggplot2") # for plotting

library(plotly) # for 3d plots
library(htmlwidgets)
library(webshot)

################
# Set the path #
################

path <- file.path("C://Users//fred1//Dropbox//Ouimet_Genest_projects//Genest_Ouimet_2024_LL_smoothing//_simulations", fsep = .Platform$file.sep)
setwd(path)

##############
# Parameters #
##############

d <- 2 # dimension of simplex
MCsim <- 10 ^ 3 # number of uniforms sampled for integral MC estimates
BB <- seq(0.01, 1, by = 0.01) # bandwidths for LSCV graphs

MM <- list("LL","NW") # list of Dirichlet kernel methods
KK <- c(7, 10, 14, 20) # indices for the mesh
JJ <- 1:3 # target regression function indices
RR <- 1:200 # replication indices

###########################
## Mesh of design points ##
###########################

# Mesh of the 2-dim simplex

mesh <- function(k) {
  w <- (k - 1 / sqrt(2)) / (k - 1)
  res <- list()
  for (i in 1:k) {
    for (j in i:k) {
      res <- append(res, list((c(w * (i - 1) + 0.5, w * (k - j) + 0.5) / (k + 1))))
    }
  }
  return(res)
}

#################################
## Target regression functions ##
#################################

# for one design point x

m <- function(j, x) { # x is a d-dim vector on the simplex
  if (j == 1) {
    # Case when j = 1
    res <- x[1] * x[2]
  } else if (j == 2) {
    # Case when j = 2
    res <- log(1 + x[1] + x[2])
  } else if (j == 3) {
    # Case when j = 3
    res <- sin(x[1]) + cos(x[2])
  } else if (j == 4) {
    # Case when j = 4
    res <- sqrt(x[1]) + sqrt(x[2])
  } else {
    # Default case if j is not 1, 2, 3, or 4
    warning("Invalid value of j. Should be 1, 2, 3, or 4.")
    res <- NULL
  }
  return(res)
}

# for a mesh of design points x_1, ..., x_n

mm <- function(j, xx) { # xx is a list of d-dim vectors on the simplex
  res <- list()
  n <- length(xx)
  for (i in 1:n) {
    res <- append(res, list(m(j, xx[[i]])))
  }
  return(res)
}

################
## Estimators ##
################

hat_m <- function(xx, b, s, j, method, y = NULL) {
  # xx is a list of d-dim vectors on the simplex, s is a d-dim vector on the simplex
  n <- length(xx)
  d <- length(xx[[1]])
  
  if (is.null(y)) {
    y <- as.numeric(mm(j, xx)) # without random noise (this is not observed)
    y <- y + 0.1 * IQR(y) # with random noise (this is observed)
  }
  
  u <- s / b + rep(1, d)
  v <- (1 - sum(s)) / b + 1
  kernel_vec <- rep(NA, n)
  for (i in 1:n) {
    kernel_vec[i] <- LaplacesDemon::ddirichlet(c(xx[[i]], 1 - sum(xx[[i]])), c(u, v), log = FALSE)
  }
  
  if (method == "LL") {
    design_mat <- matrix(1, nrow = n, ncol = d + 1)
    for (i in 1:n) {
      design_mat[i, -1] <- xx[[i]] - s
    }
    W <- diag(kernel_vec)
    
    return(solve(t(design_mat) %*% W %*% design_mat, t(design_mat) %*% W %*% y)[1])
  } else if (method == "NW") {
    return(sum(y * kernel_vec) / sum(kernel_vec))
  } else {
    stop("Invalid method. Choose either 'LL' or 'NW'.")
  }
}

###################
## LSCV Function ##
###################

LSCV <- function(xx, b, j, method, unif_sample = NULL) {
  d <- length(xx[[1]])
  
  # Generate unif_sample if not provided
  if (is.null(unif_sample)) {
    unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
  }
  
  integ_squared <- rep(NA, MCsim)
  for (i in 1:MCsim) {
    integ_squared[i] <- (hat_m(xx, b, unif_sample[i, 1:d], j, method) - m(j, unif_sample[i, 1:d])) ^ 2
  }
  
  return(1 / factorial(d) * mean(integ_squared))
}

# Function to evaluate LSCV for a range of b values and plot the results
evaluate_LSCV <- function(xx, j, method, unif_sample = NULL) {
  b_values <- BB
  LSCV_values <- numeric(length(b_values))
  
  for (i in seq_along(b_values)) {
    b <- b_values[i]
    LSCV_values[i] <- LSCV(xx, b, j, method, unif_sample)
  }
  
  # Create a data frame for plotting
  results_df <- data.frame(b = b_values, LSCV = LSCV_values)
  
  # Plot the results
  ggplot(results_df, aes(x = b, y = LSCV)) +
    geom_line(color = "blue") +
    labs(title = "LSCV Values for different bandwidths b",
         x = "b",
         y = "LSCV") +
    theme_minimal()
}

################################
## Optimal Bandwidth Function ##
################################

b_opt <- function(xx, j, method, unif_sample = NULL) {
  # Generate unif_sample if not provided
  if (is.null(unif_sample)) {
    unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
  }
  
  objective_function <- function(b) {
    return(LSCV(xx, b, j, method, unif_sample))
  }
  res <- optimize(objective_function, interval = c(min(BB), max(BB)))
  return(res$minimum)
}

####################################
## Compare LSCV values with b_opt ##
####################################

# # Parameters
# xx <- mesh(10)
# j <- 1
# method <- "LL"
# unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
# 
# # Evaluate and plot LSCV
# evaluate_LSCV(xx, j, method, unif_sample)
# 
# # Compare with b_opt
# b_opt(xx, j, method, unif_sample)

###################
## MISE Function ##
###################

ISE <- function(xx, j, method) {
  d <- length(xx[[1]])
  
  unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
  hat_b <- b_opt(xx, j, method, unif_sample)
  
  return(LSCV(xx, hat_b, j, method, unif_sample))
}

###############
## Main code ##
###############

# First, load the parallel package
library(parallel)

# Initialize a cluster with 6 processors
cl <- makeCluster(6)

# Create data frames to store the results
summary_results <- data.frame(
  k = integer(),
  j = integer(),
  method = character(),
  mean_ISE = numeric(),
  sd_ISE = numeric(),
  median_ISE = numeric(),
  IQR_ISE = numeric(),
  stringsAsFactors = FALSE
)

raw_results <- data.frame(
  k = integer(),
  j = integer(),
  method = character(),
  ISE = numeric(),
  stringsAsFactors = FALSE
)

# Capture the start time
start_time <- Sys.time()

# Function to compute ISE values for a given combination of k, j, and method
compute_ISE <- function(j, k, method) {
  xx <- mesh(k) # Generate the mesh of design points
  ISE_values <- numeric(length(RR)) # Initialize a vector to store ISE values for each replication
  for (r in RR) {
    ISE_value <- ISE(xx, j, method)
    ISE_values[r] <- ISE_value
  }
  list(ISE_values = ISE_values, k = k, j = j, method = method)
}

# Export necessary objects and functions to the cluster
clusterExport(cl, list("KK", "JJ", "MM", "RR", "mesh", "ISE", "b_opt", "LSCV", "hat_m", "m", "mm", "LaplacesDemon", "compute_ISE", "d", "MCsim", "BB"))

# Loop over all combinations of j, k, and method in parallel
results_list <- parLapply(cl, JJ, function(j) {
  lapply(KK, function(k) {
    lapply(MM, function(method) {
      compute_ISE(j, k, method)
    })
  })
})

# Stop the cluster
stopCluster(cl)

# Process the results
for (j_results in results_list) {
  for (k_results in j_results) {
    for (result in k_results) {
      ISE_values <- result$ISE_values
      k <- result$k
      j <- result$j
      method <- result$method
      
      mean_ISE <- mean(ISE_values)
      sd_ISE <- sd(ISE_values)
      median_ISE <- median(ISE_values)
      IQR_ISE <- IQR(ISE_values)
      
      # Store the summary results
      summary_results <- rbind(
        summary_results,
        data.frame(
          k = k,
          j = j,
          method = method,
          mean_ISE = mean_ISE,
          sd_ISE = sd_ISE,
          median_ISE = median_ISE,
          IQR_ISE = IQR_ISE,
          stringsAsFactors = FALSE
        )
      )
      
      # Store the raw results
      for (ISE_value in ISE_values) {
        raw_results <- rbind(
          raw_results,
          data.frame(
            k = k,
            j = j,
            method = method,
            ISE = ISE_value,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

# Calculate the duration
print(Sys.time() - start_time)

# Save the summary results to a CSV file in the specified path
summary_output_file <- file.path(path, "ISE_results.csv")
write.csv(summary_results, summary_output_file, row.names = FALSE)

print("Summary results saved to ISE_results.csv")

# Save the raw results to a CSV file in the specified path
raw_output_file <- file.path(path, "raw_ISE_results.csv")
write.csv(raw_results, raw_output_file, row.names = FALSE)

print("Raw results saved to raw_ISE_results.csv")

########################
## Plotting boxplots ##
########################

# Load the raw results
raw_results <- read.csv(file.path(path, "raw_ISE_results.csv"))

# Convert 'method' and 'k' to factors for better plotting
raw_results$method <- as.factor(raw_results$method)
raw_results$k <- as.factor(raw_results$k)  # Convert k to factor

# Define the function to convert k to n
n <- function(k) {
  return(k * (k + 1) / 2)
}

# Create a custom labeller to replace k values with corresponding n values
labeller_k_to_n <- function(k_values) {
  sapply(k_values, function(k) paste0("n = ", n(as.numeric(as.character(k)))))  # Ensure k is numeric
}

# Create a custom labeller to replace j values with "j=1", "j=2", etc.
labeller_j <- function(j_values) {
  sapply(j_values, function(j) paste0("j = ", j))
}

# Combine the two labellers
custom_labeller <- labeller(
  k = labeller_k_to_n,
  j = labeller_j
)

# Save the plot as a high-quality PDF using cairo_pdf
plot_file <- file.path(path, "ISE_boxplots.pdf")
cairo_pdf(plot_file, width = 8, height = 5)  # Adjust width and height for aspect ratio

# Plot the boxplots side by side for each j and k with increased white space between rows
p <- ggplot(raw_results, aes(x = method, y = ISE, fill = method)) +
  geom_boxplot(width = 0.25) +  # Reduce the width of the boxes
  facet_grid(j ~ k, scales = "free", labeller = custom_labeller) +
  scale_y_log10() +  # Set y-axis to logarithmic scale
  labs(x = "Method",
       y = "ISE (log scale)") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.spacing.y = unit(2, "lines"), # Increase white space between rows
        strip.text.y = element_text(angle = 0)) # Ensure j labels are displayed horizontally

print(p)

# Close the PDF device
dev.off()

print("Boxplots saved to ISE_boxplots.pdf")


#############################################
## Illustration (Section 5)                ##
#############################################

vectors <- list(
  c(77.5, 19.5, 3.0, 10.4),
  c(71.9, 24.9, 3.2, 11.7),
  c(50.7, 36.1, 13.2, 12.8),
  c(52.2, 40.9, 6.6, 13.0),
  c(70.0, 26.5, 3.5, 15.7),
  c(66.5, 32.2, 1.3, 16.3),
  c(43.1, 55.3, 1.6, 18.0),
  c(53.4, 36.8, 9.8, 18.7),
  c(15.5, 54.4, 30.1, 20.7),
  c(31.7, 41.5, 26.8, 22.1),
  c(65.7, 27.8, 6.5, 22.4),
  c(70.4, 29.0, 0.6, 24.4),
  c(17.4, 53.6, 29.0, 25.8),
  c(10.6, 69.8, 19.6, 32.5),
  c(38.2, 43.1, 18.7, 33.6),
  c(10.8, 52.7, 36.5, 36.8),
  c(18.4, 50.7, 30.9, 37.8),
  c(4.6, 47.4, 48.0, 36.9),
  c(15.6, 50.4, 34.0, 42.2),
  c(31.9, 45.1, 23.0, 47.0),
  c(9.5, 53.5, 37.0, 47.1),
  c(17.1, 48.0, 34.9, 48.4),
  c(10.5, 55.4, 34.1, 49.4),
  c(4.8, 54.7, 41.0, 49.5),
  c(2.6, 45.2, 52.2, 59.2),
  c(11.4, 52.7, 35.9, 60.1),
  c(6.7, 46.9, 46.4, 61.7),
  c(6.9, 49.7, 43.4, 62.4),
  c(4.0, 44.9, 51.1, 69.3),
  c(7.4, 51.6, 40.9, 73.6),
  c(4.8, 49.5, 45.7, 74.4),
  c(4.5, 48.5, 47.0, 78.5),
  c(6.6, 52.1, 41.3, 82.9),
  c(6.7, 47.3, 45.9, 87.7),
  c(7.4, 45.6, 46.9, 88.1),
  c(6.0, 48.9, 45.1, 90.4),
  c(6.3, 53.8, 39.9, 90.6),
  c(2.5, 48.0, 49.5, 97.7),
  c(2.0, 47.8, 50.2, 103.7)
)

# Extract xx and y from the given vectors
xx <- lapply(vectors, function(v) v[1:2] / 100)
y <- as.numeric(lapply(vectors, function(v) v[4]))

# Generate a grid of points for s
s_grid <- expand.grid(seq(0, 1, length.out = 100), seq(0, 1, length.out = 100))
s_grid <- s_grid[rowSums(s_grid) <= 1, ]

# Define the function to calculate the RSS using LOOCV for a given bandwidth b
RSS_LOOCV <- function(xx, y, b, method) {
  n <- length(xx)
  rss <- 0
  for (i in 1:n) {
    # Exclude the i-th observation
    xx_minus_i <- xx[-i]
    y_minus_i <- y[-i]
    s <- xx[[i]] # Use double brackets to access the i-th vector directly
    y_i <- y[[i]] # Use double brackets to access the i-th response value directly
    y_hat_i <- hat_m(xx_minus_i, b, s, j, method, y_minus_i)
    rss <- rss + (y_i - y_hat_i) ^ 2
  }
  return(rss / n)
}

# Define the function to find the optimal bandwidth by minimizing the RSS using LOOCV
b_opt_LOOCV <- function(xx, y, method) {
  objective_function <- function(b) {
    RSS_LOOCV(xx, y, b, method)
  }
  res <- optimize(objective_function, interval = c(min(BB), max(BB)))
  return(res$minimum)
}

# Assuming xx and y are defined as in the original code
method <- "LL"
optimal_b <- b_opt_LOOCV(xx, y, method)

print(paste("Optimal bandwidth (b):", optimal_b))

# Plot RSS as a function of b
b_values <- seq(0.05, 1, by = 0.01)
rss_values <- sapply(b_values, function(b) RSS_LOOCV(xx, y, b, method))
results_df <- data.frame(b = b_values, RSS = rss_values)

# Save the plot to a PDF file in the specified path
plot_file <- file.path(path, "rss_b.pdf")
cairo_pdf(plot_file, width = 8, height = 5)  # Adjust width and height for aspect ratio

ggplot(results_df, aes(x = b, y = RSS)) +
  geom_line(color = "blue") +
  labs(x = "b", y = "RSS") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )

dev.off()

print(paste("RSS plot saved to", plot_file))
