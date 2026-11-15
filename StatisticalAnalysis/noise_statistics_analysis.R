#!/usr/bin/env Rscript
# Analyse Pure Gaussian Noise After Quantification Pipeline
#
# This script loads the results from noise_statistics_test.m, which generated
# 100 pure Gaussian white noise spectra and quantified them via AMARES and
# integration for both 13C and 2H. 

if(!require("tidyverse")) install.packages("tidyverse")
if(!require("fitdistrplus")) install.packages("fitdistrplus")
if(!require("R.matlab")) install.packages("R.matlab")
library(tidyverse)
library(fitdistrplus)
library(R.matlab)
library(ggplot2)
library(gridExtra)
library(VGAM)
library(EnvStats)

set.seed(123456)

# Load MATLAB results
cat("Loading noise statistics test results...\n")
mat_data <- readMat("noise_statistics_test_results.mat")

# Extract the four quantification method results (each is 100 samples x 5 peaks)
carbon_amares <- mat_data$carbon.amares.amplitudes
carbon_integration <- mat_data$carbon.integration.amplitudes
deuterium_amares <- mat_data$deuterium.amares.amplitudes
deuterium_integration <- mat_data$deuterium.integration.amplitudes

noise_std <- mat_data$noise.std
n_samples <- mat_data$n.samples

cat(sprintf("Loaded %d samples with noise std = %.2f\n", n_samples, noise_std))

# Create output directory for plots
output_dir <- "noise_distribution_analysis"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Function to fit multiple distributions to a data vector
fit_multiple_distributions <- function(data_vector, column_name) {
    # Initialize results list
    fit_results <- list()

    # Remove NaN values
    data_vector <- data_vector[!is.na(data_vector)]

    # Check if we have enough data
    if (length(data_vector) < 10) {
        cat("  Skipping", column_name, "- insufficient data points\n")
        return(NULL)
    }

    # Check for sufficient variation
    data_range <- max(data_vector) - min(data_vector)
    data_sd <- sd(data_vector)

    if (is.na(data_range) || !is.finite(data_range) || data_range < 1e-10 ||
        is.na(data_sd) || !is.finite(data_sd) || data_sd < 1e-10) {
        cat("  Skipping", column_name, "- insufficient variation",
            "(range =", data_range, ", sd =", data_sd, ")\n")
        return(NULL)
    }

    # Fit distributions with error handling
    tryCatch({
        fit_results$gamma <- fitdist(data_vector, "gamma")
    }, error = function(e) {
        cat("  Gamma fit failed for", column_name, "\n")
    })

    tryCatch({
        # Use analytic estimate for starting values
        start_scale <- sqrt(mean(data_vector^2) / 2)
        fit_results$rice <- fitdist(data_vector, "rice", start = list(sigma = start_scale, vee = 1))
    }, error = function(e) {
        cat("  Rice fit failed for", column_name, "\n")
    })

    tryCatch({
        fit_results$norm <- fitdist(data_vector, "norm")
    }, error = function(e) {
        cat("  Normal fit failed for", column_name, "\n")
    })

    tryCatch({
        fit_results$chisq <- fitdist(data_vector, "chisq", start = list(df = 3))
    }, error = function(e) {
        cat("  Chi-squared fit failed for", column_name, "\n")
    })

    tryCatch({
        fit_results$lnorm <- fitdist(data_vector, "lnorm")
    }, error = function(e) {
        cat("  Log-normal fit failed for", column_name, "\n")
    })

    tryCatch({
        fit_results$weibull <- fitdist(data_vector, "weibull")
    }, error = function(e) {
        cat("  Weibull fit failed for", column_name, "\n")
    })

    return(fit_results)
}

# Function to calculate AIC and BIC for all fitted distributions
calculate_metrics <- function(fit_list, column_name) {
    if (is.null(fit_list) || length(fit_list) == 0) {
        return(NULL)
    }

    metrics <- data.frame(
        Method = character(),
        Distribution = character(),
        AIC = numeric(),
        BIC = numeric(),
        stringsAsFactors = FALSE
    )

    for (dist_name in names(fit_list)) {
        aic_val <- AIC(fit_list[[dist_name]])
        bic_val <- BIC(fit_list[[dist_name]])

        metrics <- rbind(metrics, data.frame(
            Method = column_name,
            Distribution = dist_name,
            AIC = aic_val,
            BIC = bic_val,
            stringsAsFactors = FALSE
        ))
    }

    return(metrics)
}

# Function to create Cullen-Frey plot
create_cullen_frey_plot <- function(data_vector, title_text, output_file) {
    # Remove NaN values
    data_vector <- data_vector[!is.na(data_vector)]

    if (length(data_vector) < 10) {
        cat("  Insufficient data for Cullen-Frey plot:", title_text, "\n")
        return(NULL)
    }

    # Check for sufficient variation
    data_range <- max(data_vector) - min(data_vector)
    if (is.na(data_range) || !is.finite(data_range) || data_range < 1e-10) {
        cat("  Insufficient variation for Cullen-Frey plot:", title_text,
            " (range =", data_range, ")\n")
        return(NULL)
    }

    pdf(output_file, width = 8, height = 6)
    tryCatch({
        descdist(data_vector, boot = 1000)
        title(main = title_text, cex.main = 1.2)
    }, error = function(e) {
        cat("  Cullen-Frey plot failed for", title_text, ":", e$message, "\n")
    })
    dev.off()
    cat("  Saved Cullen-Frey plot to:", output_file, "\n")
}

# Function to create diagnostic plots for fitted distributions
create_diagnostic_plots <- function(fit_list, title_text, output_file) {
    if (is.null(fit_list) || length(fit_list) == 0) {
        return(NULL)
    }

    pdf(output_file, width = 12, height = 10)

    for (dist_name in names(fit_list)) {
        # plot.fitdist creates a 2x2 layout automatically, so don't override with par(mfrow)
        plot(fit_list[[dist_name]])
        # Add overall title using mtext
        mtext(paste(title_text, "-", dist_name), side = 3, line = -2, outer = TRUE, cex = 1.2, font = 2)
    }

    dev.off()
    cat("  Saved diagnostic plots to:", output_file, "\n")
}

# Main analysis function
analyze_quantification_method <- function(amplitude_matrix, method_name, metabolite_names) {
    cat("\n=== Analyzing", method_name, "===\n")

    all_metrics <- data.frame()
    all_fits <- list()

    # Analyze each metabolite (column) across all samples (rows)
    n_metabolites <- ncol(amplitude_matrix)

    for (i in 1:n_metabolites) {
        metabolite_name <- metabolite_names[i]
        full_name <- paste0(method_name, "_", metabolite_name)

        cat("Analyzing", full_name, "\n")

        # Extract data vector for this metabolite
        data_vector <- amplitude_matrix[, i]

        # Create Cullen-Frey plot
        cullen_file <- file.path(output_dir, paste0(full_name, "_cullen_frey.pdf"))
        create_cullen_frey_plot(data_vector, full_name, cullen_file)

        # Fit distributions
        fits <- fit_multiple_distributions(data_vector, full_name)

        if (!is.null(fits) && length(fits) > 0) {
            all_fits[[full_name]] <- fits

            # Calculate metrics
            metrics <- calculate_metrics(fits, full_name)
            if (!is.null(metrics)) {
                all_metrics <- rbind(all_metrics, metrics)
            }

            # Create diagnostic plots
            diagnostic_file <- file.path(output_dir, paste0(full_name, "_diagnostics.pdf"))
            create_diagnostic_plots(fits, full_name, diagnostic_file)
        }
    }

    return(list(
        metrics = all_metrics,
        fits = all_fits
    ))
}

# Define metabolite names
metabolite_names_13C <- c("Pyruvate", "Lactate", "Bicarbonate", "Alanine", "Pyruvate-Hydrate")
metabolite_names_2H <- c("Water", "Glucose", "Glx", "Lactate", "Other")

# Analyze all four quantification methods
cat("\n====================================\n")
cat("Analyzing Pure Noise Distributions\n")
cat("====================================\n")

results_carbon_amares <- analyze_quantification_method(
    carbon_amares,
    "13C_AMARES",
    metabolite_names_13C
)

results_carbon_integration <- analyze_quantification_method(
    carbon_integration,
    "13C_Integration",
    metabolite_names_13C
)

results_deuterium_amares <- analyze_quantification_method(
    deuterium_amares,
    "2H_AMARES",
    metabolite_names_2H
)

results_deuterium_integration <- analyze_quantification_method(
    deuterium_integration,
    "2H_Integration",
    metabolite_names_2H
)

# Combine all metrics (some may be NULL if data had insufficient variation)
all_metrics <- rbind(
    results_carbon_amares$metrics,
    results_carbon_integration$metrics,
    results_deuterium_amares$metrics,
    results_deuterium_integration$metrics
)

# Check if we have any valid results
if (is.null(all_metrics) || nrow(all_metrics) == 0) {
    cat("\n=== WARNING ===\n")
    cat("No valid metrics were generated!\n")
    cat("This likely means all quantification methods produced near-zero values.\n")
    cat("For pure noise with integration: this is expected (real spectrum averages to zero).\n")
    cat("Analysis cannot proceed without valid data.\n")
    quit(save = "no", status = 1)
}

# Save metrics to CSV
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
metrics_file <- file.path(output_dir, paste0("noise_distribution_metrics_", timestamp, ".csv"))
write.csv(all_metrics, metrics_file, row.names = FALSE)
cat("\nSaved metrics to:", metrics_file, "\n")

# Find best distribution for each method based on AIC
cat("\n=== Best Fitting Distributions (by AIC) ===\n")
best_distributions <- all_metrics %>%
    group_by(Method) %>%
    slice_min(order_by = AIC, n = 1) %>%
    dplyr::select(Method, Distribution, AIC, BIC) %>%
    arrange(Method)

print(best_distributions)

# Add a column to identify quantification method type (AMARES vs Integration)
best_distributions <- best_distributions %>%
    mutate(QuantMethod = ifelse(grepl("AMARES", Method), "AMARES", "Integration"))

# Count how many times each distribution is best
cat("\n=== Distribution Frequency as Best Fit ===\n")
best_dist_summary <- best_distributions %>%
    count(Distribution) %>%
    arrange(desc(n))

print(best_dist_summary)

# Count by quantification method type
cat("\n=== Distribution Frequency by Quantification Method ===\n")
best_by_method <- best_distributions %>%
    group_by(QuantMethod, Distribution) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(QuantMethod, desc(n))

print(best_by_method)

# Create summary plots
cat("\n=== Creating Summary Plots ===\n")

# Plot AIC values by distribution
p1 <- ggplot(all_metrics, aes(x = Distribution, y = AIC, color = Method)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(title = "AIC Values for Pure Noise Quantification",
         subtitle = "Lower AIC indicates better fit",
         x = "Distribution", y = "AIC")

ggsave(file.path(output_dir, "aic_comparison.pdf"), p1, width = 10, height = 6)
cat("Saved AIC comparison plot\n")

# Plot BIC values by distribution
p2 <- ggplot(all_metrics, aes(x = Distribution, y = BIC, color = Method)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(title = "BIC Values for Pure Noise Quantification",
         subtitle = "Lower BIC indicates better fit",
         x = "Distribution", y = "BIC")

ggsave(file.path(output_dir, "bic_comparison.pdf"), p2, width = 10, height = 6)
cat("Saved BIC comparison plot\n")

# Bar plot of best distributions (overall)
p3 <- ggplot(best_dist_summary, aes(x = reorder(Distribution, n), y = n, fill = Distribution)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Frequency of Best-Fitting Distributions",
         subtitle = "Across all metabolites and quantification methods",
         x = "Distribution", y = "Number of Cases")

ggsave(file.path(output_dir, "best_distribution_counts.pdf"), p3, width = 8, height = 6)
cat("Saved best distribution counts plot\n")

# Comparison plot: AMARES vs Integration
p4 <- ggplot(best_by_method, aes(x = Distribution, y = n, fill = QuantMethod)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Best-Fitting Distributions: AMARES vs Integration",
         subtitle = "Comparison of quantification methods",
         x = "Distribution", y = "Frequency",
         fill = "Quantification Method") +
    scale_fill_manual(values = c("AMARES" = "#E69F00", "Integration" = "#56B4E9"))

ggsave(file.path(output_dir, "amares_vs_integration_comparison.pdf"), p4, width = 10, height = 6)
cat("Saved AMARES vs Integration comparison plot\n")

# Faceted plot by quantification method
p5 <- ggplot(best_by_method, aes(x = reorder(Distribution, n), y = n, fill = Distribution)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ QuantMethod, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Best-Fitting Distributions by Quantification Method",
         subtitle = "Separate panels for AMARES and Integration",
         x = "Distribution", y = "Frequency")

ggsave(file.path(output_dir, "best_distributions_faceted.pdf"), p5, width = 10, height = 6)
cat("Saved faceted distributions plot\n")

# Stacked percentage plot
best_by_method_pct <- best_by_method %>%
    group_by(QuantMethod) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

p6 <- ggplot(best_by_method_pct, aes(x = QuantMethod, y = pct, fill = Distribution)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() +
    labs(title = "Distribution of Best Fits: AMARES vs Integration",
         subtitle = "Percentage breakdown by quantification method",
         x = "Quantification Method", y = "Percentage (%)",
         fill = "Distribution") +
    scale_y_continuous(labels = function(x) paste0(x, "%"))

ggsave(file.path(output_dir, "best_distributions_percentage.pdf"), p6, width = 8, height = 6)
cat("Saved percentage breakdown plot\n")

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", output_dir, "\n")
cat("\nKey Finding: What distribution does pure Gaussian noise follow\n")
cat("             after passing through AMARES/Integration quantification?\n")
cat("See Cullen-Frey graphs and best fit tables above.\n")

cat("\n=== Important Note on Integration Methods ===\n")
if (is.null(results_carbon_integration$metrics) || nrow(results_carbon_integration$metrics) == 0) {
    cat("13C Integration: Insufficient variation detected.\n")
    cat("  Real spectrum integration of pure noise should be close to zero but not exactly zero.\n")
    cat("  If exactly zero, reduce spectral points to increase variance.\n")
}
if (is.null(results_deuterium_integration$metrics) || nrow(results_deuterium_integration$metrics) == 0) {
    cat("2H Integration: Insufficient variation detected.\n")
    cat("  Real spectrum integration of pure noise should be close to zero but not exactly zero.\n")
    cat("  If exactly zero, reduce spectral points to increase variance.\n")
}
