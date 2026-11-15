if(!require("tidyverse")) install.packages("tidyverse")
if(!require("fitdistrplus")) install.packages("fitdistrplus")
if(!require("tidyr")) install.packages("tidyr")
if(!require("gridExtra")) install.packages("gridExtra")
if(!require("purrr")) install.packages("purrr")
if(!require("VGAM")) install.packages("VGAM")
if(!require("EnvStats")) install.packages("EnvStats")



library(tidyverse)
library(fitdistrplus)
set.seed(123456) #Always the answer....
#setwd("~/you/should/fill/this/out")


#This list of distributions is ONLY used for plotting one specific sub-figure (the count of best by AICs) 
# It is entirely cosmetic -- for considering additional distributions add them to the list of results below 
# (and, as is conventional for R, you need to define ddist, rdist, pdist, etc) 
BlankDists <- data.frame(Column = "None",Distribution=c("chisq","gamma","lnorm","normal","rice","weibull"), n=0)

# Load the data
#OXSA AMARES data 
exampleNoiseData <- read_csv("ExampleAmps_for_upload_oxsa.csv")
isCensored <- F; 

#jMRUI censored amps 
#exampleNoiseData <- read_csv("ExampleAmps-jMRUI.csv")
#isCensored <- T; 

#Start with some illustrative plots about _everything_. 
# Define noise samples range -- this can easily be increased but it takes a lot longer and doesn't change much...
samplesToUse <- (exampleNoiseData %>% dplyr::select(-SpectrumIndex) %>% slice_tail(n=20) %>% pull(Sample2) %>% unlist())

#Plots -- a few quick ones 
#Uncomment to save 
#pdf("NoiseStatsBicarbonate_2_zeros_removed.pdf")
plotdist(samplesToUse, histo = TRUE, demp = TRUE)
descdist(samplesToUse, boot=1000)
#dev.off()
if(isCensored){
    #Fit ignoring censoring 
    fittedDistribution <- fitdist(samplesToUse, "gamma", optim.method="Nelder-Mead")
    print(fittedDistribution)
    #Surve2fitdistcens requires a df, containing two columns, left = NA if left censored, 
    # right = NA if right censored 
    # The same number = left = right = observation if not censored. 
fitMe <- data.frame(left=samplesToUse,right=samplesToUse); 
fitMe <- fitMe %>% mutate(right = case_when(right == 0~ NA, .default = right))
fitted <- fitdistcens(fitMe, "gamma")
print(fitted)
} else { 
    #Both empirically chosen given the below
    fittedDistribution <- fitdist(samplesToUse, "lnorm", optim.method="Nelder-Mead")
}

# Required libraries
library(dplyr)
library(tidyr)
library(fitdistrplus)
library(ggplot2)
library(gridExtra)
library(purrr)
library(VGAM)
library(EnvStats)

# We will now fit _everything_. 
#################################################################
# Function to fit multiple distributions to a data vector
fit_multiple_distributions <- function(data_vector, column_name, is_data_censored=F) {
    # Initialize results list
    fit_results <- list()

    # Fit distributions without error trapping
    # To fit more distributions, extend this list
    if(is_data_censored) {
        fitMe <- data.frame(left=data_vector,right=data_vector); 
        fitMe <- fitMe %>% mutate(right = case_when(right == 0 ~ NA, .default = right), 
                                  left = case_when(left == 0 ~ 1e-304, .default = left)) #lnorm bails if this is 0 exactly 

        start_scale <- sqrt(mean(fitMe$right^2, na.rm = TRUE) / 2) #Analytic estimate for rice and rayleigh

        fit_results$gamma    <- fitdistcens(fitMe, "gamma")
        fit_results$rice     <- fitdistcens(fitMe, "rice", start = list(sigma = start_scale, vee = 1))
        fit_results$norm     <- fitdistcens(fitMe, "norm")
        fit_results$chisq    <- fitdistcens(fitMe, "chisq", start = list(df = mean(fitMe$right, na.rm = TRUE), ncp = 0.1), lower = c(0, 0)) #bounded; df = 3, ncp = 0 fails to start with fitdistrplus 1.2.6
        fit_results$lnorm    <- fitdistcens(fitMe, "lnorm")
        fit_results$weibull  <- fitdistcens(fitMe, "weibull")
        #fit_results$rayleigh <- fitdistcens(fitMe, "rayleigh", start=list(scale=start_scale) )
        #This is horribly numerically unstable for the censored distribution but works for the uncensored one 
        #Analytically both are weibull 

        #Store the raw, zero-excluded set of data here 
        fit_results$gamma$data     <- data_vector[!data_vector==0]  
        fit_results$rice$data      <- data_vector[!data_vector==0] 
        fit_results$norm$data      <- data_vector[!data_vector==0]   
        fit_results$chisq$data     <- data_vector[!data_vector==0] 
        fit_results$lnorm$data     <- data_vector[!data_vector==0] 
        fit_results$weibull$data   <- data_vector[!data_vector==0] 
        #fit_results$rayleigh$data  <- data_vector[!data_vector==0] 


    } else {
        fit_results$gamma    <- fitdist(data_vector, "gamma")
        fit_results$rice     <- fitdist(data_vector, "rice", start = list(sigma = 50, vee = 1))
        fit_results$norm     <- fitdist(data_vector, "norm")
        
        fit_results$chisq    <- fitdist(data_vector, "chisq", start = list(df = 3))
        fit_results$lnorm    <- fitdist(data_vector, "lnorm")
        fit_results$weibull  <- fitdist(data_vector, "weibull")
        #fit_results$rayleigh <- fitdist(data_vector, "rayleigh", start=list(scale=1), method="mse")

    }
    return(fit_results)
}

# Function to calculate AIC and BIC for all fitted distributions
calculate_metrics <- function(fit_list, column_name) {
    metrics <- data.frame(
                          Column = character(),
                          Distribution = character(),
                          AIC = numeric(),
                          BIC = numeric(),
                          stringsAsFactors = FALSE
    )

    for (dist_name in names(fit_list)) {
        aic_val <- AIC(fit_list[[dist_name]])
        bic_val <- BIC(fit_list[[dist_name]])

        metrics <- rbind(metrics, data.frame(
                                             Column = column_name,
                                             Distribution = dist_name,
                                             AIC = aic_val,
                                             BIC = bic_val,
                                             stringsAsFactors = FALSE
                                             ))
    }

    return(metrics)
}

# Function to create and save diagnostic plots
create_diagnostic_plots <- function(fit_list, column_name, base_dir = "plots") {
    # Create timestamp for filenames
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    # Create directory if it doesn't exist
    if (!dir.exists(base_dir)) {
        dir.create(base_dir, recursive = TRUE)
    }

    # First create Cullen-Frey plot
    pdf(file.path(base_dir, paste0(column_name, "_cullen_frey_", timestamp, ".pdf")))
    descdist(fit_list[[1]]$data, boot = 100)
    title(paste("Cullen-Frey Plot for", column_name))
    dev.off()

    # Create individual distribution fit plots
    for (dist_name in names(fit_list)) {
        pdf_filename <- file.path(base_dir, paste0(column_name, "_", dist_name, "_", timestamp, ".pdf"))
        pdf(pdf_filename)
        plot(fit_list[[dist_name]])
        title(paste("Distribution:", dist_name, "for", column_name))
        dev.off()
        cat("Saved plot to:", pdf_filename, "\n")
    }
}

# Main function to analyze all columns
analyze_all_columns <- function(df, n_tail = 50, isCensored) {
    # Get all column names except one called SpectrumIndex
    column_names <- setdiff(names(df), "SpectrumIndex")

    # Initialize results
    all_metrics <- data.frame()
    all_fits <- list()

    # Create a timestamp directory for this run
    timestamp_dir <- format(Sys.time(), "%Y%m%d_%H%M%S")
    plots_dir <- file.path("plots", timestamp_dir)

    # Analyze each column
    for (col in column_names) {
        cat("Analyzing column:", col, "\n")

        # Extract data vector (last n_tail points)
        data_vector <- df %>%
            slice_tail(n = n_tail) %>%
            dplyr::pull(col) %>%
            unlist()

        # Skip if all values are NA or if there's no variation
        if (all(is.na(data_vector)) || sd(data_vector, na.rm = TRUE) < 1e-10) {
            cat("Skipping column", col, "due to insufficient data\n")
            next
        }

        # Fit distributions
        fits <- fit_multiple_distributions(data_vector, col, isCensored)

        all_fits[[col]] <- fits

        # Calculate metrics
        metrics <- calculate_metrics(fits, col)
        all_metrics <- rbind(all_metrics, metrics)

        # Create and save diagnostic plots
        create_diagnostic_plots(fits, col, plots_dir)
    }

    return(list(
                metrics = all_metrics,
                fits = all_fits
                ))
}

analyze_all_columns <- function(df, n_tail = 50, isCensored) {
    # Get all column names except one called SpectrumIndex
    column_names <- setdiff(names(df), "SpectrumIndex")

    # Initialize results
    all_metrics <- data.frame()
    all_fits <- list()

    # Create a timestamp directory for this run
    timestamp_dir <- format(Sys.time(), "%Y%m%d_%H%M%S")
    plots_dir <- file.path("plots", timestamp_dir)

    # Analyze each column
    for (col in column_names) {
        cat("Analyzing column:", col, "\n")

        # Extract data vector (last n_tail points)
        data_vector <- df %>%
            slice_tail(n = n_tail) %>%
            dplyr::pull(col) %>%
            unlist()

        # Skip if all values are NA or if there's no variation
        if (all(is.na(data_vector)) || sd(data_vector, na.rm = TRUE) < 1e-10) {
            cat("Skipping column", col, "due to insufficient data\n")
            next
        }

        # Fit distributions
        fits <- fit_multiple_distributions(data_vector, col, isCensored)
        all_fits[[col]] <- fits

        # Calculate metrics
        metrics <- calculate_metrics(fits, col)
        all_metrics <- rbind(all_metrics, metrics)

        # Create and save diagnostic plots
        create_diagnostic_plots(fits, col, plots_dir)
    }

    return(list(
                metrics = all_metrics,
                fits = all_fits
                ))
}

# Function to create summary plots of AIC and BIC values
create_summary_plots <- function(metrics_df, base_dir = "plots") {
    # Create timestamp for filenames
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    # Create directory if it doesn't exist
    if (!dir.exists(base_dir)) {
        dir.create(base_dir, recursive = TRUE)
    }

    # Reshape data for plotting
    metrics_long <- metrics_df %>%
        pivot_longer(cols = c(AIC, BIC), names_to = "Metric", values_to = "Value")

    # Create plot with distributions on x-axis
    p1 <- ggplot(metrics_long, aes(x = Distribution, y = Value, color = Column)) +
        geom_point(position = position_dodge(width = 0.5), size = 3) +
        facet_wrap(~ Metric, scales = "free_y") +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.text = element_text(size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "Comparison of AIC and BIC Values by Distribution",
     x = "Distribution", y = "Value")

# Save the plot
pdf_file <- file.path(base_dir, paste0("comparison_plot_", timestamp, ".pdf"))
ggsave(pdf_file, p1, width = 6, height = 4)
cat("Saved comparison plot to:", pdf_file, "\n")

# Create boxplot of values by distribution
p2 <- ggplot(metrics_long, aes(x = Distribution, y = Value, fill = Distribution)) +
    geom_boxplot() +
    facet_wrap(~ Metric, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "Distribution of AIC and BIC Values",
     x = "Distribution", y = "Value")

# Save the plot
pdf_file <- file.path(base_dir, paste0("boxplot_", timestamp, ".pdf"))
ggsave(pdf_file, p2, width = 6, height = 4)
cat("Saved boxplot to:", pdf_file, "\n")

# Create histogram comparing distributions
p3 <- ggplot(metrics_long, aes(x = Value, fill = Distribution)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ Metric, scales = "free_x") +
    theme_bw() +
    labs(title = "Density Distribution of AIC and BIC Values",
         x = "Value", y = "Density")

    # Save the plot
    pdf_file <- file.path(base_dir, paste0("density_plot_", timestamp, ".pdf"))
    ggsave(pdf_file, p3, width = 6, height = 4)
    cat("Saved density plot to:", pdf_file, "\n")

    # Find best distribution for each column - using explicit dplyr::select to fix the error
    best_fits <- metrics_df %>%
        group_by(Column) %>%
        slice_min(order_by = AIC, n = 1) %>%
        dplyr::select(Column, Distribution, AIC, BIC)

    best_fits_including_blanks <- rbind(best_fits, BlankDists) #defined above 

    # Count how many times each distribution is the best fit
    best_dist_counts <- best_fits_including_blanks %>% ungroup() %>% 
        count(Distribution) %>% mutate(n=n-1)
    #arrange(desc(n))

    # Create bar plot of best distribution counts
    p4 <- ggplot(best_dist_counts, aes(x = Distribution, y = n, fill = Distribution)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_bw() +
        theme(legend.position = "none") +
        labs(title = "Number of Columns Best Fitted by Each Distribution",
             x = "Distribution", y = "Count")

        # Save the plot
        pdf_file <- file.path(base_dir, paste0("best_dist_plot_", timestamp, ".pdf"))
        ggsave(pdf_file, p4, width = 6, height = 4)
        cat("Saved best distribution plot to:", pdf_file, "\n")

        #Print a chi-squared test 
        chisq.test(table(best_dist_counts %>% ungroup() ))

        return(list(comparison_plot = p1,
                    boxplot = p2,
                    density_plot = p3,
                    best_dist_plot = p4))
}

#Function to extract parameters from fitted distributions
extract_parameters <- function(fit_obj, dist_name) {
    params <- list()

    if (dist_name == "gamma") {
        params$shape <- fit_obj$estimate["shape"]
        params$rate <- fit_obj$estimate["rate"]
    } else if (dist_name == "rice") {
        params$sigma <- fit_obj$estimate["sigma"]
        params$vee <- fit_obj$estimate["vee"]
    } else if (dist_name == "norm") {
        params$mean <- fit_obj$estimate["mean"]
        params$sd <- fit_obj$estimate["sd"]
    } else if (dist_name == "chisq") {
        params$df <- fit_obj$estimate["df"]
        params$ncp <- fit_obj$estimate["ncp"]
    } else if (dist_name == "lnorm") {
        params$meanlog <- fit_obj$estimate["meanlog"]
        params$sdlog <- fit_obj$estimate["sdlog"]
    } else if (dist_name == "weibull") {
        params$shape <- fit_obj$estimate["shape"]
        params$scale <- fit_obj$estimate["scale"]
    } else if (dist_name == "rayleigh") {
        params$scale <- fit_obj$estimate["scale"]
    }

    return(params)
}

# Enhanced function to calculate metrics AND parameters
calculate_metrics_with_params <- function(fit_list, column_name) {
    metrics <- data.frame(
                          Column = character(),
                          Distribution = character(),
                          AIC = numeric(),
                          BIC = numeric(),
                          stringsAsFactors = FALSE
    )

    # Create a list to store parameter columns dynamically
    param_columns <- list()

    for (dist_name in names(fit_list)) {
        # Get AIC and BIC
        aic_val <- AIC(fit_list[[dist_name]])
        bic_val <- BIC(fit_list[[dist_name]])

        # Extract parameters
        params <- extract_parameters(fit_list[[dist_name]], dist_name)

        # Create a row with basic metrics
        row_data <- data.frame(
                               Column = column_name,
                               Distribution = dist_name,
                               AIC = aic_val,
                               BIC = bic_val,
                               stringsAsFactors = FALSE
        )

        # Add parameter columns
        for (param_name in names(params)) {
            col_name <- paste0(dist_name, "_", param_name)
            row_data[[col_name]] <- params[[param_name]]
        }

        # Store parameter column names
        param_columns <- c(param_columns, names(params))

        # Add to metrics dataframe
        if (nrow(metrics) == 0) {
            metrics <- row_data
        } else {
            # Ensure all columns exist in both dataframes
            missing_in_metrics <- setdiff(names(row_data), names(metrics))
            missing_in_row <- setdiff(names(metrics), names(row_data))

            for (col in missing_in_metrics) {
                metrics[[col]] <- NA
            }
            for (col in missing_in_row) {
                row_data[[col]] <- NA
            }

            metrics <- rbind(metrics, row_data[names(metrics)])
        }
    }

    return(metrics)
}

# Modified main analysis function
analyze_all_columns_with_params <- function(df, n_tail = 50, isCensored) {
    # Get all column names except one called SpectrumIndex
    column_names <- setdiff(names(df), "SpectrumIndex")

    # Initialize results
    all_metrics <- data.frame()
    all_fits <- list()

    # Create a timestamp directory for this run
    timestamp_dir <- format(Sys.time(), "%Y%m%d_%H%M%S")
    plots_dir <- file.path("plots", timestamp_dir)

    # Analyze each column
    for (col in column_names) {
        cat("Analyzing column:", col, "\n")

        # Extract data vector (last n_tail points)
        data_vector <- df %>%
            slice_tail(n = n_tail) %>%
            dplyr::pull(col) %>%
            unlist()

        # Skip if all values are NA or if there's no variation
        if (all(is.na(data_vector)) || sd(data_vector, na.rm = TRUE) < 1e-10) {
            cat("Skipping column", col, "due to insufficient data\n")
            next
        }

        # Fit distributions
        fits <- fit_multiple_distributions(data_vector, col, isCensored)
        all_fits[[col]] <- fits

        # Calculate metrics with parameters
        metrics <- calculate_metrics_with_params(fits, col)

        # Add to all metrics
        if (nrow(all_metrics) == 0) {
            all_metrics <- metrics
        } else {
            # Ensure consistent columns
            missing_in_all <- setdiff(names(metrics), names(all_metrics))
            missing_in_metrics <- setdiff(names(all_metrics), names(metrics))

            for (col_name in missing_in_all) {
                all_metrics[[col_name]] <- NA
            }
            for (col_name in missing_in_metrics) {
                metrics[[col_name]] <- NA
            }

            all_metrics <- rbind(all_metrics, metrics[names(all_metrics)])
        }

        # Create and save diagnostic plots
        create_diagnostic_plots(fits, col, plots_dir)
    }

    return(list(
                metrics = all_metrics,
                fits = all_fits
                ))
}

# Function to create a summary table of parameters
create_parameter_summary <- function(metrics_df) {
    # Find best distribution for each column based on AIC
    best_fits <- metrics_df %>%
        group_by(Column) %>%
        slice_min(order_by = AIC, n = 1) %>%
        ungroup()

    # Create a summary focusing on the best fit parameters
    param_summary <- data.frame(
                                Column = character(),
                                Best_Distribution = character(),
                                AIC = numeric(),
                                BIC = numeric(),
                                stringsAsFactors = FALSE
    )

    for (i in 1:nrow(best_fits)) {
        row <- best_fits[i, ]
        dist_name <- row$Distribution

        # Extract relevant parameter columns for this distribution
        param_cols <- grep(paste0("^", dist_name, "_"), names(row), value = TRUE)

        summary_row <- data.frame(
                                  Column = row$Column,
                                  Best_Distribution = dist_name,
                                  AIC = row$AIC,
                                  BIC = row$BIC,
                                  stringsAsFactors = FALSE
        )

        # Add parameter values with cleaner names
        for (param_col in param_cols) {
            clean_name <- sub(paste0(dist_name, "_"), "", param_col)
            summary_row[[clean_name]] <- row[[param_col]]
        }

        if (nrow(param_summary) == 0) {
            param_summary <- summary_row
        } else {
            # Handle missing columns
            missing_in_summary <- setdiff(names(summary_row), names(param_summary))
            missing_in_row <- setdiff(names(param_summary), names(summary_row))

            for (col in missing_in_summary) {
                param_summary[[col]] <- NA
            }
            for (col in missing_in_row) {
                summary_row[[col]] <- NA
            }

            param_summary <- rbind(param_summary, summary_row[names(param_summary)])
        }
    }

    return(param_summary)
}

#######################################
#######################################
#######################################
# Execute the analysis

results <- analyze_all_columns(exampleNoiseData, n_tail=50, isCensored)

# Create and save summary plots
# There is a bug in the creation of the AIC/BIC for the rayleigh distribution with censored data 
summary_plots <- create_summary_plots(results$metrics)

# Save results to CSV
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
write.csv(results$metrics, paste0("distribution_fit_metrics_", timestamp, ".csv"), row.names = FALSE)

# Print summary table of best distribution for each column - using explicit dplyr::select
best_distributions <- results$metrics %>%
    group_by(Column) %>%
    slice_min(order_by = AIC, n = 1) %>%
    dplyr::select(Column, Distribution, AIC, BIC) %>%
    arrange(Column)

print(best_distributions)

# Create a summary table of number of times each distribution is best
best_dist_summary <- best_distributions %>%
    count(Distribution) %>%
    arrange(desc(n))

print(best_dist_summary)

results_with_params <- analyze_all_columns_with_params(exampleNoiseData, n_tail = 50, isCensored)

# Create parameter summary
param_summary <- create_parameter_summary(results_with_params$metrics)

# Save results with parameters
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
write.csv(results_with_params$metrics, 
          paste0("distribution_fit_metrics_with_params_", timestamp, ".csv"), 
          row.names = FALSE)
write.csv(param_summary, 
          paste0("best_fit_parameter_summary_", timestamp, ".csv"), 
          row.names = FALSE)

# Print parameter summary
cat("\n=== Parameter Summary for Best Fitting Distributions ===\n")
print(param_summary)
