# Lia Harrington 2016 - Community Detection
# gea_paper_figures.R
#
# Usage:
# Run in command line:
#
#       Rscript Scripts/gea_paper_figures.R
#
# Output:
# Produces boxplots of various metrics as displayed in paper.  

library(cowplot)
library(dplyr)
library(readr)
library(ggplot2)

data_file <- file.path("Data", "all_iterations_data.csv")
iterations_df <- readr::read_csv("all_iterations_data.csv")

iterations_df <- dplyr::mutate(iterations_df, 
                               precision = (true_positive / (true_positive + false_positive)))

iterations_df <- dplyr::mutate(iterations_df, 
                               recall = (true_positive / (true_positive + false_negative)))

iterations_df <- dplyr::mutate(iterations_df, 
                               f1_score = (2 * ( (precision * recall) / (precision + recall))))

iterations_df <- dplyr::mutate(iterations_df, 
                               fpr = (false_positive / (false_positive + true_negative)))

iterations_df <- dplyr::mutate(iterations_df, 
                               fnr = (false_negative / (true_positive + false_negative)))

plot_f1 <- function(m, a_min, a_max, iterations_df){
  # m = number of paths
  # a_min = min percent additional
  # a_max = maximum percent additional
  # iterations_df = mutated data frame
  
  path_sub <- iterations_df$num_paths == m
  add_sub <- (iterations_df$percent_addit == a_min | iterations_df$percent_addit == a_max)
  plot_ready <- iterations_df[path_sub & add_sub, ]
  
  p <- ggplot(plot_ready, aes(method, f1_score)) + geom_boxplot(aes(fill = method)) +
    facet_grid(percent_addit ~ percent_path) + ylab("F1 Score") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    xlab(" ") + scale_fill_discrete(name = "Method")
  
  figure_name <- file.path("Paper_Figs", paste0("f1_boxplots_", m, ".png"))
  ggsave(filename = figure_name, plot = p)
}

# plot F1 for m = 2-8 using either 10% or 100% additional random genes for 
# percent pathway = .3, .475, .65, .825 and 1.0

for (m in 2:8){
  plot_f1(m, .1, 1, iterations_df)
} 
  
plot_fnr <- function(m, a_min, a_max, iterations_df){
  # plots the false negative rate
  
  path_sub <- iterations_df$num_paths == m
  add_sub <- (iterations_df$percent_addit == a_min | iterations_df$percent_addit == a_max)
  plot_ready <- iterations_df[path_sub & add_sub, ]
  
  p <- ggplot(plot_ready, aes(method, fnr)) + geom_boxplot(aes(fill = method)) +
    facet_grid(percent_addit ~ percent_path) + ylab("False Negative Rate") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    xlab(" ") + scale_fill_discrete(name = "Method")
  
  figure_name <- file.path("Paper_Figs", paste0("fnr_boxplots_", m, ".png"))
  ggsave(filename = figure_name, plot = p)
}   
  
plot_fnr(4, .1, 1, iterations_df)  
  
plot_fpr <- function(m, a_min, a_max, iterations_df){
  # plots false positive rate
  
  path_sub <- iterations_df$num_paths == m
  add_sub <- (iterations_df$percent_addit == a_min | iterations_df$percent_addit == a_max)
  plot_ready <- iterations_df[path_sub & add_sub, ]
  
  p <- ggplot(plot_ready, aes(method, fpr)) + geom_boxplot(aes(fill = method)) +
    facet_grid(percent_addit ~ percent_path) + ylab("False Positive Rate") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    xlab(" ") + scale_fill_discrete(name = "Method")
  
  figure_name <- file.path("Paper_Figs", paste0("fpr_boxplots_", m, ".png"))
  ggsave(filename = figure_name, plot = p)
}  
  
plot_fpr(4, .1, 1, iterations_df)
