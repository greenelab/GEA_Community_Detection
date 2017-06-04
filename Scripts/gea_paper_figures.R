# Lia Harrington 2016 - Community Detection
# gea_paper_figures.R
#
# Usage:
# Run in command line:
#
#       Rscript gea_paper_figures.R
#
# Output:
# Produces boxplots of various metrics as displayed in paper.  

library(cowplot)
library(dplyr)
library(ggplot2)

setwd("~/Dropbox/GEA_Community_Detection-master/Scripts")
data <- read.csv("all_iterations_data.csv")

data <- dplyr::mutate(data, precision = (true_positive / (true_positive + false_positive)))
data <- dplyr::mutate(data, recall = (true_positive / (true_positive + false_negative)))
data <- dplyr::mutate(data, f1_score = (2 * ((precision * recall)/(precision + recall))))

data <- dplyr::mutate(data, fpr = (false_positive / (false_positive + true_negative)))
data <- dplyr::mutate(data, fnr = (false_negative / (true_positive + false_negative)))


plot_f1 <- function(m, a_min, a_max, data){
  # m = number of paths
  # a_min = min percent additional
  # a_max = maximum percent additional
  # df = mutated data frame 
  
  plot_ready <- data[data$num_paths == m & (data$percent_addit == a_min | data$percent_addit == a_max), ]
  
  ggplot(plot_ready, aes(method, f1_score)) + geom_boxplot(aes(fill = method)) +
    facet_grid(percent_addit ~ percent_path)+ylab('F1 Score') +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    xlab(' ')+scale_fill_discrete(name = "Method")
  
  ggsave(sprintf("f1_boxplots%s.png", m), path = "/Users/lia/Dropbox/GEA_Community_Detection-master/Paper_Figs")
}

# plot F1 for m = 2-8 using either 10% or 100% additional random genes for percent pathway = .3, .475, .65, .825 
# and 1.0
for (m in 2:8){
  plot_f1(m, .1, 1, data)
} 
  
plot_fnr <- function(m, a_min, a_max, data){
  # plots the false negative rate
  
  plot_ready <- data[data$num_paths == m & (data$percent_addit == a_min | data$percent_addit == a_max), ]
  
  ggplot(plot_ready, aes(method, fnr)) + geom_boxplot(aes(fill = method)) +
    facet_grid(percent_addit ~ percent_path)+ylab('False Negative Rate') +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    xlab(' ')+scale_fill_discrete(name = "Method")
  
  ggsave(sprintf("fnr_boxplots%s.png", m), path = "/Users/lia/Dropbox/GEA_Community_Detection-master/Paper_Figs")
  
}   
  
plot_fnr(4, .1, 1, data)  
  
plot_fpr <- function(a, a_min, a_max, data){
  # plots false positive rate
  
  ggplot(plot_ready, aes(method, fpr)) + geom_boxplot(aes(fill = method)) +
    facet_grid(percent_addit ~ percent_path)+ylab('False Positive Rate') +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    xlab(' ')+scale_fill_discrete(name = "Method")
  
  ggsave(sprintf("fpr_boxplots%s.png",m), path = "/Users/lia/Dropbox/GEA_Community_Detection-master/Paper_Figs")
  
}  
  
plot_fpr(4, .1, 1, data)
