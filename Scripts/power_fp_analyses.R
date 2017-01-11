# Lia Harrington 2016 - Community Detection
# power_fp_analyses.R
#
# Usage:
# Run in command line:
#
#       Rscript power_fp_analyses.R
#
# Output:
# Creates paper figures of power and false positives 

data <- read.csv("full_data.csv")

library(ggplot2)
library(reshape)
library(cowplot)
library(Hmisc)

d <- data[c("iter_num", "method", "num_paths", "percent_path", "percent_addit",
           "true_positive", "false_positive", "false_negative")]

melt_df <- melt(d, measure.vars = c("true_positive", "false_positive", 
                                    "false_negative"))

colnames(melt_df)[6] <- c("metric") # renames "variable" column to "metric"

# plot of average TP, FP, and FN for each method 
ggplot(melt_df, aes(x = metric, y = value, fill = factor(method))) +
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               fun.args = list(conf.int = .95), width = .3, 
               position = position_dodge(.9)) +
  xlab("Metric") + ylab("Average Frequency") + scale_fill_hue(name = "Method")

ggsave("average_metrics.png", path = "./Paper_Figures")

# plot of false positives rate across path number for each method
ggplot(data, aes(x = num_paths, y = false_positive / (false_positive + true_negative), 
  fill = factor(method))) +
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               fun.args = list(conf.int = .95), width = .3, 
               position = position_dodge(.9)) +
  xlab("Number of Paths") + ylab("False Positives Rate") + 
  scale_fill_hue(name = "Method")

ggsave("fpr_paths.png", path = "./Paper_Figures")

# plot of power across path number for each method 
ggplot(data, aes(x = num_paths, y = true_positive / (true_positive + false_negative), 
  fill = factor(method))) +
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", width = .75) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               fun.args = list(conf.int = .95), width = .3, 
               position = position_dodge(.75)) +
  xlab("Number of Paths") + ylab("Power") + scale_fill_hue(name = "Method")

ggsave("power_paths.png", path = "./Paper_Figures")

# plot of false positives rate across percent path for each method
ggplot(data, aes(x = percent_path, y = false_positive / (false_positive + true_negative), 
  fill = factor(method))) +
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", width = .1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
              fun.args=list(conf.int = .95), width = .1, 
              position = position_dodge(.1)) +
  xlab("Percent Path") + ylab("False Positives Rate") + 
  scale_fill_hue(name="Method")

ggsave("fpr_percen_path.png", path = "./Paper_Figures")

# plot of power across percent path for each method 
ggplot(data, aes(x = percent_path, y = true_positive / (true_positive + false_negative), 
  fill = factor(method))) +
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", width = .1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               fun.args = list(conf.int = .95), width = .1, 
               position = position_dodge(.1)) +
  xlab("Percentage of Path") + ylab("Power") + scale_fill_hue(name = "Method")

ggsave("power_percen_path.png", path = "./Paper_Figures")

# plot of false positives rate across percent addit for each method
ggplot(data, aes(x = percent_addit, y = false_positive / (false_positive + true_negative), 
  fill = factor(method))) +
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", width = .1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               fun.args = list(conf.int = .95), width = .1, 
               position = position_dodge(.1)) +
  xlab("Percent Additional") + ylab("False Positives Rate") + 
  scale_fill_hue(name = "Method")

ggsave("fpr_percen_addit.png", path = "./Paper_Figures")

# plot of power across percent addit for each method 
ggplot(data, aes(x = percent_addit, y = true_positive / (true_positive + false_negative), 
  fill = factor(method))) +
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", width = .1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               fun.args = list(conf.int = .95), width = .1, 
               position = position_dodge(.1)) +
  xlab("Percentage Additional") + ylab("Power") + scale_fill_hue(name = "Method")

ggsave("power_percen_addit.png", path = "./Paper_Figures")
