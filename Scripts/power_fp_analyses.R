data = read.csv('full_data.csv')
attach(data)

library(ggplot2)
library(reshape)
library(cowplot)

d = data[,c(1:7, 9)]
melt_df <- melt(d, measure.vars =c('true_positive', 'false_positive', 
                                    'false_negative'))
colnames(melt_df)[6] = c('metric')

# > head(melt_df)
# iter_num  method num_paths percent_path percent_addit   metric         value
# 1        0 ctr_all         3          0.3           0.1 true_positive     3
# 2        1 ctr_all         3          0.3           0.1 true_positive     3
# 3        2 ctr_all         3          0.3           0.1 true_positive     3

# plot of average TP, FP, and FN for each method 
ggplot(melt_df, aes(metric, value, fill=factor(method)))+
                    geom_bar(stat = 'summary', fun.y='mean', position="dodge")+ 
                    stat_summary(fun.data=mean_cl_boot, geom='errorbar', 
                    fun.args=list(conf.int=.95),
                    width=.3, position = position_dodge(.9))+
                    xlab('Metric')+
                    ylab('Average Frequency')+
                    scale_fill_hue(name='Method')

# plot of false positives rate across path number for each method
ggplot(data, aes(num_paths, false_positive/(false_positive+true_negative), 
                    fill=factor(method)))+
                    geom_bar(stat = 'summary', fun.y='mean', position="dodge")+ 
                    stat_summary(fun.data=mean_cl_boot, geom='errorbar', 
                    fun.args=list(conf.int=.95),
                    width=.3, position = position_dodge(.9))+
                    xlab('Number of Paths')+
                    ylab('False Positives Rate')+
                    scale_fill_hue(name='Method')

# plot of power across path number for each method 
ggplot(data, aes(num_paths, true_positive/(true_positive+false_negative), 
                 fill=factor(method)))+
                 geom_bar(stat = 'summary', fun.y='mean', position="dodge", width=.75)+ 
                 stat_summary(fun.data=mean_cl_boot, geom='errorbar', 
                 fun.args=list(conf.int=.95),
                 width=.3, position = position_dodge(.75))+
                 xlab('Number of Paths')+
                 ylab('Power')+
                 scale_fill_hue(name='Method')

# plot of false positives rate across percent path for each method
ggplot(data, aes(percent_path, false_positive/(false_positive+true_negative), 
                 fill=factor(method)))+
                 geom_bar(stat = 'summary', fun.y='mean', position="dodge", width=.1)+ 
                 stat_summary(fun.data=mean_cl_boot, geom='errorbar', 
                 fun.args=list(conf.int=.95),
                 width=.1, position = position_dodge(.1))+
                 xlab('Percent Path')+
                 ylab('False Positives Rate')+
                 scale_fill_hue(name='Method')

# plot of power across percent path for each method 
ggplot(data, aes(percent_path, true_positive/(true_positive+false_negative), 
                 fill=factor(method)))+
                 geom_bar(stat = 'summary', fun.y='mean', position="dodge", width=.1)+ 
                 stat_summary(fun.data=mean_cl_boot, geom='errorbar', 
                 fun.args=list(conf.int=.95),
                 width=.1, position = position_dodge(.1))+
                 xlab('Percentage of Path')+
                 ylab('Power')+
                 scale_fill_hue(name='Method')

# plot of false positives rate across percent addit for each method
ggplot(data, aes(percent_addit, false_positive/(false_positive+true_negative), 
                 fill=factor(method)))+
                 geom_bar(stat = 'summary', fun.y='mean', position="dodge", width=.1)+ 
                 stat_summary(fun.data=mean_cl_boot, geom='errorbar', 
                 fun.args=list(conf.int=.95),
                 width=.1, position = position_dodge(.1))+
                 xlab('Percent Additional')+
                 ylab('False Positives Rate')+
                 scale_fill_hue(name='Method')

# plot of power across percent addit for each method 
ggplot(data, aes(percent_addit, true_positive/(true_positive+false_negative), 
                 fill=factor(method)))+
                 geom_bar(stat = 'summary', fun.y='mean', position="dodge", width=.1)+ 
                 stat_summary(fun.data=mean_cl_boot, geom='errorbar', 
                 fun.args=list(conf.int=.95),
                 width=.1, position = position_dodge(.1))+
                 xlab('Percentage Additional')+
                 ylab('Power')+
                 scale_fill_hue(name='Method')
